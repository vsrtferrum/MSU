package main

import (
	"fmt"
	"math"
	"math/rand"
	"msu/infra"
	"os"
	"time"

	"golang.org/x/crypto/ssh/terminal"
)

var (
	r1        = 1000.0
	r2        = 2000.0
	c1        = 0.001
	c2        = 0.0005
	x1        = 0.0
	x2        = 0.0
	dt        = 0.001
	tmax      = 10.0
	percent   = 0.2
	percent2  = 0.15
	tT        = 2.18
	separator string
)

func SystemIteration(R1, R2, C1, C2 float64, t, u float64) float64 {
	x1 = ((-(R1+R2)*x1/(C1*R1*R2))+x2/(C1*R2)+u/(C1*R1))*dt + x1
	x2 = (x1/(C2*R2)-x2/(C2*R2))*dt + x2
	return -1*x1 + u
}

func plusMin(flag rune, val float64) float64 {
	if flag == '-' {
		return -1 * val
	}
	return val
}

func RandomFloat64(min, max float64) float64 {
	rand.Seed(time.Now().UnixNano())
	return min + rand.Float64()*(max-min)
}

func numToIndex(v int) int {
	switch v {
	case 0:
		return 0
	case 1:
		return 1
	case 2:
		return 2
	case 3:
		return 3
	case 4:
		return 12
	case 5:
		return 13
	case 6:
		return 23
	default:
		return 123
	}
}

func indexToNum(v string) int {
	switch v {
	case "+++":
		return 0
	case "++-":
		return 1
	case "+-+":
		return 2
	case "+--":
		return 3
	case "-++":
		return 4
	case "-+-":
		return 5
	case "--+":
		return 6
	default:
		return 7
	}
}
func MinMaxFactor() {
	fmt.Println(separator)
	fmt.Printf("C1 от %f до %f, при изначальном %f\n", c1*(1-float64(percent)), c1*(1+float64(percent)), c1)
	fmt.Printf("C2 от %f до %f, при изначальном %f\n", c2*(1-float64(percent)), c2*(1+float64(percent)), c2)
	fmt.Printf("R1 от %f до %f, при изначальном %f\n", r1*(1-float64(percent)), r1*(1+float64(percent)), r1)
	fmt.Printf("R2 от %f до %f, при изначальном %f\n", r2*(1-float64(percent)), r2*(1+float64(percent)), r2)
	fmt.Println(separator)
}
func mid(input []float64) float64 {
	sum := 0.0
	for _, i := range input {
		sum += i
	}
	return sum / float64(len(input))
}
func removeFirstCharacter(s string) string {
	if len(s) > 0 {
		return s[1:]
	}
	return s
}

func skoAndDispersion(data []float64, mean float64) (float64, float64) {
	disp := 0.0
	for _, val := range data {
		disp += (val - mean) * (val - mean) / float64(len(data)-1)
	}
	return math.Sqrt(disp), disp

}

func LogicalAnd(a, b rune) rune {
	if (a == '-' && b == '-') || (a == '+' && b == '+') {
		return '+'
	}
	return '-'
}

func main() {
	fd := int(os.Stdout.Fd())
	width, _, _ := terminal.GetSize(fd)
	sep := make([]rune, width)
	for i := range sep {
		sep[i] = '_'
	}

	separator = string(sep)
	res := make([]float64, int(tmax/dt), int(tmax/dt))
	time := make([]float64, int(tmax/dt), int(tmax/dt))

	for i := range res {
		res[i] = SystemIteration(r1, r2, c1, c2, float64(i)*dt, 1)
		time[i] = float64(i) * dt
	}

	infra.DrawChart(time, map[string][]float64{"Переходная характеристика": res}, "Переходная_характеристика.html")

	MinMaxFactor()
	mp := make(map[string][]float64)
	ch := make(chan string, 8)
	ch <- "+"
	ch <- "-"
	for len(ch) != 0 {
		temp := <-ch
		if len(temp) != 3 {
			ch <- string(append([]rune(temp), '+'))
			ch <- string(append([]rune(temp), '-'))
		} else {
			mp[temp] = make([]float64, 0)
		}
	}

	fmt.Println("id\tВремя переходного процесса")
	ProcessTime := make(map[string]float64)
	for str := range mp {
		arr := []float64{r1, c1, c2}
		for i, val := range str {
			if val == '+' {
				arr[i] *= 1 + percent
			} else {
				arr[i] /= 1 - percent
			}
			x1, x2 = 0, 0
			ans := -1.0
			for i := range time {
				temp := SystemIteration(arr[0], r2, arr[1], arr[2], dt*float64(i), 1)
				mp[str] = append(mp[str], temp)
				if math.Abs(temp) <= float64(0.05) && ans == -1.0 {
					ans = dt * float64(i)
				} else if math.Abs(temp) > float64(0.05) {
					ans = -1.0
				}
			}
			ProcessTime[str] = ans
		}
		fmt.Printf("%s\t%f\n", str, ProcessTime[str])
	}

	fmt.Println(separator)
	infra.DrawChart(time, mp, "Переходные_характеристики8.html")
	b := make([]float64, 8)
	fmt.Printf("x0\tx1\tx2\tx3\tx1x2\tx1x3\tx2x3\tx1x2x3\ty\n")
	for iter := range mp {
		fmt.Printf("+\t%c\t%c\t%c\t%c\t%c\t%c\t%c\t%f\n", iter[0], iter[1], iter[2], LogicalAnd(rune(iter[0]), rune(iter[1])), LogicalAnd(rune(iter[0]), rune(iter[2])), LogicalAnd(rune(iter[1]), rune(iter[2])), LogicalAnd(LogicalAnd(rune(iter[0]), rune(iter[1])), rune(iter[2])), ProcessTime[iter])
		b[0] += ProcessTime[iter]
		b[1] += plusMin(rune(iter[0]), ProcessTime[iter])
		b[2] += plusMin(rune(iter[1]), ProcessTime[iter])
		b[3] += plusMin(rune(iter[2]), ProcessTime[iter])
		b[4] += plusMin(LogicalAnd(rune(iter[0]), rune(iter[1])), ProcessTime[iter])
		b[5] += plusMin(LogicalAnd(rune(iter[0]), rune(iter[2])), ProcessTime[iter])
		b[6] += plusMin(LogicalAnd(rune(iter[1]), rune(iter[2])), ProcessTime[iter])
		b[7] += plusMin(LogicalAnd(LogicalAnd(rune(iter[0]), rune(iter[1])), rune(iter[2])), ProcessTime[iter])
	}
	fmt.Println(separator)
	fmt.Println("Коэффициенты факторной математической модели:")
	for i, val := range b {
		fmt.Printf("Коэффициент b%d = %f\n", numToIndex(i), val/8)
	}
	fmt.Println(separator)
	mp2 := make(map[string][]float64)
	merger := make(map[string][]float64)
	for i := 1; i < 4; i++ {
		temp := []float64{r1, c1, c2}
		fmt.Printf("Опыт №%d\n", i)
		fmt.Printf("Вариация\tx1\t\tx2\t\tx3\t\ty\n")
		for iter := range mp {
			for i, val := range iter {
				if val == '+' {
					temp[i] = RandomFloat64(temp[i], temp[i]*(1+percent2))
				} else {
					temp[i] = RandomFloat64(temp[i]*(1-percent2), temp[i])
				}
			}
			x1, x2 = 0, 0
			ans := -1.0
			for j := range time {
				temp := SystemIteration(temp[0], r2, temp[1], temp[2], dt*float64(i), 1)
				mp2[fmt.Sprintf("%d%s", i, iter)] = append(mp2[fmt.Sprintf("%d%s", i, iter)], temp)
				if math.Abs(temp) <= float64(0.05) && ans == -1.0 {
					ans = dt * float64(j)
				} else if math.Abs(temp) > float64(0.05) {
					ans = -1.0
				}
			}
			fmt.Printf("%s\t%f\t%f\t%f\t%f\n", iter, temp[0], temp[1], temp[2], ans)
			merger[iter] = append(merger[iter], ans)
		}
	}
	fmt.Println(separator)
	infra.DrawChart(time, mp2, "Переходные_характеристики24.html")
	fmt.Printf("x0\tx1\tx2\tx3\tx1x2\tx1x3\tx2x3\tx1x2x3\ty1\t\ty2\t\ty3\t\tmidy\n")

	for iter := range merger {
		temp := merger[iter]
		fmt.Printf("+\t%c\t%c\t%c\t%c\t%c\t%c\t%c\t%f\t%f\t%f\t%f\n", iter[0], iter[1], iter[2], LogicalAnd(rune(iter[0]), rune(iter[1])), LogicalAnd(rune(iter[0]), rune(iter[2])), LogicalAnd(rune(iter[1]), rune(iter[2])), LogicalAnd(LogicalAnd(rune(iter[0]), rune(iter[1])), rune(iter[2])), temp[0], temp[1], temp[2], (temp[0]+temp[1]+temp[2])/3)
	}
	fmt.Println(separator)
	fmt.Println("Рассчёт СКО и Дисперсии")
	fmt.Println("Вариация\tСКО\t\tДисперсия")
	for iter, val := range merger {

		mean := 0.0
		for _, v := range val {
			mean += v
		}
		mean /= float64(len(val))
		sko, disp := skoAndDispersion(val, mean)
		fmt.Printf("%s\t\t%f\t%f\n", iter, sko, disp)
	}
	fmt.Println(separator)
	var mnindex, mxindex string
	mn, mx := 11.0, -1.0
	mnSKO, mxSKO := 0.0, 0.0
	mnDIFF, mxDIFF := 0.0, 0.0
	for index, val := range merger {
		temp := 0.0
		for _, v := range val {
			temp += v
		}
		mn = min(mn, temp/float64(len(val)))
		mx = max(mx, temp/float64(len(val)))
		if mx == temp/float64(len(val)) {
			mxSKO, mxDIFF = skoAndDispersion(val, temp/float64(len(val)))
			mxindex = index
		}
		if mn == temp/float64(len(val)) {
			mnSKO, mnDIFF = skoAndDispersion(val, temp/float64(len(val)))
			mnindex = index
		}
	}
	fmt.Printf("Метод вычисления максимального относительного отклонения\n\n")
	fmt.Printf("Среднее значение по максимальному времени переходного процесса: %f\n", mx)
	fmt.Printf("СКО по максимальным значениям: %f\n", mxSKO)
	fmt.Printf("Дисперсия по максимальным значениям: %f\n\n", mxDIFF)

	fmt.Printf("Среднее значение по минимальному времени переходного процесса: %f\n", mn)
	fmt.Printf("СКО по минимальному значениям: %f\n", mnSKO)
	fmt.Printf("Дисперсия по минимальному значениям: %f\n", mnDIFF)
	fmt.Println(separator)
	for i := 0; i < 3; i++ {
		fmt.Printf("Максимальное значение среди всех эксперимента №%d\n", i+1)
		fmt.Printf("Максимальное значение\ttau\t\n")
		tmn, tmx := merger[mnindex], merger[mxindex]
		fmt.Printf("%f\t\t%f\n", tmx[i], math.Abs(tmx[i]-mx)/mxSKO)
		fmt.Printf("Минимальное значение\ttau\t\n")
		fmt.Printf("%f\t\t%f\n", tmn[i], math.Abs(tmn[i]-mn)/mnSKO)
	}
	fmt.Println(separator)
	fmt.Println("Проверка гипотезы об однородности дисперссии опытов:")
	fmt.Print("\n")
	var str []string
	for iter := range merger {
		str = append(str, iter)
		fmt.Printf("\t%s\t", iter)
	}
	fmt.Print("\n")
	for iter, val := range merger {
		fmt.Printf("%s\t", iter)
		for _, iter2 := range str {
			temp, temp2 := 0.0, 0.0
			for _, v := range val {
				temp += v
			}
			for _, v2 := range merger[iter2] {
				temp2 += v2
			}
			_, disp := skoAndDispersion(val, temp/float64(len(val)))
			_, disp2 := skoAndDispersion(merger[iter2], temp2/float64(len(merger[iter2])))
			fmt.Printf("%f\t", max(disp, disp2)/min(disp, disp2))
		}
		fmt.Print("\n")
	}
	fmt.Println(separator)
	fmt.Println("Проверка однородности G-критерием Кохрена ")
	dispsum, dispmax := 0.0, 0.0
	for _, val := range merger {
		temp := 0.0
		for _, v := range val {
			temp += v
		}
		_, disp := skoAndDispersion(val, temp/float64(len(val)))
		dispsum += disp
		dispmax = max(dispmax, disp)
	}
	fmt.Printf("Gp = %f\n", dispmax/dispsum)
	fmt.Println(separator)
	fmt.Println("Дисперсия воспроизводимости")
	fmt.Printf("%f\n", dispsum/8)
	fmt.Println(separator)

	for i := range b {
		b[i] = 1 / float64(len(b))
	}
	for iter := range mp {
		fmt.Printf("+\t%c\t%c\t%c\t%c\t%c\t%c\t%c\t%f\n", iter[0], iter[1], iter[2], LogicalAnd(rune(iter[0]), rune(iter[1])), LogicalAnd(rune(iter[0]), rune(iter[2])), LogicalAnd(rune(iter[1]), rune(iter[2])), LogicalAnd(LogicalAnd(rune(iter[0]), rune(iter[1])), rune(iter[2])), ProcessTime[iter])
		b[0] += mid(merger[iter])
		b[1] += plusMin(rune(iter[0]), mid(merger[iter]))
		b[2] += plusMin(rune(iter[1]), mid(merger[iter]))
		b[3] += plusMin(rune(iter[2]), mid(merger[iter]))
		b[4] += plusMin(LogicalAnd(rune(iter[0]), rune(iter[1])), mid(merger[iter]))
		b[5] += plusMin(LogicalAnd(rune(iter[0]), rune(iter[2])), mid(merger[iter]))
		b[6] += plusMin(LogicalAnd(rune(iter[1]), rune(iter[2])), mid(merger[iter]))
		b[7] += plusMin(LogicalAnd(LogicalAnd(rune(iter[0]), rune(iter[1])), rune(iter[2])), mid(merger[iter]))
	}
	fmt.Println(separator)
	fmt.Println("Коэффициенты факторной математической модели:")
	for iter, val := range b {
		fmt.Printf("b%d = %f\n", numToIndex(iter), val)
	}
	fmt.Println(separator)
	fmt.Println("Сравнение абсолютных величин с доврительным интервалом")
	fmt.Printf("Вариация\tb_t\t\tb\n")
	ok := make(map[string]bool)
	for iter, val := range merger {
		_, disp := skoAndDispersion(val, mid(val))
		fmt.Printf("%s\t\t%f\t%f\n", iter, math.Sqrt(disp/24)*tT, math.Abs(b[indexToNum(iter)]))
		ok[iter] = math.Sqrt(disp/24)*tT < math.Abs(b[indexToNum(iter)])
	}
	fmt.Println(separator)
	fmt.Println("t-критерий Стьюдента")
	fmt.Printf("Вариация\ttp\t\ttT\n")

	for iter, val := range merger {
		_, disp := skoAndDispersion(val, mid(val))
		fmt.Printf("%s\t\t%f\t%f\n", iter, math.Abs(b[indexToNum(iter)]/math.Sqrt(disp/24)), tT)
		if ok[iter] != (math.Abs(b[indexToNum(iter)]/math.Sqrt(disp/24)) > tT) {
			panic("АдЕкВаТнОсТь НаРуШеНа")
		}
	}
	fmt.Println(separator)
	fmt.Println("Проверка гипотезы об адекватности полученной модели:")

	sum, sumd := 0.0, 0.0
	var ln float64
	for _, val := range ok {
		if val {
			ln += 1.0
		}
	}
	for iter := range mp2 {
		r := 0.0
		if ok[removeFirstCharacter(iter)] {
			for _, v := range iter {
				if v == '+' {
					r += b[indexToNum(iter)]
				} else if v == '-' {
					r -= b[indexToNum(iter)]
				}
			}

			_, disp := skoAndDispersion(merger[removeFirstCharacter(iter)], mid(merger[removeFirstCharacter(iter)]))
			sum += (mid(merger[removeFirstCharacter(iter)]) - r) * (mid(merger[removeFirstCharacter(iter)]) - r)
			sumd += disp / float64(ln)
		}
	}
	fmt.Printf("Sad^2 = %f\n", 1/sum*3.0/4.0)
	fmt.Printf("Fp = %f Ft =%f", 3.0/4.0*sumd/sumd, 4.5)
}
