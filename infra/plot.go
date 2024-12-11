package infra

import (
	"fmt"
	"os"

	"github.com/go-echarts/go-echarts/v2/charts"
	"github.com/go-echarts/go-echarts/v2/opts"
)

func generateXAxis(xValues []float64) []string {
	xAxis := make([]string, len(xValues))
	for i, v := range xValues {
		xAxis[i] = fmt.Sprintf("%.2f", v)
	}
	return xAxis
}

func DrawChart(xValues []float64, seriesData map[string][]float64, fileName string) {
	lineChart := charts.NewLine()
	lineChart.SetGlobalOptions()

	for seriesName, data := range seriesData {
		lineData := make([]opts.LineData, len(data))
		for i, v := range data {
			lineData[i] = opts.LineData{Value: v}
		}
		lineChart.AddSeries(seriesName, lineData)
	}
	lineChart.SetXAxis(generateXAxis(xValues))
	f, err := os.Create(fileName)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	if err := lineChart.Render(f); err != nil {
		panic(err)
	}
}
