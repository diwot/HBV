using HBV.InputOutput;
using HBV;
using ScottPlot;

namespace HBV
{
    public partial class MainForm : Form
    {
        public MainForm()
        {
            InitializeComponent();
        }

        private void butLoadAndPlotData_Click(object sender, EventArgs e)
        {
            CsvDataMeteoData meteoData = new CsvDataMeteoData(@"C:\Users\twidmer\Downloads\N\meteodata_24h.csv");


            List<DateTime> dateTimes = meteoData.allDateTimes;
            List<float> tground = meteoData.Tground;

            PlotTimeSeries(plotRainfall, dateTimes, tground, "Ground Temperature", "Date", "Temperature");
        }

        private void PlotTimeSeries(FormsPlot plot, IList<DateTime> times, IList<float> values, string title, string xLabel, string yLabel)
        {
            // Convert DateTime to double array of days since first date
            double[] xData = times.Select(t => (t - times[0]).TotalDays).ToArray();

            // Convert float list to double array
            double[] yData = values.Select(v => (double)v).ToArray();

            // Create plot
            plot.Plot.Clear();
            plot.Plot.AddScatter(xData, yData);

            // Customize the plot
            plot.Plot.Title(title);
            plot.Plot.XLabel(xLabel);
            plot.Plot.YLabel(yLabel);

            // Format x-axis as dates
            plot.Plot.XAxis.DateTimeFormat(true);
            plot.Plot.XAxis.TickLabelFormat("MM/dd/yyyy", true);

            // Update display
            plot.Dock = DockStyle.Fill;
            plot.Refresh();
        }

        private void butRunOptimizer_Click(object sender, EventArgs e)
        {
            PSODriver.RunOptimizer();
        }

        private void butRunModel_Click(object sender, EventArgs e)
        {
            (ElevationBandData result, CsvDataMeteoData data) = PSODriver.Run();
            PlotTimeSeries(plotRainfall, data.allDateTimes, result.r, "rainfall", "time", "mm");
            PlotTimeSeries(plotRiverDischarge, data.allDateTimes, result.qr, "river discharge", "time", "mm");
        }
    }
}
