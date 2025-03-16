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
            //CsvDataMeteoData meteoData = new CsvDataMeteoData(@"C:\Users\twidmer\Downloads\N\meteodata_24h.csv");


            //List<DateTime> dateTimes = meteoData.allDateTimes;
            //List<float> tground = meteoData.Tground;

            //PlotTimeSeries(plotRainfall, dateTimes, tground, "Ground Temperature", "Date", "Temperature");
        }

        private void PlotTimeSeries(FormsPlot plot, IList<DateTime> times, IList<float> values, string title, string xLabel, string yLabel)
        {
            PlotTimeSeries(plot, times, new List<IList<float>> { values }, title, xLabel, yLabel);
        }

        private void PlotTimeSeries(FormsPlot plot, IList<DateTime> times, List<IList<float>> values, string title, string xLabel, string yLabel, List<float> scalingFactors = null)
        {
            // Convert DateTime to double array of days since first date
            double[] xData = times.Select(t => (t - times[0]).TotalDays).ToArray();

            // Create plot
            plot.Plot.Clear();

            // Add each series
            for (int i = 0; i < values.Count; ++i) // (var series in values)
            {
                var series = values[i];
                double scaling = 1.0;
                if (scalingFactors != null)
                    scaling = scalingFactors[i];

                double[] yData = series.Select(v => scaling * (double.IsNaN((double)v) ? 0.0 : (double)v)).ToArray();
                plot.Plot.AddScatter(xData, yData);
            }

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
        CsvMeteoData meteoData;
        CsvDataElevationBands catchmentInfo;

        private PSO optimizer = null;
        private void butRunOptimizer_Click(object sender, EventArgs e)
        {
            // Load meteorological data
            meteoData = new CsvMeteoData(@"C:\Daten\meteodata_24h.csv");

            // Load catchment info
            catchmentInfo = new CsvDataElevationBands(@"C:\Daten\Elevation_bands.csv");

            (Task t, PSO optimizer) = PSODriver.RunOptimizer(meteoData,catchmentInfo);
            this.optimizer = optimizer;
            t.Start();
            tmr.Enabled = true;
        }

        private void butRunModel_Click(object sender, EventArgs e)
        {

            HBVParams[] pars = null;
            if(dlgOpen.ShowDialog()==DialogResult.OK)
                pars = ResultWriter.Load(dlgOpen.FileName).ToArray();

            // Load meteorological data
            CsvMeteoData meteoData = new CsvMeteoData(@"C:\Daten\meteodata_24h.csv");

            // Load catchment info
            var catchmentInfo = new CsvDataElevationBands(@"C:\Daten\Elevation_bands.csv");
            (ElevationBandData result, CsvMeteoData data) = PSODriver.Run(meteoData,catchmentInfo,pars);
            PlotTimeSeries(plotRainfall, data.allDateTimes, result.r, "rainfall", "time", "mm");
            PlotTimeSeries(plotRiverDischarge, data.allDateTimes, result.qr, "river discharge", "time", "mm");
        }

        private void tmr_Tick(object sender, EventArgs e)
        {
            float[] rawData = optimizer.GetCopyGlobalBest();

            (ElevationBandData result, CsvMeteoData data) = PSODriver.Run(meteoData,catchmentInfo,PSODriver.FloatToHBVParams(rawData, 8, 20));
            PlotTimeSeries(plotRainfall, data.allDateTimes, result.r, "rainfall", "time", "mm");
            PlotTimeSeries(plotRiverDischarge, data.allDateTimes, new List<IList<float>> { result.qr, data.Discharge }, "river discharge", "time", "mm", new List<float>() { 1, 0.001f });
        }

        private void butSave_Click(object sender, EventArgs e)
        {
            if (optimizer != null)
            {
                if (dlgSave.ShowDialog() == DialogResult.OK)
                {
                    float[] rawData = optimizer.GetCopyGlobalBest();

                    var r = PSODriver.FloatToHBVParams(rawData, 8, 20);

                    ResultWriter.Save(r, dlgSave.FileName);
                }
            }
            else
                MessageBox.Show("You must run an optimization before anything can be saved.");
        }
    }
}
