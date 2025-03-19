using HBV.InputOutput;
using HBV;
using ScottPlot;
using System.Diagnostics;
using ScottPlot.Drawing.Colormaps;

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
            var xData = times.Select(t => (float)(t - times[0]).TotalDays).ToList();
            PlotTimeSeries(plot, xData, values, title, xLabel, yLabel, scalingFactors);
        }

        private void PlotTimeSeries(FormsPlot plot, IList<float> times, List<IList<float>> values, string title, string xLabel, string yLabel, List<float> scalingFactors = null)
        {
            // Convert DateTime to double array of days since first date
           

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
                plot.Plot.AddScatter(times.Select(x => (double)x).ToArray(), yData);
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

            (Task t, PSO optimizer) = PSODriver.RunOptimizer(meteoData, catchmentInfo, ModelType.NormalHBV);
            this.optimizer = optimizer;
            t.Start();
            tmr.Enabled = true;
        }

        private void butRunModel_Click(object sender, EventArgs e)
        {

            HBVParams[] pars = null;
            if (dlgOpen.ShowDialog() == DialogResult.OK)
                pars = ResultWriter.Load(dlgOpen.FileName).ToArray();

            // Load meteorological data
            CsvMeteoData meteoData = new CsvMeteoData(@"C:\Daten\meteodata_24h.csv");

            // Load catchment info
            var catchmentInfo = new CsvDataElevationBands(@"C:\Daten\Elevation_bands.csv");
            (ElevationBandData result, CsvMeteoData data) = PSODriver.Run(meteoData, catchmentInfo, pars);
            PlotTimeSeries(plotRainfall, data.allDateTimes, result.r, "rainfall", "time", "mm");
            PlotTimeSeries(plotRiverDischarge, data.allDateTimes, result.qr, "river discharge", "time", "mm");
        }

        private void tmr_Tick(object sender, EventArgs e)
        {
            float[] rawData = optimizer.GetCopyGlobalBest();

            (ElevationBandData result, CsvMeteoData data) = PSODriver.Run(meteoData, catchmentInfo, PSODriver.FloatToHBVParams(rawData, 8, 20));
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

        private void butComparisonPlot_Click(object sender, EventArgs e)
        {

            string dir = Application.StartupPath;
            string pathMatlabReference = @"C:\Daten\test_discharge_lumped_ice_storage_24h.mat";

            List<IList<float>> test = new List<IList<float>>();

            using (FileStream fs = new FileStream(pathMatlabReference, FileMode.Open))
            {
                MatFileHandler.MatFileReader r = new MatFileHandler.MatFileReader(fs);
                var file = r.Read();

                var vars = file.Variables.ToList();
                var testName = vars[0].Name;
                List<float> testValue = ToFloatList(vars[0].Value.ConvertToDoubleArray());
                test.Add(testValue);
            }





            Console.WriteLine(dir);

            //C:\git\HBV\HBV\Data
            //C:\git\HBV\HBV\bin\Debug\net8.0-windows\

            //string path = dir + "../../../Data/Discharge_12h_lumped.csv";


            //List<IList<float>> test = CsvLoader.LoadCsv(path);
            List<float> x = GenerateIndexer(test[0].Count);



            string pathPars = @"C:\Daten\test_params_lumped_ice_storage_24h.csv";

            List<float> qr;
            {
                

                var csvPars = CsvLoader.LoadCsv(pathPars);
                HBVParams[] pars = new HBVParams[csvPars.Count];
                for (int i = 0;i<pars.Length;++i)
                   pars[i] = new HBVParams(csvPars[i]);
                ModelType perma = ModelType.IceStorage;

                // Load meteorological data
                CsvMeteoData meteoData = new CsvMeteoData(@"C:\Daten\meteodata_24h.csv");

                //int start = meteoData.allDateTimes.IndexOf(new DateTime(2013, 1, 1, 0, 0, 0));
                //meteoData = meteoData.GetSubRange(start, x.Count);

                // Load catchment info
                var catchmentInfo = new CsvDataElevationBands(@"C:\Daten\Elevation_bands_lumped.csv");
                (ElevationBandData result, CsvMeteoData data) = PSODriver.Run(meteoData, catchmentInfo, pars, perma);

                qr = result.qr.ToList();

                //test.Add(qr);
            }







            PlotTimeSeries(plotRainfall, x, test, "Compare", "Index", "Matlab");
        }

        private List<float> ToFloatList(double[] doubles)
        {
            List<float> result = new List<float>(doubles.Length);
            for (int i = 0; i < doubles.Length; i++)
            {
                result.Add((float)doubles[i]);
            }
            return result;
        }

        private List<float> GenerateIndexer(int count)
        {
            List<float> result = new List<float>(count);
            for (int i = 0; i < count; i++)
                result.Add(i);
            return result;
        }
    }
}
