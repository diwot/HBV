using HBV;
using HBV.InputOutput;

namespace HBV
{
    public enum ModelType
    {
        NormalHBV,IceStorage,InfiltrationBlocking,BothPermafrostAdaptions
    }
    public class PSODriver
    {
        // Conversion function that can be reused elsewhere
        public static HBVParams[] FloatToHBVParams(float[] optimParams, int numElevationBands, int numParameters)
        {
            HBVParams[] pars = new HBVParams[numElevationBands];
            for (int i = 0; i < numElevationBands; ++i)
            {
                int offset = numParameters * i;
                pars[i] = new HBVParams(
                    optimParams[offset + 0],  // FC
                    optimParams[offset + 1],  // BETA 
                    optimParams[offset + 2],  // LP
                    optimParams[offset + 3],  // ALFA
                    optimParams[offset + 4],  // KF
                    optimParams[offset + 5],  // KS
                    optimParams[offset + 6],  // PERC
                    optimParams[offset + 7],  // CFLUX
                    optimParams[offset + 8],  // TT
                    optimParams[offset + 9],  // TTI
                    optimParams[offset + 10], // CFMAX
                    optimParams[offset + 11], // FOCFMAX
                    optimParams[offset + 12], // CFR
                    optimParams[offset + 13], // WHC
                    optimParams[offset + 14], // MAXBAS
                    optimParams[offset + 15], // precip_corr
                    optimParams[offset + 16], // TTG
                    optimParams[offset + 17], // TTGI
                    optimParams[offset + 18], // corrTg
                    optimParams[offset + 19]  // infil_stop
                );
            }
            return pars;
        }

        public static float[] HBVParamsToFloat(HBVParams[] hbvParams)
        {
            int numParameters = 20;
            float[] result = new float[hbvParams.Length * numParameters];
            for (int i = 0; i < hbvParams.Length; i++)
            {
                int offset = numParameters * i;
                result[offset + 0] = hbvParams[i].FC;
                result[offset + 1] = hbvParams[i].BETA;
                result[offset + 2] = hbvParams[i].LP;
                result[offset + 3] = hbvParams[i].ALFA;
                result[offset + 4] = hbvParams[i].KF;
                result[offset + 5] = hbvParams[i].KS;
                result[offset + 6] = hbvParams[i].PERC;
                result[offset + 7] = hbvParams[i].CFLUX;
                result[offset + 8] = hbvParams[i].TT;
                result[offset + 9] = hbvParams[i].TTI;
                result[offset + 10] = hbvParams[i].CFMAX;
                result[offset + 11] = hbvParams[i].FOCFMAX;
                result[offset + 12] = hbvParams[i].CFR;
                result[offset + 13] = hbvParams[i].WHC;
                result[offset + 14] = hbvParams[i].MAXBAS;
                result[offset + 15] = hbvParams[i].precip_corr;
                result[offset + 16] = hbvParams[i].TTG;
                result[offset + 17] = hbvParams[i].TTGI;
                result[offset + 18] = hbvParams[i].corrTg;
                result[offset + 19] = hbvParams[i].infil_stop;
            }
            return result;
        }

        public static (Task, PSO) RunOptimizer(CsvMeteoData meteoDataMeasured, CsvDataElevationBands catchmentInfo, ModelType perma)
        {
            // Load meteorological data
            // CsvMeteoData meteoDataMeasured = new CsvMeteoData(@"C:\Daten\meteodata_24h.csv");

            // Load catchment info
            // var catchmentInfo = new CsvDataElevationBands(@"C:\Daten\Elevation_bands.csv");

            // Set timestep in seconds
            float seconds = 86400f; // 24h
            // Other options:
            // float seconds = 43200f; // 12h
            // float seconds = 21600f; // 6h  
            // float seconds = 10800f; // 3h
            // float seconds = 3600f;  // 1h

            // Create arrays from the needed data
            var dTime = meteoDataMeasured.allDateTimes;
            var tmean = meteoDataMeasured.Tmean;
            var pet = meteoDataMeasured.PET;
            var precip = meteoDataMeasured.Precipitation_mm_even;
            var discharge = meteoDataMeasured.Discharge.Select(x => x / 1000f).ToArray(); // Convert to m3/s
            var tground = meteoDataMeasured.Tground;

            int numSubcatchments = 8;
            int numParameters = 20;

            int numElevationBands = catchmentInfo.ElevationBottom.Count;


            HBVParams[] pars_test = new HBVParams[numElevationBands];
            for (int i = 0; i < numElevationBands; i++)
            {
                pars_test[i] = new HBVParams(500, 3, 0.5f, 1, 0.1f, 0.01f, 3, 2, 0, 2, 5, 0.5f, 0.2f, 0.5f, 1, 5, 0, 1, 1, 0);
            }

            //int perma = 1;
            int timestep = 24;
            int period = 2014;

            var global = new Global
            {
                DTime = seconds,
                Tmean = tmean.ToArray(),
                PET = pet.ToArray(),
                Precip = precip.ToArray(),
                Discharge = discharge,
                Catchmentinfo = new CatchmentInfo
                {
                    Area_m2 = catchmentInfo.AreaM2.ToArray(),
                    rel_Area = catchmentInfo.RelArea.ToArray(),
                    Elevation_bottom = catchmentInfo.ElevationBottom.ToArray(),
                    Elevation_top = catchmentInfo.ElevationTop.ToArray(),
                    Total_Area_m2 = catchmentInfo.TotalAreaM2[0]
                },
                Tground = tground.ToArray()
            };

            float test_KGE = HBVModel.KGE_HBV(pars_test, perma, seconds, global);
            ElevationBandData test_HBV = HBVModel.HBV_permafrost_semidistributed(pars_test, perma, seconds, global);



            HBVParams[] minBound = new HBVParams[numElevationBands]; // new HBVParams(0.1f, 0.01f, 0.10f, 0.10f, 0.0005f, 0.0005f, 0.01f, 0.0f, -3f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 1f, 1f, 1f, 1f);
            HBVParams[] maxBound = new HBVParams[numElevationBands]; // new HBVParams(1000f, 7f, 1f, 3f, 0.3f, 0.3f, 6f, 4f, 4f, 7f, 20f, 1f, 1f, 0.8f, 5f, 100f, 1f, 1f, 1f, 1f);

            for (int i = 0; i < numElevationBands; i++)
            {
                minBound[i] = new HBVParams(0.1f, 0.01f, 0.10f, 0.10f, 0.0005f, 0.0005f, 0.01f, 0.0f, -3f, 0f, 0f, 0f, 0f, 0f, 0.006f, 0f, 1f, 1f, 1f, 1f);
                maxBound[i] = new HBVParams(1000f, 7f, 1f, 3f, 0.3f, 0.3f, 6f, 4f, 4f, 7f, 20f, 1f, 1f, 0.8f, 5f, 100f, 1f, 1f, 1f, 1f);
            }

            Func<float[], float> objFun = delegate (float[] optimParams)
            {

               

                // Use the conversion function
                HBVParams[] pars = FloatToHBVParams(optimParams, numElevationBands, numParameters);

                return HBVModel.KGE_HBV(pars, perma, seconds, global);
            };

           

            PSO optimizer = new PSO(numSubcatchments * numParameters, 1000, objFun, HBVParamsToFloat(minBound), HBVParamsToFloat(maxBound));


            Task t = new Task(delegate ()
            {
                optimizer.IterateNTimes();
            });


            return (t, optimizer);
        }

        //public static (ElevationBandData, CsvMeteoData) Run(CsvMeteoData meteoData, CsvDataElevationBands catchmentInfo, HBVParams[] loadedParams = null)
        //{
        //    int numElevationBands = catchmentInfo.AreaM2.Count;
        //    HBVParams[] pars_test = new HBVParams[numElevationBands];
        //    for (int i = 0; i < numElevationBands; i++)
        //    {
        //        pars_test[i] = new HBVParams(500, 3, 0.5f, 1, 0.1f, 0.01f, 3, 2, 0, 2, 5, 0.5f, 0.2f, 0.5f, 1, 5, 0, 1, 1, 0);
        //    }

        //    return Run(meteoData, catchmentInfo, pars_test);
        //}

        public static (ElevationBandData, CsvMeteoData) Run(CsvMeteoData meteoData,
            CsvDataElevationBands catchmentInfo, HBVParams[] pars_test=null, ModelType perma = ModelType.NormalHBV)
        {
            int numElevationBands = catchmentInfo.AreaM2.Count;
            if (pars_test == null)
            {
                pars_test = new HBVParams[numElevationBands];
                for (int i = 0; i < numElevationBands; i++)
                {
                    pars_test[i] = new HBVParams(500, 3, 0.5f, 1, 0.1f, 0.01f, 3, 2, 0, 2, 5, 0.5f, 0.2f, 0.5f, 1, 5, 0, 1, 1, 0);
                }
            }

            // Load meteorological data
            // CsvMeteoData meteoData = new CsvMeteoData(@"C:\Daten\meteodata_24h.csv");

            // Load catchment info
            // var catchmentInfo = new CsvDataElevationBands(@"C:\Daten\Elevation_bands.csv");

            // Set timestep in seconds
            float seconds = 86400f; // 24h
            // Other options:
            // float seconds = 43200f; // 12h
            // float seconds = 21600f; // 6h  
            // float seconds = 10800f; // 3h
            // float seconds = 3600f;  // 1h

            // Create arrays from the needed data
            var dTime = meteoData.allDateTimes;
            var tmean = meteoData.Tmean;
            var pet = meteoData.PET;
            var precip = meteoData.Precipitation_mm_even;
            var discharge = meteoData.Discharge.Select(x => x / 1000f).ToArray(); // Convert to m3/s
            var tground = meteoData.Tground;

            int numSubcatchments = 8;
            int numParameters = 20;

            // int numElevationBands = catchmentInfo.ElevationBottom.Count;





            // Run the HBV model for each elevation band and combine the discharge at the end

            // Parameter bounds
            // Min bounds: [0.1, 0.01, 0.10, 0.10, 0.0005, 0.0005, 0.01, 0.0, -3, 0, 0, 0, 0, 0, 0, 0, -5, -2, 0, -5]
            // Max bounds: [1000, 7, 1, 3, 0.3, 0.3, 6, 4, 4, 7, 20, 1, 1, 0.8, 5, 100, 0, 2, 5, 2]

            // int perma = 1;
            int timestep = 24;
            int period = 2014;

            var global = new Global
            {
                DTime = seconds,
                Tmean = tmean.ToArray(),
                PET = pet.ToArray(),
                Precip = precip.ToArray(),
                Discharge = discharge,
                Catchmentinfo = new CatchmentInfo
                {
                    Area_m2 = catchmentInfo.AreaM2.ToArray(),
                    rel_Area = catchmentInfo.RelArea.ToArray(),
                    Elevation_bottom = catchmentInfo.ElevationBottom.ToArray(),
                    Elevation_top = catchmentInfo.ElevationTop.ToArray(),
                    Total_Area_m2 = catchmentInfo.TotalAreaM2[0]
                },
                Tground = tground.ToArray()
            };

            float test_KGE = HBVModel.KGE_HBV(pars_test, perma, seconds, global);
            ElevationBandData test_HBV = HBVModel.HBV_permafrost_semidistributed(pars_test, perma, seconds, global);


            return (test_HBV, meteoData);

            // Run HBV model

            //Func<float[], float> objFun = delegate (float[] optimParams)
            //{
            //    HBVModel.KGE_HBV();
            //};


        }

        public void EvaluateKGE_HBV()
        {

        }

    }
}
