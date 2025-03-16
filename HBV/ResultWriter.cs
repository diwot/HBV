namespace HBV
{
    public static class ResultWriter
    {

        public static List<HBVParams> Load(string path)
        {
            List<HBVParams> parametersList = new List<HBVParams>();

            using (StreamReader reader = new StreamReader(path))
            {
                // Read and discard header line
                reader.ReadLine();

                string line;
                while ((line = reader.ReadLine()) != null)
                {
                    string[] values = line.Split(',');

                    if (values.Length == 20) // Ensure correct number of parameters
                    {
                        HBVParams param = new HBVParams
                        {
                            FC = float.Parse(values[0]),
                            BETA = float.Parse(values[1]),
                            LP = float.Parse(values[2]),
                            ALFA = float.Parse(values[3]),
                            KF = float.Parse(values[4]),
                            KS = float.Parse(values[5]),
                            PERC = float.Parse(values[6]),
                            CFLUX = float.Parse(values[7]),
                            TT = float.Parse(values[8]),
                            TTI = float.Parse(values[9]),
                            CFMAX = float.Parse(values[10]),
                            FOCFMAX = float.Parse(values[11]),
                            CFR = float.Parse(values[12]),
                            WHC = float.Parse(values[13]),
                            MAXBAS = float.Parse(values[14]),
                            precip_corr = float.Parse(values[15]),
                            TTG = float.Parse(values[16]),
                            TTGI = float.Parse(values[17]),
                            corrTg = float.Parse(values[18]),
                            infil_stop = float.Parse(values[19])
                        };

                        parametersList.Add(param);
                    }
                }
            }

            return parametersList;
        }





        public static void Save(HBVParams[] finalParams, string path)
        {
            using (StreamWriter writer = new StreamWriter(path))
            {
                // Write header
                writer.WriteLine("FC,BETA,LP,ALFA,KF,KS,PERC,CFLUX,TT,TTI,CFMAX,FOCFMAX,CFR,WHC,MAXBAS,precip_corr,TTG,TTGI,corrTg,infil_stop");

                // Write each parameter set
                foreach (var param in finalParams)
                {
                    writer.WriteLine($"{param.FC},{param.BETA},{param.LP},{param.ALFA},{param.KF},{param.KS},{param.PERC},{param.CFLUX},{param.TT},{param.TTI},{param.CFMAX},{param.FOCFMAX},{param.CFR},{param.WHC},{param.MAXBAS},{param.precip_corr},{param.TTG},{param.TTGI},{param.corrTg},{param.infil_stop}");
                }
            }
        }
    }
}
