namespace HBV
{
    public static class ResultWriter
    {
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
