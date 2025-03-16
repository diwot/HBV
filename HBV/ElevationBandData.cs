namespace HBV
{
    public class ElevationBandData
    {
        public float[] r;      // rainfall
        public float[] ssp;    // storage snow pack
        public float[] ssw;    // storage surface water
        public float[] ssm;    // storage soil moisture
        public float[] sgw;    // storage ground water
        public float[] qd;     // direct discharge
        public float[] qc;     // capillary transport
        public float[] qin;    // infiltration
        public float[] qf;     // fast run-off
        public float[] qs;     // slow run-off
        public float[] qr;     // river discharge

        public float[] tg;     // corrected ground temperature

        public ElevationBandData() { }

        public ElevationBandData(int length)
        {
            r= new float[length];
            ssp= new float[length + 1];
            ssw= new float[length + 1];
            ssm= new float[length + 1];
            sgw= new float[length + 1];
            qd= new float[length];
            qc= new float[length];
            qin= new float[length];
            qf= new float[length];
            qs= new float[length];
            qr= new float[length];

            tg = new float[length];
        }
    }
}
