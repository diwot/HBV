namespace HBV
{
    public class HBVParams
    {
        public float FC;      // max soil moisture storage, field capacity
        public float BETA;    // shape coefficient governing fate of water input to soil moisture storage
        public float LP;      // threshold for reduction of evaporation
        public float ALFA;    // non-linearity parameter
        public float KF;      // fast run-off parameter
        public float KS;      // slow run-off parameter
        public float PERC;    // percolation, max flow from upper to lower storage
        public float CFLUX;   // rate of capillary rise
        public float TT;       // threshold temperature for snowfall
        public float TTI;      // threshold temperature for snow melt
        public float CFMAX;    // degree-day factor
        public float FOCFMAX;  // degree day factor, rate of snowmelt
        public float CFR;      // refreezing factor
        public float WHC;      // water holding capacity


        public float MAXBAS;
        public float precip_corr;

        public float TTG;      // threshold temperature for ice in ground
        public float TTGI;     // threshold temperature variation  
        public float corrTg;   // correction for ground temperature
        public float infil_stop; // threshold temperature for infiltration stop

        public HBVParams(float fC, float bETA, float lP, float aLFA, float kF, float kS, float pERC, float cFLUX, float tT, float tTI, float cFMAX, float fOCFMAX, float cFR, float wHC, float mAXBAS, float precip_corr, float tTG, float tTGI, float corrTg, float infil_stop)
        {
            FC = fC;
            BETA = bETA;
            LP = lP;
            ALFA = aLFA;
            KF = kF;
            KS = kS;
            PERC = pERC;
            CFLUX = cFLUX;
            TT = tT;
            TTI = tTI;
            CFMAX = cFMAX;
            FOCFMAX = fOCFMAX;
            CFR = cFR;
            WHC = wHC;
            MAXBAS = mAXBAS;
            this.precip_corr = precip_corr;
            TTG = tTG;
            TTGI = tTGI;
            this.corrTg = corrTg;
            this.infil_stop = infil_stop;
        }
    }
}
