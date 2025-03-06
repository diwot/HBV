namespace HBV
{
    public class HBVOptions
    {
        public float TCON;     // time constant
        public float CEVPFO;   // correction factor for evaporation in forest areas
        public float TT;       // threshold temperature for snowfall
        public float TTI;      // threshold temperature for snow melt
        public float DTTM;     // degree day factor
        public float CFMAX;    // degree-day factor
        public float FOCFMAX;  // degree day factor, rate of snowmelt
        public float CFR;      // refreezing factor
        public float WHC;      // water holding capacity
    }
}
