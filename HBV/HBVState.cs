namespace HBV
{
    public class HBVState
    {
        public float sm1;    // initial storage soil moisture
        public float sw1;    // initial storage surface water 
        public float gw1;    // initial storage ground water
        public float sp1;    // initial storage snow pack
        public float mw1;    // initial storage melt water

        public float[] sp;   // storage snow pack array
        public float[] mw;   // storage melted water array
        public float[] sm;   // storage soil moisture array
        public float[] sw;   // storage surface water array
        public float[] gw;   // storage ground water array
    }
}
