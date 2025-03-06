namespace HBV
{
    public class FloatVec
    {
        private readonly float[] data;

        public FloatVec(int size)
        {
            data = new float[size];
        }

        public FloatVec(float[] array)
        {
            data = array;
        }

        public static implicit operator float[](FloatVec vec)
        {
            return vec.data;
        }

        public float this[int index]
        {
            get { return data[index]; }
            set { data[index] = value; }
        }

        public int Length => data.Length;

        public static FloatVec operator +(FloatVec a, FloatVec b)
        {
            if (a.Length != b.Length)
                throw new ArgumentException("Vectors must be same length");

            var result = new FloatVec(a.Length);
            for (int i = 0; i < a.Length; i++)
            {
                result[i] = a[i] + b[i];
            }
            return result;
        }

        public static FloatVec operator -(FloatVec a, FloatVec b)
        {
            if (a.Length != b.Length)
                throw new ArgumentException("Vectors must be same length");

            var result = new FloatVec(a.Length);
            for (int i = 0; i < a.Length; i++)
            {
                result[i] = a[i] - b[i];
            }
            return result;
        }

        public static FloatVec operator *(FloatVec a, float scalar)
        {
            var result = new FloatVec(a.Length);
            for (int i = 0; i < a.Length; i++)
            {
                result[i] = a[i] * scalar;
            }
            return result;
        }

        public static FloatVec operator *(float scalar, FloatVec a)
        {
            return a * scalar;
        }

        public static FloatVec operator /(FloatVec a, float scalar)
        {
            if (scalar == 0)
                throw new DivideByZeroException();

            var result = new FloatVec(a.Length);
            for (int i = 0; i < a.Length; i++)
            {
                result[i] = a[i] / scalar;
            }
            return result;
        }
    }
}
