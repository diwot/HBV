using System.Runtime.CompilerServices;

namespace HBV
{
    //public enum ModelType
    //{
    //    HBV = 0,
    //    HBVPermafrost = 1
    //}



    public struct RandomSlim
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static uint WangHash(uint s) { s = (s ^ 61) ^ (s >> 16); s *= 9; s = s ^ (s >> 4); s *= 0x27d4eb2d; s = s ^ (s >> 15); return s; }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static uint RandomInt(ref uint s) { s ^= s << 13; s ^= s >> 17; s ^= s << 5; return s; }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static float RandomFloat(ref uint s) { return RandomInt(ref s) * 2.3283064365387e-10f; }

    };
}
