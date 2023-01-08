using System;

namespace JA
{
    public static class DoubleConstants
    {
        public const double ulp = 1.0 / 2251799813685248;
        public const double pi = Math.PI;
        public const double deg = pi / 180;
        public const double mm = 0.001, inch = 0.0254;

        public const double tiny = 64 * ulp;
        public const double small = 1.0/8388608;

        internal static readonly Random rng = new Random();

        #region Functions
        public static double Random(double minValue=0, double maxValue=1) => minValue + (maxValue-minValue) * rng.NextDouble();
        public static double Sqr(double x) => x*x;
        public static double Sqrt(double x) => Math.Sqrt(x);

        public static double Raise(double x, int e)
        {
            if (e == 0) return 1;
            if (e == 1) return x;
            if (e == 2) return x * x;
            if (e < 0) return 1 / Raise(x, -e);
            return Math.Pow(x, e);
        }

        public static double Sign(double magnitude, double value)
            => magnitude * Math.Sign(value);

        #endregion

    }
}
