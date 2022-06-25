using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace JA
{
    public static class SingleConstants
    {
        public const float ulp = 1f / 8388608;
        public const float pi = (float)Math.PI;
        public const float deg = pi / 180;
        public const float tiny = 4 * ulp;
        public const float small = 16* ulp;

        #region Functions
        public static float Random(float minValue=0, float maxValue=1) => minValue + (float)( (maxValue-minValue) * DoubleConstants.rng.NextDouble());
        public static float Cos(float θ) => (float)Math.Cos(θ);
        public static float Sin(float θ) => (float)Math.Sin(θ);
        public static float Tan(float θ) => (float)Math.Tan(θ);
        public static float Asin(float t) => (float)Math.Asin(t);
        public static float Acos(float t) => (float)Math.Acos(t);
        public static float Atan(float t) => (float)Math.Atan(t);
        public static float Atan(float dy, float dx) => (float)Math.Atan2(dy, dx);
        public static float Sqr(float x) => x*x;
        public static float Sqrt(float x) => (float)Math.Sqrt(x);
        public static float Pow(float x, float y) => (float)Math.Pow(x, y);
        public static float Raise(float x, int e)
        {
            if (e == 0) return 1;
            if (e == 1) return x;
            if (e == 2) return x * x;
            if (e < 0) return 1 / Raise(x, -e);
            return Pow(x, e);
        }
        #endregion

        #region Roots


        #endregion
    }
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
