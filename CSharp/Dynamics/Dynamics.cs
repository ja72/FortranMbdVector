using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace JA.Dynamics
{
    
    public static class Dynamics
    {
        public static readonly Vector3 o_ = Vector3.Zero;
        public static readonly Vector3 i_ = Vector3.UnitX;
        public static readonly Vector3 j_ = Vector3.UnitY;
        public static readonly Vector3 k_ = Vector3.UnitZ;
        public static readonly Quaternion q_eye = Quaternion.Identity;
        public static readonly Matrix3 zeros_ = Matrix3.Zero;
        public static readonly Matrix3 eye_ = Matrix3.Identity;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Dot(Vector3 a, Vector3 b)
            => Vector3.Dot(a, b);
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3 Cross(Vector3 a, Vector3 b)
            => Vector3.Cross(a, b);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Matrix3 CrossOp(this Vector3 a)
          => new Matrix3(
              0, -a.Z, a.Y,
              a.Z, 0, -a.X,
              -a.Y, a.X, 0);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Matrix3 Mmoi(this Vector3 a, double scale = 1)
        {
            //tex: ${\rm mmoi}(v) = -v\times v\times$
            double xx = scale*a.X * a.X, yy = scale*a.Y * a.Y, zz = scale*a.Z * a.Z;
            double xy = scale*a.X * a.Y, yz = scale*a.Y * a.Z, zx = scale*a.Z * a.X;
            return new Matrix3(
                yy + zz, -xy, -zx,
                -xy, xx + zz, -yz,
                -zx, -yz, xx + yy);
        }

        public static bool IsFinite(this double value)
            => !double.IsInfinity(value) && !double.IsNaN(value);

        public static bool IsFinite(this Vector3 vector)
        {
            return !double.IsInfinity(vector.X) && !double.IsNaN(vector.X)
                && !double.IsInfinity(vector.Y) && !double.IsNaN(vector.Y)
                && !double.IsInfinity(vector.Z) && !double.IsNaN(vector.Z);
        }

    }
}
