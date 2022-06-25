using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Runtime.CompilerServices;

namespace JA.Dynamics
{
    using static DoubleConstants;
    [TypeConverter(typeof(ExpandableObjectConverter))]
    public readonly struct Quaternion :
        ICollection<double>,
        System.Collections.ICollection,
        IEquatable<Quaternion>,
        IFormattable
    {
        readonly (Vector3 v, double s) data;

        #region Factory
        public Quaternion(Vector3 vector, double scalar)
        {
            this.data = (vector, scalar);
        }
        public static readonly Quaternion Zero = new Quaternion(Vector3.Zero, 0);
        public static readonly Quaternion Identity = new Quaternion(Vector3.Zero, 1);

        public static explicit operator Vector3(Quaternion quaternion)
            => quaternion.Vector;
        public static explicit operator double(Quaternion quaternion)
            => quaternion.Scalar;
        public static Quaternion FromVector(Vector3 vector) => new Quaternion(vector, 0);
        public static Quaternion FromAxisAngle(Axis axis, double angle)
        {
            switch (axis)
            {
                case Axis.X: return new Quaternion(Vector3.UnitX * Math.Sin(angle / 2), Math.Cos(angle / 2));
                case Axis.Y: return new Quaternion(Vector3.UnitY * Math.Sin(angle / 2), Math.Cos(angle / 2));
                case Axis.Z: return new Quaternion(Vector3.UnitZ * Math.Sin(angle / 2), Math.Cos(angle / 2));
                default:
                    throw new NotSupportedException($"Invalid Axis {axis}");
            }
        }
        public static Quaternion FromAxisAngle(Vector3 axis, double angle)
            => new Quaternion(Vector3.Normalize(axis) * Math.Sin(angle / 2), Math.Cos(angle / 2));
        public static Quaternion FromRotation(Matrix3 rotation)
        {
            double x = rotation.A32 - rotation.A23;
            double y = rotation.A13 - rotation.A31;
            double z = rotation.A21 - rotation.A12;
            double t = rotation.A11 + rotation.A22 + rotation.A33;
            double s = 0.5 * Math.Sqrt((x * x + y * y + z * z) / (3 - t));
            double f = 1 / (4 * s);
            Vector3 v = new Vector3(f * x, f * y, f * z);
            return new Quaternion(v, s);
        }
        public static Quaternion RotationBetweenVectors(Vector3 from, Vector3 to)
        {
            Vector3 n = Vector3.Cross(from, to);
            double nm = n.Magnitude;
            if (nm < tiny)
            {
                return Identity;
            }
            double m = Math.Sqrt(from.MagnitudeSquared * to.MagnitudeSquared);
            double a = Vector3.Dot(from, to) / m;
            double sin = Math.Sqrt((1 - a) / 2);
            double cos = Math.Sqrt((1 + a) / 2);
            return new Quaternion(n / nm * sin, cos);
        }

        #endregion

        #region Properties
        [Browsable(false)]
        public int Size { get => 4; }
        public Vector3 Vector => data.v;
        public double Scalar => data.s;
        [Browsable(false)] public double Magnitude => Math.Sqrt(MagnitudeSquared);
        [Browsable(false)] public double MagnitudeSquared => data.v.MagnitudeSquared + data.s * data.s;
        public Quaternion ToUnit() => Normalize(this);

        public void Deconstruct(out Vector3 vector, out double scalar)
        {
            vector = data.v;
            scalar = data.s;
        }

        public void Deconstruct(out double x, out double y, out double z, out double w)
        {
            x = data.v.X;
            y = data.v.Y;
            z = data.v.Z;
            w = data.s;
        }

        public Matrix3 ToRotation(bool inverse = false)
        {
            int sign = inverse ? -1 : 1;
            double x = sign * data.v.X, y = sign * data.v.Y, z = sign * data.v.Z, w = data.s;
            double xx = x * x, yy = y * y, zz = z * z;

            return new Matrix3(
                1 - 2 * (yy + zz), -2 * (w * z - x * y), -2 * (-x * z - w * y),
                -2 * (-w * z - x * y), 1 - 2 * (xx + zz), -2 * (w * x - y * z),
                -2 * (-x * z + w * y), -2 * (-w * x - y * z), 1 - 2 * (xx + yy));
        }
        public Vector3 Rotate(Vector3 vector, bool inverse = false)
        {
            int inv = inverse ? -1 : 1;
            Vector3 vxp = Vector3.Cross(data.v, vector);
            Vector3 vxvxp = Vector3.Cross(data.v, vxp);
            return vector + 2 * (inv * data.s * vxp + vxvxp);
        }

        public (Vector3 axis, double angle) GetAxisAngle()
        {
            double vm = data.v.Magnitude;
            if (vm < tiny)
            {
                return (Vector3.Zero, 0);
            }
            double sin = 2 * data.s * vm;
            double cos = data.s * data.s - vm * vm;
            double angle = Math.Atan(sin / cos);
            Vector3 axis = data.v / vm;
            return (axis, angle);
        }
        #endregion

        #region Algebra

        public static Quaternion Normalize(Quaternion quaternion)
        {
            double m2 = quaternion.MagnitudeSquared;
            if (m2 > 0 && m2 != 1)
            {
                return Multiply(quaternion, 1 / Math.Sqrt(m2));
            }
            return quaternion;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Quaternion Add(Quaternion a, Quaternion b)
            => new Quaternion(a.data.v + b.data.v, a.data.s + b.data.s);
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Quaternion Subtract(Quaternion a, Quaternion b)
            => new Quaternion(a.data.v - b.data.v, a.data.s - b.data.s);
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Quaternion Multiply(Quaternion a, double f)
            => new Quaternion(f * a.data.v, f * a.data.s);
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Quaternion Negate(Quaternion a)
            => new Quaternion(-a.data.v, -a.data.s);
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Dot(Quaternion quaternion, Quaternion other)
            => Vector3.Dot(quaternion.data.v, other.data.v) + quaternion.data.s * other.data.s;
        public static Quaternion Multiply(Quaternion a, Quaternion b)
        {
            return new Quaternion(
                a.data.s * b.data.v + b.data.s * a.data.v + Vector3.Cross(a.data.v, b.data.v),
                a.data.s * b.data.s - Vector3.Dot(a.data.v, b.data.v));
        }
        public static Quaternion Divide(Quaternion a, Quaternion b)
        {
            return a * Inverse(b);
        }
        public static Quaternion Conjugate(Quaternion a)
            => new Quaternion(-a.data.v, a.data.s);
        public static Quaternion Inverse(Quaternion a)
        {
            return Conjugate(a) / a.MagnitudeSquared;
        }
        #endregion

        #region Operators
        public static Quaternion operator +(Quaternion a, Quaternion b) => Add(a, b);
        public static Quaternion operator -(Quaternion a) => Negate(a);
        public static Quaternion operator -(Quaternion a, Quaternion b) => Subtract(a, b);
        public static Quaternion operator *(double a, Quaternion b) => Multiply(b, a);
        public static Quaternion operator *(Quaternion a, double b) => Multiply(a, b);
        public static Quaternion operator /(Quaternion a, double b) => Multiply(a, 1 / b);
        public static Quaternion operator *(Quaternion a, Quaternion b) => Multiply(a, b);
        public static Quaternion operator /(Quaternion a, Quaternion b) => Divide(a, b);

        #endregion

        #region Formatting
        public override string ToString() => ToString("g");
        public string ToString(string formatting) => ToString(formatting, null);
        public string ToString(string format, IFormatProvider provider)
        {
            return $"<{Vector.ToString(format, provider)}|{Scalar.ToString(format, provider)}>";
        }
        #endregion

        #region Equality
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public override bool Equals(object obj)
        {
            return obj is Quaternion quaternion
                && Equals(quaternion);
        }
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public bool Equals(Quaternion other)
        {
            double δs = Math.Abs(data.s - other.data.s);
            return δs < tiny && data.v.Equals(other.data.v);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public override int GetHashCode()
        {
            return 1768953197 + data.GetHashCode();
        }

        public static bool operator ==(Quaternion left, Quaternion right)
        {
            return left.Equals(right);
        }

        public static bool operator !=(Quaternion left, Quaternion right)
        {
            return !(left == right);
        }

        #endregion

        #region Collection
        [Browsable(false)] public int Count => 4;
        [Browsable(false)] public bool IsReadOnly => true;
        void ICollection<double>.Add(double item) => throw new NotSupportedException();
        void ICollection<double>.Clear() => throw new NotSupportedException();
        bool ICollection<double>.Contains(double item) => throw new NotSupportedException();
        bool ICollection<double>.Remove(double item) => throw new NotSupportedException();
        public IEnumerator<double> GetEnumerator()
        {
            yield return data.v.X;
            yield return data.v.Y;
            yield return data.v.Z;
            yield return data.s;
        }
        public double[] ToArray() => new[] { data.v.X, data.v.Y, data.v.Z, data.s };
        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator() => GetEnumerator();
        public void CopyTo(double[] array, int index) => CopyTo(array as Array, index);
        public void CopyTo(Array array, int index) => Array.Copy(ToArray(), 0, array, index, Count);
        [Browsable(false)] public object SyncRoot => null;
        [Browsable(false)] public bool IsSynchronized => false;

        public static Quaternion FromSpan(Span<double> array, int index)
        {
            if (index + 4 <= array.Length)
            {
                return new Quaternion(
                    Vector3.FromSpan(array, index),
                    array[index + 3]);
            }
            throw new ArgumentOutOfRangeException(nameof(index));
        }
        public static Quaternion FromArray(double[] array, int index)
        {
            return FromSpan(array, index);
        }

        #endregion

    }
}

