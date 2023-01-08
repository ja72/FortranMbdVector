using System;
using System.Linq;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.ComponentModel;

namespace JA.LinearAlgebra
{
    using static Math;
    using static DoubleConstants;
    using static LinearAlgebra;

    [TypeConverter(typeof(ExpandableObjectConverter))]
    public readonly struct Matrix3 :
        ICollection<double>,
        System.Collections.ICollection,
        IEquatable<Matrix3>,
        IFormattable
    {
        readonly (double a11, double a12, double a13,
            double a21, double a22, double a23,
            double a31, double a32, double a33) data;

        #region Factory

        public Matrix3(double a11, double a12, double a13,
            double a21, double a22, double a23,
            double a31, double a32, double a33)
        {
            this.data = (a11, a12, a13, a21, a22, a23, a31, a32, a33);
        }
        public static readonly Matrix3 Zero = new Matrix3(0, 0, 0, 0, 0, 0, 0, 0, 0);
        public static readonly Matrix3 Identity = new Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1);
        public static readonly Matrix3 Ones = new Matrix3(1, 1, 1, 1, 1, 1, 1, 1, 1);
        public static Matrix3 FromRows(Vector3 row1, Vector3 row2, Vector3 row3)
            => new Matrix3(
                row1.X, row1.Y, row1.Z,
                row2.X, row2.Y, row2.Z,
                row3.X, row3.Y, row3.Z);
        public static Matrix3 FromColumns(Vector3 column1, Vector3 column2, Vector3 column3)
            => new Matrix3(
                column1.X, column2.X, column3.X,
                column1.Y, column2.Y, column3.Y,
                column1.Z, column2.Z, column3.Z);

        public static Matrix3 Scalar(double a) => new Matrix3(a, 0, 0, 0, a, 0, 0, 0, a);
        public static implicit operator Matrix3(double a) => Scalar(a);
        public static Matrix3 Diagonal(double a11, double a22, double a33) => new Matrix3(a11, 0, 0, 0, a22, 0, 0, 0, a33);
        public static Matrix3 Symmetric(double a11, double a22, double a33, double a12, double a13, double a23)
            => new Matrix3(a11, a12, a13, a12, a22, a23, a13, a23, a33);
        public static Matrix3 Diagonal(Vector3 vector)
            => Diagonal(vector.X, vector.Y, vector.Z);
        public static Matrix3 RandomFull(double minValue = 0, double maxValue = 1) => new Matrix3(
                Random(minValue, maxValue), Random(minValue, maxValue), Random(minValue, maxValue),
                Random(minValue, maxValue), Random(minValue, maxValue), Random(minValue, maxValue),
                Random(minValue, maxValue), Random(minValue, maxValue), Random(minValue, maxValue));
        public static Matrix3 RandomSymm(double minValue = 0, double maxValue = 1)
            => Symmetric(Random(minValue, maxValue), Random(minValue, maxValue), Random(minValue, maxValue),
                        Random(minValue, maxValue), Random(minValue, maxValue), Random(minValue, maxValue));
        public static Matrix3 RotationAboutX(double angle)
        {
            double c = Math.Cos(angle), s = Math.Sin(angle);
            return new Matrix3(
                1, 0, 0,
                0, c, -s,
                0, s, c);
        }
        public static Matrix3 RotationAboutY(double angle)
        {
            double c = Math.Cos(angle), s = Math.Sin(angle);
            return new Matrix3(
                c, 0, s,
                0, 1, 0,
                -s, 0, c);
        }
        public static Matrix3 RotationAboutZ(double angle)
        {
            double c = Math.Cos(angle), s = Math.Sin(angle);
            return new Matrix3(
                c, -s, 0,
                s, c, 0,
                0, 0, 1);
        }
        public static Matrix3 RotationAbout(Axis axis, double angle)
        {
            switch (axis)
            {
                case Axis.X: return RotationAboutX(angle);
                case Axis.Y: return RotationAboutY(angle);
                case Axis.Z: return RotationAboutZ(angle);
                default:
                    throw new NotSupportedException($"Invalid Axis {axis}");
            }
        }
        public static Matrix3 RotationAbout(Vector3 axis, double angle)
        {
            axis = Vector3.Normalize(axis);
            var kx = axis.CrossOp();
            var neg_kxkx = LinearAlgebra.Mmoi(axis);
            return 1 + Math.Sin(angle) * kx - (1 - Math.Cos(angle)) * neg_kxkx;
        }
        public double MaxValue() => Max(GetRow(0).MaxValue(), Max(GetRow(1).MaxValue(), GetRow(2).MaxValue()));
        public double MinValue() => Min(GetRow(0).MinValue(), Min(GetRow(1).MinValue(), GetRow(2).MinValue()));

        public Matrix3 UnderflowToZero(double scale = 1)
        {
            var m = MaxValue() - MinValue();
            var s = scale * m * Ones;
            return (this + s) - s;
        }
        #endregion

        #region Properties
        public int Rows { get => 3; }
        public int Columns { get => 3; }

        public double A11 { get => data.a11; }
        public double A12 { get => data.a12; }
        public double A13 { get => data.a13; }
        public double A21 { get => data.a21; }
        public double A22 { get => data.a22; }
        public double A23 { get => data.a23; }
        public double A31 { get => data.a31; }
        public double A32 { get => data.a32; }
        public double A33 { get => data.a33; }
        public double Determinant
        {
            get => data.a11 * (data.a22 * data.a33 - data.a23 * data.a32)
                + data.a12 * (data.a23 * data.a31 - data.a21 * data.a33)
                + data.a13 * (data.a21 * data.a32 - data.a22 * data.a31);
        }
        public bool IsSingular { get => Determinant == 0; }
        public bool IsZero { get => this == Zero; }
        public bool IsIdentity { get => this == Identity; }
        public Vector3 GetRow(int row)
        {
            switch (row)
            {
                case 0: return new Vector3(data.a11, data.a12, data.a13);
                case 1: return new Vector3(data.a21, data.a22, data.a23);
                case 2: return new Vector3(data.a31, data.a32, data.a33);
                default:
                    throw new ArgumentOutOfRangeException(nameof(row));
            }
        }
        public Vector3 GetColumn(int column)
        {
            switch (column)
            {
                case 0: return new Vector3(data.a11, data.a21, data.a31);
                case 1: return new Vector3(data.a12, data.a22, data.a32);
                case 2: return new Vector3(data.a13, data.a23, data.a33);
                default:
                    throw new ArgumentOutOfRangeException(nameof(column));
            }
        }

        public Vector3 GetDiagonal() => new Vector3(data.a11, data.a22, data.a33);

        #endregion

        #region Algebra
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Matrix3 Abs(Matrix3 matrix)
            => new Matrix3(
                Math.Abs(matrix.data.a11), Math.Abs(matrix.data.a12), Math.Abs(matrix.data.a13),
                Math.Abs(matrix.data.a21), Math.Abs(matrix.data.a22), Math.Abs(matrix.data.a23),
                Math.Abs(matrix.data.a31), Math.Abs(matrix.data.a32), Math.Abs(matrix.data.a33));

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Matrix3 Add(Matrix3 matrix, Matrix3 other)
            => new Matrix3(
                matrix.data.a11 + other.data.a11, matrix.data.a12 + other.data.a12, matrix.data.a13 + other.data.a13,
                matrix.data.a21 + other.data.a21, matrix.data.a22 + other.data.a22, matrix.data.a23 + other.data.a23,
                matrix.data.a31 + other.data.a31, matrix.data.a32 + other.data.a32, matrix.data.a33 + other.data.a33);
        [MethodImpl(MethodImplOptions.AggressiveInlining)]        
        public static Matrix3 Subtract(Matrix3 matrix, Matrix3 other)
            => new Matrix3(
                matrix.data.a11 - other.data.a11, matrix.data.a12 - other.data.a12, matrix.data.a13 - other.data.a13,
                matrix.data.a21 - other.data.a21, matrix.data.a22 - other.data.a22, matrix.data.a23 - other.data.a23,
                matrix.data.a31 - other.data.a31, matrix.data.a32 - other.data.a32, matrix.data.a33 - other.data.a33);
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Matrix3 Scale(double factor, Matrix3 matrix)
            => new Matrix3(
            factor * matrix.data.a11, factor * matrix.data.a12, factor * matrix.data.a13,
            factor * matrix.data.a21, factor * matrix.data.a22, factor * matrix.data.a23,
            factor * matrix.data.a31, factor * matrix.data.a32, factor * matrix.data.a33);
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Matrix3 Negate(Matrix3 matrix)
            => new Matrix3(
            -matrix.data.a11, -matrix.data.a12, -matrix.data.a13,
            -matrix.data.a21, -matrix.data.a22, -matrix.data.a23,
            -matrix.data.a31, -matrix.data.a32, -matrix.data.a33);
        public static Vector3 Product(Matrix3 matrix, Vector3 vector)
            => new Vector3(
                Vector3.Dot(matrix.GetRow(0), vector),
                Vector3.Dot(matrix.GetRow(1), vector),
                Vector3.Dot(matrix.GetRow(2), vector));
        public static Vector3 Product(Vector3 vector, Matrix3 matrix)
            => new Vector3(
                Vector3.Dot(vector, matrix.GetColumn(0)),
                Vector3.Dot(vector, matrix.GetColumn(1)),
                Vector3.Dot(vector, matrix.GetColumn(2)));
        public static Matrix3 Product(Matrix3 matrix, Matrix3 other)
        {
            Vector3 row1 = matrix.GetRow(0);
            Vector3 row2 = matrix.GetRow(1);
            Vector3 row3 = matrix.GetRow(2);
            Vector3 column1 = other.GetColumn(0);
            Vector3 column2 = other.GetColumn(1);
            Vector3 column3 = other.GetColumn(2);

            return new Matrix3(
                           Vector3.Dot(row1, column1), Vector3.Dot(row1, column2), Vector3.Dot(row1, column3),
                           Vector3.Dot(row2, column1), Vector3.Dot(row2, column2), Vector3.Dot(row2, column3),
                           Vector3.Dot(row3, column1), Vector3.Dot(row3, column2), Vector3.Dot(row3, column3));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Matrix3 Transpose() 
            => new Matrix3(
                data.a11, data.a21, data.a31,
                data.a12, data.a22, data.a32,
                data.a13, data.a23, data.a33);

        public Matrix3 Inverse()
        {
            var id = 1 / Determinant;
            return new Matrix3(
                id*(data.a22 * data.a33 - data.a23 * data.a32), 
                id*(data.a13 * data.a32 - data.a12 * data.a33), 
                id*(data.a12 * data.a23 - data.a13 * data.a22),
                id*(data.a23 * data.a31 - data.a21 * data.a33), 
                id*(data.a11 * data.a33 - data.a13 * data.a31), 
                id*(data.a13 * data.a21 - data.a11 * data.a23),
                id*(data.a21 * data.a32 - data.a22 * data.a31), 
                id*(data.a12 * data.a31 - data.a11 * data.a32), 
                id*(data.a11 * data.a22 - data.a12 * data.a21));
        }

        public Vector3 Solve(Vector3 vector) => Inverse() * vector;
        public Matrix3 Solve(Matrix3 other) => Inverse() * other;
        #endregion

        #region Operators
        public static Matrix3 operator +(Matrix3 a, Matrix3 b) => Add(a, b);
        public static Matrix3 operator -(Matrix3 a, Matrix3 b) => Subtract(a, b);
        public static Matrix3 operator -(Matrix3 b) => Negate(b);
        public static Matrix3 operator *(double a, Matrix3 b) => Scale(a, b);
        public static Matrix3 operator *(Matrix3 a, double b) => Scale(b, a);
        public static Vector3 operator *(Matrix3 a, Vector3 b) => Product(a, b);
        public static Vector3 operator *(Vector3 a, Matrix3 b) => Product(a, b);
        public static Matrix3 operator *(Matrix3 a, Matrix3 b) => Product(a, b);
        public static Matrix3 operator /(Matrix3 a, double b) => Scale(1 / b, a);
        public static Matrix3 operator ~(Matrix3 b) => b.Transpose();
        public static Matrix3 operator !(Matrix3 b) => b.Inverse();
        public static Vector3 operator /(Vector3 a, Matrix3 b) => b.Solve(a);
        public static Matrix3 operator /(Matrix3 a, Matrix3 b) => b.Solve(a);
        #endregion

        #region Eigenvalues
        /// <summary>
        /// Calculates the three eigenvalues analytically.
        /// </summary>
        /// <remarks>
        /// Code taken from:
        /// https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/index.html
        /// </remarks>
        /// <returns>A vector containing the three eigenvalues.</returns>
        public Vector3 GetEigenValues()
        {
            //      Determine coefficients of characteristic polynomial. We write
            //      | A   D   F  |
            // A =  | D*  B   E  |
            //      | F*  E*  C  |

            var de = data.a12 * data.a23;
            var dd = data.a12 * data.a12;
            var ee = data.a23 * data.a23;
            var ff = data.a13 * data.a13;
            var m = data.a11 + data.a22 + data.a33;
            var c1 = (data.a11 * data.a22 + data.a11 * data.a33 + data.a22 * data.a33) - (dd + ee + ff);
            var c0 = data.a33 * dd + data.a11 * ee + data.a22 * ff - data.a11 * data.a22 * data.a33 - 2.0 * data.a13 * de;

            var p = m * m - 3.0 * c1;
            var q = m * (p - (3.0 / 2.0) * c1) - (27.0 / 2.0) * c0;
            var sqrt_p =  Math.Sqrt(Math.Abs(p));

            var sqrt_z =  Math.Sqrt(Math.Abs(27.0 * (0.25 * c1 * c1 * (p - c1) + c0 * (q + 27.0 / 4.0 * c0))));
            var phi = (1 / 3.0) *  Math.Atan2(sqrt_z, q);

            var c = sqrt_p *  Math.Cos(phi);
            var s = sqrt_p *  Math.Abs(Math.Sin(phi))/ Math.Sqrt(3);

            var w = (1 / 3.0) * (m - c);

            // sort the eigenvalues
            if (c >= s)
            {
                return new Vector3(
                    w - s,
                    w + s,
                    w + c);
            }
            else if (c >= -s)
            {
                return new Vector3(
                    w - s,
                    w + c,
                    w + s);
            }
            else
            {
                return new Vector3(
                    w + c,
                    w - s,
                    w + s);
            }
        }

        public Matrix3 GetDiagonalEigenvalues()
            => Diagonal(GetEigenValues());

        public Matrix3 GetEigenVectors() => GetEigenVectors(GetEigenValues());
        public Matrix3 GetEigenVectors(Vector3 eigenValues)
        {
            Vector3 ev1 = GetEigenVector(eigenValues.X).ToUnit();
            Vector3 ev2 = GetEigenVector(eigenValues.Y).ToUnit();
            Vector3 ev3 = GetEigenVector(eigenValues.Z).ToUnit();

            return FromColumns(ev1, ev2, ev3);
        }
        Vector3 GetEigenVector(double w)
        {
            return new Vector3(
                  data.a12 * (data.a23 - data.a33 + w) - data.a13 * (data.a22 - data.a23 - w)
                + data.a22 * (data.a33 - w) - data.a23 * data.a23 - w * (data.a33 - w),
                -data.a11 * (data.a23 - data.a33 + w) + data.a12 * (data.a13 - data.a33 + w)
                - data.a13 * data.a13 + data.a13 * data.a23 + w * (data.a23 - data.a33 + w),
                  data.a11 * (data.a22 - data.a23 - w) - data.a12 * data.a12 + data.a12 * (data.a13 + data.a23)
                + data.a13 * (w - data.a22) - w * (data.a22 - data.a23 - w));
        }

        #endregion


        #region Formatting
        public override string ToString() => ToString("g");
        public string ToString(string formatting) => ToString(formatting, null);
        public string ToString(string format, IFormatProvider provider)
        {
            return $"[{data.a11.ToString(format, provider)},{data.a12.ToString(format, provider)},{data.a13.ToString(format, provider)}|" + Environment.NewLine +
                $"{data.a21.ToString(format, provider)},{data.a22.ToString(format, provider)},{data.a23.ToString(format, provider)}|" + Environment.NewLine +
                $"{data.a31.ToString(format, provider)},{data.a32.ToString(format, provider)},{data.a33.ToString(format, provider)}]";
        }
        #endregion

        #region Equality
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public override bool Equals(object obj)
        {
            return obj is Matrix3 matrix &&
                   Equals(matrix);
        }
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public bool Equals(Matrix3 other)
        {
            other = Abs(this - other);
            return other.data.a11 < tiny
                && other.data.a12 < tiny
                && other.data.a13 < tiny
                && other.data.a21 < tiny
                && other.data.a22 < tiny
                && other.data.a23 < tiny
                && other.data.a31 < tiny
                && other.data.a32 < tiny
                && other.data.a33 < tiny;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public override int GetHashCode()
        {
            return 1768953197 + data.GetHashCode();
        }

        public static bool operator ==(Matrix3 left, Matrix3 right)
        {
            return left.Equals(right);
        }

        public static bool operator !=(Matrix3 left, Matrix3 right)
        {
            return !(left == right);
        }
        #endregion

        #region Collection
        public int Count => 9;
        public bool IsReadOnly => true;
        void ICollection<double>.Add(double item) => throw new NotSupportedException();
        void ICollection<double>.Clear() => throw new NotSupportedException();
        bool ICollection<double>.Contains(double item) => throw new NotSupportedException();
        bool ICollection<double>.Remove(double item) => throw new NotSupportedException();
        public IEnumerator<double> GetEnumerator()
        {
            yield return data.a11;
            yield return data.a12;
            yield return data.a13;
            yield return data.a21;
            yield return data.a22;
            yield return data.a23;
            yield return data.a31;
            yield return data.a32;
            yield return data.a33;
        }
        public double[] ToArray() => new[] { 
            data.a11,data.a12,data.a13,
            data.a21,data.a22,data.a23,
            data.a31,data.a32,data.a33 };
        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator() => GetEnumerator();
        public void CopyTo(double[] array, int index) => CopyTo(array as Array, index);
        public void CopyTo(Array array, int index) => Array.Copy(ToArray(), 0, array, index, Count);
        public object SyncRoot => null;
        public bool IsSynchronized => false;
        #endregion

    }
}

