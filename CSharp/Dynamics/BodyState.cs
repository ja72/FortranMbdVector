using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace JA.Dynamics
{
    using JA.LinearAlgebra;
    using static DoubleConstants;
    using static LinearAlgebra.LinearAlgebra;

    public unsafe struct BodyState :
        System.Collections.ICollection,
        ICollection<double>,
        IEquatable<BodyState>
    {
        public const int StateSize = 13;
        public const int PositionOffset = 0;
        public const int OrientationOffset = 3;
        public const int MomnetumOffset = 7;
        public const int AngularOffset = 10;

        fixed double data[StateSize];

        #region Factory
        public BodyState(double[] y)
        {
            for (int i = 0; i < StateSize; i++)
            {
                data[i] = y[i];
            }
        }
        public BodyState(Vector3 position, Quaternion orientation, Vector3 momentum, Vector3 angular)
        {
            data[PositionOffset] = position.X;
            data[PositionOffset + 1] = position.Y;
            data[PositionOffset + 2] = position.Z;
            data[OrientationOffset] = orientation.Vector.X;
            data[OrientationOffset + 1] = orientation.Vector.Y;
            data[OrientationOffset + 2] = orientation.Vector.Z;
            data[OrientationOffset + 3] = orientation.Scalar;
            data[MomnetumOffset] = momentum.X;
            data[MomnetumOffset + 1] = momentum.Y;
            data[MomnetumOffset + 2] = momentum.Z;
            data[AngularOffset] = angular.X;
            data[AngularOffset + 1] = angular.Y;
            data[AngularOffset + 2] = angular.Z;
        }
        #endregion        

        #region Algebra
        public static BodyState Scale(BodyState a, double alpha)
        {
            BodyState c = default;
            for (int i = 0; i < StateSize; i++)
            {
                c.data[i] = alpha * a.data[i];
            }
            return c;
        }
        public static BodyState AddScale(BodyState a, BodyState b, double alpha, double beta)
        {
            BodyState c = default;
            for (int i = 0; i < StateSize; i++)
            {
                c.data[i] = alpha * a.data[i] + beta * b.data[i];
            }
            return c;
        }
        public void AddScale(BodyState b, double factor)
        {
            for (int i = 0; i < StateSize; i++)
            {
                data[i] += factor * b.data[i];
            }
        }
        #endregion

        #region Operators   
        public static BodyState operator +(BodyState a, BodyState b) => AddScale(a,b,1,1);
        public static BodyState operator -(BodyState a) => Scale(a,-1);
        public static BodyState operator -(BodyState a, BodyState b) => AddScale(a,b,1,-1);
        public static BodyState operator *(double a, BodyState b) => Scale(b,a);
        public static BodyState operator *(BodyState a, double b) => Scale(a,b);
        public static BodyState operator /(BodyState a, double b) => Scale(a, 1 / b);

        public static bool operator ==(BodyState left, BodyState right)
        {
            return left.Equals(right);
        }

        public static bool operator !=(BodyState left, BodyState right)
        {
            return !(left == right);
        }
        #endregion

        #region Components

        public void NormalizeOrientation()
        {
            var span = AsSpan();
            var orientation = Quaternion.FromSpan(span, OrientationOffset);
            orientation = orientation.ToUnit();
            data[OrientationOffset] = orientation.Vector.X;
            data[OrientationOffset + 1] = orientation.Vector.Y;
            data[OrientationOffset + 2] = orientation.Vector.Z;
            data[OrientationOffset + 3] = orientation.Scalar;
        }

        public void Deconstruct(out Vector3 position, out Quaternion orientation, out Vector3 momentum, out Vector3 angular)
        {
            var span = AsSpan();
            position = Vector3.FromSpan(span, PositionOffset);
            orientation = Quaternion.FromSpan(span, OrientationOffset);
            momentum = Vector3.FromSpan(span, MomnetumOffset);
            angular = Vector3.FromSpan(span, AngularOffset);
        }
        #endregion

        #region Collections
        public ref double this[int index] => ref data[index];
        public ref double GetPinnableRef() => ref data[0];
        [Browsable(false)]public int Count { get => StateSize; }
        [Browsable(false)]public bool IsReadOnly => true;
        void ICollection<double>.Add(double item) => throw new NotSupportedException();
        void ICollection<double>.Clear() => throw new NotSupportedException();
        bool ICollection<double>.Contains(double item) => throw new NotSupportedException();
        bool ICollection<double>.Remove(double item) => throw new NotSupportedException();
        public IEnumerator<double> GetEnumerator()
        {
            foreach (var x in ToArray())
            {
                yield return x;
            }
        }
        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator() => GetEnumerator();
        public Span<double> AsSpan()
        {
            fixed (double* ptr = data)
            {
                return new Span<double>(ptr, StateSize);
            }
        }
        public double[] ToArray() => AsSpan().ToArray();
        public void CopyTo(Array array, int index)
        {
            Array.Copy(ToArray(), 0, array, index, StateSize);
        }
        public void CopyTo(double[] array, int index)
        {
            for (int i = 0; i < StateSize; i++)
            {
                array[index+i] = data[i];
            }
        }
        [Browsable(false)]public object SyncRoot => null;
        [Browsable(false)]public bool IsSynchronized => false;

        #endregion

        #region Equality
        public override bool Equals(object obj)
        {
            return obj is BodyState state && Equals(state);
        }

        public bool Equals(BodyState other)
        {
            for (int i = 0; i < StateSize; i++)
            {
                if (data[i] != other.data[i])
                {
                    return false;
                }
            }
            return true;
        }

        public override int GetHashCode()
        {
            int hashCode = 1264045569;
            for (int i = 0; i < StateSize; i++)
            {
                hashCode = hashCode * -1521134295 + data[i].GetHashCode();

            }
            return hashCode;
        }


        #endregion

    }
}
