using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using JA.Dynamics;

namespace JA.Geomtery
{
    using JA.LinearAlgebra;

    using static Math;
    using static DoubleConstants;
    using static LinearAlgebra.LinearAlgebra;

    public enum KnownShape
    {
        Sphere,
        Cylinder
    }
    public static class Shapes
    {
        public static Shape New(KnownShape known, params double[] dims)
        {
            switch (known)
            {
                case KnownShape.Sphere:
                    return Sphere(dims[0]);
                case KnownShape.Cylinder:
                    return CylinderAlongZ(dims[0], dims[1]);
                default:
                    throw new NotSupportedException($"{known} is not supported.");
            }
        }
        public static SphereShape Sphere(double radius)
            => new SphereShape(radius);
        public static CylinderShape CylinderAlongZ(double length, double radius)
            => new CylinderShape(length, radius);
    }
    public abstract class Shape
    {
        protected Shape(KnownShape known, double volume, Vector3 center, (double v_1, double v_2, double v_3) vmmoi)
        {
            Known = known;
            Volume = volume;
            Center = center;
            Vmmoi = vmmoi;
        }
        public KnownShape Known { get;  }
        public double Volume { get; }
        public Vector3 Center { get; }
        public (double v_1, double v_2, double v_3) Vmmoi { get; }
        public abstract Vector3 GetNearestPoint(Vector3 direction, Matrix3 rotation);
    }

    public class SphereShape : Shape
    {
        static double CalcVol(double radius) => 4 * pi * radius * radius * radius / 3;
        static (double v_1, double v_2, double v_3) CalcVmmoi(double radius)
        {
            double v = 2 * radius * radius / 5;
            return (v, v, v);
        }
        public SphereShape(double radius) 
            : base(KnownShape.Sphere, 
                  CalcVol(radius), 
                  Vector3.Zero, 
                  CalcVmmoi(radius))
        {
            Radius = radius;
        }

        public double Radius { get; }
        public double Diameter { get => 2*Radius; }

        public override Vector3 GetNearestPoint(Vector3 direction, Matrix3 rotation)
        {
            //tex: Point $\vec{\rm pos} = -r\, \vec{e}$
            return -direction * Radius;
        }
    }
    public class CylinderShape : Shape
    {
        public CylinderShape(double length, double radius)
            : base(KnownShape.Cylinder, CalcVol(length,radius), Vector3.Zero, CalcVmmoi(length,radius))
        {
            Length = length;
            Radius = radius;
        }

        static double CalcVol(double l, double r) => l*pi*r*r;
        static (double v_1, double v_2, double v_3) CalcVmmoi(double l, double r)
        {
            double v1 = l*l/12 + r*r/4;
            double v2 = r * r / 2;
            return (v1, v1, v2);
        }

        public double Length { get; }
        public double Radius { get; }
        public double Diameter { get => 2*Radius; }

        public override Vector3 GetNearestPoint(Vector3 direction, Matrix3 rotation)
        {
            //tex: Point $\vec{\rm pos} = \mathbf{R}\, \pmatrix{
            //  \frac{d}{2} \cos \varphi \\ 
            //  \frac{d}{2} \sin \varphi \\ 
            //  \pm \frac{\ell}{2} }$
            double d = Diameter;
            double l = Length;
            double alp = -Sign(1.0, Dot(direction, rotation.GetColumn(2)));
            double x = Dot(direction, rotation.GetColumn(0));
            double y = Dot(direction, rotation.GetColumn(1));
            if (Abs(x) <= tiny && Abs(y) <= tiny)
            {
                Vector3 point = new Vector3(0, 0, (l / 2) * alp);
                return rotation * point;
            }
            else
            {
                double φ = Atan2(y, x) - pi;
                Vector3 point = new Vector3((d / 2) * Cos(φ), (d / 2) * Sin(φ), (l / 2) * alp);
                return rotation * point;
            }
        }
    }

}
