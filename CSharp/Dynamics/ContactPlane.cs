using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace JA.Dynamics
{
    using JA.LinearAlgebra;

    using static Math;
    using static DoubleConstants;
    using static LinearAlgebra.LinearAlgebra;

    public class ContactPlane
    {
        private Vector3 normal;

        public ContactPlane(Vector3 position, Vector3 normal, double epsilon = 0, double mu = 0)
        {
            Position = position;
            Normal = normal.ToUnit();
            Epsilon = epsilon;
            Mu = mu;
            Reset();
        }

        public Vector3 Position { get; set; }
        public Vector3 Normal { get => normal; set => normal = value.ToUnit(); }
        public bool IsOk { get => Normal.Magnitude > 0; }
        public double Epsilon { get; set; }
        public double Mu { get; set; }
        public Vector3 Slip { get; private set; }
        public double Delta { get; private set; }
        public double Vimp { get; private set; }
        public double Vslip { get; private set; }
        public double Jn { get; private set; }
        public double Js { get; private set; }
        public bool Active { get; private set; }

        public void Reset()
        {
            Slip = Vector3.Zero;
            Delta = 0;
            Vimp = 0;
            Vslip = 0;
            Jn = 0;
            Js = 0;
            Active = false;
        }

        public BodyState Calculate(RigidBody rb, BodyState Y)
        {
            (Vector3 pos_b, Quaternion ori, Vector3 p, Vector3 L_b) = Y;
            Matrix3 R = ori.ToRotation();
            Matrix3 Rt = ori.ToRotation(inverse: true);
            Vector3 c = R * rb.Cg;
            var (I_1, I_2, I_3) = rb.MMoi;
            Matrix3 I_inv = R * Matrix3.Diagonal(1 / I_1, 1 / I_2, 1 / I_3) * Rt;
            Vector3 omega = I_inv * (L_b - (c ^ p));
            Vector3 vee_b = p / rb.Mass - (omega ^ c);

            // Contact point on surface of shape
            Vector3 r_A = pos_b + c + R * rb.Shape.GetNearestPoint(Rt * Normal);
            Delta = Dot(Normal, r_A - Position);
            // Clip contact point on plane
            r_A -= Min(Delta, 0) * Normal;

            // Get contact speed
            Vector3 v_A = vee_b + (omega ^ (r_A - pos_b));
            Vimp = Dot(Normal, v_A);

            Active = Delta <= 0 && Vimp < 0;

            if (Active)
            {
                // Location of cg relative to contact
                Vector3 c_A = pos_b + c - r_A;
                Matrix3 cx = c_A.CrossOp();
                // Move inverse inertia to contact point
                Matrix3 M_inv = (1 / rb.Mass) - (cx * I_inv * cx);
                double m_eff = 1 / Dot(Normal, M_inv * Normal);
                Jn = -(1 + Epsilon) * m_eff * Vimp;

                Slip = v_A - Vimp * Normal;
                Vslip = Slip.Magnitude;
                if (Abs(Vslip) < tiny)
                {
                    Slip = o_;
                    Vslip = 0.0;
                    Js = 0.0;
                }
                else
                {
                    // Find direction friction
                    Slip = -Slip / Vslip;
                    double m_slip = 1 / Dot(Slip, M_inv * Slip);
                    Js = m_slip * Vslip;
                    Js = Min(Mu * Jn, Js);
                }

                Vector3 Jimp = Jn * Normal + Js * Slip;
                p += Jimp;
                L_b += (r_A - pos_b) ^ Jimp;

                return new BodyState(
                    pos_b,
                    ori,
                    p,
                    L_b);
            }
            else
            {
                Slip = o_;
                Vslip = 0.0;
                Jn = 0.0;
                Js = 0.0;

                return Y;
            }
        }
    }
}
