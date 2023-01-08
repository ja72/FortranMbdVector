using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace JA.Dynamics
{
    using JA.Geomtery;
    using JA.LinearAlgebra;

    using static Math;
    using static DoubleConstants;
    using static LinearAlgebra.LinearAlgebra;

    public class RigidBody
    {
        public RigidBody(Shape shape, double mass = 0, double density = 0)
        {

            if (mass == 0 && density == 0)
            {
                density = 1000.0;
            }
            if (mass == 0)
            {
                mass = density * shape.Volume;
            }
            Shape = shape;
            Mass = mass;
            Cg = shape.Center;
            var (v_1, v_2, v_3) = shape.Vmmoi;
            MMoi = (mass * v_1, mass * v_2, mass * v_3);

            InitialPosition = Vector3.Zero;
            InitialOrientation = Quaternion.Identity;
            InitialVelocity = Vector3.Zero;
            InitialOmega = Vector3.Zero;
        }
        public static bool NormalizeOrientation { get; set; } = true;
        public Shape Shape { get; }
        public double Mass { get; }
        public Vector3 Cg { get; }
        public (double I_1, double I_2, double I_3) MMoi { get; }
        public Vector3 InitialPosition { get; set; }
        public Quaternion InitialOrientation { get; set; }
        public Vector3 InitialVelocity { get; set; }
        public Vector3 InitialOmega { get; set; }

        public double EstMinTimeStep(int n_steps, double endTime)
        {
            double h = endTime / n_steps;
            double omg = InitialOmega.Magnitude;
            if (omg > tiny)
            {
                h = Min(h,  0.25 * deg / omg);
            }
            return h;
        }
        public double EstMaxTimeStep(int n_steps, double endTime)
        {
            double h = endTime / n_steps;
            double omg = InitialOmega.Magnitude;
            if (omg > tiny)
            {
                h = Min(h,  15 * deg / omg);
            }
            return h;
        }

        public void SetPoseAtCg(Vector3 positionOfCg, Quaternion orientation)
        {
            Vector3 c = orientation.Rotate(Cg);
            InitialPosition = positionOfCg - c;
            InitialOrientation = orientation;
        }
        public void SetMotionAtCg(Vector3 velocityAtCg, Vector3 omega)
        {
            Vector3 c = InitialOrientation.Rotate(Cg);
            InitialVelocity = velocityAtCg + (c ^ omega);
            InitialOmega = omega;
        }

        public BodyState GetStateFromMotion(Vector3 positionAtRef, Quaternion orientation, Vector3 veocityAtRef, Vector3 omega)
        {
            //tex: Get momentum at reference point from motion
            // $$\begin{aligned}\vec{p} & =m\left(\vec{v}_{b}+\vec{\omega}\times\vec{c}\right)\\
            //\vec{L}_{b} & ={\rm I}_{c}\vec{\omega}+\vec{c}\times\vec{p}
            //\end{aligned}$$

            Matrix3 R = orientation.ToRotation();
            Vector3 c = R * Cg;
            var (I_1, I_2, I_3) = MMoi;
            Matrix3 I_c = R * Matrix3.Diagonal(I_1, I_2, I_3) * R.Transpose();
            Vector3 momentum = Mass * (veocityAtRef - (c ^ omega));
            Vector3 angularAtRef = I_c * omega + (c ^ momentum);

            return new BodyState(
                positionAtRef,
                orientation,
                momentum,
                angularAtRef);
        }
        public BodyState GetInitialState()
        {
            return GetStateFromMotion(
                InitialPosition,
                InitialOrientation,
                InitialVelocity,
                InitialOmega);
        }

        public (Vector3 velocity, Vector3 omega) GetMotionFromState(BodyState Y)
        {
            //tex: Get motion from momentum
            // $$\begin{aligned}\vec{v}_{b} & =\tfrac{1}{m}\vec{p}-\vec{\omega}\times\vec{c}\\
            //\vec{\omega} & ={\rm I}_{c}^{-1}\left(\vec{L}_{b}-\vec{c}\times\vec{p}\right)
            //\end{aligned}$$

            (_, Quaternion orientation, Vector3 p, Vector3 L_b) = Y;
            Matrix3 R = orientation.ToRotation();
            Vector3 c = R * Cg;
            var (I_1, I_2, I_3) = MMoi;
            Matrix3 I_inv = R * Matrix3.Diagonal(1 / I_1, 1 / I_2, 1 / I_3) * R.Transpose();
            Vector3 omega = I_inv * (L_b - (c ^ p));
            Vector3 vee_b = p / Mass - (omega ^ c);
            return (vee_b, omega);
        }

        public BodyState GetRate(double time, BodyState Y, Vector3 gravity)
        {
            //tex: Get time derivative of state
            // $$\frac{{\rm d}}{{\rm d}t}\begin{Bmatrix}\vec{r}_{b}\\
            //q\\
            //\vec{p}\\
            //\vec{L}_{b}
            //\end{Bmatrix}=\begin{Bmatrix}\vec{v}_{b}\\
            //\tfrac{1}{2}\,\omega\otimes q\\
            //\vec{F}\\
            //\vec{\tau}_{b}-\vec{v}_{b}\times\vec{p}
            //\end{Bmatrix}$$

            if (NormalizeOrientation)
            {
                Y.NormalizeOrientation();
            }
            (_, Quaternion orientation, Vector3 p, Vector3 L_b) = Y;
            Matrix3 R = orientation.ToRotation();
            Vector3 c = R * Cg;
            var (I_1, I_2, I_3) = MMoi;
            Matrix3 I_inv = R * Matrix3.Diagonal(1 / I_1, 1 / I_2, 1 / I_3) * R.Transpose();
            Vector3 omega = I_inv * (L_b - (c ^ p));
            Vector3 vee_b = p / Mass - (omega ^ c);
            Quaternion qp = 0.5 * Quaternion.FromVector(omega) * orientation;
            Vector3 F = Mass * gravity;
            Vector3 τ_b = c ^ F;
            Vector3 dLdt = τ_b - (vee_b ^ p);

            return new BodyState(
                vee_b,
                qp,
                F,
                dLdt);
        }

    }

}