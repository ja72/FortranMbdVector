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

    public class Simulation
    {
        public event EventHandler<SimulationEventArgs> StepTaken;

        readonly List<Results> _history;

        public Simulation(RigidBody body, ContactPlane contact, Vector3 gravity)
        {
            Body = body;
            Contact = contact;
            _history = new List<Results>();
            Gravity = gravity;
            Reset();
        }

        public Vector3 Gravity { get; set; }
        public RigidBody Body { get; }
        public ContactPlane Contact { get; }
        public bool HasContact { get => Contact != null; }
        public double Time { get; private set; }
        public int Step { get; private set; }
        public IReadOnlyList<Results> History => _history;
        public BodyState Current { get; private set; }

        public void Reset()
        {
            Step = 0;
            Time = 0;
            _history.Clear();
            Current = Body.GetInitialState();
            if (HasContact)
            {
                Current = Contact.Calculate(Body, Current);
            }
            _history.Add(new Results(Body, Step, Time, Current));

            StepTaken?.Invoke(this, new SimulationEventArgs(this));
        }
        BodyState Integrate(double timeStep) => Integrate(Time, Current, timeStep);
        BodyState Integrate(double time, BodyState Y, double timeStep)
        {
            BodyState k0 = Body.GetRate(time, Y, Gravity);
            BodyState Y1 = BodyState.AddScale(Y, k0, 1, timeStep / 2);
            BodyState k1 = Body.GetRate(time+timeStep/2, Y1, Gravity);
            BodyState Y2 = BodyState.AddScale(Y, k1, 1, timeStep / 2);
            BodyState k2 = Body.GetRate(time+timeStep/2, Y2, Gravity);
            BodyState Y3 = BodyState.AddScale(Y, k2, 1, timeStep);
            BodyState k3 = Body.GetRate(time+timeStep, Y3, Gravity);

            BodyState Y_next = Y;

            Y_next.AddScale(k0, timeStep / 6);
            Y_next.AddScale(k1, timeStep / 3);
            Y_next.AddScale(k2, timeStep / 3);
            Y_next.AddScale(k3, timeStep / 6);

            if (RigidBody.NormalizeOrientation)
            {
                Y_next.NormalizeOrientation();
            }

            return Y_next;
        }


        public void Run(double endTime, int n_steps)
        {
            StepTaken?.Invoke(this, new SimulationEventArgs(this));
            double h = (endTime - Time) / n_steps;
            n_steps += Step;
            while (Step<n_steps)
            {
                double h_next = Min(h, endTime - Time);
                var _next = Integrate(h_next);
                if (HasContact && Contact.IsOk)
                {
                    Current = Contact.Calculate(Body, _next);
                }
                else
                {
                    Current = _next;
                }
                Time += h_next;
                Step++;
                _history.Add(new Results(Body, Step, Time, Current));

                StepTaken?.Invoke(this, new SimulationEventArgs(this));
            }
        }

        #region Formatting
        public static string GetTableHeader()
        {
            const string fmt = "{0,9} {1,9}|{2,9}{3,9}{4,9}|{5,9}{6,9}{7,9}{8,9}|{9,9}{10,9}{11,9}|{12,10}{13,10}{14,10}|{15,8}";
            var line1 = new object[] { "#", "time",
                "pos-x", "pos-y", "pos-z",
                "ori-x", "ori-y", "ori-z", "ori-w",
                "mom-x", "mom-y", "mom-z",
                "ang-x", "ang-y", "ang-z",
                "J" };
            var line2 = new object[] { "", "[s]",
                "[m]", "[m]", "[m]",
                "[]", "[]", "[]", "[]",
                "[kg m/s]", "[kg m/s]", "[kg m/s]",
                "[g m^2/s]", "[g m^2/s]", "[g m^2/s]",
                "[N s]" };
            var line3 = new object[] { "---",  "---",
                "---", "---", "---",
                "---", "---", "---", "---",
                "---", "---", "---",
                "---", "---", "---",
                "---" };

            var sb = new StringBuilder();
            sb.AppendFormat(fmt, line1);
            sb.AppendLine();
            sb.AppendFormat(fmt, line2);
            sb.AppendLine();
            sb.AppendFormat(fmt, line3);

            return sb.ToString();
        }

        public string GetTableRow()
        {
            const string fmt = "{0,9:g} {1,9:F3}|{2,9:F5}{3,9:F5}{4,9:F5}|{5,9:F4}{6,9:F4}{7,9:F4}{8,9:F4}|{9,9:F5}{10,9:F5}{11,9:F5}|{12,10:F4}{13,10:F4}{14,10:F4}|{15,8:F4}";

            int i = this.Step;
            double time = this.Time;
            double Jn = this.Contact.Jn;
            (Vector3 position, Quaternion orientation, Vector3 momentum, Vector3 angular) = this.Current;
            var line = new object[] {
                i, time,
                position.X, position.Y, position.Z,
                orientation.Vector.X, orientation.Vector.Y, orientation.Vector.Z, orientation.Scalar,
                momentum.X, momentum.Y, momentum.Z,
                1000*angular.X, 1000*angular.Y, 1000*angular.Z,
                Jn
            };
            return string.Format(fmt, line);
        }

        public static string GetCsvHeader()
        {
            var line = new object[] { "#", "time [s]",
                "pos-x [m]", "pos-y [m]", "pos-z [m]",
                "ori-x []", "ori-y []", "ori-z []", "ori-w []",
                "mom-x [kg m/s]", "mom-y [kg m/s]", "mom-z [kg m/s]",
                "ang-x [g m^2/s]", "ang-y [g m^2/s]", "ang-z [g m^2/s]",
                "J [N s]" };
            return string.Join(",", line);
        }
        public string GetCsvRow()
        {
            int i = this.Step;
            double time = this.Time;
            double Jn = this.Contact.Jn;
            (Vector3 position, Quaternion orientation, Vector3 momentum, Vector3 angular) = this.Current;
            var line = new object[] {
                i, time,
                position.X, position.Y, position.Z,
                orientation.Vector.X, orientation.Vector.Y, orientation.Vector.Z, orientation.Scalar,
                momentum.X, momentum.Y, momentum.Z,
                1000*angular.X, 1000*angular.Y, 1000*angular.Z,
                Jn
            };
            return string.Join(",", line);
        }

        #endregion
    }

    public class SimulationEventArgs : EventArgs
    {
        public SimulationEventArgs(Simulation simulation)
        {
            Simulation = simulation;
        }
        public Simulation Simulation { get; }
    }
}
