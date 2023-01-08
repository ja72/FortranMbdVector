using System;

namespace JA.Dynamics
{
    using System.Text;
    using JA.LinearAlgebra;

    public readonly struct Results
    {
        public Results(Simulation simulation) 
            : this(simulation.Body, simulation.Step, simulation.Time, simulation.Current)
        {
        }
        public Results(RigidBody body, int step, double time, BodyState state) : this()
        {
            Body = body;
            Step = step;
            Time = time;

            (Vector3 position, Quaternion orientation, Vector3 momentum, Vector3 angular) = state;
            (Vector3 velocity, Vector3 omega) = body.GetMotionFromState(state);

            Position = position;
            Orientation = orientation;
            Momentum = momentum;
            Angular = angular;
            Velocity = velocity;
            Omega = omega;
        }

        public RigidBody Body { get; }
        public int Step { get; }
        public double Time { get; }
        public Vector3 Position { get; }
        public Quaternion Orientation { get; }
        public Vector3 Momentum { get; }
        public Vector3 Angular { get; }
        public Vector3 Velocity {get; }
        public Vector3 Omega { get; }

        #region Formatting

        public static string GetTableHeader()
        {
            const string fmt = "{0,9} {1,9}|{2,9}{3,9}{4,9}|{5,9}{6,9}{7,9}{8,9}|{9,9}{10,9}{11,9}|{12,10}{13,10}{14,10}";
            var line1 = new object[] { "#", "time",
                "pos-x", "pos-y", "pos-z",
                "ori-x", "ori-y", "ori-z", "ori-w",
                "mom-x", "mom-y", "mom-z",
                "ang-x", "ang-y", "ang-z",
                 };
            var line2 = new object[] { "", "[s]",
                "[m]", "[m]", "[m]",
                "[]", "[]", "[]", "[]",
                "[kg m/s]", "[kg m/s]", "[kg m/s]",
                "[g m^2/s]", "[g m^2/s]", "[g m^2/s]",
                };
            var line3 = new object[] { "---",  "---",
                "---", "---", "---",
                "---", "---", "---", "---",
                "---", "---", "---",
                "---", "---", "---",
                };

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
            const string fmt = "{0,9:g} {1,9:F3}|{2,9:F5}{3,9:F5}{4,9:F5}|{5,9:F4}{6,9:F4}{7,9:F4}{8,9:F4}|{9,9:F5}{10,9:F5}{11,9:F5}|{12,10:F4}{13,10:F4}{14,10:F4}";

            int i = Step;
            double time = Time;
            (Vector3 position, Quaternion orientation, Vector3 momentum, Vector3 angular) = (Position, Orientation, Momentum, Angular);
            var line = new object[] {
                i, time,
                position.X, position.Y, position.Z,
                orientation.Vector.X, orientation.Vector.Y, orientation.Vector.Z, orientation.Scalar,
                momentum.X, momentum.Y, momentum.Z,
                1000*angular.X, 1000*angular.Y, 1000*angular.Z,                
            };
            return string.Format(fmt, line);
        }

        public static string GetCSVHead()
        {
            var objs = new object[] {
                "Step",
                "Time",
                Vector3.GetCSVHead("r"),
                Quaternion.GetCSVHead("q"),
                Vector3.GetCSVHead("mom"),
                Vector3.GetCSVHead("ang"),
                Vector3.GetCSVHead("vee"),
                Vector3.GetCSVHead("omg"),
            };
            return string.Join(",", objs);

        }
        public string ToCSVRow()
        {
            var objs = new object[] {
                Step,
                Time,
                Position.ToCSVRow(),
                Orientation.ToCSVRow(),
                Momentum.ToCSVRow(),
                Angular.ToCSVRow(),
                Velocity.ToCSVRow(),
                Omega.ToCSVRow(),
            };
            return string.Join(",", objs);
        } 
        #endregion
    }
}
