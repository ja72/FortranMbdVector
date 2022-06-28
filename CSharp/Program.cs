using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using System.IO;

namespace JA
{
    using JA.LinearAlgebra;
    using JA.Dynamics;
    using JA.Geomtery;

    using static DoubleConstants;
    using static LinearAlgebra.LinearAlgebra;

    internal class Program
    {
        public const int ReportCount = 36;
        public static readonly Vector3 Gravity = -10 * k_;

        static void Main(string[] args)
        {
            //MbdSolver(9*36);
            MbdSolver(360 * 1000);

        }

        public static IReadOnlyList<Results> MbdSolver(int n_steps, double endTime = 2.0, double angularResolutionDegrees = 0.5)
        {
            const double L = 90 * mm;
            const double D = 60 * mm;
            const double ε = 0.87;
            const double μ = 0.15;
            Console.WriteLine($" CSharp {n_steps,12} steps.");
            Console.WriteLine();

            RigidBody rb = new RigidBody(Shapes.CylinderAlongZ(L, D / 2), mass: 2.0);

            Vector3 pos = L / 2 * k_;
            Quaternion ori = q_eye;
            Vector3 omg = 5 * j_;
            Vector3 vee = (-1.0) * k_;

            rb.SetPoseAtCg(pos, ori);
            rb.SetMotionAtCg(vee, omg);

            ContactPlane floor = new ContactPlane(o_, k_, ε, μ);

            Simulation mbd = new Simulation(rb, floor, Gravity);

            double h = rb.EstMaxTimeStep(n_steps, endTime, angularResolutionDegrees);
            endTime = n_steps * h;

            Console.WriteLine(GetTableHeader());

            bool output = false;
            StringWriter io = null; 
            if (output)
            {
                io = new StringWriter();
                io.WriteLine(GetCsvHeader());
            }
            mbd.StepTaken += (s, ev) =>
            {
                var sim = ev.Simulation;
                if (sim.Step == n_steps || (sim.Step % (n_steps / ReportCount) == 0))
                {
                    Console.WriteLine(GetTableRow(sim));
                }
                io?.WriteLine(GetCsvRow(sim));
            };

            var sw = Stopwatch.StartNew();
            mbd.Run(endTime, n_steps);
            sw.Stop();

            if (output)
            {
                File.WriteAllText("results-cs.csv", io.ToString());

                io.Close();
            }
            float cpu = (float)sw.Elapsed.TotalSeconds;
            float ops = n_steps / cpu;

            Console.WriteLine();
            Console.WriteLine($" steps={n_steps} time={cpu} kops={ops / 1000}");

            return mbd.History;
        }


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

        public static string GetTableRow(Simulation simulation)
        {
            const string fmt = "{0,9:g} {1,9:F3}|{2,9:F5}{3,9:F5}{4,9:F5}|{5,9:F4}{6,9:F4}{7,9:F4}{8,9:F4}|{9,9:F5}{10,9:F5}{11,9:F5}|{12,10:F4}{13,10:F4}{14,10:F4}|{15,8:F4}";

            int i = simulation.Step;
            double time = simulation.Time;
            double Jn = simulation.Contact.Jn;
            (Vector3 position, Quaternion orientation, Vector3 momentum, Vector3 angular) = simulation.Current;
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
        public static string GetCsvRow(Simulation simulation)
        {
            int i = simulation.Step;
            double time = simulation.Time;
            double Jn = simulation.Contact.Jn;
            (Vector3 position, Quaternion orientation, Vector3 momentum, Vector3 angular) = simulation.Current;
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

    }
}
