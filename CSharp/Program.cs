using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using System.IO;
    using System.Runtime.InteropServices;

namespace JA
{
    using JA.LinearAlgebra;
    using JA.Dynamics;
    using JA.Geomtery;

    using static DoubleConstants;
    using static LinearAlgebra.LinearAlgebra;

    internal static class Fortran
    {
        [DllImport("FortranMbdVector.exe", CharSet = CharSet.Ansi, CallingConvention = CallingConvention.Cdecl)]
        //internal static extern void sim_export([In] vec3 gravity, int n_steps, double t_end, string f_export, ref int f_export_length);
        internal static extern void sim_export([In] ref double gravity);
    }

    internal class Program
    {
        public const int ReportCount = 36;
        public static readonly Vector3 Gravity = -10 * k_;

        static void Main(string[] args)
        {
            //MbdSolver(9*36);
            //MbdSolver(360 * 2250, 2.0);

            //var results = FreeSolver(360 * 2250, 2.0);

            FortranSolver(360, 2.0);
        }

        public static void FortranSolver(int n_steps, double endTime)
        {
            string f_export = "results-cs-f.csv";
            int f_export_length = f_export.Length;
            double[] gravity = new double[] { 0, 0, -10 };
            //Fortran.sim_export(gravity, n_steps, endTime, f_export, ref f_export_length);
            Fortran.sim_export(ref gravity[0]);
        }

        public static IReadOnlyList<Results> MbdSolver(int n_steps, double endTime)
        {
            const double ε = 0.87;
            const double μ = 0.15;
            Console.WriteLine($" CSharp {n_steps,12} steps.");
            Console.WriteLine();

            const double L = 90 * mm;
            const double D = 60 * mm;

            RigidBody rb = new RigidBody(Shapes.CylinderAlongZ(L, D / 2), mass: 2.0);

            Vector3 pos = L / 2 * k_;
            Quaternion ori = q_eye;
            Vector3 omg = 5 * j_;
            Vector3 vee = (-1.0) * k_;

            rb.SetPoseAtCg(pos, ori);
            rb.SetMotionAtCg(vee, omg);

            ContactPlane floor = new ContactPlane(o_, k_, ε, μ);

            Simulation mbd = new Simulation(rb, floor, Gravity);

            double h = rb.EstMinTimeStep(n_steps, endTime);
            endTime = n_steps * h;

            Console.WriteLine(Simulation.GetTableHeader());

            bool output = false;
            StringWriter io = null;
            if (output)
            {
                io = new StringWriter();
                io.WriteLine(Simulation.GetCsvHeader());
            }
            mbd.StepTaken += (s, ev) =>
            {
                var sim = ev.Simulation;
                if (sim.Step == n_steps || (sim.Step % (n_steps / ReportCount) == 0))
                {
                    Console.WriteLine(sim.GetTableRow());
                }
                io?.WriteLine(sim.GetCsvRow());
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

        public static IReadOnlyList<Results> FreeSolver(int n_steps, double endTime)
        {
            const double L = 90 * mm;
            const double D = 60 * mm;

            RigidBody rb = new RigidBody(Shapes.CylinderAlongZ(L, D / 2), mass: 2.0);

            Vector3 pos = L / 2 * k_;
            Quaternion ori = q_eye;
            Vector3 omg = 5 * j_ - 1 * k_;
            Vector3 vee = (0.0) * k_;

            rb.SetPoseAtCg(pos, ori);
            rb.SetMotionAtCg(vee, omg);
            Simulation mbd = new Simulation(rb, null, Vector3.Zero);
            mbd.Run(endTime, 360);
            Console.WriteLine(Results.GetTableHeader());
            foreach (var item in mbd.History)
            {
                if (item.Step % 15 == 0)
                {
                    Console.WriteLine(item.GetTableRow());
                }
            }

            Console.WriteLine("Checking Accuracy of Simulation.");

            bool output = true;
            StringWriter io = null;
            if (output)
            {
                io = new StringWriter();
                io.WriteLine(Results.GetCSVHead());
            }

            List<Results> results = new List<Results>();
            Console.WriteLine($"End Time = {endTime}");
            Console.WriteLine();
            Console.WriteLine($"{"Steps",9} {"Angle",24} {"Omega",24}");
            Console.WriteLine($"{"",9} {"[deg]",24} {"[rad/s]",24}");
            Quaternion q_expected = new Quaternion(0.45393475, -0.830921987, -0.064323953, 0.315236931);
            Vector3 ω_expected = new Vector3(0.291136175, 5.089648811, -0.103511894);
            do
            {
                mbd.Reset();
                mbd.Run(endTime, n_steps);
                Results result = new Results(mbd);
                results.Add(result);
                if (output)
                {
                    io.WriteLine(result.ToCSVRow());
                }
                var dq = result.Orientation * q_expected.Inverse();
                var dw = result.Omega - ω_expected;
                Console.WriteLine($"{result.Step,9} {dq.GetAxisAngle().angle / deg,24} {dw.Magnitude,24}");
                n_steps /= 2;
            } while (n_steps>1);

            if (output)
            {
                File.WriteAllText("step-doe-cs.csv", io.ToString());

                io.Close();
            }

            return results.AsReadOnly();
        }


    }
}
