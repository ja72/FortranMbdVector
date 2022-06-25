using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace JA.Dynamics
{
    using static Math;
    using static DoubleConstants;
    using static Dynamics;

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
            Current = Contact.Calculate(Body, Current);
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

            Y_next.NormalizeOrientation();

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
                Current = Contact.Calculate(Body, _next);
                Time += h_next;
                Step++;
                _history.Add(new Results(Body, Step, Time, Current));

                StepTaken?.Invoke(this, new SimulationEventArgs(this));
            }
        }
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
