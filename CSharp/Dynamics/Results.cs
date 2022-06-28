using System;

namespace JA.Dynamics
{
    using JA.LinearAlgebra;

    public readonly struct Results
    {
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
    }
}
