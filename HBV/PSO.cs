using System.Diagnostics;

namespace HBV
{
    public class PSO
    {
        private readonly int numParticles;
        private readonly int dimensions;
        private readonly float[][] positions;        // [particle_idx * dimensions + dim_idx]
        private readonly float[][] velocities;       // [particle_idx * dimensions + dim_idx] 
        private readonly float[][] personalBest;     // [particle_idx * dimensions + dim_idx]
        private readonly float[] personalBestVal;  // [particle_idx]
        private readonly float[] globalBest;       // [dimensions]

        private readonly float[] searchSpaceLowerLimit;
        private readonly float[] searchSpaceUpperLimit;



        public float[] GetCopyGlobalBest()
        {
            float[] copy = new float[dimensions];
            lock (globalBest)
            {
                Array.Copy(globalBest, copy, dimensions);
            }
            return copy;
        }


        private float globalBestVal;
        private readonly Func<float[], float> objectiveFunction;


        private const float W = 0.729f;      // Inertia weight
        private const float C1 = 1.49445f;   // Cognitive parameter
        private const float C2 = 1.49445f;   // Social parameter

        public PSO(int dimensions, int numParticles, Func<float[], float> objectiveFunction, float[] searchSpaceLowerLimit, float[] searchSpaceUpperLimit)
        {
            if(searchSpaceLowerLimit.Length != dimensions || searchSpaceUpperLimit.Length != dimensions)
            {
                throw new ArgumentException("Search space limits must match number of dimensions");
            }

            this.numParticles = numParticles;
            this.dimensions = dimensions;
            this.objectiveFunction = objectiveFunction;

            positions = Alloc(numParticles, dimensions); // new float[numParticles * dimensions];
            velocities = Alloc(numParticles, dimensions); //new float[numParticles * dimensions];
            personalBest = Alloc(numParticles, dimensions); //new float[numParticles * dimensions];
            personalBestVal = new float[numParticles];
            globalBest = new float[dimensions];

            this.searchSpaceLowerLimit = searchSpaceLowerLimit;
            this.searchSpaceUpperLimit = searchSpaceUpperLimit;

            // Initialize with worst possible value
            globalBestVal = float.MaxValue;

            Initialize();
        }

        private float[][] Alloc(int firstDim, int secondDim)
        {
            float[][] result = new float[firstDim][];
            for (int i = 0; i < firstDim; i++)
            {
                result[i] = new float[secondDim];
            }
            return result;
        }

        public void Initialize()
        {
            //for (int i = 0; i < numParticles; i++)
            Parallel.For(0, numParticles, i =>
            {
                uint seed = (uint)(0 * numParticles * dimensions + i * dimensions + 1);

                for (int d = 0; d < dimensions; d++)
                {
                    //int idx = i * dimensions + d;
                    positions[i][d] = (RandomSlim.RandomFloat(ref seed) * (searchSpaceUpperLimit[d] - searchSpaceLowerLimit[d]) + searchSpaceLowerLimit[d]);
                    velocities[i][d] = (RandomSlim.RandomFloat(ref seed) * (searchSpaceUpperLimit[d] - searchSpaceLowerLimit[d]) * 0.1f);
                    personalBest[i][d] = positions[i][d];
                }

                //float[] particlePosition = new float[dimensions];
                //Array.Copy(positions, i * dimensions, particlePosition, 0, dimensions);
                personalBestVal[i] = objectiveFunction(positions[i]);

                if (float.IsNaN(personalBestVal[i]))
                    throw new Exception();

                if (personalBestVal[i] < globalBestVal)
                {
                    lock (globalBest)
                    {
                        if (personalBestVal[i] < globalBestVal)
                        {
                            globalBestVal = personalBestVal[i];
                            Array.Copy(positions[i], globalBest, dimensions);
                        }
                    }
                }
            });
        }

        public void IterateNTimes(int n = 10000)
        {
            //ParallelOptions options = new ParallelOptions
            //{
            //    MaxDegreeOfParallelism = Environment.ProcessorCount*4 // Limit to 4 threads
            //};

            Stopwatch watch = new Stopwatch();
            watch.Start();


            int total = n * numParticles;

            int atomicCounter = 0;

            //for (int i = 0; i < total; i++)
            Parallel.For(0, total, j =>
            {
                int i = Interlocked.Add(ref atomicCounter, 1) - 1;

                int particleIndex = i % numParticles;
                int iterationIndex = i / numParticles;
                Iterate(particleIndex, iterationIndex);

                if (particleIndex == 0)
                {
                    string time = (watch.ElapsedMilliseconds / 1000.0).ToString("0.000");
                    Console.WriteLine($"PSO Iteration {iterationIndex}, elapsed seconds {time}, current best score "+globalBestVal.ToString("0.000000"));
                    watch.Restart();
                }
            });
            
        }

        public void Iterate(int particleIndex, int iterationIndex)
        {
            int i = particleIndex;

            uint seed = (uint)(iterationIndex * numParticles * dimensions + i * dimensions + 1);

            for (int d = 0; d < dimensions; d++)
            {
                //int idx = i * dimensions + d;
                float r1 = RandomSlim.RandomFloat(ref seed);
                float r2 = RandomSlim.RandomFloat(ref seed);

                velocities[i][d] = W * velocities[i][d] +
                                C1 * r1 * (personalBest[i][d] - positions[i][d]) +
                                C2 * r2 * (globalBest[d] - positions[i][d]);

                positions[i][d] += velocities[i][d];

                positions[i][d] = Clamp(positions[i][d], searchSpaceLowerLimit[d], searchSpaceUpperLimit[d]);
            }

            //float[] particlePosition = new float[dimensions];
            //Array.Copy(positions, i * dimensions, particlePosition, 0, dimensions);
            float value = objectiveFunction(positions[i]);

            if (value < personalBestVal[i])
            {
                personalBestVal[i] = value;
                Array.Copy(positions[i], personalBest[i], dimensions);

                if (value < globalBestVal)
                {
                    lock (globalBest)
                    {
                        if (value < globalBestVal)
                        {
                            globalBestVal = value;
                            Array.Copy(positions[i], globalBest, dimensions);
                        }
                    }
                }
            }

        }

        private float Clamp(float value, float lower, float upper)
        {
            if (value < lower)
            {
                return lower;
            }
            else if (value > upper)
            {
                return upper;
            }
            return value;
        }

        //public float[] GetGlobalBest() => globalBest;
        //public float GetGlobalBestValue() => globalBestVal;
    }
}
