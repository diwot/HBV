using HBV;
using HBV.InputOutput;
using System.Runtime.CompilerServices;

namespace HBV
{
    public class FloatVec
    {
        private readonly float[] data;

        public FloatVec(int size)
        {
            data = new float[size];
        }

        public FloatVec(float[] array)
        {
            data = array;
        }

        public static implicit operator float[](FloatVec vec)
        {
            return vec.data;
        }

        public float this[int index]
        {
            get { return data[index]; }
            set { data[index] = value; }
        }

        public int Length => data.Length;

        public static FloatVec operator +(FloatVec a, FloatVec b)
        {
            if (a.Length != b.Length)
                throw new ArgumentException("Vectors must be same length");

            var result = new FloatVec(a.Length);
            for (int i = 0; i < a.Length; i++)
            {
                result[i] = a[i] + b[i];
            }
            return result;
        }

        public static FloatVec operator -(FloatVec a, FloatVec b)
        {
            if (a.Length != b.Length)
                throw new ArgumentException("Vectors must be same length");

            var result = new FloatVec(a.Length);
            for (int i = 0; i < a.Length; i++)
            {
                result[i] = a[i] - b[i];
            }
            return result;
        }

        public static FloatVec operator *(FloatVec a, float scalar)
        {
            var result = new FloatVec(a.Length);
            for (int i = 0; i < a.Length; i++)
            {
                result[i] = a[i] * scalar;
            }
            return result;
        }

        public static FloatVec operator *(float scalar, FloatVec a)
        {
            return a * scalar;
        }

        public static FloatVec operator /(FloatVec a, float scalar)
        {
            if (scalar == 0)
                throw new DivideByZeroException();

            var result = new FloatVec(a.Length);
            for (int i = 0; i < a.Length; i++)
            {
                result[i] = a[i] / scalar;
            }
            return result;
        }
    }



    public class HBVState
    {
        public float sm1;    // initial storage soil moisture
        public float sw1;    // initial storage surface water 
        public float gw1;    // initial storage ground water
        public float sp1;    // initial storage snow pack
        public float mw1;    // initial storage melt water

        public float[] sp;   // storage snow pack array
        public float[] mw;   // storage melted water array
        public float[] sm;   // storage soil moisture array
        public float[] sw;   // storage surface water array
        public float[] gw;   // storage ground water array
    }



    public class HBVOptions
    {
        public float TCON;     // time constant
        public float CEVPFO;   // correction factor for evaporation in forest areas
        public float TT;       // threshold temperature for snowfall
        public float TTI;      // threshold temperature for snow melt
        public float DTTM;     // degree day factor
        public float CFMAX;    // degree-day factor
        public float FOCFMAX;  // degree day factor, rate of snowmelt
        public float CFR;      // refreezing factor
        public float WHC;      // water holding capacity
    }


    public class MeteoData
    {
        public int ElevationBand;
        public float[] TairElev;
        public float[] PrecipElev;
        public float[] PetElev;
    }

    public class ElevationBandData
    {
        public float[] r;      // rainfall
        public float[] ssp;    // storage snow pack
        public float[] ssw;    // storage surface water
        public float[] ssm;    // storage soil moisture
        public float[] sgw;    // storage ground water
        public float[] qd;     // direct discharge
        public float[] qc;     // capillary transport
        public float[] qin;    // infiltration
        public float[] qf;     // fast run-off
        public float[] qs;     // slow run-off
        public float[] qr;     // river discharge

        public ElevationBandData() { }

        public ElevationBandData(int length)
        {
            r= new float[length];
            ssp= new float[length + 1];
            ssw= new float[length + 1];
            ssm= new float[length + 1];
            sgw= new float[length + 1];
            qd= new float[length];
            qc= new float[length];
            qin= new float[length];
            qf= new float[length];
            qs= new float[length];
            qr= new float[length];
        }
    }

    public class CatchmentInfo
    {
        public float[] Area_m2;
        public float[] rel_Area;
        public float[] Elevation_bottom;
        public float[] Elevation_top;
        public float Total_Area_m2;
    }

    public class Global
    {
        public float DTime;
        public float[] Tmean;
        public float[] PET;
        public float[] Precip;
        public float[] Discharge;
        public CatchmentInfo Catchmentinfo;
    }

    public class HBVModel
    {
        private static float Corr(float[] x, float[] y)
        {
            // Pearson correlation coefficient
            float meanX = x.Average();
            float meanY = y.Average();

            float sum1 = 0f;
            float sum2 = 0f;
            float sum3 = 0f;

            for (int i = 0; i < x.Length; i++)
            {
                float xDiff = x[i] - meanX;
                float yDiff = y[i] - meanY;
                sum1 += xDiff * yDiff;
                sum2 += xDiff * xDiff;
                sum3 += yDiff * yDiff;
            }

            return sum1 / (float)Math.Sqrt(sum2 * sum3);
        }

        private static float StdDev(float[] values)
        {
            float mean = values.Average();
            float sumSquaredDiffs = values.Sum(x => (x - mean) * (x - mean));
            return MathF.Sqrt(sumSquaredDiffs / values.Length);
        }

        public static float KGE_HBV(HBVParams[] pars, int permafrostOption, float seconds, Global global)
        {
            float[] runoffMeasured = global.Discharge;

            var outputHBV = HBV_permafrost_semidistributed(pars, permafrostOption, seconds, global);
            float[] runoffModelled = outputHBV.qr;

            // Get valid indices (where runoffMeasured is not NaN)
            var validIndices = new List<int>();
            for (int i = 0; i < runoffMeasured.Length; i++)
            {
                if (!float.IsNaN(runoffMeasured[i]))
                {
                    validIndices.Add(i);
                }
            }

            float[] qObs = validIndices.Select(i => runoffMeasured[i]).ToArray();
            float[] qSim = validIndices.Select(i => runoffModelled[i]).ToArray();

            // Compute components of KGE
            float r = Corr(qObs, qSim); // Pearson correlation coefficient
            float beta = qSim.Average() / qObs.Average(); // Bias ratio
            float alpha = StdDev(qSim) / StdDev(qObs); // Variability ratio Kling et al., 2012

            // Compute KGE score
            float KGE = 1f - (float)Math.Sqrt((r - 1f) * (r - 1f) +
                                             (beta - 1f) * (beta - 1f) +
                                             (alpha - 1f) * (alpha - 1f));

            return -KGE;
        }



        public static ElevationBandData HBV_permafrost_semidistributed(HBVParams[] pars, int permafrostOption, 
            float seconds, Global global)
        {
            CatchmentInfo Catchmentinfo = global.Catchmentinfo;

            var tm = global.Tmean;

            // Constants
            const float t_lapse = -0.8f;  // temperature lapse rate in °C/100m
            const float P_LAPSE = 1f;     // precipitation lapse rate in mm/100m
            const float PET_LAPSE = -0.6f; // pet lapse rate in °C/100m

            float[] h_elev_band = new float[Catchmentinfo.Elevation_bottom.Length];
            string[] elev_band = new string[h_elev_band.Length];
            for (int i = 0; i < h_elev_band.Length; i++)
            {
                elev_band[i] = $"eb_{i + 1}";
            }
            float[] area = Catchmentinfo.Area_m2;   // Area of each elevation band in m2
            float[] rel_area = new FloatVec(Catchmentinfo.rel_Area) / 100f; // relative area of each elevation band (sum should equal 1)

            for (int i = 0; i < h_elev_band.Length; i++)
            {
                h_elev_band[i] = (Catchmentinfo.Elevation_bottom[i] + Catchmentinfo.Elevation_top[i]) / 2f;
            }
            float h_meteo = 10f; // Height of Meteostation above sea level
            float A_total = Catchmentinfo.Total_Area_m2; // catchment area in m2
            float[] A_elev_band = area;

            if (pars.Length != elev_band.Length)
            {
                throw new ArgumentException("Number of elevation bands does not match number of parameter sets");
            }

            MeteoData[] meteoData = new MeteoData[elev_band.Length];
            ElevationBandData[] elevationBandData = new ElevationBandData[elev_band.Length];

            // Process each elevation band
            for (int i = 0; i < elev_band.Length; i++)
            {
                var pars_elev = pars[i];

                float MAXBAS = pars_elev.MAXBAS; // pars[i,15];
                float precip_corr = pars_elev.precip_corr; // pars[i,14];

                // Calculate meteorological adjustments
                float diff_h = h_elev_band[i] - h_meteo;

                // Calculate adjusted meteorological variables following MATLAB code
                float[] tair_elev = new float[tm.Length];
                float[] precip_elev = new float[global.Precip.Length];
                float[] pet_elev = new float[global.PET.Length];

                // Temperature adjustment
                for (int j = 0; j < tm.Length; j++)
                {
                    tair_elev[j] = tm[j] + (t_lapse / 100f * diff_h);
                }

                // Precipitation adjustment
                for (int j = 0; j < global.Precip.Length; j++)
                {
                    float po = global.Precip[j] * precip_corr;
                    precip_elev[j] = po;
                    if (po != 0)
                    {
                        precip_elev[j] = Math.Max(0, po + (P_LAPSE / 100f * diff_h));
                    }
                }

                // PET adjustment
                for (int j = 0; j < global.PET.Length; j++)
                {
                    float etpo = global.PET[j];
                    pet_elev[j] = etpo;
                    if (etpo != 0)
                    {
                        pet_elev[j] = Math.Max(0, etpo + (PET_LAPSE / 100f * diff_h));
                    }
                }

                // Store results in MeteoData
                meteoData[i] = new MeteoData
                {
                    ElevationBand = i,
                    TairElev = tair_elev,
                    PrecipElev = precip_elev,
                    PetElev = pet_elev
                };

                ElevationBandData outputHBV = new ElevationBandData(meteoData[i].PetElev.Length); // default;
                switch (permafrostOption)
                {
                    case 1:
                        Core.HBV_no_permafrost_optimized(meteoData[i].TairElev, meteoData[i].PrecipElev,
                            meteoData[i].PetElev, pars_elev, Catchmentinfo.Area_m2[i], seconds, outputHBV);

                        //outputHBV = Core.HBV_no_permafrost(meteoData[i].TairElev, meteoData[i].PrecipElev,
                        //    meteoData[i].PetElev, pars_elev, Catchmentinfo.Area_m2[i], seconds);
                        break;
                        //case 2:
                        //    outputHBV = HBV_permafrost_ice_storage(meteoData.TairElev, meteoData.PrecipElev,
                        //        meteoData.PetElev, pars_elev, Catchmentinfo.Area_m2[i], seconds);
                        //    break;
                        //case 3:
                        //    outputHBV = HBV_permafrost_infiltration_blocking(meteoData.TairElev, meteoData.PrecipElev,
                        //        meteoData.PetElev, pars_elev, Catchmentinfo.Area_m2[i], seconds);
                        //    break;
                        //default:
                        //    outputHBV = HBV_permafrost_combined(meteoData.TairElev, meteoData.PrecipElev,
                        //        meteoData.PetElev, pars_elev, Catchmentinfo.Area_m2[i], seconds);
                        //   break;
                }

                // Store results
                elevationBandData[i] = new ElevationBandData
                {
                    r = outputHBV.r,
                    ssp = outputHBV.ssp,
                    ssw = outputHBV.ssw,
                    ssm = outputHBV.ssm,
                    sgw = outputHBV.sgw,
                    qd = outputHBV.qd,
                    qc = outputHBV.qc,
                    qin = outputHBV.qin,
                    qf = outputHBV.qf,
                    qs = outputHBV.qs,
                    qr = Routing(outputHBV.qr, MAXBAS)
                };
            }

            return CombineElevationBandResults(elevationBandData);
        }

        private static float[] Routing(float[] Q, float MAXBAS)
        {
            // Integration step
            float step = 0.005f;

            // Create arrays for y and h
            int numSteps = (int)(MAXBAS / step) + 1;
            float[] y = new float[numSteps];
            float[] h = new float[numSteps];

            // Fill y array with values from 0 to MAXBAS
            for (int i = 0; i < numSteps; i++)
            {
                y[i] = i * step;
            }

            if (MAXBAS == 0.0)
                throw new Exception("division by zero");

            // Construct triangular weighting function
            for (int i = 0; i < numSteps; i++)
            {
                if (y[i] < MAXBAS / 2)
                {
                    h[i] = step * (y[i] * 4 / (MAXBAS * MAXBAS));
                }
                else
                {
                    h[i] = step * (4 / MAXBAS - y[i] * 4 / (MAXBAS * MAXBAS));
                }
            }

            // Adjust for non-integer MAXBAS
            int numSegments = (int)Math.Ceiling(MAXBAS) + 1;
            float[] I = new float[numSegments];

            if (MAXBAS % 1 > 0)
            {
                for (int i = 0; i < numSegments - 1; i++)
                {
                    I[i] = i * (numSteps - 1) / (MAXBAS);
                }
                I[numSegments - 1] = numSteps - 1;
            }
            else
            {
                for (int i = 0; i < numSegments; i++)
                {
                    I[i] = i * (numSteps - 1) / MAXBAS;
                }
            }

            // Calculate MAXBAS weights
            float[] MAXBAS_w = new float[I.Length - 1];
            for (int k = 1; k < I.Length; k++)
            {
                int start = (int)Math.Floor(I[k - 1]);
                int end = (int)Math.Floor(I[k]);

                float sum = 0;
                for (int j = start; j <= end; j++)
                {
                    sum += h[j];
                }
                MAXBAS_w[k - 1] = sum;
            }

            // Normalize weights
            float weightSum = MAXBAS_w.Sum();

            if (weightSum == 0.0)
                throw new Exception("division by zero");

            for (int i = 0; i < MAXBAS_w.Length; i++)
            {
                MAXBAS_w[i] /= weightSum;
            }

            // Perform convolution
            float[] result = new float[Q.Length];
            int halfWindow = MAXBAS_w.Length / 2;

            for (int i = 0; i < Q.Length; i++)
            {
                float sum = 0;
                for (int j = 0; j < MAXBAS_w.Length; j++)
                {
                    int idx = i - halfWindow + j;
                    if (idx >= 0 && idx < Q.Length)
                    {
                        sum += Q[idx] * MAXBAS_w[j];
                    }
                }
                result[i] = sum;
            }

            return result;
        }

        private static ElevationBandData CombineElevationBandResults(ElevationBandData[] elevationBandData)
        {
            // Based on MATLAB reference:
            // Sum qr across all elevation bands for each timestep
            int timeSteps = elevationBandData[0].qr.Length;
            float[] totalDischarge = new float[timeSteps];

            for (int t = 0; t < timeSteps; t++)
            {
                for (int band = 0; band < elevationBandData.Length; band++)
                {
                    totalDischarge[t] += elevationBandData[band].qr[t];
                }
            }

            // Return combined results
            return new ElevationBandData
            {
                qr = totalDischarge,
                // Copy other fields from first band or combine as needed
                r = elevationBandData[0].r,
                ssp = elevationBandData[0].ssp,
                ssw = elevationBandData[0].ssw,
                ssm = elevationBandData[0].ssm,
                sgw = elevationBandData[0].sgw,
                qd = elevationBandData[0].qd,
                qc = elevationBandData[0].qc,
                qin = elevationBandData[0].qin,
                qf = elevationBandData[0].qf,
                qs = elevationBandData[0].qs
            };
        }

        //public class DataCache
        //{
        //    public float[] p ;
        //    public float[] ta ;
        //    public float[] ss ;
        //    public float[] r ;
        //    public float[] sm ;
        //    public float[] sr ;
        //    public float[] in_;
        //    public float[] qd ;
        //    public float[] qin;
        //    public float[] etp;
        //    public float[] eta;
        //    public float[] qc ;
        //    public float[] qf ;
        //    public float[] qs ;
        //    public float[] qt ;
        //    public float[] qr ;

        //    public HBVState state;
        //    public HBVOptions O;
        //}



        //public static ElevationBandData HBV_no_permafrost_optimized2(float[] tm, float[] po, float[] etpo, HBVParams P, float A, float seconds)
        //{
        //    const float FFO = 0.0f;  // percent of forest
        //    const float FFI = 1.00f; // not forested
        //    int dtot = etpo.Length; // number of timesteps

        //    // Pre-compute constants used in loops
        //    float ttMinusHalfTti = P.TT - P.TTI * 0.5f;
        //    float ttPlusHalfTti = P.TT + P.TTI * 0.5f;
        //    float ttiInv = 1.0f / P.TTI;
        //    float tconInv = 1.0f / seconds;
        //    float fcInv = 1.0f / P.FC;
        //    float lpFcInv = 1.0f / (P.LP * P.FC);
        //    float areaScaling = A / (seconds * 1000);
        //    float cfmaxScaled = P.CFMAX * (FFO * P.FOCFMAX + FFI);
        //    float cfrScaled = P.CFR * cfmaxScaled;
        //    float evapoScaling = FFO * 1.15f + FFI;
        //    float onePlusAlfa = 1 + P.ALFA;

        //    var O = new HBVOptions
        //    {
        //        TCON = seconds,
        //        CEVPFO = 1.15f,
        //        TT = P.TT,
        //        TTI = P.TTI,
        //        DTTM = 0.0f,
        //        CFMAX = P.CFMAX,
        //        FOCFMAX = P.FOCFMAX,
        //        CFR = P.CFR,
        //        WHC = P.WHC
        //    };

        //    // Preallocate all arrays at once
        //    var state = new HBVState
        //    {
        //        sm1 = 50f,
        //        sw1 = 1f,
        //        gw1 = 40f,
        //        sp1 = 50f,
        //        mw1 = 0f,
        //        sp = new float[dtot + 1],
        //        mw = new float[dtot + 1],
        //        sm = new float[dtot + 1],
        //        sw = new float[dtot + 1],
        //        gw = new float[dtot + 1]
        //    };

        //    float[] r = new float[dtot];
        //    float[] qd = new float[dtot];
        //    float[] qc = new float[dtot];
        //    float[] qin = new float[dtot];
        //    float[] qf = new float[dtot];
        //    float[] qs = new float[dtot];
        //    float[] qr = new float[dtot];

        //    // Set initial conditions
        //    state.sp[0] = state.sp1;
        //    state.mw[0] = state.mw1;
        //    state.sm[0] = Math.Min(state.sm1, P.FC);
        //    state.sw[0] = state.sw1;
        //    state.gw[0] = state.gw1;
        //    state.sp[1] = state.sp[0];
        //    state.mw[1] = state.mw[0];
        //    state.sm[1] = state.sm[0];
        //    state.sw[1] = state.sw[0];
        //    state.gw[1] = state.gw[0];

        //    // Initial values
        //    qf[0] = P.KF * MathF.Pow(state.sw[0], onePlusAlfa);
        //    qs[0] = P.KS * state.gw[0];
        //    qr[0] = (qf[0] + qs[0]) * areaScaling;

        //    // Main simulation loop
        //    for (int t = 0; t < dtot; t++)
        //    {
        //        float ta_t = tm[t];
        //        float p_t = po[t];
        //        float ss_t, sm_t, sr_t, in_t;

        //        // Snow pack balance components
        //        if (ta_t < ttMinusHalfTti)
        //        {
        //            ss_t = p_t * (FFO + FFI);
        //            r[t] = 0;
        //        }
        //        else if (ta_t < ttPlusHalfTti)
        //        {
        //            ss_t = p_t * (ttPlusHalfTti - ta_t) * ttiInv * (FFO + FFI);
        //            r[t] = p_t * (ta_t - ttMinusHalfTti) * ttiInv;
        //        }
        //        else
        //        {
        //            ss_t = 0;
        //            r[t] = p_t;
        //        }

        //        // Snow melt and refreezing
        //        sm_t = Math.Min(Math.Max(P.CFMAX * (ta_t - ttPlusHalfTti), 0), state.sp[t]);
        //        sr_t = Math.Min(Math.Max(cfrScaled * (ttPlusHalfTti - ta_t), 0), state.mw[t]);

        //        // Meltwater and rain reservoir
        //        in_t = Math.Max(state.mw[t] + sm_t + r[t] - sr_t - O.WHC * state.sp[t], 0);

        //        // Soil moisture balance components
        //        qd[t] = Math.Max(in_t + state.sm[t] - P.FC, 0);

        //        // Amount of infiltration
        //        float smFcRatio = state.sm[t] * fcInv;
        //        qin[t] = Math.Max(MathF.Pow(smFcRatio, P.BETA) * (in_t - qd[t]), 0);

        //        // Evaporation
        //        float etp_t = etpo[t] * evapoScaling;
        //        float eta_t = Math.Min(etp_t, etp_t * state.sm[t] * lpFcInv);

        //        // Water balance components
        //        qc[t] = P.CFLUX * (1 - smFcRatio);
        //        qf[t] = Math.Min(P.KF * MathF.Pow(state.sw[t], onePlusAlfa), state.sw[t]);
        //        qs[t] = P.KS * state.gw[t];
        //        qr[t] = (qs[t] + qf[t]) * areaScaling;

        //        // Update states
        //        state.sp[t + 1] = state.sp[t] + ss_t + sr_t - sm_t;
        //        state.mw[t + 1] = state.mw[t] + r[t] - sr_t + sm_t - in_t;

        //        float swNext = state.sw[t] + Math.Max(qd[t] + qin[t] - P.PERC, 0) - 
        //                     P.KF * MathF.Pow(state.sw[t], onePlusAlfa) - Math.Min(state.sw[t], qc[t]);
        //        state.sw[t + 1] = Math.Max(swNext, 0);

        //        if (state.sw[t + 1] == 0)
        //        {
        //            qc[t] = Math.Max(state.sw[t] + Math.Max(qd[t] + qin[t] - P.PERC, 0) - 
        //                    P.KF * MathF.Pow(state.sw[t], onePlusAlfa), 0);
        //        }
        //        else
        //        {
        //            qc[t] = Math.Min(state.sw[t], qc[t]);
        //        }

        //        state.sm[t + 1] = Math.Max(state.sm[t] + in_t - qd[t] - qin[t] + qc[t] - eta_t, 0);
        //        state.gw[t + 1] = (1 - P.KS) * state.gw[t] + Math.Min(qd[t] + qin[t], P.PERC);
        //    }

        //    return new ElevationBandData
        //    {
        //        r = r,
        //        ssp = state.sp,
        //        ssw = state.sw,
        //        ssm = state.sm,
        //        sgw = state.gw,
        //        qd = qd,
        //        qc = qc,
        //        qin = qin,
        //        qf = qf,
        //        qs = qs,
        //        qr = qr
        //    };
        //}


         }




    public class HBVParams
    {
        public float FC;      // max soil moisture storage, field capacity
        public float BETA;    // shape coefficient governing fate of water input to soil moisture storage
        public float LP;      // threshold for reduction of evaporation
        public float ALFA;    // non-linearity parameter
        public float KF;      // fast run-off parameter
        public float KS;      // slow run-off parameter
        public float PERC;    // percolation, max flow from upper to lower storage
        public float CFLUX;   // rate of capillary rise
        public float TT;       // threshold temperature for snowfall
        public float TTI;      // threshold temperature for snow melt
        public float CFMAX;    // degree-day factor
        public float FOCFMAX;  // degree day factor, rate of snowmelt
        public float CFR;      // refreezing factor
        public float WHC;      // water holding capacity


        public float MAXBAS;
        public float precip_corr;

        public float TTG;      // threshold temperature for ice in ground
        public float TTGI;     // threshold temperature variation  
        public float corrTg;   // correction for ground temperature
        public float infil_stop; // threshold temperature for infiltration stop

        public HBVParams(float fC, float bETA, float lP, float aLFA, float kF, float kS, float pERC, float cFLUX, float tT, float tTI, float cFMAX, float fOCFMAX, float cFR, float wHC, float mAXBAS, float precip_corr, float tTG, float tTGI, float corrTg, float infil_stop)
        {
            FC = fC;
            BETA = bETA;
            LP = lP;
            ALFA = aLFA;
            KF = kF;
            KS = kS;
            PERC = pERC;
            CFLUX = cFLUX;
            TT = tT;
            TTI = tTI;
            CFMAX = cFMAX;
            FOCFMAX = fOCFMAX;
            CFR = cFR;
            WHC = wHC;
            MAXBAS = mAXBAS;
            this.precip_corr = precip_corr;
            TTG = tTG;
            TTGI = tTGI;
            this.corrTg = corrTg;
            this.infil_stop = infil_stop;
        }
    }

    //public enum ModelType
    //{
    //    HBV = 0,
    //    HBVPermafrost = 1
    //}



    public struct RandomSlim
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static uint WangHash(uint s) { s = (s ^ 61) ^ (s >> 16); s *= 9; s = s ^ (s >> 4); s *= 0x27d4eb2d; s = s ^ (s >> 15); return s; }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static uint RandomInt(ref uint s) { s ^= s << 13; s ^= s >> 17; s ^= s << 5; return s; }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static float RandomFloat(ref uint s) { return RandomInt(ref s) * 2.3283064365387e-10f; }

    };


    public class PSO
    {
        private readonly int numParticles;
        private readonly int dimensions;
        private readonly float[] positions;        // [particle_idx * dimensions + dim_idx]
        private readonly float[] velocities;       // [particle_idx * dimensions + dim_idx] 
        private readonly float[] personalBest;     // [particle_idx * dimensions + dim_idx]
        private readonly float[] personalBestVal;  // [particle_idx]
        private readonly float[] globalBest;       // [dimensions]

        private readonly float[] searchSpaceLowerLimit;
        private readonly float[] searchSpaceUpperLimit;


        private float globalBestVal;
        private readonly Func<float[], float> objectiveFunction;


        private const float W = 0.729f;      // Inertia weight
        private const float C1 = 1.49445f;   // Cognitive parameter
        private const float C2 = 1.49445f;   // Social parameter

        public PSO(int dimensions, Func<float[], float> objectiveFunction, float[] searchSpaceLowerLimit, float[] searchSpaceUpperLimit, int particleFactor = 10)
        {
            if(searchSpaceLowerLimit.Length != dimensions || searchSpaceUpperLimit.Length != dimensions)
            {
                throw new ArgumentException("Search space limits must match number of dimensions");
            }

            this.numParticles = dimensions * particleFactor;
            this.dimensions = dimensions;
            this.objectiveFunction = objectiveFunction;

            positions = new float[numParticles * dimensions];
            velocities = new float[numParticles * dimensions];
            personalBest = new float[numParticles * dimensions];
            personalBestVal = new float[numParticles];
            globalBest = new float[dimensions];

            this.searchSpaceLowerLimit = searchSpaceLowerLimit;
            this.searchSpaceUpperLimit = searchSpaceUpperLimit;

            // Initialize with worst possible value
            globalBestVal = float.MaxValue;

            Initialize();
        }

        public void Initialize()
        {
            //for (int i = 0; i < numParticles; i++)
            Parallel.For(0, numParticles, i =>
            {
                uint seed = (uint)(0 * numParticles * dimensions + i * dimensions + 1);

                for (int d = 0; d < dimensions; d++)
                {
                    int idx = i * dimensions + d;
                    positions[idx] = (RandomSlim.RandomFloat(ref seed) * (searchSpaceUpperLimit[d] - searchSpaceLowerLimit[d]) + searchSpaceLowerLimit[d]);
                    velocities[idx] = (RandomSlim.RandomFloat(ref seed) * (searchSpaceUpperLimit[d] - searchSpaceLowerLimit[d]) * 0.1f);
                    personalBest[idx] = positions[idx];
                }

                float[] particlePosition = new float[dimensions];
                Array.Copy(positions, i * dimensions, particlePosition, 0, dimensions);
                personalBestVal[i] = objectiveFunction(particlePosition);

                if (float.IsNaN(personalBestVal[i]))
                    throw new Exception();

                if (personalBestVal[i] < globalBestVal)
                {
                    globalBestVal = personalBestVal[i];
                    Array.Copy(particlePosition, globalBest, dimensions);
                }
            });
        }

        public void IterateNTimes(int n = 10000)
        {
            for (int i = 0; i < n; i++)
            {
                Iterate(i);

                Console.WriteLine($"PSO Iteration {i}");
            }
        }

        public void Iterate(int iterationIndex)
        {
            //for (int i = 0; i < numParticles; i++)
            Parallel.For(0, numParticles, i =>
            {
                uint seed = (uint)(iterationIndex * numParticles * dimensions + i * dimensions + 1);

                for (int d = 0; d < dimensions; d++)
                {
                    int idx = i * dimensions + d;
                    float r1 = RandomSlim.RandomFloat(ref seed);
                    float r2 = RandomSlim.RandomFloat(ref seed);

                    velocities[idx] = W * velocities[idx] +
                                    C1 * r1 * (personalBest[idx] - positions[idx]) +
                                    C2 * r2 * (globalBest[d] - positions[idx]);

                    positions[idx] += velocities[idx];

                    positions[idx] = Clamp(positions[idx], searchSpaceLowerLimit[d], searchSpaceUpperLimit[d]);
                }

                float[] particlePosition = new float[dimensions];
                Array.Copy(positions, i * dimensions, particlePosition, 0, dimensions);
                float value = objectiveFunction(particlePosition);

                if (value < personalBestVal[i])
                {
                    personalBestVal[i] = value;
                    Array.Copy(positions, i * dimensions, personalBest, i * dimensions, dimensions);

                    if (value < globalBestVal)
                    {
                        lock (globalBest)
                        {
                            if (value < globalBestVal)
                            {
                                globalBestVal = value;
                                Array.Copy(particlePosition, globalBest, dimensions);
                            }
                        }
                    }
                }
            });
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

        public float[] GetGlobalBest() => globalBest;
        public float GetGlobalBestValue() => globalBestVal;
    }


    public class PSODriver
    {
        public static void RunOptimizer()
        {
            // Load meteorological data
            CsvDataMeteoData meteoDataMeasured = new CsvDataMeteoData(@"C:\Users\twidmer\Downloads\N\meteodata_24h.csv");

            // Load catchment info
            var catchmentInfo = new CsvDataElevationBands(@"C:\Users\twidmer\Downloads\N\Elevation_bands.csv");

            // Set timestep in seconds
            float seconds = 86400f; // 24h
            // Other options:
            // float seconds = 43200f; // 12h
            // float seconds = 21600f; // 6h  
            // float seconds = 10800f; // 3h
            // float seconds = 3600f;  // 1h

            // Create arrays from the needed data
            var dTime = meteoDataMeasured.allDateTimes;
            var tmean = meteoDataMeasured.Tmean;
            var pet = meteoDataMeasured.PET;
            var precip = meteoDataMeasured.Precipitation_mm_even;
            var discharge = meteoDataMeasured.Discharge.Select(x => x / 1000f).ToArray(); // Convert to m3/s
            var tground = meteoDataMeasured.Tground;

            int numSubcatchments = 8;
            int numParameters = 20;

            int numElevationBands = catchmentInfo.ElevationBottom.Count;


            HBVParams[] pars_test = new HBVParams[numElevationBands];
            for (int i = 0; i < numElevationBands; i++)
            {
                pars_test[i] = new HBVParams(500, 3, 0.5f, 1, 0.1f, 0.01f, 3, 2, 0, 2, 5, 0.5f, 0.2f, 0.5f, 1, 5, 0, 1, 1, 0);
            }

            int perma = 1;
            int timestep = 24;
            int period = 2014;

            var global = new Global
            {
                DTime = seconds,
                Tmean = tmean.ToArray(),
                PET = pet.ToArray(),
                Precip = precip.ToArray(),
                Discharge = discharge,
                Catchmentinfo = new CatchmentInfo
                {
                    Area_m2 = catchmentInfo.AreaM2.ToArray(),
                    rel_Area = catchmentInfo.RelArea.ToArray(),
                    Elevation_bottom = catchmentInfo.ElevationBottom.ToArray(),
                    Elevation_top = catchmentInfo.ElevationTop.ToArray(),
                    Total_Area_m2 = catchmentInfo.TotalAreaM2[0]
                }
            };

            float test_KGE = HBVModel.KGE_HBV(pars_test, perma, seconds, global);
            ElevationBandData test_HBV = HBVModel.HBV_permafrost_semidistributed(pars_test, perma, seconds, global);



            HBVParams[] minBound = new HBVParams[numElevationBands]; // new HBVParams(0.1f, 0.01f, 0.10f, 0.10f, 0.0005f, 0.0005f, 0.01f, 0.0f, -3f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 1f, 1f, 1f, 1f);
            HBVParams[] maxBound = new HBVParams[numElevationBands]; // new HBVParams(1000f, 7f, 1f, 3f, 0.3f, 0.3f, 6f, 4f, 4f, 7f, 20f, 1f, 1f, 0.8f, 5f, 100f, 1f, 1f, 1f, 1f);

            for (int i = 0; i < numElevationBands; i++)
            {
                minBound[i] = new HBVParams(0.1f, 0.01f, 0.10f, 0.10f, 0.0005f, 0.0005f, 0.01f, 0.0f, -3f, 0f, 0f, 0f, 0f, 0f, 0.006f, 0f, 1f, 1f, 1f, 1f);
                maxBound[i] = new HBVParams(1000f, 7f, 1f, 3f, 0.3f, 0.3f, 6f, 4f, 4f, 7f, 20f, 1f, 1f, 0.8f, 5f, 100f, 1f, 1f, 1f, 1f);
            }

            Func<float[], float> objFun = delegate (float[] optimParams)
            {

                HBVParams[] pars = new HBVParams[numElevationBands];
                for (int i = 0; i < numElevationBands; ++i)
                {
                    int offset = numParameters * i;
                    pars[i] = new HBVParams(
                        optimParams[offset + 0],  // FC
                        optimParams[offset + 1],  // BETA 
                        optimParams[offset + 2],  // LP
                        optimParams[offset + 3],  // ALFA
                        optimParams[offset + 4],  // KF
                        optimParams[offset + 5],  // KS
                        optimParams[offset + 6],  // PERC
                        optimParams[offset + 7],  // CFLUX
                        optimParams[offset + 8],  // TT
                        optimParams[offset + 9],  // TTI
                        optimParams[offset + 10], // CFMAX
                        optimParams[offset + 11], // FOCFMAX
                        optimParams[offset + 12], // CFR
                        optimParams[offset + 13], // WHC
                        optimParams[offset + 14], // MAXBAS
                        optimParams[offset + 15], // precip_corr
                        optimParams[offset + 16], // TTG
                        optimParams[offset + 17], // TTGI
                        optimParams[offset + 18], // corrTg
                        optimParams[offset + 19]  // infil_stop
                    );
                }

                return HBVModel.KGE_HBV(pars, perma, seconds, global);
            };

            // Helper function to convert HBVParams[] to float[]
            Func<HBVParams[], float[]> paramsToFloat = delegate (HBVParams[] hbvParams)
            {
                float[] result = new float[hbvParams.Length * numParameters];
                for (int i = 0; i < hbvParams.Length; i++)
                {
                    int offset = numParameters * i;
                    result[offset + 0] = hbvParams[i].FC;
                    result[offset + 1] = hbvParams[i].BETA;
                    result[offset + 2] = hbvParams[i].LP;
                    result[offset + 3] = hbvParams[i].ALFA;
                    result[offset + 4] = hbvParams[i].KF;
                    result[offset + 5] = hbvParams[i].KS;
                    result[offset + 6] = hbvParams[i].PERC;
                    result[offset + 7] = hbvParams[i].CFLUX;
                    result[offset + 8] = hbvParams[i].TT;
                    result[offset + 9] = hbvParams[i].TTI;
                    result[offset + 10] = hbvParams[i].CFMAX;
                    result[offset + 11] = hbvParams[i].FOCFMAX;
                    result[offset + 12] = hbvParams[i].CFR;
                    result[offset + 13] = hbvParams[i].WHC;
                    result[offset + 14] = hbvParams[i].MAXBAS;
                    result[offset + 15] = hbvParams[i].precip_corr;
                    result[offset + 16] = hbvParams[i].TTG;
                    result[offset + 17] = hbvParams[i].TTGI;
                    result[offset + 18] = hbvParams[i].corrTg;
                    result[offset + 19] = hbvParams[i].infil_stop;
                }
                return result;
            };


            PSO optimizer = new PSO(numSubcatchments * numParameters, objFun, paramsToFloat(minBound), paramsToFloat(maxBound));

            optimizer.IterateNTimes();



            // Alternative bounds:
            // Min: [0.1, 0.01, 0.10, 0.10, 0.0005, 0.0005, 0.01, 0.0, -3, 0, 0, 0, 0, 0, 0, 0, -5, -2, 0, 1]
            // Max: [1000, 7, 1, 3, 0.3, 0.3, 6, 4, 4, 7, 20, 1, 1, 0.8, 5, 100, 0, 2, 5, 1]

            // Alternative bounds 2:
            // Min: [0.1, 0.01, 0.10, 0.10, 0.0005, 0.0005, 0.01, 0.0, -3, 0, 0, 0, 0, 0, 0, 0, -5, -2, 0, -5]
            // Max: [1000, 7, 1, 3, 0.3, 0.3, 6, 4, 4, 7, 20, 1, 1, 0.8, 5, 100, 0, 2, 5, 2]




        }

        public static (ElevationBandData, CsvDataMeteoData) Run()
        {
            // Load meteorological data
            CsvDataMeteoData meteoData = new CsvDataMeteoData(@"C:\Users\twidmer\Downloads\N\meteodata_24h.csv");

            // Load catchment info
            var catchmentInfo = new CsvDataElevationBands(@"C:\Users\twidmer\Downloads\N\Elevation_bands_lumped.csv");

            // Set timestep in seconds
            float seconds = 86400f; // 24h
            // Other options:
            // float seconds = 43200f; // 12h
            // float seconds = 21600f; // 6h  
            // float seconds = 10800f; // 3h
            // float seconds = 3600f;  // 1h

            // Create arrays from the needed data
            var dTime = meteoData.allDateTimes;
            var tmean = meteoData.Tmean;
            var pet = meteoData.PET;
            var precip = meteoData.Precipitation_mm_even;
            var discharge = meteoData.Discharge.Select(x => x / 1000f).ToArray(); // Convert to m3/s
            var tground = meteoData.Tground;

            int numSubcatchments = 8;
            int numParameters = 20;

            int numElevationBands = catchmentInfo.ElevationBottom.Count;


            HBVParams[] pars_test = new HBVParams[numElevationBands];
            for (int i = 0; i < numElevationBands; i++)
            {
                pars_test[i] = new HBVParams(500, 3, 0.5f, 1, 0.1f, 0.01f, 3, 2, 0, 2, 5, 0.5f, 0.2f, 0.5f, 1, 5, 0, 1, 1, 0);
            }



            // Run the HBV model for each elevation band and combine the discharge at the end

            // Parameter bounds
            // Min bounds: [0.1, 0.01, 0.10, 0.10, 0.0005, 0.0005, 0.01, 0.0, -3, 0, 0, 0, 0, 0, 0, 0, -5, -2, 0, -5]
            // Max bounds: [1000, 7, 1, 3, 0.3, 0.3, 6, 4, 4, 7, 20, 1, 1, 0.8, 5, 100, 0, 2, 5, 2]

            int perma = 1;
            int timestep = 24;
            int period = 2014;

            var global = new Global
            {
                DTime = seconds,
                Tmean = tmean.ToArray(),
                PET = pet.ToArray(),
                Precip = precip.ToArray(),
                Discharge = discharge,
                Catchmentinfo = new CatchmentInfo
                {
                    Area_m2 = catchmentInfo.AreaM2.ToArray(),
                    rel_Area = catchmentInfo.RelArea.ToArray(),
                    Elevation_bottom = catchmentInfo.ElevationBottom.ToArray(),
                    Elevation_top = catchmentInfo.ElevationTop.ToArray(),
                    Total_Area_m2 = catchmentInfo.TotalAreaM2[0]
                }
            };

            float test_KGE = HBVModel.KGE_HBV(pars_test, perma, seconds, global);
            ElevationBandData test_HBV = HBVModel.HBV_permafrost_semidistributed(pars_test, perma, seconds, global);


            return (test_HBV, meteoData);

            // Run HBV model

            //Func<float[], float> objFun = delegate (float[] optimParams)
            //{
            //    HBVModel.KGE_HBV();
            //};


        }

        public void EvaluateKGE_HBV()
        {

        }

    }
}
