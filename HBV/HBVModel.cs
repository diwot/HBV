namespace HBV
{
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

            float[] qObs = new float[validIndices.Count];
            float[] qSim = new float[validIndices.Count];
            for (int i = 0; i < validIndices.Count; i++)
            {
                qObs[i] = runoffMeasured[validIndices[i]];
                qSim[i] = runoffModelled[validIndices[i]];
            }

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
                    PetElev = pet_elev,

                    Tground = global.Tground
                };

                ElevationBandData outputHBV = new ElevationBandData(meteoData[i].PetElev.Length); // default;
                switch (permafrostOption)
                {
                    case 1:
                        NoPermafrost.HBV_no_permafrost_optimized(meteoData[i].TairElev, meteoData[i].PrecipElev,
                            meteoData[i].PetElev, pars_elev, Catchmentinfo.Area_m2[i], seconds, outputHBV);

                        //outputHBV = Core.HBV_no_permafrost(meteoData[i].TairElev, meteoData[i].PrecipElev,
                        //    meteoData[i].PetElev, pars_elev, Catchmentinfo.Area_m2[i], seconds);
                        break;
                    case 2:
                        PermafrostIceStorage.HBV_permafrost_ice_storage(meteoData[i].Tground, meteoData[i].TairElev, meteoData[i].PrecipElev,
                            meteoData[i].PetElev, pars_elev, Catchmentinfo.Area_m2[i], seconds, outputHBV);
                        break;
                    //case 3:
                    //    outputHBV = NoPermafrost.HBV_permafrost_infiltration_blocking(meteoData[i].TairElev, meteoData[i].PrecipElev,
                    //        meteoData[i].PetElev, pars_elev, Catchmentinfo.Area_m2[i], seconds);
                    //    break;
                    //default:
                    //    outputHBV = NoPermafrost.HBV_permafrost_combined(meteoData[i].TairElev, meteoData[i].PrecipElev,
                    //        meteoData[i].PetElev, pars_elev, Catchmentinfo.Area_m2[i], seconds);
                    //    break;
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
}
