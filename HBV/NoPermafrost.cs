using HBV;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HBV
{
    public static class NoPermafrost
    {
        public static void HBV_no_permafrost_optimized(float[] tm, float[] po, float[] etpo, 
            HBVParams P, float A, float seconds, ElevationBandData result)
        {
            float FFO = 0.0f;  // percent of forest
            float FFI = 1.00f; // not forested
            int dtot = etpo.Length; // number of timesteps

            //var O = new HBVOptions
            //{
                var O_TCON = seconds;
                var O_CEVPFO = 1.15f;
                var O_TT = P.TT;
                var O_TTI = P.TTI;
                var O_DTTM = 0.0f;
                var O_CFMAX = P.CFMAX;
                var O_FOCFMAX = P.FOCFMAX;
                var O_CFR = P.CFR;
                var O_WHC = P.WHC;
            //};

            //var state = new HBVState
            //{
                var state_sm1 = 50f;
                var state_sw1 = 1f;
                var state_gw1 = 40f;
                var state_sp1 = 50f;
                var state_mw1 = 0f;
                
                var state_sp = result.ssp; // new float[dtot + 1];
                //var state_mw = result.ssp; // new float[dtot + 1];
                var state_sm = result.ssm; // new float[dtot + 1];
                var state_sw = result.ssw; // new float[dtot + 1];
                var state_gw = result.sgw; // new float[dtot + 1]
           // }; //

            // Only keep arrays needed for output
            float[] r = result.r; // new float[dtot];
            float[] qd = result.qd; // new float[dtot];
            float[] qc = result.qc; // new float[dtot];
            float[] qin = result.qin; // new float[dtot];
            float[] qf = result.qf; //new float[dtot];
            float[] qs = result.qs; //new float[dtot];
            float[] qr = result.qr; //new float[dtot];

            // Set initial conditions
            state_sp[0] = state_sp1;
            float mw_t = state_mw1;
            state_sm[0] = Math.Min(state_sm1, P.FC);
            state_sw[0] = state_sw1;
            state_gw[0] = state_gw1;
            state_sp[1] = state_sp[0];
            //state_mw[1] = state_mw[0];
            state_sm[1] = state_sm[0];
            state_sw[1] = state_sw[0];
            state_gw[1] = state_gw[0];

            //if (O_TCON == 0.0)
            //    throw new Exception("Division by zero below -> leads to NaN");

            qf[0] = P.KF * MathF.Pow(state_sw[0], 1 + P.ALFA);
            qs[0] = P.KS * state_gw[0];
            qr[0] = (qf[0] + qs[0]) * A / (O_TCON * 1000);

            // Main simulation loop
            for (int t = 0; t < dtot; t++)
            {
                float p_t, ta_t, ss_t, sm_t, sr_t, in_t, etp_t, eta_t;

                // Snow pack balance components
                p_t = po[t];
                ta_t = tm[t];

                if (ta_t < (O_TT - O_TTI * 0.5f))
                {
                    ss_t = p_t * (FFO + FFI);
                    r[t] = 0;
                }
                else if (ta_t >= (O_TT - O_TTI * 0.5f) && ta_t < (O_TT + O_TTI * 0.5f))
                {
                    //if (O_TTI == 0.0)
                    //    throw new Exception("Division by zero below -> leads to NaN");

                    ss_t = p_t * ((O_TT + O_TTI * 0.5f) - ta_t) / O_TTI * (FFO + FFI);
                    r[t] = p_t * (ta_t - (O_TT - O_TTI * 0.5f)) / O_TTI;
                }
                else
                {
                    ss_t = 0;
                    r[t] = p_t;
                }

                // Snow melt, but not more than available snow storage
                sm_t = Math.Min(Math.Max(O_CFMAX * (ta_t - (O_TT + O_DTTM)), 0), state_sp[t]);

                // Refreezing of water in snowpack 
                sr_t = Math.Min(Math.Max(O_CFR * (FFO * O_FOCFMAX + FFI) * O_CFMAX * ((O_TT + O_DTTM) - ta_t), 0), mw_t);

                // Meltwater and rain reservoir
                in_t = Math.Max(mw_t + sm_t + r[t] - sr_t - O_WHC * state_sp[t], 0);

                // Soil moisture balance components (mm/d)
                qd[t] = Math.Max((in_t + state_sm[t] - P.FC), 0); // Direct runoff->overland flow

                //if (P.FC == 0.0)
                //    throw new Exception("Division by zero below -> leads to NaN");

                // Amount of infiltration
                qin[t] = MathF.Pow(state_sm[t] / P.FC, P.BETA) * (in_t - qd[t]);
                qin[t] = Math.Max(qin[t], 0); // Blocking of overinfiltration

                // Potential evaporation
                etp_t = etpo[t] * (FFO * O_CEVPFO + FFI);

                //if (P.LP == 0.0 || P.FC == 0.0)
                //    throw new Exception("Division by zero below -> leads to NaN");

                // Actual evaporation
                eta_t = Math.Min(etp_t, (etp_t * state_sm[t] / (P.LP * P.FC)));

                // Refilling from groundwater
                qc[t] = P.CFLUX * (P.FC - state_sm[t]) / P.FC;

                // Surface water balance components (mm/d)
                qf[t] = Math.Min(P.KF * MathF.Pow(state_sw[t], (1 + P.ALFA)), state_sw[t]);

                // Ground water balance components (mm/d)
                qs[t] = P.KS * state_gw[t];

                //if (O_TCON == 0.0)
                //    throw new Exception("Division by zero below -> leads to NaN");

                // Total discharge (m3/s)
                qr[t] = (qs[t] + qf[t]) * A / (O_TCON * 1000);

                // Snow pack balance (mm)
                state_sp[t + 1] = state_sp[t] + ss_t + sr_t - sm_t;

                // Meltwater in snowpack
                mw_t = mw_t + r[t] - sr_t + sm_t - in_t;

                // Surface water balance (mm)
                state_sw[t + 1] = Math.Max(state_sw[t] + Math.Max((qd[t] + qin[t] - P.PERC), 0) -
                    P.KF * MathF.Pow(state_sw[t], (1 + P.ALFA)) - Math.Min(state_sw[t], qc[t]), 0);

                // Soil moisture balance (mm)
                if (state_sw[t + 1] == 0)
                {
                    qc[t] = state_sw[t] + Math.Max((qd[t] + qin[t] - P.PERC), 0) - P.KF * MathF.Pow(state_sw[t], (1 + P.ALFA));
                    qc[t] = Math.Max(qc[t], 0);
                }
                else
                {
                    qc[t] = Math.Min(state_sw[t], qc[t]);
                }

                state_sm[t + 1] = state_sm[t] + in_t - qd[t] - qin[t] + qc[t] - eta_t;
                state_sm[t + 1] = Math.Max(state_sm[t + 1], 0);

                // Ground water balance (mm)
                state_gw[t + 1] = (1 - P.KS) * state_gw[t] + Math.Min((qd[t] + qin[t]), P.PERC);
            }

            //return new ElevationBandData
            //{
            //    r = r,
            //    ssp = state.sp,
            //    ssw = state.sw,
            //    ssm = state.sm,
            //    sgw = state.gw,
            //    qd = qd,
            //    qc = qc,
            //    qin = qin,
            //    qf = qf,
            //    qs = qs,
            //    qr = qr
            //};
        }



        public static ElevationBandData HBV_no_permafrost(float[] tm, float[] po, float[] etpo, HBVParams P, float A, float seconds)
        {
            float FFO = 0.0f;  // percent of forest
            float FFI = 1.00f; // not forested
            int dtot = etpo.Length; // number of timesteps

            var O = new HBVOptions
            {
                TCON = seconds,
                CEVPFO = 1.15f,
                TT = P.TT,
                TTI = P.TTI,
                DTTM = 0.0f,
                CFMAX = P.CFMAX,
                FOCFMAX = P.FOCFMAX,
                CFR = P.CFR,
                WHC = P.WHC
            };

            var state = new HBVState
            {
                sm1 = 50f,
                sw1 = 1f,
                gw1 = 40f,
                sp1 = 50f,
                mw1 = 0f,

                sp = new float[dtot + 1],
                mw = new float[dtot + 1],
                sm = new float[dtot + 1],
                sw = new float[dtot + 1],
                gw = new float[dtot + 1]
            };

            // Initialize arrays
            float[] p = new float[dtot];
            float[] ta = new float[dtot];
            float[] ss = new float[dtot];
            float[] r = new float[dtot];
            float[] sm = new float[dtot];
            float[] sr = new float[dtot];
            float[] in_ = new float[dtot];
            float[] qd = new float[dtot];
            float[] qin = new float[dtot];
            float[] etp = new float[dtot];
            float[] eta = new float[dtot];
            float[] qc = new float[dtot];
            float[] qf = new float[dtot];
            float[] qs = new float[dtot];
            float[] qt = new float[dtot];
            float[] qr = new float[dtot];

            // Set initial conditions
            state.sp[0] = state.sp1;
            state.mw[0] = state.mw1;
            state.sm[0] = Math.Min(state.sm1, P.FC);
            state.sw[0] = state.sw1;
            state.gw[0] = state.gw1;
            state.sp[1] = state.sp[0];
            state.mw[1] = state.mw[0];
            state.sm[1] = state.sm[0];
            state.sw[1] = state.sw[0];
            state.gw[1] = state.gw[0];

            if (O.TCON == 0.0)
                throw new Exception("Division by zero below -> leads to NaN");

            qf[0] = P.KF * MathF.Pow(state.sw[0], 1 + P.ALFA);
            qs[0] = P.KS * state.gw[0];
            qt[0] = (qf[0] + qs[0]) * A / (O.TCON * 1000);

            // Main simulation loop
            for (int t = 0; t < dtot; t++)
            {
                // Snow pack balance components
                p[t] = po[t];
                ta[t] = tm[t];

                if (ta[t] < (O.TT - O.TTI / 2))
                {
                    ss[t] = p[t] * (FFO + FFI);
                    r[t] = 0;
                }
                else if (ta[t] >= (O.TT - O.TTI / 2) && ta[t] < (O.TT + O.TTI / 2))
                {
                    if (O.TTI == 0.0)
                        throw new Exception("Division by zero below -> leads to NaN");

                    ss[t] = p[t] * ((O.TT + O.TTI / 2) - ta[t]) / O.TTI * (FFO + FFI);
                    r[t] = p[t] * (ta[t] - (O.TT - O.TTI / 2)) / O.TTI;
                }
                else
                {
                    ss[t] = 0;
                    r[t] = p[t];
                }

                // Snow melt, but not more than available snow storage
                sm[t] = Math.Min(Math.Max(O.CFMAX * (ta[t] - (O.TT + O.DTTM)), 0), state.sp[t]);

                // Refreezing of water in snowpack 
                sr[t] = Math.Min(Math.Max(O.CFR * (FFO * O.FOCFMAX + FFI) * O.CFMAX * ((O.TT + O.DTTM) - ta[t]), 0), state.mw[t]);

                // Meltwater and rain reservoir
                in_[t] = Math.Max(state.mw[t] + sm[t] + r[t] - sr[t] - O.WHC * state.sp[t], 0);

                // Soil moisture balance components (mm/d)
                qd[t] = Math.Max((in_[t] + state.sm[t] - P.FC), 0); // Direct runoff->overland flow

                if (P.FC == 0.0)
                    throw new Exception("Division by zero below -> leads to NaN");

                // Amount of infiltration
                qin[t] = MathF.Pow(state.sm[t] / P.FC, P.BETA) * (in_[t] - qd[t]);
                qin[t] = Math.Max(qin[t], 0); // Blocking of overinfiltration

                // Potential evaporation
                etp[t] = etpo[t] * (FFO * O.CEVPFO + FFI);

                if (P.LP == 0.0 || P.FC == 0.0)
                    throw new Exception("Division by zero below -> leads to NaN");

                // Actual evaporation
                eta[t] = Math.Min(etp[t], (etp[t] * state.sm[t] / (P.LP * P.FC)));

                // Refilling from groundwater
                qc[t] = P.CFLUX * (P.FC - state.sm[t]) / P.FC;

                // Surface water balance components (mm/d)
                qf[t] = Math.Min(P.KF * MathF.Pow(state.sw[t], (1 + P.ALFA)), state.sw[t]);

                // Ground water balance components (mm/d)
                qs[t] = P.KS * state.gw[t];

                if (O.TCON == 0.0)
                    throw new Exception("Division by zero below -> leads to NaN");

                // Total discharge (m3/s)
                qt[t] = (qs[t] + qf[t]) * A / (O.TCON * 1000);

                // Snow pack balance (mm)
                state.sp[t + 1] = state.sp[t] + ss[t] + sr[t] - sm[t];

                // Meltwater in snowpack
                state.mw[t + 1] = state.mw[t] + r[t] - sr[t] + sm[t] - in_[t];

                // Surface water balance (mm)
                state.sw[t + 1] = Math.Max(state.sw[t] + Math.Max((qd[t] + qin[t] - P.PERC), 0) -
                    P.KF * MathF.Pow(state.sw[t], (1 + P.ALFA)) - Math.Min(state.sw[t], qc[t]), 0);

                // Soil moisture balance (mm)
                if (state.sw[t + 1] == 0)
                {
                    qc[t] = state.sw[t] + Math.Max((qd[t] + qin[t] - P.PERC), 0) - P.KF * MathF.Pow(state.sw[t], (1 + P.ALFA));
                    qc[t] = Math.Max(qc[t], 0);
                }
                else
                {
                    qc[t] = Math.Min(state.sw[t], qc[t]);
                }

                state.sm[t + 1] = state.sm[t] + in_[t] - qd[t] - qin[t] + qc[t] - eta[t];
                state.sm[t + 1] = Math.Max(state.sm[t + 1], 0);

                // Ground water balance (mm)
                state.gw[t + 1] = (1 - P.KS) * state.gw[t] + Math.Min((qd[t] + qin[t]), P.PERC);

                // River discharge (m3/s)
                qr[t] = qt[t];
            }

            return new ElevationBandData
            {
                r = r,
                ssp = state.sp,
                ssw = state.sw,
                ssm = state.sm,
                sgw = state.gw,
                qd = qd,
                qc = qc,
                qin = qin,
                qf = qf,
                qs = qs,
                qr = qr
            };
        }

    }
}
