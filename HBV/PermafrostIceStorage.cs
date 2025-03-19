using ScottPlot.Drawing.Colormaps;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HBV
{
    public class PermafrostIceStorage
    {
        public static void HBV_permafrost_ice_storage(float[] Tground, float[] tm, float[] po, float[] etpo, HBVParams P, float A, float seconds, ElevationBandData result)
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

            var O_TTG = P.TTG;
            var O_TTGI = P.TTGI;
            var corrTg = P.corrTg;

            //};

            //var state = new HBVState
            //{
            var state_sm1 = 50f;
            var state_sw1 = 1f;
            var state_gw1 = 40f;
            var state_sp1 = 50f;
            var state_mw1 = 0f;

            var i_mw1 = 10.0f;
            var i_sw1 = 50.0f;
            var i_gw1 = 50.0f;

            var state_sp = result.ssp; // new float[dtot + 1];
                                       //var state_mw = result.ssp; // new float[dtot + 1];
            var state_sm = result.ssm; // new float[dtot + 1];
            var state_sw = result.ssw; // new float[dtot + 1];
            var state_gw = result.sgw; // new float[dtot + 1]
                                       // }; //

            // Only keep arrays needed for output
            float[] r = result.r; // new float[dtot];
            float[] tg = result.tg;
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

            float[] i_mw = new float[dtot + 1];
            float[] i_sw = new float[dtot + 1];
            float[] i_gw = new float[dtot + 1];
            //i_mw[0] = i_mw1;
            //i_sw[0] = i_sw1;
            //i_gw[0] = i_gw1;
            i_mw[1] = i_mw1;
            i_sw[1] = i_sw1;
            i_gw[1] = i_gw1;

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
                tg[t] = corrTg * Tground[t];

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

                // Ice storage balance (mm)
                if (Tground[t] < (O_TTG - O_TTGI / 2)) // freezing of water in soil
                {
                    i_mw[t + 1] =   i_mw[t] + (O_CFR * O_CFMAX * mw_t); // freezing of water in mw
                    i_sw[t + 1] =   i_sw[t] + (O_CFR * O_CFMAX * state_sw[t]); // freezing of water in sw
                    i_gw[t + 1] = i_gw[t] + (O_CFR * O_CFMAX * state_gw[t]); // freezing of water in gw
                }
                else // melting of ice in soil
                {
                    i_mw[t + 1] =   Math.Max(i_mw[t] - (O_CFMAX * i_mw1), 0); // melting of ice in mw
                    i_sw[t + 1] =   Math.Max(i_sw[t] - (O_CFMAX * i_sw1), 0); // melting of ice in sw
                    i_gw[t + 1] = Math.Max(i_gw[t] - (O_CFMAX * i_gw1), 0); // melting of ice in gw
                }

                // Meltwater in snowpack (meltwater + rain - refreeze + snowmelt - infiltration)
                mw_t = mw_t + r[t] - sr_t + sm_t - in_t - (i_mw[t + 1] - i_mw[t]);

                // Surface water balance (mm)
                state_sw[t + 1] = Math.Max(state_sw[t] + Math.Max((qd[t] + qin[t] - P.PERC), 0) -
                    P.KF * MathF.Pow(state_sw[t], (1 + P.ALFA)) - Math.Min(state_sw[t], qc[t]), 0);
                state_sw[t + 1] = Math.Max(state_sw[t + 1] - (i_sw[t + 1] - i_sw[t]), 0);

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
                state_gw[t + 1] = (1 - P.KS) * state_gw[t] + Math.Min((qd[t] + qin[t]), P.PERC) - (i_gw[t + 1] - i_gw[t]);
                state_gw[t + 1] = Math.Max(state_gw[t + 1], 0);
            }
        }
    }
}
