/*
 * helicopter.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "helicopter".
 *
 * Model version              : 11.9
 * Simulink Coder version : 9.4 (R2020b) 29-Jul-2020
 * C source code generated on : Fri Mar 31 14:08:43 2023
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "helicopter.h"
#include "helicopter_private.h"
#include "helicopter_dt.h"

/* Block signals (default storage) */
B_helicopter_T helicopter_B;

/* Continuous states */
X_helicopter_T helicopter_X;

/* Block states (default storage) */
DW_helicopter_T helicopter_DW;

/* Real-time model */
static RT_MODEL_helicopter_T helicopter_M_;
RT_MODEL_helicopter_T *const helicopter_M = &helicopter_M_;

/*
 * Writes out MAT-file header.  Returns success or failure.
 * Returns:
 *      0 - success
 *      1 - failure
 */
int_T rt_WriteMat4FileHeader(FILE *fp, int32_T m, int32_T n, const char *name)
{
  typedef enum { ELITTLE_ENDIAN, EBIG_ENDIAN } ByteOrder;

  int16_T one = 1;
  ByteOrder byteOrder = (*((int8_T *)&one)==1) ? ELITTLE_ENDIAN : EBIG_ENDIAN;
  int32_T type = (byteOrder == ELITTLE_ENDIAN) ? 0: 1000;
  int32_T imagf = 0;
  int32_T name_len = (int32_T)strlen(name) + 1;
  if ((fwrite(&type, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&m, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&n, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&imagf, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&name_len, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(name, sizeof(char), name_len, fp) == 0)) {
    return(1);
  } else {
    return(0);
  }
}                                      /* end rt_WriteMat4FileHeader */

/*
 * This function updates continuous states using the ODE1 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE1_IntgData *id = (ODE1_IntgData *)rtsiGetSolverData(si);
  real_T *f0 = id->f[0];
  int_T i;
  int_T nXc = 4;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);
  rtsiSetdX(si, f0);
  helicopter_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; ++i) {
    x[i] += h * f0[i];
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helicopter_output(void)
{
  /* local block i/o variables */
  real_T rtb_Clock;
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T rtb_TmpSignalConversionAtToFile[8];
  real_T Clock_tmp;
  real_T rtb_Backgain;
  real_T rtb_Derivative;
  real_T rtb_Frontgain;
  real_T *lastU;
  int8_T rtAction;
  if (rtmIsMajorTimeStep(helicopter_M)) {
    /* set solver stop time */
    if (!(helicopter_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helicopter_M->solverInfo,
                            ((helicopter_M->Timing.clockTickH0 + 1) *
        helicopter_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helicopter_M->solverInfo,
                            ((helicopter_M->Timing.clockTick0 + 1) *
        helicopter_M->Timing.stepSize0 + helicopter_M->Timing.clockTickH0 *
        helicopter_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helicopter_M)) {
    helicopter_M->Timing.t[0] = rtsiGetT(&helicopter_M->solverInfo);
  }

  /* Reset subsysRan breadcrumbs */
  srClearBC(helicopter_DW.IfActionSubsystem_SubsysRanBC);
  if (rtmIsMajorTimeStep(helicopter_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S7>/HIL Read Encoder Timebase' */

    /* S-Function Block: helicopter/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder(helicopter_DW.HILReadEncoderTimebase_Task,
        1, &helicopter_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          helicopter_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 =
          helicopter_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_HILReadEncoderTimebase_o3 =
          helicopter_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }
  }

  /* Clock: '<Root>/Clock' incorporates:
   *  Clock: '<S6>/Clock'
   *  Derivative: '<S7>/Derivative'
   */
  Clock_tmp = helicopter_M->Timing.t[0];

  /* Clock: '<Root>/Clock' */
  helicopter_B.Clock = Clock_tmp;
  if (rtmIsMajorTimeStep(helicopter_M)) {
    /* ToFile: '<Root>/To File2' */
    {
      if (!(++helicopter_DW.ToFile2_IWORK.Decimation % 1) &&
          (helicopter_DW.ToFile2_IWORK.Count * (1 + 1)) + 1 < 100000000 ) {
        FILE *fp = (FILE *) helicopter_DW.ToFile2_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[1 + 1];
          helicopter_DW.ToFile2_IWORK.Decimation = 0;
          u[0] = helicopter_M->Timing.t[1];
          u[1] = helicopter_B.Clock;
          if (fwrite(u, sizeof(real_T), 1 + 1, fp) != 1 + 1) {
            rtmSetErrorStatus(helicopter_M, "Error writing to MAT-file Time.mat");
            return;
          }

          if (((++helicopter_DW.ToFile2_IWORK.Count) * (1 + 1))+1 >= 100000000)
          {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file Time.mat.\n");
          }
        }
      }
    }
  }

  /* FromWorkspace: '<S8>/From Workspace' */
  {
    real_T *pDataValues = (real_T *) helicopter_DW.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *) helicopter_DW.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex = helicopter_DW.FromWorkspace_IWORK.PrevIndex;
    real_T t = helicopter_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[140]) {
      currTimeIndex = 139;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    helicopter_DW.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          rtb_Clock = pDataValues[currTimeIndex];
        } else {
          rtb_Clock = pDataValues[currTimeIndex + 1];
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_Clock = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 141;
      }
    }
  }

  /* FromWorkspace: '<S8>/From Workspace1' */
  {
    real_T *pDataValues = (real_T *) helicopter_DW.FromWorkspace1_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *) helicopter_DW.FromWorkspace1_PWORK.TimePtr;
    int_T currTimeIndex = helicopter_DW.FromWorkspace1_IWORK.PrevIndex;
    real_T t = helicopter_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[140]) {
      currTimeIndex = 139;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    helicopter_DW.FromWorkspace1_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 4; ++elIdx) {
              (&helicopter_B.FromWorkspace1[0])[elIdx] =
                pDataValues[currTimeIndex];
              pDataValues += 141;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 4; ++elIdx) {
              (&helicopter_B.FromWorkspace1[0])[elIdx] =
                pDataValues[currTimeIndex + 1];
              pDataValues += 141;
            }
          }
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;

        {
          int_T elIdx;
          for (elIdx = 0; elIdx < 4; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&helicopter_B.FromWorkspace1[0])[elIdx] = (real_T) rtInterpolate(d1,
              d2, f1, f2);
            pDataValues += 141;
          }
        }
      }
    }
  }

  if (rtmIsMajorTimeStep(helicopter_M)) {
    /* Gain: '<S7>/Travel: Count to rad' incorporates:
     *  Gain: '<S7>/Travel_gain'
     */
    helicopter_B.TravelCounttorad = helicopter_P.travel_gain *
      rtb_HILReadEncoderTimebase_o1 * helicopter_P.TravelCounttorad_Gain;

    /* Gain: '<S16>/Gain' */
    helicopter_B.Gain = helicopter_P.Gain_Gain * helicopter_B.TravelCounttorad;

    /* Sum: '<Root>/Sum4' incorporates:
     *  Constant: '<Root>/Constant'
     */
    helicopter_B.travel = helicopter_P.Constant_Value + helicopter_B.Gain;

    /* Gain: '<S7>/Pitch: Count to rad' */
    helicopter_B.PitchCounttorad = helicopter_P.PitchCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S13>/Gain' */
    helicopter_B.Gain_i = helicopter_P.Gain_Gain_a *
      helicopter_B.PitchCounttorad;
  }

  /* Gain: '<S17>/Gain' incorporates:
   *  TransferFcn: '<S7>/Travel: Transfer Fcn'
   */
  helicopter_B.Gain_d = (helicopter_P.TravelTransferFcn_C *
    helicopter_X.TravelTransferFcn_CSTATE + helicopter_P.TravelTransferFcn_D *
    helicopter_B.TravelCounttorad) * helicopter_P.Gain_Gain_l;

  /* Gain: '<S14>/Gain' incorporates:
   *  TransferFcn: '<S7>/Pitch: Transfer Fcn'
   */
  helicopter_B.Gain_b = (helicopter_P.PitchTransferFcn_C *
    helicopter_X.PitchTransferFcn_CSTATE + helicopter_P.PitchTransferFcn_D *
    helicopter_B.PitchCounttorad) * helicopter_P.Gain_Gain_ae;

  /* Gain: '<S5>/Gain1' */
  helicopter_B.Gain1[0] = helicopter_P.Gain1_Gain * helicopter_B.travel;
  helicopter_B.Gain1[1] = helicopter_P.Gain1_Gain * helicopter_B.Gain_d;
  helicopter_B.Gain1[2] = helicopter_P.Gain1_Gain * helicopter_B.Gain_i;
  helicopter_B.Gain1[3] = helicopter_P.Gain1_Gain * helicopter_B.Gain_b;

  /* Sum: '<S8>/Sum' */
  helicopter_B.delta_x[0] = helicopter_B.Gain1[0] - helicopter_B.FromWorkspace1
    [0];
  helicopter_B.delta_x[1] = helicopter_B.Gain1[1] - helicopter_B.FromWorkspace1
    [1];
  helicopter_B.delta_x[2] = helicopter_B.Gain1[2] - helicopter_B.FromWorkspace1
    [2];
  helicopter_B.delta_x[3] = helicopter_B.Gain1[3] - helicopter_B.FromWorkspace1
    [3];

  /* Sum: '<S8>/Sum1' incorporates:
   *  Gain: '<S8>/Gain'
   */
  helicopter_B.Sum1 = rtb_Clock - (((helicopter_P.K[0] * helicopter_B.delta_x[0]
    + helicopter_P.K[1] * helicopter_B.delta_x[1]) + helicopter_P.K[2] *
    helicopter_B.delta_x[2]) + helicopter_P.K[3] * helicopter_B.delta_x[3]);
  if (rtmIsMajorTimeStep(helicopter_M)) {
    /* Gain: '<S3>/Gain1' */
    helicopter_B.Gain1_g = helicopter_P.Gain1_Gain_n * helicopter_B.Gain_i;

    /* Constant: '<Root>/elevation_ref' */
    helicopter_B.elev_ref = helicopter_P.elevation_ref_Value;

    /* Gain: '<S7>/Elevation: Count to rad' incorporates:
     *  Gain: '<S7>/Elevation_gain'
     */
    helicopter_B.ElevationCounttorad = helicopter_P.elevation_gain *
      rtb_HILReadEncoderTimebase_o3 * helicopter_P.ElevationCounttorad_Gain;

    /* Gain: '<S11>/Gain' */
    helicopter_B.Gain_e = helicopter_P.Gain_Gain_lv *
      helicopter_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    helicopter_B.elev = helicopter_B.Gain_e +
      helicopter_P.elavation_offsetdeg_Value;

    /* Gain: '<S4>/Gain1' */
    helicopter_B.Gain1_p = helicopter_P.Gain1_Gain_a * helicopter_B.elev;
  }

  /* Gain: '<S12>/Gain' incorporates:
   *  TransferFcn: '<S7>/Elevation: Transfer Fcn'
   */
  helicopter_B.Gain_dg = (helicopter_P.ElevationTransferFcn_C *
    helicopter_X.ElevationTransferFcn_CSTATE +
    helicopter_P.ElevationTransferFcn_D * helicopter_B.ElevationCounttorad) *
    helicopter_P.Gain_Gain_n;

  /* Gain: '<S2>/Gain1' */
  helicopter_B.Gain1_e[0] = helicopter_P.Gain1_Gain_f * helicopter_B.travel;
  helicopter_B.Gain1_e[1] = helicopter_P.Gain1_Gain_f * helicopter_B.Gain_d;
  helicopter_B.Gain1_e[2] = helicopter_P.Gain1_Gain_f * helicopter_B.Gain_i;
  helicopter_B.Gain1_e[3] = helicopter_P.Gain1_Gain_f * helicopter_B.Gain_b;
  helicopter_B.Gain1_e[4] = helicopter_P.Gain1_Gain_f * helicopter_B.elev;
  helicopter_B.Gain1_e[5] = helicopter_P.Gain1_Gain_f * helicopter_B.Gain_dg;
  if (rtmIsMajorTimeStep(helicopter_M)) {
    /* ToFile: '<Root>/To File' */
    {
      if (!(++helicopter_DW.ToFile_IWORK.Decimation % 1) &&
          (helicopter_DW.ToFile_IWORK.Count * (6 + 1)) + 1 < 100000000 ) {
        FILE *fp = (FILE *) helicopter_DW.ToFile_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[6 + 1];
          helicopter_DW.ToFile_IWORK.Decimation = 0;
          u[0] = helicopter_M->Timing.t[1];
          u[1] = helicopter_B.Gain1_e[0];
          u[2] = helicopter_B.Gain1_e[1];
          u[3] = helicopter_B.Gain1_e[2];
          u[4] = helicopter_B.Gain1_e[3];
          u[5] = helicopter_B.Gain1_e[4];
          u[6] = helicopter_B.Gain1_e[5];
          if (fwrite(u, sizeof(real_T), 6 + 1, fp) != 6 + 1) {
            rtmSetErrorStatus(helicopter_M,
                              "Error writing to MAT-file observed_10_3.mat");
            return;
          }

          if (((++helicopter_DW.ToFile_IWORK.Count) * (6 + 1))+1 >= 100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file observed_10_3.mat.\n");
          }
        }
      }
    }
  }

  /* Sum: '<Root>/Sum1' incorporates:
   *  Constant: '<Root>/Vd_bias'
   *  Gain: '<S9>/K_pd'
   *  Gain: '<S9>/K_pp'
   *  Sum: '<S9>/Sum2'
   *  Sum: '<S9>/Sum3'
   */
  rtb_Frontgain = ((helicopter_B.Sum1 - helicopter_B.Gain1_e[2]) *
                   helicopter_P.K_pp - helicopter_P.K_pd * helicopter_B.Gain1_e
                   [3]) + helicopter_P.Vd_ff;

  /* Integrator: '<S6>/Integrator' */
  /* Limited  Integrator  */
  if (helicopter_X.Integrator_CSTATE >= helicopter_P.Integrator_UpperSat) {
    helicopter_X.Integrator_CSTATE = helicopter_P.Integrator_UpperSat;
  } else {
    if (helicopter_X.Integrator_CSTATE <= helicopter_P.Integrator_LowerSat) {
      helicopter_X.Integrator_CSTATE = helicopter_P.Integrator_LowerSat;
    }
  }

  /* Clock: '<S6>/Clock' incorporates:
   *  Integrator: '<S6>/Integrator'
   */
  rtb_Clock = helicopter_X.Integrator_CSTATE;

  /* Sum: '<S6>/Sum' */
  rtb_Derivative = helicopter_B.elev_ref - helicopter_B.Gain1_e[4];

  /* Sum: '<Root>/Sum2' incorporates:
   *  Constant: '<Root>/Vs_bias'
   *  Gain: '<S6>/K_ed'
   *  Gain: '<S6>/K_ep'
   *  Sum: '<S6>/Sum1'
   */
  rtb_Backgain = ((helicopter_P.K_ep * rtb_Derivative + rtb_Clock) -
                  helicopter_P.K_ed * helicopter_B.Gain1_e[5]) +
    helicopter_P.Vs_ff;

  /* Clock: '<S6>/Clock' */
  rtb_Clock = Clock_tmp;

  /* If: '<S6>/If' incorporates:
   *  Gain: '<S6>/K_ei'
   *  Inport: '<S10>/In1'
   */
  if (rtmIsMajorTimeStep(helicopter_M)) {
    rtAction = (int8_T)!(rtb_Clock >= 2.0);
    helicopter_DW.If_ActiveSubsystem = rtAction;
  } else {
    rtAction = helicopter_DW.If_ActiveSubsystem;
  }

  if (rtAction == 0) {
    /* Outputs for IfAction SubSystem: '<S6>/If Action Subsystem' incorporates:
     *  ActionPort: '<S10>/Action Port'
     */
    helicopter_B.In1 = helicopter_P.K_ei * rtb_Derivative;
    if (rtmIsMajorTimeStep(helicopter_M)) {
      srUpdateBC(helicopter_DW.IfActionSubsystem_SubsysRanBC);
    }

    /* End of Outputs for SubSystem: '<S6>/If Action Subsystem' */
  }

  /* End of If: '<S6>/If' */
  if (rtmIsMajorTimeStep(helicopter_M)) {
  }

  /* Derivative: '<S7>/Derivative' */
  if ((helicopter_DW.TimeStampA >= Clock_tmp) && (helicopter_DW.TimeStampB >=
       Clock_tmp)) {
    rtb_Derivative = 0.0;
  } else {
    rtb_Derivative = helicopter_DW.TimeStampA;
    lastU = &helicopter_DW.LastUAtTimeA;
    if (helicopter_DW.TimeStampA < helicopter_DW.TimeStampB) {
      if (helicopter_DW.TimeStampB < Clock_tmp) {
        rtb_Derivative = helicopter_DW.TimeStampB;
        lastU = &helicopter_DW.LastUAtTimeB;
      }
    } else {
      if (helicopter_DW.TimeStampA >= Clock_tmp) {
        rtb_Derivative = helicopter_DW.TimeStampB;
        lastU = &helicopter_DW.LastUAtTimeB;
      }
    }

    rtb_Derivative = (helicopter_B.PitchCounttorad - *lastU) / (Clock_tmp -
      rtb_Derivative);
  }

  /* Gain: '<S15>/Gain' */
  helicopter_B.Gain_l = helicopter_P.Gain_Gain_a1 * rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter_M)) {
  }

  /* Gain: '<S1>/Back gain' incorporates:
   *  Sum: '<S1>/Subtract'
   */
  Clock_tmp = (rtb_Backgain - rtb_Frontgain) * helicopter_P.Backgain_Gain;

  /* Saturate: '<S7>/Back motor: Saturation' */
  if (Clock_tmp > helicopter_P.BackmotorSaturation_UpperSat) {
    /* Saturate: '<S7>/Back motor: Saturation' */
    helicopter_B.BackmotorSaturation = helicopter_P.BackmotorSaturation_UpperSat;
  } else if (Clock_tmp < helicopter_P.BackmotorSaturation_LowerSat) {
    /* Saturate: '<S7>/Back motor: Saturation' */
    helicopter_B.BackmotorSaturation = helicopter_P.BackmotorSaturation_LowerSat;
  } else {
    /* Saturate: '<S7>/Back motor: Saturation' */
    helicopter_B.BackmotorSaturation = Clock_tmp;
  }

  /* End of Saturate: '<S7>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter_M)) {
  }

  /* Gain: '<S1>/Front gain' incorporates:
   *  Sum: '<S1>/Add'
   */
  Clock_tmp = (rtb_Frontgain + rtb_Backgain) * helicopter_P.Frontgain_Gain;

  /* Saturate: '<S7>/Front motor: Saturation' */
  if (Clock_tmp > helicopter_P.FrontmotorSaturation_UpperSat) {
    /* Saturate: '<S7>/Front motor: Saturation' */
    helicopter_B.FrontmotorSaturation =
      helicopter_P.FrontmotorSaturation_UpperSat;
  } else if (Clock_tmp < helicopter_P.FrontmotorSaturation_LowerSat) {
    /* Saturate: '<S7>/Front motor: Saturation' */
    helicopter_B.FrontmotorSaturation =
      helicopter_P.FrontmotorSaturation_LowerSat;
  } else {
    /* Saturate: '<S7>/Front motor: Saturation' */
    helicopter_B.FrontmotorSaturation = Clock_tmp;
  }

  /* End of Saturate: '<S7>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter_M)) {
    /* S-Function (hil_write_analog_block): '<S7>/HIL Write Analog' */

    /* S-Function Block: helicopter/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helicopter_DW.HILWriteAnalog_Buffer[0] = helicopter_B.FrontmotorSaturation;
      helicopter_DW.HILWriteAnalog_Buffer[1] = helicopter_B.BackmotorSaturation;
      result = hil_write_analog(helicopter_DW.HILInitialize_Card,
        helicopter_P.HILWriteAnalog_channels, 2,
        &helicopter_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
      }
    }

    /* SignalConversion generated from: '<S8>/To File' */
    rtb_TmpSignalConversionAtToFile[0] = helicopter_B.FromWorkspace1[0];
    rtb_TmpSignalConversionAtToFile[4] = helicopter_B.Gain1[0];
    rtb_TmpSignalConversionAtToFile[1] = helicopter_B.FromWorkspace1[1];
    rtb_TmpSignalConversionAtToFile[5] = helicopter_B.Gain1[1];
    rtb_TmpSignalConversionAtToFile[2] = helicopter_B.FromWorkspace1[2];
    rtb_TmpSignalConversionAtToFile[6] = helicopter_B.Gain1[2];
    rtb_TmpSignalConversionAtToFile[3] = helicopter_B.FromWorkspace1[3];
    rtb_TmpSignalConversionAtToFile[7] = helicopter_B.Gain1[3];

    /* ToFile: '<S8>/To File' */
    {
      if (!(++helicopter_DW.ToFile_IWORK_j.Decimation % 1) &&
          (helicopter_DW.ToFile_IWORK_j.Count * (8 + 1)) + 1 < 100000000 ) {
        FILE *fp = (FILE *) helicopter_DW.ToFile_PWORK_e.FilePtr;
        if (fp != (NULL)) {
          real_T u[8 + 1];
          helicopter_DW.ToFile_IWORK_j.Decimation = 0;
          u[0] = helicopter_M->Timing.t[1];
          u[1] = rtb_TmpSignalConversionAtToFile[0];
          u[2] = rtb_TmpSignalConversionAtToFile[1];
          u[3] = rtb_TmpSignalConversionAtToFile[2];
          u[4] = rtb_TmpSignalConversionAtToFile[3];
          u[5] = rtb_TmpSignalConversionAtToFile[4];
          u[6] = rtb_TmpSignalConversionAtToFile[5];
          u[7] = rtb_TmpSignalConversionAtToFile[6];
          u[8] = rtb_TmpSignalConversionAtToFile[7];
          if (fwrite(u, sizeof(real_T), 8 + 1, fp) != 8 + 1) {
            rtmSetErrorStatus(helicopter_M,
                              "Error writing to MAT-file q_r_equals1.mat");
            return;
          }

          if (((++helicopter_DW.ToFile_IWORK_j.Count) * (8 + 1))+1 >= 100000000)
          {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file q_r_equals1.mat.\n");
          }
        }
      }
    }
  }
}

/* Model update function */
void helicopter_update(void)
{
  real_T *lastU;

  /* Update for Derivative: '<S7>/Derivative' */
  if (helicopter_DW.TimeStampA == (rtInf)) {
    helicopter_DW.TimeStampA = helicopter_M->Timing.t[0];
    lastU = &helicopter_DW.LastUAtTimeA;
  } else if (helicopter_DW.TimeStampB == (rtInf)) {
    helicopter_DW.TimeStampB = helicopter_M->Timing.t[0];
    lastU = &helicopter_DW.LastUAtTimeB;
  } else if (helicopter_DW.TimeStampA < helicopter_DW.TimeStampB) {
    helicopter_DW.TimeStampA = helicopter_M->Timing.t[0];
    lastU = &helicopter_DW.LastUAtTimeA;
  } else {
    helicopter_DW.TimeStampB = helicopter_M->Timing.t[0];
    lastU = &helicopter_DW.LastUAtTimeB;
  }

  *lastU = helicopter_B.PitchCounttorad;

  /* End of Update for Derivative: '<S7>/Derivative' */
  if (rtmIsMajorTimeStep(helicopter_M)) {
    rt_ertODEUpdateContinuousStates(&helicopter_M->solverInfo);
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++helicopter_M->Timing.clockTick0)) {
    ++helicopter_M->Timing.clockTickH0;
  }

  helicopter_M->Timing.t[0] = rtsiGetSolverStopTime(&helicopter_M->solverInfo);

  {
    /* Update absolute timer for sample time: [0.002s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++helicopter_M->Timing.clockTick1)) {
      ++helicopter_M->Timing.clockTickH1;
    }

    helicopter_M->Timing.t[1] = helicopter_M->Timing.clockTick1 *
      helicopter_M->Timing.stepSize1 + helicopter_M->Timing.clockTickH1 *
      helicopter_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void helicopter_derivatives(void)
{
  XDot_helicopter_T *_rtXdot;
  boolean_T lsat;
  boolean_T usat;
  _rtXdot = ((XDot_helicopter_T *) helicopter_M->derivs);

  /* Derivatives for TransferFcn: '<S7>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter_P.TravelTransferFcn_A *
    helicopter_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S7>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter_P.PitchTransferFcn_A *
    helicopter_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S7>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE += helicopter_P.ElevationTransferFcn_A *
    helicopter_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE += helicopter_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S6>/Integrator' */
  lsat = (helicopter_X.Integrator_CSTATE <= helicopter_P.Integrator_LowerSat);
  usat = (helicopter_X.Integrator_CSTATE >= helicopter_P.Integrator_UpperSat);
  if (((!lsat) && (!usat)) || (lsat && (helicopter_B.In1 > 0.0)) || (usat &&
       (helicopter_B.In1 < 0.0))) {
    _rtXdot->Integrator_CSTATE = helicopter_B.In1;
  } else {
    /* in saturation */
    _rtXdot->Integrator_CSTATE = 0.0;
  }

  /* End of Derivatives for Integrator: '<S6>/Integrator' */
}

/* Model initialize function */
void helicopter_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &helicopter_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options(helicopter_DW.HILInitialize_Card,
      "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(helicopter_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_M, _rt_error_message);
      return;
    }

    if ((helicopter_P.HILInitialize_AIPStart && !is_switching) ||
        (helicopter_P.HILInitialize_AIPEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums = &helicopter_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] = (helicopter_P.HILInitialize_AILow);
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &helicopter_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] = helicopter_P.HILInitialize_AIHigh;
        }
      }

      result = hil_set_analog_input_ranges(helicopter_DW.HILInitialize_Card,
        helicopter_P.HILInitialize_AIChannels, 8U,
        &helicopter_DW.HILInitialize_AIMinimums[0],
        &helicopter_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_P.HILInitialize_AOPStart && !is_switching) ||
        (helicopter_P.HILInitialize_AOPEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums = &helicopter_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] = (helicopter_P.HILInitialize_AOLow);
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &helicopter_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] = helicopter_P.HILInitialize_AOHigh;
        }
      }

      result = hil_set_analog_output_ranges(helicopter_DW.HILInitialize_Card,
        helicopter_P.HILInitialize_AOChannels, 8U,
        &helicopter_DW.HILInitialize_AOMinimums[0],
        &helicopter_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_P.HILInitialize_AOStart && !is_switching) ||
        (helicopter_P.HILInitialize_AOEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter_P.HILInitialize_AOInitial;
        }
      }

      result = hil_write_analog(helicopter_DW.HILInitialize_Card,
        helicopter_P.HILInitialize_AOChannels, 8U,
        &helicopter_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
        return;
      }
    }

    if (helicopter_P.HILInitialize_AOReset) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter_P.HILInitialize_AOWatchdog;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (helicopter_DW.HILInitialize_Card, helicopter_P.HILInitialize_AOChannels,
         8U, &helicopter_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_P.HILInitialize_EIPStart && !is_switching) ||
        (helicopter_P.HILInitialize_EIPEnter && is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &helicopter_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = helicopter_P.HILInitialize_EIQuadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode(helicopter_DW.HILInitialize_Card,
        helicopter_P.HILInitialize_EIChannels, 8U, (t_encoder_quadrature_mode *)
        &helicopter_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_P.HILInitialize_EIStart && !is_switching) ||
        (helicopter_P.HILInitialize_EIEnter && is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &helicopter_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] = helicopter_P.HILInitialize_EIInitial;
        }
      }

      result = hil_set_encoder_counts(helicopter_DW.HILInitialize_Card,
        helicopter_P.HILInitialize_EIChannels, 8U,
        &helicopter_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_P.HILInitialize_POPStart && !is_switching) ||
        (helicopter_P.HILInitialize_POPEnter && is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues = &helicopter_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicopter_P.HILInitialize_POModes;
        }
      }

      result = hil_set_pwm_mode(helicopter_DW.HILInitialize_Card,
        helicopter_P.HILInitialize_POChannels, 8U, (t_pwm_mode *)
        &helicopter_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_POChannels =
          helicopter_P.HILInitialize_POChannels;
        int32_T *dw_POModeValues = &helicopter_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE ||
              dw_POModeValues[i1] == PWM_RAW_MODE) {
            helicopter_DW.HILInitialize_POSortedChans[num_duty_cycle_modes] =
              (p_HILInitialize_POChannels[i1]);
            helicopter_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes] =
              helicopter_P.HILInitialize_POFrequency;
            num_duty_cycle_modes++;
          } else {
            helicopter_DW.HILInitialize_POSortedChans[7U - num_frequency_modes] =
              (p_HILInitialize_POChannels[i1]);
            helicopter_DW.HILInitialize_POSortedFreqs[7U - num_frequency_modes] =
              helicopter_P.HILInitialize_POFrequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(helicopter_DW.HILInitialize_Card,
          &helicopter_DW.HILInitialize_POSortedChans[0], num_duty_cycle_modes,
          &helicopter_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(helicopter_DW.HILInitialize_Card,
          &helicopter_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &helicopter_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues = &helicopter_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicopter_P.HILInitialize_POConfiguration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues = &helicopter_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] = helicopter_P.HILInitialize_POAlignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &helicopter_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] = helicopter_P.HILInitialize_POPolarity;
        }
      }

      result = hil_set_pwm_configuration(helicopter_DW.HILInitialize_Card,
        helicopter_P.HILInitialize_POChannels, 8U,
        (t_pwm_configuration *) &helicopter_DW.HILInitialize_POModeValues[0],
        (t_pwm_alignment *) &helicopter_DW.HILInitialize_POAlignValues[0],
        (t_pwm_polarity *) &helicopter_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs = &helicopter_DW.HILInitialize_POSortedFreqs[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] = helicopter_P.HILInitialize_POLeading;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &helicopter_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter_P.HILInitialize_POTrailing;
        }
      }

      result = hil_set_pwm_deadband(helicopter_DW.HILInitialize_Card,
        helicopter_P.HILInitialize_POChannels, 8U,
        &helicopter_DW.HILInitialize_POSortedFreqs[0],
        &helicopter_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_P.HILInitialize_POStart && !is_switching) ||
        (helicopter_P.HILInitialize_POEnter && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter_P.HILInitialize_POInitial;
        }
      }

      result = hil_write_pwm(helicopter_DW.HILInitialize_Card,
        helicopter_P.HILInitialize_POChannels, 8U,
        &helicopter_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
        return;
      }
    }

    if (helicopter_P.HILInitialize_POReset) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter_P.HILInitialize_POWatchdog;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (helicopter_DW.HILInitialize_Card, helicopter_P.HILInitialize_POChannels,
         8U, &helicopter_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S7>/HIL Read Encoder Timebase' */

  /* S-Function Block: helicopter/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader(helicopter_DW.HILInitialize_Card,
      helicopter_P.HILReadEncoderTimebase_SamplesI,
      helicopter_P.HILReadEncoderTimebase_Channels, 3,
      &helicopter_DW.HILReadEncoderTimebase_Task);
    if (result >= 0) {
      result = hil_task_set_buffer_overflow_mode
        (helicopter_DW.HILReadEncoderTimebase_Task, (t_buffer_overflow_mode)
         (helicopter_P.HILReadEncoderTimebase_Overflow - 1));
    }

    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_M, _rt_error_message);
    }
  }

  /* Start for ToFile: '<Root>/To File2' */
  {
    FILE *fp = (NULL);
    char fileName[509] = "Time.mat";
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helicopter_M, "Error creating .mat file Time.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp, 1 + 1, 0, "Time")) {
      rtmSetErrorStatus(helicopter_M,
                        "Error writing mat file header to file Time.mat");
      return;
    }

    helicopter_DW.ToFile2_IWORK.Count = 0;
    helicopter_DW.ToFile2_IWORK.Decimation = -1;
    helicopter_DW.ToFile2_PWORK.FilePtr = fp;
  }

  /* Start for FromWorkspace: '<S8>/From Workspace' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829859, 0.38860612234601766, 0.10951766463249335,
      -0.11003765145431615, -0.27691054168983442, -0.39790612952526,
      -0.47963605905675011, -0.52359877533084864, -0.52359877559826973,
      -0.52359877559827184, -0.52359877559525625, -0.5034343602334812,
      -0.46496995806905961, -0.42077303003309041, -0.37347825137113466,
      -0.32522954476016258, -0.27772317725680296, -0.23225508114983789,
      -0.18977001056679055, -0.15091060290302505, -0.11606484361633929,
      -0.085410818093897833, -0.058957969249168296, -0.036584363690165767,
      -0.018069704889241978, -0.0031240223552705393, 0.0085878843546138839,
      0.017426054250776479, 0.023758751406675049, 0.027949299291127216,
      0.03034600016893263, 0.031274762590088767, 0.031034060580304,
      0.029891860211338828, 0.028084169832659822, 0.025814896844602364,
      0.023256724330032452, 0.020552753342749086, 0.017818689720132674,
      0.015145386795530458, 0.012601586456742364, 0.010236729994402727,
      0.00808373667752782, 0.0061616717254202591, 0.0044782462003949419,
      0.0030321093234589869, 0.001814908903052137, 0.00081310811936197158,
      9.5570245348675087E-6, -0.00061517496659357906, -0.0010816933623657432,
      -0.0014108831691479473, -0.0016232044676299795, -0.0017381581198929741,
      -0.0017739001476768657, -0.0017469835241546994, -0.0016722069316523758,
      -0.0015625512955923515, -0.0014291864719376823, -0.0012815322244053018,
      -0.0011273594817623467, -0.000972919736394906, -0.00082309227180032973,
      -0.00068154064252001767, -0.00055087144231957819, -0.00043278986339023717,
      -0.00032824785874985274, -0.00023758186734545994, -0.00016063804836563556,
      -9.68838045870557E-5, -4.55050645613575E-5, -5.4893532657196431E-6,
      2.4304874244651664E-5, 4.5091814633080318E-5, 5.8113724249464482E-5,
      6.4606854902171662E-5, 6.5776445238086012E-5, 6.2779129562096081E-5,
      5.6710935114057115E-5, 4.859878450969024E-5, 3.9393060361381238E-5,
      2.995829078344947E-5, 2.1058387470906936E-5, 1.3332252028819269E-5,
      7.25543154300734E-6, 3.0850963395057107E-6, 7.9190416257812757E-7, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helicopter_DW.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter_DW.FromWorkspace_PWORK.DataPtr = (void *) pDataValues0;
    helicopter_DW.FromWorkspace_IWORK.PrevIndex = 0;
  }

  /* Start for FromWorkspace: '<S8>/From Workspace1' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0 } ;

    static real_T pDataValues0[] = { 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1378421413625808, 3.1262155534522376,
      3.10330930000849, 3.0666274151650832, 3.0144539223697704, 2.94565627709518,
      2.85950776327925, 2.7555515879714023, 2.6335051105171137,
      2.4931956060994906, 2.334518576158215, 2.1583782673061158,
      1.9678000585014415, 1.7674109958869935, 1.5625809553885213,
      1.3586838044741918, 1.1605918373780086, 0.97235388559446667,
      0.79694979478239891, 0.63643262362857911, 0.49215933336189333,
      0.36485709441233694, 0.25463979483044608, 0.16109616202899482,
      0.083400574803900243, 0.020423656748079128, -0.029167441590430717,
      -0.066822903969511238, -0.0940390640172371, -0.11230155616148145,
      -0.12304153581961431, -0.12760344279547234, -0.12722271682343744,
      -0.12301190960382118, -0.11595372706812129, -0.10689965404656832,
      -0.096572947115304619, -0.085574921723409178, -0.074393600547682945,
      -0.063413927297430783, -0.052928880644284908, -0.043150944333506851,
      -0.034223500301872239, -0.02623181081720723, -0.019213342946752358,
      -0.013167263934834849, -0.0080629998988050941, -0.0038478030673174597,
      -0.00045331552627555529, 0.0021988489649454289, 0.0041924581198616428,
      0.0056128839341651514, 0.0065441176560935939, 0.0070664922857152168,
      0.0072550250650869064, 0.00717829283458278, 0.0068977561501107907,
      0.0064674529784764755, 0.005933989061613674, 0.00533675915543249,
      0.0047083409026493122, 0.0040750107669852613, 0.0034573389686893179,
      0.0028708275268994549, 0.0023265621896184441, 0.0018318551163349604,
      0.0013908606226841701, 0.0010051510697664814, 0.00067424408850566043,
      0.00039607579366873183, 0.00016741749637961824, -1.5764280054477936E-5,
      -0.00015800288847750221, -0.0002640788176064482, -0.00033881524283117575,
      -0.00038691737178839263, -0.00041285075154368273, -0.00042075359457995727,
      -0.000414378251453702, -0.00039705716697264982, -0.00037168896779048367,
      -0.00034074070922520696, -0.00030626273215289665, -0.000269913023132571,
      -0.00023298841456133491, -0.00019646039111510013, -0.00016101367185642179,
      -0.00012708610401591719, -9.4908725853876536E-5, -6.45451237302267E-5,
      -3.5929413424288044E-5, -8.9023066244445362E-6, 1.6755233402276745E-5,
      4.1291304275298244E-5, 6.4956082728002542E-5, 8.7979585544560287E-5,
      0.0001105555787483951, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, -0.015002048903878, -0.046506351640397819,
      -0.091625013783775022, -0.14672753939691041, -0.20869397116250202,
      -0.27519058111832578, -0.34459405525957371, -0.4158247012369547,
      -0.48818590980886728, -0.56123801767189563, -0.63470811975886954,
      -0.70456123540313487, -0.76231283523884974, -0.80155625041600664,
      -0.81932016204519575, -0.81558860361519658, -0.79236786840187745,
      -0.75295180713247389, -0.70161636323804444, -0.642068684660706,
      -0.57709316104270769, -0.50920895584185244, -0.44086919827306431,
      -0.37417453121468969, -0.31078234884614853, -0.25190767219852694,
      -0.19836439337914524, -0.15062184952143698, -0.10886464019548273,
      -0.073049968558182077, -0.042959918625599378, -0.018247627912883521,
      0.0015229038869867745, 0.016843228881422569, 0.028232730144190418,
      0.036216292085951965, 0.041306827725282651, 0.043992101566076976,
      0.044725284702394126, 0.043918693001197925, 0.041940186612312878,
      0.039111745245561408, 0.035709776128128416, 0.031966757933770365,
      0.028073871481921388, 0.024184316047724877, 0.020417056145219648,
      0.016860787325750293, 0.013577950164931163, 0.010608657965180377,
      0.0079744366184999627, 0.0056817032585330567, 0.0037249348881154419,
      0.0020894985184511071, 0.00075413111769133336, -0.00030692892219717712,
      -0.0011221467374137966, -0.0017212126866203033, -0.0021338556673479929,
      -0.0023889196247665425, -0.00251367301108498, -0.002533320542593945,
      -0.0024706871932980556, -0.0023460457672039273, -0.002177061349284855,
      -0.0019788282931703, -0.0017639779746761834, -0.0015428382115734575,
      -0.0013236279250899308, -0.0011126731791773407, -0.00091463318916472846,
      -0.00073272710574210547, -0.0005689544337469266, -0.00042430371648433384,
      -0.00029894570092755613, -0.00019240851579769408, -0.00010373351897131799,
      -3.1611372151375904E-5, 2.5501372521759298E-5, 6.9284337935007908E-5,
      0.00010147279673138174, 0.00012379303425946403, 0.00013791190828809448,
      0.00014539883608049629, 0.00014769843428317239, 0.00014611209379397023,
      0.00014178687703287154, 0.00013571027134947928, 0.0001287095126496527,
      0.00012145440848885078, 0.00011446284122104207, 0.00010810842720103505,
      0.00010263016010766058, 9.8144283492925032E-5, 9.4659113810727054E-5,
      9.20940112671617E-5, 9.030397281539556E-5, 8.9110874038378314E-5, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.10602875206930784,
      0.2226603793434293, 0.31888147183731452, 0.38944360632143116,
      0.4379550737757536, 0.46997264229815694, 0.49051724876992886,
      0.50343100140899366, 0.51142138585518537, 0.516304398574798,
      0.51925862127101052, 0.49369514260411862, 0.40816625102502213,
      0.27735747107777331, 0.12554854255532832, -0.026373229438175494,
      -0.16411528845645573, -0.27857766738513789, -0.362819311450757,
      -0.4208602502585716, -0.45922218563117462, -0.47977963643710364,
      -0.48299930543731057, -0.47137243454263411, -0.44803173400532159,
      -0.41610372877300006, -0.37842344550628659, -0.33742606603963138,
      -0.29512400752269558, -0.25312442038415944, -0.21266498097863495,
      -0.17465703257617921, -0.13973056785971222, -0.10827820593816484,
      -0.08049664499348852, -0.0564247667243048, -0.035977961704848072,
      -0.01897848995636775, -0.0051818583767153648, 0.0057006821054954582,
      0.013983327509792698, 0.019990343349923223, 0.02404381838027192,
      0.026454223049266035, 0.027513434685478955, 0.027489892292598461,
      0.026625554207334368, 0.0251343497620774, 0.023201839258737977,
      0.020985823207948284, 0.018617670409738718, 0.016204163754237166,
      0.01382969150783564, 0.011558639647707536, 0.0094378667867968113,
      0.0074991671099748913, 0.0057616481708270273, 0.0042339693313325766,
      0.0029164029902768451, 0.00180269463583127, 0.00088170928900110912,
      0.0001388612487848162, -0.00044266859100561273, -0.00088091799708456264,
      -0.0011943173278429198, -0.0014010355302781496, -0.0015184799953429362,
      -0.0015629313842647763, -0.001549294584623051, -0.0014909475746790468,
      -0.0013996710154913217, -0.0012856427222851519, -0.0011574826974115604,
      -0.0010223360237741064, -0.00088598257730188035, -0.00075296413520287153,
      -0.00062672101046423734, -0.00050973178876656533, -0.00040365106676487628,
      -0.00030944127777021713, -0.00022749573342695228, -0.00015775091428027466,
      -9.9786809308777258E-5, -5.2914746211718544E-5, -1.6252681847461758E-5,
      1.1211648811682018E-5, 3.0568980114997757E-5, 4.2947128098269616E-5,
      4.9478688788351377E-5, 5.1276305284897461E-5, 4.941372717570669E-5,
      4.4910571430523127E-5, 3.87182995665869E-5, 3.1704462702397684E-5,
      2.4631848240530552E-5, 1.8129165106239853E-5, 1.2651308116784499E-5,
      8.4323664838459322E-6, 5.4485824269034921E-6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.42411500823370818, 0.466526509059192,
      0.384884369970716, 0.28224853795481525, 0.1940458698442025,
      0.12807027409865138, 0.082178425889304046, 0.051655010567624865,
      0.031961537785506429, 0.019532050874372563, 0.011816890780859821,
      -0.10225391466372982, -0.34211556634416723, -0.52323511982304327,
      -0.60723571408286625, -0.60768708787951142, -0.55096823613402035,
      -0.45784951580900241, -0.33696657616759879, -0.23216375525542965,
      -0.15344774147362727, -0.082229803218585043, -0.012878675993616461,
      0.046507483576612475, 0.093362802144847767, 0.12771202093432249,
      0.15072113305825491, 0.16398951787829988, 0.16920823408626276,
      0.1679983485524339, 0.1618377576222218, 0.15203179361761979,
      0.13970585886237327, 0.12580944768309119, 0.11112624377401992,
      0.0962875130794896, 0.081787220081070366, 0.067997886998175841,
      0.055186526321878725, 0.043530161922794151, 0.033130581613803368,
      0.024028063355149173, 0.016213900124826562, 0.0096416186716963819,
      0.004236846544652389, -9.4169570413849332E-5, -0.0034573523455627537,
      -0.0059648177844745494, -0.0077300420127662112, -0.008864064203450173,
      -0.0094726111895629961, -0.0096540266212576317, -0.0094978889852741378,
      -0.0090842074417527615, -0.0084830914431835071, -0.0077547987070700278,
      -0.0069500757565846335, -0.0061107153581038036, -0.0052702653643232136,
      -0.0044548334177871166, -0.0036839413874207769, -0.0029713921608254955,
      -0.0023261193589782171, -0.0017529976245864206, -0.0012535973234416623,
      -0.000826872809307852, -0.00046977786003194551, -0.00017780555601100751,
      5.4547198236506066E-5, 0.00023338803922496476, 0.00036510623716786404,
      0.00045611317275738865, 0.00051264009964396373, 0.0005405866945664352,
      0.00054541378608049835, 0.00053207376852782058, 0.00050497249888472456,
      0.00046795688682453373, 0.00042432288797252654, 0.00037683915597405443,
      0.00032778217737786466, 0.00027897927658785028, 0.00023185641987378682,
      0.00018748825240816832, 0.00014664825742278172, 0.00010985732264056078,
      7.7429325204311011E-5, 4.9512591927819521E-5, 2.6126242741866694E-5,
      7.19046596709461E-6, -7.4503124280890768E-6, -1.8012622987570466E-5,
      -2.4769087460053166E-5, -2.8055347461253622E-5, -2.8290457844519193E-5,
      -2.6010732537513016E-5, -2.1911427956716194E-5, -1.6875766531471516E-5,
      -1.1935136227784256E-5, -8.0237304746936312E-6, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    } ;

    helicopter_DW.FromWorkspace1_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter_DW.FromWorkspace1_PWORK.DataPtr = (void *) pDataValues0;
    helicopter_DW.FromWorkspace1_IWORK.PrevIndex = 0;
  }

  /* Start for ToFile: '<Root>/To File' */
  {
    FILE *fp = (NULL);
    char fileName[509] = "observed_10_3.mat";
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helicopter_M,
                        "Error creating .mat file observed_10_3.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp, 6 + 1, 0, "observed_10_3")) {
      rtmSetErrorStatus(helicopter_M,
                        "Error writing mat file header to file observed_10_3.mat");
      return;
    }

    helicopter_DW.ToFile_IWORK.Count = 0;
    helicopter_DW.ToFile_IWORK.Decimation = -1;
    helicopter_DW.ToFile_PWORK.FilePtr = fp;
  }

  /* Start for If: '<S6>/If' */
  helicopter_DW.If_ActiveSubsystem = -1;

  /* Start for ToFile: '<S8>/To File' */
  {
    FILE *fp = (NULL);
    char fileName[509] = "q_r_equals1.mat";
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helicopter_M, "Error creating .mat file q_r_equals1.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp, 8 + 1, 0, "ans")) {
      rtmSetErrorStatus(helicopter_M,
                        "Error writing mat file header to file q_r_equals1.mat");
      return;
    }

    helicopter_DW.ToFile_IWORK_j.Count = 0;
    helicopter_DW.ToFile_IWORK_j.Decimation = -1;
    helicopter_DW.ToFile_PWORK_e.FilePtr = fp;
  }

  /* InitializeConditions for TransferFcn: '<S7>/Travel: Transfer Fcn' */
  helicopter_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S7>/Pitch: Transfer Fcn' */
  helicopter_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S7>/Elevation: Transfer Fcn' */
  helicopter_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S6>/Integrator' */
  helicopter_X.Integrator_CSTATE = helicopter_P.Integrator_IC;

  /* InitializeConditions for Derivative: '<S7>/Derivative' */
  helicopter_DW.TimeStampA = (rtInf);
  helicopter_DW.TimeStampB = (rtInf);
}

/* Model terminate function */
void helicopter_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(helicopter_DW.HILInitialize_Card);
    hil_monitor_stop_all(helicopter_DW.HILInitialize_Card);
    is_switching = false;
    if ((helicopter_P.HILInitialize_AOTerminate && !is_switching) ||
        (helicopter_P.HILInitialize_AOExit && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter_P.HILInitialize_AOFinal;
        }
      }

      num_final_analog_outputs = 8U;
    } else {
      num_final_analog_outputs = 0;
    }

    if ((helicopter_P.HILInitialize_POTerminate && !is_switching) ||
        (helicopter_P.HILInitialize_POExit && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter_P.HILInitialize_POFinal;
        }
      }

      num_final_pwm_outputs = 8U;
    } else {
      num_final_pwm_outputs = 0;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(helicopter_DW.HILInitialize_Card
                         , helicopter_P.HILInitialize_AOChannels,
                         num_final_analog_outputs
                         , helicopter_P.HILInitialize_POChannels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &helicopter_DW.HILInitialize_AOVoltages[0]
                         , &helicopter_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog(helicopter_DW.HILInitialize_Card,
            helicopter_P.HILInitialize_AOChannels, num_final_analog_outputs,
            &helicopter_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(helicopter_DW.HILInitialize_Card,
            helicopter_P.HILInitialize_POChannels, num_final_pwm_outputs,
            &helicopter_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(helicopter_DW.HILInitialize_Card);
    hil_monitor_delete_all(helicopter_DW.HILInitialize_Card);
    hil_close(helicopter_DW.HILInitialize_Card);
    helicopter_DW.HILInitialize_Card = NULL;
  }

  /* Terminate for ToFile: '<Root>/To File2' */
  {
    FILE *fp = (FILE *) helicopter_DW.ToFile2_PWORK.FilePtr;
    if (fp != (NULL)) {
      char fileName[509] = "Time.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helicopter_M, "Error closing MAT-file Time.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helicopter_M, "Error reopening MAT-file Time.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 1 + 1, helicopter_DW.ToFile2_IWORK.Count,
           "Time")) {
        rtmSetErrorStatus(helicopter_M,
                          "Error writing header for Time to MAT-file Time.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helicopter_M, "Error closing MAT-file Time.mat");
        return;
      }

      helicopter_DW.ToFile2_PWORK.FilePtr = (NULL);
    }
  }

  /* Terminate for ToFile: '<Root>/To File' */
  {
    FILE *fp = (FILE *) helicopter_DW.ToFile_PWORK.FilePtr;
    if (fp != (NULL)) {
      char fileName[509] = "observed_10_3.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helicopter_M,
                          "Error closing MAT-file observed_10_3.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helicopter_M,
                          "Error reopening MAT-file observed_10_3.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 6 + 1, helicopter_DW.ToFile_IWORK.Count,
           "observed_10_3")) {
        rtmSetErrorStatus(helicopter_M,
                          "Error writing header for observed_10_3 to MAT-file observed_10_3.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helicopter_M,
                          "Error closing MAT-file observed_10_3.mat");
        return;
      }

      helicopter_DW.ToFile_PWORK.FilePtr = (NULL);
    }
  }

  /* Terminate for ToFile: '<S8>/To File' */
  {
    FILE *fp = (FILE *) helicopter_DW.ToFile_PWORK_e.FilePtr;
    if (fp != (NULL)) {
      char fileName[509] = "q_r_equals1.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helicopter_M, "Error closing MAT-file q_r_equals1.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helicopter_M,
                          "Error reopening MAT-file q_r_equals1.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 8 + 1, helicopter_DW.ToFile_IWORK_j.Count,
           "ans")) {
        rtmSetErrorStatus(helicopter_M,
                          "Error writing header for ans to MAT-file q_r_equals1.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helicopter_M, "Error closing MAT-file q_r_equals1.mat");
        return;
      }

      helicopter_DW.ToFile_PWORK_e.FilePtr = (NULL);
    }
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/

/* Solver interface called by GRT_Main */
#ifndef USE_GENERATED_SOLVER

void rt_ODECreateIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEDestroyIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEUpdateContinuousStates(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

#endif

void MdlOutputs(int_T tid)
{
  helicopter_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  helicopter_update();
  UNUSED_PARAMETER(tid);
}

void MdlInitializeSizes(void)
{
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
}

void MdlStart(void)
{
  helicopter_initialize();
}

void MdlTerminate(void)
{
  helicopter_terminate();
}

/* Registration function */
RT_MODEL_helicopter_T *helicopter(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helicopter_P.Integrator_UpperSat = rtInf;
  helicopter_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helicopter_M, 0,
                sizeof(RT_MODEL_helicopter_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helicopter_M->solverInfo,
                          &helicopter_M->Timing.simTimeStep);
    rtsiSetTPtr(&helicopter_M->solverInfo, &rtmGetTPtr(helicopter_M));
    rtsiSetStepSizePtr(&helicopter_M->solverInfo,
                       &helicopter_M->Timing.stepSize0);
    rtsiSetdXPtr(&helicopter_M->solverInfo, &helicopter_M->derivs);
    rtsiSetContStatesPtr(&helicopter_M->solverInfo, (real_T **)
                         &helicopter_M->contStates);
    rtsiSetNumContStatesPtr(&helicopter_M->solverInfo,
      &helicopter_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&helicopter_M->solverInfo,
      &helicopter_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&helicopter_M->solverInfo,
      &helicopter_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&helicopter_M->solverInfo,
      &helicopter_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&helicopter_M->solverInfo, (&rtmGetErrorStatus
      (helicopter_M)));
    rtsiSetRTModelPtr(&helicopter_M->solverInfo, helicopter_M);
  }

  rtsiSetSimTimeStep(&helicopter_M->solverInfo, MAJOR_TIME_STEP);
  helicopter_M->intgData.f[0] = helicopter_M->odeF[0];
  helicopter_M->contStates = ((real_T *) &helicopter_X);
  rtsiSetSolverData(&helicopter_M->solverInfo, (void *)&helicopter_M->intgData);
  rtsiSetSolverName(&helicopter_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helicopter_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helicopter_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helicopter_M->Timing.sampleTimes = (&helicopter_M->Timing.sampleTimesArray[0]);
    helicopter_M->Timing.offsetTimes = (&helicopter_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helicopter_M->Timing.sampleTimes[0] = (0.0);
    helicopter_M->Timing.sampleTimes[1] = (0.002);

    /* task offsets */
    helicopter_M->Timing.offsetTimes[0] = (0.0);
    helicopter_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helicopter_M, &helicopter_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helicopter_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helicopter_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helicopter_M, -1);
  helicopter_M->Timing.stepSize0 = 0.002;
  helicopter_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  helicopter_M->Sizes.checksums[0] = (3184592147U);
  helicopter_M->Sizes.checksums[1] = (2081142745U);
  helicopter_M->Sizes.checksums[2] = (4223715458U);
  helicopter_M->Sizes.checksums[3] = (2053704134U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[2];
    helicopter_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    systemRan[1] = (sysRanDType *)&helicopter_DW.IfActionSubsystem_SubsysRanBC;
    rteiSetModelMappingInfoPtr(helicopter_M->extModeInfo,
      &helicopter_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helicopter_M->extModeInfo, helicopter_M->Sizes.checksums);
    rteiSetTPtr(helicopter_M->extModeInfo, rtmGetTPtr(helicopter_M));
  }

  helicopter_M->solverInfoPtr = (&helicopter_M->solverInfo);
  helicopter_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&helicopter_M->solverInfo, 0.002);
  rtsiSetSolverMode(&helicopter_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helicopter_M->blockIO = ((void *) &helicopter_B);

  {
    int32_T i;
    for (i = 0; i < 6; i++) {
      helicopter_B.Gain1_e[i] = 0.0;
    }

    helicopter_B.Clock = 0.0;
    helicopter_B.FromWorkspace1[0] = 0.0;
    helicopter_B.FromWorkspace1[1] = 0.0;
    helicopter_B.FromWorkspace1[2] = 0.0;
    helicopter_B.FromWorkspace1[3] = 0.0;
    helicopter_B.TravelCounttorad = 0.0;
    helicopter_B.Gain = 0.0;
    helicopter_B.travel = 0.0;
    helicopter_B.Gain_d = 0.0;
    helicopter_B.PitchCounttorad = 0.0;
    helicopter_B.Gain_i = 0.0;
    helicopter_B.Gain_b = 0.0;
    helicopter_B.Gain1[0] = 0.0;
    helicopter_B.Gain1[1] = 0.0;
    helicopter_B.Gain1[2] = 0.0;
    helicopter_B.Gain1[3] = 0.0;
    helicopter_B.delta_x[0] = 0.0;
    helicopter_B.delta_x[1] = 0.0;
    helicopter_B.delta_x[2] = 0.0;
    helicopter_B.delta_x[3] = 0.0;
    helicopter_B.Sum1 = 0.0;
    helicopter_B.Gain1_g = 0.0;
    helicopter_B.elev_ref = 0.0;
    helicopter_B.ElevationCounttorad = 0.0;
    helicopter_B.Gain_e = 0.0;
    helicopter_B.elev = 0.0;
    helicopter_B.Gain1_p = 0.0;
    helicopter_B.Gain_dg = 0.0;
    helicopter_B.Gain_l = 0.0;
    helicopter_B.BackmotorSaturation = 0.0;
    helicopter_B.FrontmotorSaturation = 0.0;
    helicopter_B.In1 = 0.0;
  }

  /* parameters */
  helicopter_M->defaultParam = ((real_T *)&helicopter_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helicopter_X;
    helicopter_M->contStates = (x);
    (void) memset((void *)&helicopter_X, 0,
                  sizeof(X_helicopter_T));
  }

  /* states (dwork) */
  helicopter_M->dwork = ((void *) &helicopter_DW);
  (void) memset((void *)&helicopter_DW, 0,
                sizeof(DW_helicopter_T));

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  helicopter_DW.TimeStampA = 0.0;
  helicopter_DW.LastUAtTimeA = 0.0;
  helicopter_DW.TimeStampB = 0.0;
  helicopter_DW.LastUAtTimeB = 0.0;
  helicopter_DW.HILWriteAnalog_Buffer[0] = 0.0;
  helicopter_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helicopter_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 16;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.BTransTable = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.PTransTable = &rtPTransTable;
  }

  /* Initialize Sizes */
  helicopter_M->Sizes.numContStates = (4);/* Number of continuous states */
  helicopter_M->Sizes.numPeriodicContStates = (0);
                                      /* Number of periodic continuous states */
  helicopter_M->Sizes.numY = (0);      /* Number of model outputs */
  helicopter_M->Sizes.numU = (0);      /* Number of model inputs */
  helicopter_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helicopter_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helicopter_M->Sizes.numBlocks = (77);/* Number of blocks */
  helicopter_M->Sizes.numBlockIO = (24);/* Number of block outputs */
  helicopter_M->Sizes.numBlockPrms = (152);/* Sum of parameter "widths" */
  return helicopter_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
