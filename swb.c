/* Functions for the I cells - WB model.                                  */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "nr.h"
#include "sa.h"
#include "sn.h"
#include "swb.h"

/* This function substitutes the functions for the GA model */
void substiture_function_WB(func_cell_model_I *fcmI, fl_st fl)
{
  fcmI->read_cell_par = &read_cell_par_WB;
  fcmI->write_cell_par = &write_cell_par_WB;
  fcmI->steady_state_var = &steady_state_var_WB;
  fcmI->update_cell = &update_cell_WB;
}

/* This function reads the data for the I cells */
void read_cell_par_WB(cell_par *cpar, fl_st fl)
{
  fscanf(fl.tmp, "gNa=%lf gKdr=%lf gL=%lf\n", &cpar->gNa, &cpar->gKdr, 
  &cpar->gL);

  I_app_read(cpar, fl);

  fscanf(fl.tmp, "Cm=%lf VNa=%lf VK=%lf VL=%lf phi=%lf\n", &cpar->Cm,
  &cpar->VNa, &cpar->VK, &cpar->VL, &cpar->phi);

  fscanf(fl.tmp, "non=%d rho=%lf Vinc1=%lf Vinc2=%lf\n", &cpar->non,
  &cpar->rho, &cpar->Vinc1, &cpar->Vinc2);

  /* Defining the structure of the cell's compartments and variables */
  cpar->nceq = 3;    /* number of cellular variables: V, h, n   */
  cpar->nseq = 2;    /* number of synaptic variables: s         */
  cpar->neq = cpar->nceq + cpar->nseq;
  fprintf(fl.out, "I: nceq=%d nseq=%d neq=%d\n", cpar->nceq, cpar->nseq,
  cpar->neq);
  fflush(fl.out);
}

/* This function writes the data for for the I cells */
void write_cell_par_WB(cell_par *cpar, run_par runpar, fl_st fl)
{
  int itIa;

  fprintf(fl.out, "I CELL: WB model\n");

  fprintf(fl.out, "gNa=%lf gKdr=%lf gL=%lf\n", cpar->gNa, cpar->gKdr,
  cpar->gL);

  I_app_write(cpar, runpar.nt, fl);
  
  fprintf(fl.out, "Cm=%lf VNa=%lf VK=%lf VL=%lf phi=%lf\n", cpar->Cm,
  cpar->VNa, cpar->VK, cpar->VL, cpar->phi);

  fprintf(fl.out, "non=%d rho=%lf Vinc1=%lf Vinc2=%lf\n", cpar->non,
  cpar->rho, cpar->Vinc1, cpar->Vinc2);
  fflush(fl.out);
}

/* This function calculates the steady-state variables for a specific V */
void steady_state_var_WB(net_par *netpar, double *Varb, fl_st fl)
{
  double alphahs, betahs, Hinfs, tauHs;
  double alphans, betans, Ninfs, tauNs;
  double Vc, s_inf;
  int nceq;

  nceq = netpar->I.nceq;

  Vc = Varb[1];                                                     /* V */

  alphahs = 0.07*exp(-(Vc+58.0)/20.0);
  betahs  = 1.0/(1.0+exp(-(Vc+28.0)/10.0));
  Varb[2] = alphahs/(alphahs+betahs);                               /* h */

  alphans = 0.01*(Vc+34.0)/(1.0-exp(-(Vc+34.0)/10.0));
  betans  = 0.125*exp(-(Vc+44.0)/80.0);
  Varb[3] = alphans/(alphans+betans);                               /* n */

  Varb[nceq+1] = 1.0;                                               /* tA */

  s_inf   = Gammaf(Vc, netpar->A.ths, netpar->A.sigs);              /* sA */
  Varb[nceq+2] = netpar->A.kf * s_inf / (netpar->A.kf * s_inf + netpar->A.kr); 
}

/* This function updates the variables of an inhibitory cell   */
void update_cell_WB(cell_par *parI, syn_GABAA_par *synA, double Iapp,
     double *Varc, double *kout, int ion, int it, fl_st fl)
{
  double Vc, hc, nc, tc, sc, s_inf;
  double alphams, betams, Minfs, alphahs, betahs, Hinfs, tauHs;
  double alphans, betans, Ninfs, tauNs;
  double INa, IKdr, ncsq;

  Vc = Varc[1];
  hc = Varc[2];
  nc = Varc[3];
  tc = Varc[parI->nceq+1];
  sc = Varc[parI->nceq+2];

  alphams = 0.1*(Vc+35.0)/(1.0-exp(-(Vc+35.0)/10.0));
  betams  = 4.0*exp(-(Vc+60.0)/18.0);
  Minfs   = alphams/(alphams+betams);
  
  alphahs = 0.07*exp(-(Vc+58.0)/20.0);
  betahs  = 1.0/(1.0+exp(-(Vc+28.0)/10.0));
  Hinfs   = alphahs/(alphahs+betahs);
  tauHs   = 1.0/(alphahs+betahs);

  alphans = 0.01*(Vc+34.0)/(1.0-exp(-(Vc+34.0)/10.0));
  betans  = 0.125*exp(-(Vc+44.0)/80.0);
  Ninfs   = alphans/(alphans+betans);
  tauNs   = 1.0/(alphans+betans);

  INa  = parI->gNa * Minfs * Minfs * Minfs * hc * (Vc - parI->VNa);
  ncsq = nc*nc;
  IKdr = parI->gKdr * ncsq * ncsq * (Vc - parI->VK);

  kout[1] = (-parI->gL * (Vc - parI->VL) - INa - IKdr + Iapp)     /* V */
            / parI->Cm;        

  kout[2] = parI->phi * (Hinfs - hc) / tauHs;                     /* h */

  kout[3] = parI->phi * (Ninfs - nc) / tauNs;                     /* n */

  s_inf   = Gammaf(Vc, synA->ths, synA->sigs);                    /* t */
  kout[parI->nceq+1] = -synA->kt * s_inf * tc + synA->kv * (1.0 - tc);

                                                                  /* s */
  kout[parI->nceq+2] = synA->kf * s_inf * tc * (1.0 - sc) - synA->kr * sc;  
} 
