/* Functions for the E cells - GA model.                                  */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "nr.h"
#include "sa.h"
#include "sn.h"
#include "sgy.h"

/* This function substitutes the functions for the GY model */
void substiture_function_GY(func_cell_model_E *fcmE, fl_st fl)
{
  fcmE->read_cell_par = &read_cell_par_GY;
  fcmE->write_cell_par = &write_cell_par_GY;
  fcmE->steady_state_var = &steady_state_var_GY;
  fcmE->update_cell = &update_cell_GY;
}

/* This function reads the data for the E cells */
void read_cell_par_GY(cell_par *cpar, fl_st fl)
{
  fscanf(fl.tmp, "gNa=%lf gNaP=%lf gKdr=%lf gA=%lf gZ=%lf gL=%lf\n",
  &cpar->gNa, &cpar->gNaP, &cpar->gKdr, &cpar->gA, &cpar->gZ, &cpar->gL);

  I_app_read(cpar, fl);

  fscanf(fl.tmp, "Cm=%lf VNa=%lf VK=%lf VL=%lf phi=%lf\n", &cpar->Cm,
  &cpar->VNa, &cpar->VK, &cpar->VL, &cpar->phi);

  fscanf(fl.tmp, "thetam=%lf sigmam=%lf thetap=%lf sigmap=%lf\n",
  &cpar->thetam, &cpar->sigmam, &cpar->thetap, &cpar->sigmap);
  fscanf(fl.tmp, "thetah=%lf sigmah=%lf thetan=%lf sigman=%lf\n",
  &cpar->thetah, &cpar->sigmah, &cpar->thetan, &cpar->sigman);
  fscanf(fl.tmp, "thetaa=%lf sigmaa=%lf thetab=%lf sigmab=%lf tauB=%lf\n",
  &cpar->thetaa, &cpar->sigmaa, &cpar->thetab, &cpar->sigmab, &cpar->tauB);
  fscanf(fl.tmp, "thetaz=%lf sigmaz=%lf tauZ=%lf\n",
  &cpar->thetaz, &cpar->sigmaz, &cpar->tauZ);

  fscanf(fl.tmp, "non=%d rho=%lf Vinc1=%lf Vinc2=%lf\n", &cpar->non,
  &cpar->rho, &cpar->Vinc1, &cpar->Vinc2);

  /* Defining the structure of the cell's compartments and variables */
  cpar->nceq = 5;    /* number of cellular variables: V, h, n, b, z */
  cpar->nseq = 4;    /* number of synaptic variables: T, sA, sN     */
  cpar->neq = cpar->nceq + cpar->nseq;
  fprintf(fl.out, "E: nceq=%d nseq=%d neq=%d \n", cpar->nceq, cpar->nseq,
  cpar->neq);
  fflush(fl.out);
}

/* This function writes the data for for the E cells */
void write_cell_par_GY(cell_par *cpar, run_par runpar, fl_st fl)
{
  fprintf(fl.out, "E CELL: GY model\n");

  fprintf(fl.out, "gNa=%lf gNaP=%lf gKdr=%lf gA=%lf gZ=%lf gL=%lf\n",
  cpar->gNa, cpar->gNaP, cpar->gKdr, cpar->gA, cpar->gZ, cpar->gL);

  I_app_write(cpar, runpar.nt, fl);
  
  fprintf(fl.out, "Cm=%lf VNa=%lf VK=%lf VL=%lf phi=%lf\n", cpar->Cm,
  cpar->VNa, cpar->VK, cpar->VL, cpar->phi);

  fprintf(fl.out, "thetam=%lf sigmam=%lf thetap=%lf sigmap=%lf\n",
  cpar->thetam, cpar->sigmam, cpar->thetap, cpar->sigmap);
  fprintf(fl.out, "thetah=%lf sigmah=%lf thetan=%lf sigman=%lf\n",
  cpar->thetah, cpar->sigmah, cpar->thetan, cpar->sigman);
  fprintf(fl.out, "thetaa=%lf sigmaa=%lf thetab=%lf sigmab=%lf tauB=%lf\n",
  cpar->thetaa, cpar->sigmaa, cpar->thetab, cpar->sigmab, cpar->tauB);
  fprintf(fl.out, "thetaz=%lf sigmaz=%lf tauZ=%lf\n",
  cpar->thetaz, cpar->sigmaz, cpar->tauZ);

  fprintf(fl.out, "non=%d rho=%lf Vinc1=%lf Vinc2=%lf\n", cpar->non,
  cpar->rho, cpar->Vinc1, cpar->Vinc2);
  fflush(fl.out);
}

/* This function calculates the steady-state variables for a specific V */
void steady_state_var_GY(net_par *netpar, double *Varb, fl_st fl)
{
  double sAMPAinf, sNMDAinf, Vc;
  int nceq;

  nceq = netpar->E.nceq;

  Vc = Varb[1];                                                        /* V */

  Varb[2] = Gammaf(Vc, netpar->E.thetah, netpar->E.sigmah);            /* h */

  Varb[3] = Gammaf(Vc, netpar->E.thetan, netpar->E.sigman);            /* n */

  Varb[4] = Gammaf(Vc, netpar->E.thetab, netpar->E.sigmab);            /* b */

  Varb[5] = Gammaf(Vc, netpar->E.thetaz, netpar->E.sigmaz);            /* z */

  Varb[nceq+1] = 1.0;                                                  /* tG */

  sAMPAinf = Gammaf(Vc, netpar->P.ths, netpar->P.sigs);
  Varb[nceq+2] = netpar->P.kf * sAMPAinf * Varb[nceq+1] /              /* sP */
  (netpar->P.kf * sAMPAinf * Varb[nceq+1] + (1.0 / netpar->P.tAMPA));

  /* Fix  this!!! */
  sNMDAinf= Gammaf(Vc, netpar->N.ths, netpar->N.sigs);                 /* sN */

  Varb[nceq+3] = netpar->N.kx * sNMDAinf /  
  (netpar->N.kx * sNMDAinf + ((1.0 - sNMDAinf) / netpar->N.tsrNMDA) );


  Varb[nceq+4] = netpar->N.kf * Varb[nceq+3] * Varb[nceq+1] /
  (netpar->N.kf * Varb[nceq+3] * Varb[nceq+1] * sNMDAinf + 
  (1.0 / netpar->N.tsdNMDA) );
}

/* This function updates the variables of an excitatory cell   */
void update_cell_GY(cell_par *parE, syn_Glu_par *synG, syn_AMPA_par *synP,
     syn_NMDA_par *synN, double Iapp, double *Varc, double *kout, int ion,
     int it, fl_st fl)
{
  double Vc, hc, nc, bc, zc, tGlu, sAMPA, xNMDA, sNMDA;
  double Minf, Pinf, Hinf, tauH, Ninf, tauN, Ainf, Binf, Zinf;
  double sAMPAinf, sNMDAinf;

  Vc = Varc[1];
  hc = Varc[2];
  nc = Varc[3];
  bc = Varc[4];
  zc = Varc[5];
  tGlu = Varc[parE->nceq+1];
  sAMPA = Varc[parE->nceq+2];
  xNMDA = Varc[parE->nceq+3];
  sNMDA = Varc[parE->nceq+4];

  Minf = Gammaf(Vc, parE->thetam, parE->sigmam);
  Pinf = Gammaf(Vc, parE->thetap, parE->sigmap);
  Hinf = Gammaf(Vc, parE->thetah, parE->sigmah);
  tauH = 1.0+7.5*Gammaf(Vc,-40.5,-6.0);
  Ninf = Gammaf(Vc, parE->thetan, parE->sigman);
  tauN = 1.0+5.0*Gammaf(Vc,-27.0,-15.0);
  Ainf = Gammaf(Vc, parE->thetaa, parE->sigmaa);
  Binf = Gammaf(Vc, parE->thetab, parE->sigmab);
  Zinf = Gammaf(Vc, parE->thetaz, parE->sigmaz);

  sAMPAinf = Gammaf(Vc, synP->ths, synP->sigs);
  sNMDAinf = Gammaf(Vc, synN->ths, synN->sigs);
  
  kout[1] = (-parE->gL * (Vc-parE->VL) - (parE->gNa * pow(Minf,3.0) * hc
               + parE->gNaP * Pinf) * (Vc-parE->VNa) - (parE->gKdr * pow(nc,
               4.0) + parE->gA * pow(Ainf,3.0) * bc + parE->gZ * zc)
               * (Vc-parE->VK) + Iapp) / parE->Cm;

  kout[2] = parE->phi * (Hinf - hc) / tauH;

  kout[3] = parE->phi * (Ninf - nc) / tauN;

  kout[4] = (Binf - bc) / parE->tauB;

  kout[5] = (Zinf - zc) / parE->tauZ;

  kout[parE->nceq+1] = -synG->kt * sAMPAinf * tGlu + synG->kv * (1.0-  tGlu);

  kout[parE->nceq+2] = synP->kf * sAMPAinf * tGlu * (1.0 - sAMPA) -
                       sAMPA / synP->tAMPA;

  /* kout[parE->nceq+3] = synN->kf * sNMDAinf * tGlu * (1.0 - sNMDA) -
     synN->kr * sNMDA; */

  kout[parE->nceq+3] = synN->kx * sNMDAinf * (1.0 - xNMDA) -
     (1.0 - sNMDAinf) * xNMDA / synN->tsrNMDA;

  kout[parE->nceq+4] = synN->kf * xNMDA * tGlu * (1.0 - sNMDA) -
     sNMDA / synN->tsdNMDA;
}
