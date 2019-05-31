/* This program simulates a network with excitatory coupling.               */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "nr.h"
#include "sa.h"
#include "sn.h"
#include "sgy.h"
#include "swb.h"

/* This function simulates the system for one parameter set. */
/* It looks for a firing state/ by reducing I_app.            */
/* The calculation is done only for geometry=a.              */
void vary_I_one_par(fl_st fl, int* rng_ptr, int sm, double *Vstore,
     avr_val *av)
{
  func_cell_model *fcm;
  net_par *netpar;
  run_par *runpar;
  Varb_ar Varbar;
  spike_time spktm;
  bur_pattern_str burps;
  vary_I_par vIpar;
  double Ip, gNp, tsdNp;
  double gNstore, tsdNstore, gAstore, gNaPstore;
  int iI, itIa, ieq, ion;

  read_input(&fcm, &netpar, &runpar, &burps, &vIpar, sm, fl);

  Varbar.E = dmatrix(1, netpar->E.non, 1, netpar->E.neq);
  spktm.E.tm_spike_a = dvector(1, netpar->E.non);
  spktm.E.tm_spike_b = dvector(1, netpar->E.non);
  spktm.E.t_per_sum  = dvector(1, netpar->E.non);
  spktm.E.t_per2_sum = dvector(1, netpar->E.non);
  spktm.E.cv         = dvector(1, netpar->E.non); 
  spktm.E.n_spike    = ivector(1, netpar->E.non);
  spktm.E.spike_per_tchi = dvector(1, netpar->E.non);

  Varbar.I = dmatrix(1, netpar->I.non, 1, netpar->I.neq);
  spktm.I.tm_spike_a = dvector(1, netpar->I.non);
  spktm.I.tm_spike_b = dvector(1, netpar->I.non);
  spktm.I.t_per_sum  = dvector(1, netpar->I.non);
  spktm.I.t_per2_sum = dvector(1, netpar->I.non);
  spktm.I.cv         = dvector(1, netpar->I.non); 
  spktm.I.n_spike    = ivector(1, netpar->I.non);
  spktm.I.spike_per_tchi = dvector(1, netpar->I.non);

  if (burps.anburp == 'y') bur_ps_create(&burps, fl);

  in_con(fcm, netpar, runpar, Varbar, rng_ptr, fl);

  gNstore = netpar->EE.gNMDA;
  gAstore = netpar->EE.gAMPA;
  tsdNstore = netpar->N.tsdNMDA;
  gNaPstore = netpar->E.gNaP;

  for (iI=0; iI<=vIpar.nI; iI++)
  {
    for (ion=burps.kon_s_min; ion<=burps.kon_s_max; ion++)
    {
      burps.nisi[ion] = 0;
      burps.previous_spike_time[ion] = 0.0;
    }

    if (iI == 0)
    {
      Ip = vIpar.Istart;
      if (vIpar.difstart == 'y')
      {
        netpar->EE.gNMDA = vIpar.gNstart;
        netpar->EE.gAMPA = vIpar.gAstart;
        netpar->N.tsdNMDA = vIpar.tsdNstart;
        netpar->E.gNaP = vIpar.gNaPstart;
      }
    }
    else
    {
      Ip = vIpar.Istart + iI * (vIpar.Iend - vIpar.Istart) / vIpar.nI;
      if (vIpar.difstart == 'y')
      {
        netpar->EE.gNMDA = vIpar.gNstart + iI * 
          (gNstore - vIpar.gNstart) / vIpar.nI;
        netpar->EE.gAMPA = vIpar.gAstart + iI * 
          (gAstore - vIpar.gAstart) / vIpar.nI;
        netpar->N.tsdNMDA = vIpar.tsdNstart + iI *
          (tsdNstore - vIpar.tsdNstart) / vIpar.nI;
        netpar->E.gNaP = vIpar.gNaPstart + iI * 
          (gNaPstore - vIpar.gNaPstart) / vIpar.nI;
      }
    }
    fprintf(fl.out, "iI=%d Ip=%lf gN=%lf tN=%lf gA=%lf gNaP=%lf", iI, Ip, 
    netpar->EE.gNMDA, netpar->N.tsdNMDA, netpar->EE.gAMPA, netpar->E.gNaP);
    printf("iI=%d Ip=%lf gN=%lf tN=%lf gA=%lf gNaP=%lf ", iI, Ip,
    netpar->EE.gNMDA, netpar->N.tsdNMDA, netpar->EE.gAMPA, netpar->E.gNaP);


    fprintf(fl.out, " Iapp=");
    for (itIa=0; itIa<=netpar->E.n_Iapp-1; itIa++)
    {
      netpar->E.Iapp[itIa] = Ip;
      fprintf(fl.out, " %lf",  netpar->E.Iapp[itIa]);
    }
    fprintf(fl.out, " Iapp_a=");
    for (itIa=0; itIa<=netpar->E.n_Iapp-1; itIa++)
    {
      netpar->E.Iapp_a[itIa] = Ip;
      fprintf(fl.out, " %lf", netpar->E.Iapp_a[itIa]);
    }
    fprintf(fl.out, "\n");
    
    n_run(fcm, netpar, *runpar, Varbar, &spktm, &burps, av, fl);
    fprintf(fl.col, "  \n");

    if (burps.anburp == 'y')
    {
      analyze_burst_pattern_all(runpar, &burps, av, fl);

      fprintf(fl.out, "iI=%d spikes_in_burst=%d sp_bhv=%c\n", iI,
      av->spikes_in_burst, av->sp_bhv);
    }
  }

  if (netpar->E.neq >= Meq)
  {
    printf("neq=%d > Meq=%d\n", netpar->E.neq, Meq);
    exit(0);
  }

  for (ieq=1; ieq<=netpar->E.neq; ieq++) Vstore[ieq] = Varbar.E[1][ieq];

  if (burps.anburp == 'y') bur_ps_free(&burps, fl);

  free_dmatrix(Varbar.E, 1, netpar->E.non, 1, netpar->E.neq);
  free_dvector(spktm.E.tm_spike_a, 1, netpar->E.non);
  free_dvector(spktm.E.tm_spike_b, 1, netpar->E.non);
  free_dvector(spktm.E.t_per_sum , 1, netpar->E.non); 
  free_dvector(spktm.E.t_per2_sum, 1, netpar->E.non); 
  free_dvector(spktm.E.cv        , 1, netpar->E.non);
  free_ivector(spktm.E.n_spike   , 1, netpar->E.non);
  free_dvector(spktm.E.spike_per_tchi, 1, netpar->E.non);

  free_dmatrix(Varbar.I, 1, netpar->I.non, 1, netpar->I.neq);
  free_dvector(spktm.I.tm_spike_a, 1, netpar->I.non);
  free_dvector(spktm.I.tm_spike_b, 1, netpar->I.non);
  free_dvector(spktm.I.t_per_sum , 1, netpar->I.non);
  free_dvector(spktm.I.t_per2_sum, 1, netpar->I.non);
  free_dvector(spktm.I.cv        , 1, netpar->I.non);
  free_ivector(spktm.I.n_spike   , 1, netpar->I.non);
  free_dvector(spktm.I.spike_per_tchi, 1, netpar->I.non);
  free_ivector(runpar->E.nwritear, 1, runpar->I.nwrite);
  free_ivector(runpar->I.nwritear, 1, runpar->I.nwrite); 
}

/* This function simulates the system for one parameter set, starting */
/* from Vstore.                                                       */
void one_par_Vstore(fl_st fl, int* rng_ptr, int sm, double *Vstore_in,
     double *Vstore_out, avr_val *av)
{
  func_cell_model *fcm;
  net_par *netpar;
  run_par *runpar;
  Varb_ar Varbar;
  spike_time spktm;
  bur_pattern_str burps;
  vary_I_par vIpar;
  double  Ip;
  int iI, itIa, ieq;

  read_input(&fcm, &netpar, &runpar, &burps, &vIpar, sm, fl);

  Varbar.E = dmatrix(1, netpar->E.non, 1, netpar->E.neq);
  spktm.E.tm_spike_a = dvector(1, netpar->E.non);
  spktm.E.tm_spike_b = dvector(1, netpar->E.non);
  spktm.E.t_per_sum  = dvector(1, netpar->E.non);
  spktm.E.t_per2_sum = dvector(1, netpar->E.non);
  spktm.E.cv         = dvector(1, netpar->E.non); 
  spktm.E.n_spike    = ivector(1, netpar->E.non);
  spktm.E.spike_per_tchi = dvector(1, netpar->E.non);

  Varbar.I = dmatrix(1, netpar->I.non, 1, netpar->I.neq);
  spktm.I.tm_spike_a = dvector(1, netpar->I.non);
  spktm.I.tm_spike_b = dvector(1, netpar->I.non);
  spktm.I.t_per_sum  = dvector(1, netpar->I.non);
  spktm.I.t_per2_sum = dvector(1, netpar->I.non);
  spktm.I.cv         = dvector(1, netpar->I.non); 
  spktm.I.n_spike    = ivector(1, netpar->I.non);
  spktm.I.spike_per_tchi = dvector(1, netpar->I.non);

  if (burps.anburp == 'y') bur_ps_create(&burps, fl);

  if (netpar->E.neq >= Meq)
  {
    printf("neq=%d > Meq=%d\n", netpar->E.neq, Meq);
    exit(0);
  }

  for (ieq=1; ieq<=netpar->E.neq; ieq++) Varbar.E[1][ieq] = Vstore_in[ieq];

  n_run(fcm, netpar, *runpar, Varbar, &spktm, &burps, av, fl);
  fprintf(fl.col, "  \n");
  
  if (burps.anburp == 'y')
  {
    analyze_burst_pattern_all(runpar, &burps, av, fl);
  
    fprintf(fl.out, "spikes_in_burst=%d sp_bhv=%c\n", av->spikes_in_burst,
    av->sp_bhv);
  }

  for (ieq=1; ieq<=netpar->E.neq; ieq++) Vstore_out[ieq] = Varbar.E[1][ieq];

  if (burps.anburp == 'y') bur_ps_free(&burps, fl);

  free_dmatrix(Varbar.E, 1, netpar->E.non, 1, netpar->E.neq);
  free_dvector(spktm.E.tm_spike_a, 1, netpar->E.non);
  free_dvector(spktm.E.tm_spike_b, 1, netpar->E.non);
  free_dvector(spktm.E.t_per_sum , 1, netpar->E.non); 
  free_dvector(spktm.E.t_per2_sum, 1, netpar->E.non); 
  free_dvector(spktm.E.cv        , 1, netpar->E.non);
  free_ivector(spktm.E.n_spike   , 1, netpar->E.non);
  free_dvector(spktm.E.spike_per_tchi, 1, netpar->E.non);


  free_dmatrix(Varbar.I, 1, netpar->I.non, 1, netpar->I.neq);
  free_dvector(spktm.I.tm_spike_a, 1, netpar->I.non);
  free_dvector(spktm.I.tm_spike_b, 1, netpar->I.non);
  free_dvector(spktm.I.t_per_sum , 1, netpar->I.non);
  free_dvector(spktm.I.t_per2_sum, 1, netpar->I.non);
  free_dvector(spktm.I.cv        , 1, netpar->I.non);
  free_ivector(spktm.I.n_spike   , 1, netpar->I.non);
  free_dvector(spktm.I.spike_per_tchi, 1, netpar->I.non);
  free_ivector(runpar->E.nwritear, 1, runpar->I.nwrite);
  free_ivector(runpar->I.nwritear, 1, runpar->I.nwrite); 
}

/* This function simulates the system for one parameter set */
void one_par(fl_st fl, int* rng_ptr, int sm, avr_val *av)
{
  func_cell_model *fcm;
  net_par *netpar;
  run_par *runpar;
  /* syn_coup syncoup; */
  Varb_ar Varbar;
  spike_time spktm;
  bur_pattern_str burps;
  vary_I_par vIpar;

  read_input(&fcm, &netpar, &runpar, &burps, &vIpar, sm, fl);

  Varbar.E = dmatrix(1, netpar->E.non, 1, netpar->E.neq);
  spktm.E.tm_spike_a = dvector(1, netpar->E.non);
  spktm.E.tm_spike_b = dvector(1, netpar->E.non);
  spktm.E.t_per_sum  = dvector(1, netpar->E.non);
  spktm.E.t_per2_sum = dvector(1, netpar->E.non);
  spktm.E.cv         = dvector(1, netpar->E.non); 
  spktm.E.n_spike    = ivector(1, netpar->E.non);
  spktm.E.spike_per_tchi = dvector(1, netpar->E.non);

  Varbar.I = dmatrix(1, netpar->I.non, 1, netpar->I.neq);
  spktm.I.tm_spike_a = dvector(1, netpar->I.non);
  spktm.I.tm_spike_b = dvector(1, netpar->I.non);
  spktm.I.t_per_sum  = dvector(1, netpar->I.non);
  spktm.I.t_per2_sum = dvector(1, netpar->I.non);
  spktm.I.cv         = dvector(1, netpar->I.non); 
  spktm.I.n_spike    = ivector(1, netpar->I.non);
  spktm.I.spike_per_tchi = dvector(1, netpar->I.non);

  switch (netpar->S.geometry)
  {
    case 's':
      /*
      determine_sparse_par(&netpar->I, &netpar->S, &netpar->II, *runpar, fl);
      sparse_coup(&netpar->I, &netpar->II, *runpar, rng_ptr, &syncoup, fl);
      */
      break;
    case 'o':
      g_od_subs(netpar, runpar, fl);
      break;
    case 'a':
      break;
    default:
      printf("Wrong geometry!!!\n");
      exit(0);
      break;
  }

  if (burps.anburp == 'y') bur_ps_create(&burps, fl);

  in_con(fcm, netpar, runpar, Varbar, rng_ptr, fl);
  n_run(fcm, netpar, *runpar, Varbar, &spktm, &burps, av, fl);

  free_dmatrix(Varbar.E, 1, netpar->E.non, 1, netpar->E.neq);
  free_dvector(spktm.E.tm_spike_a, 1, netpar->E.non);
  free_dvector(spktm.E.tm_spike_b, 1, netpar->E.non);
  free_dvector(spktm.E.t_per_sum , 1, netpar->E.non); 
  free_dvector(spktm.E.t_per2_sum, 1, netpar->E.non); 
  free_dvector(spktm.E.cv        , 1, netpar->E.non);
  free_ivector(spktm.E.n_spike   , 1, netpar->E.non);
  free_dvector(spktm.E.spike_per_tchi, 1, netpar->E.non);


  free_dmatrix(Varbar.I, 1, netpar->I.non, 1, netpar->I.neq);
  free_dvector(spktm.I.tm_spike_a, 1, netpar->I.non);
  free_dvector(spktm.I.tm_spike_b, 1, netpar->I.non);
  free_dvector(spktm.I.t_per_sum , 1, netpar->I.non);
  free_dvector(spktm.I.t_per2_sum, 1, netpar->I.non);
  free_dvector(spktm.I.cv        , 1, netpar->I.non);
  free_ivector(spktm.I.n_spike   , 1, netpar->I.non);
  free_dvector(spktm.I.spike_per_tchi, 1, netpar->I.non);
  free_ivector(runpar->E.nwritear, 1, runpar->I.nwrite);
  free_ivector(runpar->I.nwritear, 1, runpar->I.nwrite); 

  switch (netpar->S.geometry)
  {
  case 's':
    /* sparse_syn_free(netpar->E.non, &syncoup, fl); */
      break;
    case 'o':
      free_g_od_subs(netpar, runpar, fl);
      break;
    case 'a':
      break;
    default:
      printf("Wrong geometry!!!\n");
      exit(0);
      break;
  }

  if (burps.anburp == 'y')
  {
    analyze_burst_pattern_all(runpar, &burps, av, fl);

    fprintf(fl.out, "spikes_in_burst=%d sp_bhv=%c\n", av->spikes_in_burst,
    av->sp_bhv);

    bur_ps_free(&burps, fl);
  }
}

/* This function reads the data */
void read_input(func_cell_model **fcm_in, net_par **netpar_in,
     run_par **runpar_in, bur_pattern_str *burps, vary_I_par *vIpar, int sm, 
     fl_st fl)
{
  static net_par netpar;
  static run_par runpar;
  static func_cell_model fcm;
  int iwrite;
  char model_type[3];

  runpar.epsilon = 1.0e-10;

  rewind(fl.tmp);

  /* Reading the data */
  fscanf(fl.tmp, "E_CELL: %c%c model\n", &model_type[0], &model_type[1]);
  model_type[2] = '\n';

  if (strncmp(model_type, "GY", 2) == 0)  /* Golomb-Yaari model  */
    substiture_function_GY(&fcm.E, fl);
  else
  {
    printf("wrong model_type=%c%c!\n", model_type[0], model_type[1]);
    exit(0);
  }

  fcm.E.read_cell_par(&netpar.E, fl);
  
  fscanf(fl.tmp, "I CELL: %c%c model\n", &model_type[0], &model_type[1]);
  model_type[2] = '\n';

  if (strncmp(model_type, "WB", 2) == 0)       /* Wang-Buzsaki model */
    substiture_function_WB(&fcm.I, fl);

  fcm.I.read_cell_par(&netpar.I, fl);
  
  fscanf(fl.tmp, "SYNAPSE\n");
  fscanf(fl.tmp, "geometry=%c fshape=%c concur=%c\n", &netpar.S.geometry,
  &netpar.S.fshape, &netpar.S.concur);

  fscanf(fl.tmp, "Glu: kt=%lf kv=%lf Vrev=%lf\n", &netpar.G.kt, &netpar.G.kv,
  &netpar.G.Vrev);

  fscanf(fl.tmp, "AMPA: ths=%lf sigs=%lf kf=%lf tAMPA=%lf ivar=%d\n",
  &netpar.P.ths, &netpar.P.sigs, &netpar.P.kf, &netpar.P.tAMPA,
  &netpar.P.ivar);

  fscanf(fl.tmp, "NMDA: kx=%lf tsrNMDA=%lf kf=%lf tsdNMDA=%lf ivar=%d\n",
  &netpar.N.kx, &netpar.N.tsrNMDA, &netpar.N.kf, &netpar.N.tsdNMDA,
  &netpar.N.ivar);
  fscanf(fl.tmp, "ths=%lf sigs=%lf thetanp=%lf sigmanp=%lf zeromag=%c\n",
  &netpar.N.ths, &netpar.N.sigs, &netpar.N.thetanp, &netpar.N.sigmanp, 
  &netpar.N.zeromag);

  fscanf(fl.tmp, "GABAA: ths=%lf sigs=%lf kt=%lf kv=%lf kf=%lf kr=%lf "
  "Vrev=%lf\n", &netpar.A.ths, &netpar.A.sigs, &netpar.A.kt, &netpar.A.kv,
  &netpar.A.kf, &netpar.A.kr, &netpar.A.Vrev);

  fscanf(fl.tmp, "EE: gAMPA=%lf gNMDA=%lf Min=%lf sig=%lf\n", &netpar.EE.gAMPA,
  &netpar.EE.gNMDA, &netpar.EE.Min, &netpar.EE.sig);

  fscanf(fl.tmp, "IE: gGABAA=%lf Min=%lf sig=%lf\n", &netpar.IE.gGABAA,
  &netpar.IE.Min, &netpar.IE.sig);

  fscanf(fl.tmp, "EI: gAMPA=%lf gNMDA=%lf Min=%lf sig=%lf Nalpha=%c alpha=%lf"
  "\n", &netpar.EI.gAMPA, &netpar.EI.gNMDA, &netpar.EI.Min, &netpar.EI.sig,
  &netpar.EI.Nalpha, &netpar.EI.alpha);

  fscanf(fl.tmp, "II: gGABAA=%lf Min=%lf sig=%lf gel=%lf Mel=%lf\n",
  &netpar.II.gGABAA, &netpar.II.Min, &netpar.II.sig, &netpar.II.gel,
  &netpar.II.Mel);

  fscanf(fl.tmp, "GENERAL\n");
  fscanf(fl.tmp, "deltat=%lf nt=%d\n", &runpar.deltat, &runpar.nt);

  fscanf(fl.tmp, "E: nwrite=%d nwritear=", &runpar.E.nwrite);
  runpar.E.nwritear = ivector(1, runpar.E.nwrite);
  fscanf(fl.tmp, "%d", &runpar.E.nwritear[1]);
  for (iwrite=2; iwrite<=runpar.E.nwrite; iwrite++) 
    fscanf(fl.tmp, " %d", &runpar.E.nwritear[iwrite]);
  fscanf(fl.tmp, "\n");

  fscanf(fl.tmp, "I: nwrite=%d nwritear=", &runpar.I.nwrite);
  runpar.I.nwritear = ivector(1, runpar.I.nwrite);
  fscanf(fl.tmp, "%d", &runpar.I.nwritear[1]);
  for (iwrite=2; iwrite<=runpar.I.nwrite; iwrite++) 
    fscanf(fl.tmp, " %d", &runpar.I.nwritear[iwrite]);
  fscanf(fl.tmp, "\n");

  fscanf(fl.tmp, "twrite=%d tmcol=%d tmtrj=%d imtrj=%d tchi=%d traster=%d\n", 
  &runpar.twrite, &runpar.tmcol, &runpar.tmtrj, &runpar.imtrj, &runpar.tchi,
  &runpar.traster);
  fscanf(fl.tmp, "chirange=%lf\n", &runpar.chirange);
  fscanf(fl.tmp, "method=%c incond=%c smforce=%c firestop=%c\n",
  &runpar.method, &runpar.incond, &runpar.smforce, &runpar.firestop);
  fscanf(fl.tmp, "anburp=%c Ttransient=%lf tolerance=%lf delta=%lf mapval=%d"
  "\n", &burps->anburp, &burps->Ttransient, &burps->tolerance, &burps->delta,
  &burps->mapval);
  fscanf(fl.tmp, "Istart=%lf Iend=%lf nI=%d difstart=%c\n", &vIpar->Istart,
  &vIpar->Iend, &vIpar->nI, &vIpar->difstart);
  fscanf(fl.tmp, "gNstart=%lf tsdNstart=%lf gAstart=%lf gNaPstart=%lf\n",
  &vIpar->gNstart, &vIpar->tsdNstart, &vIpar->gAstart, &vIpar->gNaPstart);


  burps->misi = 10000;
  if (netpar.E.non == 1)
  {
    burps->kon_s_min = 1;
    burps->kon_s_max = 1;
  }
  else if (netpar.E.non <= 3)
  {
    burps->kon_s_min = 1;
    burps->kon_s_max = netpar.E.non;
  }
  else
  {
    burps->kon_s_min = (netpar.E.non / 4) + 1;
    burps->kon_s_max = 3 * netpar.E.non / 4;
    if (burps->kon_s_min >= burps->kon_s_max)
    {
      printf("kon_s_min=%d >= kon_s_max=%d\n", burps->kon_s_min,
      burps->kon_s_max);
      exit(0);
    }
  }

  /* printing flag */

  runpar.sm = sm;    /* 0 - print as if there is only one parameter set */
                     /* 1 - no such printing                            */
  if (runpar.smforce == 'p')      /* always print    */
    runpar.sm = 0;
  else if (runpar.smforce == 'a') /* always print (regardless adj_str)  */
    runpar.sm = 0;
  else if (runpar.smforce == 'n') /* always no print */
    runpar.sm = 1;
  else if (runpar.smforce != 'l') /* leave as is     */
  {
    printf("smforce should be either p or n or l or a !!! smforce=%c\n",
    runpar.smforce);
    exit(0);
  }
  fprintf(fl.out, "sm=%d\n", sm);

  /* updating data */

  if (!runpar.sm)
  {
    if (runpar.imtrj > netpar.E.non)
    {
      printf("imtrj=%d > non=%d !\n", runpar.imtrj, netpar.E.non);
      exit(0);
    }
  }

  runpar.spike_threshold = 0.0;

  *fcm_in = &fcm;
  *netpar_in = &netpar;
  *runpar_in = &runpar;
  runpar.E.nwrite = min(runpar.E.nwrite, netpar.E.non);

  if (runpar.tchi < 0) runpar.tchi = 0;
  if (runpar.tchi > runpar.nt) runpar.tchi = runpar.nt;

  for (iwrite=1; iwrite<=runpar.E.nwrite; iwrite++) 
  {
    if (runpar.E.nwritear[iwrite] > netpar.E.non) 
      runpar.E.nwritear[iwrite] = netpar.E.non;
  }

  if (netpar.S.geometry == 'o')
  {
    runpar.chi_range_min = (int) ((netpar.E.non / 2) -
      (runpar.chirange - runpar.epsilon) * netpar.E.rho);
    runpar.chi_range_max = (int) ((netpar.E.non / 2) +
      (runpar.chirange + runpar.epsilon) * netpar.E.rho);
  }
  else if (netpar.S.geometry == 'a')
  {
    runpar.chi_range_min = 1;
    runpar.chi_range_max = netpar.E.non;
  }
  else
  {
    printf("wrong geometry=%c\n", netpar.S.geometry);
    exit(0);
  }

  if (netpar.EI.Nalpha == 'y')
  {
    netpar.EI.gNMDA = netpar.EI.alpha * netpar.EI.gAMPA;
  }
  else if (netpar.EI.Nalpha != 'n')
  {
    printf("Nalpha=%c should be either y or n !!!\n", netpar.EI.Nalpha);
    exit(0);
  }

  /* writing the data */
  fcm.E.write_cell_par(&netpar.E, runpar, fl);
  fcm.I.write_cell_par(&netpar.I, runpar, fl);

  fprintf(fl.out, "SYNAPSE\n");
  fprintf(fl.out, "geometry=%c fshape=%c concur=%c\n", netpar.S.geometry,
  netpar.S.fshape, netpar.S.concur);

  fprintf(fl.out, "Glu: kt=%lf kv=%lf Vrev=%lf\n", netpar.G.kt, netpar.G.kv,
  netpar.G.Vrev);

  fprintf(fl.out, "AMPA: ths=%lf sigs=%lf kf=%lf tAMPA=%lf ivar=%d\n",
  netpar.P.ths, netpar.P.sigs, netpar.P.kf, netpar.P.tAMPA, netpar.P.ivar);

  fprintf(fl.out, "NMDA: kx=%lf tsrNMDA=%lf kf=%lf tsdNMDA=%lf ivar=%d\n",
  netpar.N.kx, netpar.N.tsrNMDA, netpar.N.kf, netpar.N.tsdNMDA, netpar.N.ivar);
  fprintf(fl.out, "ths=%lf sigs=%lf thetanp=%lf sigmanp=%lf zeromag=%c\n",
  netpar.N.ths, netpar.N.sigs, netpar.N.thetanp, netpar.N.sigmanp, 
  netpar.N.zeromag);

  fprintf(fl.out, "GABAA: ths=%lf sigs=%lf kt=%lf kv=%lf\nkf=%lf kr=%lf "
  "Vrev=%lf\n", netpar.A.ths, netpar.A.sigs, netpar.A.kt, netpar.A.kv,
  netpar.A.kf, netpar.A.kr, netpar.A.Vrev);

  fprintf(fl.out, "EE: gAMPA=%lf gNMDA=%lf Min=%lf sig=%lf\n", netpar.EE.gAMPA,
  netpar.EE.gNMDA, netpar.EE.Min, netpar.EE.sig);


  fprintf(fl.out, "IE: gGABAA=%lf Min=%lf sig=%lf\n", netpar.IE.gGABAA,
  netpar.IE.Min, netpar.IE.sig);

  fprintf(fl.out, "EI: gAMPA=%lf gNMDA=%lf Min=%lf sig=%lf\n", netpar.EI.gAMPA,
  netpar.EI.gNMDA, netpar.EI.Min, netpar.EI.sig);

  fprintf(fl.out, "II: gGABAA=%lf Min=%lf sig=%lf gel=%lf Mel=%lf\n",
  netpar.II.gGABAA, netpar.II.Min, netpar.II.sig, netpar.II.gel,
  netpar.II.Mel);

  fprintf(fl.out, "GENERAL\n");
  fprintf(fl.out, "deltat=%lf nt=%d\n", runpar.deltat, runpar.nt);

  fprintf(fl.out, "E: nwrite=%d nwritear=", runpar.E.nwrite);
  fprintf(fl.out, "%d", runpar.E.nwritear[1]);
  for (iwrite=2; iwrite<=runpar.E.nwrite; iwrite++) 
    fprintf(fl.out, " %d", runpar.E.nwritear[iwrite]);
  fprintf(fl.out, "\n");

  fprintf(fl.out, "I: nwrite=%d nwritear=", runpar.I.nwrite);
  fprintf(fl.out, "%d", runpar.I.nwritear[1]);
  for (iwrite=2; iwrite<=runpar.I.nwrite; iwrite++) 
    fprintf(fl.out, " %d", runpar.I.nwritear[iwrite]);
  fprintf(fl.out, "\n");

  fprintf(fl.out, "twrite=%d tmcol=%d tmtrj=%d imtrj=%d tchi=%d traster=%d\n", 
  runpar.twrite, runpar.tmcol, runpar.tmtrj, runpar.imtrj, runpar.tchi,
  runpar.traster);
  fprintf(fl.out, "chirange=%lf chi_range_min=%d chi_range_max=%d\n",
  runpar.chirange, runpar.chi_range_min, runpar.chi_range_max);
  fprintf(fl.out, "method=%c incond=%c smforce=%c firestop=%c\n",
  runpar.method, runpar.incond, runpar.smforce, runpar.firestop);
  fprintf(fl.out, "anburp=%c Ttransient=%lf tolerance=%lf delta=%lf mapval=%d"
  "\n", burps->anburp, burps->Ttransient, burps->tolerance, burps->delta,
  burps->mapval);
  fprintf(fl.out, "Istart=%lf Iend=%lf nI=%d difstart=%c\n", vIpar->Istart,
  vIpar->Iend, vIpar->nI, vIpar->difstart);
  fprintf(fl.out, "gNstart=%lf tsdNstart=%lf gAstart=%lf gNaPstart=%lf\n",
  vIpar->gNstart, vIpar->tsdNstart, vIpar->gAstart, vIpar->gNaPstart);
  fprintf(fl.out, "kon_s_min=%d kon_s_max=%d\n", burps->kon_s_min,
  burps->kon_s_max);
  fprintf(fl.out, "\n");
  fflush(fl.out);
}

/* This function generates the arrays in burps */
void bur_ps_create(bur_pattern_str *burps, fl_st fl)
{
  int ion;

  burps->isiar = dmatrix(burps->kon_s_min, burps->kon_s_max, 1, burps->misi);
  burps->nisi = ivector(burps->kon_s_min, burps->kon_s_max);
  burps->previous_spike_time = dvector(burps->kon_s_min, burps->kon_s_max);

  for (ion=burps->kon_s_min; ion<=burps->kon_s_max; ion++)
  {
      burps->nisi[ion] = 0;
      burps->previous_spike_time[ion] = 0.0;
  }
}

/* This function frees the arrays in burps */
void bur_ps_free(bur_pattern_str *burps, fl_st fl)
{
  free_dmatrix(burps->isiar, burps->kon_s_min, burps->kon_s_max, 1,
  burps->misi);
  free_dvector(burps->previous_spike_time, burps->kon_s_min, burps->kon_s_max);
}

void I_app_read(cell_par *cpar, fl_st fl)
{
  int itIa;

  fscanf(fl.tmp, "n_Iapp=%d nin_Iapp_a=%d nt_Iapp=", &cpar->n_Iapp,
  &cpar->nin_Iapp_a);
  for (itIa=0; itIa<=cpar->n_Iapp-2; itIa++) fscanf(fl.tmp, " %d", 
  &cpar->nt_Iapp[itIa]);
  fscanf(fl.tmp, " Iapp=");
  for (itIa=0; itIa<=cpar->n_Iapp-1; itIa++) fscanf(fl.tmp, " %lf", 
  &cpar->Iapp[itIa]);
  fscanf(fl.tmp, " Iapp_a=");
  for (itIa=0; itIa<=cpar->n_Iapp-1; itIa++) fscanf(fl.tmp, " %lf", 
  &cpar->Iapp_a[itIa]);
  fscanf(fl.tmp, "\n");
}

void I_app_write(cell_par *cpar, int nt, fl_st fl)
{
  int itIa;

  for (itIa=1; itIa<=cpar->n_Iapp-2; itIa++) cpar->nt_Iapp[itIa] += 
  cpar->nt_Iapp[itIa-1];
  cpar->nt_Iapp[cpar->n_Iapp-1] = nt;
  cpar->i_Iapp = 0;

  fprintf(fl.out, "n_Iapp=%d nin_Iapp_a=%d nt_Iapp=", cpar->n_Iapp,
  cpar->nin_Iapp_a);
  for (itIa=0; itIa<=cpar->n_Iapp-1; itIa++) fprintf(fl.out, " %d", 
  cpar->nt_Iapp[itIa]);
  fprintf(fl.out, "\n");
  fprintf(fl.out, "Iapp=");
  for (itIa=0; itIa<=cpar->n_Iapp-1; itIa++) fprintf(fl.out, " %lf", 
  cpar->Iapp[itIa]);
  fprintf(fl.out, " Iapp_a=");
  for (itIa=0; itIa<=cpar->n_Iapp-1; itIa++) fprintf(fl.out, " %lf", 
  cpar->Iapp_a[itIa]);
  fprintf(fl.out, "\n");
}

/* This function substitutes the synaptic coupling term for 1-d coupling */
void g_od_subs(net_par *netpar, run_par *runpar, fl_st fl)
{
  int ion;
  int non_E, non_I;

  non_E = netpar->E.non;
  non_I = netpar->I.non;

  netpar->EE.gvec = dvector(1-(non_E+1)/2 , non_E/2);
  g_syn_od_cal(netpar->EE.gvec, netpar->S.fshape, netpar->EE.sig, 
  netpar->E.non, netpar->E.rho, "EE", fl);

  netpar->IE.gvec = dvector(1-(non_I+1)/2 , non_I/2);
  g_syn_od_cal(netpar->IE.gvec, netpar->S.fshape, netpar->IE.sig, 
  netpar->I.non, netpar->I.rho, "IE", fl);

  netpar->EI.gvec = dvector(1-(non_E+1)/2 , non_E/2);
  g_syn_od_cal(netpar->EI.gvec, netpar->S.fshape, netpar->EI.sig, 
  netpar->E.non, netpar->E.rho, "EI", fl);

  netpar->II.gvec = dvector(1-(non_I+1)/2 , non_I/2);
  g_syn_od_cal(netpar->II.gvec, netpar->S.fshape, netpar->II.sig, 
  netpar->I.non, netpar->I.rho, "II", fl);

  if (!runpar->sm)
  {
    fprintf(fl.out, "\n EE    IE    EI    II\n");
    for (ion=1-(non_E+1)/2; ion<=non_E/2; ion++)
    {
      fprintf(fl.out, "%4d %10.5lf %10.5lf %10.5lf %10.5lf\n", ion,
      netpar->EE.gvec[ion], netpar->IE.gvec[ion], netpar->EI.gvec[ion],
      netpar->II.gvec[ion]);
    }
  }
}

/* This function calculates the synaptic coupling term for one 1-d case */
void g_syn_od_cal(double *gvec, char fshape, double sig, int non, double rho,
     char *syn_type, fl_st fl)
{
  double econst, sume, nor_const;
  double lneigr, lneigl;
  int ion, npre;
  char scaling;

  lneigr = lneigl = sig * rho;
  npre = non;
  scaling = 'a';
  fprintf(fl.out, "%s lneigr=%lf lneigl=%lf npre=%d scaling=%c\n", syn_type,
  lneigr, lneigl, npre, scaling);

  for (ion=1-(npre+1)/2; ion<=npre/2; ion++)
    gvec[ion]=0.0;

  if (npre <= 1) 
  {
    gvec[0] = 1.0;
  }
  else
  {
    if (fshape == 'e') 
    {
      econst = 2.0 * tanh(0.5/lneigl) * tanh(0.5/lneigr) / 
        ( tanh(0.5/lneigl) + tanh(0.5/lneigr) );
      nor_const = 0.0;
    
      /* Synaptic coupling before normalization */
      for (ion=1-npre/2; ion<=0; ion++) 
        gvec[ion]=econst*exp(-fabs(1.0*ion/lneigl));
      for (ion=1; ion<=npre/2-1; ion++) 
        gvec[ion]=econst*exp(-fabs(1.0*ion/lneigr));
      gvec[npre/2]=0.0;

      sume=0.0;
      /* Calculating the normalization constant */
      if (scaling == 'a')           /* all */
      {
        for (ion=1-npre/2; ion<=npre/2-1; ion++) sume=sume+gvec[ion];
        nor_const = 1.0 / sume;
      }
      else if (scaling == 'l')      /* left */
      {
        for (ion=1-npre/2; ion<=-1; ion++) sume=sume+gvec[ion];
        sume=sume+(gvec[0]/2.0);
        nor_const = (0.5)/sume;
      }
      else if (scaling == 'r')      /* right */
      {
        for (ion=1; ion<=npre/2-1; ion++) sume=sume+gvec[ion];
        sume=sume+(gvec[0]/2.0);
        nor_const = (0.5)/sume;
      }
      else
      {
       printf ("scaling should be a, l or r");
       exit(1);
      }
      /* Normalizing synaptic efficacies */
        for (ion=1-npre/2; ion<=npre/2-1; ion++) gvec[ion]=nor_const*gvec[ion];
    }
    else if (fshape == 's')
    {
      for (ion = -((int) lneigl); ion <= ((int) lneigr); ion++)
        gvec[ion] = 1.0 / (lneigl+lneigr+1.0);
    }
  }
}

/* This function frees the synaptic coupling term for 'o' geometry */
void free_g_od_subs(net_par *netpar, run_par *runpar, fl_st fl)
{
  int non_E;

  non_E = netpar->E.non;

  free_dvector(netpar->EE.gvec, 1-(non_E+1)/2 , non_E/2);
}

/* This function substitutes the initial conditions */
void in_con(func_cell_model *fcm, net_par *netpar, run_par *runpar,
     Varb_ar Varbar, int *rng_ptr, fl_st fl)
{
  double xx;
  int ion, ieq;
  char line[Mline];

  /* neuronal variables */

  if (runpar->incond == 'r')
  {
    fscanf(fl.tmp, "\n");
    fscanf(fl.tmp, "INITIAL CONDITIONS\n");
    fscanf(fl.tmp, "E\n");
    fgets(line, Mline, fl.tmp);
                /* "V     h     n     b     z     T     sP     xN     sN\n" */
    for (ion=1; ion<=netpar->E.non; ion++)  
    {
      for (ieq=1; ieq<=netpar->E.neq; ieq++)
        fscanf(fl.tmp, "%lf", &(Varbar.E[ion][ieq]));
      fscanf(fl.tmp, "\n");
    }

    fscanf(fl.tmp, "I\n");
    fgets(line, Mline, fl.tmp);   /* "V     h     n     T     sA\n" */
    for (ion=1; ion<=netpar->I.non; ion++)  
    {
      for (ieq=1; ieq<=netpar->I.neq; ieq++)
        fscanf(fl.tmp, "%lf", &(Varbar.I[ion][ieq]));
      fscanf(fl.tmp, "\n");
    }
  }
  else if (runpar->incond == 'b')
  {
    for (ion=1; ion<=netpar->E.non; ion++)
    {
      Varbar.E[ion][1] = netpar->E.Vinc1;
      fcm->E.steady_state_var(netpar, Varbar.E[ion], fl);
    }

    for (ion=1; ion<=netpar->I.non; ion++)
    {
      Varbar.I[ion][1] = netpar->I.Vinc1;
      fcm->I.steady_state_var(netpar, Varbar.I[ion], fl);
    }
  }
  else
  {
    printf("Wrong incond:%c\n", runpar->incond);
    exit(1);
  }

  if (!runpar->sm)
  {
    fprintf(fl.out, "\nInitial conditions\n");
    fprintf(fl.out, "E\n");
    for (ion=1; ion<=netpar->E.non; ion++)
    {
      fprintf(fl.out, "%5d", ion);
      fprintf(fl.out, " %10.5lf", Varbar.E[ion][1]);
      for (ieq=2; ieq<=netpar->E.neq; ieq++) 
        fprintf(fl.out, " %8.5lf", Varbar.E[ion][ieq]);
      fprintf(fl.out, "\n");
    }

    fprintf(fl.out, "I\n");
    for (ion=1; ion<=netpar->I.non; ion++)
    {
      fprintf(fl.out, "%5d", ion);
      fprintf(fl.out, " %10.5lf", Varbar.I[ion][1]);
      for (ieq=2; ieq<=netpar->I.neq; ieq++) 
        fprintf(fl.out, " %8.5lf", Varbar.I[ion][ieq]);
      fprintf(fl.out, "\n");
    }
  }
  fflush(fl.out);
}

/* This function solves the differential equations */
void n_run(func_cell_model *fcm, net_par *netpar, run_par runpar,
     Varb_ar Varbar, spike_time *spktm, bur_pattern_str *burps, avr_val *av,
     fl_st fl)
{
  syn_input synput;
  pop_aver pop;
  spike_in_time_step spk_ts;
  Varb_ar k0, k1, k2, k3, k4, Varc;
  array_convolve arc;
  Vold Vo;
  v_cal vcal;
  cell_end_already_fire cells_already_fire;
  double time;
  int it, ion, ieq;

  varb_create(&k0, netpar->E.non, netpar->E.neq, netpar->I.non, netpar->I.neq);
  varb_create(&k1, netpar->E.non, netpar->E.neq, netpar->I.non, netpar->I.neq);
  varb_create(&k2, netpar->E.non, netpar->E.neq, netpar->I.non, netpar->I.neq);
  varb_create(&k3, netpar->E.non, netpar->E.neq, netpar->I.non, netpar->I.neq);
  varb_create(&k4, netpar->E.non, netpar->E.neq, netpar->I.non, netpar->I.neq);
  varb_create(&Varc, netpar->E.non, netpar->E.neq, netpar->I.non, 
  netpar->I.neq);

  synp_create(&synput, netpar->E.non, netpar->I.non);
  /* synput.sc = syncoup; */

  spk_ts_create(&spk_ts, netpar->E.non, netpar->I.non);

  if (netpar->S.geometry == 'o')
  {
    create_one_d_coup_space(&arc, netpar->E.non, runpar, fl);
    vel_ar_determine_par(&vcal, netpar, fl);
    vel_ar_create(&vcal.E, netpar->E.non, netpar->E.rho, 'E', runpar, fl);
  }

  spike_detect_ar_create(&Vo.E, netpar->E.non);
  spike_detect_ar_create(&Vo.I, netpar->I.non);

  cells_already_fire.E = 0;
  it=0;
  time=0.0;

  stat_pop(Varbar.E, &pop.E, netpar->E.non, 'E', it, runpar, fl);
  stat_pop(Varbar.I, &pop.I, netpar->I.non, 'I', it, runpar, fl);

  if (((!runpar.sm) || (runpar.smforce == 'a')) && (runpar.tmcol >= runpar.nt))
     pr_fct(netpar, runpar, pop.E.Vpop, pop.I.Vpop, Varbar, time, it, fl);

  /* Main running loop */ 
  
  while ((++it) <= runpar.nt &&
         !check_fire_stop(runpar, cells_already_fire, fl))
  {
    time=it*runpar.deltat;
    update_i_Iapp(netpar, it, fl);

    if (runpar.method == 'r' || runpar.method == 'o')
                                                     /* Runge-Kutta-4 method */
    {
      one_integration_step(fcm,netpar,runpar,Varbar,k0,k1,0.0,              
      it,Varc,synput,&arc,fl);
      one_integration_step(fcm,netpar,runpar,Varbar,k1,k2,runpar.deltat/2.0,
      it,Varc,synput,&arc,fl);
      one_integration_step(fcm,netpar,runpar,Varbar,k2,k3,runpar.deltat/2.0,
      it,Varc,synput,&arc,fl);
      one_integration_step(fcm,netpar,runpar,Varbar,k3,k4,runpar.deltat,
      it,Varc,synput,&arc,fl);

      for (ion=1; ion<=netpar->E.non; ion++)
      {
        for (ieq=1; ieq<=netpar->E.neq; ieq++)
          Varbar.E[ion][ieq] = Varbar.E[ion][ieq]+(runpar.deltat/6.0)*
          (k1.E[ion][ieq]+2.0*(k2.E[ion][ieq]+
          k3.E[ion][ieq])+k4.E[ion][ieq]);
      }

     for (ion=1; ion<=netpar->I.non; ion++)
      {
        for (ieq=1; ieq<=netpar->I.neq; ieq++)
          Varbar.I[ion][ieq] = Varbar.I[ion][ieq]+(runpar.deltat/6.0)*
          (k1.I[ion][ieq]+2.0*(k2.I[ion][ieq]+
          k3.I[ion][ieq])+k4.I[ion][ieq]);
      }
    }
    else if (runpar.method =='t' || runpar.method == 'w')  
                                                     /* Runge-Kutta-2 method */
    {
      one_integration_step(fcm,netpar,runpar,Varbar,k0,k1,0.0,              
      it,Varc,synput,&arc,fl);
      one_integration_step(fcm,netpar,runpar,Varbar,k1,k2,runpar.deltat/2.0,
      it,Varc,synput,&arc,fl);

      for (ion=1; ion<=netpar->E.non; ion++)
      {
        for (ieq=1; ieq<=netpar->E.neq; ieq++)
          Varbar.E[ion][ieq] = Varbar.E[ion][ieq]+runpar.deltat*k2.E[ion][ieq];
      }

      for (ion=1; ion<=netpar->I.non; ion++)
      {
        for (ieq=1; ieq<=netpar->I.neq; ieq++)
          Varbar.I[ion][ieq] = Varbar.I[ion][ieq]+runpar.deltat*k2.I[ion][ieq];
      }
    }
    else if (runpar.method =='e')  /* Eulear method */
    {
      one_integration_step(fcm,netpar,runpar,Varbar,k0,k1,0.0,              
      it,Varc,synput,&arc,fl);

      for (ion=1; ion<=netpar->E.non; ion++)
      {
        for (ieq=1; ieq<=netpar->E.neq; ieq++)
          Varbar.E[ion][ieq]=Varbar.E[ion][ieq]+runpar.deltat*k1.E[ion][ieq];
      }

      for (ion=1; ion<=netpar->I.non; ion++)
      {
        for (ieq=1; ieq<=netpar->I.neq; ieq++)
          Varbar.I[ion][ieq]=Varbar.I[ion][ieq]+runpar.deltat*k1.I[ion][ieq];
      }
    }
    else
    {
      printf("wrong method!\n");
      exit(0);
    }

    stat_pop(Varbar.E, &pop.E, netpar->E.non, 'E', it, runpar, fl);
    stat_pop(Varbar.I, &pop.I, netpar->I.non, 'I', it, runpar, fl);

    if (((!runpar.sm) || (runpar.smforce == 'a')) && !(it%runpar.twrite) &&
      (it >= runpar.nt-runpar.tmcol))
     pr_fct(netpar, runpar, pop.E.Vpop, pop.I.Vpop, Varbar, time, it, fl);
    
    spike_detect_peak(it, time, netpar->E, netpar->S.geometry, runpar, 
    Varbar.E, &spktm->E, &spk_ts.E, Vo.E, &vcal.E, 'E', burps, fl);
    
    spike_detect_peak(it, time, netpar->I, netpar->S.geometry, runpar,
    Varbar.I, &spktm->I, &spk_ts.I, Vo.I, &vcal.I, 'I', burps, fl);
    
    if (netpar->S.geometry == 'o')
    {
      cells_already_fire.E = check_fire_all(&vcal.E, 'E', fl);
    }
  }


  fprintf(fl.out, "chiE=%lf chiI=%lf\n", pop.E.chi, pop.I.chi);
  av->chiE = pop.E.chi;

  if (netpar->S.geometry == 'o')
  {
    vel_cal(&vcal, &vcal.E, netpar->E.non, netpar->E.rho, 'E', runpar, fl);
    av->velE = vcal.E.vel;
    av->vel_endE = vcal.E.vel_end;
  }

  varb_free(&k0, netpar->E.non, netpar->E.neq, netpar->I.non, netpar->I.neq);
  varb_free(&k1, netpar->E.non, netpar->E.neq, netpar->I.non, netpar->I.neq);
  varb_free(&k2, netpar->E.non, netpar->E.neq, netpar->I.non, netpar->I.neq);
  varb_free(&k3, netpar->E.non, netpar->E.neq, netpar->I.non, netpar->I.neq);
  varb_free(&k4, netpar->E.non, netpar->E.neq, netpar->I.non, netpar->I.neq);
  varb_free(&Varc, netpar->E.non, netpar->E.neq, netpar->I.non, 
  netpar->I.neq);

  synp_free(&synput, netpar->E.non, netpar->I.non);

  spk_ts_free(&spk_ts, netpar->E.non, netpar->I.non);

  if (netpar->S.geometry == 'o')
  {
    free_one_d_coup_space(&arc, netpar->E.non, runpar, fl);
    vel_ar_free(&vcal.E, netpar->E.non);
  }

  spike_detect_ar_free(&Vo.E, netpar->E.non);
  spike_detect_ar_free(&Vo.I, netpar->I.non);
}

void varb_create(Varb_ar *kk, int nonE, int neqE, int nonI, int neqI)
{
  int ion, ieq;

  kk->E = dmatrix(1, nonE, 1, neqE);
  kk->I = dmatrix(1, nonI, 1, neqI);

  for (ion=1; ion<=nonE; ion++)
  {
    for (ieq=1; ieq<=neqE; ieq++)
      kk->E[ion][ieq]=0.0; 
  }

  for (ion=1; ion<=nonI; ion++)
  {
    for (ieq=1; ieq<=neqI; ieq++)
      kk->I[ion][ieq]=0.0; 
  }
}

void varb_free(Varb_ar *kk, int nonE, int neqE, int nonI, int neqI)
{
  free_dmatrix(kk->E, 1, nonE, 1, neqE);
  free_dmatrix(kk->I, 1, nonI, 1, neqI);
}

void spk_ts_create(spike_in_time_step *spk_ts, int nonE, int nonI)
{
  spk_ts->E.ion   = ivector(1, nonE);
  spk_ts->E.tpeak = dvector(1, nonE);
  spk_ts->E.n_fire = 0;

  spk_ts->I.ion   = ivector(1, nonI);
  spk_ts->I.tpeak = dvector(1, nonI);
  spk_ts->I.n_fire = 0;
}

void spk_ts_free(spike_in_time_step *spk_ts, int nonE, int nonI)
{
  free_ivector(spk_ts->E.ion  , 1, nonE);
  free_dvector(spk_ts->E.tpeak, 1, nonE);

  free_ivector(spk_ts->I.ion  , 1, nonI);
  free_dvector(spk_ts->I.tpeak, 1, nonI);
}

void synp_create(syn_input *synput, int nonE, int nonI)
{
     synput->EE_P = dvector(1, nonE);
     synput->EE_N = dvector(1, nonE);
     synput->IE_A = dvector(1, nonE);
     synput->EI_P = dvector(1, nonI);
     synput->EI_N = dvector(1, nonI);
     synput->II_A = dvector(1, nonI);
}

void synp_free(syn_input *synput, int nonE, int nonI)
{
     free_dvector(synput->EE_P, 1, nonE);
     free_dvector(synput->EE_N, 1, nonE);
     free_dvector(synput->IE_A, 1, nonE);
     free_dvector(synput->EI_P, 1, nonI);
     free_dvector(synput->EI_N, 1, nonI);
     free_dvector(synput->II_A, 1, nonI);
}

/* This function updates i_Iapp if necessary    */
void update_i_Iapp(net_par *netpar, int it, fl_st fl)
{
  if (it == 1)
  {
    fprintf(fl.out, "E it=%d i_Iapp=%d Iapp_a=%lf Iapp=%lf\n", it,
    netpar->E.i_Iapp, netpar->E.Iapp_a[netpar->E.i_Iapp],
    netpar->E.Iapp[netpar->E.i_Iapp]);

    fprintf(fl.out, "I it=%d i_Iapp=%d Iapp_a=%lf Iapp=%lf\n", it,
    netpar->I.i_Iapp, netpar->I.Iapp_a[netpar->I.i_Iapp],
    netpar->I.Iapp[netpar->I.i_Iapp]);
  }

  if (it > netpar->E.nt_Iapp[netpar->E.i_Iapp])
  {
    netpar->E.i_Iapp++;
    fprintf(fl.out, "E it=%d i_Iapp=%d Iapp_a=%lf Iapp=%lf\n", it,
    netpar->E.i_Iapp, netpar->E.Iapp_a[netpar->E.i_Iapp],
    netpar->E.Iapp[netpar->E.i_Iapp]);
  }

  if (it > netpar->I.nt_Iapp[netpar->I.i_Iapp])
  {
    netpar->I.i_Iapp++;
    fprintf(fl.out, "I it=%d i_Iapp=%d Iapp_a=%lf Iapp=%lf\n", it,
    netpar->I.i_Iapp, netpar->I.Iapp_a[netpar->I.i_Iapp],
    netpar->I.Iapp[netpar->I.i_Iapp]);
  }
}

/* This function computes one integration step */
void one_integration_step(func_cell_model *fcm, net_par *netpar,
     run_par runpar, Varb_ar Varbar, Varb_ar kin, Varb_ar kout, double delt,
     int it, Varb_ar Varc, syn_input synput, array_convolve *arc, fl_st fl)
{
  cell_par *parE = &netpar->E;
  cell_par *parI = &netpar->I;
  syn_Glu_par   *synG = &netpar->G;
  syn_AMPA_par  *synP = &netpar->P;
  syn_NMDA_par  *synN = &netpar->N;
  syn_GABAA_par *synA = &netpar->A;

  Varb_ar Varo;         /* Variables for one_integration_step */
  double Vc, epost_NMDA, syn_term_P, syn_term_N, syn_term_A;
  double Iapp;
  int ion, ieq;

  /* Runge-Kutta input variables */

  if (delt > runpar.epsilon)
  {
    for (ion=1; ion<=netpar->E.non; ion++)
    {
      for (ieq=1; ieq<=netpar->E.neq; ieq++)
        Varc.E[ion][ieq] = Varbar.E[ion][ieq] + delt*kin.E[ion][ieq];
    }
    Varo.E = Varc.E;

    for (ion=1; ion<=netpar->I.non; ion++)
    {
      for (ieq=1; ieq<=netpar->I.neq; ieq++)
        Varc.I[ion][ieq] = Varbar.I[ion][ieq] + delt*kin.I[ion][ieq];
    }
    Varo.I = Varc.I;
  }
  else
  {
    Varo.E = Varbar.E;
    Varo.I = Varbar.I;
  }

  /* Calculating the input field of the chemical and/or electrical coupling */
  if ( runpar.method == 'r' || runpar.method == 't' || runpar.method == 'e' ||
      (runpar.method == 'o' && delt < runpar.epsilon) ||
      (runpar.method == 'w' && delt < runpar.epsilon) )
    syn_field_cal(netpar, runpar, Varo, synput, arc, fl);

  /* Updating E cells */
  for (ion=1; ion<=netpar->E.non; ion++)
  {
    /* Intrinsic properties */
    Iapp = ion <= parE->nin_Iapp_a ? parE->Iapp_a[parE->i_Iapp] : 
                                     parE->Iapp[parE->i_Iapp];
    fcm->E.update_cell(parE, synG, synP, synN, Iapp, Varo.E[ion], kout.E[ion],
    ion, it, fl);

    /* Synaptic coupling */
    Vc = Varo.E[ion][1];

    if (netpar->N.zeromag == 'n')
       epost_NMDA = Gammaf(Vc, netpar->N.thetanp, netpar->N.sigmanp);
    else if (netpar->N.zeromag == 'y')
      epost_NMDA = 1.0;
    else
    {
      printf("wrong zeromag=%c\n", netpar->N.zeromag);
      exit(0);
    }

    syn_term_P = netpar->EE.gAMPA * synput.EE_P[ion] * (Vc - netpar->G.Vrev);
    syn_term_N = netpar->EE.gNMDA * synput.EE_N[ion] * epost_NMDA * 
                 (Vc - netpar->G.Vrev);
    syn_term_A = netpar->IE.gGABAA * synput.IE_A[ion] * (Vc - netpar->A.Vrev);
   kout.E[ion][1] -= (syn_term_P + syn_term_N + syn_term_A);
  }

  /* Updating I cells */
  for (ion=1; ion<=netpar->I.non; ion++)
  {
    /* Intrinsic properties */
    Iapp = ion <= parI->nin_Iapp_a ? parI->Iapp_a[parI->i_Iapp] : 
                                     parI->Iapp[parI->i_Iapp];
    fcm->I.update_cell(parI, synA, Iapp, Varo.I[ion], kout.I[ion], ion, it,
    fl);

    /* Synaptic coupling */
    Vc = Varo.I[ion][1];

    if (netpar->N.zeromag == 'n')
      epost_NMDA = Gammaf(Vc, netpar->N.thetanp, netpar->N.sigmanp);
    else if (netpar->N.zeromag == 'y')
      epost_NMDA = 1.0;

    syn_term_P = netpar->EI.gAMPA * synput.EI_P[ion] * (Vc - netpar->G.Vrev);
    syn_term_N = netpar->EI.gNMDA * synput.EI_N[ion] * epost_NMDA * 
                 (Vc - netpar->G.Vrev);
    syn_term_A = netpar->II.gGABAA * synput.II_A[ion] * (Vc - netpar->A.Vrev);
    kout.I[ion][1] -= (syn_term_P + syn_term_N + syn_term_A);
  }
}

/* This function calculates the synaptic fields of all neurons */
void syn_field_cal(net_par *netpar, run_par runpar, Varb_ar Varbar,
     syn_input synput, array_convolve *arc, fl_st fl)
{
  int ieq;

  if (netpar->S.geometry == 'a')     /* all-to-all coupling */
  {
    /* EE AMPA */
    ieq = netpar->E.nceq + netpar->P.ivar;
    all_coup_cal(Varbar.E, ieq, synput.EE_P, netpar->E.non, netpar->E.non,
    "EEP", fl);

    /* EE NMDA */
    ieq = netpar->E.nceq + netpar->N.ivar;
    all_coup_cal(Varbar.E, ieq, synput.EE_N, netpar->E.non, netpar->E.non,
    "EEN", fl);

    /* IE GABAA */
    ieq = netpar->I.nceq + 2;
    all_coup_cal(Varbar.I, ieq, synput.IE_A, netpar->E.non, netpar->I.non,
    "IEA", fl);

    /* EI AMPA */
    ieq = netpar->E.nceq + netpar->P.ivar;
    all_coup_cal(Varbar.E, ieq, synput.EI_P, netpar->I.non, netpar->E.non,
    "EIP", fl);

    /* EI NMDA */
    ieq = netpar->E.nceq + netpar->N.ivar;
    all_coup_cal(Varbar.E, ieq, synput.EI_N, netpar->I.non, netpar->E.non,
    "EIN", fl);

    /* II GABAA */
    ieq = netpar->I.nceq + 2;
    all_coup_cal(Varbar.I, ieq, synput.II_A, netpar->I.non, netpar->I.non,
    "IIA", fl);
  }
  else if (netpar->S.geometry == 'o')     /* 1-d connectivity */
  {
    /* EE AMPA */
    ieq = netpar->E.nceq + netpar->P.ivar;
    one_d_coup_cal(Varbar.E, ieq, netpar->EE.gvec, synput.EE_P, netpar->E.non, 
    netpar->E.non, arc, "EEP", fl);

    /* EE NMDA */
    ieq = netpar->E.nceq + netpar->N.ivar;
    one_d_coup_cal(Varbar.E, ieq, netpar->EE.gvec, synput.EE_N, netpar->E.non,
    netpar->E.non, arc, "EEN", fl);

    /* IE GABAA */
    ieq = netpar->I.nceq + 2;
    one_d_coup_cal(Varbar.I, ieq, netpar->IE.gvec, synput.IE_A, netpar->E.non,
    netpar->I.non, arc, "IEA", fl);

    /* EI AMPA */
    ieq = netpar->E.nceq + netpar->P.ivar;
    one_d_coup_cal(Varbar.E, ieq, netpar->EI.gvec, synput.EI_P, netpar->I.non,
    netpar->E.non, arc, "EIP", fl);

    /* EI NMDA */
    ieq = netpar->E.nceq + netpar->N.ivar;
    one_d_coup_cal(Varbar.E, ieq, netpar->EI.gvec, synput.EI_N, netpar->I.non,
    netpar->E.non, arc, "EIN", fl);

    /* II GABAA */
    ieq = netpar->I.nceq + 2;
    one_d_coup_cal(Varbar.I, ieq, netpar->II.gvec, synput.II_A, netpar->I.non,
    netpar->I.non, arc, "IIA", fl);
  }
}

/* This function calculates the synaptic fields of all neurons for       */
/* all-to-all coupling                                                   */
void all_coup_cal (double **VecVar, int ieq, double *syn_field, 
     int npost, int npre, char *coup_form, fl_st fl)
{
  double syn_all;
  int ipost, ipre;

  syn_all = 0.0;
  for (ipre=1; ipre<=npre; ipre++)
    syn_all = syn_all + VecVar[ipre][ieq];
  syn_all = syn_all / npre;
  
  for (ipost=1; ipost<=npost; ipost++)
    syn_field[ipost] = syn_all;
}

/* This function creates the space for the one_d_coup arays */
void create_one_d_coup_space(array_convolve *arc, int npre, run_par runpar,
     fl_st fl)
{
  int ion;
  arc->svecin = dvector(1, 2*npre);
  arc->gwrap  = dvector(1, 2*npre);
  arc->ans    = dvector(1, 4*npre);

  for (ion=1; ion<=2*npre; ion++) arc->svecin[ion] = 0.0;
  for (ion=1; ion<=2*npre; ion++) arc->gwrap[ion]  = 0.0;
  for (ion=1; ion<=4*npre; ion++) arc->ans[ion]    = 0.0;
}

/* This function frees the space for the one_d_coup auxiliary arays */
void free_one_d_coup_space(array_convolve *arc, int npre, run_par runpar,
     fl_st fl)
{
  free_dvector(arc->svecin ,1, 2*npre);
  free_dvector(arc->gwrap  ,1, 2*npre);
  free_dvector(arc->ans    ,1, 4*npre);
}

/* This function calculates the synaptic fields of all neurons for       */
/* 1-d architecture with open boundary conditions.                       */
void one_d_coup_cal(double **VecVar, int ieq, double *gvec, double *sing,
     int npost, int npre, array_convolve *arc, char *coup_form, fl_st fl)
{
  int ipre,ipost;

  for (ipre=1; ipre<=npre; ipre++) arc->svecin[ipre] = VecVar[ipre][ieq];
  for (ipre=1; ipre<=npre; ipre++) arc->svecin[npre+ipre] = 0.0;
  for (ipre=1; ipre<=npre/2; ipre++) arc->gwrap[ipre] = gvec[1-ipre];
  for (ipre=npre+npre/2+1; ipre<=2*npre; ipre++) arc->gwrap[ipre]=
     gvec[2*npre+1-ipre];
  for (ipre=npre/2+1; ipre<=npre+npre/2; ipre++) arc->gwrap[ipre]=0.0;

  convlv(arc->svecin, 2*npre, arc->gwrap, 2*npre, 1, arc->ans);

  for (ipost=1; ipost<=npost; ipost++) sing[ipost]=arc->ans[ipost];
}

/* This function calculates population statistics */
void stat_pop(double **Varbar, pop_aver_cell *popc, int non, char cell_type,
     int it, run_par runpar, fl_st fl)
{
  double sigma_Vpop_sq, *sigma_V_sq, pop_av_sigma_V_sq;
  int ion;

  popc->Vpop = 0.0;
  for (ion=runpar.chi_range_min; ion<=runpar.chi_range_max; ion++)
  {
    popc->Vpop += Varbar[ion][1];
  }
  popc->Vpop /= runpar.chi_range_max - runpar.chi_range_min + 1;

  if (it == runpar.nt - runpar.tchi + 1)
  {
    popc->V_avt    = dvector(runpar.chi_range_min, runpar.chi_range_max);
    popc->V_sq_avt = dvector(runpar.chi_range_min, runpar.chi_range_max);
    popc->Vpop_avt    = 0.0;
    popc->Vpop_sq_avt = 0.0;
    for (ion=runpar.chi_range_min; ion<=runpar.chi_range_max; ion++)
    {
      popc->V_avt[ion]    = 0.0;
      popc->V_sq_avt[ion] = 0.0;
    }
  }
  
  if (it >= runpar.nt - runpar.tchi + 1)
  {
    for (ion=runpar.chi_range_min; ion<=runpar.chi_range_max; ion++)
    {
      popc->V_avt[ion]    += Varbar[ion][1];
      popc->V_sq_avt[ion] += pow(Varbar[ion][1], 2.0);
    }
    popc->Vpop_avt    += popc->Vpop;
    popc->Vpop_sq_avt += pow(popc->Vpop, 2.0);
  }

  if (it >= runpar.nt)
  {
    sigma_V_sq = dvector(runpar.chi_range_min, runpar.chi_range_max);

    popc->Vpop_avt    /= runpar.tchi;
    popc->Vpop_sq_avt /= runpar.tchi;
    sigma_Vpop_sq = popc->Vpop_sq_avt - pow(popc->Vpop_avt, 2.0);

    pop_av_sigma_V_sq = 0.0;
    for (ion=runpar.chi_range_min; ion<=runpar.chi_range_max; ion++)
    {
      popc->V_avt[ion]    /= runpar.tchi;
      popc->V_sq_avt[ion] /= runpar.tchi;
      sigma_V_sq[ion] = popc->V_sq_avt[ion] - pow(popc->V_avt[ion], 2.0);
      pop_av_sigma_V_sq += sigma_V_sq[ion];
    }
    pop_av_sigma_V_sq /= runpar.chi_range_max - runpar.chi_range_min + 1;

    if (pop_av_sigma_V_sq < 0.0)
      popc->chi = -9999.0;
    else if (sigma_Vpop_sq < 0.0)
      popc->chi = -9999.0;
    else
    {
      popc->chi = sqrt(sigma_Vpop_sq / pop_av_sigma_V_sq);
    }
    fprintf(fl.out, "Vpop_avt=%lf Vpop_sq_avt=%lf sigma_Vpop_sq=%lf\n",
    popc->Vpop_avt, popc->Vpop_sq_avt, sigma_Vpop_sq);
    fprintf(fl.out, "pop_av_sigma_V_sq=%lf\n", pop_av_sigma_V_sq);


    free_dvector(popc->V_avt  , runpar.chi_range_min, runpar.chi_range_max);
    free_dvector(popc->V_sq_avt, runpar.chi_range_min, runpar.chi_range_max);
    free_dvector(sigma_V_sq  , runpar.chi_range_min, runpar.chi_range_max);
  }
}

/* This function writes the data on fl.col and fl.trj */
void pr_fct(net_par *netpar, run_par runpar, double VpopE, double VpopI,
     Varb_ar Varbar, double time, int it, fl_st fl)
{
  int ion, ieq;

  fprintf(fl.col, "%12.5lf %12.5lf", time, VpopE);

  for (ion=1; ion<=runpar.E.nwrite; ion++)
  {
    fprintf(fl.col, "%12.5lf", Varbar.E[runpar.E.nwritear[ion]][1]);
    /* if (ion%10 == 0) fprintf(fl.col, "\n"); */
  }

  fprintf(fl.col, "%12.5lf", VpopI);

  for (ion=1; ion<=runpar.I.nwrite; ion++)
  {
    fprintf(fl.col, "%12.5lf", Varbar.I[runpar.I.nwritear[ion]][1]);
    /* if (ion%10 == 0) fprintf(fl.col, "\n"); */
  }

  fprintf(fl.col, "\n");

  fflush(fl.col);
  
  if (it >= runpar.nt - runpar.tmtrj)
  {
    fprintf(fl.trj, "%12.5lf", time);

    for (ieq=1; ieq<=netpar->E.neq; ieq++)
      fprintf(fl.trj, "%12.5lf", Varbar.E[runpar.imtrj][ieq]);

    for (ieq=1; ieq<=netpar->I.neq; ieq++)
      fprintf(fl.trj, "%12.5lf", Varbar.I[runpar.imtrj][ieq]);

    fprintf(fl.trj, "\n");
  }
}

void spike_detect_ar_create(Vold_cell *Voc, int non)
{
  int ion;

  Voc->Voldma = dvector(1, non);
  Voc->Voldmb = dvector(1, non);
  Voc->after_max_vol = ivector(1, non);

  for (ion=1; ion<=non; ion++)
  {
    Voc->Voldma[ion] = -80.0;
    Voc->Voldmb[ion] = -80.0;
    Voc->after_max_vol[ion] = 0;
  }
}

void spike_detect_ar_free(Vold_cell *Voc, int non)
{
  free_dvector(Voc->Voldma, 1, non);
  free_dvector(Voc->Voldmb, 1, non);
  free_ivector(Voc->after_max_vol, 1, non);
}

/* This function detects the spike time */
void spike_detect_peak(int it, double time, cell_par par, char geometry,
     run_par runpar, double **Varbar, spike_time_cell *spktmc,
     spike_in_time_type *spk_tsc, Vold_cell Voc, v_cal_cell *vcalc,
     char cell_type, bur_pattern_str *burps, fl_st fl)
{
  double Volt_thresh=-30.0, after_min_vol=-33.0, xb, xc, tpeak, Vpeak;
  double V0, V1, V2;
  double t_per;
  int ion;

  if (it == 1)
  {

    for (ion=1; ion<=par.non; ion++)
    {
      spktmc->n_spike[ion] = 0;
      spktmc->tm_spike_a[ion] = 0;
      spktmc->t_per_sum[ion]=0.0;
      spktmc->t_per2_sum[ion]=0.0;
    }
  }

  spk_tsc->n_fire = 0;

 for (ion=1; ion<=par.non; ion++)
  {
    V0 = Varbar[ion][1];
    if ((!Voc.after_max_vol[ion]) && it>=10)
    {
      V1 = Voc.Voldma[ion];
      V2 = Voc.Voldmb[ion];
      if (V1 >= V0 && V1 >= V2 && V1 > Volt_thresh)  /* detect a spike */
      {
        xb = V2 - V0;
        xc = V0 - 2.0 * V1 + V2;
        if (fabs(xc) < runpar.epsilon)  
	{
          tpeak = time - runpar.deltat;
          Vpeak = V1;
	}
        else        
        {
          tpeak = time - runpar.deltat + 0.5 * (xb / xc) * runpar.deltat;
          Vpeak = V1 - 0.125 * xb * xb / xc;
	}
        Voc.after_max_vol[ion] = 1;

        spk_tsc->n_fire++;
        spk_tsc->ion[spk_tsc->n_fire] = ion;
        spk_tsc->tpeak[spk_tsc->n_fire] = tpeak;

        if (!runpar.sm)
	{
          fprintf(fl.ras, "%14.8e %d %lf %c\n", tpeak, ion, Vpeak, cell_type);
          fflush(fl.ras);
	}

        /* Storing the new ISI in burps */
        if (cell_type == 'E')
	{
          if ((burps->anburp == 'y') && (ion >= burps->kon_s_min) && 
              (ion <= burps->kon_s_max) )
	  {
            if (burps->previous_spike_time[ion] > burps->Ttransient)
            {
              burps->nisi[ion]++;
              if (burps->nisi[ion] > burps->misi)
	      {
                printf("ion=%d burps->nisi=%d > burps->misi=%d !\n", ion,
	        burps->nisi[ion], burps->misi);
                exit(0);
              }
              burps->isiar[ion][burps->nisi[ion]] = tpeak - 
                burps->previous_spike_time[ion];
	    }
            burps->previous_spike_time[ion] = tpeak;
          }
	}
        else if (cell_type != 'I')
	{
          printf("cell_type=%c should be E or I\n", cell_type);
          exit(0);
	}

        /* If this is the first spike of the cell, stote it in vcalc  */
        if (cell_type == 'E')
	{
          if (geometry == 'o')
	  {
            if (!vcalc->yspike[ion])
	    {
              vcalc->tspike[ion] = tpeak;
              vcalc->yspike[ion] = 1;
            }
	  }
        }

        if (it > runpar.nt - runpar.tchi + 1)
        /* Storing data for time-interval statistics */
	{
          if (++spktmc->n_spike[ion] >= 1)
          {
            spktmc->tm_spike_b[ion] = spktmc->tm_spike_a[ion];
            spktmc->tm_spike_a[ion] = tpeak;
            if (spktmc->tm_spike_b[ion] > (runpar.nt - runpar.tchi + 1) *
              runpar.deltat)
	    {
              t_per = spktmc->tm_spike_a[ion] - spktmc->tm_spike_b[ion];
              spktmc->t_per_sum[ion] += t_per;
              spktmc->t_per2_sum[ion] += t_per * t_per;
              /* if (!runpar.sm) 
                fprintf(fl.out, "%d %d %lf %lf %lf %lf %lf\n", ion, 
                spktmc->n_spike[ion], spktmc->tm_spike_b[ion], 
                spktmc->tm_spike_a[ion], t_per,
                spktmc->t_per_sum[ion], spktmc->t_per2_sum[ion]); */
	    }
          }
          else
          {
            spktmc->tm_spike_a[ion] = tpeak;
          }
	}
      }
    }
    else
    {
      if (V0 < after_min_vol) Voc.after_max_vol[ion] = 0;
    }

    /* updating the voltages */
    Voc.Voldmb[ion] = Voc.Voldma[ion];
    Voc.Voldma[ion] = Varbar[ion][1];
  }
}

void vel_ar_determine_par(v_cal *vcal, net_par *netpar, fl_st fl)
{
  double percent_omit_min, percent_omit_max;
  int iminE, imaxE;

  percent_omit_min = 0.3;
  percent_omit_max = 0.1;

  iminE = netpar->E.nin_Iapp_a + (int) (percent_omit_min * netpar->E.non) + 1;

  vcal->ion_min = iminE;

  imaxE = netpar->E.non - (int) (percent_omit_max * netpar->E.non) - 1;

  vcal->ion_max = imaxE;

  if (vcal->ion_max - vcal->ion_min > 3)
  {
    vcal->cal_possible = 1;
  }
  else
  {
    vcal->cal_possible = 0;
  }

  fprintf(fl.out, "iminE=%d ion_min=%d\n", iminE, vcal->ion_min);
  fprintf(fl.out, "imaxE=%d ion_max=%d cal_possible=%d\n", imaxE,
  vcal->ion_max, vcal->cal_possible);
  fprintf(fl.out, "\n");
  fflush(fl.out);
}

void vel_ar_create(v_cal_cell *vcalc, int non, double rho, char cell_type,
     run_par runpar, fl_st fl)
{
  int ion;

  vcalc->yspike = ivector(1, non);
  vcalc->tspike = dvector(1, non);

  for (ion=1; ion<=non; ion++)
  {
    vcalc->yspike[ion] = 0;
    vcalc->tspike[ion] = 0.0;
  }

  vcalc->ion_end_2 = non - (int) (rho + runpar.epsilon);
  vcalc->ion_end_1 = vcalc->ion_end_2 - 1;
  fprintf(fl.out, "%c ion_end_1=%d ion_end_2=%d\n", cell_type, 
  vcalc->ion_end_1, vcalc->ion_end_2);

  vcalc->ion_stop_min = (int) (0.9 * non);
  vcalc->ion_stop_max = non;
  fprintf(fl.out, "%c ion_stop_min=%d ion_stop_max=%d\n", cell_type,
  vcalc->ion_stop_min, vcalc->ion_stop_max);
}

void vel_ar_free(v_cal_cell *vcalc, int non)
{
  free_ivector(vcalc->yspike, 1, non);
  free_dvector(vcalc->tspike, 1, non);
}

void vel_cal(v_cal *vcal, v_cal_cell *vcalc, int non, int rho, char cell_type,
     run_par runpar, fl_st fl)
{
  double tav, iav, t2av, tiav, denom;
  int nspike, ion;

  if(!vcal->cal_possible)
  {
    fprintf(fl.out, "cal_possible=%d No calculation of velocity!\n",
    vcal->cal_possible);
    vcalc->vel = -999.0;
    vcalc->tav = -999.0;
    return;
  }

  nspike=0;
  tav = 0.0;
  iav = 0.0;
  t2av = 0.0;
  tiav = 0.0;

  for (ion=vcal->ion_min; ion<=vcal->ion_max; ion++)
  {
    if (vcalc->yspike[ion])
    {
      nspike++;
      tav = tav + vcalc->tspike[ion];
      iav = iav + ion;
      t2av = t2av + pow(vcalc->tspike[ion], 2.0);
      tiav = tiav + vcalc->tspike[ion] * ion;
      /* fprintf(fl.out, "ion=%d nspike=%d tspike=%lf\n", ion, nspike,
       vcalc->tspike[ion]); */
    }
  }

  if (nspike <= 3)
  {
    fprintf(fl.out, "nspike=%d No calculation of velocity!\n", nspike);
    vcalc->vel = -998.0;
    vcalc->tav = -998.0;
    return;
  }

  tav /= nspike;
  iav /= nspike;
  t2av /= nspike;
  tiav /= nspike;

  denom = t2av - tav*tav;
  if (denom > 0.0)
  {
    vcalc->vel = ((tiav - tav * iav) / denom) / rho;
    vcalc->tav = tav;
  }
  else 
  {
    vcalc->vel = -997.0;
    vcalc->tav = -997.0;
  }

  if (vcalc->tspike[vcalc->ion_end_1] > runpar.epsilon && 
      vcalc->tspike[vcalc->ion_end_2] > runpar.epsilon)
  {
    vcalc->vel_end = 1.0 * (vcalc->ion_end_2 - vcalc->ion_end_1) /
    ((vcalc->tspike[vcalc->ion_end_2] -vcalc->tspike[vcalc->ion_end_1]) * rho);
  }
  else
  {
    vcalc->vel_end = -996.0;
  }

  fprintf(fl.out, "%c vel=%lf tav=%lf vel_end=%lf\n", cell_type, vcalc->vel,
  vcalc->tav, vcalc->vel_end);
}

int check_fire_all(v_cal_cell *vcalc, char cell_type, fl_st fl)
{
  int ion, mul;

  mul=1;

  for (ion=vcalc->ion_stop_min; ion<=vcalc->ion_stop_max; ion++)
  {
    mul = mul * vcalc->yspike[ion];
  }
  return mul;
}

int check_fire_stop(run_par runpar, cell_end_already_fire cells_already_fire,
    fl_st fl)
{
  int stop_cal;

  stop_cal = 0;

  if (runpar.firestop == 'y')
  {
    if (cells_already_fire.E && cells_already_fire.I) stop_cal = 1;
  }

  return stop_cal;
}

/* Linear interpolation */
double lininter(double x1, double x2, double xc, double y1, double y2)
{
  double linter;

  linter = ((xc-x1)*y2+(x2-xc)*y1) / (x2-x1) ;
  return(linter);

}

/* This function analyses the bursting pattern for all the nuerons     */
/* and finds the number of spikes within a burst                       */
void analyze_burst_pattern_all(run_par *runpar, bur_pattern_str *burps, 
     avr_val *av, fl_st fl)
{
  double freq;
  int ion, iisi;

  if (!runpar->sm)
  {
  fprintf(fl.out, "imtrj=%d\n", runpar->imtrj);
  for (iisi=1; iisi<=burps->nisi[runpar->imtrj]; iisi++)
    fprintf(fl.out, "i=%d bur=%lf\n", iisi, burps->isiar[runpar->imtrj][iisi]);
  }

  burps->n_spk_in_bur_all = 0.0;
  burps->n_spk_in_bur_a_all = 0.0;
  burps->freq_all = 0.0;
  burps->duty_cycle_all = 0.0;

  for (ion=burps->kon_s_min; ion<=burps->kon_s_max; ion++)
  {
    if (!runpar->sm) fprintf(fl.out, "\nion=%d\n", ion);
    analyze_burst_pattern_one(ion, runpar, burps, fl);

    if (ion == runpar->imtrj) burps->xmnmx_imtrj = burps->xmnmx;

    burps->n_spk_in_bur_all += burps->n_spk_in_bur;
    burps->n_spk_in_bur_a_all += burps->n_spk_in_bur_a;

    if (fabs(burps->T_per_av) > runpar->epsilon)
    {
      freq = 1000.0 /  burps->T_per_av;
    }
    else
    {
      freq = -999.999;
      if (!runpar->sm) fprintf(fl.out, "ion=%d T_per_av=%lf freq=%lf\n", ion, 
      burps->T_per_av, freq);
    }
    burps->freq_all += freq;
    burps->duty_cycle_all += burps->duty_cycle;    
  }

  burps->n_spk_in_bur_all /= burps->kon_s_max - burps->kon_s_min + 1;
  burps->n_spk_in_bur_a_all /= burps->kon_s_max - burps->kon_s_min + 1;
  burps->freq_all /= burps->kon_s_max - burps->kon_s_min + 1;
  burps->duty_cycle_all /= burps->kon_s_max - burps->kon_s_min + 1;

  fprintf(fl.out, "\nn_spk_in_bur_all=%lf n_spk_in_bur_a_all=%lf\nfreq_all=%lf"
  " duty_cycle_all=%lf\n", burps->n_spk_in_bur_all, burps->n_spk_in_bur_a_all,
  burps->freq_all, burps->duty_cycle_all);
  av->xmnmx_imtrj = burps->xmnmx_imtrj;
  av->n_spk_in_bur_all = burps->n_spk_in_bur_all;
  av->n_spk_in_bur_a_all = burps->n_spk_in_bur_a_all;
  av->freq_all = burps->freq_all;
  av->duty_cycle_all = burps->duty_cycle_all;
  av->spikes_in_burst = determine_value(av, burps->mapval, fl);

  printf("av->n_spk_in_bur_a=%lf\n", av->n_spk_in_bur_a_all);
  if (av->n_spk_in_bur_a_all <= 0)
    av->sp_bhv = 'q';                /* quiescennt */
  else if (fabs(av->n_spk_in_bur_a_all - 1) <= runpar->epsilon)
    av->sp_bhv = 's';                /* spiking    */
  else if (av->n_spk_in_bur_a_all > 1 + runpar->epsilon)
    av->sp_bhv = 'b';                /* bursting   */

  fprintf(fl.out, "%lf %lf %lf %lf %c", av->xmnmx_imtrj, 
  av->n_spk_in_bur_a_all, av->freq_all, av->duty_cycle_all, av->sp_bhv);
  fprintf(fl.out, "\n");
}

/* This function analyses the bursting pattern for all the nuerons     */
/* and finds the number of spikes within a burst                       */
void analyze_burst_pattern_one(int ion, run_par *runpar,
     bur_pattern_str *burps, fl_st fl)
{
  double xmin, xmax, xmed, aver_large_isi, aver_large_isi_a;
  double xratio;
  int iisi, iisi_max, iisi_min, iisi_first, iisi_last;
  int num_large_isi;
  int nlarge, nsmall;

  if (burps->nisi[ion] == 0)
  {
    burps->n_spk_in_bur = 0;
    burps->n_spk_in_bur_a = 0;
    burps->xmnmx = -999.999;
    burps->T_per_av = -999.999;
    burps->duty_cycle = -999.999;
    if (!runpar->sm) fprintf(fl.out, "n_spk_in_bur=%lf n_spk_in_bur_a=%lf\n",
    burps->n_spk_in_bur, burps->n_spk_in_bur_a);
    return;
  }
  else if (burps->nisi[ion] == 1)
  {
    burps->n_spk_in_bur = -1;
    burps->n_spk_in_bur_a = -1;
    burps->xmnmx = -999.999;
    burps->T_per_av = -999.999;
    burps->duty_cycle = -999.999;
    if (!runpar->sm) fprintf(fl.out, "n_spk_in_bur=%lf n_spk_in_bur_a=%lf\n",
    burps->n_spk_in_bur, burps->n_spk_in_bur_a);
    return;
  }

  /* Find the maximal ISI */
  iisi_max = 1;
  iisi_min = 1;
  for (iisi=2; iisi<=burps->nisi[ion]; iisi++)
  {
    if (burps->isiar[ion][iisi] > burps->isiar[ion][iisi_max])
      iisi_max = iisi;
    if (burps->isiar[ion][iisi] < burps->isiar[ion][iisi_min])
      iisi_min = iisi;
  }
  if (!runpar->sm)
    fprintf(fl.out, "iisi_max=%d %lf iisi_min=%d %lf\n", iisi_max,
    burps->isiar[ion][iisi_max], iisi_min, burps->isiar[ion][iisi_min]);

  /* Find the first and "large" ISIs */
  iisi_first = 1;
  while (burps->isiar[ion][iisi_first] <
        (1.0 - burps->tolerance) * burps->isiar[ion][iisi_max])
  {
    iisi_first++;
    if (iisi_first > burps->nisi[ion])
    {
      printf("iisi_first=%d > burps->nisi=%d\n", iisi_first, burps->nisi);
      exit(0);
    }
  }

  iisi_last = burps->nisi[ion];
  while (burps->isiar[ion][iisi_last] <
         (1.0 - burps->tolerance) * burps->isiar[ion][iisi_max])
  {
    iisi_last--;
    if (iisi_last < 1)
    {
      printf("iisi_last=%d < 1\n", iisi_last);
      exit(0);
    }
  }

  if (!runpar->sm)
    fprintf(fl.out, "iisi_first=%d iisi_last=%d\n", iisi_first, iisi_last);

  /* Find numbers of intervals between large ISIs and the average large ISI*/
  num_large_isi = 1;
  aver_large_isi = 0;
  for (iisi=iisi_first+1; iisi<=iisi_last-1; iisi++)
  {
    if (burps->isiar[ion][iisi] >=
        (1.0 - burps->tolerance) * burps->isiar[ion][iisi_max])
    {
      num_large_isi++;
      aver_large_isi += burps->isiar[ion][iisi];
    }
  }

  burps->n_spk_in_bur = ((double) (iisi_last - iisi_first)) / 
                        ((double) num_large_isi);
  if (!runpar->sm)
    fprintf(fl.out, "num_large_isi=%d n_spk_in_bur=%lf\n", num_large_isi,
    burps->n_spk_in_bur);

  /* Find the ratio of the small and large ISI numbers */
  xmin = burps->isiar[ion][iisi_min];
  xmax = burps->isiar[ion][iisi_max];
  /* burps->xmnmx = 2.0 * (xmax - xmin) / (xmax + xmin); */
  burps->xmnmx = 2.0 * xmin / (xmax + xmin);
  if ( burps->xmnmx > 1.0 - burps->delta)
  {
   burps->n_spk_in_bur_a = 1;
  }
  else
  {    xmed = (xmax + xmin) / 2.0;
    nlarge = 2;
    nsmall = 0;
    aver_large_isi_a = burps->isiar[ion][iisi_first] +
                       burps->isiar[ion][iisi_last];
/* xx */ if (ion == runpar->imtrj && !runpar->sm) fprintf(fl.out, "i nl=%d aver_large_isi_a=%lf\n", nlarge, aver_large_isi_a);

    /* xx */ if (ion == runpar->imtrj && !runpar->sm) fprintf(fl.out, "xmin=%lf xmax=%lf xmed=%lf\n", xmin, xmax, xmed);
    for (iisi=iisi_first+1; iisi<=iisi_last-1; iisi++)
    {
      if (burps->isiar[ion][iisi] >= xmed)
      {
        nlarge++;
        aver_large_isi_a += burps->isiar[ion][iisi];
/* xx */ if (ion == runpar->imtrj && !runpar->sm) fprintf(fl.out, "+ nl=%d aver_large_isi_a=%lf\n", nlarge, aver_large_isi_a);
      }
      else
      {
        nsmall++;
      }
/* xx */ if (ion == runpar->imtrj && !runpar->sm) fprintf(fl.out, "%d %lf %d %d\n", iisi, burps->isiar[ion][iisi], nlarge, nsmall);
    }

    burps->n_spk_in_bur_a = ((double) (nsmall+nlarge-1)) /
                            ((double) (nlarge-1));
    if (!runpar->sm) fprintf(fl.out, "nlarge=%d nsmall=%d ", nlarge, nsmall);
  }  

  /* Find the average time period */

  /* xx */ if (ion == runpar->imtrj && !runpar->sm)
    fprintf(fl.out, "b aver_large_isi_a=%lf\n", aver_large_isi_a);
  aver_large_isi_a /= nlarge;
  /* xx */ if (ion == runpar->imtrj && !runpar->sm)
    fprintf(fl.out, "a aver_large_isi_a=%lf\n", aver_large_isi_a);
  burps->T_per_av = 0.0;
  for (iisi=iisi_first+1; iisi<=iisi_last; iisi++)
  {
    burps->T_per_av += burps->isiar[ion][iisi];
    if (ion == runpar->imtrj && !runpar->sm)
      fprintf(fl.out, "i=%d b=%lf Tp=%lf\n", iisi, burps->isiar[ion][iisi],
      burps->T_per_av);
  }

  if (burps->n_spk_in_bur_a > 1.0 + runpar->epsilon)
  {
    burps->T_per_av = burps->T_per_av / (nlarge-1);
    burps->duty_cycle = -aver_large_isi_a / burps->T_per_av + 1.0;
  }
  else
  {
    burps->T_per_av = burps->T_per_av / (iisi_last - iisi_first);
    burps->duty_cycle = 1.0;
  }
  if (!runpar->sm) 
    fprintf(fl.out, "T_per_av=%lf aver_large_isi_a=%lf duty_cycle=%lf\n",
    burps->T_per_av, aver_large_isi_a, burps->duty_cycle);

  xratio = burps->xmnmx / (2.0 - burps->xmnmx);
  if (!runpar->sm)
    fprintf(fl.out, "xmnmx=%lf xratio=%lf, n_spk_in_bur_a=%lf\n", burps->xmnmx,
    xratio, burps->n_spk_in_bur_a);
}

/* This function determines the integer value, which is a function of */
/* the firing pattern, and defines region in a map.                   */
int determine_value(avr_val *av, int mapval, fl_st fl)
{
  double epsl=1.0e-6;
  int dval;

  if (mapval == 1)
  {
    dval = av->spike_num;
    if (dval > 10) dval = 10;
    return(dval);
  }
  else if (mapval == 2)
  {
    dval = (int) (av->n_spk_in_bur_all + 0.5);
    if (dval > 10) dval = 10;
     return(dval);
  }
  else if (mapval == 3)
  {
    dval = (int) (av->n_spk_in_bur_a_all + 0.5);
    if (dval > 10) dval = 10;
    return(dval);
  }
  else if (mapval == 4)
  {
    dval = (int) (av->n_spk_in_bur_all + epsl);

    if (fabs(dval - av->n_spk_in_bur_all) > epsl) /* chaos */
    {
      dval = -1;
    }

    if (dval > 10) dval = 10;
     return(dval);
  }
  else
  {
    printf("wrong mapval=%d\n", mapval);
    exit(0);
  }
}
