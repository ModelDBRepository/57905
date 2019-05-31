#ifndef min
#define min(x,y) ( ((x) < (y)) ? (x) : (y) )
#define max(x,y) ( ((x) > (y)) ? (x) : (y) )
#endif

#ifndef Pi
#define Pi 3.1415926535897931
#endif

#ifndef NR_END
#define NR_END 1
#endif

#ifndef FREE_ARG
#define FREE_ARG char*
#endif

#ifndef div_nz
#define div_nz(x,y) ( (fabs(y) < (runpar.epsilon)) ? (0) : (x/y) )
#endif

#ifndef Gammaf
#define Gammaf(VV, theta,sigma) ( 1.0/(1.0+exp(-(VV-(theta))/(sigma))) )
#endif

#ifndef Meq
#define Meq 101
#endif

/* Structure Declaration */

typedef struct fl_st{
  FILE *in, *tmp, *avr, *out, *col, *trj, *ras;
} fl_st;


typedef struct avr_val{
  double chiE, velE, vel_endE;
  double n_spk_in_bur, n_spk_in_bur_a, T_per_av, duty_cycle, xmnmx;
  double n_spk_in_bur_all, n_spk_in_bur_a_all, freq_all, duty_cycle_all;
  double xmnmx_imtrj;
  int spike_num, spikes_in_burst;
  char sp_bhv;
} avr_val;

typedef struct col_write{
  int nwrite, *nwritear;
} col_write;

typedef struct run_par{
  col_write E, I;
  double epsilon, deltat, spike_threshold, chirange;
  int nt, twrite, tmcol, tmtrj, imtrj, tchi, traster;
  int chi_range_min, chi_range_max;
  int sm;
  char method, incond, smforce, firestop;
} run_par;

/* Function Declaration */

/* lfibrng6a.c */
int **init_rng_s_dbl(int ngen, int length, int seed);
double get_rn_dbl(int *genptr);
