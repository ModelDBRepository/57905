#ifndef Mline
#define Mline 120
#endif

#define MIapp 10

typedef struct cell_par{
  double gNa, gNaP, gKdr, gKdr1, gKdr3, gA, gZ, gL;
  double Iapp[MIapp], Iapp_a[MIapp];
  int n_Iapp, nin_Iapp_a, nt_Iapp[MIapp], i_Iapp;
  double Cm, VNa, VK, VL;
  double phi;
  double thetam, sigmam, thetap, sigmap;
  double thetah, sigmah, thetan, sigman;
  double thetaa, sigmaa, thetab, sigmab, tauB;
  double thetaz, sigmaz, tauZ;
  double vsh;
  double rho, Vinc1, Vinc2;
  int non, nceq, nseq, neq;
} cell_par;

typedef struct syn_par_all{
  char geometry, fshape, concur;
} syn_par_all;

typedef struct syn_Glu_par{
  double kt, kv, Vrev;
} syn_Glu_par;

typedef struct syn_AMPA_par{
  double ths, sigs, kf, tAMPA; /* kr */
  int ivar;
} syn_AMPA_par;

typedef struct syn_NMDA_par{
  double ths, sigs, kx, tsrNMDA, kf, tsdNMDA, thetanp, sigmanp; /* ky, kr */
  int ivar;
  char zeromag;
} syn_NMDA_par;

typedef struct syn_GABAA_par{
  double ths, sigs, kt, kv, kf, kr, Vrev;
} syn_GABAA_par;

typedef struct syn_coup_par{
  double gAMPA, gNMDA, gGABAA, Min, sig, rho_in;
  double gel, Mel, rho_el;
  double *gvec;
  double alpha;
  char Nalpha;
} syn_coup_par;

typedef struct net_par{
  cell_par E, I;           /* E,I cells */
  syn_par_all S;           /* Synapses - general parameters */
  syn_Glu_par G;
  syn_AMPA_par P;
  syn_NMDA_par N;
  syn_GABAA_par A;
  syn_coup_par EE;         /* E -> E connections */
  syn_coup_par IE;         /* I -> E connections */
  syn_coup_par EI;         /* E -> I connections */
  syn_coup_par II;         /* I -> I connections */
} net_par;

typedef struct Varb_ar{
  double **E, **I;
} Varb_ar;

typedef struct pop_aver_cell{
  double Vpop, chi;
  double Vpop_avt, Vpop_sq_avt, *V_avt, *V_sq_avt;
} pop_aver_cell;

typedef struct pop_aver{
  pop_aver_cell E, I;
} pop_aver;

typedef struct spike_time_cell{
  double *tm_spike_a, *tm_spike_b;
  double *t_per_sum, *t_per2_sum, *cv;
  double *spike_per_tchi;
  int *n_spike;
} spike_time_cell;

typedef struct spike_time{
  spike_time_cell E, I;
} spike_time;

typedef struct spike_in_time_type{
  double *tpeak;
  int n_fire, *ion;
} spike_in_time_type;

typedef struct spike_in_time_step{
  spike_in_time_type E, I;
} spike_in_time_step;

typedef struct Vold_cell{
  double *Voldma, *Voldmb;
  int *after_max_vol;
} Vold_cell;

typedef struct Vold{
  Vold_cell E, I;
} Vold;

typedef struct syn_coup{
  int **elect_mat, *m_elect, **inhib_mat, *m_inhib;
} syn_coup;

typedef struct syn_input{
  double *EE_P, *EE_N;
  double *IE_A;
  double *EI_P, *EI_N;
  double *II_A;
} syn_input;

typedef struct array_convolve{
  double *svecin, *gwrap, *ans;
} array_convolve;

typedef struct v_cal_cell{
  int *yspike;
  int ion_end_1, ion_end_2;
  int ion_stop_min, ion_stop_max;
  double *tspike, vel, tav, vel_end;
} v_cal_cell;

typedef struct v_cal{
  int ion_min, ion_max;
  int cal_possible;
  v_cal_cell E, I;
} v_cal;

typedef struct cell_end_already_fire{
  int E, I;
} cell_end_already_fire;

typedef struct func_cell_model_E
{
  void (*read_cell_par)(cell_par *cpar, fl_st fl);
  void (*write_cell_par)(cell_par *cpar, run_par runpar, fl_st fl);
  void (*steady_state_var)(net_par *netpar, double *Varb, fl_st fl);
  void (*update_cell)(cell_par *parE, syn_Glu_par *synG, syn_AMPA_par *synP,
     syn_NMDA_par *synN, double Iapp, double *Varc, double *kout, int ion,
     int it, fl_st fl);
} func_cell_model_E;

typedef struct func_cell_model_I
{
void (*read_cell_par)(cell_par *cpar, fl_st fl);
void (*write_cell_par)(cell_par *cpar, run_par runpar, fl_st fl);
void (*steady_state_var)(net_par *netpar, double *Varb, fl_st fl);
void (*update_cell)(cell_par *parI, syn_GABAA_par *synA, double Iapp,
     double *Varc, double *kout, int ion, int it, fl_st fl);
} func_cell_model_I;

typedef struct func_cell_model
{
  func_cell_model_E E;
  func_cell_model_I I;
} func_cell_model;

typedef struct bur_pattern_str{
  double **isiar, n_spk_in_bur, n_spk_in_bur_a;
  double *previous_spike_time;
  double T_per_av, duty_cycle, xmnmx;
  double Ttransient, tolerance, delta;
  double xmnmx_imtrj;
  double n_spk_in_bur_all, n_spk_in_bur_a_all, freq_all, duty_cycle_all;
  int misi, *nisi;
  int kon_s_min, kon_s_max;
  int mapval;
  char anburp;
} bur_pattern_str;

typedef struct vary_I_par{
  double Istart, Iend;
  double gNstart, tsdNstart, gAstart, gNaPstart;
  int nI;
  char difstart;
} vary_I_par;

/*
typedef struct vary_I_par{ 
  double Istart, Iend;
  int nI;
} vary_I_apr;
*/

/* Function Declaration */

void read_input(func_cell_model **fcm_in, net_par **netpar_in,
     run_par **runpar_in, bur_pattern_str *burps, vary_I_par *vIpar, int sm, 
     fl_st fl);
void bur_ps_create(bur_pattern_str *burps, fl_st fl);
void bur_ps_free(bur_pattern_str *burps, fl_st fl);
void I_app_read(cell_par *cpar, fl_st fl);
void I_app_write(cell_par *cpar, int nt, fl_st fl);
void g_od_subs(net_par *netpar, run_par *runpar, fl_st fl);
void g_syn_od_cal(double *gvec, char fshape, double sig, int non, double rho,
     char *syn_type, fl_st fl);
void free_g_od_subs(net_par *netpar, run_par *runpar, fl_st fl);
void in_con(func_cell_model *fcm, net_par *netpar, run_par *runpar,
     Varb_ar Varbar, int *rng_ptr, fl_st fl);
void n_run(func_cell_model *fcm, net_par *netpar, run_par runpar,
     Varb_ar Varbar, spike_time *spktm, bur_pattern_str *burps, avr_val *av,
     fl_st fl);
void varb_create(Varb_ar *kk, int nonE, int neqE, int nonI, int neqI);
void varb_free(Varb_ar *kk, int nonE, int neqE, int nonI, int neqI);
void spk_ts_create(spike_in_time_step *spk_ts, int nonE, int nonI);
void spk_ts_free(spike_in_time_step *spk_ts, int nonE, int nonI);
void synp_create(syn_input *synput, int nonE, int nonI);
void synp_free(syn_input *synput, int nonE, int nonI);
void update_i_Iapp(net_par *netpar, int it, fl_st fl);
void one_integration_step(func_cell_model *fcm, net_par *netpar,
     run_par runpar, Varb_ar Varbar, Varb_ar kin, Varb_ar kout, double delt,
     int it, Varb_ar Varc, syn_input synput, array_convolve *arc, fl_st fl);
void syn_field_cal(net_par *netpar, run_par runpar, Varb_ar Varbar,
     syn_input synput, array_convolve *arc, fl_st fl);
void all_coup_cal (double **VecVar, int ieq, double *syn_field, 
     int npost, int npre, char *coup_form, fl_st fl);
void create_one_d_coup_space(array_convolve *arc, int npre, run_par runpar,
     fl_st fl);
void free_one_d_coup_space(array_convolve *arc, int npre, run_par runpar,
     fl_st fl);
void one_d_coup_cal(double **VecVar, int ieq, double *gvec, double *sing,
     int npost, int npre, array_convolve *arc, char *coup_form, fl_st fl);
void stat_pop(double **Varbar, pop_aver_cell *popc, int non, char cell_type,
     int it, run_par runpar, fl_st fl);
void pr_fct(net_par *netpar, run_par runpar, double VpopE, double VpopI,
     Varb_ar Varbar, double time, int it, fl_st fl);
void spike_detect_ar_create(Vold_cell *Voc, int non);
void spike_detect_ar_free(Vold_cell *Voc, int non);
void spike_detect_peak(int it, double time, cell_par par, char geometry,
     run_par runpar, double **Varbar, spike_time_cell *spktmc,
     spike_in_time_type *spk_tsc, Vold_cell Voc, v_cal_cell *vcalc,
     char cell_type, bur_pattern_str *burps, fl_st fl);
void vel_ar_determine_par(v_cal *vcal, net_par *netpar, fl_st fl);
void vel_ar_create(v_cal_cell *vcalc, int non, double rho, char cell_type,
     run_par runpar, fl_st fl);
void vel_ar_free(v_cal_cell *vcalc, int non);
void vel_cal(v_cal *vcal, v_cal_cell *vcalc, int non, int rho, char cell_type,
     run_par runpar, fl_st fl);
int check_fire_all(v_cal_cell *vcalc, char cell_type, fl_st fl);
int check_fire_stop(run_par runpar, cell_end_already_fire cells_already_fire,
    fl_st fl);
double lininter(double x1, double x2, double xc, double y1, double y2);
void analyze_burst_pattern_all(run_par *runpar, bur_pattern_str *burps, 
     avr_val *av, fl_st fl);
void analyze_burst_pattern_one(int ion, run_par *runpar,
     bur_pattern_str *burps, fl_st fl);
