void substiture_function_GY(func_cell_model_E *fcmE, fl_st fl);
void read_cell_par_GY(cell_par *cpar, fl_st fl);
void write_cell_par_GY(cell_par *cpar, run_par runpar, fl_st fl);
void steady_state_var_GY(net_par *netpar, double *Varb, fl_st fl);
void update_cell_GY(cell_par *parE, syn_Glu_par *synG, syn_AMPA_par *synP,
     syn_NMDA_par *synN, double Iapp, double *Varc, double *kout, int ion,
     int it, fl_st fl);
