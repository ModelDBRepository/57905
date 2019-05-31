void substiture_function_WB(func_cell_model_I *fcmE, fl_st fl);
void read_cell_par_WB(cell_par *cpar, fl_st fl);
void write_cell_par_WB(cell_par *cpar, run_par runpar, fl_st fl);
void steady_state_var_WB(net_par *netpar, double *Varb, fl_st fl);
void update_cell_WB(cell_par *parI, syn_GABAA_par *synA, double Iapp,
     double *Varc, double *kout, int ion, int it, fl_st fl);
