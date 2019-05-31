#define Mpar 500
#define Mline 2000
#define Mword 50
#define Mdatraw 200
#define Mdatcol 1000

/* Structure Declaration */

typedef struct scan_val{
  double eps;
  char scan_type, bhv;
  double parmina, parmaxa, para_ar[Mpar];
  int  seed, npara, ipara, npta, nrepeata;
  char par1a[Mword], par2a[Mword];

  double parmin, parmax, par_ar[Mpar];
  int npar, ipar, npt, nrepeat, irepeat;
  char par1[Mword], par2[Mword];
} scan_val;


typedef struct dat_str{
  char datar[Mdatraw][Mdatcol];
  int n_datar;
} dat_str;

/* Function Declaration */
void read_first_input_line(scan_val *sval, fl_st fl);
void read_first_input_lines_f(scan_val *sval, fl_st fl);
void read_parameters_from_line(char *line, char *par1, char *par2, int *npt,
     double *parmin, double *parmax, int *npar, double *par_ar, fl_st fl);

void update_file_old(scan_val *sval, int skip_lines, fl_st fl);
int process_line_old(scan_val *sval, char line[], FILE *ftmp);
int process_line_no_colon_old(scan_val *sval, char line[], FILE *ftmp);
int process_line_yes_colon_old(scan_val *sval, char line[], FILE *ftmp);

void vary_I_one_par(fl_st fl, int* rng_ptr, int sm, double *Vstore,
     avr_val *av);
void one_par_Vstore(fl_st fl, int* rng_ptr, int sm, double *Vstore_in,
     double *Vstore_out, avr_val *av);
void one_par(fl_st fl, int *rng_ptr, int sm, avr_val *av);
void write_avr(scan_val sval, avr_val av, int savr, fl_st fl);
void border_find(scan_val *sval, int *rng_ptr, int sm, avr_val *av, fl_st fl);
int border_cal(double x1, char n1, double x2, char n2, scan_val *sval,
    int *rng_ptr, int sm, avr_val *av, fl_st fl);
void update_and_run_two_par(scan_val *sval, fl_st fl);
void read_file_in(dat_str *datstr, char scan_type, fl_st fl);
void update_file(char *par1, char *par2, int npt, double par_ar, int ipar,
     dat_str *datstr, fl_st fl);
int process_line(char *par1, char *par2, int npt, double par_ar, int ipar,
    char line[], FILE *ftmp);
void border_find_ic(scan_val *sval, int *rng_ptr, int sm, avr_val *av,
     fl_st fl);
void locate_border_behavior(int ipar_yes, double *Vstore_in,
     double *Vstore_out, int ipar_no, scan_val *sval, int *rng_ptr, int sm,
     avr_val *av, fl_st fl);
