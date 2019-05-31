/* This program calculate the phase reduction of a conductance-based neuron */
/* model.                                                                   */
/* In this file, the loop over the parameter is executed, and the           */
/* parameters of each run are determined.                                   */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "sa.h"
#include "sp.h"

main(int argc, char*argv[])
{
  scan_val sval;
  avr_val av;
  fl_st fl;
  int **ran_gen, *ran_one, ngen, length;
  int sm, savr;
  char suffix[3]="a6", file_name[36], fnm[8];

  strcpy(fnm,"sp.");

  if (argc >= 2) strcpy(suffix,argv[1]);
  /* printf("argc=%d argv[0]=%s fnm=%s suffix=%s\n", argc, argv[0], fnm, 
     suffix); */

  fl.in  = fopen(strcat(strcat(strcpy(file_name,fnm),"n."   ),suffix), "r");
  fl.tmp = fopen(strcat(strcat(strcpy(file_name,fnm),"tmp." ),suffix), "w+");
  fl.out = fopen(strcat(strcat(strcpy(file_name,fnm),"out." ),suffix), "w");
  fl.col = fopen(strcat(strcat(strcpy(file_name,fnm),"col." ),suffix), "w");
  fl.trj = fopen(strcat(strcat(strcpy(file_name,fnm),"trj." ),suffix), "w");
  fl.ras = fopen(strcat(strcat(strcpy(file_name,fnm),"ras." ),suffix), "w");

  fscanf(fl.in, "scan=%c\n", &sval.scan_type);

  if (sval.scan_type == 'n')
  {
    /* Simulating one parameter set */
    savr = sm = 0;
    fscanf(fl.in, "seed=%d\n", &sval.seed);
    fprintf(fl.out, " seed=%d\n", sval.seed);
    fl.avr = fl.out;
    update_file_old(&sval, 2, fl);

    ngen = 2;
    length = 17;

    one_par(fl, ran_one, sm, &av);
    write_avr(sval, av, savr, fl);
  }
  else if ((sval.scan_type == 'e') || (sval.scan_type == 'u'))
  {
    /* Simulating several parameter sets */
    savr = sm = 1;
    read_first_input_line(&sval, fl);
    fscanf(fl.in, "seed=%d\n", &sval.seed);
    fprintf(fl.out, "seed=%d\n", sval.seed);
    fl.avr = fopen(strcat(strcat(strcpy(file_name,fnm),"avr." ),suffix), "w");

    for (sval.ipar=0; sval.ipar <= sval.npar; sval.ipar++)
    {
      for (sval.irepeat=1; sval.irepeat <= sval.nrepeat; sval.irepeat++)
      {
        printf("ipar=%d npar=%d irepeat=%d nrepeat=%d\n", sval.ipar, 
        sval.npar, sval.irepeat, sval.nrepeat);
        if (sval.ipar > 0) fprintf(fl.out, "\n");

        update_file_old(&sval, 3, fl);
        one_par(fl, ran_one, sm, &av);
        write_avr(sval, av, savr, fl);
        fflush(fl.avr);
      }
    }

    fclose(fl.avr);
  }
  else if (sval.scan_type == 'f')
  {
    /* Finding borders between regimes of different firing patterns */
    savr = sm = 1;
    read_first_input_lines_f(&sval, fl);

    for (sval.ipara=0; sval.ipara <= sval.npara; sval.ipara++)
    {
      printf("ipara=%d npara=%d para=%lf\n", sval.ipara, sval.npara,
      sval.para_ar[sval.ipara]);
      border_find(&sval, ran_one, sm, &av, fl);
    }
   }
  else if (sval.scan_type == 'p')
  {
    /* Finding borders between regimes of different firing patterns, */
    /* considering the initial conditions and bistability.           */
    savr = sm = 1;
    read_first_input_lines_f(&sval, fl);

    for (sval.ipara=0; sval.ipara <= sval.npara; sval.ipara++)
    {
      printf("ipara=%d npara=%d para=%lf\n", sval.ipara, sval.npara,
      sval.para_ar[sval.ipara]);
      border_find_ic(&sval, ran_one, sm, &av, fl);
    }
   }

  fclose(fl.in);
  fclose(fl.tmp);
  fclose(fl.out);
  fclose(fl.col);
  fclose(fl.trj);
  fclose(fl.ras);
}

/* This function reads and processes the first input line */
void read_first_input_line(scan_val *sval, fl_st fl)
{
  int ipar;
  char line[Mline], *p1, par2all[Mword], *p2, *p3;

  if (fgets(line, Mline, fl.in) == NULL)
  {
    printf("empty input file!!!\n");
    exit(0);
  }
  
  p1 = strtok(line, " ");
  sscanf(p1, " %s", sval->par1);
  p1 = strtok(NULL, " ");
  sscanf(p1, " %s", par2all);

  if (strstr(par2all, "#") == NULL)
  {
    strcpy(sval->par2, par2all);
    sval->npt = 1;
  }
  else
  {
    p3 = &sval->par2[0];
    for (p2 = &par2all[0]; *p2 != '#'; p2++) *p3++ = *p2;
    *p3 = '\0';
    p2++;
   sscanf(p2, "%d", &sval->npt);
  }

  fprintf(fl.out, " par1=%s par2=%s npt=%d", sval->par1, sval->par2,
  sval->npt);

  if (sval->scan_type == 'e')
  {
    p1 = strtok(NULL, " ");
    sscanf(p1, "parmin=%lf", &sval->parmin);
    p1 = strtok(NULL, " ");
    sscanf(p1, "parmax=%lf", &sval->parmax);
    p1 = strtok(NULL, " ");
    sscanf(p1, "npar=%d", &sval->npar);
    p1 = strtok(NULL, " ");
    sscanf(p1, "nrepeat=%d", &sval->nrepeat);
    fprintf(fl.out, " parmin=%lf parmax=%lf\n npar=%d nrepeat=%d\n", 
    sval->parmin, sval->parmax, sval->npar, sval->nrepeat);

    if (sval->npar == 0)
    {
      sval->par_ar[0] = sval->parmin;
    }
    else
    {
      for (ipar=0; ipar<=sval->npar; ipar++)
      {
        sval->par_ar[ipar] = sval->parmin + (sval->parmax - sval->parmin) *
        ipar / sval->npar;
      } 
    }
  }
  else if (sval->scan_type == 'u')
  {
    ipar=-1;
    while((p1 = strtok(NULL, " ")) != NULL)
    {
      sscanf(p1, "%lf", &sval->par_ar[++ipar]);
    }
    sval->npar=ipar;
  }
  else
  {
    printf("wrong scan_type!!!\n");
    exit(0);
  }

  for (ipar=0; ipar<=sval->npar; ipar++)
    fprintf(fl.out, "ipar=%d par_ar=%lf\n", ipar, sval->par_ar[ipar]);

  return;
}

/* This function reads and processes the first input lines for category 'f' */
void read_first_input_lines_f(scan_val *sval, fl_st fl)
{
  char line[Mline];

  if (fgets(line, Mline, fl.in) == NULL)
  {
    printf("empty input file!!!\n");
    exit(0);
  }

  read_parameters_from_line(line, sval->par1a, sval->par2a, &sval->npta,
  &sval->parmina, &sval->parmaxa, &sval->npara, sval->para_ar, fl);

  if (fgets(line, Mline, fl.in) == NULL)
  {
    printf("empty input file!!!\n");
    exit(0);
  }

  read_parameters_from_line(line, sval->par1, sval->par2, &sval->npt,
  &sval->parmin, &sval->parmax, &sval->npar, sval->par_ar, fl);

  if (sval->scan_type == 'f')
  {
    fscanf(fl.in, "eps=%lf\n", &sval->eps);
    fprintf(fl.out, "eps=%lf\n", sval->eps);
  }
  else if (sval->scan_type == 'p')
  {
    fscanf(fl.in, "eps=%lf bhv=%c\n", &sval->eps, &sval->bhv);
    fprintf(fl.out, "eps=%lf bhv=%c\n", sval->eps, sval->bhv);
  }
}

/* This function reads parameters from the first input lines for category  */
/* 'f'                                                                     */
void read_parameters_from_line(char *line, char *par1, char *par2, int *npt,
     double *parmin, double *parmax, int *npar, double *par_ar, fl_st fl)
{
  int ipar;
  char *p1, par2all[Mword], *p2, *p3;

  p1 = strtok(line, " ");
  sscanf(p1, " %s", par1);
  p1 = strtok(NULL, " ");
  sscanf(p1, " %s", par2all);

  if (strstr(par2all, "#") == NULL)
  {
    strcpy(par2, par2all);
    *npt = 1;
  }
  else
  {
    p3 = &par2[0];
    for (p2 = &par2all[0]; *p2 != '#'; p2++) *p3++ = *p2;
    *p3 = '\0';
    p2++;
   sscanf(p2, "%d", npt);
  }

  fprintf(fl.out, "par1=%s par2=%s npt=%d",par1, par2, *npt);

  p1 = strtok(NULL, " ");
  sscanf(p1, "parmin=%lf", parmin);
  p1 = strtok(NULL, " ");
  sscanf(p1, "parmax=%lf", parmax);
  p1 = strtok(NULL, " ");
  sscanf(p1, "npar=%d", npar);
  fprintf(fl.out, " parmin=%lf parmax=%lf npar=%d\n", *parmin, *parmax,
  *npar);

  if (*npar == 0)
  {
    par_ar[0] = *parmin;
  }
  else
  {
    for (ipar=0; ipar<=*npar; ipar++)
    {
      par_ar[ipar] = *parmin + (*parmax - *parmin) * ipar / *npar;
    } 
  }
  for (ipar=0; ipar<=*npar; ipar++)
    fprintf(fl.out, "ipar=%d par_ar=%lf\n", ipar, par_ar[ipar]);
}

/* This function updates the input file and write the new parameter   */
/* value(s)                                                           */
void update_file_old(scan_val *sval, int skip_lines, fl_st fl)
{
  int nget, nchange;
  char line[Mline], iline;

  rewind(fl.in);
  rewind(fl.tmp);
  for (iline=1; iline<=skip_lines; iline++) fgets(line, Mline, fl.in);

  /* no scanning */
  if (sval->scan_type == 'n')
  {
    while (fgets(line, Mline, fl.in) != NULL) fputs(line, fl.tmp);
    rewind(fl.tmp);
    return;
  }

  /* scanning - multiplying all the occurrences of the specific parameter */
  /* value                                                                */

  if (strcmp(sval->par1, "ALL") == 0)
  {
    nchange=0;
    while (nget = (fgets(line, Mline, fl.in)) != NULL)
    {
      if (process_line_old(sval, line, fl.tmp)) nchange++;
    }
    rewind(fl.tmp);
    return;
  }

  /* scanning - changing the specific parameter value */

  while (nget = (fgets(line, Mline, fl.in)) != NULL)
  {
    if (strncmp(line, sval->par1, strlen(sval->par1)) != 0) 
      fputs(line, fl.tmp);
    else
      break;
  }
  fputs(line, fl.tmp);
  
  if (nget == 0) 
  {
    printf("par1 not found!!!\n");
    exit(0);
  }

  while (nget = (fgets(line, Mline, fl.in)) != NULL)
  {
    if (process_line_old(sval, line, fl.tmp))
    {
       break;
    }
  }

  /* checking for end of file */

  if (fgets(line, Mline, fl.in) != NULL) 
    fputs(line, fl.tmp);
  else
  {
    printf("match not found!!!\n");
   exit(0);
  }

  while (fgets(line, Mline, fl.in) != NULL) fputs(line, fl.tmp);
  rewind(fl.tmp);
}

/* This function processes one line and relplace a value by sval.par2 or */
/* multiplies it by sval->par2                                           */
int process_line_old(scan_val *sval, char line[], FILE *ftmp)
{
  int ret_val;

  if (strstr(sval->par2, ":") == NULL)
  {
    /* printf("par2=%s =NU\n", sval->par2); */
    ret_val = process_line_no_colon_old(sval, line, ftmp);
  }
  else
  {
    /* printf("par2=%s !=NU\n", sval->par2); */
    ret_val = process_line_yes_colon_old(sval, line, ftmp);
  }

  return(ret_val);
}

/* This function processes one line and relplace a value by sval.par2 or */
/* multiplies it by sval->par2 .                                         */
/* There is no colon (:) in sval->par2 .                                 */
int process_line_no_colon_old(scan_val *sval, char line[], FILE *ftmp)
{
  double par_ref;
  int  il, il_end, im, ipt, iw;
  int cond1, cond2, cond3;
  char *pline, word[Mword];

  il_end = -1;
  while (line[il_end+1] != '\n' && il_end < Mline-2) il_end++;
  pline = line;
  il = -1;

  while (++il <= il_end)
  {
    /* Condition: matched pattern, '=' at the end, ' ' or beginning of  */
    /* line at the beginning                                            */
    cond1 = strncmp(pline+il, sval->par2, strlen(sval->par2)) == 0;
    cond2 = line[il + strlen(sval->par2)] == '=';
    if (il == 0)
      cond3 = 1;
    else
      cond3 = line[il - 1] == ' ';
    if (cond1 && cond2 && cond3) break;
  }

  if (il >= il_end-1)
  /* par2 does not appear in line */
  {
    fputs(line, ftmp);
    return(0);
  }
  else
  /* par2 appears in line */
  {
    for (im=0; im<il; im++) fputc(line[im], ftmp);
    fprintf(ftmp, "%s=", sval->par2);

    while (line[il-1] != '=') il++;
    ipt=0;
    while (ipt < sval->npt-1)
    {
      putc(line[++il], ftmp);
      if (line[il-1] != ' ' && line[il] == ' ') ipt++;
    }

    while (line[il] == ' ') il++;
    iw=-1;
    while (line[il] != ' ')
    {
      word[++iw] = line[il];
      il++;
    }
    word[++iw] = '\0';

    if (strcmp(sval->par1, "ALL") != 0)
      fprintf(ftmp, "%lf", sval->par_ar[sval->ipar]);
    else
    {
      sscanf(word, "%lf", &par_ref);
      fprintf(ftmp, "%lf", sval->par_ar[sval->ipar] * par_ref);
    }

    for (im=il; im<=il_end; im++) fputc(line[im], ftmp);
    fputc('\n', ftmp);
    return(1);
  }
}

/* This function processes one line and relplace a value by sval.par2 or */
/* multiplies it by sval->par2 .                                         */
/* There is a colon (:) in sval->par2 .                                  */
int process_line_yes_colon_old(scan_val *sval, char line[], FILE *ftmp)
{
  double par_ref;
  int  il, il_end, im, ipt, iw;
  int cond1, cond2, cond3;
  int len2an;
  char line_cp[Mword], par2a[Mword], par2an[Mword], par2b[Mword];
  char *pline, word[Mword], *p1;

  il_end = -1;
  while (line[il_end+1] != '\n' && il_end < Mline-2) il_end++;
  pline = line;

  strcpy(line_cp, sval->par2);
  p1 = strtok(line_cp, ":");
  sscanf(p1, " %s", par2a);
  strcat(strcpy(par2an, par2a),":");
  p1 = strtok(NULL, " ");
  sscanf(p1, " %s", par2b);
  len2an = strlen(par2an);

  /* Checking for par2an */

  if (strncmp(line, par2an, len2an))
  {
    fputs(line, ftmp);
    return(0);
  }

  il = len2an;

  while (++il <= il_end)
  {
    /* Condition: matched pattern, '=' at the end, ' ' or beginning of  */
    /* line at the beginning                                            */
    cond1 = strncmp(pline+il, par2b, strlen(par2b)) == 0;
    cond2 = line[il + strlen(par2b)] == '=';
    if (il == 0)
      cond3 = 1;
    else
      cond3 = line[il - 1] == ' ';
    if (cond1 && cond2 && cond3) break;
  }


  if (il >= il_end-1)
  /* par2 does not appear in line */
  {
    fputs(line, ftmp);
    return(0);
  }
  else
  /* par2 appears in line */
  {
    for (im=0; im<il; im++) fputc(line[im], ftmp);
    fprintf(ftmp, "%s=", par2b);

    while (line[il-1] != '=') il++;
    ipt=0;
    while (ipt < sval->npt-1)
    {
      putc(line[++il], ftmp);
      if (line[il-1] != ' ' && line[il] == ' ') ipt++;
    }

    while (line[il] == ' ') il++;
    iw=-1;
    while (line[il] != ' ')
    {
      word[++iw] = line[il];
      il++;
    }
    word[++iw] = '\0';

    if (strcmp(sval->par1, "ALL") != 0)
      fprintf(ftmp, "%lf", sval->par_ar[sval->ipar]);
    else
    {
      sscanf(word, "%lf", &par_ref);
      fprintf(ftmp, "%lf", sval->par_ar[sval->ipar] * par_ref);
    }

    for (im=il; im<=il_end; im++) fputc(line[im], ftmp);
    fputc('\n', ftmp);
    return(1);
  }
}

/* This function writes population-average results on fl.avr */
void write_avr(scan_val sval, avr_val av, int savr, fl_st fl)
{

  if (savr >= 1) 
  { 
    fprintf(fl.avr, "%9.5lf", sval.par_ar[sval.ipar]);
    fprintf(fl.avr, " %lf %lf %lf %lf %c", av.xmnmx_imtrj, 
    av.n_spk_in_bur_a_all, av.freq_all, av.duty_cycle_all, av.sp_bhv);
    fprintf(fl.avr, " %lf", av.chiE);
    fprintf(fl.avr, "\n");
    fflush(fl.avr);
  }
}

/* This function finds, for a specific value of one parameter para,      */
/* the values of the parameter par for which the firing pattern varies   */
void border_find(scan_val *sval, int *rng_ptr, int sm, avr_val *av, fl_st fl)
{
  double x1, x2;
  int ipar;
  char snl, snr, n1, n2;

  sval->ipar = 0;
  printf("ipar=%d npar=%d par=%lf\n", sval->ipar, sval->npar,
  sval->par_ar[sval->ipar]);

  update_and_run_two_par(sval, fl);
  one_par(fl, rng_ptr, sm, av);
  snl = av->sp_bhv;
  fprintf(fl.out, "ipar=%d snl=%c\n\n", sval->ipar, snl);

  for (ipar=1; ipar<=sval->npar; ipar++)
  {
    sval->ipar = ipar;
    printf("ipar=%d npar=%d par=%lf\n", sval->ipar, sval->npar,
    sval->par_ar[sval->ipar]);

    update_and_run_two_par(sval, fl);
    one_par(fl, rng_ptr, sm, av);
    snr = av->sp_bhv;
    fprintf(fl.out, "ipar=%d snr=%c\n\n", sval->ipar, snr);

    if (snr != snl)
    {
      x1 = sval->par_ar[sval->ipar - 1];
      x2 = sval->par_ar[sval->ipar];
      n1 = snl;
      n2 = snr;
      border_cal(x1, n1, x2, n2, sval, rng_ptr, sm, av, fl);
    }

    snl = snr;
  }
}

/* This function finds the border points recursively */
int border_cal(double x1, char n1, double x2, char n2, scan_val *sval,
    int *rng_ptr, int sm, avr_val *av, fl_st fl)
{
  double x3;
  char n3;

  if (n1 == n2)
  {
    return (0);
  }
  else if (fabs(x1 - x2) < sval->eps)
  {
    fprintf(fl.out, "%lf\n", x1 + 0.5 * (x2 - x1));
    fprintf(fl.ras, "%lf %lf %c %c\n", sval->para_ar[sval->ipara],
    x1 + 0.5 * (x2 - x1), n1, n2);
    fflush(fl.ras);
    printf("%lf\n", x1 + 0.5 * (x2 - x1));
    return (0);
  }
  else
  {
    x3 = x1 + (x2 - x1) / 2.0;

    sval->ipar = sval->npar+1;
    sval->par_ar[sval->ipar] = x3;

    printf("ipar=%d npar=%d par=%lf\n", sval->ipar, sval->npar,
    sval->par_ar[sval->ipar]);

    update_and_run_two_par(sval, fl);
    one_par(fl, rng_ptr, sm, av);
    n3 = av->sp_bhv;

    if (n3 != n1)
    {
      border_cal(x1, n1, x3, n3, sval, rng_ptr, sm, av, fl);
    }

    if (n3 != n2)
    {
      border_cal(x3, n3, x2, n2, sval, rng_ptr, sm, av, fl);
    }
  
    return (1);
  }
}

/* This function modifies two parameters according to sval and then */
/* runs the simulation */
void update_and_run_two_par(scan_val *sval, fl_st fl)
{
  dat_str datstr;
  int idatar;

  read_file_in(&datstr, sval->scan_type, fl);
  update_file(sval->par1a, sval->par2a, sval->npta,
    sval->para_ar[sval->ipara], sval->ipara, &datstr, fl);
  update_file(sval->par1, sval->par2, sval->npt,
    sval->par_ar[sval->ipar], sval->ipar, &datstr, fl);
  for (idatar=1; idatar<=datstr.n_datar; idatar++)
    fprintf(fl.tmp, "%s", datstr.datar[idatar]);
}

/* This function reads the input file into an array of strings        */
void read_file_in(dat_str *datstr, char scan_type, fl_st fl)
{
  char line[Mline];
  int idat, nscan, iscan;

  rewind(fl.in);

  nscan = 1;
  if (scan_type == 'e' || scan_type == 'u') 
    nscan = 2;
  else if (scan_type == 't')
    nscan = 3;
  else if (scan_type == 'f' || scan_type == 'p')
    nscan = 5;

  for (iscan=1; iscan<=nscan; iscan++) fgets(line, Mline, fl.in);

  idat = 0;
  while (fgets(datstr->datar[++idat], Mline, fl.in) != NULL) { }
  datstr->n_datar = idat - 1;
}

/* This function updates the input file and write the new parameter   */
/* value(s)                                                           */
void update_file(char *par1, char *par2, int npt, double par_ar, int ipar,
     dat_str *datstr, fl_st fl)
{
  int nget, nget1, nchange, idatar;

  if (strcmp(par1, "ALL") == 0)
  {
    /* scanning - multiplying all the occurrences of the specific */
    /* parameter value                                            */

    nchange=0;
    for (idatar=1; idatar<=datstr->n_datar; idatar++)
    {
      if (process_line(par1, par2, npt, par_ar, ipar, datstr->datar[idatar],
      fl.tmp)) nchange++;
    }
  }
  else
  {
    /* scanning - changing the specific parameter value */

    nget=0;
    for (idatar=1; idatar<=datstr->n_datar; idatar++)
    {
      nget++;
      if (strncmp(datstr->datar[idatar], par1, strlen(par1)) == 0) 
        break;
    }

    if (nget == datstr->n_datar)
    {
      printf("par1=%s not found!!!\n", par1);
      exit(0);
    }
    nget1 = nget+1;
    for (idatar=nget1; idatar<=datstr->n_datar; idatar++)
    {
      nget++;
      if (process_line(par1, par2, npt, par_ar, ipar, datstr->datar[idatar],
      fl.tmp)) break;
    }

    /* checking for end of file */

    if (nget >= datstr->n_datar) 
    {
      printf("match not found!!!\n");
     exit(0);
    }
  }

  rewind(fl.tmp);
  return;
}

/* This function processes one line and relplace a value by sval.par2 or */
/* multiplies it by sval->par2                                           */
int process_line(char *par1, char *par2, int npt, double par_ar, int ipar,
    char line[], FILE *ftmp)
{
  double par_ref;
  int  il, il_end, im, ipt, iw, len;
  int cond1, cond2, cond3;
  char *pline, word[Mword], newline[Mdatcol], auxline[Mdatcol];
  
  for (il=0; il<=Mdatcol-1; il++) newline[il] = '\0';

  il_end = -1;
  /* 
     while (line[il_end+1] != '\n' && il_end < Mline-2) il_end++; */
  len = strlen(line);
  il_end = len - 2;
  pline = line;
  il = -1;

  while (++il <= il_end)
  {
    /* Condition: matched pattern, '=' at the end, ' ' or beginning of  */
    /* line at the beginning                                            */
    cond1 = strncmp(pline+il, par2, strlen(par2)) == 0;
    cond2 = line[il + strlen(par2)] == '=';
    if (il == 0)
      cond3 = 1;
    else
      cond3 = line[il - 1] == ' ';
    if (cond1 && cond2 && cond3) break;
  }
  
  if (il >= il_end-1)
  /* par2 does not appear in line */
  {
    return(0);
  }
  else
  /* par2 appears in line */
  {
    strncpy(newline, line, il);
    for(im=0; im<strlen(par2); im++)
      newline[il+im] = par2[im];
    newline[il+strlen(par2)] = '=';
    newline[il+strlen(par2)+1] = '\0';

    while (line[il-1] != '=') il++;
    len=strlen(newline);
    ipt=0;
    while (ipt < npt-1)
    {
      sprintf(auxline, "%c", line[il++]);
      strcat(newline, auxline);
      if (line[il-1] != ' ' && line[il] == ' ') 
      {
        ipt++;
        sprintf(auxline, "%c", ' ');
        strcat(newline, auxline);
      }
    }

    while (line[il] == ' ') il++;
    iw=-1;
    while ((line[il] != ' ') && (line[il] != '\n'))
    {
      word[++iw] = line[il];
      il++;
    }

    word[++iw] = '\0';

    if (strcmp(par1, "ALL") != 0)
    {
      sprintf(auxline, "%lf", par_ar);
      strcat(newline, auxline);
    }
    else
    {
      sscanf(word, "%lf", &par_ref);
      sprintf(auxline, "%lf", par_ar * par_ref);
      strcat(newline, auxline);
    }

    len=strlen(newline);
    for (im=il; im<=il_end; im++) 
    {
      if (line[im] != '\n') newline[len+im-il] = line[im];
    }

    /* len=strlen(newline); */
    newline[len+il_end+1-il] = '\n';
    newline[len+il_end+2-il] = '\0';
    strcpy(line, newline);
    return(1);
  }
}

/* This function finds, for a specific value of one parameter para,      */
/* the values of the parameter par for which the firing pattern varies   */
void border_find_ic(scan_val *sval, int *rng_ptr, int sm, avr_val *av,
     fl_st fl)
{
  double x1, x2;
  double Vstore1[Meq], Vstore2[Meq], Vstore3[Meq];
  int ipar, ieq;
  char snl, snr, n1, n2;

  sval->ipar = 0;
  printf("ipar=%d npar=%d par=%lf\n", sval->ipar, sval->npar,
  sval->par_ar[sval->ipar]);

  update_and_run_two_par(sval, fl);
  vary_I_one_par(fl, rng_ptr, sm, Vstore1, av);
  snl = av->sp_bhv;
  fprintf(fl.out, "ipar=%d snl=%c\n\n", sval->ipar, snl);
  fprintf(fl.out, "Vstore1=\n");
  for (ieq=1; ieq<=9; ieq++) fprintf(fl.out, " %lf\n", Vstore1[ieq]);
  fprintf(fl.out, "\n");

  for (ipar=1; ipar<=sval->npar; ipar++)
  {
    sval->ipar = ipar;
    printf("ipar=%d npar=%d par=%lf\n", sval->ipar, sval->npar,
    sval->par_ar[sval->ipar]);

    update_and_run_two_par(sval, fl);
    vary_I_one_par(fl, rng_ptr, sm, Vstore2, av);
    snr = av->sp_bhv;
    fprintf(fl.out, "ipar=%d snr=%c\n\n", sval->ipar, snr);
    fprintf(fl.out, "Vstore2=\n");
    for (ieq=1; ieq<=9; ieq++) fprintf(fl.out, " %lf\n", Vstore2[ieq]);
    fprintf(fl.out, " %lf\n");
    for (ieq=1; ieq<=Meq; ieq++) Vstore3[ieq] = Vstore2[ieq];

    fprintf(fl.out, "p1=%lf bhv1=%c p2=%lf bhv2=%c\n",
      sval->par_ar[sval->ipar-1], snl, sval->par_ar[sval->ipar], snr);

    if ((snl == sval->bhv) && (snr!= sval->bhv))
    {
      locate_border_behavior(sval->ipar-1, Vstore1, Vstore2, sval->ipar, sval,
      rng_ptr, sm, av, fl);
    }
    else if ((snl != sval->bhv) && (snr== sval->bhv))
    {
      locate_border_behavior(sval->ipar, Vstore2, Vstore1, sval->ipar-1, sval,
      rng_ptr,  sm, av, fl);
    }
    else
    {
      printf("same bhv snr=%c snl=%c\n", snr, snl);
    }

    for (ieq=1; ieq<=Meq; ieq++) Vstore1[ieq] = Vstore3[ieq];
    snl = snr;
  }
}

/* This function finds the border of a specific behavior type */
void locate_border_behavior(int ipar_yes, double *Vstore_in,
     double *Vstore_out, int ipar_no, scan_val *sval, int *rng_ptr, int sm,
     avr_val *av, fl_st fl)
{
  double x1, x2, x3, deltax, epsilon, *Vs1, *Vs2, *Vs3, sign;
  char n1, n3;
  int ieq, nstep;

  epsilon = 1.0e-10;

  fprintf(fl.out, "\nipy=%d py=%lf ipn=%d pn=%lf\n", ipar_yes,
  sval->par_ar[ipar_yes], ipar_no, sval->par_ar[ipar_no]);
  fprintf(fl.out, "Vstore_in=\n");
  for (ieq=1; ieq<=9; ieq++) fprintf(fl.out, " %lf\n", Vstore_in[ieq]);
  fprintf(fl.out, " %lf\n");

  x1 = sval->par_ar[ipar_yes];
  n3 =n1 = sval->bhv;
  x2 = sval->par_ar[ipar_no];
  Vs1 = Vstore_in;
  Vs2 = Vstore_out;

  deltax = (x2 - x1);

  if (deltax > 0)
    sign = 1.0;
  else
    sign = -1.0;

  while (fabs(deltax) > sval->eps / 2.0)
  {
    deltax /= 2.0;
    nstep = 1;

    /* advancing the parameters to continue the desired state */
    while ( ( sign * (x3 = x1 + deltax) <= sign * (x2 - epsilon) ) &&
            ((nstep == 1) || (n3 == sval->bhv)) )
    {
      printf("in loop x1=%lf x2=%lf x3=%lf deltax=%lf nstep=%d\n", x1, x2, x3, 
      deltax, nstep);
      printf("V=%lf", Vs1[1]);
      for (ieq=2; ieq<=9; ieq++) printf(" %lf", Vs1[ieq]);
      printf("\n");
      fprintf(fl.out, "V=%20.14lf", Vs1[1]);
      for (ieq=2; ieq<=9; ieq++) fprintf(fl.out, " %17.14lf", Vs1[ieq]);
      fprintf(fl.out, "\n");
      sval->ipar = sval->npar+1;
      sval->par_ar[sval->ipar] = x3;
      update_and_run_two_par(sval, fl);
      one_par_Vstore(fl, rng_ptr, sm, Vs1, Vs2, av);
      n3 = av->sp_bhv;
      printf("in loop n3=%c\n", n3);
      if (n3 == sval->bhv)
      {
        x1 = x3;
        Vs3 = Vs1;
        Vs1 = Vs2;
        Vs2 = Vs3;
      }
      nstep++;
    }
  }
    
  if (sign * x3 <= sign * (x2 - epsilon))
  {
    fprintf(fl.out, "border: x3=%lf\n", x3);
  }

  fprintf(fl.out, "x1=%lf\n", x1);
  printf("x1=%lf\n", x1);

  if (deltax > 0.0)
    fprintf(fl.ras, "%lf %lf %c %c\n", sval->para_ar[sval->ipara],
    x1, sval->bhv, n3);
  else
    fprintf(fl.ras, "%lf %lf %c %c\n", sval->para_ar[sval->ipara],
    x1, n3, sval->bhv);

    fflush(fl.ras);
}
