#include "raytrace.h"
#include "levmar.h"
#include "poly.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>
# define M_PI   3.14159265358979323846  /* pi */


static lens_element_t lenses[50];
static int lenses_cnt = 0;
static const float zoom = 0.0f;
static int max_degree = 4;

static inline float ap(float x, int n, float p_dist, float p_rad, int dim)
{ // sample incoming pupil
  return 2.0f*(x/(n-1.0f)-.5f); // just return [-1,1] arbitrarily
}

typedef struct tmpopt_t
{
  float *sample_in;
  poly_system_t *poly;
  float *reference;
  float *last_valid_param;
  int fit_idx;
}
tmpopt_t;

void print_coeffs(float *param, int param_cnt)
{
  int row = 8;
  for(int p=0;p<param_cnt/row;p++)
  {
    fprintf(stderr, "p[%d] = ", p*row);
    for(int q=0;q<row;q++)
      if(p*row+q < param_cnt)
        if(param[p] != 0.0f) fprintf(stderr, "%f ", param[p]);
    fprintf(stderr, "\n");
  }
}
void eval_poly(float *param, float *sample, int param_cnt, int sample_cnt, void *data)
{
  tmpopt_t *tmp = (tmpopt_t *)data;

  // recover broken parameters
  for(int k=0;k<param_cnt;k++)
  {
    if(param[k] == param[k]) tmp->last_valid_param[k] = param[k];
    else param[k] = tmp->last_valid_param[k];
  }

  // fill params into poly
  //poly_system_set_coeffs(tmp->poly, max_degree, param);
  poly_set_coeffs(tmp->poly->poly + tmp->fit_idx, max_degree, param);

#pragma omp parallel for default(shared)
  for(int s=0;s<sample_cnt;s++)
  {
    const float *in = tmp->sample_in + 5*s;
    float out[5];
    poly_system_evaluate(tmp->poly, in, out, max_degree);
    // system is really 5x4, don't output wavelength:
    int k = tmp->fit_idx;
    //for(int k=0;k<4;k++)
    {
      sample[s] = out[k];
      assert(sample[s] == sample[s]);
    }
  }
}

// our 'variables' are the coefficients of the polynomial. The derivatives hence,
// equal the evaluation of the terms with coefficient = 1.
void eval_jac(float *param, float *j, int param_cnt, int sample_cnt, void *data)
{
  tmpopt_t *tmp = (tmpopt_t *)data;

  //no need to set coefficients - we set them to one anyway

  int idx = 0;
  for(int s=0;s<sample_cnt;s++)
  {
    const float *in = tmp->sample_in + 5*s;
    for(int i=0;i<param_cnt;i++)
    {
      poly_term_t term = tmp->poly->poly[tmp->fit_idx].term[i];
      term.coeff = 1;
      j[idx++] = poly_term_evaluate(&term, in);
      assert(j[idx-1] == j[idx-1]);
    }
  }
}

int main(int argc, char *arg[])
{
  if(argc < 2)
  {
    fprintf(stderr, "usage: %s lensfile\n", arg[0]);
    exit(1);
  }
  char *lensfilename = arg[1];
  lenses_cnt = lens_configuration(lenses, lensfilename, sizeof(lenses));
  const float p_dist = lens_get_thickness(lenses + lenses_cnt-1, zoom);
  const float p_rad = lenses[lenses_cnt-1].housing_radius;

  int user_degree = 4;
  max_degree = 9;
  if(argc > 2) user_degree = atol(arg[2]);
  if(user_degree < 1) user_degree = 1;
  if(user_degree > 9) user_degree = 9;

  int min_degree = 3;
  if(argc > 3) min_degree = atol(arg[3]);
  if(min_degree < 1) min_degree = 1;
  if(min_degree > user_degree) min_degree = user_degree;

  int min_degree_aperture = 1;
  if(argc > 4) min_degree_aperture = atol(arg[4]);
  if(min_degree_aperture < 1) min_degree_aperture = 1;
  if(min_degree_aperture > user_degree) min_degree_aperture = user_degree;

  int pass2 = 0;
  for(int i = 1; i < argc && pass2 == 0; i++)
    if(atol(arg[i]) == -2)
      pass2 = 1;

  char fitfile[2048], apfitfile[2048];
  snprintf(fitfile, 2048, "%s.fit", lensfilename);
  snprintf(apfitfile, 2048, "%s_ap.fit", lensfilename);

  // load generic 1233-coefficient degree 9 polynomial template with all zero coeffs:
  poly_system_t poly, poly_ap;
  if(!pass2)
  {
    if(poly_system_read(&poly, "degree9-sorted.poly") || poly_system_read(&poly_ap, "degree9-sorted.poly"))
    {
      fprintf(stderr, "[fit] could not read `degree9.poly' template!\n");
      exit(1);
    }
  }
  else
  {
    if(!poly_system_read(&poly, fitfile) && !poly_system_read(&poly_ap, apfitfile))
      user_degree = min_degree = min_degree_aperture = max_degree = 9;
    else
    {
      fprintf(stderr, "[fit] could not read fitted polynomials %s and %s!\n", fitfile, apfitfile);
      exit(1);
    }
  }
  // already storing sorted:
  // poly_system_sort(&poly);
  // poly_system_sort(&poly_ap);

  const int coeff_size = max(poly_system_get_coeffs(&poly, max_degree, 0),
    poly_system_get_coeffs(&poly_ap, max_degree, 0));
  float *coeff = (float *)malloc(sizeof(float)*coeff_size);

  const int sample_cnt = 100000;
  float *sample = (float *)malloc(sample_cnt*sizeof(float)*5);
  float *sample_in = (float *)malloc(sample_cnt*sizeof(float)*5);
  const int oversample = 10; // only do this x coeff count many ray tracing samples

  int valid = 0;
  for(int i=0;i<sample_cnt;i++)
  {
    const float u = drand48(), v = drand48(), w = drand48(), x = drand48(), y = drand48();
    float ray_in[] = {
      p_rad * 4.0f * (x-0.5),
      p_rad * 4.0f * (y-0.5),
      p_rad/p_dist * cosf(2.0f*M_PI*u)*sqrtf(v),
      p_rad/p_dist * sinf(2.0f*M_PI*u)*sqrtf(v),
      0.4 + 0.3*w};
    ray_in[2] -= ray_in[0] / p_dist;
    ray_in[3] -= ray_in[1] / p_dist;
    float out[5];
    int error = evaluate(lenses, lenses_cnt, zoom, ray_in, out);
    if(!error)
    {
      for(int k=0;k<5;k++) sample_in[5*valid + k] = ray_in[k];
      for(int k=0;k<5;k++) sample[valid+k*sample_cnt] = out[k];
      valid++;
    }
    // only need to be able to determine the dimensionality of our problem, not much more:
    if(valid > oversample*coeff_size) break;
  }
  fprintf(stderr, "[ sensor->outer pp ] optimising %d coeffs by %d/%d valid sample points\n", coeff_size, valid, sample_cnt);

  // levmar setup
  tmpopt_t tmp;
  tmp.sample_in = sample_in;
  float opts[LM_OPTS_SZ], info[LM_INFO_SZ];
  // opts[0]=LM_INIT_MU; opts[1]=1E-3; opts[2]=1E-5; opts[3]=1E-7; // terminates way to early
  opts[0]=LM_INIT_MU; opts[1]=1E-7; opts[2]=1E-7; opts[3]=1E-12; // known to go through
  /// opts[0]=1e-2f; opts[1]=1E-7; opts[2]=1E-7; opts[3]=1E-12; // known to go through
  // opts[0]=LM_INIT_MU; opts[1]=1E-7; opts[2]=1E-8; opts[3]=1E-13; // goes through, some nans
  // opts[0] = LM_INIT_MU; opts[1]=1E-8; opts[2]=1E-8; opts[3]=1E-15;
  opts[4] = LM_DIFF_DELTA;
  tmp.reference = sample;
  tmp.poly = &poly;
  tmp.last_valid_param = (float *)malloc(sizeof(float)*coeff_size);

  // ===================================================================================================
  // evaluate poly sensor -> outer pupil
  // ===================================================================================================
  float last_error[5] = {FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX};
  poly_system_t poly_backup = {0};

  poly_system_copy(&poly, &poly_backup);
  for(max_degree=min_degree;max_degree <= user_degree; max_degree++)
  {
    int sumCoeffs = 0;
    float errorSum = 0.0f;
    memset(tmp.last_valid_param, 0, sizeof(float)*coeff_size);
    for(int j = 0; j < 5; j++)
    {
      tmp.fit_idx = j;
      //const int degree_coeff_size = poly_system_get_coeffs(&poly, max_degree, 0);
      const int degree_coeff_size = poly_get_coeffs(poly.poly + j, max_degree, 0);

      //restore coefficients from backup poly as initial guess
      memset(coeff + sumCoeffs, 0, sizeof(float)*degree_coeff_size);
      //poly_get_coeffs(poly_backup.poly+j, max_degree, coeff + sumCoeffs);

      // optimize taylor polynomial a bit
      //slevmar_dif(eval_poly, coeff + sumCoeffs, sample+j*sample_cnt, degree_coeff_size, valid, 1000, opts, info, NULL, NULL, &tmp);
      slevmar_der(eval_poly, eval_jac, coeff + sumCoeffs, sample+j*sample_cnt, degree_coeff_size, valid, 1000, opts, info, NULL, NULL, &tmp);
      sumCoeffs += degree_coeff_size;
      if(info[1] < last_error[j])
      {
        last_error[j] = info[1];
        poly_destroy(poly_backup.poly+j);
        poly_copy(poly.poly+j, poly_backup.poly+j);
      }
//      else
//      {
//        poly_destroy(poly.poly+j);
//        poly_copy(poly_backup.poly+j, poly.poly+j);
//        poly_get_coeffs(poly_backup.poly+j, max_degree, coeff+sumCoeffs);
//      }
      fprintf(stderr, "error: %g\n", info[1]);
      errorSum += max(0.0f, info[1]);
    }
    fprintf(stderr, "degree %d has %d samples, fitting error %g\n", max_degree, sumCoeffs, errorSum);
  }

  // write optimised poly
  poly_system_simplify(&poly_backup);
  fprintf(stderr, "output poly has %d coeffs.\n", poly_system_get_coeffs(&poly_backup, user_degree, 0));
  poly_system_write(&poly_backup, fitfile);

  // ===================================================================================================
  // evaluate_aperture poly sensor -> aperture
  // ===================================================================================================
  memset(sample, 0, sizeof(float)*sample_cnt*5);
  for(int i=0;i<valid;i++)
  {
    float *ray_in = sample_in+5*i;
    float out[5] = {0.0f};
    int error = evaluate_aperture(lenses, lenses_cnt, zoom, ray_in, out);
    (void)error; // silence non-debug build warning
    assert(error == 0);
    for(int k=0;k<5;k++)
      sample[i+k*sample_cnt] = out[k];
  }
  fprintf(stderr, "[ sensor->aperture ] optimising %d coeffs by %d/%d valid sample points\n", coeff_size, valid, sample_cnt);

  tmp.poly = &poly_ap;

  for(int i = 0; i < 5; i++) last_error[i] = FLT_MAX;
  poly_system_copy(&poly_ap, &poly_backup);

  for(max_degree=min_degree_aperture;max_degree <= user_degree; max_degree++)
  {
    int sumCoeffs = 0;
    float errorSum = 0.0f;
    memset(tmp.last_valid_param, 0, sizeof(float)*coeff_size);
    for(int j = 0; j < 5; j++)
    {
      tmp.fit_idx = j;
      //const int degree_coeff_size = poly_system_get_coeffs(&poly, max_degree, 0);
      const int degree_coeff_size = poly_get_coeffs(poly_ap.poly + j, max_degree, 0);

      //restore coefficients from backup poly as initial guess
      memset(coeff + sumCoeffs, 0, sizeof(float)*degree_coeff_size);
      //poly_get_coeffs(poly_backup.poly+j, max_degree, coeff + sumCoeffs);

      // optimize taylor polynomial a bit
      //slevmar_dif(eval_poly, coeff + sumCoeffs, sample+j*sample_cnt, degree_coeff_size, valid, 1000, opts, info, NULL, NULL, &tmp);
      slevmar_der(eval_poly, eval_jac, coeff + sumCoeffs, sample+j*sample_cnt, degree_coeff_size, valid, 1000, opts, info, NULL, NULL, &tmp);
      sumCoeffs += degree_coeff_size;
      if(info[1] < last_error[j])
      {
        last_error[j] = info[1];
        poly_destroy(poly_backup.poly+j);
        poly_copy(poly_ap.poly+j, poly_backup.poly+j);
      }
//      else
//      {
//        poly_destroy(poly_ap.poly+j);
//        poly_copy(poly_backup.poly+j, poly_ap.poly+j);
//        poly_get_coeffs(poly_backup.poly+j, max_degree, coeff+sumCoeffs);
//      }
      fprintf(stderr, "error: %g\n", info[1]);
      errorSum += max(0.0f, info[1]);
    }
    fprintf(stderr, "degree %d has %d samples, fitting error %g\n", max_degree, sumCoeffs, errorSum);
  }

  poly_system_simplify(&poly_backup);
  // TODO: this totally doesn't throw away useless coeffs, the fitter will make ineffective ones != 0!
  fprintf(stderr, "output aperture poly has %d coeffs.\n", poly_system_get_coeffs(&poly_backup, user_degree, 0));
  poly_system_write(&poly_backup, apfitfile);
  exit(0);
}
