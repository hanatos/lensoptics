#include "raytrace.h"
#include "poly.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>

static lens_element_t lenses[50];
static int lenses_cnt = 0;
static float zoom = 0;
static int aspheric_elements = 1;

int main(int argc, char *arg[])
{
  if(argc < 2)
  {
    fprintf(stderr, "usage: %s lensfile\n", arg[0]);
    exit(1);
  }

  char *lensfilename = arg[1];

  poly_system_t poly, poly_ap;
  char fitf[1024], afitf[1024];
  snprintf(fitf,  sizeof(fitf),  "%s.fit", lensfilename);
  snprintf(afitf, sizeof(afitf), "%s_ap.fit", lensfilename);
  if(poly_system_read(&poly, fitf) || poly_system_read(&poly_ap, afitf))
  {
    fprintf(stderr, "[simplify] could not read poly fits for `%s'!\n", lensfilename);
    exit(1);
  }

  lenses_cnt = lens_configuration(lenses, lensfilename, sizeof(lenses));
  const float p_dist = lens_get_thickness(lenses + lenses_cnt-1, zoom);
  const float p_rad = lenses[lenses_cnt-1].housing_radius;

  int max_degree = 14;

  const int coeff_size = poly_system_get_coeffs(&poly, max_degree, 0);
  float *coeff = (float *)malloc(sizeof(float)*coeff_size);

  const int sample_cnt = 100000;
  float *sample_in = (float *)malloc(sample_cnt*sizeof(float)*5);
  const int oversample = 10; // only do this x coeff count many ray tracing samples

  // only rays that pass the apperture and reach the outer pupil are used for
  // calculating the polynomials error
  int valid = 0;
  while(1)
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
    int error = evaluate(lenses, lenses_cnt, zoom, ray_in, out, aspheric_elements);
    if(!error)
    {
      for(int k=0;k<5;k++) sample_in[5*valid + k] = ray_in[k];
      valid++;
    }
    // only need to be able to determine the dimensionality of our problem, not much more:
    if(valid > max(1000, oversample*coeff_size)) break;
  }

  // calculate error introduced by removing a term
  float *err = (float *)malloc(coeff_size*sizeof(float));
  int sum_terms = 0;

  for(int variable = 0; variable < poly_num_vars; variable++)
  {
    for(int term = 0; term < poly.poly[variable].num_terms; term++, sum_terms++)
    {
      err[sum_terms] = 0;
      poly_term_t t = poly.poly[variable].term[term];
      for(int j = 0; j < valid; j++)
      {
        float error = poly_term_evaluate(&t, sample_in+5*j);
        err[sum_terms] += error*error;
      }
    }
  }

  poly_system_get_coeffs(&poly, max_degree, coeff);
  // remove terms with error < eps
  for(int i=0;i<coeff_size;i++)
  {
    if(err[i]/valid < 1e-6)
      coeff[i] = 0;
  }
  poly_system_set_coeffs(&poly, max_degree, coeff);

  // write optimised poly
  char fitfile[2048];
  snprintf(fitfile, 2048, "%s.fit", lensfilename);
  poly_system_simplify(&poly);
  fprintf(stderr, "output poly has %d coeffs.\n", poly_system_get_coeffs(&poly, max_degree, 0));
  //poly_print(&poly, 0, stderr);
  poly_system_write(&poly, fitfile);
  free(err);
  free(coeff);

  //Do the same for the aperture polynomial
  const int ap_coeff_size = poly_system_get_coeffs(&poly_ap, max_degree, 0);
  coeff = (float *)malloc(sizeof(float)*ap_coeff_size);
  // calculate error introduced by removing a term
  err = (float *)malloc(ap_coeff_size*sizeof(float));
  sum_terms = 0;

  for(int variable = 0; variable < poly_num_vars; variable++)
  {
    for(int term = 0; term < poly_ap.poly[variable].num_terms; term++, sum_terms++)
    {
      err[sum_terms] = 0;
      poly_term_t t = poly_ap.poly[variable].term[term];
      for(int j = 0; j < valid; j++)
      {
        float error = poly_term_evaluate(&t, sample_in+5*j);
        err[sum_terms] += error*error;
      }
    }
  }

  poly_system_get_coeffs(&poly_ap, max_degree, coeff);
  // remove terms with error < eps
  for(int i=0;i<ap_coeff_size;i++)
  {
    if(err[i]/valid < 1e-6)
      coeff[i] = 0;
  }
  poly_system_set_coeffs(&poly_ap, max_degree, coeff);

  // write optimised poly
  snprintf(fitfile, 2048, "%s_ap.fit", lensfilename);
  poly_system_simplify(&poly_ap);
  fprintf(stderr, "output poly has %d coeffs.\n", poly_system_get_coeffs(&poly_ap, max_degree, 0));
  //poly_print(&poly_ap, 0, stderr);
  poly_system_write(&poly_ap, fitfile);

  free(err);
  free(coeff);
  free(sample_in);

  exit(0);
}
