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

  int max_degree = 9;

  const int coeff_size = poly_system_get_coeffs(&poly, max_degree, 0);
  float *coeff = (float *)malloc(sizeof(float)*coeff_size);

  const int sample_cnt = 100000;
  float *sample_in = (float *)malloc(sample_cnt*sizeof(float)*5);
  const int oversample = 10; // only do this x coeff count many ray tracing samples

  // only rays that pass the apperture and reach the outer pupil are used for
  // calculating the polynomials error
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
    int error = evaluate(lenses, lenses_cnt, zoom, ray_in, out, aspheric_elements);
    if(!error)
    {
      for(int k=0;k<5;k++) sample_in[5*valid + k] = ray_in[k];
      valid++;
    }
    // only need to be able to determine the dimensionality of our problem, not much more:
    if(valid > max(1000, oversample*coeff_size)) break;
  }

  // calculate reference evaluation of polynomial
  fprintf(stderr, "[ sensor->outer pp ] optimising %d coeffs by %d/%d valid sample points\n", coeff_size, valid, sample_cnt);
  float *sample_ref = (float *)malloc(valid*sizeof(float)*5);
  for(int i=0;i<valid;i++)
    poly_system_evaluate(&poly, sample_in+5*i, sample_ref+5*i, max_degree);

  // get (fitted) coefficients of original polynomial
  poly_system_get_coeffs(&poly, max_degree, coeff);
  float *err = (float *)malloc(coeff_size*sizeof(float));
  for(int i=0;i<coeff_size;i++)
  {
    //set one of the coefficients to zero at a time
    float curr_coeff = coeff[i];
    coeff[i] = 0;
    poly_system_set_coeffs(&poly, max_degree, coeff);

    // evaluate the same rays as above and calculate the error introduced by
    // removing term i
    err[i] = 0;
    for(int j=0;j<valid;j++)
    {
      float tmp[5];
      poly_system_evaluate(&poly, sample_in+5*j, tmp, max_degree);
      for(int k=0;k<5;k++)
        err[i] += (tmp[k]-sample_ref[5*j+k])*(tmp[k]-sample_ref[5*j+k]);
    }

    //restore coefficient and iterate
    coeff[i] = curr_coeff;
    poly_system_set_coeffs(&poly, max_degree, coeff);
  }

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

  //Do the same for the aperture
  const int ap_coeff_size = poly_system_get_coeffs(&poly_ap, max_degree, 0);
  err = (float *)malloc(ap_coeff_size*sizeof(float));
  coeff = (float *)malloc(ap_coeff_size*sizeof(float));
  fprintf(stderr, "[ sensor->aperture ] optimising %d coeffs by %d/%d valid sample points\n", ap_coeff_size, valid, sample_cnt);
  for(int i=0;i<valid;i++)
    poly_system_evaluate(&poly_ap, sample_in+5*i, sample_ref+5*i, max_degree);

  poly_system_get_coeffs(&poly_ap, max_degree, coeff);
  for(int i=0;i<ap_coeff_size;i++)
  {
    //set one of the coefficients to zero at a time
    float curr_coeff = coeff[i];
    coeff[i] = 0;
    poly_system_set_coeffs(&poly_ap, max_degree, coeff);

    //evaluate rays and calculate error
    err[i] = 0;
    for(int j=0;j<valid;j++)
    {
      float tmp[5];
      poly_system_evaluate(&poly_ap, sample_in+5*j, tmp, max_degree);
      for(int k=0;k<5;k++)
        err[i] += (tmp[k]-sample_ref[5*j+k])*(tmp[k]-sample_ref[5*j+k]);
    }

    //restore coefficient
    coeff[i] = curr_coeff;
    poly_system_set_coeffs(&poly_ap, max_degree, coeff);
  }
  //poly_print(&poly_ap, 0, stderr);
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
  free(sample_ref);

  exit(0);
}
