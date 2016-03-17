#include <Eigen/Dense>
#include "raytrace.h"
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
static int aspheric_elements = 1;
static const float precision[5] = {1e-6, 1e-6, 1e-8, 1e-8, 1e-9};

static inline float ap(float x, int n, float p_dist, float p_rad, int dim)
{ // sample incoming pupil
  return 2.0f*(x/(n-1.0f)-.5f); // just return [-1,1] arbitrarily
}

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
  max_degree = 15;
  if(argc > 2) user_degree = atol(arg[2]);
  if(user_degree < 1) user_degree = 1;
  if(user_degree > 15) user_degree = 15;

  int min_degree = 3;
  if(argc > 3) min_degree = atol(arg[3]);
  if(min_degree < 1) min_degree = 1;
  if(min_degree > user_degree) min_degree = user_degree;

  int min_degree_aperture = 1;
  if(argc > 4) min_degree_aperture = atol(arg[4]);
  if(min_degree_aperture < 1) min_degree_aperture = 1;
  if(min_degree_aperture > user_degree) min_degree_aperture = user_degree;

  int max_coeffs = 10000;
  if(argc > 5) max_coeffs = atol(arg[5]);

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
    if(poly_system_read(&poly, "sorted.poly") || poly_system_read(&poly_ap, "sorted.poly"))
    {
      fprintf(stderr, "[fit] could not read `sorted.poly' template!\n");
      exit(1);
    }
  }
  else
  {
    if(!poly_system_read(&poly, fitfile) && !poly_system_read(&poly_ap, apfitfile))
      user_degree = min_degree = min_degree_aperture = max_degree = 15;
    else
    {
      fprintf(stderr, "[fit] could not read fitted polynomials %s and %s!\n", fitfile, apfitfile);
      exit(1);
    }
  }
  // already storing sorted:
  // poly_system_sort(&poly);
  // poly_system_sort(&poly_ap);

  const int coeff_size = max(poly_system_get_coeffs(&poly, user_degree, 0),
    poly_system_get_coeffs(&poly_ap, user_degree, 0));
  float *coeff = (float *)malloc(sizeof(float)*coeff_size);

  const int sample_cnt = 1000000;
  float *sample = (float *)malloc(sample_cnt*sizeof(float)*5);
  float *sample_in = (float *)malloc(sample_cnt*sizeof(float)*5);
  const int oversample = 4; // only do this x coeff count many ray tracing samples
//  #define BUCKET
  #ifdef BUCKET
  const int bucket_cnt = 10;
  int bucket[bucket_cnt];
  for(int i = 0; i < bucket_cnt; i++) bucket[i] = 0;
  #endif

  int valid = 0;
  while(1)
  {
    const float u = drand48(), v = drand48(), w = drand48(), x = drand48(), y = drand48();
    float ray_in[] = {
      // p_rad * 4.0f * (x-0.5f),
      // p_rad * 4.0f * (y-0.5f),
      35.0f/2.0f - x*35.0f, // 35mm film, isotropic
      35.0f/2.0f - y*35.0f,
      p_rad/p_dist * cosf(2.0f*M_PI*u)*sqrtf(v),
      p_rad/p_dist * sinf(2.0f*M_PI*u)*sqrtf(v),
      0.4f + 0.3f*w};
    ray_in[2] -= ray_in[0] / p_dist;
    ray_in[3] -= ray_in[1] / p_dist;
    float out[5];
    int error = evaluate(lenses, lenses_cnt, zoom, ray_in, out, aspheric_elements);
    #ifdef BUCKET
    int bucket_num = bucket_cnt * (out[0]*out[0]+out[1]*out[1]) / (lenses[0].housing_radius*lenses[0].housing_radius);
    if(!error && bucket[bucket_num] <= oversample*coeff_size/bucket_cnt)
    #else
    if(!error)
    #endif
    {
      for(int k=0;k<5;k++) sample_in[5*valid + k] = ray_in[k];
      for(int k=0;k<5;k++) sample[valid+k*sample_cnt] = out[k];
      valid++;
      #ifdef BUCKET
      bucket[bucket_num]++;
      #endif
    }
    // only need to be able to determine the dimensionality of our problem, not much more:
    if(valid >= oversample*coeff_size) break;
    if(valid >= sample_cnt) break;
  }
  fprintf(stderr, "[ sensor->outer pp ] optimising %d coeffs by %d/%d valid sample points\n", coeff_size, valid, sample_cnt);
  const char *outvname[] = {"x", "y", "dx", "dy", "transmittance"};

  // ===================================================================================================
  // evaluate poly sensor -> outer pupil
  // ===================================================================================================
  float last_error[5] = {FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX};
  poly_system_t poly_backup = {0};
  poly_system_copy(&poly, &poly_backup);

  for(max_degree=min_degree;max_degree <= user_degree; max_degree++)
  {
    int sumCoeffs = 0;
    int maxSumCoeffs = 0;
    float errorSum = 0.0f;
    for(int j = 0; j < 5; j++)
    {
      //const int degree_coeff_size = poly_system_get_coeffs(&poly, max_degree, 0);
      const int degree_coeff_size = poly_get_coeffs(poly.poly + j, max_degree, 0);
      int coeff_cnt = degree_coeff_size;
      const int degree_num_samples = std::min(valid, oversample*degree_coeff_size);

      // optimize taylor polynomial a bit
      Eigen::MatrixXd A(degree_coeff_size, degree_num_samples);
      #pragma omp parallel for
      for(int y = 0; y < degree_num_samples; y++)
      {
        for(int x = 0; x < degree_coeff_size; x++)
        {
          poly_term_t term = poly.poly[j].term[x];
          term.coeff = 1;
          A(x,y) = poly_term_evaluate(&term, sample_in+5*y);
        }
      }
      //cout<<A<<endl;
      Eigen::VectorXd b(degree_num_samples);
      for(int y = 0; y < degree_num_samples; y++)
        b(y) = sample[j*sample_cnt+y];
      //cout<<b<<endl;

      Eigen::VectorXd result = Eigen::ArrayXd::Zero(degree_coeff_size);
      Eigen::VectorXd residual = b;
      Eigen::VectorXd factor(degree_coeff_size);
      Eigen::VectorXd used(degree_coeff_size);
      int permutation[degree_coeff_size];
      for(int i = 0; i < degree_coeff_size; i++) permutation[i] = i;
      for(int i = 0; i < degree_coeff_size; i++) used[i] = 0.0;
      for(int i = 0; i < degree_coeff_size; i++)
        factor(i) = 1 / (A.row(i).norm() * (poly_term_get_degree(poly.poly[j].term+i)+1.0));

      coeff_cnt = 0;
      for(int i = 0; i < degree_coeff_size; i++)
      {
        int maxidx = 0;
        Eigen::VectorXd prod = (Eigen::ArrayXd(A * residual) * Eigen::ArrayXd(factor) * (1-Eigen::ArrayXd(used))).abs();
        prod.maxCoeff(&maxidx);
        permutation[coeff_cnt] = maxidx;
        coeff_cnt++;

        Eigen::MatrixXd tmp2(degree_num_samples, coeff_cnt);
        for(int k = 0; k < coeff_cnt; k++)
          tmp2.col(k) = A.row(permutation[k]).transpose();
        result = (tmp2.transpose()*tmp2).ldlt().solve(tmp2.transpose()*b);
        residual = b-tmp2*result;
        used(maxidx) = 1.0;
        if(j < 4) // don't limit transmittance precision
        if(residual.squaredNorm() < precision[j] * degree_num_samples)
          break;
        if(i > max_coeffs) break; // force sparsity
      }

      Eigen::VectorXd coeffs = Eigen::ArrayXd::Zero(degree_coeff_size);
      for(int k = 0; k < coeff_cnt; k++)
        coeffs(permutation[k]) = result(k);
      result = coeffs;
      float error = residual.squaredNorm() / degree_num_samples;

      sumCoeffs += coeff_cnt;
      maxSumCoeffs += degree_coeff_size;
      if(error < last_error[j])
      {
        last_error[j] = error;
        for(int i = 0; i < degree_coeff_size; i++) coeff[i] = result[i];
        poly_set_coeffs(poly.poly + j, max_degree, coeff);
        poly_destroy(poly_backup.poly+j);
        poly_copy(poly.poly+j, poly_backup.poly+j);
      }
      fprintf(stderr, "%s: %.4f ", outvname[j], error);
      errorSum += error;
    }

    if(pass2) fprintf(stderr, "\n%d coeffs, fitting error %g\n", sumCoeffs, errorSum);
    else fprintf(stderr, "\ndegree %d has %d/%d coeffs, fitting error %g\n", max_degree, sumCoeffs, maxSumCoeffs, errorSum);
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
    int error = evaluate_aperture(lenses, lenses_cnt, zoom, ray_in, out, aspheric_elements);
    (void)error; // silence non-debug build warning
    assert(error == 0);
    for(int k=0;k<5;k++)
      sample[i+k*sample_cnt] = out[k];
  }
  fprintf(stderr, "[ sensor->aperture ] optimising %d coeffs by %d/%d valid sample points\n", coeff_size, valid, sample_cnt);
  for(int i = 0; i < 5; i++) last_error[i] = FLT_MAX;
  poly_system_copy(&poly_ap, &poly_backup);

  for(max_degree=min_degree_aperture;max_degree <= user_degree; max_degree++)
  {
    int sumCoeffs = 0;
    int maxSumCoeffs = 0;
    float errorSum = 0.0f;
    for(int j = 0; j < 5; j++)
    {
      //const int degree_coeff_size = poly_system_get_coeffs(&poly, max_degree, 0);
      const int degree_coeff_size = poly_get_coeffs(poly_ap.poly + j, max_degree, 0);
      int coeff_cnt = degree_coeff_size;
      const int degree_num_samples = std::min(valid, oversample*degree_coeff_size);

      // optimize taylor polynomial a bit
      Eigen::MatrixXd A(degree_coeff_size, degree_num_samples);
      #pragma omp parallel for
      for(int y = 0; y < degree_num_samples; y++)
      {
        for(int x = 0; x < degree_coeff_size; x++)
        {
          poly_term_t term = poly_ap.poly[j].term[x];
          term.coeff = 1;
          A(x,y) = poly_term_evaluate(&term, sample_in+5*y);
        }
      }
      //cout<<A<<endl;
      Eigen::VectorXd b(degree_num_samples);
      for(int y = 0; y < degree_num_samples; y++)
        b(y) = sample[j*sample_cnt+y];
      //cout<<b<<endl;
      Eigen::VectorXd result = Eigen::ArrayXd::Zero(degree_coeff_size);
      Eigen::VectorXd residual = b;
      Eigen::VectorXd factor(degree_coeff_size);
      Eigen::VectorXd used(degree_coeff_size);
      int permutation[degree_coeff_size];
      for(int i = 0; i < degree_coeff_size; i++) permutation[i] = i;
      for(int i = 0; i < degree_coeff_size; i++) used[i] = 0.0;
      for(int i = 0; i < degree_coeff_size; i++)
        factor(i) = 1 / (A.row(i).norm() * (poly_term_get_degree(poly_ap.poly[j].term+i)+1.0));

      coeff_cnt = 0;
      for(int i = 0; i < degree_coeff_size; i++)
      {
        int maxidx = 0;
        Eigen::VectorXd prod = (Eigen::ArrayXd(A * residual) * Eigen::ArrayXd(factor) * (1-Eigen::ArrayXd(used))).abs();
        prod.maxCoeff(&maxidx);
        permutation[coeff_cnt] = maxidx;
        coeff_cnt++;

        Eigen::MatrixXd tmp2(degree_num_samples, coeff_cnt);
        for(int k = 0; k < coeff_cnt; k++)
          tmp2.col(k) = A.row(permutation[k]).transpose();
        result = (tmp2.transpose()*tmp2).ldlt().solve(tmp2.transpose()*b);
        residual = b-tmp2*result;
        used(maxidx) = 1.0;
        float max_err = residual.squaredNorm();
        if(max_err < precision[j]*degree_num_samples)
          break;
        if(i > max_coeffs) break; // force sparsity
      }

      Eigen::VectorXd coeffs = Eigen::ArrayXd::Zero(degree_coeff_size);
      for(int k = 0; k < coeff_cnt; k++)
        coeffs(permutation[k]) = result(k);
      result = coeffs;
      float error = residual.squaredNorm() / degree_num_samples;

      sumCoeffs += coeff_cnt;
      maxSumCoeffs += degree_coeff_size;
      if(error < last_error[j])
      {
        last_error[j] = error;
        for(int i = 0; i < degree_coeff_size; i++) coeff[i] = result[i];
        poly_set_coeffs(poly_ap.poly + j, max_degree, coeff);
        poly_destroy(poly_backup.poly+j);
        poly_copy(poly_ap.poly+j, poly_backup.poly+j);
      }
      fprintf(stderr, "%s: %.4f ", outvname[j], error);
      errorSum += error;
    }
    if(pass2) fprintf(stderr, "\n%d coeffs, fitting error %g\n", sumCoeffs, errorSum);
    else fprintf(stderr, "\ndegree %d has %d/%d coeffs, fitting error %g\n", max_degree, sumCoeffs, maxSumCoeffs, errorSum);
  }

  poly_system_simplify(&poly_backup);
  // TODO: this totally doesn't throw away useless coeffs, the fitter will make ineffective ones != 0!
  fprintf(stderr, "output aperture poly has %d coeffs.\n", poly_system_get_coeffs(&poly_backup, user_degree, 0));
  poly_system_write(&poly_backup, apfitfile);
  exit(0);
}
