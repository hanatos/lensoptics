#include "raytrace.h"

static lens_element_t lenses[50];
static int lenses_cnt = 0;
static const float zoom = 0.0f;

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
  const int sample_cnt = 500000;
  FILE *f = fopen("fresnel.dat", "wb");
  //XXX better find a way to set color-range in gnuplot than to add two points
  fprintf(f, "%g, %g, %g\n", 0.f, 0.f, 1.f);
  fprintf(f, "%g, %g, %g\n", 0.f, 0.f, 0.f);
  for(int i=0;i<sample_cnt;i++)
  {
    const float u = drand48(), v = drand48(), w = drand48(), x = drand48(), y = drand48();
    float ray_in[] = {
      (x-0.5)*36.0f,
      (y-0.5)*24.0f,
      p_rad/p_dist * cosf(2.0f*M_PI*u)*sqrtf(v),
      p_rad/p_dist * sinf(2.0f*M_PI*u)*sqrtf(v),
      0.4 + 0.3*w};
    ray_in[2] -= ray_in[0] / p_dist;
    ray_in[3] -= ray_in[1] / p_dist;
    float out[5] = {0.0f, 0.0f, 0.0f, 0.0f, ray_in[4]};
    int error = evaluate(lenses, lenses_cnt, zoom, ray_in, out);
    if(!error)
    {
      float intensity = out[4];
      fprintf(f, "%g, %g, %g\n", out[0], out[1], intensity);
    }
    else
      --i;
  }
  fclose(f);
  exit(0);
}
