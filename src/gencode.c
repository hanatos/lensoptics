#include <math.h>
#include <string.h>
#include "gencode.h"
#include "lenssystem.h"

static lens_element_t lenses[50];
static int lenses_cnt = 0;
static float zoom = 0;

int main(int argc, char **argv)
{
  FILE *f = 0;

  const char *lensfilename = "lenses/ItohZoom.fx";
  if(argc > 1)
    lensfilename = argv[1];

  poly_system_t poly, poly_ap;
  char fitf[1024], afitf[1024];
  snprintf(fitf,  sizeof(fitf),  "%s.fit", lensfilename);
  snprintf(afitf, sizeof(afitf), "%s_ap.fit", lensfilename);
  if(poly_system_read(&poly, fitf) || poly_system_read(&poly_ap, afitf))
  {
    fprintf(stderr, "[gencode] could not read poly fits for `%s'!\n", lensfilename);
    exit(1);
  }

  const char *varnames[poly_num_vars] =
  {
    "x",
    "y",
    "dx",
    "dy",
    "lambda"
  };

  f = fopen("pt_evaluate.h", "wb");
  print_poly_system_code(f, &poly, varnames);
  fclose(f);

  f = fopen("pt_evaluate_jacobian.h", "wb");
  print_jacobian(f, &poly, varnames);
  fclose(f);

  f = fopen("pt_evaluate_aperture.h", "wb");
  print_poly_system_code(f, &poly_ap, varnames);
  fclose(f);

  f = fopen("pt_evaluate_aperture_jacobian.h", "wb");
  print_jacobian(f, &poly_ap, varnames);
  fclose(f);

  f = fopen("pt_sample_ap.h", "wb");
  print_solve_omega(f, &poly_ap, varnames);
  fclose(f);

  f = fopen("lt_sample_ap.h", "wb");
  print_connect(f, &poly, &poly_ap, varnames);
  fclose(f);

  char name[512];
  lens_canonicalize_name(lensfilename, name);
  lenses_cnt = lens_configuration(lenses, lensfilename, sizeof(lenses));
  float lens_length = 0;
  for(int i = 0; i < lenses_cnt; i++) lens_length += lens_get_thickness(lenses + i, zoom);
  const float aperture_housing_radius = lens_get_aperture_radius(lenses, lenses_cnt);
  const float aperture_pos = lens_get_aperture_pos(lenses, lenses_cnt, zoom);
  const float bfl = lens_get_thickness(lenses + lenses_cnt-1, zoom);

  f = fopen("init.h", "wb");
  fprintf(f, "static const char *lens_name = \"%s\"; // descriptive name of the lens\n", name);
  fprintf(f, "static const float lens_outer_pupil_radius = %f; // scene facing radius in mm\n", lenses[0].housing_radius);
  fprintf(f, "static const float lens_inner_pupil_radius = %f; // sensor facing radius in mm\n", lenses[lenses_cnt-1].housing_radius);
  fprintf(f, "static const float lens_length = %f; // overall lens length in mm\n", lens_length);
  fprintf(f, "static const float lens_focal_length = %f; // approximate lens focal length in mm (BFL)\n", bfl);
  fprintf(f, "static const float lens_aperture_pos = %f; // distance aperture -> outer pupil in mm\n", aperture_pos);
  fprintf(f, "static const float lens_aperture_housing_radius = %f; // lens housing radius at the aperture\n", aperture_housing_radius);
  fprintf(f, "static const float lens_outer_pupil_curvature_radius = %f; // radius of curvature of the outer pupil\n", lenses[0].lens_radius);
  fclose(f);

  exit(0);
}
