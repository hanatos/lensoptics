#include <string.h>
#include <math.h>
#include "poly.h"

int main(int argc, char *arg[])
{
  if(argc < 2)
  {
    fprintf(stderr, "usage: %s input.poly\n", arg[0]);
    exit(-1);
  }
  
  poly_system_t s;
  poly_system_read(&s, arg[1]);
  
  for(int i = 0; i < poly_num_vars; i++)
  {
    static const char *vn[poly_num_vars] = {"x0", "x1", "x2", "x3", "x4"};
    fprintf(stdout, "%s = ", vn[i]);
    for(int t=0;t<p->num_terms;t++)
    {
      fprintf(f, " + %g ", p->term[t].coeff);
      for(int k=0;k<poly_num_vars;k++)
        if(p->term[t].exp[k] == 1)
          fprintf(f, "*%s", vname ? vname[k] : vn[k]);
        else if(p->term[t].exp[k] > 0)
          fprintf(f, "*lens_ipow(%s, %d)", vname ? vname[k] : vn[k], p->term[t].exp[k]);
    }
    fprintf(stdout, "\n\n");
  }
}
