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

  FILE *f = stdout;

  for(int i = 0; i < poly_num_vars; i++)
  {
    poly_t *p = s.poly+i;
    poly_term_t *term = p->term;
    for(int t=0;t<p->num_terms;t++,term++)
    {
      fprintf(f, "%g,%d,%d,%d,%d,%d\n", term->coeff, term->exp[0], term->exp[1],
        term->exp[2], term->exp[3], term->exp[4]);
    }
    fprintf(f, "\n");
  }
}
