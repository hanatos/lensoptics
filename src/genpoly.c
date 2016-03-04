#include <math.h>
#include <string.h>
#include "poly.h"

static int degree = 9;

int main(int argc, char *arg[])
{
  poly_system_t s;
  //Generate array {0,0,0,0,0}, {0,0,0,0,1}, {0,0,0,0,2} ... {0,0,0,1,0} 
  // ... {degree, degree, ...}
  //this is essentially the degree-bit representation of a number with
  //poly_num_vars digits (with base degree)
  int numTerms = (int)powf(1.0f+degree, 1.0f*poly_num_vars);
  fprintf(stderr, "polynomial of degree %d has up to %d terms\n", degree, numTerms);
  //int *expMatrix = new int[numTerms*poly_num_vars];
  poly_term_t *expMatrix = malloc((numTerms+1)*poly_num_vars*sizeof(poly_term_t));
  for(int i = 0; i <= numTerms; i++)
  {
    expMatrix[i].coeff = 0.0f;
    int tmp = i;
    for(int j = poly_num_vars-1; j >= 0; j--)
    {
      expMatrix[i].exp[j] = tmp%(degree+1);
      tmp /= (degree+1);
    }
  }
  int validTerms = 0;
  for(int i=0;i<numTerms;i++)
  {
    //remove terms with sum exp > degree
    int *exp = expMatrix[i].exp;
    if(exp[0]+exp[1]+exp[2]+exp[3]+exp[4] <= degree)
      validTerms++;
    //  fprintf(stderr, "Term: %d %d %d %d %d\n", exp[0], exp[1], exp[2], exp[3], exp[4]);
  }
  for(int i=0;i<poly_num_vars;i++)
  {
    s.poly[i].num_terms = validTerms;
    s.poly[i].term = (poly_term_t *)malloc(sizeof(poly_term_t)*validTerms);
    int termsWritten = 0;
    for(int j=0;j<numTerms;j++)
    {
      int *exp = expMatrix[j].exp;
      if(exp[0]+exp[1]+exp[2]+exp[3]+exp[4] <= degree)
        memcpy(s.poly[i].term+(termsWritten++), expMatrix+j, sizeof(poly_term_t));
    }
  }
  fprintf(stderr, "polynomial of degree %d has %d terms\n", degree, validTerms);
  
  free(expMatrix);
  poly_system_sort(&s);
  poly_system_write(&s, "degree9-complete-sorted.poly");
}