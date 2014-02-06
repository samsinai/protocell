#include <math.h>
#include "random.h"

using namespace std;

/*
 * This is ran2 from Numerical Recipes in C++
 * Press et al. 2002
 */
double Random::Uniform(int& dummy)
{
  int j,k;
  static int dummy2 = 123456789;
  static int iy = 0;

  static int iv[NTAB];

  double temp;

  if (dummy <= 0)
    {
      dummy = (dummy == 0) ? 1 : -dummy;
      dummy2 = dummy;
      
      for (j = NTAB + 7; j >= 0; j --)
	{
	  k = (dummy)/IQ1;
	  dummy = IA1*(dummy - k*IQ1) - k*IR1;
	  if (dummy < 0)
	    {
	      dummy += IM1;
	    }
	  if (j < NTAB)
	    {
	      *(iv + j) = dummy;
	    }
	}
      iy = *iv;
    }

  k = (dummy)/IQ1;
  dummy = IA1*(dummy - k*IQ1) - k*IR1;

  if (dummy < 0)
    {
      dummy += IM1;
    }

  k = dummy/IQ2;
  dummy2 = IA2*(dummy2 - k*IQ2) - k*IR2;

  if (dummy2 < 0)
    {
      dummy2 += IM2;
    }

  j = iy/NDIV;
  iy = *(iv + j) - dummy2;
  *(iv + j) = dummy;

  if (iy < 1)
    {
      iy += IMM1;
    }
  
  if ((temp = AM*iy) > RNMX)
    {
      return RNMX;
    }
  else
    {
      return temp;
    }
}

// int Random::UniformInt(int& dummy, int min, int max)
// {
//   return round(min + max*Uniform(dummy));
// }


double Random::Gaussian(int& dummy)
{
  static int iset = 0;
  static double gset;
  double r1, r2, rsq, fac;
  
  if (dummy < 0) iset = 0;
  if (iset == 0)
    {
      do
	{
	  r1 = 2.0*Uniform(dummy) - 1;
	  r2 = 2.0*Uniform(dummy) - 1;	  
	  rsq = r1*r1 + r2*r2;
	}
      while(rsq >= 1.0 || rsq == 0.0);

      fac = sqrt(-2.0*log(rsq)/rsq);
      gset = r1*fac;
      iset = 1;
      return r2*fac;
    }
  else
    {
      iset = 0;
      return gset;
    }    
}
