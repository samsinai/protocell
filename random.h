#ifndef	RANDOM
#define	RANDOM

#include <math.h>

static const int IM1 = 2147483563;
static const int IM2 = 2147483399;
static const int IA1 = 40014;
static const int IA2 = 40692;
static const int IQ1 = 53668;
static const int IQ2 = 52774;
static const int IR1 = 12211;
static const int IR2 = 3791;
static const int NTAB = 32;
static const int IMM1 = IM1-1;
static const int NDIV = (1+IMM1/NTAB);
static const double REPS = 3.0e-16;
static const double RNMX = (1.0-REPS);
static const double AM = (1.0/double(IM1));

// for Gaussian
static const double SQRT_TWO_PI = sqrt(2*M_PI);
static const double INV_SQRT_TWO_PI = 1/sqrt(2*M_PI);

class Random
{
 public:
  /*
   * double between zero and one from uniform distribution
   */
  static double Uniform(int& dummy);
  /* 
   *double between min and max from uniform distribution
   */
  static double Uniform(int& dummy, int min, int max);
  /*
   * int between min and max from uniform distribution
   */
  // static int UniformInt(int& dummy, int min, int max);

  /*
   * double from Gaussian distribution with mean zero and
   * scale one.
   */
  static double Gaussian(int&);
  /*
   * double from Gaussian distribution with mean
   *as supplied and scale one.
   */
  static double Gaussian(int&, double&);
  /*
   * double from Gaussian distribution with mean and
   * scale as supplied.
   */
  static double Gaussian(int&, double&, double&);

/*   template <class T> static void RandPerm(int& dummy, T* array, int len) */
/*     { */
/*       T temp; */
/*       int r; */
/*       for (int i = 0; i < len; i ++) */
/* 	{ */
/* 	  temp = array[i]; */
/* 	  r = UniformInt(dummy, 0, len-1); */
	  
/* 	  // swap */
/* 	  array[i] = array[r]; */
/* 	  array[r] = temp;	   */
/* 	} */
/*     } */

}; // RANDOM_H
  
#endif
