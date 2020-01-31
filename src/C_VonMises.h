#include <math.h>
#include <stdlib.h>
static inline float parcels_vonmisesvariate(float mu, float kappa)
/* Circular data distribution. */
/* Returns a float between 0 and 2*pi */
/* mu is the mean angle, expressed in radians between 0 and 2*pi, and  */
/* kappa is the concentration parameter, which must be greater than or */
/* equal to zero.  If kappa is equal to zero, this distribution reduces */
/* to a uniform random angle over the range 0 to 2*pi. */
/* Based upon an algorithm published in: Fisher, N.I.,*/
/* Statistical Analysis of Circular Data", Cambridge University Press, 1993.*/
{
  float u1, u2, u3, r, s, z, d, f, q, theta;
  
  if (kappa <= 1e-6){
    return (2.0 * M_PI * (float)rand()/(float)(RAND_MAX));
  }
    
  s = 0.5 / kappa;
  r = s + sqrt(1.0 + s * s);
  
  do {
    u1 = (float)rand()/(float)(RAND_MAX);
    z = cos(M_PI * u1);
    
    d = z / (r + z);
    u2 = (float)rand()/(float)(RAND_MAX);
  }  while ( ( u2 >= (1.0 - d * d) ) && ( u2 > (1.0 - d) * exp(d) ) );
        
  q = 1.0 / r;
  f = (q + z) / (1.0 + q * z);
  u3 = (float)rand()/(float)(RAND_MAX);
  
  if (u3 > 0.5){
    theta = mu + acos(f), 2*M_PI ;
  }
  else {
    theta = mu - acos(f), 2*M_PI;
  }
  if (theta < 0){
    theta = 2*M_PI+theta;
  }
  
  return theta;
}


#Cette fonction est tirée directement du livre. Elle marche mais elle est beaucoup trop lente.
static inline float parcels_vonmisesvariate(float mu, float kappa)
/* Circular data distribution. */
/* mu is the mean angle, expressed in radians between 0 and 2*pi, and  */
/* kappa is the concentration parameter, which must be greater than or */
/* equal to zero.  If kappa is equal to zero, this distribution reduces */
/* to a uniform random angle over the range 0 to 2*pi. */
/* Based upon an algorithm published in: Fisher, N.I.,*/
/* Statistical Analysis of Circular Data", Cambridge University Press, 1993.*/
{
  float u1, u2, u3, a, b, c, r, z, f, theta;
  
  if (kappa <= 1e-6){
    return (2.0 * M_PI * (float)rand()/(float)(RAND_MAX));
  }
  
  a = 1 + sqrt(1 + 4 * kappa * kappa);
  b = (a - sqrt(2 * a)) / 2 * kappa;
  r = (1 + b * b) / (2 * b);
  
  
  do {
    u1 = (float)rand()/(float)(RAND_MAX);
    z = cos(M_PI * u1);
    f = (1 + r * z) / (r + z);
    c = kappa * (r - f);
    u2 = (float)rand()/(float)(RAND_MAX);
  }  while ( ( (c * (2.0 - c) - u2) <= 0.0 ) && ( (log(c / u2) + 1 - c) < 0.0 ) );
  
  
    
  u3 = (float)rand()/(float)(RAND_MAX);
  if (u3 > 0.5){
    theta = mu + acos(f);
  }
  else {
    theta = mu - acos(f);
  }
  if (theta < 0){
    theta = 2*M_PI+theta;
  }
  
  return theta;
}


int main()
{   float a; 
    int i = 1;
    
    while (i <= 10000)
    {   a=parcels_vonmisesvariate(M_PI, 4);
        printf("%f", a);
        printf(", ");
        ++i;
    }
}


#ifdef __cplusplus
}
#endif
#endif