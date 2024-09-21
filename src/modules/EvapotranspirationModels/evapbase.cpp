#include "EvapotranspirationModels/evapbase.hpp"

double evapT_base::delta(double& t) // Slope of sat vap p vs t, kPa/DEGREE_CELSIUS
{
  if (t > 0.0)
    return(2504.0*exp(17.27 * t/(t+237.3)) / sqr(t+237.3));
  else
    return(3549.0*exp( 21.88 * t/(t+265.5)) / sqr(t+265.5));
}

double evapT_base::lambda(double& t) // Latent heat of vaporization (mJ/(kg DEGREE_CELSIUS))
{
   return( 2.501 - 0.002361 * t );
}

double evapT_base::gamma(double& Pa, double& t) // Psychrometric constant (kPa/DEGREE_CELSIUS)
{
   return( 0.00163 * Pa / lambda(t)); // lambda (mJ/(kg DEGREE_CELSIUS))
}
