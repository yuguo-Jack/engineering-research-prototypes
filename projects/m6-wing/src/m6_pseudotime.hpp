#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace m6 {
// Consistent physical or nondimensional units. The normal momentum flux is
// the arithmetic face average of rho*velocity, dotted with the area vector.
inline double coupled_time_convective_radius(double momentum_flux,double rho,double reference2,double acoustic2,double area2){
 if(!(rho>0&&reference2>0&&acoustic2>0&&area2>0))throw std::runtime_error("invalid coupled time radius inputs");
 double un=std::abs(momentum_flux/rho),half=.5*(1-reference2/acoustic2)*un;
 return un-half+std::sqrt(reference2*area2+half*half);
}
inline double coupled_time_diffusivity(double rho,double cp,double mu,double mut,double conductivity,double prandtl,double schmidt,double molecular_diffusivity){
 if(!(rho>0&&cp>0&&mu>0&&mut>=0&&conductivity>0&&prandtl>0&&schmidt>0&&molecular_diffusivity>=0))throw std::runtime_error("invalid coupled time diffusivity inputs");
 return std::max({(mu+mut)/rho,(conductivity/cp+mut/prandtl)/rho,molecular_diffusivity+mut/(schmidt*rho)});
}
inline double coupled_time_viscous_radius(double diffusivity,double area,double normal_distance){
 if(!(diffusivity>=0&&area>0&&normal_distance>0))throw std::runtime_error("invalid coupled viscous time radius inputs");
 return 2*diffusivity*area/normal_distance;
}
}
