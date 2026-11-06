#pragma once
#include <cmath>
#include <stdexcept>

namespace m6 {
// The steady density-based outer iteration still supplies a local pseudo-time
// to the scalar absolute system. Density belongs to the post-flow scalar stage;
// dt is retained from the pre-flow coupled assembly.
inline double scalar_time_coefficient(double rho,double volume,double dt){
 if(!(rho>0&&volume>0&&dt>0&&std::isfinite(rho)&&std::isfinite(volume)&&std::isfinite(dt)))throw std::runtime_error("invalid scalar pseudo-time inputs");
 return (rho*volume)/dt;
}
}
