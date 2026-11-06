#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace m6 {
struct SmoothWallThermal {
 double tplus=0,slope=0,conductance=0,viscous_source=0,adiabatic_jump=0,virtual_fraction=0,fraction_y_derivative=0,fraction_pr_derivative=0;
};
// Smooth stationary enhanced wall; zero blowing and pressure-gradient terms,
// Jaya thermal correction enabled, non-SA model. Inputs use consistent units.
// The thermal slope is a separate blend, not d(Tplus)/d(yplus).
inline SmoothWallThermal smooth_wall_thermal(double yplus,double ustar,double rho,double cp,double pr,double prt,double tangent_speed,bool rough_wall_enabled=true){
 for(double x:{yplus,ustar,rho,cp,pr,prt,tangent_speed})if(!std::isfinite(x))throw std::runtime_error("nonfinite thermal wall input");
 if(!(yplus>0&&ustar>0&&rho>0&&cp>0&&pr>0&&prt>0&&tangent_speed>=0))throw std::runtime_error("invalid thermal wall input");
 constexpr double kappa=.4187,E=9.793;
 double y=std::max(yplus,.001),chi=rough_wall_enabled?std::min(1.,9.783206999999999*yplus):1.;
 // The rough-wall switch also caps the smooth-wall low-y logarithm/blend.
 double blend=std::exp(chi-1),py=pr*y;
 double gamma=-.01*blend*py*py*py*py/(1+(5/blend)*pr*pr*pr*y);
 double laminar=std::exp(gamma),turbulent=gamma==0?0:std::exp(1/gamma),ratio=pr/prt;
 double jaya_factor=1+.28*std::exp(-.007*ratio),p=9.24*(std::pow(ratio,.75)-1)*jaya_factor;
 double recovery=std::abs(ratio-1)<1e-10?6.93*jaya_factor:p/(ratio-1);
 double log_velocity=std::log((E/chi)*y)/kappa;
 SmoothWallThermal result;
 result.tplus=pr*y*laminar+prt*turbulent*(log_velocity+p);
 if(!(result.tplus>0&&std::isfinite(result.tplus)))throw std::runtime_error("invalid thermal wall Tplus");
 result.slope=pr*laminar+prt*turbulent/(kappa*y);
 result.conductance=rho*cp*ustar/result.tplus;
 double speed2=tangent_speed*tangent_speed;
 result.viscous_source=(-.5*rho*ustar/result.tplus)*(pr*laminar*speed2+prt*turbulent*(speed2+(ratio-1)*recovery*ustar*recovery*ustar));
 result.adiabatic_jump=-result.viscous_source/result.conductance;
 result.virtual_fraction=std::min(1.,yplus*result.slope/result.tplus);
 if(result.virtual_fraction<1){
  double yp=yplus>.001?1.:0.,chip=rough_wall_enabled&&9.783206999999999*yplus<1?9.783206999999999:0.;
  double denominator_term=(5/blend)*pr*pr*pr*y,denominator_fraction=denominator_term/(1+denominator_term);
  double gy=gamma*(chip+4*yp/y-denominator_fraction*(yp/y-chip)),gp=gamma*(4-3*denominator_fraction)/pr;
  double ly=laminar*gy,lp=laminar*gp,ty=turbulent==0?0:-turbulent*gy/(gamma*gamma),tp=turbulent==0?0:-turbulent*gp/(gamma*gamma);
  double log_y=(yp/y-chip/chi)/kappa;
  double p_pr=9.24/prt*(.75*std::pow(ratio,-.25)*jaya_factor+(std::pow(ratio,.75)-1)*.28*std::exp(-.007*ratio)*(-.007));
  double t_y=pr*yp*laminar+pr*y*ly+prt*ty*(log_velocity+p)+prt*turbulent*log_y;
  double t_pr=y*laminar+pr*y*lp+prt*tp*(log_velocity+p)+prt*turbulent*p_pr;
  double s_y=pr*ly+prt/(kappa*y)*ty-prt*turbulent/(kappa*y*y)*yp,s_pr=laminar+pr*lp+prt/(kappa*y)*tp;
  result.fraction_y_derivative=(result.slope+yplus*s_y-result.virtual_fraction*t_y)/result.tplus;
  result.fraction_pr_derivative=(yplus*s_pr-result.virtual_fraction*t_pr)/result.tplus;
 }
 return result;
}
// Adiabatic heat balance in the observed turbulent-wall branch. Conductivity
// epsilon has conductivity units: 1e-8 W/(m K) in the captured physical calls;
// a nondimensional caller must convert it with the conductivity scale.
inline double smooth_adiabatic_wall_temperature(const SmoothWallThermal& thermal,double cell_temperature,double distance,double molecular_conductivity,double conductivity_epsilon){
 if(!(cell_temperature>0&&distance>0&&molecular_conductivity>0&&conductivity_epsilon>=0))throw std::runtime_error("invalid adiabatic thermal wall input");
 double effective_conductivity=std::max(molecular_conductivity,distance*thermal.conductance);
 return cell_temperature-thermal.viscous_source*distance/(effective_conductivity+conductivity_epsilon);
}
}
