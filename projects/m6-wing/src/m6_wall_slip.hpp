#pragma once
#include "m6_physics.hpp"
namespace m6 {
struct WallSlipFactor {double value=0,derivative=0;};
inline double wall_friction_velocity_from_k(double rho,double mu,double k,double tangent_speed,double distance){
 if(!(rho>0&&mu>0&&k>0&&tangent_speed>=0&&distance>0))throw std::runtime_error("invalid wall friction-velocity input");
 return std::sqrt((mu/rho)*tangent_speed/distance+.3*k);
}
// Smooth stationary wall, zero pressure-gradient/blowing/thermal parameters.
// The velocity law AND its slope are independently blended; the latter is
// not the derivative of the blended velocity law. eps is part of each blend.
inline WallSlipFactor wall_slip_factor(double y){
 constexpr double eps=1e-10,kappa=.4187,E=9.793;
 double yl=std::max(y,.001),yt=std::max(y,.2),a=yl+eps,b=std::log(E*yt)/kappa+eps;
 double p=std::pow(a,-4)+std::pow(b,-4),uplus=std::pow(p,-.25);
 double da=y>.001?1.:0.,db=y>.2?1/(kappa*y):0.;
 double uplus_prime=std::pow(p,-1.25)*(std::pow(a,-5)*da+std::pow(b,-5)*db);
 double c=1+eps,d=1/(kappa*yt)+eps,pd=std::pow(c,-4)+std::pow(d,-4),slope=std::pow(pd,-.25);
 double dd=y>.2?-1/(kappa*y*y):0.,slope_prime=std::pow(pd,-1.25)*std::pow(d,-5)*dd;
 double den=uplus+eps;return {1-y*slope/den,-(slope+y*slope_prime)/den+y*slope*uplus_prime/(den*den)};
}
inline D3 wall_virtual_velocity_from_star(double ustar,double distance,double rho,double mu,const D3&velocity,const D3&wall,const D3&tangent){
 if(!(ustar>=0&&distance>0&&rho>0&&mu>0))throw std::runtime_error("invalid wall slip input");
 double y=ustar*distance*rho/mu;auto factor=wall_slip_factor(y);D3 result{};
 for(int j=0;j<3;++j){result[j]=wall[j]+factor.value*tangent[j];if(y>11.224708)result[j]+=velocity[j]-wall[j]-tangent[j];}
 return result;
}
struct WallVirtualState {D3 velocity{};std::array<double,12>jacobian{};double yplus=0;};
// Primitive derivatives with respect to u,v,w,rho; k and molecular mu frozen.
inline WallVirtualState stationary_wall_virtual_state(double rho,double mu,double k,const D3&velocity,const D3&normal,double distance){
 if(!(rho>0&&mu>0&&k>0&&distance>0))throw std::runtime_error("invalid stationary wall gradient input");
 double un=dot3(velocity,normal);D3 tangent{};for(int j=0;j<3;++j)tangent[j]=velocity[j]-un*normal[j];
 double speed=length3(tangent),nu=mu/rho,ustar=wall_friction_velocity_from_k(rho,mu,k,speed,distance),y=ustar*distance/nu;
 auto factor=wall_slip_factor(y);WallVirtualState result;result.yplus=y;
 double dy_ds=1/(2*ustar),dy_drho=distance*ustar/mu-speed/(2*rho*ustar);bool restore_normal=y>11.224708;
 for(int i=0;i<3;++i){result.velocity[i]=factor.value*tangent[i]+(restore_normal?un*normal[i]:0.);
  for(int j=0;j<3;++j){double projection=(i==j?1.:0.)-normal[i]*normal[j];result.jacobian[4*i+j]=factor.value*projection+(speed>1e-30?factor.derivative*dy_ds*tangent[i]*tangent[j]/speed:0.)+(restore_normal?normal[i]*normal[j]:0.);}
  result.jacobian[4*i+3]=factor.derivative*dy_drho*tangent[i];
 }return result;
}
}
