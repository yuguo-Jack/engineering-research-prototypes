#pragma once
#include "m6_physics.hpp"
namespace m6 {
struct ViscousFlux {D3 traction{};double energy=0;};
// Inputs are the already corrected physical gradients and effective transport.
// Works in either consistently dimensional or nondimensional units.
inline ViscousFlux viscous_flux(const D3&normal,const std::array<D3,4>&g,double mu,double conductivity,const D3&velocity){
 ViscousFlux result;double div=g[0][0]+g[1][1]+g[2][2];
 for(int j=0;j<3;++j)for(int k=0;k<3;++k)result.traction[j]+=mu*(g[j][k]+g[k][j]-(j==k?2*div/3:0))*normal[k];
 result.energy=conductivity*dot3(g[3],normal);for(int j=0;j<3;++j)result.energy+=result.traction[j]*velocity[j];
 return result;
}
// Normal symmetry callback: zero tangential normal derivatives, mirror the
// cell's normal velocity, and evaluate viscous work at projected CELL velocity.
// The normal thermal/species flux is zero. Inputs use either consistent units.
inline ViscousFlux symmetry_viscous_flux(const D3&normal,double distance,const std::array<D3,3>&physical,double mu,const D3&velocity){
 if(!(distance>0))throw std::runtime_error("nonpositive symmetry normal distance");
 std::array<D3,4>g{};D3 tangent=velocity;double un=dot3(velocity,normal);
 for(int k=0;k<3;++k){tangent[k]-=un*normal[k];double correction=-un*normal[k]/distance-dot3(physical[k],normal);for(int j=0;j<3;++j)g[k][j]=physical[k][j]+correction*normal[j];}
 return viscous_flux(normal,g,mu,0.,tangent);
}
}
