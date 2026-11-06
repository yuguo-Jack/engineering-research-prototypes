#pragma once
#include "m6_physics.hpp"

namespace m6 {
struct FarfieldEvaluation{Q5 flux{};W7 boundary{};int branch=0;};
// Independent implementation for the supplied Case controls: characteristic
// subsonic far field, M4/P5 near tangency, positive tangent flux as an outlet.
// Acoustics, density and cached cell total enthalpy are explicit inputs so that
// the same compiled kernel can be checked against live native observations.
inline FarfieldEvaluation characteristic_farfield(const Gas&gas,const W7&cell,const W7&face,const W7&infinity,const D3&normal,double cell_rho,double cell_a2,double face_rho,double cell_total_enthalpy){
 FarfieldEvaluation result;auto&b=result.boundary;double rho=0,volume=0,H=0;
 D3 inf_velocity={infinity[1],infinity[2],infinity[3]};double un_inf=dot3(inf_velocity,normal),un_face=0;
 for(int j=0;j<3;++j)un_face+=face[j+1]*normal[j];double ai=std::sqrt(cell_a2),ainf=gas.sound(infinity),rho_inf=gas.rho(infinity);
 double direction_cosine=un_inf/length3(inf_velocity);bool tangent=std::abs(direction_cosine)<.05;
 if(tangent){
  double a=.5*(ai+ainf),ml=un_face/a,mr=un_inf/a;
  auto mp=[](double m){return std::abs(m)>=1?.5*(m+std::abs(m)):.25*(m+1)*(m+1)+.125*(m*m-1)*(m*m-1);};
  auto mm=[](double m){return std::abs(m)>=1?.5*(m-std::abs(m)):-.25*(m-1)*(m-1)-.125*(m*m-1)*(m*m-1);};
  auto pp=[](double m){return std::abs(m)>=1?(m>=0?1.:0.):.25*(m+1)*(m+1)*(2-m)+.1875*m*(m*m-1)*(m*m-1);};
  auto pm=[](double m){return std::abs(m)>=1?(m<=0?1.:0.):.25*(m-1)*(m-1)*(2+m)-.1875*m*(m*m-1)*(m*m-1);};
  volume=a*(mp(ml)+mm(mr));bool inflow=volume<=0;b=inflow?infinity:face;
  rho=inflow?rho_inf:face_rho;H=inflow?gas.enthalpy(infinity):cell_total_enthalpy;
  if(!inflow)volume=un_face;b[0]=pp(ml)*face[0]+pm(mr)*infinity[0];result.branch=inflow?5:6;
 }else if(std::abs(un_face)>=ai){
  bool inflow=un_face<=0;b=inflow?infinity:face;rho=inflow?rho_inf:face_rho;volume=inflow?un_inf:un_face;H=gas.enthalpy(b);result.branch=inflow?3:4;
 }else{
  double cp=gas.cp(infinity[4]),gamma=cp/(cp-gas.rg),gm1=gamma-1;
  double minus=un_inf-2*ainf/gm1,plus=un_face+2*ai/gm1;
  volume=.5*(plus+minus);double sound=.25*gm1*(plus-minus);bool inflow=volume<=0;
  double entropy=inflow?infinity[0]/std::pow(rho_inf,gamma):cell[0]/std::pow(cell_rho,gamma);
  rho=std::pow(sound*sound/(gamma*entropy),1/gm1);b=inflow?infinity:face;b[0]=entropy*std::pow(rho,gamma);
  double change=volume-(inflow?un_inf:un_face);for(int j=0;j<3;++j)b[j+1]+=change*normal[j];
  b[4]=inflow?b[0]/(gas.rg*rho):cell[4];H=gas.enthalpy(b);result.branch=inflow?1:2;
 }
 double mass=rho*volume;result.flux={mass,mass*b[1]+b[0]*normal[0],mass*b[2]+b[0]*normal[1],mass*b[3]+b[0]*normal[2],mass*H};
 if(!(rho>0&&b[0]>0&&b[4]>0))throw std::runtime_error("nonpositive characteristic far-field state");
 for(double x:result.flux)if(!std::isfinite(x))throw std::runtime_error("nonfinite far-field flux");return result;
}
// Boundary-only centered numerical linearization. Face reconstruction increments
// remain frozen; cell acoustics and cached cell H follow the perturbed cell.
// Only 1632 original far-field faces use this path, not the interior face loop.
inline B5 characteristic_farfield_jacobian(const Gas&gas,const W7&cell,const W7&face,const W7&infinity,const D3&normal,double face_rho,bool independent_density){
 auto original=gas.conserved(cell);B5 matrix{};double original_rho=gas.rho(cell);
 for(int column=0;column<5;++column){double eps=1e-6*std::max(std::abs(original[column]),.01);std::array<Q5,2>values{};
  for(int side=0;side<2;++side){auto q=original;q[column]+=(side==0?eps:-eps);auto changed=gas.primitive(q,cell[4],cell[5],cell[6]);auto reconstructed=face;
   for(int j=0;j<5;++j)reconstructed[j]+=changed[j]-cell[j];double rho=gas.rho(changed),rf=independent_density?face_rho+rho-original_rho:gas.rho(reconstructed);
   values[side]=characteristic_farfield(gas,changed,reconstructed,infinity,normal,rho,std::pow(gas.sound(changed),2),rf,gas.enthalpy(changed)).flux;
  }
  for(int row=0;row<5;++row)matrix[5*row+column]=(values[0][row]-values[1][row])/(2*eps);
 }return matrix;
}
}
