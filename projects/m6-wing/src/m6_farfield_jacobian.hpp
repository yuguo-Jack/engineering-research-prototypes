#pragma once
#include "m6_farfield.hpp"

namespace m6 {
// Approximate implicit characteristic boundary linearization, separately
// verified on actual normal assemblies for branches 1, 2, 5 and 6. The flux
// itself still comes from characteristic_farfield. Units are ND, per area.
inline B5 approximate_farfield_primitive_jacobian(const Gas&gas,const W7&cell,const W7&face,
 const W7&infinity,const D3&normal,double cell_rho,double cell_a2,double face_rho,
 double cell_cp,double cell_rg,double cell_total_enthalpy,const FarfieldEvaluation&evaluated){
 const int branch=evaluated.branch;
 if(branch!=1&&branch!=2&&branch!=5&&branch!=6)throw std::runtime_error("unobserved approximate farfield branch");
 const auto&b=evaluated.boundary;const double mass=evaluated.flux[0];
 const double cp_inf=gas.cp(infinity[4]),gamma=cp_inf/(cp_inf-gas.rg),gm1=gamma-1;
 const double cc=std::sqrt(cell_a2),a_inf=gas.sound(infinity);
 D3 velocity={b[1],b[2],b[3]},inf_velocity={infinity[1],infinity[2],infinity[3]},face_velocity={face[1],face[2],face[3]};
 const double H=branch==5?gas.enthalpy(infinity):branch==6?cell_total_enthalpy:gas.enthalpy(b);
 Q5 dm{},dp{},dh{};std::array<Q5,3>du{};
 if(branch==1||branch==2){
  const double un=dot3(velocity,normal),uf=dot3(face_velocity,normal),ui=dot3(inf_velocity,normal);
  const double bc=.25*gm1*(uf+2*cc/gm1-ui+2*a_inf/gm1);
  const double rho=std::abs(un)>1e-30?mass/un:b[0]/(gas.rg*b[4]);
  for(int j=0;j<3;++j){dm[j+1]=.5*rho*(1+un/bc)*normal[j];dp[j+1]=.5*gamma*b[0]/bc*normal[j];du[j][4]=normal[j]*cc/(2*gm1*cell[4]);}
  if(branch==1){
   dm[4]=rho*cc/(2*gm1*cell[4])*(1+un/bc);
   dp[4]=gamma*b[0]*cc/(2*gm1*bc*cell[4]);
   for(int i=0;i<3;++i)for(int j=0;j<3;++j)du[i][j+1]=.5*normal[i]*normal[j];
  }else{
   dm[0]=mass/cell[0];dm[4]=rho/(2*gm1*cell[4])*(cc+un*(cc/bc-2*gamma));
   dp[0]=b[0]/cell[0];dp[4]=gamma*b[0]/(gm1*cell[4])*(cc/(2*bc)-1);
   for(int i=0;i<3;++i)for(int j=0;j<3;++j)du[i][j+1]=(i==j?1.:0.)-.5*normal[i]*normal[j];
  }
  for(int j=0;j<5;++j)for(int i=0;i<3;++i)dh[j]+=velocity[i]*du[i][j];
  // The observed outlet implicit operator uses this characteristic thermal
  // sensitivity although the residual's outlet temperature equals cell T.
  for(int j=0;j<3;++j)dh[j+1]+=cp_inf*gm1*b[4]/(2*bc)*normal[j];
  dh[4]+=cp_inf*b[4]*cc/(2*bc*cell[4]);
 }else{
  const double a=.5*(cc+a_inf),uf=dot3(face_velocity,normal),mach=uf/a;
  double p_plus=mach>0?1.:0.,p_prime=0,m_prime=mach>0?1.:0.;
  if(std::abs(mach)<1){p_plus=.25*(mach+1)*(mach+1)*(2-mach)+.1875*mach*(mach*mach-1)*(mach*mach-1);p_prime=.9375*(mach*mach-1)*(mach*mach-1);m_prime=.5*(1+mach*mach*mach);}
  dp[0]=p_plus;for(int j=0;j<3;++j)dp[j+1]=face[0]*p_prime/a*normal[j];
  if(branch==5){for(int j=0;j<3;++j)dm[j+1]=gas.rho(infinity)*m_prime*normal[j];}
  else{
   dm[0]=uf/(cell_rg*cell[4]);dm[4]=-uf*cell_rho/cell[4];
   for(int j=0;j<3;++j){dm[j+1]=face_rho*normal[j];du[j][j+1]=1;dh[j+1]=face[j+1];}dh[4]=cell_cp;
  }
 }
 B5 result{};for(int col=0;col<5;++col){
  result[col]=dm[col];for(int i=0;i<3;++i)result[5*(i+1)+col]=velocity[i]*dm[col]+mass*du[i][col]+normal[i]*dp[col];
  result[20+col]=H*dm[col]+mass*dh[col];
 }
 // Retain the observed coordinate coupling, independently established by
 // before/after actual blocks and an explicit disabled-term control.
 if(branch==2)result[15]=velocity[0]*dm[0]+normal[2]*dp[0];
 return result;
}

inline B5 normal_farfield_jacobian(const Gas&gas,const W7&cell,const W7&face,
 const W7&infinity,const D3&normal,double face_rho){
 const double rho=gas.rho(cell),a2=std::pow(gas.sound(cell),2);
 auto evaluated=characteristic_farfield(gas,cell,face,infinity,normal,rho,a2,face_rho,gas.enthalpy(cell));
 // Supersonic boundary branches were not exercised by the retained normal
 // observation. Preserve the existing finite-difference path for them.
 if(evaluated.branch!=1&&evaluated.branch!=2&&evaluated.branch!=5&&evaluated.branch!=6)return characteristic_farfield_jacobian(gas,cell,face,infinity,normal,face_rho,true);
 auto primitive=approximate_farfield_primitive_jacobian(gas,cell,face,infinity,normal,rho,a2,face_rho,gas.cp(cell[4]),gas.rg,gas.enthalpy(cell),evaluated);
 auto dp=gas.pressure_jacobian(cell);auto dw=gas.primitive_jacobian(cell);B5 result{};
 for(int row=0;row<5;++row)for(int col=0;col<5;++col)for(int k=0;k<5;++k)result[5*row+col]+=primitive[5*row+k]*(k==0?dp[col]:dw[5*k+col]);
 return result;
}
}
