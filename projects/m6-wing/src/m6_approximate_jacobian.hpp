#pragma once
#include "m6_flux_jacobian.hpp"

namespace m6 {
// The normal AUSM residual and its approximate implicit operator are different
// formulas. This stationary, ideal-mixture five-equation approximation was
// checked against actual assembled face increments, not callback replay.
// All inputs and output blocks are nondimensional, per unit face area.
struct FaceCellThermodynamics { double rho,cp,rg; };

inline FluxJacobians approximate_coupled_primitive_jacobian(
 const W7&cell_l,const W7&cell_r,const W7&face_l,const W7&face_r,
 const D3&normal,double mass,const std::array<double,4>&face_thermo,
 double a2_l,double a2_r,double ref_l,double ref_r,
 const std::array<FaceCellThermodynamics,2>&cell_thermo,bool stage_p=true){
 const bool left_upwind=mass>0;
 const double rho=face_thermo[left_upwind?0:1],volume_flux=mass/rho;
 const double reference2=.5*(ref_l+ref_r);
 const double alpha=stage_p?.5*(1-.5*(ref_l/a2_l+ref_r/a2_r)):.5;
 if(!(rho>0&&reference2>0))throw std::runtime_error("nonpositive approximate face linearization inputs");
 const double wave=std::sqrt(alpha*alpha*volume_flux*volume_flux+reference2);
 const double plus=std::abs((1-alpha)*volume_flux+wave);
 const double minus=std::abs((1-alpha)*volume_flux-wave);
 const double difference=(plus-minus)/(2*wave),average=.5*(plus+minus),absolute=std::abs(volume_flux);
 const double pressure_coefficient=(alpha*(2*absolute-volume_flux*difference)+average-absolute)/reference2;
 const double velocity_jump=(face_l[1]-face_r[1])*normal[0]+(face_l[2]-face_r[2])*normal[1]+(face_l[3]-face_r[3])*normal[2];
 const double dissipation_mass=.5*((face_l[0]-face_r[0])*pressure_coefficient+velocity_jump*difference*rho);
 const auto&up=left_upwind?face_l:face_r;
 const double up_h=face_thermo[left_upwind?2:3];
 FluxJacobians result;
 for(int side=0;side<2;++side){
  auto&d=side?result.right:result.left;
  const bool own_side=left_upwind==(side==0);
  const auto&cell=side?cell_r:cell_l;const auto&t=cell_thermo[side];
  const double dm_p=own_side?volume_flux/(t.rg*cell[4]):0;
  const double dm_t=own_side?-volume_flux*t.rho/cell[4]:0;
  D3 dm_u{};for(int j=0;j<3;++j)dm_u[j]=.5*rho*normal[j];
  d[0]=dm_p;d[4]=dm_t;for(int j=0;j<3;++j)d[1+j]=dm_u[j];
  for(int i=0;i<3;++i){
   d[5*(i+1)]=up[i+1]*dm_p+.5*normal[i];
   for(int j=0;j<3;++j)d[5*(i+1)+j+1]=up[i+1]*dm_u[j];
   if(own_side)d[6*(i+1)]+=rho*volume_flux;
   d[5*(i+1)+4]=up[i+1]*dm_t;
  }
  d[20]=up_h*dm_p;d[24]=up_h*dm_t+(own_side?rho*volume_flux*t.cp:0);
  for(int j=0;j<3;++j)d[21+j]=up_h*dm_u[j]+(own_side?rho*volume_flux*up[j+1]:0);
 }
 const bool positive_dissipation=dissipation_mass>0;
 const auto&d_up=positive_dissipation?face_l:face_r;
 const double d_h=face_thermo[positive_dissipation?2:3];
 const double velocity_factor=rho*difference;
 const double pressure_velocity_factor=rho*(average-absolute+alpha*volume_flux*difference);
 B5 dissipation{};dissipation[0]=pressure_coefficient;
 for(int j=0;j<3;++j)dissipation[1+j]=velocity_factor*normal[j];
 for(int i=0;i<3;++i){
  dissipation[5*(i+1)]=d_up[i+1]*pressure_coefficient+normal[i]*difference;
  for(int j=0;j<3;++j)dissipation[5*(i+1)+j+1]=(d_up[i+1]*velocity_factor+pressure_velocity_factor*normal[i])*normal[j];
 }
 dissipation[20]=d_h*pressure_coefficient;
 for(int j=0;j<3;++j)dissipation[21+j]=d_h*velocity_factor*normal[j];
 for(int j=0;j<25;++j){result.left[j]+=.5*dissipation[j];result.right[j]-=.5*dissipation[j];}
 auto&d=positive_dissipation?result.left:result.right;
 const auto&cell=positive_dissipation?cell_l:cell_r;
 for(int j=0;j<3;++j){d[6*(j+1)]+=dissipation_mass;d[21+j]+=dissipation_mass*cell[j+1];}
 d[24]+=dissipation_mass*cell_thermo[positive_dissipation?0:1].cp;
 return result;
}

inline FluxJacobians conservative_face_jacobians(const Gas&gas,const W7&cell_l,const W7&cell_r,const FluxJacobians&primitive){
 FluxJacobians result;
 for(int side=0;side<2;++side){
  const auto&cell=side?cell_r:cell_l;const auto&source=side?primitive.right:primitive.left;
  auto&target=side?result.right:result.left;auto dp=gas.pressure_jacobian(cell);auto dw=gas.primitive_jacobian(cell);
  for(int row=0;row<5;++row)for(int col=0;col<5;++col)for(int k=0;k<5;++k)target[5*row+col]+=source[5*row+k]*(k==0?dp[col]:dw[5*k+col]);
 }
 return result;
}

inline FluxJacobians coupled_approximate_jacobian(const Gas&gas,const W7&cell_l,const W7&cell_r,
 const W7&face_l,const W7&face_r,const D3&normal,double mass,
 double a2_l,double a2_r,double ref_l,double ref_r,const std::array<double,4>*face_thermo=nullptr){
 const std::array<double,4>thermo=face_thermo?*face_thermo:std::array<double,4>{gas.rho(face_l),gas.rho(face_r),gas.enthalpy(face_l),gas.enthalpy(face_r)};
 const std::array<FaceCellThermodynamics,2>cell_thermo={FaceCellThermodynamics{gas.rho(cell_l),gas.cp(cell_l[4]),gas.rg},FaceCellThermodynamics{gas.rho(cell_r),gas.cp(cell_r[4]),gas.rg}};
 auto primitive=approximate_coupled_primitive_jacobian(cell_l,cell_r,face_l,face_r,normal,mass,thermo,a2_l,a2_r,ref_l,ref_r,cell_thermo);
 return conservative_face_jacobians(gas,cell_l,cell_r,primitive);
}
}
