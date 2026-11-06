#pragma once
// Pure derivatives captured from a fixed cell W.  These are not reconstructed
// face states, SST coefficients, boundary caches, or an approximate derivative.
#include "m6_approximate_jacobian.hpp"

namespace m6 {
struct CellConservativeJacobian {
 Q5 pressure;
 B5 primitive;
#ifdef M6_CELL_JACOBIAN_THERMO_CACHE
 double density;
 double heat_capacity;
#endif
};
inline CellConservativeJacobian capture_cell_conservative_jacobian(const Gas&gas,const W7&cell){
 // Keep both original evaluations, even though pressure_jacobian internally
 // repeats primitive_jacobian.  No algebra or arithmetic order is changed.
 return {gas.pressure_jacobian(cell),gas.primitive_jacobian(cell)
#ifdef M6_CELL_JACOBIAN_THERMO_CACHE
  ,gas.rho(cell),gas.cp(cell[4])
#endif
 };
}
inline FluxJacobians conservative_face_jacobians_cached(
 const CellConservativeJacobian&left,const CellConservativeJacobian&right,
 const FluxJacobians&primitive){
 FluxJacobians result;
 for(int side=0;side<2;++side){
  const auto&source=side?primitive.right:primitive.left;
  auto&target=side?result.right:result.left;const auto&cached=side?right:left;
  const auto&dp=cached.pressure;const auto&dw=cached.primitive;
  // Identical row/column/k traversal and accumulation to the original function.
  for(int row=0;row<5;++row)for(int col=0;col<5;++col)for(int k=0;k<5;++k)target[5*row+col]+=source[5*row+k]*(k==0?dp[col]:dw[5*k+col]);
 }
 return result;
}
inline FluxJacobians coupled_approximate_jacobian_cached(const Gas&gas,const W7&cell_l,const W7&cell_r,
 const W7&face_l,const W7&face_r,const D3&normal,double mass,
 double a2_l,double a2_r,double ref_l,double ref_r,
 const CellConservativeJacobian&cached_l,const CellConservativeJacobian&cached_r,
 const std::array<double,4>*face_thermo=nullptr){
 const std::array<double,4>thermo=face_thermo?*face_thermo:std::array<double,4>{gas.rho(face_l),gas.rho(face_r),gas.enthalpy(face_l),gas.enthalpy(face_r)};
#ifdef M6_CELL_JACOBIAN_THERMO_CACHE
 // The existing owner captures these values from the same immutable cell W
 // once per assemble(matrix=true), including synchronized ghost cells. Face
 // reconstruction, pressure/primitive Jacobians and product order are intact.
 const std::array<FaceCellThermodynamics,2>cell_thermo={FaceCellThermodynamics{cached_l.density,cached_l.heat_capacity,gas.rg},FaceCellThermodynamics{cached_r.density,cached_r.heat_capacity,gas.rg}};
#else
 const std::array<FaceCellThermodynamics,2>cell_thermo={FaceCellThermodynamics{gas.rho(cell_l),gas.cp(cell_l[4]),gas.rg},FaceCellThermodynamics{gas.rho(cell_r),gas.cp(cell_r[4]),gas.rg}};
#endif
 auto primitive=approximate_coupled_primitive_jacobian(cell_l,cell_r,face_l,face_r,normal,mass,thermo,a2_l,a2_r,ref_l,ref_r,cell_thermo);
 return conservative_face_jacobians_cached(cached_l,cached_r,primitive);
}
}
