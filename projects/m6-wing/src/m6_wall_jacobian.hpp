#pragma once
#include "m6_physics.hpp"

namespace m6 {
// The normal stationary-wall implicit block freezes reconstructed face density
// and reference speed. It is an approximate solver matrix, separate from the
// derivative of the physical pressure p_face + rho_face * Uref * u_normal.
inline Q5 normal_wall_pressure_jacobian(const Gas&gas,const W7&cell,const D3&normal,double reference2,double face_density){
 auto primitive=gas.primitive_jacobian(cell);auto result=gas.pressure_jacobian(cell);
 const double coefficient=face_density*std::sqrt(reference2);
 for(int j=0;j<5;++j){double normal_velocity_change=0;for(int k=0;k<3;++k)normal_velocity_change+=normal[k]*primitive[5*(k+1)+j];result[j]+=coefficient*normal_velocity_change;}
 return result;
}
}
