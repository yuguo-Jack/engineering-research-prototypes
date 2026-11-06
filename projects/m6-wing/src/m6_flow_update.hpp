#pragma once
#include "m6_physics.hpp"
#include <cstdint>
#include <fstream>

namespace m6 {
// Physical limits parsed from the original Case; separate from numerical
// forcing tolerances and from the legacy fractional update safeguard.
struct FlowUpdateLimits {
 std::array<double,5>physical{}; // min/max pressure, min/max temperature, min density
 void load(const std::string&path){
  std::ifstream input(path,std::ios::binary);char magic[8]{};uint32_t version=0;
  input.read(magic,8);input.read(reinterpret_cast<char*>(&version),4);input.read(reinterpret_cast<char*>(physical.data()),40);
  if(!input||std::string(magic,8)!="M6FLIM01"||version!=1||input.peek()!=std::char_traits<char>::eof())throw std::runtime_error("invalid original flow-update limits input");
  for(double value:physical)if(!std::isfinite(value)||value<=0)throw std::runtime_error("nonpositive original flow-update limit");
  if(physical[0]>=physical[1]||physical[2]>=physical[3])throw std::runtime_error("unordered original flow-update limits");
 }
};
inline bool within_case_flow_update_limits(const Gas&gas,const W7&state,const FlowUpdateLimits&limits){
 for(int j=0;j<5;++j)if(!std::isfinite(state[j]))return false;
 const double pressure=state[0]*gas.rho0*gas.a0*gas.a0,temperature=state[4]*gas.T0,rho=gas.rho(state)*gas.rho0;
 return pressure>=limits.physical[0]&&pressure<=limits.physical[1]&&temperature>=limits.physical[2]&&temperature<=limits.physical[3]&&std::isfinite(rho)&&rho>=limits.physical[4];
}
// Linear change of variables for an independently solved conserved correction.
// This is distinct from the nonlinear Q -> W inversion after a finite update.
inline Q5 linearized_primitive_increment(const Gas&gas,const W7&state,const Q5&dq){
 const double rho=gas.rho(state),dr=dq[0];auto q=gas.conserved(state);Q5 dw{};
 double de=(dq[4]-q[4]/rho*dr)/rho;
 for(int j=1;j<=3;++j){dw[j]=(dq[j]-state[j]*dr)/rho;de-=state[j]*dw[j];}
 dw[4]=de/(gas.cp(state[4])-gas.rg);dw[0]=gas.rg*(state[4]*dr+rho*dw[4]);return dw;
}
inline W7 primitive_flow_trial(const W7&state,const Q5&dw,double relaxation){
 W7 result=state;for(int j=0;j<5;++j)result[j]+=relaxation*dw[j];return result;
}
}
