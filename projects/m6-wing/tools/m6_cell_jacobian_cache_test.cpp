#include "../src/m6_cell_jacobian_cache.hpp"
#include <cstring>
#include <iomanip>
#include <iostream>
#include <vector>

template<size_t N>void equal_bits(const std::array<double,N>&a,const std::array<double,N>&b,const char*what){
 for(size_t i=0;i<N;++i){if(!std::isfinite(a[i])||!std::isfinite(b[i]))throw std::runtime_error(std::string(what)+" nonfinite");
  if(std::memcmp(&a[i],&b[i],sizeof(double))!=0)throw std::runtime_error(std::string(what)+" bitwise mismatch at "+std::to_string(i));}
}

int main(int argc,char**argv){try{
 std::array<double,48>meta{};meta[0]=1.1361339060496805;meta[1]=352.517926811456;meta[2]=298.15;meta[3]=299.124701264874;meta[4]=101325;meta[5]=1.72e-5;meta[6]=.0454;meta[7]=1059.4853515625;
 std::array<double,5>low={985.1700244,.4572524926,-.0010998299666,1.5381987004e-6,-6.675425469e-10},high={901.5941633,.4449195144,-.00015867025232,2.7184049518e-8,-1.7782013004e-12};for(int i=0;i<5;++i){meta[8+i]=low[i];meta[13+i]=high[i];}
 m6::Gas gas(meta);if(argc==2)gas.load_material(argv[1]);else if(argc!=1)throw std::runtime_error("usage: m6_cell_jacobian_cache_test [original material-polynomials.bin]");
 std::vector<double>temperatures={200.,249.,298.15,299.999,300.,300.001,700.,999.999,1000.,1000.001,1200.,4800.,5000.,5000.001,5500.};
 size_t states=0,pairs=0;double checksum=0;
 for(double temperature:temperatures)for(int sample=0;sample<19;++sample){
  m6::W7 left={.35+.025*sample,.04*sample-.31,.019*sample-.17,.05-.011*sample,temperature/meta[2],1e-4,1.};
  m6::W7 right=left;right[0]*=1.013;right[1]-=.017;right[2]+=.009;right[4]*=1.0003;
  auto lc=m6::capture_cell_conservative_jacobian(gas,left),rc=m6::capture_cell_conservative_jacobian(gas,right);
#ifdef M6_CELL_JACOBIAN_THERMO_CACHE
  equal_bits(std::array<double,2>{lc.density,lc.heat_capacity},std::array<double,2>{gas.rho(left),gas.cp(left[4])},"left captured cell thermodynamics");
  equal_bits(std::array<double,2>{rc.density,rc.heat_capacity},std::array<double,2>{gas.rho(right),gas.cp(right[4])},"right captured cell thermodynamics");
#endif
  equal_bits(lc.primitive,gas.primitive_jacobian(left),"left primitive");equal_bits(lc.pressure,gas.pressure_jacobian(left),"left pressure");
  equal_bits(rc.primitive,gas.primitive_jacobian(right),"right primitive");equal_bits(rc.pressure,gas.pressure_jacobian(right),"right pressure");states+=2;
  m6::FluxJacobians seeded;for(int j=0;j<25;++j){seeded.left[j]=std::sin(.37*(sample+j+1));seeded.right[j]=std::cos(.29*(sample+2*j+1));}
  auto old_transform=m6::conservative_face_jacobians(gas,left,right,seeded),cached_transform=m6::conservative_face_jacobians_cached(lc,rc,seeded);
  equal_bits(old_transform.left,cached_transform.left,"left conservative transform");equal_bits(old_transform.right,cached_transform.right,"right conservative transform");
  for(double mass:{-.37,0.,.41})for(m6::D3 normal:std::array<m6::D3,3>{{{1,0,0},{0,-1,0},{.3,-.4,std::sqrt(.75)}}}){
   auto fl=left,fr=right;fl[0]*=1.0002;fl[2]+=.003;fr[3]-=.004;fr[4]*=.9997;
   const double al=std::pow(gas.sound(left),2),ar=std::pow(gas.sound(right),2),rl=.21*al,rr=.22*ar;
   std::array<double,4>thermo={gas.rho(fl)+.001,gas.rho(fr)+.002,gas.enthalpy(fl)+.004,gas.enthalpy(fr)-.003};
   for(bool independent:{false,true}){
    const auto*input=independent?&thermo:nullptr;
    auto old=m6::coupled_approximate_jacobian(gas,left,right,fl,fr,normal,mass,al,ar,rl,rr,input);
    auto cached=m6::coupled_approximate_jacobian_cached(gas,left,right,fl,fr,normal,mass,al,ar,rl,rr,lc,rc,input);
    equal_bits(old.left,cached.left,"left coupled approximate block");equal_bits(old.right,cached.right,"right coupled approximate block");checksum+=cached.left[0]+cached.right[24];++pairs;
   }
  }
 }
 std::cout<<std::setprecision(17)<<"M6_CELL_JACOBIAN_CACHE_BITWISE states="<<states<<" face_pairs="<<pairs<<" cache_bytes_per_cell="<<sizeof(m6::CellConservativeJacobian)<<" material="<<(argc==2?"original-case":"legacy-polynomial")<<" checksum="<<checksum<<" passed\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
