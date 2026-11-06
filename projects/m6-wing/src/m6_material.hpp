#pragma once
// Original-Case species polynomials. Float property evaluation and double
// integration are independent implementations of the observed public format.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>

namespace m6 {
struct MaterialPiece {
 float minimum=0,maximum=0,minimum_cp=0,maximum_cp=0;
 std::array<float,5>coefficient{};
};
struct MaterialSpecies {
 double fraction=0,molecular_weight=0;std::array<MaterialPiece,2>piece;
};
class MixtureMaterial {
 std::array<MaterialSpecies,3>species;
 template<class T>static void read(std::ifstream&f,T*value,size_t count){
  f.read(reinterpret_cast<char*>(value),static_cast<std::streamsize>(sizeof(T)*count));
  if(!f)throw std::runtime_error("truncated original material input");
 }
 // An observable float intermediate prevents a compiler from fusing the
 // multiply and add on a later AVX/FMA Linux build.
 static float unfused(float a,float b,float c){volatile float product=a*b;return product+c;}
 static float evaluate(const MaterialPiece&p,float temperature){
  float value=p.coefficient[4];for(int k=3;k>=0;--k)value=unfused(value,temperature,p.coefficient[k]);return value;
 }
 static double integral(const MaterialPiece&p,double a,double b){
  double lower=static_cast<double>(p.coefficient[4])/5,upper=lower;
  for(int k=3;k>=0;--k){double coefficient=static_cast<double>(p.coefficient[k])/(k+1);lower=lower*a+coefficient;upper=upper*b+coefficient;}
  return upper*b-lower*a;
 }
 static double species_h(const MaterialSpecies&s,double temperature){
  const auto&first=s.piece.front();const auto&last=s.piece.back();
  double result=static_cast<double>(first.minimum_cp)*std::min(temperature,static_cast<double>(first.minimum));
  for(const auto&p:s.piece)if(temperature>p.minimum)result+=integral(p,p.minimum,std::min(temperature,static_cast<double>(p.maximum)));
  if(temperature>last.maximum)result+=(temperature-last.maximum)*static_cast<double>(last.maximum_cp);return result;
 }
 public:
 explicit MixtureMaterial(const std::string&path){
  std::ifstream input(path,std::ios::binary);char magic[8];std::array<int32_t,4>header{};read(input,magic,8);read(input,header.data(),header.size());
  if(std::string(magic,8)!="M6MCPP01"||header!=std::array<int32_t,4>{1,3,2,5})throw std::runtime_error("unsupported original material format");
  double total=0;
  for(auto&s:species){
   read(input,&s.fraction,1);read(input,&s.molecular_weight,1);
   if(!(s.fraction>=0&&std::isfinite(s.fraction)&&s.molecular_weight>0&&std::isfinite(s.molecular_weight)))throw std::runtime_error("invalid original species values");total+=s.fraction;
   for(auto&p:s.piece){std::array<float,4>bounds{};read(input,bounds.data(),bounds.size());p.minimum=bounds[0];p.maximum=bounds[1];p.minimum_cp=bounds[2];p.maximum_cp=bounds[3];read(input,p.coefficient.data(),p.coefficient.size());
    if(!(p.minimum>0&&p.maximum>p.minimum&&p.minimum_cp>0&&p.maximum_cp>0))throw std::runtime_error("invalid material polynomial range");
    for(float value:p.coefficient)if(!std::isfinite(value))throw std::runtime_error("nonfinite material coefficient");
    if(p.minimum_cp!=evaluate(p,p.minimum)||p.maximum_cp!=evaluate(p,p.maximum))throw std::runtime_error("material endpoint precision mismatch");
   }
   if(s.piece[0].maximum!=s.piece[1].minimum)throw std::runtime_error("noncontiguous material ranges");
  }
  if(std::abs(total-1)>1e-12||input.peek()!=std::char_traits<char>::eof())throw std::runtime_error("invalid material mass fractions or trailing bytes");
 }
 double cp(double temperature,const double*fractions=nullptr)const{
  float t=static_cast<float>(temperature);double result=0;
  for(size_t i=0;i<species.size();++i){const auto&s=species[i];float value;
   if(t<=s.piece.front().minimum)value=s.piece.front().minimum_cp;
   else if(t>=s.piece.back().maximum)value=s.piece.back().maximum_cp;
   else{const auto&p=s.piece[t<s.piece[0].maximum?0:1];value=evaluate(p,t);}
   result+=(fractions?fractions[i]:s.fraction)*static_cast<double>(value);
  }return result;
 }
 double h(double temperature,const double*fractions=nullptr)const{
  double result=0;for(size_t i=0;i<species.size();++i)result+=(fractions?fractions[i]:species[i].fraction)*species_h(species[i],temperature);return result;
 }
 double cp_derivative(double temperature)const{
  float t=static_cast<float>(temperature);double result=0;
  for(const auto&s:species){if(t<s.piece.front().minimum||t>=s.piece.back().maximum)continue;
   const auto&p=s.piece[t<s.piece[0].maximum?0:1];float value=4*p.coefficient[4];
   for(int k=3;k>=1;--k)value=unfused(value,t,static_cast<float>(k)*p.coefficient[k]);result+=s.fraction*static_cast<double>(value);
  }return result;
 }
};
}
