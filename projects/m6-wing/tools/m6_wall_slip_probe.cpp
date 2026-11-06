// Same compiled wall virtual-velocity law used by the independent candidate.
#include "../src/m6_wall_slip.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
int main(int argc,char**argv){try{
 if(argc==2&&std::string(argv[1])=="--self-test"){
  m6::D3 normal={.3,-.4,std::sqrt(.75)},velocity={100,20,-10};double worst=0;int checks=0;
  for(double distance:{1e-16,1e-10,1e-6,1e-4,1e-2}){
   auto base=m6::stationary_wall_virtual_state(1.2,1.72e-5,16,velocity,normal,distance);
   for(int column=0;column<4;++column){auto plus=velocity,minus=velocity;double rp=1.2,rm=1.2,eps=1e-6*std::max(1.,column<3?std::abs(velocity[column]):1.2);
    if(column<3){plus[column]+=eps;minus[column]-=eps;}else{rp+=eps;rm-=eps;}
    auto p=m6::stationary_wall_virtual_state(rp,1.72e-5,16,plus,normal,distance),m=m6::stationary_wall_virtual_state(rm,1.72e-5,16,minus,normal,distance);
    for(int row=0;row<3;++row){double fd=(p.velocity[row]-m.velocity[row])/(2*eps),analytic=base.jacobian[4*row+column];worst=std::max(worst,std::abs(fd-analytic)/std::max(1.,std::abs(fd)+std::abs(analytic)));++checks;}
   }
  }
  if(worst>3e-7)throw std::runtime_error("wall virtual-velocity analytic derivative check failed");std::cout<<"M6_WALL_SLIP_JACOBIAN checks="<<checks<<" maximum_scaled_error="<<worst<<'\n';return 0;
 }
 if(argc!=3)throw std::runtime_error("m6_wall_slip_probe INPUT.csv OUTPUT.csv");std::ifstream input(argv[1]);std::ofstream output(argv[2]);if(!input||!output)throw std::runtime_error("cannot open wall slip files");
 std::string line,token;std::getline(input,line);std::istringstream header(line);std::map<std::string,size_t>columns;size_t count=0;while(std::getline(header,token,','))columns[token]=count++;
 output<<"epoch,sequence,iteration,face,cell,us0,us1,us2,yplus,ustar_from_k,uks0,uks1,uks2\n"<<std::setprecision(17);size_t rows=0;
 while(std::getline(input,line)){if(line.empty())continue;std::istringstream row(line);std::vector<double>v;while(std::getline(row,token,','))v.push_back(std::stod(token));if(v.size()!=columns.size())throw std::runtime_error("wrong wall slip input width");for(double x:v)if(!std::isfinite(x))throw std::runtime_error("nonfinite wall slip input");auto get=[&](const std::string&key){return v.at(columns.at(key));};
  for(auto key:{"ks","vwPlus","alpha","beta","gamma","les"})if(get(key)!=0)throw std::runtime_error("wall slip probe only covers the observed smooth zero-gradient branch");
  m6::D3 velocity={get("u"),get("v"),get("w")},wall{},tangent{};for(int j=0;j<3;++j){wall[j]=get("uw"+std::to_string(j));tangent[j]=get("Vp"+std::to_string(j));}
  auto result=m6::wall_virtual_velocity_from_star(get("uStar"),get("yp"),get("rho"),get("mu"),velocity,wall,tangent);
  double from_k=m6::wall_friction_velocity_from_k(get("rho"),get("mu"),get("k"),get("up"),get("yp"));auto k_result=m6::wall_virtual_velocity_from_star(from_k,get("yp"),get("rho"),get("mu"),velocity,wall,tangent);
  for(auto key:{"epoch","sequence","iteration","face","cell"})output<<get(key)<<',';for(double x:result)output<<x<<',';output<<get("uStar")*get("yp")*get("rho")/get("mu")<<','<<from_k;for(double x:k_result)output<<','<<x;output<<'\n';++rows;
 }if(!rows)throw std::runtime_error("no wall slip observations");std::cout<<"M6_WALL_SLIP_PROBE rows="<<rows<<" integration_steps=0\n";return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}}
