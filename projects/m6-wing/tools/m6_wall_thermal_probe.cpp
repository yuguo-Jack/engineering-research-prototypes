#include <limits>
#include "../src/m6_wall_thermal.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
int main(int argc,char**argv){try{
 if(argc==2&&std::string(argv[1])=="--self-test"){
  int checks=0;double derivative_error=0;
  for(double y:{.0001,.01,.4,1.,5.,30.,1000.})for(double pr:{.1,.4,.85,2.,10.}){
   auto a=m6::smooth_wall_thermal(y,8,1.2,1060,pr,.85,100),b=m6::smooth_wall_thermal(y,8,2.4,1060,pr,.85,100),c=m6::smooth_wall_thermal(y,8,1.2,2120,pr,.85,100);
   if(!(a.virtual_fraction>=0&&a.virtual_fraction<=1&&std::isfinite(a.adiabatic_jump)))throw std::runtime_error("thermal wall physical bounds failed");
   if(std::abs(b.adiabatic_jump-a.adiabatic_jump)>1e-14||std::abs(c.adiabatic_jump*2-a.adiabatic_jump)>1e-14)throw std::runtime_error("thermal wall density/Cp scaling failed");
   auto lo=m6::smooth_wall_thermal(y,8,1.2,1060,.85*(1-1e-7),.85,100),hi=m6::smooth_wall_thermal(y,8,1.2,1060,.85*(1+1e-7),.85,100);
   if(std::abs(hi.adiabatic_jump-lo.adiabatic_jump)>1e-5)throw std::runtime_error("thermal Pr/Prt unity limit failed");++checks;
   for(int parameter=0;parameter<2;++parameter){double step=1e-6*(parameter==0?y:pr);auto plus=m6::smooth_wall_thermal(y+(parameter==0?step:0),8,1.2,1060,pr+(parameter==1?step:0),.85,100),minus=m6::smooth_wall_thermal(y-(parameter==0?step:0),8,1.2,1060,pr-(parameter==1?step:0),.85,100);double fd=(plus.virtual_fraction-minus.virtual_fraction)/(2*step),analytic=parameter==0?a.fraction_y_derivative:a.fraction_pr_derivative;derivative_error=std::max(derivative_error,std::abs(fd-analytic)/std::max(1.,std::abs(fd)+std::abs(analytic)));}
  }
  if(derivative_error>3e-7)throw std::runtime_error("thermal virtual-temperature derivative check failed");std::cout<<"M6_WALL_THERMAL_CORE physical_cases="<<checks<<" derivative_max_scaled_error="<<derivative_error<<'\n';return 0;
 }
 if(argc!=3)throw std::runtime_error("m6_wall_thermal_probe INPUT.csv OUTPUT.csv");
 std::ifstream input(argv[1]);std::ofstream output(argv[2]);if(!input||!output)throw std::runtime_error("cannot open thermal probe files");
 std::string line,token;std::getline(input,line);std::istringstream header(line);std::map<std::string,size_t>columns;size_t count=0;while(std::getline(header,token,','))columns[token]=count++;
 if(columns.count("Tw")!=columns.count("Tc"))throw std::runtime_error("both wall and cell temperatures are required");
 output<<"sample,sequence,iteration,a,b,adiabatic_jump,tplus,slope,virtual_fraction,virtual_T,adiabatic_wall_T,heat_source,heat_conductance\n"<<std::setprecision(17);size_t rows=0;
 while(std::getline(input,line)){if(line.empty())continue;std::istringstream row(line);std::vector<double>v;while(std::getline(row,token,','))v.push_back(std::stod(token));if(v.size()!=columns.size())throw std::runtime_error("wrong thermal input width");for(double x:v)if(!std::isfinite(x))throw std::runtime_error("nonfinite thermal input");auto get=[&](const std::string&key){return v.at(columns.at(key));};
  auto result=m6::smooth_wall_thermal(get("yplus"),get("ustar"),get("rho"),get("cp"),get("Pr"),get("Prt"),get("up"),get("rough_wall")!=0);
  double virtual_temperature=columns.count("Tw")?get("Tc")+(get("Tw")-get("Tc"))*result.virtual_fraction:0;
  double wall_temperature=0,source=0,conductance=0;if(columns.count("distance")){wall_temperature=m6::smooth_adiabatic_wall_temperature(result,get("Tc"),get("distance"),get("conductivity"),get("conductivity_epsilon"));double effective=std::max(get("conductivity"),get("distance")*result.conductance);source=-result.viscous_source*effective/(effective+get("conductivity_epsilon"));conductance=effective/get("distance");}
  output<<get("sample")<<','<<get("sequence")<<','<<get("iteration")<<','<<result.viscous_source<<','<<result.conductance<<','<<result.adiabatic_jump<<','<<result.tplus<<','<<result.slope<<','<<result.virtual_fraction<<','<<virtual_temperature<<','<<wall_temperature<<','<<source<<','<<conductance<<'\n';++rows;
 }
 if(!rows||!output)throw std::runtime_error("no thermal observations or output failure");std::cout<<"M6_WALL_THERMAL_PROBE rows="<<rows<<" integration_steps=0\n";return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}}
