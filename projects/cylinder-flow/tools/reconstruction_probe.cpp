// Diagnostic operator comparison only: no flow iteration or solution is produced.
#define main cylinder_driver_main
#include "../src/cylinder.cpp"
#undef main
#include <sstream>

int main(int argc,char**argv){try{
 if(argc<4)throw std::runtime_error("usage: reconstruction_probe mesh.bin native-cell-probe.csv output.csv [limiter]");
 Solver solver(argv[1]);solver.limiter=argc>4?argv[4]:"barth";
 std::ifstream input(argv[2]);std::string line;std::getline(input,line);std::vector<std::vector<double>>native(solver.mesh.n);
 while(std::getline(input,line)){std::stringstream stream(line);std::vector<double>row;std::string field;while(std::getline(stream,field,','))row.push_back(std::stod(field));int id=static_cast<int>(row[1])-1;
  if(id<0||id>=solver.mesh.n)throw std::runtime_error("invalid reference cell id");
  native[id]=row;V w{row[4]/solver.mesh.rho0,row[5]/solver.mesh.a0,row[6]/solver.mesh.a0,row[7]/(solver.mesh.rho0*solver.mesh.a0*solver.mesh.a0)};solver.q[id]=solver.mesh.gas.conserved(w);
 }
 for(const auto&row:native)if(row.empty())throw std::runtime_error("incomplete diagnostic field");
 solver.gradients();solver.update_reference_velocity();std::ofstream out(argv[3]);out<<"cell,x,y,p_gx,p_gy,p_gx_native,p_gy_native,u_gx,u_gy,u_gx_native,u_gy_native,v_gx,v_gy,v_gx_native,v_gy_native,uref2,uref2_native\n"<<std::setprecision(16);
 double maxerr=0;
 for(int i=0;i<solver.mesh.n;++i){const auto&r=native[i];double ps=solver.mesh.rho0*solver.mesh.a0*solver.mesh.a0;out<<i+1<<','<<solver.mesh.center[i][0]<<','<<solver.mesh.center[i][1]<<','<<solver.gx[i][3]*ps<<','<<solver.gy[i][3]*ps<<','<<r[11]<<','<<r[12]<<','<<solver.gx[i][1]*solver.mesh.a0<<','<<solver.gy[i][1]*solver.mesh.a0<<','<<r[13]<<','<<r[14]<<','<<solver.gx[i][2]*solver.mesh.a0<<','<<solver.gy[i][2]*solver.mesh.a0<<','<<r[15]<<','<<r[16]<<','<<solver.uref2[i]*solver.mesh.a0*solver.mesh.a0<<','<<r[17]<<'\n';
  maxerr=std::max(maxerr,std::abs(solver.uref2[i]*solver.mesh.a0*solver.mesh.a0/r[17]-1));}
 std::cout<<"Reference-velocity operator max relative difference="<<maxerr<<"; diagnostic only, no flow solve\n";return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
