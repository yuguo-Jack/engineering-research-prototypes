// Independent 3-D finite-volume RANS/SST solver; standard C++17 + MPI only.
#define NOMINMAX
#include "m6_linear.hpp"
#include "m6_acceleration.hpp"
#include "m6_flux_jacobian.hpp"
#include "m6_approximate_jacobian.hpp"
#ifdef M6_CELL_JACOBIAN_CACHE
#include "m6_cell_jacobian_cache.hpp"
#endif
#include "m6_farfield.hpp"
#include "m6_farfield_jacobian.hpp"
#include "m6_wall_jacobian.hpp"
#include "m6_viscous.hpp"
#include "m6_wall_slip.hpp"
#include "m6_wall_thermal.hpp"
#include "m6_sst_transport.hpp"
#include "m6_pseudotime.hpp"
#include "m6_scalar_time.hpp"
#include "m6_scalar_update.hpp"
#include "m6_flow_update.hpp"
#include <memory>
#ifdef _WIN32
#include <windows.h>
#else
#include <sched.h>
#include <unistd.h>
#endif
#include <chrono>
#include <deque>
#include <map>
#include <unordered_map>

namespace m6 {
using Gradient=std::array<double,21>;
using Properties=std::array<double,8>; // mut,F1,beta,gamma,Pk,CD,strain2,vorticity
void apply_cpu_affinity(const MPIWorld&world,const std::string&cpu_list){
 if(!cpu_list.empty()){std::vector<int>cpus;std::istringstream s(cpu_list);std::string token;while(std::getline(s,token,','))cpus.push_back(std::stoi(token));if(static_cast<int>(cpus.size())!=world.size)throw std::runtime_error("one CPU index is required per MPI rank");int cpu=cpus[world.rank];
#ifdef _WIN32
   if(cpu<0||cpu>=64||!SetProcessAffinityMask(GetCurrentProcess(),static_cast<DWORD_PTR>(1ULL<<cpu)))throw std::runtime_error("cannot set rank CPU affinity");
#else
   cpu_set_t mask;CPU_ZERO(&mask);CPU_SET(cpu,&mask);if(sched_setaffinity(0,sizeof(mask),&mask)!=0)throw std::runtime_error("cannot set rank CPU affinity");
#endif
   for(int rank=0;rank<world.size;++rank){MPI_Barrier(MPI_COMM_WORLD);if(rank==world.rank)std::cout<<"M6_AFFINITY rank="<<world.rank<<" logical_cpu="<<cpu<<std::endl;}
  }
}
struct Options {
 std::string input,output="m6-output",restart,cpu_list;int iterations=20,order=2,linear_steps=8,jacobian_period=1,checkpoint=50,anderson_depth=0,anderson_start=10,linear_audit_period=0;
 bool recycle_linear_direction=false;
 bool export_linear_operator=false;
 bool turbulence_k_first=false;
 bool characteristic_farfield=false;
 bool capture_window_start=false;
 bool wall_slip_gradients=false;
 bool coupled_wall_reconstruction=false;
 bool symmetry_normal_lsq=false;
 bool enhanced_thermal_wall=false;
 bool symmetry_face_gradients=false;
 bool audit_turbulence_linear=false;
 bool coupled_farfield_reconstruction=false;
 bool symmetry_viscous_flux=false;
 bool farfield_zero_viscous_flux=false;
 bool normal_sst_transport=false;
 bool lagged_sst_mass_flux=false;
 bool normal_coupled_timestep=false;
 bool normal_scalar_step=false;
 bool cache_flow_diffusion=false;
 bool case_material_polynomials=false;
 bool normal_approximate_jacobian=false;
 bool normal_farfield_jacobian=false;
 bool normal_wall_jacobian=false;
 bool normal_symmetry_jacobian=false;
 int turbulence_linear_max_steps=0;
 double turbulence_linear_tolerance=.1;
 bool primitive_flow_update=false;
 bool case_flow_update_limits=false;
 bool skip_unused_face_work=false;
 // Coupled native reports retain the reference velocity from Update_Time_Step.
 // This switch changes force reporting only; the next assembly keeps its fresh cache.
 bool lagged_force_reference_velocity=false;
 bool coupled_force_report_stage=false;
 double diagnostic_directional_step=1e-4;
 // coeff-kw-mu-dok=.001 is read from the supplied M6 Case, not a tuned value.
 double cfl=50,linear_tolerance=.1,flow_relax=1,turb_relax=.6,nonorth_limit=-1,omega_laminar_factor=.001,cross_diffusion_limit=1e6;bool validate=false,wall_production=true,tke_pressure=false,implicit_gradients=false,check_viscous=false,cell_omega_wall=true,steady_turbulence=false,sst_vorticity=false,viscous_gauss=true,reconstruct_enthalpy=false,harmonic_diffusion=true,harmonic_scalars=false,edge_correction=false,edge_implicit=false,coupled_ausm=false,turbulence_gauss=false,local_preconditioning=false,exact_flux_jacobian=false,linearize_negative_cross=true,independent_face_thermo=false,boundary_face_ray=false,batched_orthogonalization=false,check_linear=false;
};
struct WindowSnapshot {std::vector<W7>owned;std::vector<double>wall_temperature;std::vector<std::array<double,3>>scalar_boundary;D3 force{};double seconds=0;int iteration=-1;};
struct Solver {
 Mesh mesh;const MPIWorld&world;Gas gas;Options options;Sparse<5>flow;Sparse<1>scalar;FlowUpdateLimits flow_update_limits;
 std::vector<W7>w,lo,hi;Field<5>q,residual;std::vector<Gradient>grad;std::vector<std::array<double,9>>inverse_lsq;
 std::vector<double>scalar_dt;std::array<long long,2>scalar_negative_repairs{},scalar_negative_source_rows{};
 std::vector<std::array<double,3>>flow_diffusion_cell_cache;
 std::vector<std::array<double,8>>flow_diffusion_face_cache;
 bool flow_diffusion_cache_active=false;double flow_diffusion_cache_seconds=0;int flow_diffusion_cache_builds=0;
 std::vector<Properties>prop;std::vector<double>pseudo,massflux,wall_weight,wall_pk;std::vector<std::array<double,12>>physical_gradient;std::vector<std::array<double,6>>scalar_gradient;std::vector<std::array<double,1>>reference_velocity;std::vector<double>total_enthalpy;std::vector<std::array<double,4>>enthalpy_gradient;
 Field<2>thermodynamics,thermo_limiter;Field<6>thermo_gradient;
 W7 case_farfield{};
 std::unordered_map<int,int>wall_gradient_slot;std::vector<std::array<double,12>>wall_virtual_jacobian;
 std::vector<double>wall_temperature;std::vector<std::array<double,5>>wall_temperature_jacobian;
 uint64_t global_thermal_wall_faces=0;
 std::unique_ptr<SSTTransportCache>sst_transport;
 D3 force{},pressure_force{},viscous_force{};std::array<double,7>norms{};
 double linear_estimate=0,linear_true_residual=-1,maximum_linear_true_residual=0;int linear_audits=0,linear_target_missed=0;
 std::array<double,2>turbulence_linear_estimate{},turbulence_linear_true{{-1,-1}},turbulence_linear_maximum{};std::array<int,2>turbulence_linear_steps{},turbulence_linear_audits{},turbulence_linear_misses{};
 Field<5>recycled_flow_direction;KrylovStart linear_start;bool recycle_flow_direction=false;int recycled_attempts=0,recycled_accepted=0;double recycled_initial_ratio_sum=0;
 double gradient_time=0,assembly_time=0,linear_time=0,turbulence_time=0,acceleration_time=0;uint64_t limited_updates=0;int iteration_offset=0,aa_accepted=0,aa_rejected=0;std::array<double,2>scalar_merit{};
 uint64_t window_capture_peak_payload_bytes=0;double window_capture_time=0;
 uint64_t skipped_postflow_assemblies=0;
 Solver(const MPIWorld&mpi,const Options&o):mesh(o.input,mpi),world(mpi),gas(mesh.meta),options(o),flow(mesh,mpi),scalar(mesh,mpi),w(mesh.initial),lo(mesh.nt),hi(mesh.nt),q(mesh.n),residual(mesh.n),grad(mesh.nt),inverse_lsq(mesh.n),prop(mesh.nt),pseudo(mesh.n),massflux(mesh.faces.size()),wall_weight(mesh.n),wall_pk(mesh.n){
  if(options.coupled_wall_reconstruction&&(!options.wall_slip_gradients||!options.boundary_face_ray||!options.viscous_gauss))throw std::runtime_error("coupled wall reconstruction requires wall-slip-gradients, boundary-face-ray and physical Gauss gradients");
  if(options.enhanced_thermal_wall&&!options.coupled_wall_reconstruction)throw std::runtime_error("enhanced thermal wall requires coupled-wall-reconstruction");
  if(options.enhanced_thermal_wall&&options.anderson_depth>0)throw std::runtime_error("enhanced thermal cache has not been validated with Anderson state trials");
  if(options.symmetry_face_gradients&&(!options.symmetry_normal_lsq||!options.viscous_gauss))throw std::runtime_error("symmetry face gradients require normal symmetry LSQ and physical Gauss gradients");
  if(options.coupled_farfield_reconstruction&&(!options.characteristic_farfield||!options.boundary_face_ray||!options.coupled_wall_reconstruction))throw std::runtime_error("coupled farfield reconstruction requires characteristic farfield, face-ray geometry and the coupled reconstruction limiter");
  if(options.symmetry_viscous_flux&&!options.viscous_gauss)throw std::runtime_error("normal symmetry viscous flux requires physical Gauss gradients");
  if(options.farfield_zero_viscous_flux&&!options.characteristic_farfield)throw std::runtime_error("observed zero farfield viscous flux requires characteristic farfield");
  if(options.normal_sst_transport&&(!options.turbulence_k_first||!options.steady_turbulence||!options.cell_omega_wall||!options.turbulence_gauss||!options.coupled_ausm||!options.harmonic_scalars||!options.characteristic_farfield))throw std::runtime_error("normal SST requires k-first steady scalar solves, cell omega walls, Gauss reconstruction, coupled fluxes, harmonic scalar diffusion and characteristic farfield");
  if(options.normal_sst_transport&&options.anderson_depth>0)throw std::runtime_error("normal SST boundary caches have not been validated with Anderson state trials");
  if(options.lagged_sst_mass_flux&&!options.normal_sst_transport)throw std::runtime_error("lagged SST mass flux requires normal SST transport");
  if(options.normal_coupled_timestep&&!options.coupled_ausm)throw std::runtime_error("normal coupled timestep requires coupled flux/reference-velocity storage");
  if(options.lagged_force_reference_velocity&&!options.coupled_ausm)throw std::runtime_error("lagged force reference velocity requires coupled fluxes");
  if(options.coupled_force_report_stage&&(!options.normal_sst_transport||options.anderson_depth>0||options.lagged_force_reference_velocity))throw std::runtime_error("coupled force report stage requires normal SST, no Anderson trials, and no separate lagged-force override");
  if(options.normal_scalar_step&&(!options.normal_sst_transport||!options.lagged_sst_mass_flux||!options.normal_coupled_timestep))throw std::runtime_error("normal scalar step requires normal SST transport and retained coupled mass/time inputs");
  if(options.cache_flow_diffusion&&!options.implicit_gradients)throw std::runtime_error("flow diffusion cache requires the implicit full-gradient operator");
  if(options.normal_approximate_jacobian&&(!options.coupled_ausm||!options.exact_flux_jacobian||!options.independent_face_thermo||!options.local_preconditioning||options.reconstruct_enthalpy))throw std::runtime_error("normal approximate Jacobian requires coupled AUSM, exact boundary derivatives, independent face thermodynamics and local preconditioning");
  if(options.normal_farfield_jacobian&&(!options.characteristic_farfield||!options.exact_flux_jacobian||!options.independent_face_thermo))throw std::runtime_error("normal farfield Jacobian requires characteristic farfield, exact boundary mode and independent face density");
  if(options.normal_wall_jacobian&&(!options.coupled_ausm||!options.exact_flux_jacobian||!options.independent_face_thermo))throw std::runtime_error("normal wall Jacobian requires coupled fluxes, exact boundary mode and independent face density");
  if(options.normal_symmetry_jacobian&&(!options.coupled_ausm||!options.exact_flux_jacobian||!options.independent_face_thermo))throw std::runtime_error("normal symmetry Jacobian requires coupled fluxes, exact boundary mode and independent face density");
  recycle_flow_direction=o.recycle_linear_direction;
  if(options.case_material_polynomials)gas.load_material((std::filesystem::path(o.input)/"material-polynomials.bin").string());
  if(options.case_flow_update_limits){if(!options.primitive_flow_update)throw std::runtime_error("Case flow-update limits require primitive flow updates");flow_update_limits.load((std::filesystem::path(o.input)/"flow-update-limits.bin").string());}
  case_farfield=mesh.initial[0];
  if(options.characteristic_farfield){
   if(mesh.meta[47]!=1)throw std::runtime_error("characteristic far-field needs separately parsed original Case profiles; prepare a new partition");
   case_farfield={mesh.meta[44]/(gas.rho0*gas.a0*gas.a0),0,0,0,mesh.meta[43]/gas.T0,mesh.meta[45]/(gas.a0*gas.a0),mesh.meta[46]/gas.a0};
   D3 direction={mesh.meta[40],mesh.meta[41],mesh.meta[42]};double speed=mesh.meta[39]*gas.sound(case_farfield)/length3(direction);for(int j=0;j<3;++j)case_farfield[j+1]=speed*direction[j];
  }
  for(int i=0;i<mesh.n;++i)q[i]=gas.conserved(w[i]);
  for(auto&f:mesh.faces)if(f.r>=0){f.lr=flow.index(f.l,f.r);if(f.r<mesh.n)f.rl=flow.index(f.r,f.l);}
  if(options.wall_slip_gradients)for(auto&f:mesh.faces)if(f.type==3){wall_gradient_slot.emplace(f.id,static_cast<int>(wall_virtual_jacobian.size()));wall_virtual_jacobian.emplace_back();}
  if(options.enhanced_thermal_wall){wall_temperature.resize(wall_virtual_jacobian.size());wall_temperature_jacobian.resize(wall_temperature.size());for(auto&f:mesh.faces)if(f.type==3)wall_temperature[wall_gradient_slot.at(f.id)]=w[f.l][4];global_thermal_wall_faces=static_cast<uint64_t>(world.sum(static_cast<double>(wall_temperature.size())));}
  if(options.turbulence_gauss)scalar_gradient.resize(mesh.nt);if(options.coupled_ausm)reference_velocity.resize(mesh.nt);if(options.viscous_gauss)physical_gradient.resize(mesh.nt);if(options.reconstruct_enthalpy){total_enthalpy.resize(mesh.nt);enthalpy_gradient.resize(mesh.nt);}if(options.independent_face_thermo){thermodynamics.resize(mesh.nt);thermo_limiter.resize(mesh.nt);thermo_gradient.resize(mesh.nt);}setup_lsq();if(options.normal_sst_transport)sst_transport=std::make_unique<SSTTransportCache>(mesh,gas,world,o.input);if(!o.restart.empty())load(o.restart);
 }
 static double now(){return MPI_Wtime();}
 W7 freestream()const{return options.characteristic_farfield?case_farfield:mesh.initial[0];}
 Wall wall(const Face&f)const{
  const auto&v=w[f.l];double un=v[1]*f.normal[0]+v[2]*f.normal[1]+v[3]*f.normal[2];double u2=0;for(int j=0;j<3;++j)u2+=std::pow(v[j+1]-un*f.normal[j],2);
  return wall_model(gas.rho(v),gas.mu,v[5],v[6],std::sqrt(u2),f.dn);
 }
 SmoothWallThermal thermal_wall(const Face&f)const{
  const auto&v=w[f.l];D3 velocity={v[1],v[2],v[3]};double un=dot3(velocity,f.normal);for(int j=0;j<3;++j)velocity[j]-=un*f.normal[j];double speed=length3(velocity),rho=gas.rho(v),star=wall_friction_velocity_from_k(rho,gas.mu,v[5],speed,f.dn);
  return smooth_wall_thermal(star*f.dn*rho/gas.mu,star,rho,gas.cp(v[4]),gas.mu*gas.cp(v[4])/gas.conductivity,mesh.meta[26],speed);
 }
 void update_wall_temperature(){
  if(!options.enhanced_thermal_wall)return;double epsilon=1e-8*gas.T0/(gas.rho0*std::pow(gas.a0,3));
  for(auto&f:mesh.faces)if(f.type==3)wall_temperature[wall_gradient_slot.at(f.id)]=smooth_adiabatic_wall_temperature(thermal_wall(f),w[f.l][4],f.dn,gas.conductivity,epsilon);
 }
 W7 virtual_neighbor(const Face&f)const{
  W7 b=w[f.l];if(f.type==3){for(int j=1;j<=3;++j)b[j]=-b[j];auto model=wall(f);b[6]=options.cell_omega_wall?b[6]+2*f.dn*model.omega_outward_derivative:std::max(1e-12,2*model.omega-b[6]);}
  else if(f.type==7){double un=b[1]*f.normal[0]+b[2]*f.normal[1]+b[3]*f.normal[2];for(int j=0;j<3;++j)b[j+1]-=2*un*f.normal[j];}
  else if(f.type==9){auto inf=freestream();double un=b[1]*f.normal[0]+b[2]*f.normal[1]+b[3]*f.normal[2];
   if(un<0)b=inf;else{double dp=inf[0]-b[0],r=gas.rho(b),a=gas.sound(b);b[0]=inf[0];for(int j=0;j<3;++j)b[j+1]-=dp/(r*a)*f.normal[j];}
  }return b;
 }
 D3 neighbor_dx(const Face&f)const{
  if(f.r>=0)return sub3(mesh.center[f.r],mesh.center[f.l]);
  if(f.type==7&&options.symmetry_normal_lsq)return {2*f.dn*f.normal[0],2*f.dn*f.normal[1],2*f.dn*f.normal[2]};
  if(options.boundary_face_ray){auto dx=sub3(f.center,mesh.center[f.l]);for(double&v:dx)v*=2;return dx;}
  return {2*f.dn*f.normal[0],2*f.dn*f.normal[1],2*f.dn*f.normal[2]};
 }
 void setup_lsq(){
  std::vector<std::array<double,9>>a(mesh.n);
  for(auto&f:mesh.faces){D3 dx=neighbor_dx(f);double weight=1/dot3(dx,dx);for(int j=0;j<3;++j)for(int k=0;k<3;++k){double t=weight*dx[j]*dx[k];a[f.l][3*j+k]+=t;if(f.r>=0&&f.r<mesh.n)a[f.r][3*j+k]+=t;}}
  for(int i=0;i<mesh.n;++i)inverse_lsq[i]=binverse<3>(a[i]);
 }
 void update_reference_velocity(){
  if(!options.coupled_ausm)return;
  for(int i=0;i<mesh.n;++i){double v2=w[i][1]*w[i][1]+w[i][2]*w[i][2]+w[i][3]*w[i][3],nu=(gas.mu+prop[i][0])/gas.rho(w[i]);reference_velocity[i][0]=std::max(v2,std::pow(4*nu/std::cbrt(mesh.volume[i]),2));}
  for(auto&f:mesh.faces)if(f.r>=0){double rl=gas.rho(w[f.l]),rr=gas.rho(w[f.r]);double jump=4*std::abs(w[f.l][0]-w[f.r][0])/(rl+rr),nu=std::max((gas.mu+prop[f.l][0])/rl,(gas.mu+prop[f.r][0])/rr);jump=std::max(jump,std::pow(4*nu/f.dn,2));reference_velocity[f.l][0]=std::max(reference_velocity[f.l][0],jump);if(f.r<mesh.n)reference_velocity[f.r][0]=std::max(reference_velocity[f.r][0],jump);}
  for(int i=0;i<mesh.n;++i){double a2=std::pow(gas.sound(w[i]),2),v=reference_velocity[i][0];reference_velocity[i][0]=std::min(a2,std::max(1e-20*a2,v/(.5-std::min(v/a2,.25))));}mesh.exchange(reference_velocity);
 }
 void reconstruct_normal_scalar(int slot,const std::vector<W7>*frozen=nullptr){auto generated=sst_transport->reconstruct(slot,w,frozen);for(int i=0;i<mesh.nt;++i){lo[i][5+slot]=generated.limiter[i][0];for(int j=0;j<3;++j){grad[i][15+3*slot+j]=generated.raw[i][j];scalar_gradient[i][3*slot+j]=generated.physical[i][j];}}}
 void gradients(const std::vector<W7>*frozen_limiter=nullptr,bool defer_normal_k=false){
  if(defer_normal_k&&(!options.normal_sst_transport||frozen_limiter))throw std::invalid_argument("deferred normal-k reconstruction requires unfrozen normal SST caller context");
  double start=now();mesh.exchange(w);std::fill(grad.begin(),grad.end(),Gradient{});lo=w;hi=w;
#ifdef M6_SKIP_UNUSED_NORMAL_LSQ
  // Normal SST replaces all scalar gradients/limiters before their first live
  // consumer. Keep every frozen-limiter path and every halo exchange unchanged.
  // hi[5,6] is dead local scratch, not a persisted or subsequently read output.
  const int early_variables=options.normal_sst_transport&&frozen_limiter==nullptr?5:7;
#else
  constexpr int early_variables=7;
#endif
  for(auto&f:mesh.faces){W7 b=f.r>=0?w[f.r]:virtual_neighbor(f);D3 dx=neighbor_dx(f);double weight=1/dot3(dx,dx);
   if(f.type==9&&options.coupled_farfield_reconstruction)for(int k=0;k<5;++k)b[k]=w[f.l][k];
   if(f.type==3&&options.coupled_wall_reconstruction){const auto&v=w[f.l];auto boundary=stationary_wall_virtual_state(gas.rho(v),gas.mu,v[5],D3{v[1],v[2],v[3]},f.normal,f.dn);for(int k=0;k<3;++k)b[k+1]=2*boundary.velocity[k]-v[k+1];}
   if(f.type==3&&options.enhanced_thermal_wall)b[4]=2*wall_temperature[wall_gradient_slot.at(f.id)]-w[f.l][4];
   for(int k=0;k<early_variables;++k){double delta=b[k]-w[f.l][k];if(!((f.type==3&&options.coupled_wall_reconstruction||f.type==9&&options.coupled_farfield_reconstruction)&&k<5)){lo[f.l][k]=std::min(lo[f.l][k],b[k]);hi[f.l][k]=std::max(hi[f.l][k],b[k]);}
    for(int j=0;j<3;++j)grad[f.l][3*k+j]+=weight*dx[j]*delta;
    if(f.r>=0&&f.r<mesh.n){lo[f.r][k]=std::min(lo[f.r][k],w[f.l][k]);hi[f.r][k]=std::max(hi[f.r][k],w[f.l][k]);for(int j=0;j<3;++j)grad[f.r][3*k+j]+=weight*dx[j]*delta;}
   }
  }
  for(int i=0;i<mesh.n;++i)for(int k=0;k<early_variables;++k){D3 v={grad[i][3*k],grad[i][3*k+1],grad[i][3*k+2]};for(int j=0;j<3;++j){grad[i][3*k+j]=0;for(int l=0;l<3;++l)grad[i][3*k+j]+=inverse_lsq[i][3*j+l]*v[l];}}
  if(options.coupled_wall_reconstruction){
   // One current-gradient wall extrapolation in face order. This is not a
   // converged Neumann solve and does not remove the wall from the matrix.
   for(auto&f:mesh.faces)if(f.type==3){auto ray=sub3(f.center,mesh.center[f.l]);double delta=0;for(int j=0;j<3;++j)delta+=grad[f.l][j]*ray[j];for(int j=0;j<3;++j)for(int l=0;l<3;++l)grad[f.l][j]+=inverse_lsq[f.l][3*j+l]*ray[l]/dot3(ray,ray)*delta;}
   // Wall extrema use normal extrapolation and tangential velocity, whereas
   // the weighted LSQ matrix above uses the face-centre connector.
   for(auto&f:mesh.faces)if(f.type==3){W7 boundary=w[f.l];for(int k=0;k<5;++k)for(int j=0;j<3;++j)boundary[k]+=2*f.dn*grad[f.l][3*k+j]*f.normal[j];double un=0;for(int j=0;j<3;++j)un+=boundary[j+1]*f.normal[j];for(int j=0;j<3;++j)boundary[j+1]-=un*f.normal[j];for(int k=0;k<5;++k){lo[f.l][k]=std::min(lo[f.l][k],boundary[k]);hi[f.l][k]=std::max(hi[f.l][k],boundary[k]);}}
  }
  if(options.coupled_farfield_reconstruction)for(auto&f:mesh.faces)if(f.type==9){for(int k=0;k<5;++k){double boundary=w[f.l][k];for(int j=0;j<3;++j)boundary+=2*f.dn*grad[f.l][3*k+j]*f.normal[j];lo[f.l][k]=std::min(lo[f.l][k],boundary);hi[f.l][k]=std::max(hi[f.l][k],boundary);}}
  if(options.turbulence_gauss&&!options.normal_sst_transport){for(int i=0;i<mesh.n;++i)for(int j=15;j<21;++j)grad[i][j]=0.;for(auto&f:mesh.faces){auto right=f.r>=0?w[f.r]:virtual_neighbor(f);for(int k=5;k<7;++k){double value=f.r>=0?2*w[f.l][k]*right[k]/(w[f.l][k]+right[k]):.5*(w[f.l][k]+right[k]);for(int j=0;j<3;++j){double area=f.normal[j]*f.area;grad[f.l][3*k+j]+=(value-w[f.l][k])*area;if(f.r>=0&&f.r<mesh.n)grad[f.r][3*k+j]-=(value-w[f.r][k])*area;}}}for(int i=0;i<mesh.n;++i)for(int j=15;j<21;++j)grad[i][j]/=mesh.volume[i];}
  // Limit the actual face increments; a smooth rational limiter is continuous
  // near extrema and never changes an exactly linear admissible increment.
  std::vector<W7>phi(mesh.n);for(auto&v:phi)v.fill(1);
  for(auto&f:mesh.faces)for(int side=0;side<2;++side){int i=side==0?f.l:f.r;if(i<0||i>=mesh.n)continue;D3 dx=sub3(f.center,mesh.center[i]);
   for(int k=0;k<early_variables;++k){double delta=0;for(int j=0;j<3;++j)delta+=grad[i][3*k+j]*dx[j];double d=std::abs(delta),bound=delta>0?hi[i][k]-w[i][k]:w[i][k]-lo[i][k];
    if(bound>=2*d)continue;double t=bound*bound+4*d*d;double fraction=bound*t/(bound*bound*bound+d*t+1e-100);if(options.coupled_wall_reconstruction&&k<5){double scale=k==0?gas.rho0*gas.a0*gas.a0:k==4?gas.T0:gas.a0,epsilon=1e-32/(scale*scale*scale);fraction=(bound*t+epsilon)/(bound*bound*bound+d*t+epsilon);}phi[i][k]=std::min(phi[i][k],fraction);
   }
  }
  // Store reconstruction gradients separately from raw viscous/SST gradients.
  lo.resize(mesh.nt);for(int i=0;i<mesh.n;++i)lo[i]=frozen_limiter?(*frozen_limiter)[i]:phi[i];mesh.exchange(lo);mesh.exchange(grad);
  if(options.reconstruct_enthalpy){std::fill(enthalpy_gradient.begin(),enthalpy_gradient.end(),std::array<double,4>{});for(int i=0;i<mesh.nt;++i)total_enthalpy[i]=gas.enthalpy(w[i]);std::vector<double>hlo=total_enthalpy,hhi=total_enthalpy;
   for(auto&f:mesh.faces){double right=f.r>=0?total_enthalpy[f.r]:gas.enthalpy(virtual_neighbor(f));D3 dx=neighbor_dx(f);double weight=1/dot3(dx,dx),delta=right-total_enthalpy[f.l];hlo[f.l]=std::min(hlo[f.l],right);hhi[f.l]=std::max(hhi[f.l],right);
    for(int j=0;j<3;++j)enthalpy_gradient[f.l][j]+=weight*dx[j]*delta;if(f.r>=0&&f.r<mesh.n){hlo[f.r]=std::min(hlo[f.r],total_enthalpy[f.l]);hhi[f.r]=std::max(hhi[f.r],total_enthalpy[f.l]);for(int j=0;j<3;++j)enthalpy_gradient[f.r][j]+=weight*dx[j]*delta;}
   }
   for(int i=0;i<mesh.n;++i){D3 b={enthalpy_gradient[i][0],enthalpy_gradient[i][1],enthalpy_gradient[i][2]};for(int j=0;j<3;++j){enthalpy_gradient[i][j]=0;for(int k=0;k<3;++k)enthalpy_gradient[i][j]+=inverse_lsq[i][3*j+k]*b[k];}enthalpy_gradient[i][3]=1;}
   for(auto&f:mesh.faces)for(int side=0;side<2;++side){int i=side==0?f.l:f.r;if(i<0||i>=mesh.n)continue;D3 dx=sub3(f.center,mesh.center[i]);double delta=0;for(int j=0;j<3;++j)delta+=enthalpy_gradient[i][j]*dx[j];double d=std::abs(delta),bound=delta>0?hhi[i]-total_enthalpy[i]:total_enthalpy[i]-hlo[i];if(bound>=2*d)continue;double t=bound*bound+4*d*d;enthalpy_gradient[i][3]=std::min(enthalpy_gradient[i][3],bound*t/(bound*bound*bound+d*t+1e-100));}
   mesh.exchange(enthalpy_gradient);
  }
  if(options.independent_face_thermo)thermodynamic_gradients();
  if(options.viscous_gauss){std::fill(physical_gradient.begin(),physical_gradient.end(),std::array<double,12>{});
   for(auto&f:mesh.faces){std::array<double,4>value{};if(f.r>=0){auto l=reconstruct(f.l,f),r=reconstruct(f.r,f);for(int k=0;k<4;++k)value[k]=.5*(l[k+1]+r[k+1]);}
    else{auto b=virtual_neighbor(f);for(int k=0;k<4;++k)value[k]=.5*(w[f.l][k+1]+b[k+1]);}
    if(f.type==9&&options.coupled_farfield_reconstruction){auto boundary=reconstruct(f.l,f);for(int k=0;k<4;++k)value[k]=boundary[k+1];}
    if(f.type==7&&options.symmetry_face_gradients){auto boundary=reconstruct(f.l,f);double un=0;for(int j=0;j<3;++j)un+=boundary[j+1]*f.normal[j];for(int j=0;j<3;++j)value[j]=boundary[j+1]-un*f.normal[j];}
    if(f.type==3&&options.wall_slip_gradients){const auto&v=w[f.l];D3 velocity={v[1],v[2],v[3]};auto boundary=stationary_wall_virtual_state(gas.rho(v),gas.mu,v[5],velocity,f.normal,f.dn);for(int k=0;k<3;++k)value[k]=boundary.velocity[k];wall_virtual_jacobian[wall_gradient_slot.at(f.id)]=boundary.jacobian;}
    if(f.type==3&&options.enhanced_thermal_wall){const auto&v=w[f.l];int slot=wall_gradient_slot.at(f.id);auto thermal=thermal_wall(f);double jump=wall_temperature[slot]-v[4],rho=gas.rho(v);D3 tangent={v[1],v[2],v[3]};double un=dot3(tangent,f.normal);for(int j=0;j<3;++j)tangent[j]-=un*f.normal[j];double speed=length3(tangent),star=wall_friction_velocity_from_k(rho,gas.mu,v[5],speed,f.dn);auto&jacobian=wall_temperature_jacobian[slot];for(int j=0;j<3;++j)jacobian[j]=speed>1e-30?jump*thermal.fraction_y_derivative*tangent[j]/(2*star*speed):0.;jacobian[3]=jump*thermal.fraction_y_derivative*(f.dn*star/gas.mu-speed/(2*rho*star));jacobian[4]=1-thermal.virtual_fraction+jump*thermal.fraction_pr_derivative*gas.mu/gas.conductivity*gas.cp_derivative(v[4]);value[3]=v[4]+jump*thermal.virtual_fraction;}
    for(int k=0;k<4;++k)for(int j=0;j<3;++j){double a=f.normal[j]*f.area;physical_gradient[f.l][3*k+j]+=(value[k]-w[f.l][k+1])*a;if(f.r>=0&&f.r<mesh.n)physical_gradient[f.r][3*k+j]-=(value[k]-w[f.r][k+1])*a;}
   }
   for(int i=0;i<mesh.n;++i)for(double&v:physical_gradient[i])v/=mesh.volume[i];mesh.exchange(physical_gradient);
  }
  if(options.normal_sst_transport){if(!defer_normal_k)reconstruct_normal_scalar(0,frozen_limiter);reconstruct_normal_scalar(1,frozen_limiter);}else if(options.turbulence_gauss){std::fill(scalar_gradient.begin(),scalar_gradient.end(),std::array<double,6>{});for(auto&f:mesh.faces){auto b=f.r>=0?w[f.r]:virtual_neighbor(f);for(int k=5;k<7;++k){double value=f.r>=0?.5*(reconstruct_scalar(f.l,f,k)+reconstruct_scalar(f.r,f,k)):.5*(w[f.l][k]+b[k]);for(int j=0;j<3;++j){double area=f.area*f.normal[j];scalar_gradient[f.l][3*(k-5)+j]+=(value-w[f.l][k])*area;if(f.r>=0&&f.r<mesh.n)scalar_gradient[f.r][3*(k-5)+j]-=(value-w[f.r][k])*area;}}}for(int i=0;i<mesh.n;++i)for(double&v:scalar_gradient[i])v/=mesh.volume[i];mesh.exchange(scalar_gradient);}
  gradient_time+=now()-start;
 }
 W7 reconstruct(int cell,const Face&f)const{
  W7 b=w[cell];if(options.order==1)return b;D3 dx=sub3(f.center,mesh.center[cell]);
  for(int k=0;k<7;++k){double delta=0;for(int j=0;j<3;++j)delta+=grad[cell][3*k+j]*dx[j];b[k]+=lo[cell][k]*delta;}
  if(options.reconstruct_enthalpy){double h=total_enthalpy[cell];for(int j=0;j<3;++j)h+=enthalpy_gradient[cell][3]*enthalpy_gradient[cell][j]*dx[j];h-=.5*(b[1]*b[1]+b[2]*b[2]+b[3]*b[3]);b[4]=gas.temperature_from_h(h,b[4]);}
  if(!(b[0]>0&&b[4]>0)){std::ostringstream s;s<<"nonpositive mean-flow face state face="<<f.id<<" cell="<<mesh.global[cell]<<" p="<<b[0]<<" T="<<b[4];throw std::runtime_error(s.str());}return b;
 }
 double reconstruct_scalar(int cell,const Face&f,int variable)const{double value=w[cell][variable];if(options.order==2){D3 dx=sub3(f.center,mesh.center[cell]);double delta=0;for(int j=0;j<3;++j)delta+=grad[cell][3*variable+j]*dx[j];value+=lo[cell][variable]*delta;}if(!(value>0)){std::ostringstream s;s<<"nonpositive scalar face state face="<<f.id<<" cell="<<mesh.global[cell]<<" variable="<<variable<<" value="<<value;throw std::runtime_error(s.str());}return value;}
 // Density and total enthalpy are reconstructed independently for the
 // inviscid flux. Temperature retains its own reconstruction for diffusion.
 void thermodynamic_gradients(){
  std::fill(thermo_gradient.begin(),thermo_gradient.end(),std::array<double,6>{});
  for(int i=0;i<mesh.nt;++i)thermodynamics[i]={gas.rho(w[i]),gas.enthalpy(w[i])};
  auto low=thermodynamics,high=thermodynamics;
  for(auto&f:mesh.faces){std::array<double,2>right;if(f.r>=0)right=thermodynamics[f.r];else{auto b=virtual_neighbor(f);right={gas.rho(b),gas.enthalpy(b)};if(f.type==9&&options.coupled_farfield_reconstruction)right=thermodynamics[f.l];}auto dx=neighbor_dx(f);double weight=1/dot3(dx,dx);
   for(int k=0;k<2;++k){double delta=right[k]-thermodynamics[f.l][k];low[f.l][k]=std::min(low[f.l][k],right[k]);high[f.l][k]=std::max(high[f.l][k],right[k]);for(int j=0;j<3;++j)thermo_gradient[f.l][3*k+j]+=weight*dx[j]*delta;
    if(f.r>=0&&f.r<mesh.n){low[f.r][k]=std::min(low[f.r][k],thermodynamics[f.l][k]);high[f.r][k]=std::max(high[f.r][k],thermodynamics[f.l][k]);for(int j=0;j<3;++j)thermo_gradient[f.r][3*k+j]+=weight*dx[j]*delta;}
   }
  }
  for(int i=0;i<mesh.n;++i){thermo_limiter[i]={1,1};for(int k=0;k<2;++k){D3 v={thermo_gradient[i][3*k],thermo_gradient[i][3*k+1],thermo_gradient[i][3*k+2]};for(int j=0;j<3;++j){thermo_gradient[i][3*k+j]=0;for(int l=0;l<3;++l)thermo_gradient[i][3*k+j]+=inverse_lsq[i][3*j+l]*v[l];}}}
  // No face storage exists for these auxiliary scalars. The scalar
  // limiter extends extrema by a mirror-point normal extrapolation.
  for(auto&f:mesh.faces)if(f.r<0&&f.type!=7){for(int k=0;k<2;++k){double delta=0;for(int j=0;j<3;++j)delta+=2*f.dn*thermo_gradient[f.l][3*k+j]*f.normal[j];double ghost=thermodynamics[f.l][k]+delta;low[f.l][k]=std::min(low[f.l][k],ghost);high[f.l][k]=std::max(high[f.l][k],ghost);}}
  for(auto&f:mesh.faces)for(int side=0;side<2;++side){int i=side==0?f.l:f.r;if(i<0||i>=mesh.n)continue;auto dx=sub3(f.center,mesh.center[i]);for(int k=0;k<2;++k){double delta=0;for(int j=0;j<3;++j)delta+=thermo_gradient[i][3*k+j]*dx[j];double d=std::abs(delta),bound=delta>0?high[i][k]-thermodynamics[i][k]:thermodynamics[i][k]-low[i][k];if(bound>=2*d)continue;double t=bound*bound+4*d*d,scale=k==0?gas.rho0:gas.a0*gas.a0,epsilon=1e-32/(scale*scale*scale);thermo_limiter[i][k]=std::min(thermo_limiter[i][k],(bound*t+epsilon)/(bound*bound*bound+d*t+epsilon));}}
  mesh.exchange(thermo_gradient);mesh.exchange(thermo_limiter);
 }
 std::array<double,2>reconstruct_thermodynamics(int cell,const Face&f)const{
  auto v=thermodynamics[cell];if(options.order==2){auto dx=sub3(f.center,mesh.center[cell]);for(int k=0;k<2;++k)for(int j=0;j<3;++j)v[k]+=thermo_limiter[cell][k]*thermo_gradient[cell][3*k+j]*dx[j];}
  if(!(v[0]>0&&std::isfinite(v[1])))throw std::runtime_error("invalid independent density/enthalpy face reconstruction");return v;
 }
 std::array<double,4>face_thermodynamics(const Face&f)const{auto l=reconstruct_thermodynamics(f.l,f),r=reconstruct_thermodynamics(f.r,f);return {l[0],r[0],l[1],r[1]};}
 double boundary_face_density(const Face&f,const W7&face)const{return options.independent_face_thermo?reconstruct_thermodynamics(f.l,f)[0]:gas.rho(face);}
 void properties(bool keep_transport=false,bool keep_production=false){
  std::fill(wall_weight.begin(),wall_weight.end(),0);std::fill(wall_pk.begin(),wall_pk.end(),0);
  for(auto&f:mesh.faces)if(f.type==3){auto x=wall(f);wall_weight[f.l]+=f.area;wall_pk[f.l]+=f.area*x.production;}
  const double beta1=mesh.meta[27],beta2=mesh.meta[28],sw1=1/mesh.meta[31],sw2=1/mesh.meta[32],betaStar=mesh.meta[34],kap=mesh.meta[35];
   for(int i=0;i<mesh.n;++i){const auto&v=w[i];auto g=grad[i];if(options.viscous_gauss)for(int j=0;j<9;++j)g[3+j]=physical_gradient[i][j];if(options.turbulence_gauss)for(int j=0;j<6;++j)g[15+j]=scalar_gradient[i][j];double r=gas.rho(v),k=std::max(v[5],options.normal_scalar_step?1e-30:1e-18),omega=std::max(v[6],options.normal_scalar_step?1e-30:1e-15),d=mesh.distance[i],nu=gas.mu/r;
   double dotkw=0;for(int j=0;j<3;++j)dotkw+=g[15+j]*g[18+j];double cross=2*r*sw2/omega*dotkw,cd=std::max(cross,1e-30);
   double arg1=std::min(std::max(std::sqrt(k)/(betaStar*omega*d),500*nu/(d*d*omega)),4*r*sw2*k/(cd*d*d));
   double f1=std::tanh(std::pow(std::min(arg1,100.),4)),arg2=std::max(2*std::sqrt(k)/(betaStar*omega*d),500*nu/(d*d*omega)),f2=std::tanh(std::pow(std::min(arg2,100.),2));
   double strain2=0;for(int j=0;j<3;++j)for(int l=0;l<3;++l)strain2+=.5*std::pow(g[3*(j+1)+l]+g[3*(l+1)+j],2);
   D3 vort={g[3*3+1]-g[3*2+2],g[3*1+2]-g[3*3+0],g[3*2+0]-g[3*1+1]};double om=length3(vort);
   double eddy_measure=options.sst_vorticity?om:std::sqrt(strain2);
   double mut=r*mesh.meta[33]*k/std::max(mesh.meta[33]*omega,eddy_measure*f2);mut=std::max(0.,std::min(mut,1e5*gas.mu));if(options.normal_scalar_step&&v[6]*gas.a0<1e-20)mut=1e-10/(gas.rho0*gas.a0);if(keep_transport)mut=prop[i][0];
   double pk=std::max(0.,std::min(mut*strain2,mesh.meta[36]*betaStar*r*k*omega));
   if(options.wall_production&&wall_weight[i]>0)pk=wall_pk[i]/wall_weight[i];if(keep_production)pk=prop[i][4];
   double beta=f1*beta1+(1-f1)*beta2,gamma=f1*(beta1/betaStar-sw1*kap*kap/std::sqrt(betaStar))+(1-f1)*(beta2/betaStar-sw2*kap*kap/std::sqrt(betaStar));
   prop[i]={mut,f1,beta,gamma,pk,(1-f1)*cross,strain2,om};
  }mesh.exchange(prop);if(!keep_transport)update_reference_velocity();
 }
 D3 face_gradient(const Face&f,int k)const{
  D3 g{};int l=f.l,r=f.r;for(int j=0;j<3;++j){auto derivative=[&](int i){if(options.turbulence_gauss&&k>=5)return scalar_gradient[i][3*(k-5)+j];return options.viscous_gauss&&k>=1&&k<=4?physical_gradient[i][3*(k-1)+j]:grad[i][3*k+j];};g[j]=r>=0?.5*(derivative(l)+derivative(r)):derivative(l);}
  if(r>=0){D3 dx=sub3(mesh.center[r],mesh.center[l]);g=correct_face_gradient(g,dx,f.normal,w[r][k]-w[l][k],options.edge_correction,options.nonorth_limit);}
  else if(f.type==7||f.type==3){double target=0;if(f.type==3&&k>=1&&k<=3)target=-w[l][k]/f.dn;double correction=target-dot3(g,f.normal);for(int j=0;j<3;++j)g[j]+=correction*f.normal[j];}
  return g;
 }
 std::array<double,2>face_transport(const Face&f)const{
  double ml=gas.mu+(f.type==3?0:prop[f.l][0]),kl=gas.conductivity+gas.cp(w[f.l][4])*(f.type==3?0:prop[f.l][0])/mesh.meta[26];
  if(f.r<0)return {ml,kl};double mr=gas.mu+prop[f.r][0],kr=gas.conductivity+gas.cp(w[f.r][4])*prop[f.r][0]/mesh.meta[26];
  if(options.harmonic_diffusion)return {2*ml*mr/(ml+mr),2*kl*kr/(kl+kr)};return {.5*(ml+mr),.5*(kl+kr)};
 }
 ViscousFlux symmetry_viscous(const Face&f,double mu)const{
  std::array<D3,3>g{};D3 velocity{};for(int k=0;k<3;++k){velocity[k]=w[f.l][k+1];for(int j=0;j<3;++j)g[k][j]=physical_gradient[f.l][3*k+j];}
  return m6::symmetry_viscous_flux(f.normal,f.dn,g,mu,velocity);
 }
 Q5 flux(const Face&f){
  auto left=reconstruct(f.l,f);Q5 result{};double pwall=0;
  if(f.type==9&&options.characteristic_farfield){result=m6::characteristic_farfield(gas,w[f.l],left,freestream(),f.normal,gas.rho(w[f.l]),std::pow(gas.sound(w[f.l]),2),boundary_face_density(f,left),gas.enthalpy(w[f.l])).flux;}
  else if(f.type==3||f.type==7){auto right=left;double un=left[1]*f.normal[0]+left[2]*f.normal[1]+left[3]*f.normal[2];for(int j=0;j<3;++j)right[j+1]-=2*un*f.normal[j];if(options.coupled_ausm){double pressure=coupled_wall_pressure(left,f.normal,reference_velocity[f.l][0],boundary_face_density(f,left));for(int j=0;j<3;++j)result[j+1]=pressure*f.normal[j];}else result=gas.ausm(left,right,f.normal);
   result[0]=0;result[4]=0;for(int j=0;j<3;++j)pwall+=result[j+1]*f.normal[j];
  }else{auto right=f.r>=0?reconstruct(f.r,f):freestream();if(options.coupled_ausm&&f.r>=0){auto thermo=options.independent_face_thermo?face_thermodynamics(f):std::array<double,4>{};result=gas.coupled_ausm(left,right,f.normal,std::pow(gas.sound(w[f.l]),2),std::pow(gas.sound(w[f.r]),2),reference_velocity[f.l][0],reference_velocity[f.r][0],options.independent_face_thermo?&thermo:nullptr);}else result=gas.ausm(left,right,f.normal);}
  D3 traction{};double energy_viscous=0;
  if(f.type==3){auto wm=wall(f);double un=w[f.l][1]*f.normal[0]+w[f.l][2]*f.normal[1]+w[f.l][3]*f.normal[2];
   for(int j=0;j<3;++j)traction[j]=-wm.velocity_coefficient*(w[f.l][j+1]-un*f.normal[j]);
  }else if(f.type==7&&options.symmetry_viscous_flux){auto viscous=symmetry_viscous(f,face_transport(f)[0]);traction=viscous.traction;energy_viscous=viscous.energy;
  }else if(f.type!=7&&!(f.type==9&&options.farfield_zero_viscous_flux)){std::array<D3,4>g={face_gradient(f,1),face_gradient(f,2),face_gradient(f,3),face_gradient(f,4)};auto transport=face_transport(f);D3 velocity{};for(int j=0;j<3;++j)velocity[j]=f.r>=0?.5*(w[f.l][j+1]+w[f.r][j+1]):w[f.l][j+1];
   auto viscous=m6::viscous_flux(f.normal,g,transport[0],transport[1],velocity);traction=viscous.traction;energy_viscous=viscous.energy;
  }
  for(int j=0;j<3;++j)result[j+1]-=traction[j];result[4]-=energy_viscous;
  if(options.tke_pressure){double kpress=2./3*gas.rho(left)*left[5];for(int j=0;j<3;++j)result[j+1]+=kpress*f.normal[j];}
  if(f.type==3){double scale=gas.rho0*gas.a0*gas.a0*f.area;for(int j=0;j<3;++j){pressure_force[j]+=(pwall-gas.pinf/(gas.rho0*gas.a0*gas.a0))*f.normal[j]*scale;viscous_force[j]-=traction[j]*scale;}}
  return result;
 }
 D3 diffusion_direction(const Face&f)const{D3 s=f.normal;for(double&v:s)v/=f.dn;if(options.edge_correction&&options.edge_implicit&&f.r>=0){auto dx=sub3(mesh.center[f.r],mesh.center[f.l]);double distance2=dot3(dx,dx);if(f.dn/std::sqrt(distance2)<.2)for(int j=0;j<3;++j)s[j]=dx[j]/distance2;}return s;}
 B5 viscous_jacobian(const Face&f,const W7&cell,double mu,double conduct,const D3&u)const{
  auto d=gas.primitive_jacobian(cell);B5 k{};auto s=diffusion_direction(f);double sn=dot3(s,f.normal);
  for(int col=0;col<5;++col){double un=0,us=0;for(int j=0;j<3;++j){un+=f.normal[j]*d[5*(j+1)+col];us+=s[j]*d[5*(j+1)+col];}
   for(int j=0;j<3;++j){double derivative=mu*f.area*(sn*d[5*(j+1)+col]+s[j]*un-2./3*f.normal[j]*us);k[5*(j+1)+col]=derivative;k[20+col]+=u[j]*derivative;}
   k[20+col]+=conduct*f.area*sn*d[20+col];
  }return k;
 }
#ifdef M6_CELL_JACOBIAN_CACHE
 B5 viscous_jacobian_cached(const Face&f,const B5&d,double mu,double conduct,const D3&u)const{
  B5 k{};auto s=diffusion_direction(f);double sn=dot3(s,f.normal);
  for(int col=0;col<5;++col){double un=0,us=0;for(int j=0;j<3;++j){un+=f.normal[j]*d[5*(j+1)+col];us+=s[j]*d[5*(j+1)+col];}
   for(int j=0;j<3;++j){double derivative=mu*f.area*(sn*d[5*(j+1)+col]+s[j]*un-2./3*f.normal[j]*us);k[5*(j+1)+col]=derivative;k[20+col]+=u[j]*derivative;}
   k[20+col]+=conduct*f.area*sn*d[20+col];
  }return k;
 }
#endif
 B5 symmetry_viscous_jacobian(const Face&f,double mu)const{
  B5 result{};auto d=gas.primitive_jacobian(w[f.l]);auto base=symmetry_viscous(f,mu);D3 velocity={w[f.l][1],w[f.l][2],w[f.l][3]};double un=dot3(velocity,f.normal);for(int j=0;j<3;++j)velocity[j]-=un*f.normal[j];
  for(int col=0;col<5;++col){D3 delta{};for(int j=0;j<3;++j)delta[j]=d[5*(j+1)+col];double dn=dot3(delta,f.normal);D3 traction{};for(int j=0;j<3;++j){traction[j]=-4./3*mu*dn/f.dn*f.normal[j];result[5*(j+1)+col]=-traction[j]*f.area;delta[j]-=dn*f.normal[j];}result[20+col]=-(dot3(traction,velocity)+dot3(base.traction,delta))*f.area;}
  return result;
 }
 template<size_t N>std::array<double,N>linear_boundary(const Face&f,const std::vector<std::array<double,N>>&values,int variable)const{
  auto b=values[f.l];
  if constexpr(N==4){if(f.type==3)for(int j=0;j<3;++j)b[j]=-b[j];else if(f.type==7){double un=0;for(int j=0;j<3;++j)un+=b[j]*f.normal[j];for(int j=0;j<3;++j)b[j]-=2*un*f.normal[j];}else if(f.type==9){double un=0;for(int j=0;j<3;++j)un+=w[f.l][j+1]*f.normal[j];if(un<0)b.fill(0);}}
  else if constexpr(N==1){if(f.type==3&&variable==6&&!options.cell_omega_wall)b[0]=2*wall(f).omega-w[f.l][6]>1e-12?-b[0]:0.;else if(f.type==9){double un=0;for(int j=0;j<3;++j)un+=w[f.l][j+1]*f.normal[j];if(un<0)b[0]=0;}}
  return b;
 }
 template<size_t N>std::vector<std::array<double,3*N>>linear_gradients(std::vector<std::array<double,N>>&values,int variable=-1,const Field<5>*conserved_delta=nullptr)const{
  mesh.exchange(values);std::vector<std::array<double,3*N>>g(mesh.nt);
  if constexpr(N==1){if(options.turbulence_gauss){for(size_t fi=0;fi<mesh.faces.size();++fi){auto&f=mesh.faces[fi];auto b=f.r>=0?values[f.r]:linear_boundary(f,values,variable);double face=.5*(values[f.l][0]+b[0]);if(f.r>=0){double a=options.normal_sst_transport?sst_transport->reference[variable-5][f.l]:w[f.l][variable],right=options.normal_sst_transport?sst_transport->reference[variable-5][f.r]:w[f.r][variable],den=(a+right)*(a+right);face=2*(right*right*values[f.l][0]+a*a*b[0])/den;}else if(options.normal_sst_transport)face=sst_transport->allocated[variable-5][fi]?0.:values[f.l][0];for(int j=0;j<3;++j){double area=f.area*f.normal[j];g[f.l][j]+=(face-values[f.l][0])*area;if(f.r>=0&&f.r<mesh.n)g[f.r][j]-=(face-values[f.r][0])*area;}}for(int i=0;i<mesh.n;++i)for(double&v:g[i])v/=mesh.volume[i];mesh.exchange(g);return g;}}
  for(auto&f:mesh.faces){auto b=f.r>=0?values[f.r]:linear_boundary(f,values,variable);D3 dx=neighbor_dx(f);double weight=1/dot3(dx,dx);
   if constexpr(N==4){if(f.type==9&&options.coupled_farfield_reconstruction)b=values[f.l];}
   if constexpr(N==4){if(f.type==3&&options.coupled_wall_reconstruction){if(!conserved_delta)throw std::runtime_error("coupled wall LSQ derivative requires density direction");const auto&jacobian=wall_virtual_jacobian[wall_gradient_slot.at(f.id)];for(int k=0;k<3;++k){double virtual_delta=jacobian[4*k+3]*(*conserved_delta)[f.l][0];for(int j=0;j<3;++j)virtual_delta+=jacobian[4*k+j]*values[f.l][j];b[k]=2*virtual_delta-values[f.l][k];}}}
   if constexpr(N==4){if(f.type==3&&options.enhanced_thermal_wall)b[3]=-values[f.l][3];}
   for(size_t k=0;k<N;++k){double delta=b[k]-values[f.l][k];for(int j=0;j<3;++j){g[f.l][3*k+j]+=weight*dx[j]*delta;if(f.r>=0&&f.r<mesh.n)g[f.r][3*k+j]+=weight*dx[j]*delta;}}
  }
  for(int i=0;i<mesh.n;++i)for(size_t k=0;k<N;++k){D3 v={g[i][3*k],g[i][3*k+1],g[i][3*k+2]};for(int j=0;j<3;++j){g[i][3*k+j]=0;for(int l=0;l<3;++l)g[i][3*k+j]+=inverse_lsq[i][3*j+l]*v[l];}}
  mesh.exchange(g);return g;
 }
 template<size_t N>std::vector<std::array<double,3*N>>linear_physical_gradients(std::vector<std::array<double,N>>&values,int variable=-1,const Field<5>*conserved_delta=nullptr)const{
  if constexpr(N==4){if(options.wall_slip_gradients&&!conserved_delta)throw std::runtime_error("wall virtual-velocity derivative needs conservative density direction");}
  auto raw=linear_gradients(values,variable,conserved_delta);bool gauss=N==4?options.viscous_gauss:options.turbulence_gauss;if(!gauss)return raw;
  std::vector<std::array<double,3*N>>g(mesh.nt);int first=N==4?1:variable;
  for(auto&f:mesh.faces){auto boundary=f.r>=0?values[f.r]:linear_boundary(f,values,variable);auto dl=sub3(f.center,mesh.center[f.l]);D3 dr{},symmetry_face{};if(f.r>=0)dr=sub3(f.center,mesh.center[f.r]);if constexpr(N==4){if(f.type==7&&options.symmetry_face_gradients){for(int k=0;k<3;++k){symmetry_face[k]=values[f.l][k];if(options.order==2)for(int j=0;j<3;++j)symmetry_face[k]+=lo[f.l][k+1]*raw[f.l][3*k+j]*dl[j];}double un=dot3(symmetry_face,f.normal);for(int k=0;k<3;++k)symmetry_face[k]-=un*f.normal[k];}}for(size_t k=0;k<N;++k){double left=values[f.l][k],right=boundary[k];if(f.r>=0&&options.order==2){for(int j=0;j<3;++j){left+=lo[f.l][first+k]*raw[f.l][3*k+j]*dl[j];right+=lo[f.r][first+k]*raw[f.r][3*k+j]*dr[j];}}double face=.5*(left+right);
   if constexpr(N==1){if(options.normal_sst_transport&&f.r<0){bool stored=sst_transport->allocated[variable-5][static_cast<size_t>(&f-mesh.faces.data())]!=0;face=stored?0.:values[f.l][0];if(!stored&&options.order==2)for(int j=0;j<3;++j)face+=lo[f.l][variable]*raw[f.l][j]*dl[j];}}
   if constexpr(N==4){if(f.type==9&&options.coupled_farfield_reconstruction){face=values[f.l][k];if(options.order==2)for(int j=0;j<3;++j)face+=lo[f.l][k+1]*raw[f.l][3*k+j]*dl[j];}}
   if constexpr(N==4){if(f.type==7&&options.symmetry_face_gradients&&k<3)face=symmetry_face[k];}
   if constexpr(N==4){if(f.type==3&&options.wall_slip_gradients&&k<3){const auto&jacobian=wall_virtual_jacobian[wall_gradient_slot.at(f.id)];face=jacobian[4*k+3]*(*conserved_delta)[f.l][0];for(int j=0;j<3;++j)face+=jacobian[4*k+j]*values[f.l][j];}}
   if constexpr(N==4){if(f.type==3&&options.enhanced_thermal_wall&&k==3){const auto&jacobian=wall_temperature_jacobian[wall_gradient_slot.at(f.id)];face=jacobian[3]*(*conserved_delta)[f.l][0]+jacobian[4]*values[f.l][3];for(int j=0;j<3;++j)face+=jacobian[j]*values[f.l][j];}}
   for(int j=0;j<3;++j){double area=f.area*f.normal[j];g[f.l][3*k+j]+=(face-values[f.l][k])*area;if(f.r>=0&&f.r<mesh.n)g[f.r][3*k+j]-=(face-values[f.r][k])*area;}}}
  for(int i=0;i<mesh.n;++i)for(double&v:g[i])v/=mesh.volume[i];mesh.exchange(g);return g;
 }
 std::array<double,3>flow_diffusion_cell_coefficients(int i)const{
  double rho=gas.rho(w[i]);auto conserved=gas.conserved(w[i]);return {rho,conserved[4]/rho,gas.cp(w[i][4])-gas.rg};
 }
 std::array<double,8>flow_diffusion_face_coefficients(const Face&f)const{
  std::array<double,8>result{};if(f.r<0&&!(f.type==7&&options.symmetry_viscous_flux))return result;
  auto transport=face_transport(f);result[0]=transport[0];result[1]=transport[1];
  if(f.r<0){D3 velocity={w[f.l][1],w[f.l][2],w[f.l][3]};double un=dot3(velocity,f.normal);for(int j=0;j<3;++j)result[2+j]=velocity[j]-un*f.normal[j];return result;}
  for(int j=0;j<3;++j)result[2+j]=.5*(w[f.l][j+1]+w[f.r][j+1]);
  std::array<D3,3>original={face_gradient(f,1),face_gradient(f,2),face_gradient(f,3)};double div=original[0][0]+original[1][1]+original[2][2];D3 tau{};
  for(int j=0;j<3;++j)for(int k=0;k<3;++k)tau[j]+=transport[0]*(original[j][k]+original[k][j]-(j==k?2*div/3:0))*f.normal[k];
  for(int j=0;j<3;++j)result[5+j]=tau[j];return result;
 }
 void prepare_flow_diffusion_cache(){
  flow_diffusion_cache_active=false;if(!options.cache_flow_diffusion)return;double started=now();
  flow_diffusion_cell_cache.resize(mesh.n);flow_diffusion_face_cache.resize(mesh.faces.size());
  for(int i=0;i<mesh.n;++i)flow_diffusion_cell_cache[i]=flow_diffusion_cell_coefficients(i);
  for(size_t i=0;i<mesh.faces.size();++i)flow_diffusion_face_cache[i]=flow_diffusion_face_coefficients(mesh.faces[i]);
  flow_diffusion_cache_active=true;++flow_diffusion_cache_builds;flow_diffusion_cache_seconds+=now()-started;
 }
 void extra_flow_diffusion(const Field<5>&x,Field<5>&y)const{
  std::vector<std::array<double,4>>dw(mesh.nt);
  for(int i=0;i<mesh.n;++i){auto coefficients=flow_diffusion_cache_active?flow_diffusion_cell_cache[i]:flow_diffusion_cell_coefficients(i);double r=coefficients[0];for(int j=0;j<3;++j)dw[i][j]=(x[i][j+1]-w[i][j+1]*x[i][0])/r;
   double de=(x[i][4]-coefficients[1]*x[i][0])/r;for(int j=0;j<3;++j)de-=w[i][j+1]*dw[i][j];dw[i][3]=de/coefficients[2];
  }
  auto gd=linear_physical_gradients(dw,-1,&x);
  if(options.symmetry_viscous_flux)for(auto&f:mesh.faces)if(f.type==7){auto coefficients=flow_diffusion_cache_active?flow_diffusion_face_cache[static_cast<size_t>(&f-mesh.faces.data())]:flow_diffusion_face_coefficients(f);std::array<D3,3>g{};for(int k=0;k<3;++k)for(int j=0;j<3;++j)g[k][j]=gd[f.l][3*k+j];auto derivative=m6::symmetry_viscous_flux(f.normal,f.dn,g,coefficients[0],D3{});D3 velocity={coefficients[2],coefficients[3],coefficients[4]};for(int j=0;j<3;++j)y[f.l][j+1]-=derivative.traction[j]*f.area;y[f.l][4]-=dot3(derivative.traction,velocity)*f.area;}
  for(auto&f:mesh.faces)if(f.r>=0){D3 dx=sub3(mesh.center[f.r],mesh.center[f.l]);std::array<D3,4>g{};
   for(int k=0;k<4;++k){for(int j=0;j<3;++j)g[k][j]=.5*(gd[f.l][3*k+j]+gd[f.r][3*k+j]);double correction=-dot3(g[k],dx)/f.dn;for(int j=0;j<3;++j)g[k][j]+=correction*f.normal[j];}
   auto coefficients=flow_diffusion_cache_active?flow_diffusion_face_cache[static_cast<size_t>(&f-mesh.faces.data())]:flow_diffusion_face_coefficients(f);double mu=coefficients[0],div=g[0][0]+g[1][1]+g[2][2];D3 traction{};
   for(int j=0;j<3;++j)for(int k=0;k<3;++k)traction[j]+=mu*(g[j][k]+g[k][j]-(j==k?2*div/3:0))*f.normal[k];
   double conductivity=coefficients[1],energy=conductivity*dot3(g[3],f.normal);
   for(int j=0;j<3;++j)energy+=traction[j]*coefficients[2+j];
   // Viscous work also differentiates its face velocity. This term matters when
   // the initial thin-layer shear is large; the sparse normal block omits it.
   for(int j=0;j<3;++j)energy+=coefficients[5+j]*.5*(dw[f.l][j]+dw[f.r][j]);
   for(int j=0;j<3;++j){double value=-traction[j]*f.area;y[f.l][j+1]+=value;if(f.r<mesh.n)y[f.r][j+1]-=value;}y[f.l][4]-=energy*f.area;if(f.r<mesh.n)y[f.r][4]+=energy*f.area;
  }
 }
 double scalar_transport(const Face&f,int variable)const{
  double s1=1/mesh.meta[variable==5?29:31],s2=1/mesh.meta[variable==5?30:32],blend=.5*(prop[f.l][1]+prop[f.r][1]),sigma=blend*s1+(1-blend)*s2;
  if(!options.harmonic_scalars)return gas.mu+sigma*.5*(prop[f.l][0]+prop[f.r][0]);double dl=gas.mu+(prop[f.l][1]*s1+(1-prop[f.l][1])*s2)*prop[f.l][0],dr=gas.mu+(prop[f.r][1]*s1+(1-prop[f.r][1])*s2)*prop[f.r][0];return 2*dl*dr/(dl+dr);
 }
 void extra_scalar_diffusion(const Field<1>&x,Field<1>&y,int variable)const{
  auto values=x;auto gd=linear_physical_gradients(values,variable);
  for(auto&f:mesh.faces)if(f.r>=0){D3 dx=sub3(mesh.center[f.r],mesh.center[f.l]),g{};for(int j=0;j<3;++j)g[j]=.5*(gd[f.l][j]+gd[f.r][j]);double derivative=dot3(g,f.normal)-dot3(g,dx)/f.dn;
   double value=-scalar_transport(f,variable)*f.area*derivative;
   if(!(variable==6&&options.cell_omega_wall&&wall_weight[f.l]>0))y[f.l][0]+=value;if(f.r<mesh.n&&!(variable==6&&options.cell_omega_wall&&wall_weight[f.r]>0))y[f.r][0]-=value;
  }
 }
 Field<5>interior_viscous_residual(const std::vector<std::array<double,2>>&transport){Field<5>r(mesh.n);for(size_t fi=0;fi<mesh.faces.size();++fi){auto&f=mesh.faces[fi];if(f.r<0){if(f.type==7&&options.symmetry_viscous_flux){auto v=symmetry_viscous(f,transport[fi][0]);for(int j=0;j<3;++j)r[f.l][j+1]-=v.traction[j]*f.area;r[f.l][4]-=v.energy*f.area;}continue;}std::array<D3,3>g={face_gradient(f,1),face_gradient(f,2),face_gradient(f,3)};double div=g[0][0]+g[1][1]+g[2][2],energy=transport[fi][1]*dot3(face_gradient(f,4),f.normal);D3 traction{};for(int j=0;j<3;++j){for(int k=0;k<3;++k)traction[j]+=transport[fi][0]*(g[j][k]+g[k][j]-(j==k?2*div/3:0))*f.normal[k];energy+=traction[j]*.5*(w[f.l][j+1]+w[f.r][j+1]);double value=-traction[j]*f.area;r[f.l][j+1]+=value;if(f.r<mesh.n)r[f.r][j+1]-=value;}r[f.l][4]-=energy*f.area;if(f.r<mesh.n)r[f.r][4]+=energy*f.area;}return r;}
 void check_scalar_diffusion_operator(int variable){
  gradients();properties();auto original=w,limiter=lo;Field<1>direction(mesh.nt),analytic(mesh.n);double eps=options.diagnostic_directional_step;std::vector<double>transport(mesh.faces.size());
  auto constrained=[&](int i){return variable==6&&options.cell_omega_wall&&wall_weight[i]>0;};
  for(int i=0;i<mesh.n;++i)if(!constrained(i)){auto c=mesh.center[i];direction[i][0]=w[i][variable]*.02*std::sin(c[0]+.7*c[1]+2*c[2]);}mesh.exchange(direction);
  for(size_t fi=0;fi<mesh.faces.size();++fi){auto&f=mesh.faces[fi];if(f.r<0)continue;transport[fi]=scalar_transport(f,variable);double value=transport[fi]*f.area/f.dn*(direction[f.l][0]-direction[f.r][0]);if(!constrained(f.l))analytic[f.l][0]+=value;if(f.r<mesh.n&&!constrained(f.r))analytic[f.r][0]-=value;}
  extra_scalar_diffusion(direction,analytic,variable);std::array<Field<1>,2>sample;
  for(int sign=0;sign<2;++sign){w=original;double e=sign==0?eps:-eps;for(int i=0;i<mesh.n;++i)w[i][variable]+=e*direction[i][0];gradients(&limiter);sample[sign].resize(mesh.n);for(size_t fi=0;fi<mesh.faces.size();++fi){auto&f=mesh.faces[fi];if(f.r<0)continue;double value=-transport[fi]*f.area*dot3(face_gradient(f,variable),f.normal);if(!constrained(f.l))sample[sign][f.l][0]+=value;if(f.r<mesh.n&&!constrained(f.r))sample[sign][f.r][0]-=value;}}
  double error=0,reference=0;for(int i=0;i<mesh.n;++i){double fd=(sample[0][i][0]-sample[1][i][0])/(2*eps);error+=std::pow(fd-analytic[i][0],2);reference+=fd*fd;}double relative=std::sqrt(world.sum(error)/world.sum(reference));w=original;gradients();properties();if(world.rank==0)std::cout<<"M6_SCALAR_DIFFUSION_JACOBIAN variable="<<variable<<" relative_L2="<<relative<<std::endl;if(relative>3e-6)throw std::runtime_error("scalar full-gradient directional check failed");
 }
 void check_viscous_operator(){
  gradients();properties();auto original=w,limiter=lo;Field<5>x(mesh.nt);double eps=options.diagnostic_directional_step;std::vector<std::array<double,2>>transport(mesh.faces.size());for(size_t fi=0;fi<mesh.faces.size();++fi)transport[fi]=face_transport(mesh.faces[fi]);
  auto direction=[&](int i){const auto&c=mesh.center[i];return std::array<double,4>{.1*std::sin(c[0]+.3*c[2]),.1*std::cos(c[0]+c[1]),.1*std::sin(.7*c[2]+c[1]),.05*std::sin(c[0]+c[1]+c[2])};};
  for(int i=0;i<mesh.n;++i){auto d=direction(i);auto p=w[i],m=w[i];for(int j=0;j<3;++j){p[j+1]+=eps*d[j];m[j+1]-=eps*d[j];}p[4]+=eps*d[3];m[4]-=eps*d[3];auto qp=gas.conserved(p),qm=gas.conserved(m);for(int j=0;j<5;++j)x[i][j]=(qp[j]-qm[j])/(2*eps);}
  mesh.exchange(x);Field<5>analytic(mesh.n);for(auto&f:mesh.faces)if(f.r>=0){auto coefficients=face_transport(f);double mu=coefficients[0],conduct=coefficients[1];D3 u{};for(int j=0;j<3;++j)u[j]=.5*(w[f.l][j+1]+w[f.r][j+1]);auto kl=viscous_jacobian(f,w[f.l],mu,conduct,u),kr=viscous_jacobian(f,w[f.r],mu,conduct,u);auto al=bmultiply<5>(kl,x[f.l]),ar=bmultiply<5>(kr,x[f.r]);for(int j=0;j<5;++j){double t=al[j]-ar[j];analytic[f.l][j]+=t;if(f.r<mesh.n)analytic[f.r][j]-=t;}}
  if(options.symmetry_viscous_flux)for(auto&f:mesh.faces)if(f.type==7){auto block=symmetry_viscous_jacobian(f,face_transport(f)[0]);auto value=bmultiply<5>(block,x[f.l]);for(int j=0;j<5;++j)analytic[f.l][j]+=value[j];}
  prepare_flow_diffusion_cache();extra_flow_diffusion(x,analytic);flow_diffusion_cache_active=false;std::array<Field<5>,2>samples;
  for(int sign=0;sign<2;++sign){w=original;double e=sign==0?eps:-eps;for(int i=0;i<mesh.n;++i){auto d=direction(i);for(int j=0;j<3;++j)w[i][j+1]+=e*d[j];w[i][4]+=e*d[3];}gradients(&limiter);samples[sign]=interior_viscous_residual(transport);}
  double error=0,reference=0;for(int i=0;i<mesh.n;++i)for(int j=0;j<5;++j){double fd=(samples[0][i][j]-samples[1][i][j])/(2*eps);error+=std::pow(fd-analytic[i][j],2);reference+=fd*fd;}
  double relative=std::sqrt(world.sum(error)/world.sum(reference));w=original;gradients();if(world.rank==0)std::cout<<"M6_VISCOUS_JACOBIAN relative_L2="<<relative<<std::endl;if(relative>3e-6)throw std::runtime_error("full viscous-gradient Jacobian directional check failed");
  check_scalar_diffusion_operator(5);check_scalar_diffusion_operator(6);
 }
 void assemble(bool matrix){
  double started=now();std::fill(residual.begin(),residual.end(),Q5{});std::fill(pseudo.begin(),pseudo.end(),0);pressure_force={};viscous_force={};if(matrix)flow.clear();
  std::vector<double>time_diffusivity;if(options.normal_coupled_timestep){time_diffusivity.resize(mesh.nt);for(int i=0;i<mesh.nt;++i)time_diffusivity[i]=coupled_time_diffusivity(gas.rho(w[i]),gas.cp(w[i][4]),gas.mu,prop[i][0],gas.conductivity,mesh.meta[26],.7,2.88e-5/gas.a0);}
#ifdef M6_CELL_JACOBIAN_CACHE
  // Scope contract: W[0..4], Gas/material and mesh indexing are immutable
  // throughout this assemble invocation. Include synchronized ghost cells.
  // Rebuild only for matrix assembly; no cache survives the function return.
  std::vector<CellConservativeJacobian>cell_jacobians;
  if(matrix){cell_jacobians.reserve(mesh.nt);for(int i=0;i<mesh.nt;++i)cell_jacobians.push_back(capture_cell_conservative_jacobian(gas,w[i]));}
#endif
  for(size_t fi=0;fi<mesh.faces.size();++fi){auto&f=mesh.faces[fi];int l=f.l,r=f.r;auto fl=flux(f);massflux[fi]=fl[0]*f.area;
   for(int j=0;j<5;++j){residual[l][j]+=fl[j]*f.area;if(r>=0&&r<mesh.n)residual[r][j]-=fl[j]*f.area;}
   // Exact/normal convective blocks replace the old physical-flux Jacobians.
   // Normal pseudo-time also uses its own spectrum; retain the old work only
   // where a fallback matrix or the selected pseudo-time actually consumes it.
   bool fallback_jacobian=r>=0?!options.exact_flux_jacobian:(f.type!=3&&f.type!=7&&!(f.type==9&&options.characteristic_farfield));
   bool need_radius=(!options.normal_coupled_timestep&&!options.local_preconditioning)||(matrix&&fallback_jacobian);
   double radius=0.;if(!options.skip_unused_face_work||need_radius){double un=w[l][1]*f.normal[0]+w[l][2]*f.normal[1]+w[l][3]*f.normal[2];radius=std::abs(un)+gas.sound(w[l]);
    if(r>=0){double unr=w[r][1]*f.normal[0]+w[r][2]*f.normal[1]+w[r][3]*f.normal[2];radius=std::max(radius,std::abs(unr)+gas.sound(w[r]));}}
   std::array<double,2>transport{};if(!options.skip_unused_face_work||matrix||!options.normal_coupled_timestep)transport=face_transport(f);double mu=transport[0],conduct=transport[1];
   double time_radius=radius;if(options.local_preconditioning&&(!options.skip_unused_face_work||!options.normal_coupled_timestep)){time_radius=gas.preconditioned_radius(w[l],f.normal,reference_velocity[l][0]);if(r>=0)time_radius=std::max(time_radius,gas.preconditioned_radius(w[r],f.normal,reference_velocity[r][0]));}
   // Each owned cell uses its own density. A cut face can be locally reversed,
   // so sharing the left-cell kinematic diffusivity makes pseudo-time rank dependent.
   if(options.normal_coupled_timestep){
    double rl=gas.rho(w[l]),rr=r>=0?gas.rho(w[r]):rl,momentum=0.;for(int j=0;j<3;++j)momentum+=(r>=0?.5*(rl*w[l][j+1]+rr*w[r][j+1]):rl*w[l][j+1])*f.normal[j]*f.area;
    double distance=r>=0?dot3(sub3(mesh.center[r],mesh.center[l]),f.normal):2*dot3(sub3(f.center,mesh.center[l]),f.normal);
    double visc=coupled_time_viscous_radius(r>=0?std::max(time_diffusivity[l],time_diffusivity[r]):time_diffusivity[l],f.area,distance);
    pseudo[l]+=(coupled_time_convective_radius(momentum,rl,reference_velocity[l][0],std::pow(gas.sound(w[l]),2),f.area*f.area)+visc)/(2*options.cfl);
    if(r>=0&&r<mesh.n)pseudo[r]+=(coupled_time_convective_radius(momentum,rr,reference_velocity[r][0],std::pow(gas.sound(w[r]),2),f.area*f.area)+visc)/(2*options.cfl);
   }else{double viscous_weight=4*mu*f.area*dot3(diffusion_direction(f),f.normal),convective_weight=time_radius*f.area;pseudo[l]+=(convective_weight+viscous_weight/gas.rho(w[l]))/options.cfl;if(r>=0&&r<mesh.n)pseudo[r]+=(convective_weight+viscous_weight/gas.rho(w[r]))/options.cfl;}
   if(!matrix)continue;
   B5 al{},ar{};if(!options.skip_unused_face_work||fallback_jacobian)al=gas.jacobian(w[l],f.normal);D3 uf={w[l][1],w[l][2],w[l][3]};if(r>=0)for(int j=0;j<3;++j)uf[j]=.5*(uf[j]+w[r][j+1]);
#ifdef M6_CELL_JACOBIAN_CACHE
   B5 kl=viscous_jacobian_cached(f,cell_jacobians[l].primitive,mu,conduct,uf),kr{};
#else
   B5 kl=viscous_jacobian(f,w[l],mu,conduct,uf),kr{};
#endif
   if(f.type==9&&options.farfield_zero_viscous_flux)kl={};
   if(r>=0){if(!options.skip_unused_face_work||!options.exact_flux_jacobian)ar=gas.jacobian(w[r],f.normal);
#ifdef M6_CELL_JACOBIAN_CACHE
    kr=viscous_jacobian_cached(f,cell_jacobians[r].primitive,mu,conduct,uf);
#else
    kr=viscous_jacobian(f,w[r],mu,conduct,uf);
#endif
    if(options.exact_flux_jacobian){auto thermo=options.independent_face_thermo?face_thermodynamics(f):std::array<double,4>{};auto left=reconstruct(l,f),right=reconstruct(r,f);double a2_l=std::pow(gas.sound(w[l]),2),a2_r=std::pow(gas.sound(w[r]),2);auto jac=options.normal_approximate_jacobian?
#ifdef M6_CELL_JACOBIAN_CACHE
     coupled_approximate_jacobian_cached(gas,w[l],w[r],left,right,f.normal,fl[0],a2_l,a2_r,reference_velocity[l][0],reference_velocity[r][0],cell_jacobians[l],cell_jacobians[r],&thermo)
#else
     coupled_approximate_jacobian(gas,w[l],w[r],left,right,f.normal,fl[0],a2_l,a2_r,reference_velocity[l][0],reference_velocity[r][0],&thermo)
#endif
     :coupled_flux_jacobian(gas,w[l],w[r],left,right,f.normal,a2_l,a2_r,reference_velocity[l][0],reference_velocity[r][0],options.independent_face_thermo?&thermo:nullptr);for(int j=0;j<25;++j){al[j]=jac.left[j]*f.area+kl[j];ar[j]=jac.right[j]*f.area-kr[j];}}
    else{for(int j=0;j<25;++j){al[j]*=.5*f.area;ar[j]*=.5*f.area;al[j]+=kl[j];ar[j]-=kr[j];}for(int j=0;j<5;++j){al[6*j]+=.5*radius*f.area;ar[6*j]-=.5*radius*f.area;}}
    for(int j=0;j<25;++j){flow.a[flow.diag[l]][j]+=al[j];flow.a[f.lr][j]+=ar[j];if(r<mesh.n){flow.a[f.rl][j]-=al[j];flow.a[flow.diag[r]][j]-=ar[j];}}
   }else if(f.type==3||f.type==7){
    B5 b{};double rhol=gas.rho(w[l]),c=gas.sound(w[l]);auto d=gas.primitive_jacobian(w[l]);
    auto face=reconstruct(l,f);auto exact=options.exact_flux_jacobian?(((f.type==3&&options.normal_wall_jacobian)||(f.type==7&&options.normal_symmetry_jacobian))?normal_wall_pressure_jacobian(gas,w[l],f.normal,reference_velocity[l][0],boundary_face_density(f,face)):coupled_wall_pressure_jacobian(gas,w[l],face,f.normal,reference_velocity[l][0],boundary_face_density(f,face),options.independent_face_thermo)):Q5{};double gamma=gas.cp(w[l][4])/(gas.cp(w[l][4])-gas.rg);for(int j=0;j<5;++j){double dp=(gamma-1)*((j==4?1.:0.)-(j>=1&&j<=3?w[l][j]:0.)+(j==0?.5*(w[l][1]*w[l][1]+w[l][2]*w[l][2]+w[l][3]*w[l][3]):0.));double dun=0;for(int k=0;k<3;++k)dun+=f.normal[k]*d[5*(k+1)+j];double pressure_derivative=options.exact_flux_jacobian?exact[j]:dp+rhol*c*dun;for(int k=0;k<3;++k)b[5*(k+1)+j]=pressure_derivative*f.normal[k]*f.area;}
    if(f.type==3){double coeff=wall(f).velocity_coefficient*f.area;for(int j=0;j<5;++j){double dn=0;for(int k=0;k<3;++k)dn+=f.normal[k]*d[5*(k+1)+j];for(int k=0;k<3;++k)b[5*(k+1)+j]+=coeff*(d[5*(k+1)+j]-f.normal[k]*dn);}}
    if(f.type==7&&options.symmetry_viscous_flux){auto viscous=symmetry_viscous_jacobian(f,mu);for(int j=0;j<25;++j)b[j]+=viscous[j];}
    for(int j=0;j<25;++j)flow.a[flow.diag[l]][j]+=b[j];
   }else if(f.type==9&&options.characteristic_farfield){auto face=reconstruct(l,f);double density=boundary_face_density(f,face);auto b=options.normal_farfield_jacobian?normal_farfield_jacobian(gas,w[l],face,freestream(),f.normal,density):characteristic_farfield_jacobian(gas,w[l],face,freestream(),f.normal,density,options.independent_face_thermo);for(int j=0;j<25;++j)flow.a[flow.diag[l]][j]+=b[j]*f.area+kl[j];}
   else{for(int j=0;j<25;++j)flow.a[flow.diag[l]][j]+=.5*f.area*al[j]+kl[j];for(int j=0;j<5;++j)flow.a[flow.diag[l]][6*j]+=.5*radius*f.area;}
  }
  for(int i=0;i<mesh.n;++i)if(matrix){if(options.local_preconditioning){auto p=gas.preconditioning_matrix(w[i],reference_velocity[i][0]);for(int j=0;j<25;++j)flow.a[flow.diag[i]][j]+=pseudo[i]*p[j];}else for(int j=0;j<5;++j)flow.a[flow.diag[i]][6*j]+=pseudo[i];}
  pressure_force=world.sum(pressure_force);viscous_force=world.sum(viscous_force);for(int j=0;j<3;++j)force[j]=pressure_force[j]+viscous_force[j];assembly_time+=now()-started;
 }
 int meanflow(bool factor,bool audit=false){
  double started=now();prepare_flow_diffusion_cache();if(factor)flow.factor();Field<5>b=residual,delta;for(auto&v:b)for(double&t:v)t=-t;std::function<void(const Field<5>&,Field<5>&)>extra;if(options.implicit_gradients)extra=[this](const Field<5>&x,Field<5>&y){extra_flow_diffusion(x,y);};bool attempt=recycle_flow_direction&&!recycled_flow_direction.empty();int steps=flow.gmres(b,delta,options.linear_steps,options.linear_tolerance,extra,options.batched_orthogonalization,&linear_estimate,attempt?&recycled_flow_direction:nullptr,&linear_start);if(attempt){++recycled_attempts;recycled_initial_ratio_sum+=linear_start.relative_residual;if(linear_start.recycled)++recycled_accepted;}
  linear_true_residual=-1;if(audit){Field<5>applied;flow.multiply(delta,applied);if(extra)extra(delta,applied);for(int i=0;i<mesh.n;++i)for(int j=0;j<5;++j)applied[i][j]-=b[i][j];double bn=std::sqrt(flow.dot(b,b));linear_true_residual=std::sqrt(flow.dot(applied,applied))/std::max(bn,1e-30);if(!std::isfinite(linear_true_residual)||!std::isfinite(linear_estimate))throw std::runtime_error("nonfinite audited linear residual");maximum_linear_true_residual=std::max(maximum_linear_true_residual,linear_true_residual);++linear_audits;if(linear_true_residual>options.linear_tolerance*(1+1e-8))++linear_target_missed;}
  flow_diffusion_cache_active=false;if(recycle_flow_direction)recycled_flow_direction.assign(delta.begin(),delta.begin()+mesh.n);
  for(int i=0;i<mesh.n;++i){double scale=options.flow_relax;bool accepted=false;Q5 primitive_delta{};if(options.primitive_flow_update)primitive_delta=linearized_primitive_increment(gas,w[i],delta[i]);
   for(int trial=0;trial<18;++trial){Q5 candidate=q[i];for(int j=0;j<5;++j)candidate[j]+=scale*delta[i][j];
    try{W7 v;if(options.primitive_flow_update){v=primitive_flow_trial(w[i],primitive_delta,scale);bool finite=true;for(int j=0;j<5;++j)finite=finite&&std::isfinite(v[j]);if(!finite||!(v[0]>0&&v[4]>0))throw std::runtime_error("inadmissible primitive flow trial");candidate=gas.conserved(v);}else v=gas.primitive(candidate,w[i][4],w[i][5],w[i][6]);if(options.case_flow_update_limits?within_case_flow_update_limits(gas,v,flow_update_limits):(v[4]>=.1&&v[4]<15&&v[0]>.01*freestream()[0]&&v[4]>=.8*w[i][4]&&v[4]<=1.2*w[i][4]&&candidate[0]>=.5*q[i][0]&&candidate[0]<=2*q[i][0])){q[i]=candidate;w[i]=v;accepted=true;break;}}
    catch(const std::runtime_error&){}scale*=.5;
   }if(!accepted){std::ostringstream s;s<<"mean-flow positivity line search exhausted cell="<<mesh.global[i]<<" w=";for(double v:w[i])s<<v<<',';s<<" delta=";for(double v:delta[i])s<<v<<',';throw std::runtime_error(s.str());}if(scale<options.flow_relax)++limited_updates;
  }linear_time+=now()-started;return steps;
 }
 void export_linear_operator_observation(const std::function<void()>&diagnostic_property_update={},bool already_assembled=false){
  if(!already_assembled){gradients();properties();if(diagnostic_property_update)diagnostic_property_update();assemble(true);}Field<5>direction(mesh.nt),applied,rhs=residual;
  for(int i=0;i<mesh.n;++i)for(int j=0;j<5;++j){direction[i][j]=(.1+.02*j)*std::sin(.013*(mesh.global[i]+1)+.4*j);rhs[i][j]=-rhs[i][j];}
  flow.multiply(direction,applied);if(options.implicit_gradients)extra_flow_diffusion(direction,applied);
  std::array<double,10>moments{};for(int i=0;i<mesh.n;++i)for(int j=0;j<5;++j){if(!std::isfinite(rhs[i][j])||!std::isfinite(applied[i][j]))throw std::runtime_error("nonfinite frozen operator observation");moments[j]+=rhs[i][j]*rhs[i][j];moments[j+5]+=applied[i][j]*applied[i][j];}moments=world.sum(moments);
  std::ostringstream name;name<<"linear-operator-rank-"<<std::setfill('0')<<std::setw(3)<<world.rank<<".bin";std::ofstream file(std::filesystem::path(options.output)/name.str(),std::ios::binary);int32_t header[4]={mesh.n,world.rank,world.size,static_cast<int32_t>(mesh.global_cells)};file.write("M6OPR001",8);file.write(reinterpret_cast<const char*>(header),sizeof(header));file.write(reinterpret_cast<const char*>(mesh.global.data()),mesh.n*sizeof(int));file.write(reinterpret_cast<const char*>(rhs.data()),mesh.n*sizeof(Q5));file.write(reinterpret_cast<const char*>(applied.data()),mesh.n*sizeof(Q5));if(!file)throw std::runtime_error("frozen operator observation write failed");
  if(world.rank==0){std::ofstream out(std::filesystem::path(options.output)/"operator-observation.json");out<<std::setprecision(16)<<std::boolalpha<<"{\n\"diagnostic_only\":true,\n\"integration_steps\":0,\n\"format\":\"M6OPR001\",\n\"cells\":"<<mesh.global_cells<<",\n\"ranks\":"<<world.size<<",\n\"restart_iteration\":"<<iteration_offset<<",\n\"implicit_viscous_gradients\":"<<options.implicit_gradients<<",\n\"rhs_component_L2\":[";for(int j=0;j<5;++j){if(j)out<<',';out<<std::sqrt(moments[j]);}out<<"],\n\"operator_application_component_L2\":[";for(int j=0;j<5;++j){if(j)out<<',';out<<std::sqrt(moments[j+5]);}out<<"]\n}\n";std::cout<<"M6_FROZEN_OPERATOR_OBSERVATION ranks="<<world.size<<" integration_steps=0"<<std::endl;}
 }
 struct LinearPair{std::array<double,2>seconds{},true_residual{};std::array<int,2>steps{};double solution_difference=0;};
 template<int N>LinearPair compare_linear_paths(Sparse<N>&matrix,const Field<N>&b,int maxsteps,double tolerance,const std::function<void(const Field<N>&,Field<N>&)>&extra,bool reverse){
  LinearPair result;std::array<Field<N>,2>solution;
  double bn=std::sqrt(matrix.dot(b,b));
  for(int position=0;position<2;++position){int mode=reverse?1-position:position;MPI_Barrier(MPI_COMM_WORLD);double start=now();result.steps[mode]=matrix.gmres(b,solution[mode],maxsteps,tolerance,extra,mode==1);result.seconds[mode]=world.max(now()-start);
   Field<N>applied;matrix.multiply(solution[mode],applied);if(extra)extra(solution[mode],applied);for(int i=0;i<mesh.n;++i)for(int j=0;j<N;++j)applied[i][j]-=b[i][j];result.true_residual[mode]=std::sqrt(matrix.dot(applied,applied))/std::max(bn,1e-30);
  }
  std::array<double,2*N>difference{};for(int i=0;i<mesh.n;++i)for(int j=0;j<N;++j){difference[j]+=std::pow(solution[0][i][j]-solution[1][i][j],2);difference[N+j]+=solution[0][i][j]*solution[0][i][j];}difference=world.sum(difference);
  for(int j=0;j<N;++j)result.solution_difference=std::max(result.solution_difference,std::sqrt(difference[j]/std::max(difference[N+j],1e-60)));
  if(result.steps[0]!=result.steps[1]||result.solution_difference>1e-7||std::abs(result.true_residual[0]-result.true_residual[1])>1e-8*std::max(1.,result.true_residual[0]))throw std::runtime_error("batched CGS2 failed real-matrix equivalence/true-residual checks");
  return result;
 }
 void check_linear_solver(){
  gradients();properties();assemble(true);flow.factor();Field<5>b=residual;for(auto&row:b)for(double&x:row)x=-x;
  std::function<void(const Field<5>&,Field<5>&)>extra;if(options.implicit_gradients)extra=[this](const Field<5>&x,Field<5>&y){extra_flow_diffusion(x,y);};
  std::array<LinearPair,3>flow_results;for(int repeat=0;repeat<3;++repeat)flow_results[repeat]=compare_linear_paths<5>(flow,b,options.linear_steps,options.linear_tolerance,extra,repeat%2!=0);
  std::array<LinearPair,2>scalar_results;
  for(int variable=5;variable<=6;++variable){turbulence(variable,false);scalar.factor();Field<1>known(mesh.nt),rhs;for(int i=0;i<mesh.n;++i)known[i][0]=1+.1*std::sin(.017*(mesh.global[i]+1));scalar.multiply(known,rhs);
   std::function<void(const Field<1>&,Field<1>&)>scalar_extra;if(options.implicit_gradients){scalar_extra=[this,variable](const Field<1>&x,Field<1>&y){extra_scalar_diffusion(x,y,variable);};scalar_extra(known,rhs);}
   scalar_results[variable-5]=compare_linear_paths<1>(scalar,rhs,std::max(4,options.linear_steps/2),.1,scalar_extra,variable==6);
  }
  if(world.rank==0){std::ofstream out(std::filesystem::path(options.output)/"linear-solver-validation.json");out<<std::setprecision(16)<<std::boolalpha<<"{\n\"diagnostic_only\":true,\n\"integration_steps\":0,\n\"cells\":"<<mesh.global_cells<<",\n\"ranks\":"<<world.size<<",\n\"restart_iteration\":"<<iteration_offset<<",\n\"implicit_viscous_gradients\":"<<options.implicit_gradients<<",\n\"timing_scope\":\"frozen factored real operator: solve only, excludes assembly, factorization and true-residual verification; research timing\",\n\"flow_pairs\":[";
   auto write=[&](const LinearPair&r){out<<"{\"mgs2_seconds\":"<<r.seconds[0]<<",\"cgs2_seconds\":"<<r.seconds[1]<<",\"mgs2_true_relative_residual\":"<<r.true_residual[0]<<",\"cgs2_true_relative_residual\":"<<r.true_residual[1]<<",\"mgs2_steps\":"<<r.steps[0]<<",\"cgs2_steps\":"<<r.steps[1]<<",\"max_component_relative_solution_difference\":"<<r.solution_difference<<'}';};
   for(int i=0;i<3;++i){if(i)out<<',';write(flow_results[i]);}out<<"],\n\"scalar_pairs\":[";for(int i=0;i<2;++i){if(i)out<<',';write(scalar_results[i]);}out<<"],\n\"status\":\"equivalence checks passed; not a full-CFD performance acceptance\"\n}\n";
   std::cout<<"M6_BATCHED_CGS2_VALIDATION real flow/k/omega matrices passed, integration_steps=0"<<std::endl;
  }
 }
 int turbulence(int variable,bool advance=true,bool audit=false,Field<1>*assembly_rhs=nullptr,bool normal_scalar_diagnostic=false){
  normal_scalar_diagnostic=normal_scalar_diagnostic||options.normal_sst_transport;
  double started=now();std::vector<double>constrained(mesh.n,0.);
  if(variable==6&&options.cell_omega_wall){
   for(auto&f:mesh.faces)if(f.type==3)constrained[f.l]=std::max(constrained[f.l],wall(f).omega);
   if(advance)for(int i=0;i<mesh.n;++i)if(constrained[i]>0)w[i][6]=.8*constrained[i]+.2*w[i][6];
   // TRB_CELL_VALUE constrains the adjacent cell; the face value is extrapolated.
   if(advance){if(options.normal_sst_transport){mesh.exchange(w);sst_transport->refresh_wall_omega(w);properties(true,true);}else{gradients();properties();assemble(false);}}
  }
  std::vector<double>wall_k_omega;if(normal_scalar_diagnostic&&variable==5){wall_k_omega.assign(mesh.n,0.);for(auto&f:mesh.faces)if(f.type==3)wall_k_omega[f.l]=std::max(wall_k_omega[f.l],wall(f).omega);}
  scalar.clear();Field<1>r(mesh.n),delta;
  for(size_t fi=0;fi<mesh.faces.size();++fi){auto&f=mesh.faces[fi];int l=f.l,rr=f.r;double mf=massflux[fi],rl=0,df=0,cl=0,cr=0;
   if(rr>=0){double value=mf>=0?reconstruct_scalar(l,f,variable):reconstruct_scalar(rr,f,variable);double diffusion=scalar_transport(f,variable);auto g=face_gradient(f,variable);rl=mf*value-diffusion*f.area*dot3(g,f.normal);df=diffusion*f.area*dot3(diffusion_direction(f),f.normal);cl=std::max(mf,0.)+df;cr=std::min(mf,0.)-df;
   }else if(f.type==3){if(variable==6&&!options.cell_omega_wall){double boundary=wall(f).omega;df=gas.mu*f.area/f.dn;rl=df*(w[l][variable]-boundary);cl=df;}}
   else if(f.type==9){rl=mf*(mf>=0?(normal_scalar_diagnostic?w[l][variable]:reconstruct_scalar(l,f,variable)):freestream()[variable]);cl=std::max(mf,0.);}
   r[l][0]+=rl;scalar.a[scalar.diag[l]][0]+=cl;
   if(rr>=0){scalar.a[f.lr][0]+=cr;if(rr<mesh.n){r[rr][0]-=rl;scalar.a[f.rl][0]-=cl;scalar.a[scalar.diag[rr]][0]-=cr;}}
  }
  double residual_norm=0;
  for(int i=0;i<mesh.n;++i){double rho=gas.rho(w[i]),vol=mesh.volume[i],source=0,negative_derivative=0;
   if(variable==5){double omega=normal_scalar_diagnostic&&wall_k_omega[i]>0?wall_k_omega[i]:w[i][6];source=prop[i][4]-mesh.meta[34]*rho*w[i][5]*omega;negative_derivative=mesh.meta[34]*rho*omega;}
   else{auto terms=sst_omega_source(rho,prop[i][0],gas.mu,prop[i][4],w[i][6],prop[i][2],prop[i][3],prop[i][5],options.omega_laminar_factor,options.cross_diffusion_limit,options.linearize_negative_cross);source=terms.value;negative_derivative=terms.negative_derivative;}
   // Primitive transport subtracts phi times the mass residual. This maintains
   // a uniform passive species exactly even during nonzero continuity residual.
   r[i][0]-=source*vol+w[i][variable]*residual[i][0];
   if(constrained[i]>0){for(int p=scalar.row[i];p<scalar.row[i+1];++p)scalar.a[p][0]=0.;scalar.a[scalar.diag[i]][0]=1.;r[i][0]=advance?0:constrained[i]-w[i][6];continue;}
   residual_norm+=std::abs(r[i][0]);
    double temporal=options.steady_turbulence?0:rho*pseudo[i];if(options.normal_scalar_step)temporal=0.;scalar.a[scalar.diag[i]][0]+=temporal+negative_derivative*vol-residual[i][0];
   double offsum=0;for(int p=scalar.row[i];p<scalar.row[i+1];++p)if(p!=scalar.diag[i])offsum+=std::abs(scalar.a[p][0]);
    scalar.a[scalar.diag[i]][0]=std::max(scalar.a[scalar.diag[i]][0],offsum+temporal+1e-30);if(options.steady_turbulence&&!options.normal_scalar_step)scalar.a[scalar.diag[i]][0]/=options.turb_relax;r[i][0]=-r[i][0];
  }
   // A_spatial*old+r is the deferred absolute source. Preserve the physical
   // increment RHS while folding negative source into the diagonal, then add
   // the retained pre-flow pseudo-time. Normal scalar AMG uses the face graph.
   if(options.normal_scalar_step){
    if(scalar_dt.size()!=static_cast<size_t>(mesh.n))throw std::runtime_error("missing pre-flow scalar local times");
    int negative_sources=0;for(int i=0;i<mesh.n;++i){if(constrained[i]>0)continue;double source=r[i][0];for(int p=scalar.row[i];p<scalar.row[i+1];++p)source+=scalar.a[p][0]*w[scalar.col[p]][variable];
     if(source<0)++negative_sources;auto adjusted=linearize_negative_scalar_source(w[i][variable],scalar.a[scalar.diag[i]][0],source);
     scalar.a[scalar.diag[i]][0]=adjusted.diagonal+scalar_time_coefficient(gas.rho(w[i]),mesh.volume[i],scalar_dt[i]);
    }scalar_negative_source_rows[variable-5]+=static_cast<long long>(world.sum(static_cast<double>(negative_sources)));
   }
  if(assembly_rhs)*assembly_rhs=r;
  if(!advance){double sum=0,floor=.01*mesh.initial[0][variable];for(int i=0;i<mesh.n;++i){double t=.1*r[i][0]/(scalar.a[scalar.diag[i]][0]*std::max(w[i][variable],floor));sum+=t*t;}scalar_merit[variable-5]=world.sum(sum);return 0;}
  norms[variable]=world.sum(residual_norm);scalar.factor();std::function<void(const Field<1>&,Field<1>&)>extra;if(options.implicit_gradients&&!options.normal_scalar_step)extra=[this,variable](const Field<1>&x,Field<1>&y){extra_scalar_diffusion(x,y,variable);};int slot=variable-5;int steps=scalar.gmres(r,delta,options.turbulence_linear_max_steps?options.turbulence_linear_max_steps:std::max(4,options.linear_steps/2),options.turbulence_linear_tolerance,extra,options.batched_orthogonalization,&turbulence_linear_estimate[slot]);turbulence_linear_steps[slot]=steps;turbulence_linear_true[slot]=-1;
  if(audit&&options.audit_turbulence_linear){Field<1>applied;scalar.multiply(delta,applied);if(extra)extra(delta,applied);for(int i=0;i<mesh.n;++i)applied[i][0]-=r[i][0];double relative=std::sqrt(scalar.dot(applied,applied))/std::max(std::sqrt(scalar.dot(r,r)),1e-30);if(!std::isfinite(relative)||!std::isfinite(turbulence_linear_estimate[slot]))throw std::runtime_error("nonfinite audited turbulence linear residual");turbulence_linear_true[slot]=relative;turbulence_linear_maximum[slot]=std::max(turbulence_linear_maximum[slot],relative);++turbulence_linear_audits[slot];if(relative>options.turbulence_linear_tolerance*(1+1e-8))++turbulence_linear_misses[slot];}
   if(options.normal_scalar_step){
    Field<1>absolute(mesh.nt);int negative=0;for(int i=0;i<mesh.n;++i){absolute[i][0]=w[i][variable]+delta[i][0];if(!std::isfinite(absolute[i][0]))throw std::runtime_error("nonfinite absolute scalar solution");if(constrained[i]<=0&&absolute[i][0]<0)++negative;}
    int total_negative=static_cast<int>(world.sum(static_cast<double>(negative)));scalar_negative_repairs[slot]+=total_negative;if(total_negative)mesh.exchange(absolute);
    double scale=variable==5?gas.a0*gas.a0:gas.a0,minimum=variable==5?1e-14:1e-20;
    for(int i=0;i<mesh.n;++i){if(constrained[i]>0)continue;double solution=absolute[i][0];
     // An independent, partition-stable positive-neighbor snapshot repair.
     // Native sequential in-place repair has not been runtime validated.
     if(solution<0){double sum=0.;int count=0;for(int p=scalar.row[i];p<scalar.row[i+1];++p){int j=scalar.col[p];if(j!=i&&absolute[j][0]>0){sum+=absolute[j][0];++count;}}solution=count?sum/count:minimum/scale;}
     w[i][variable]=clipped_relaxed_scalar(w[i][variable]*scale,solution*scale,minimum,options.turb_relax)/scale;
    }
   }else{for(int i=0;i<mesh.n;++i){if(constrained[i]>0)continue;double value=w[i][variable],change=(options.steady_turbulence?1:options.turb_relax)*delta[i][0];double lower=-.8*value,upper=variable==5?std::max(4*value,1e-5):std::max(1e6*value,1.);
     change=std::max(lower,std::min(change,upper));w[i][variable]=std::max(variable==5?1e-14:1e-10,value+change);
   }}mesh.exchange(w);turbulence_time+=now()-started;return steps;
 }
 bool may_skip_discarded_postflow_assembly(bool production_step)const{
#if defined(M6_SKIP_DISCARDED_POSTFLOW_ASSEMBLY) && M6_SKIP_DISCARDED_POSTFLOW_ASSEMBLY
  // Direct callers and diagnostics retain every assembly output. The proof is
  // only for the ordinary run() step, with fresh mean-flow assembly each time.
  return production_step&&options.normal_sst_transport&&options.normal_scalar_step&&options.lagged_sst_mass_flux
   &&options.normal_coupled_timestep&&options.anderson_depth==0&&options.jacobian_period==1
   &&!options.validate&&!options.check_viscous&&!options.check_linear&&!options.export_linear_operator;
#else
  (void)production_step;return false;
#endif
 }
 void validate_discarded_postflow_assembly()const{
  // Retain exactly the domain checks reachable from assemble(false), in the
  // original cell/face order. Its fluxes, spectra, five residuals and force
  // reductions have no live consumer on the guarded run() path. Reuse the
  // original throwing kernels rather than introduce independent thresholds.
  std::vector<double>time_diffusivity(mesh.nt);
  for(int i=0;i<mesh.nt;++i)time_diffusivity[i]=coupled_time_diffusivity(gas.rho(w[i]),gas.cp(w[i][4]),gas.mu,prop[i][0],gas.conductivity,mesh.meta[26],.7,2.88e-5/gas.a0);
  for(const auto&f:mesh.faces){int l=f.l,r=f.r;auto left=reconstruct(l,f);
   if(f.type==9&&options.characteristic_farfield){(void)m6::characteristic_farfield(gas,w[l],left,freestream(),f.normal,gas.rho(w[l]),std::pow(gas.sound(w[l]),2),boundary_face_density(f,left),gas.enthalpy(w[l]));}
   else if(f.type==3||f.type==7){if(options.coupled_ausm)(void)boundary_face_density(f,left);}
   else{if(r>=0)(void)reconstruct(r,f);if(options.coupled_ausm&&r>=0&&options.independent_face_thermo)(void)face_thermodynamics(f);}
   // Ordinary wall_model/viscous_flux/face_transport do not throw or mutate.
   if(f.type==7&&options.symmetry_viscous_flux)(void)symmetry_viscous(f,face_transport(f)[0]);
   double rl=gas.rho(w[l]),rr=r>=0?gas.rho(w[r]):rl,momentum=0.;
   for(int j=0;j<3;++j)momentum+=(r>=0?.5*(rl*w[l][j+1]+rr*w[r][j+1]):rl*w[l][j+1])*f.normal[j]*f.area;
   double distance=r>=0?dot3(sub3(mesh.center[r],mesh.center[l]),f.normal):2*dot3(sub3(f.center,mesh.center[l]),f.normal);
   (void)coupled_time_viscous_radius(r>=0?std::max(time_diffusivity[l],time_diffusivity[r]):time_diffusivity[l],f.area,distance);
   (void)coupled_time_convective_radius(momentum,rl,reference_velocity[l][0],std::pow(gas.sound(w[l]),2),f.area*f.area);
   if(r>=0&&r<mesh.n)(void)coupled_time_convective_radius(momentum,rr,reference_velocity[r][0],std::pow(gas.sound(w[r]),2),f.area*f.area);
  }
 }
 void advance_turbulence(bool audit=false,bool production_step=false){
  if(options.normal_sst_transport){
   // Normal F_FLUX is written by the mean-flow inviscid assembly, then held
   // through the primitive update and both scalar solves. Regeneration of
   // other post-flow coefficients must not silently replace that cache.
   const bool skip_postflow=may_skip_discarded_postflow_assembly(production_step);
   std::vector<double>transport_mass;if(options.lagged_sst_mass_flux&&!skip_postflow)transport_mass=massflux;
    if(options.normal_scalar_step){scalar_dt.resize(mesh.n);for(int i=0;i<mesh.n;++i)scalar_dt[i]=mesh.volume[i]/pseudo[i];}
   // Refresh post-flow gradients, but hold transport caches through k/omega.
   // The normal-k outputs at the end of gradients are overwritten after the
   // farfield cache refresh, before any consumer. Defer only this call site;
   // slot1, all other gradient callers, arithmetic and boundary order stay.
#if defined(M6_DEFER_POSTFLOW_K_RECONSTRUCTION)
   gradients(nullptr,true);
#else
   gradients();
#endif
   sst_transport->refresh_farfield(0,w,freestream());reconstruct_normal_scalar(0);properties(true,false);
   if(skip_postflow){double started=now();validate_discarded_postflow_assembly();assembly_time+=now()-started;++skipped_postflow_assemblies;}
   else assemble(false);
   if(options.lagged_sst_mass_flux&&!skip_postflow)massflux.swap(transport_mass);
   for(auto&r:residual)r[0]=0.;for(size_t fi=0;fi<mesh.faces.size();++fi){auto&f=mesh.faces[fi];double mf=massflux[fi];residual[f.l][0]+=mf;if(f.r>=0&&f.r<mesh.n)residual[f.r][0]-=mf;}
   sst_transport->refresh_wall_k(w);turbulence(5,true,audit);reconstruct_normal_scalar(0);
   sst_transport->refresh_farfield(1,w,freestream());reconstruct_normal_scalar(1);turbulence(6,true,audit);return;
  }
  if(options.turbulence_k_first){
   // Refresh the post-flow state before k. The existing omega solve then
   // imposes its adjacent-cell wall profile using the newly updated k.
   gradients();properties();assemble(false);turbulence(5,true,audit);turbulence(6,true,audit);
  }else{turbulence(6,true,audit);properties();turbulence(5,true,audit);}
 }
 std::vector<W7>acceleration_state()const{
  std::vector<W7>y(mesh.n);for(int i=0;i<mesh.n;++i){for(int j=0;j<5;++j)y[i][j]=q[i][j];y[i][5]=.1*std::log(w[i][5]/mesh.initial[i][5]);y[i][6]=.1*std::log(w[i][6]/mesh.initial[i][6]);}return y;
 }
 bool set_acceleration_state(const std::vector<W7>&y,const std::vector<W7>&base,const std::vector<W7>&original,double alpha){
  bool valid=true;for(int i=0;i<mesh.n;++i){Q5 trial{};for(int j=0;j<5;++j)trial[j]=base[i][j]+alpha*(y[i][j]-base[i][j]);double lk=base[i][5]+alpha*(y[i][5]-base[i][5]),omega_log=base[i][6]+alpha*(y[i][6]-base[i][6]);
   if(!std::isfinite(lk)||!std::isfinite(omega_log)||std::abs(lk-base[i][5])>.23||std::abs(omega_log-base[i][6])>.23){valid=false;continue;}
   double k=mesh.initial[i][5]*std::exp(10*lk),omega=mesh.initial[i][6]*std::exp(10*omega_log);
   try{auto v=gas.primitive(trial,original[i][4],k,omega);if(!std::isfinite(v[0])||v[4]<.8*original[i][4]||v[4]>1.2*original[i][4]||v[0]<=.01*freestream()[0]){valid=false;continue;}q[i]=trial;w[i]=v;}catch(const std::runtime_error&){valid=false;}
  }
  if(world.min(valid?1.:0.)<1)return false;
  if(options.cell_omega_wall){std::vector<double>target(mesh.n);for(auto&f:mesh.faces)if(f.type==3)target[f.l]=std::max(target[f.l],wall(f).omega);for(int i=0;i<mesh.n;++i)if(target[i]>0)w[i][6]=target[i];}
  return true;
 }
 std::array<double,2>merit(){
  double sum=0;for(int i=0;i<mesh.n;++i)for(int j=0;j<5;++j){double t=residual[i][j]/std::max(pseudo[i],1e-30);sum+=t*t;}double flow_merit=world.sum(sum);turbulence(5,false);turbulence(6,false);return {std::sqrt((flow_merit+scalar_merit[0]+scalar_merit[1])/mesh.global_cells),std::sqrt(flow_merit/mesh.global_cells)};
 }
 void accelerate(Anderson7&accel,const std::vector<W7>&before){
  double start=now();auto image=acceleration_state(),trial=accel.propose(before,image);auto basew=w;auto baseq=q;
  assemble(false);auto baseline=merit();bool accepted=false;
  for(double alpha:{1.,.5,.25}){w=basew;q=baseq;if(!set_acceleration_state(trial,image,basew,alpha))continue;gradients();properties();assemble(false);auto m=merit();if(m[0]<.99*baseline[0]&&m[1]<=1.05*baseline[1]){accepted=true;break;}}
  if(accepted)++aa_accepted;else{++aa_rejected;w=std::move(basew);q=std::move(baseq);gradients();properties();assemble(false);}
  acceleration_time+=now()-start;
 }
 void save_values(const std::string&directory,int iteration,const std::vector<W7>&values,const std::vector<double>*thermal_values=nullptr,const std::vector<std::array<double,3>>*scalar_values=nullptr)const{
  if(values.size()<static_cast<size_t>(mesh.n))throw std::runtime_error("M6 checkpoint has incomplete owned state");
  std::filesystem::create_directories(directory);std::ostringstream name;name<<"state-rank-"<<std::setfill('0')<<std::setw(3)<<world.rank<<".bin";
  std::ofstream f(std::filesystem::path(directory)/name.str(),std::ios::binary);f.write("M6STATE1",8);int32_t header[3]={mesh.n,iteration,world.size};f.write(reinterpret_cast<char*>(header),sizeof(header));f.write(reinterpret_cast<const char*>(mesh.global.data()),mesh.n*sizeof(int));f.write(reinterpret_cast<const char*>(values.data()),mesh.n*sizeof(W7));if(!f)throw std::runtime_error("M6 checkpoint write failed");
  if(options.enhanced_thermal_wall){if(!thermal_values||thermal_values->size()!=wall_temperature.size())throw std::runtime_error("incomplete M6 thermal wall cache snapshot");std::ostringstream thermal_name;thermal_name<<"wall-thermal-rank-"<<std::setfill('0')<<std::setw(3)<<world.rank<<".bin";std::ofstream cache(std::filesystem::path(directory)/thermal_name.str(),std::ios::binary);cache.write("M6WALL01",8);int32_t h[3]={static_cast<int32_t>(thermal_values->size()),iteration,world.size};cache.write(reinterpret_cast<char*>(h),sizeof(h));std::vector<int>ids(thermal_values->size());for(auto&face:mesh.faces)if(face.type==3)ids[wall_gradient_slot.at(face.id)]=face.id;cache.write(reinterpret_cast<char*>(ids.data()),ids.size()*sizeof(int));cache.write(reinterpret_cast<const char*>(thermal_values->data()),thermal_values->size()*sizeof(double));if(!cache)throw std::runtime_error("M6 thermal wall cache write failed");}
  if(options.normal_sst_transport){if(!scalar_values)throw std::runtime_error("normal SST checkpoint requires its boundary cache");sst_transport->save(directory,iteration,*scalar_values);}
 }
 void save(const std::string&directory,int iteration)const{auto scalar_snapshot=sst_transport?sst_transport->snapshot():std::vector<std::array<double,3>>{};save_values(directory,iteration,w,&wall_temperature,&scalar_snapshot);}
 void save_first_window(int iteration,int start,double seconds,double start_seconds,double ratio,const WindowSnapshot*first)const{
  auto directory=std::filesystem::path(options.output)/"first-fy-window";save(directory.string(),iteration);
  auto metadata=[&](const std::filesystem::path&path,int state_iteration,const D3&state_force,double state_seconds,const char*role){
   if(world.rank!=0)return;std::ofstream f(path/"window-acceptance.json");
   f<<std::setprecision(16)<<std::boolalpha<<"{\n\"case\":\"M6\",\n\"status\":\"first force-window observation; accuracy not accepted\",\n\"state_role\":\""<<role<<"\",\n\"cells\":"<<mesh.global_cells<<",\n\"ranks\":"<<world.size<<",\n\"iterations\":"<<state_iteration
    <<",\n\"restart_iteration_offset\":"<<iteration_offset<<",\n\"from_original_initial_field\":"<<(iteration_offset==0)
    <<",\n\"iteration_seconds_at_state\":"<<state_seconds<<",\n\"iteration_seconds_at_start\":"<<start_seconds<<",\n\"iteration_seconds_at_confirmation\":"<<seconds
    <<",\n\"fx\":"<<state_force[0]<<",\n\"fy\":"<<state_force[1]<<",\n\"fz\":"<<state_force[2]
    <<",\n\"first_100_step_fy_window_start\":"<<start<<",\n\"first_100_step_fy_window_confirmation\":"<<iteration<<",\n\"fy_window_ratio\":"<<ratio
    <<",\n\"window_start_field_captured\":"<<(first!=nullptr)<<",\n\"official_accuracy_tolerances_available\":false\n}\n";
   if(!f)throw std::runtime_error("first window metadata write failed");
  };
  metadata(directory,iteration,force,seconds,"confirmation");
  if(first){if(first->iteration!=start)throw std::runtime_error("window ring does not contain the required start iteration");auto path=directory/"start";save_values(path.string(),start,first->owned,&first->wall_temperature,&first->scalar_boundary);metadata(path,start,first->force,first->seconds,"window_start");}
 }
 void load(const std::string&directory){
  std::ostringstream name;name<<"state-rank-"<<std::setfill('0')<<std::setw(3)<<world.rank<<".bin";std::ifstream f(std::filesystem::path(directory)/name.str(),std::ios::binary);if(!f)throw std::runtime_error("M6 restart does not exist");
  char magic[8];Mesh::read(f,magic,8);int32_t header[3];Mesh::read(f,header,3);if(std::string(magic,8)!="M6STATE1"||header[0]!=mesh.n||header[2]!=world.size)throw std::runtime_error("M6 restart partition mismatch");
  iteration_offset=header[1];std::vector<int>ids(mesh.n);Mesh::read(f,ids.data(),mesh.n);if(!std::equal(ids.begin(),ids.end(),mesh.global.begin()))throw std::runtime_error("M6 restart ordering mismatch");Mesh::read(f,w.data(),mesh.n);for(int i=0;i<mesh.n;++i)q[i]=gas.conserved(w[i]);mesh.exchange(w);
  if(options.enhanced_thermal_wall){std::ostringstream cache_name;cache_name<<"wall-thermal-rank-"<<std::setfill('0')<<std::setw(3)<<world.rank<<".bin";std::ifstream cache(std::filesystem::path(directory)/cache_name.str(),std::ios::binary);if(!cache)throw std::runtime_error("enhanced thermal restart requires its saved wall cache; no reinitialization is inferred");char tag[8];int32_t h[3];Mesh::read(cache,tag,8);Mesh::read(cache,h,3);if(std::string(tag,8)!="M6WALL01"||h[0]!=wall_temperature.size()||h[1]!=iteration_offset||h[2]!=world.size)throw std::runtime_error("M6 thermal restart identity mismatch");std::vector<int>face_ids(wall_temperature.size());Mesh::read(cache,face_ids.data(),face_ids.size());for(auto&face:mesh.faces)if(face.type==3&&face_ids[wall_gradient_slot.at(face.id)]!=face.id)throw std::runtime_error("M6 thermal restart face ordering mismatch");Mesh::read(cache,wall_temperature.data(),wall_temperature.size());for(double t:wall_temperature)if(!(t>0&&std::isfinite(t)))throw std::runtime_error("invalid thermal restart cache value");}
  if(options.normal_sst_transport)sst_transport->load(directory,iteration_offset);
 }
 void wall_output()const{
  std::ostringstream name;name<<"wall-rank-"<<std::setfill('0')<<std::setw(3)<<world.rank<<".csv";std::ofstream f(std::filesystem::path(options.output)/name.str());f<<"face,cell,x,y,z,cell_pressure,reconstructed_pressure,k,omega,wall_omega,yplus\n"<<std::setprecision(16);
  for(auto&face:mesh.faces)if(face.type==3){auto a=wall(face);auto b=reconstruct(face.l,face);f<<face.id<<','<<mesh.global[face.l]<<','<<face.center[0]<<','<<face.center[1]<<','<<face.center[2]<<','<<w[face.l][0]*gas.rho0*gas.a0*gas.a0<<','<<b[0]*gas.rho0*gas.a0*gas.a0<<','<<w[face.l][5]*gas.a0*gas.a0<<','<<w[face.l][6]*gas.a0<<','<<a.omega*gas.a0<<','<<a.yplus<<'\n';}
 }
 void forces_only(){pressure_force={};viscous_force={};for(auto&face:mesh.faces)if(face.type==3)flux(face);pressure_force=world.sum(pressure_force);viscous_force=world.sum(viscous_force);for(int j=0;j<3;++j)force[j]=pressure_force[j]+viscous_force[j];}
 void forces_with_reference_cache(std::vector<std::array<double,1>>&cache){reference_velocity.swap(cache);forces_only();reference_velocity.swap(cache);}
 void finish_iteration_report(Anderson7&accelerator,const std::vector<W7>&before,bool use_aa,std::vector<std::array<double,1>>&force_reference_velocity){
  // Normal SST refreshes scalar reconstruction, while retaining the primary
  // convective RG and reference velocity used earlier in this flow update.
  // Native fast wall reports read these caches with the updated cell state.
  if(options.coupled_force_report_stage)forces_only();
  gradients();properties();
  if(use_aa)accelerate(accelerator,before);
  if(!options.coupled_force_report_stage){if(options.lagged_force_reference_velocity)forces_with_reference_cache(force_reference_velocity);else forces_only();}
 }
 void write_summary(int done,double elapsed,double scaled,int window_start,double window_start_seconds,double window_confirmation_seconds,double ratio)const{
  if(world.rank!=0)return;
  std::ofstream f(std::filesystem::path(options.output)/"summary.json");
  f<<std::setprecision(16)<<std::boolalpha
   <<"{\n\"case\":\"M6\",\n\"status\":\"development candidate; accuracy not yet accepted\""
   <<",\n\"ranks\":"<<world.size<<",\n\"cells\":"<<mesh.global_cells<<",\n\"iterations\":"<<done
   <<",\n\"order\":"<<options.order<<",\n\"cfl\":"<<options.cfl
   <<",\n\"scaled_flow_residual\":"<<scaled<<",\n\"iteration_seconds\":"<<elapsed
   <<",\n\"residual_normalization_basis\":\"maximum L1 in the first five updates of this run; reset on restart\""
   <<",\n\"gradient_seconds_rank0\":"<<gradient_time<<",\n\"assembly_seconds_rank0\":"<<assembly_time
   <<",\n\"linear_seconds_rank0\":"<<linear_time<<",\n\"turbulence_seconds_rank0\":"<<turbulence_time
   <<",\n\"stage_timers_overlap\":true"
   <<",\n\"fx\":"<<force[0]<<",\n\"fy\":"<<force[1]<<",\n\"fz\":"<<force[2]
   <<",\n\"first_100_step_fy_window_start\":"<<window_start
   <<",\n\"first_100_step_fy_window_confirmation\":"<<(window_start<0?-1:window_start+99)
   <<",\n\"first_100_step_fy_window_start_seconds\":"<<window_start_seconds
   <<",\n\"first_100_step_fy_window_confirmation_seconds\":"<<window_confirmation_seconds
   <<",\n\"fy_window_ratio\":"<<ratio<<",\n\"fy_window_history_first_iteration\":"<<iteration_offset+1
   <<",\n\"capture_window_start\":"<<options.capture_window_start<<",\n\"window_capture_peak_global_payload_bytes\":"<<window_capture_peak_payload_bytes<<",\n\"window_capture_seconds_rank0\":"<<window_capture_time
   <<",\n\"species_uniform_invariant\":true"
   <<",\n\"anderson_depth\":"<<options.anderson_depth<<",\n\"anderson_accepted\":"<<aa_accepted
   <<",\n\"anderson_rejected\":"<<aa_rejected<<",\n\"acceleration_seconds\":"<<acceleration_time
   <<",\n\"nonorthogonal_limit\":"<<options.nonorth_limit
   <<",\n\"cell_omega_profile\":"<<options.cell_omega_wall
   <<",\n\"implicit_viscous_gradients\":"<<options.implicit_gradients
   <<",\n\"viscous_gauss_gradients\":"<<options.viscous_gauss
   <<",\n\"wall_virtual_slip_gradients\":"<<options.wall_slip_gradients
   <<",\n\"coupled_wall_reconstruction\":"<<options.coupled_wall_reconstruction
   <<",\n\"symmetry_normal_lsq\":"<<options.symmetry_normal_lsq
   <<",\n\"enhanced_thermal_wall\":"<<options.enhanced_thermal_wall
   <<",\n\"symmetry_reconstructed_face_gradients\":"<<options.symmetry_face_gradients
   <<",\n\"coupled_farfield_reconstruction\":"<<options.coupled_farfield_reconstruction
   <<",\n\"symmetry_viscous_flux\":"<<options.symmetry_viscous_flux
   <<",\n\"farfield_zero_viscous_flux\":"<<options.farfield_zero_viscous_flux<<",\n\"normal_sst_transport\":"<<options.normal_sst_transport
   <<",\n\"lagged_sst_mass_flux\":"<<options.lagged_sst_mass_flux
   <<",\n\"normal_coupled_timestep\":"<<options.normal_coupled_timestep
   <<",\n\"case_material_polynomials\":"<<options.case_material_polynomials
   <<",\n\"flow_diffusion_cache\":"<<options.cache_flow_diffusion<<",\n\"flow_diffusion_cache_builds\":"<<flow_diffusion_cache_builds
   <<",\n\"flow_diffusion_cache_seconds_rank0\":"<<flow_diffusion_cache_seconds
   <<",\n\"flow_diffusion_cache_payload_bytes_rank0\":"<<(flow_diffusion_cell_cache.size()*sizeof(flow_diffusion_cell_cache[0])+flow_diffusion_face_cache.size()*sizeof(flow_diffusion_face_cache[0]))
   <<",\n\"normal_scalar_step\":"<<options.normal_scalar_step<<",\n\"scalar_negative_repairs\":["<<scalar_negative_repairs[0]<<","<<scalar_negative_repairs[1]<<"]"
   <<",\n\"skipped_discarded_postflow_assemblies\":"<<skipped_postflow_assemblies
#if defined(M6_SKIP_DISCARDED_POSTFLOW_ASSEMBLY) && M6_SKIP_DISCARDED_POSTFLOW_ASSEMBLY
   <<",\n\"compiled_skip_discarded_postflow_assembly\":true"
#else
   <<",\n\"compiled_skip_discarded_postflow_assembly\":false"
#endif
   <<",\n\"case_flow_update_limits\":"<<options.case_flow_update_limits
   <<",\n\"scalar_negative_source_rows\":["<<scalar_negative_source_rows[0]<<","<<scalar_negative_source_rows[1]<<"]"
#ifdef M6_AVX2_BUILD
   <<",\n\"compiled_for_avx2\":true"
#else
   <<",\n\"compiled_for_avx2\":false"
#endif
   <<",\n\"turbulence_gauss_gradients\":"<<options.turbulence_gauss
   <<",\n\"sst_vorticity_measure\":"<<options.sst_vorticity
   <<",\n\"harmonic_meanflow_diffusion\":"<<options.harmonic_diffusion
   <<",\n\"harmonic_scalar_diffusion\":"<<options.harmonic_scalars
   <<",\n\"omega_production_laminar_factor\":"<<options.omega_laminar_factor
   <<",\n\"omega_cross_diffusion_limit_to_destruction\":"<<options.cross_diffusion_limit
   <<",\n\"omega_linearize_negative_cross_diffusion\":"<<options.linearize_negative_cross
   <<",\n\"coupled_ausm\":"<<options.coupled_ausm
   <<",\n\"viscous_reference_velocity_floor\":"<<options.coupled_ausm
   <<",\n\"local_preconditioning\":"<<options.local_preconditioning
   <<",\n\"exact_local_flux_jacobian\":"<<options.exact_flux_jacobian
   <<",\n\"normal_approximate_interior_jacobian\":"<<options.normal_approximate_jacobian
   <<",\n\"normal_farfield_jacobian\":"<<options.normal_farfield_jacobian
   <<",\n\"normal_wall_jacobian\":"<<options.normal_wall_jacobian
   <<",\n\"normal_symmetry_jacobian\":"<<options.normal_symmetry_jacobian
   <<",\n\"turbulence_linear_tolerance\":"<<options.turbulence_linear_tolerance
   <<",\n\"turbulence_linear_max_steps\":"<<(options.turbulence_linear_max_steps?options.turbulence_linear_max_steps:std::max(4,options.linear_steps/2))
   <<",\n\"primitive_flow_update\":"<<options.primitive_flow_update
   <<",\n\"skip_unused_face_work\":"<<options.skip_unused_face_work
   <<",\n\"lagged_force_reference_velocity\":"<<options.lagged_force_reference_velocity
   <<",\n\"coupled_force_report_stage\":"<<options.coupled_force_report_stage
   <<",\n\"force_reference_velocity_stage\":\""<<((options.lagged_force_reference_velocity||options.coupled_force_report_stage)?"step-entry-before-meanflow-update":"current-state-after-turbulence")<<"\""
   <<",\n\"force_reconstruction_cache_stage\":\""<<(options.coupled_force_report_stage?"post-meanflow-before-final-sst-refresh":"post-sst-refreshed")<<"\""
   <<",\n\"steady_turbulence\":"<<options.steady_turbulence
   <<",\n\"turbulence_k_first\":"<<options.turbulence_k_first
   <<",\n\"characteristic_farfield_inviscid\":"<<options.characteristic_farfield
   <<",\n\"reconstruct_enthalpy\":"<<options.reconstruct_enthalpy
   <<",\n\"independent_face_thermodynamics\":"<<options.independent_face_thermo
   <<",\n\"boundary_independent_face_density\":"<<options.independent_face_thermo
   <<",\n\"lsq_boundary_face_ray\":"<<options.boundary_face_ray
   <<",\n\"batched_cgs2_orthogonalization\":"<<options.batched_orthogonalization
   <<",\n\"edge_correction\":"<<options.edge_correction<<",\n\"edge_implicit\":"<<options.edge_implicit
   <<",\n\"flow_relaxation\":"<<options.flow_relax<<",\n\"turbulence_relaxation\":"<<options.turb_relax
   <<",\n\"linear_steps\":"<<options.linear_steps<<",\n\"linear_tolerance\":"<<options.linear_tolerance
   <<",\n\"linear_audit_period\":"<<options.linear_audit_period<<",\n\"linear_audits\":"<<linear_audits<<",\n\"linear_audits_above_target\":"<<linear_target_missed<<",\n\"maximum_audited_true_relative_linear_residual\":"<<maximum_linear_true_residual
   <<",\n\"turbulence_linear_audit_enabled\":"<<options.audit_turbulence_linear<<",\n\"turbulence_linear_audits\":["<<turbulence_linear_audits[0]<<','<<turbulence_linear_audits[1]<<"],\n\"turbulence_linear_audits_above_target\":["<<turbulence_linear_misses[0]<<','<<turbulence_linear_misses[1]<<"],\n\"maximum_turbulence_true_linear_residual\":["<<turbulence_linear_maximum[0]<<','<<turbulence_linear_maximum[1]<<']'
   <<",\n\"recycle_flow_direction\":"<<recycle_flow_direction<<",\n\"recycle_attempts\":"<<recycled_attempts<<",\n\"recycle_accepted_starts\":"<<recycled_accepted<<",\n\"recycle_mean_initial_residual_ratio\":"<<(recycled_attempts?recycled_initial_ratio_sum/recycled_attempts:1.)
   <<",\n\"jacobian_period\":"<<options.jacobian_period
   <<",\n\"restart_iteration_offset\":"<<iteration_offset<<",\n\"updates_this_run\":"<<done-iteration_offset
   <<",\n\"vendor_libraries_linked\":false\n}\n";
 }
 void validate(){
  Field<1>x(mesh.nt);for(int i=0;i<mesh.n;++i)x[i][0]=static_cast<double>(mesh.global[i])+1;mesh.exchange(x);double halo_error=0;for(int i=0;i<mesh.nt;++i)halo_error=std::max(halo_error,std::abs(x[i][0]-(mesh.global[i]+1.)));
  double closure=0;std::vector<D3>sum(mesh.n);for(auto&f:mesh.faces)for(int j=0;j<3;++j){sum[f.l][j]+=f.normal[j]*f.area;if(f.r>=0&&f.r<mesh.n)sum[f.r][j]-=f.normal[j]*f.area;}
  for(auto&v:sum)for(double a:v)closure=std::max(closure,std::abs(a));halo_error=world.max(halo_error);closure=world.max(closure);
  if(halo_error!=0||closure>1e-10)throw std::runtime_error("M6 halo/geometry validation failed");
  if(world.rank==0)std::cout<<"M6_VALIDATION ranks="<<world.size<<" cells="<<mesh.global_cells<<" halo_error="<<halo_error<<" area_closure="<<closure<<std::endl;
 }
 void run(){
  std::filesystem::create_directories(options.output);std::ofstream history;if(world.rank==0){history.open(std::filesystem::path(options.output)/"history.csv");history<<"iteration,seconds,res_rho,res_mx,res_my,res_mz,res_energy,res_k,res_omega,scaled_flow,fx,fy,fz,pressure_fx,pressure_fy,pressure_fz,viscous_fx,viscous_fy,viscous_fz,gmres,limited_updates,fy_window_ratio,aa_accepted,aa_rejected\n"<<std::setprecision(16);}
  double start=now();std::array<double,5>base{};std::deque<double>fy,fy_seconds;int window_start=-1;double window_start_seconds=-1,window_confirmation_seconds=-1,scaled=1,ratio=1e100;int done=iteration_offset;Anderson7 accelerator(world,options.anderson_depth);
  std::vector<WindowSnapshot>window_snapshots(options.capture_window_start?100:0);
  std::ofstream linear_history;if(options.linear_audit_period>0&&world.rank==0){linear_history.open(std::filesystem::path(options.output)/"linear-audit.csv");linear_history<<"iteration,steps,estimated_relative_residual,true_relative_residual,target,recycled_start_used,recycled_start_relative_residual,recycle_coefficient,projected_trial_relative_residual\n"<<std::setprecision(16);}
  std::ofstream turbulence_linear_history;if(options.audit_turbulence_linear&&world.rank==0){turbulence_linear_history.open(std::filesystem::path(options.output)/"turbulence-linear-audit.csv");turbulence_linear_history<<"iteration,variable,steps,estimated_relative_residual,true_relative_residual,target\n"<<std::setprecision(16);}
  gradients();properties();std::vector<std::array<double,1>>force_reference_velocity;
  for(int local_it=1;local_it<=options.iterations;++local_it){int it=iteration_offset+local_it;bool use_aa=options.anderson_depth>0&&local_it>=options.anderson_start;auto before=use_aa?acceleration_state():std::vector<W7>{};bool matrix=(local_it-1)%options.jacobian_period==0;update_wall_temperature();assemble(matrix);
   if(options.lagged_force_reference_velocity)force_reference_velocity=reference_velocity;
   for(int j=0;j<5;++j){double t=0;for(auto&v:residual)t+=std::abs(v[j]);norms[j]=world.sum(t);if(local_it<=5)base[j]=std::max(base[j],norms[j]);}
   scaled=0;for(int j=0;j<5;++j)scaled=std::max(scaled,norms[j]/std::max(base[j],1e-30));
   bool audit=options.linear_audit_period>0&&local_it%options.linear_audit_period==0;int linear=meanflow(matrix,audit);if(audit&&world.rank==0){linear_history<<it<<','<<linear<<','<<linear_estimate<<','<<linear_true_residual<<','<<options.linear_tolerance<<','<<(linear_start.recycled?1:0)<<','<<linear_start.relative_residual<<','<<linear_start.coefficient<<','<<linear_start.trial_relative_residual<<'\n';linear_history.flush();}advance_turbulence(audit,true);if(audit&&options.audit_turbulence_linear&&world.rank==0){for(int variable=5;variable<=6;++variable){int slot=variable-5;turbulence_linear_history<<it<<','<<variable<<','<<turbulence_linear_steps[slot]<<','<<turbulence_linear_estimate[slot]<<','<<turbulence_linear_true[slot]<<','<<options.turbulence_linear_tolerance<<'\n';}turbulence_linear_history.flush();}done=it;finish_iteration_report(accelerator,before,use_aa,force_reference_velocity);
   // W7 is this mode's persisted primary state. Use the same map at every
   // step boundary as at restart, so a restart does not introduce a one-off
   // conservative/primitive rounding remap. No extra per-cell snapshot array.
   if(options.enhanced_thermal_wall||options.normal_sst_transport)for(int i=0;i<mesh.n;++i)q[i]=gas.conserved(w[i]);
   if(!window_snapshots.empty()){double capture_start=now();auto&entry=window_snapshots[(local_it-1)%100];entry.owned.assign(w.begin(),w.begin()+mesh.n);entry.wall_temperature=wall_temperature;if(sst_transport)entry.scalar_boundary=sst_transport->snapshot();entry.force=force;entry.iteration=it;window_capture_time+=now()-capture_start;window_capture_peak_payload_bytes=static_cast<uint64_t>(std::min(local_it,100))*(mesh.global_cells*sizeof(W7)+global_thermal_wall_faces*sizeof(double)+(sst_transport?sst_transport->global_faces*sizeof(std::array<double,3>):0));}
   double elapsed=world.max(now()-start),limited=world.sum(static_cast<double>(limited_updates));
   if(!window_snapshots.empty())window_snapshots[(local_it-1)%100].seconds=elapsed;
   fy.push_back(force[1]);fy_seconds.push_back(elapsed);if(fy.size()>100){fy.pop_front();fy_seconds.pop_front();}if(fy.size()==100){auto bounds=std::minmax_element(fy.begin(),fy.end());double avg=std::accumulate(fy.begin(),fy.end(),0.)/100;ratio=(*bounds.second-*bounds.first)/std::max(std::abs(avg),1e-30);if(ratio<=.01&&window_start<0){window_start=it-99;window_start_seconds=fy_seconds.front();window_confirmation_seconds=elapsed;const WindowSnapshot*first=window_snapshots.empty()?nullptr:&window_snapshots[(local_it-100)%100];save_first_window(it,window_start,elapsed,window_start_seconds,ratio,first);std::vector<WindowSnapshot>().swap(window_snapshots);}}
   if(world.rank==0){history<<it<<','<<elapsed;for(double v:norms)history<<','<<v;history<<','<<scaled;for(double v:force)history<<','<<v;for(double v:pressure_force)history<<','<<v;for(double v:viscous_force)history<<','<<v;history<<','<<linear<<','<<limited<<','<<ratio<<','<<aa_accepted<<','<<aa_rejected<<'\n';history.flush();std::cout<<"M6_ITER "<<it<<" scaled="<<scaled<<" Fx="<<force[0]<<" Fy="<<force[1]<<" Fz="<<force[2]<<" time="<<elapsed<<" limited="<<limited<<std::endl;}
   if(options.checkpoint>0&&it%options.checkpoint==0)save(options.output+"/checkpoint",it);
  }
  double elapsed=world.max(now()-start);save(options.output,done);wall_output();
  write_summary(done,elapsed,scaled,window_start,window_start_seconds,window_confirmation_seconds,ratio);
 }
};

void self_test(){
 std::array<double,48>meta{};meta[0]=1.1361339060496805;meta[1]=352.517926811456;meta[2]=298.15;meta[3]=299.124701264874;meta[4]=101325;meta[5]=1.72e-5;meta[6]=.0454;meta[7]=1059.4853515625;
 std::array<double,5>low={985.1700244,.4572524926,-.0010998299666,1.5381987004e-6,-6.675425469e-10},high={901.5941633,.4449195144,-.00015867025232,2.7184049518e-8,-1.7782013004e-12};for(int i=0;i<5;++i){meta[8+i]=low[i];meta[13+i]=high[i];}
 Gas gas(meta);double roundtrip=0,fluxerror=0,jacerror=0;D3 n={.3,-.4,std::sqrt(.75)};
 double primitive_increment_error=0;
 for(double temperature:{250.,500.,1200.,4800.}){
  W7 state={.71,.5,.04,-.08,temperature/meta[2],1e-4,1.};auto q=gas.conserved(state);Q5 direction={.013,-.021,.017,.008,.031};auto analytic=linearized_primitive_increment(gas,state,direction);double eps=1e-5;auto plus=q,minus=q;for(int j=0;j<5;++j){plus[j]+=eps*direction[j];minus[j]-=eps*direction[j];}auto wp=gas.primitive(plus,state[4],state[5],state[6]),wm=gas.primitive(minus,state[4],state[5],state[6]);
  for(int j=0;j<5;++j)primitive_increment_error=std::max(primitive_increment_error,std::abs((wp[j]-wm[j])/(2*eps)-analytic[j])/std::max(1e-5,std::abs(analytic[j])));
  auto trial=primitive_flow_trial(state,analytic,.75);if(trial[5]!=state[5]||trial[6]!=state[6])throw std::runtime_error("primitive meanflow trial changed scalar fields");
 }if(primitive_increment_error>1e-6)throw std::runtime_error("primitive change-of-variables directional check failed");std::cout<<"M6_PRIMITIVE_INCREMENT_TEST relative_error="<<primitive_increment_error<<std::endl;
 // Every homogeneous far-field branch must preserve a uniform physical flux.
 // This also exercises supersonic branches not present in the M6 observations.
 for(double mach:{.6,1.4})for(double temperature:{298.15,500.})for(D3 normal:std::array<D3,5>{{{1,0,0},{-1,0,0},{0,1,0},{.01,0,std::sqrt(.9999)},{.3,-.4,std::sqrt(.75)}}}){
  W7 state={meta[4]/(meta[0]*meta[1]*meta[1]),0,0,0,temperature/meta[2],1e-4,1};state[1]=mach*gas.sound(state);
  auto observed=characteristic_farfield(gas,state,state,state,normal,gas.rho(state),std::pow(gas.sound(state),2),gas.rho(state),gas.enthalpy(state));auto expected=gas.physical(state,normal);
  for(int j=0;j<5;++j)if(std::abs(observed.flux[j]-expected[j])>1e-11*std::max(1.,std::abs(expected[j])))throw std::runtime_error("far-field uniform-flow invariant failed");
 }
 for(double T:{250.,298.15,500.,1200.,4800.}){W7 w={.71,.5,.04,-.08,T/meta[2],1e-4,1};auto q=gas.conserved(w);auto recovered=gas.primitive(q,.99*w[4],w[5],w[6]);for(int j=0;j<5;++j)roundtrip=std::max(roundtrip,std::abs(recovered[j]-w[j]));auto a=gas.physical(w,n),b=gas.ausm(w,w,n);for(int j=0;j<5;++j)fluxerror=std::max(fluxerror,std::abs(a[j]-b[j]));
  auto jac=gas.jacobian(w,n);for(int j=0;j<5;++j){double eps=1e-6*std::max(std::abs(q[j]),.01);auto p=q,m=q;p[j]+=eps;m[j]-=eps;auto fp=gas.physical(gas.primitive(p,w[4],w[5],w[6]),n),fm=gas.physical(gas.primitive(m,w[4],w[5],w[6]),n);for(int k=0;k<5;++k)jacerror=std::max(jacerror,std::abs((fp[k]-fm[k])/(2*eps)-jac[5*k+j])/std::max(1.,std::abs(jac[5*k+j])));}
 }
 for(double beta:{1.,.2,.01}){W7 v={.71,.5,.04,-.08,1.,1e-4,1.};double a2=std::pow(gas.sound(v),2);auto a=gas.physical(v,n),b=gas.coupled_ausm(v,v,n,a2,a2,.5*beta*a2,.5*beta*a2);for(int j=0;j<5;++j)fluxerror=std::max(fluxerror,std::abs(a[j]-b[j]));}
 W7 left={.72,.6,.03,.02,1,1e-4,1},right={.68,.55,-.01,.03,.98,1e-4,1};auto f=gas.ausm(left,right,n);auto cp=gas.coupled_ausm(left,right,n,1.1,.9,.1,.2);D3 reverse={-n[0],-n[1],-n[2]};auto cm=gas.coupled_ausm(right,left,reverse,.9,1.1,.2,.1);for(int j=0;j<5;++j)fluxerror=std::max(fluxerror,std::abs(cp[j]+cm[j]));for(double&t:n)t=-t;auto reverse_flux=gas.ausm(right,left,n);for(int j=0;j<5;++j)fluxerror=std::max(fluxerror,std::abs(f[j]+reverse_flux[j]));
 {W7 cell={.71,.5,.04,-.08,1.,1e-4,1.};double a2=std::pow(gas.sound(cell),2);auto dp=gas.pressure_jacobian(cell);Q5 v={1,cell[1],cell[2],cell[3],gas.enthalpy(cell)};double dot=0;for(int j=0;j<5;++j)dot+=dp[j]*v[j];if(std::abs(dot/a2-1)>1e-12)throw std::runtime_error("preconditioning acoustic eigenvalue test failed");for(double beta:{1.,.1,.01}){auto p=gas.preconditioning_matrix(cell,beta*a2);Q5 product{};for(int i=0;i<5;++i)for(int j=0;j<5;++j)product[i]+=p[5*i+j]*v[j];for(int i=0;i<5;++i)if(std::abs(product[i]-v[i]/beta)>1e-9)throw std::runtime_error("preconditioning rank-one test failed");}}
 double flux_jacobian_error=0;
 for(double speed:{-2.3,-.4,.4,2.3})for(double beta:{1.,.2,.002}){
  W7 l={.72,speed,.06,-.08,1.2,1e-4,1.},r={.68,.9*speed,-.01,.03,1.1,1e-4,1.};auto fl=l,fr=r;fl[0]+=.002;fl[1]+=.013;fl[4]+=.011;fr[0]-=.001;fr[2]+=.007;fr[4]-=.012;
  D3 normal={.8,-.36,.48};double a2l=std::pow(gas.sound(l),2),a2r=std::pow(gas.sound(r),2),ref=.25*beta*(a2l+a2r);auto jac=coupled_flux_jacobian(gas,l,r,fl,fr,normal,a2l,a2r,ref,ref);
  for(int side=0;side<2;++side){auto cell=side==0?l:r,face=side==0?fl:fr;auto qcell=gas.conserved(cell);for(int col=0;col<5;++col){double eps=1e-6*std::max(std::abs(qcell[col]),.01);auto qp=qcell,qm=qcell;qp[col]+=eps;qm[col]-=eps;auto wp=gas.primitive(qp,cell[4],cell[5],cell[6]),wm=gas.primitive(qm,cell[4],cell[5],cell[6]);for(int k=0;k<5;++k){wp[k]+=face[k]-cell[k];wm[k]+=face[k]-cell[k];}auto fp=side==0?gas.coupled_ausm(wp,fr,normal,a2l,a2r,ref,ref):gas.coupled_ausm(fl,wp,normal,a2l,a2r,ref,ref);auto fm=side==0?gas.coupled_ausm(wm,fr,normal,a2l,a2r,ref,ref):gas.coupled_ausm(fl,wm,normal,a2l,a2r,ref,ref);for(int row=0;row<5;++row){double exact=(side==0?jac.left:jac.right)[5*row+col];flux_jacobian_error=std::max(flux_jacobian_error,std::abs((fp[row]-fm[row])/(2*eps)-exact)/std::max(1.,std::abs(exact)));}}}
 }
 // Hold independent density/H reconstruction increments fixed while varying
 // conservative cell states; this is distinct from the EOS face-state chain.
 for(double speed:{-2.3,-.4,.4,2.3})for(double beta:{1.,.2,.002}){
  W7 l={.72,speed,.06,-.08,1.2,1e-4,1.},r={.68,.9*speed,-.01,.03,1.1,1e-4,1.};auto fl=l,fr=r;fl[0]+=.002;fl[1]+=.013;fl[4]+=.011;fr[0]-=.001;fr[2]+=.007;fr[4]-=.012;
  D3 independent_normal={.8,-.36,.48};double a2l=std::pow(gas.sound(l),2),a2r=std::pow(gas.sound(r),2),ref=.25*beta*(a2l+a2r);std::array<double,4>thermo={gas.rho(l)+.013,gas.rho(r)-.007,gas.enthalpy(l)-.017,gas.enthalpy(r)+.023};auto jac=coupled_flux_jacobian(gas,l,r,fl,fr,independent_normal,a2l,a2r,ref,ref,&thermo);
  for(int side=0;side<2;++side){auto cell=side==0?l:r,face=side==0?fl:fr;auto qc=gas.conserved(cell);
   for(int col=0;col<5;++col){double eps=1e-6*std::max(std::abs(qc[col]),.01);auto qp=qc,qm=qc;qp[col]+=eps;qm[col]-=eps;auto wp=gas.primitive(qp,cell[4],cell[5],cell[6]),wm=gas.primitive(qm,cell[4],cell[5],cell[6]);auto tp=thermo,tm=thermo;
    tp[side]+=gas.rho(wp)-gas.rho(cell);tm[side]+=gas.rho(wm)-gas.rho(cell);tp[side+2]+=gas.enthalpy(wp)-gas.enthalpy(cell);tm[side+2]+=gas.enthalpy(wm)-gas.enthalpy(cell);
    for(int k=0;k<5;++k){wp[k]+=face[k]-cell[k];wm[k]+=face[k]-cell[k];}
    auto fp=side==0?gas.coupled_ausm(wp,fr,independent_normal,a2l,a2r,ref,ref,&tp):gas.coupled_ausm(fl,wp,independent_normal,a2l,a2r,ref,ref,&tp);auto fm=side==0?gas.coupled_ausm(wm,fr,independent_normal,a2l,a2r,ref,ref,&tm):gas.coupled_ausm(fl,wm,independent_normal,a2l,a2r,ref,ref,&tm);
    for(int row=0;row<5;++row){double exact=(side==0?jac.left:jac.right)[5*row+col];flux_jacobian_error=std::max(flux_jacobian_error,std::abs((fp[row]-fm[row])/(2*eps)-exact)/std::max(1.,std::abs(exact)));}
   }
  }
  auto forward=gas.coupled_ausm(fl,fr,independent_normal,a2l,a2r,ref,ref,&thermo);D3 independent_reverse={-independent_normal[0],-independent_normal[1],-independent_normal[2]};std::array<double,4>swapped={thermo[1],thermo[0],thermo[3],thermo[2]};auto backward=gas.coupled_ausm(fr,fl,independent_reverse,a2r,a2l,ref,ref,&swapped);for(int j=0;j<5;++j)fluxerror=std::max(fluxerror,std::abs(forward[j]+backward[j]));
 }
 if(flux_jacobian_error>1e-6)throw std::runtime_error("coupled flux automatic-differentiation check failed");
 std::cout<<"M6_AUSM_JACOBIAN_TEST max_scaled_error="<<flux_jacobian_error<<std::endl;
 double wall_jacobian_error=0;
 for(bool independent:{false,true})for(double normal_velocity:{-.1,.03,.4}){W7 cell={.71,.5,normal_velocity,-.08,1.,1e-4,1.},face=cell;face[0]*=1.02;face[1]+=.02;face[2]-=.01;face[4]*=.98;D3 normal={.3,-.4,std::sqrt(.75)};double reference2=.07,density=independent?1.03*gas.rho(cell):gas.rho(face);auto exact=coupled_wall_pressure_jacobian(gas,cell,face,normal,reference2,density,independent);auto conserved=gas.conserved(cell);
  for(int j=0;j<5;++j){double epsilon=1e-6*std::max(std::abs(conserved[j]),.01);auto qp=conserved,qm=conserved;qp[j]+=epsilon;qm[j]-=epsilon;auto plus=gas.primitive(qp,cell[4],cell[5],cell[6]),minus=gas.primitive(qm,cell[4],cell[5],cell[6]);double rp=density+gas.rho(plus)-gas.rho(cell),rm=density+gas.rho(minus)-gas.rho(cell);for(int k=0;k<5;++k){plus[k]+=face[k]-cell[k];minus[k]+=face[k]-cell[k];}if(!independent){rp=gas.rho(plus);rm=gas.rho(minus);}double fd=(coupled_wall_pressure(plus,normal,reference2,rp)-coupled_wall_pressure(minus,normal,reference2,rm))/(2*epsilon);wall_jacobian_error=std::max(wall_jacobian_error,std::abs(fd-exact[j])/std::max(1.,std::abs(exact[j])));}
 }if(wall_jacobian_error>1e-6)throw std::runtime_error("wall pressure frozen-reconstruction Jacobian check failed");std::cout<<"M6_WALL_PRESSURE_JACOBIAN_TEST max_scaled_error="<<wall_jacobian_error<<std::endl;
 D3 manufactured={2.,3.,4.},connector={1.,1e-7,0.},normal={0.,1.,0.};for(bool edge:{false,true}){auto g=correct_face_gradient(manufactured,connector,normal,dot3(manufactured,connector),edge,-1);for(int j=0;j<3;++j)if(std::abs(g[j]-manufactured[j])>1e-12)throw std::runtime_error("gradient linear-exactness test failed");}auto bounded=correct_face_gradient(manufactured,connector,normal,dot3(manufactured,connector)+.2,true,-1);if(length3(bounded)>6.)throw std::runtime_error("tangent-face gradient amplification test failed");
 auto wall=wall_model(1.1,1.72e-5,1e-8,1e6,.0001,1e-6);double lam=2*1.72e-5/(1.1*.072*1e-12),werror=std::abs(wall.omega/lam-1);
 if(roundtrip>1e-10||fluxerror>1e-12||jacerror>1e-6||werror>1e-5)throw std::runtime_error("M6 thermodynamic/flux/Jacobian/wall self-test failed");
 std::cout<<"M6_CORE_TEST roundtrip="<<roundtrip<<" flux="<<fluxerror<<" jacobian="<<jacerror<<" viscous_wall_limit="<<werror<<std::endl;
}
void acceleration_self_test(const MPIWorld&world){Anderson7 aa(world,3);std::vector<W7>y(1);for(int it=0;it<60;++it){auto image=y;for(int k=0;k<7;++k){double a=.3+.08*k;image[0][k]=a*y[0][k]+(1-a)*(k+1);}y=aa.propose(y,image);}double error=0;for(int k=0;k<7;++k)error=std::max(error,std::abs(y[0][k]-(k+1)));if(error>1e-9)throw std::runtime_error("M6 Anderson linear fixed-point verification failed");if(world.rank==0)std::cout<<"M6_ANDERSON_TEST error="<<error<<std::endl;}
}
#ifndef M6_NO_MAIN
int main(int argc,char**argv){
 m6::MPIWorld world(argc,argv);try{
  m6::Options o;bool self=false;
  for(int i=1;i<argc;++i){std::string a=argv[i];auto value=[&](){if(i+1>=argc)throw std::runtime_error("missing argument for "+a);return std::string(argv[++i]);};
   if(a=="--linear-audit-period"){o.linear_audit_period=std::stoi(value());if(o.linear_audit_period<0)throw std::runtime_error("linear audit period must be nonnegative");continue;}
   if(a=="--recycle-linear-direction"){o.recycle_linear_direction=true;continue;}
   if(a=="--cache-flow-diffusion"){o.cache_flow_diffusion=true;continue;}
   if(a=="--case-material-polynomials"){o.case_material_polynomials=true;continue;}
   if(a=="--normal-approximate-jacobian"){o.normal_approximate_jacobian=true;continue;}
   if(a=="--normal-farfield-jacobian"){o.normal_farfield_jacobian=true;continue;}
   if(a=="--normal-wall-jacobian"){o.normal_wall_jacobian=true;continue;}
   if(a=="--normal-symmetry-jacobian"){o.normal_symmetry_jacobian=true;continue;}
   if(a=="--turbulence-linear-tolerance"){o.turbulence_linear_tolerance=std::stod(value());continue;}
   if(a=="--turbulence-linear-steps"){o.turbulence_linear_max_steps=std::stoi(value());continue;}
   if(a=="--primitive-flow-update"){o.primitive_flow_update=true;continue;}
   if(a=="--case-flow-update-limits"){o.case_flow_update_limits=true;continue;}
   if(a=="--skip-unused-face-work"){o.skip_unused_face_work=true;continue;}
    if(a=="--lagged-force-reference-velocity"){o.lagged_force_reference_velocity=true;continue;}
    if(a=="--coupled-force-report-stage"){o.coupled_force_report_stage=true;continue;}
   if(a=="--export-linear-operator"){o.export_linear_operator=true;continue;}
   if(a=="--turbulence-k-first"){o.turbulence_k_first=true;continue;}
   if(a=="--characteristic-farfield"){o.characteristic_farfield=true;continue;}
   if(a=="--capture-window-start"){o.capture_window_start=true;continue;}
   if(a=="--wall-slip-gradients"){o.wall_slip_gradients=true;continue;}
   if(a=="--coupled-wall-reconstruction"){o.coupled_wall_reconstruction=true;continue;}
   if(a=="--symmetry-normal-lsq"){o.symmetry_normal_lsq=true;continue;}
   if(a=="--enhanced-thermal-wall"){o.enhanced_thermal_wall=true;continue;}
   if(a=="--symmetry-face-gradients"){o.symmetry_face_gradients=true;continue;}
   if(a=="--audit-turbulence-linear"){o.audit_turbulence_linear=true;continue;}
   if(a=="--coupled-farfield-reconstruction"){o.coupled_farfield_reconstruction=true;continue;}
   if(a=="--symmetry-viscous-flux"){o.symmetry_viscous_flux=true;continue;}
   if(a=="--farfield-zero-viscous-flux"){o.farfield_zero_viscous_flux=true;continue;}
   if(a=="--normal-sst-transport"){o.normal_sst_transport=true;continue;}
   if(a=="--lagged-sst-mass-flux"){o.lagged_sst_mass_flux=true;continue;}
   if(a=="--normal-coupled-timestep"){o.normal_coupled_timestep=true;continue;}
   if(a=="--normal-scalar-step"){o.normal_scalar_step=true;continue;}
   if(a=="--diagnostic-directional-step"){o.diagnostic_directional_step=std::stod(value());if(!(o.diagnostic_directional_step>0&&o.diagnostic_directional_step<=.01))throw std::runtime_error("diagnostic step must be positive and at most .01");continue;}
   if(a=="--input")o.input=value();else if(a=="--output")o.output=value();else if(a=="--restart")o.restart=value();else if(a=="--rank-cpus")o.cpu_list=value();else if(a=="--limit-nonorth")o.nonorth_limit=std::stod(value());else if(a=="--omega-laminar-factor")o.omega_laminar_factor=std::stod(value());else if(a=="--cross-diffusion-limit")o.cross_diffusion_limit=std::stod(value());else if(a=="--explicit-negative-cross-diffusion")o.linearize_negative_cross=false;else if(a=="--anderson-depth")o.anderson_depth=std::stoi(value());else if(a=="--anderson-start")o.anderson_start=std::stoi(value());else if(a=="--iterations")o.iterations=std::stoi(value());else if(a=="--order")o.order=std::stoi(value());else if(a=="--linear-steps")o.linear_steps=std::stoi(value());else if(a=="--linear-tolerance")o.linear_tolerance=std::stod(value());else if(a=="--jacobian-period")o.jacobian_period=std::stoi(value());else if(a=="--checkpoint")o.checkpoint=std::stoi(value());else if(a=="--flow-relax")o.flow_relax=std::stod(value());else if(a=="--turb-relax")o.turb_relax=std::stod(value());else if(a=="--validate")o.validate=true;else if(a=="--check-viscous-jacobian")o.check_viscous=true;else if(a=="--implicit-viscous-gradients")o.implicit_gradients=true;else if(a=="--cell-omega-wall")o.cell_omega_wall=true;else if(a=="--steady-turbulence")o.steady_turbulence=true;else if(a=="--sst-vorticity")o.sst_vorticity=true;else if(a=="--lsq-viscous-gradients")o.viscous_gauss=false;else if(a=="--reconstruct-enthalpy")o.reconstruct_enthalpy=true;else if(a=="--independent-face-thermodynamics")o.independent_face_thermo=true;else if(a=="--boundary-face-ray")o.boundary_face_ray=true;else if(a=="--batched-orthogonalization")o.batched_orthogonalization=true;else if(a=="--check-linear-solver")o.check_linear=true;else if(a=="--arithmetic-diffusion")o.harmonic_diffusion=false;else if(a=="--harmonic-scalar-diffusion")o.harmonic_scalars=true;else if(a=="--edge-correction")o.edge_correction=true;else if(a=="--edge-implicit")o.edge_implicit=true;else if(a=="--coupled-ausm")o.coupled_ausm=true;else if(a=="--gauss-turbulence")o.turbulence_gauss=true;else if(a=="--local-preconditioning")o.local_preconditioning=true;else if(a=="--exact-flux-jacobian")o.exact_flux_jacobian=true;else if(a=="--face-omega-wall")o.cell_omega_wall=false;else if(a=="--standard-wall-production")o.wall_production=false;else if(a=="--tke-pressure")throw std::runtime_error("unvalidated isotropic-stress experiment is disabled");else if(a=="--self-test")self=true;else if(a=="--help"){if(world.rank==0)std::cout<<"m6 --input PARTITION_DIRECTORY --output NEW_RUN --iterations N [--validate] [--order 1|2] [--linear-steps N] [--jacobian-period N] [--restart CHECKPOINT_DIRECTORY]\n";return 0;}else throw std::runtime_error("unknown argument "+a);
  }
  if(self){m6::acceleration_self_test(world);if(world.rank==0)m6::self_test();return 0;}if(o.input.empty())throw std::runtime_error("--input is required");if(!std::isfinite(o.turbulence_linear_tolerance)||o.turbulence_linear_tolerance<=0||o.turbulence_linear_tolerance>.1||o.turbulence_linear_max_steps<0||(o.turbulence_linear_max_steps>0&&o.turbulence_linear_max_steps<4)||o.turbulence_linear_max_steps>256)throw std::runtime_error("invalid turbulence linear controls; use tolerance in (0,.1] and steps zero/default or 4..256");if(o.iterations<0||o.linear_steps<1||o.jacobian_period<1||(o.order!=1&&o.order!=2)||o.anderson_depth<0||o.anderson_depth>8)throw std::runtime_error("invalid numerical controls");
  m6::apply_cpu_affinity(world,o.cpu_list);
  if(!std::isfinite(o.cross_diffusion_limit)||o.cross_diffusion_limit<0)throw std::runtime_error("cross diffusion limit must be finite and nonnegative");
  if(!std::isfinite(o.omega_laminar_factor)||o.omega_laminar_factor<0)throw std::runtime_error("omega laminar factor must be finite and nonnegative");
  if(o.check_linear&&o.check_viscous)throw std::runtime_error("choose either the linear solver or diffusion directional validation");
  if(o.export_linear_operator&&(o.check_linear||o.check_viscous||o.validate))throw std::runtime_error("choose one validation workflow");
  if(o.recycle_linear_direction&&(o.check_linear||o.check_viscous||o.export_linear_operator||o.validate))throw std::runtime_error("recycled directions require an integration run; use --linear-audit-period to verify the actual operator");
  if(o.local_preconditioning&&!o.coupled_ausm)throw std::runtime_error("local preconditioning requires --coupled-ausm");
  if(o.wall_slip_gradients&&!o.viscous_gauss)throw std::runtime_error("wall virtual velocity currently requires physical Gauss gradients");
  if(o.audit_turbulence_linear&&o.linear_audit_period<=0)throw std::runtime_error("turbulence linear audit requires a positive linear-audit-period");
  if(o.independent_face_thermo&&(!o.coupled_ausm||o.reconstruct_enthalpy))throw std::runtime_error("independent density/enthalpy reconstruction requires coupled AUSM and separate primitive temperature");
  if(o.exact_flux_jacobian&&(!o.coupled_ausm||o.reconstruct_enthalpy))throw std::runtime_error("exact flux Jacobian requires --coupled-ausm and primitive reconstruction");
  if((o.implicit_gradients||o.check_viscous)&&(o.nonorth_limit>=0||o.edge_correction||o.reconstruct_enthalpy))throw std::runtime_error("full-gradient Jacobian requires unlimited normal correction and primitive reconstruction");
  m6::Solver solver(world,o);solver.validate();if(o.export_linear_operator)solver.export_linear_operator_observation();else if(o.check_linear)solver.check_linear_solver();else if(o.check_viscous)solver.check_viscous_operator();else if(!o.validate)solver.run();return 0;
 }catch(const std::exception&e){std::cerr<<"M6_ERROR rank="<<world.rank<<" "<<e.what()<<std::endl;MPI_Abort(MPI_COMM_WORLD,2);return 2;}
}
#endif
