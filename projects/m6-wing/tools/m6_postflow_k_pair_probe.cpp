// Direct comparison at the actual post-flow caller boundary, no time advance.
// Two Solvers share the same input and implementation. Only gradients' explicit
// caller context differs; both retain the subsequent real farfield/k refresh.
#define M6_NO_MAIN
#include "../src/m6.cpp"
#include <cstring>
#include "m6_cartesian_fixture.hpp"

namespace {
template<class T>uint64_t same(const std::vector<T>&a,const std::vector<T>&b,const char*name){
 if(a.size()!=b.size())throw std::runtime_error(std::string(name)+" lengths differ");
 auto bytes=a.size()*sizeof(T);
 if(bytes&&std::memcmp(a.data(),b.data(),bytes))throw std::runtime_error(std::string(name)+" live bytes differ");
 return bytes;
}
uint64_t common(const m6::Solver&a,const m6::Solver&b){
 uint64_t bytes=0;
#define SAME(member) bytes+=same(a.member,b.member,#member)
 SAME(w);SAME(q);SAME(hi);SAME(physical_gradient);SAME(thermodynamics);SAME(thermo_gradient);SAME(thermo_limiter);
 SAME(total_enthalpy);SAME(enthalpy_gradient);SAME(wall_virtual_jacobian);SAME(wall_temperature);SAME(wall_temperature_jacobian);
 for(int slot=0;slot<2;++slot){bytes+=same(a.sst_transport->boundary[slot],b.sst_transport->boundary[slot],"boundary");
  bytes+=same(a.sst_transport->allocated[slot],b.sst_transport->allocated[slot],"allocated");}
 return bytes;
}
uint64_t full(const m6::Solver&a,const m6::Solver&b){
 uint64_t bytes=common(a,b);
 SAME(lo);SAME(grad);SAME(scalar_gradient);SAME(prop);SAME(pseudo);SAME(massflux);SAME(reference_velocity);
 SAME(wall_weight);SAME(wall_pk);SAME(residual);SAME(flow.a);SAME(scalar.a);
#undef SAME
 for(int slot=0;slot<2;++slot)bytes+=same(a.sst_transport->reference[slot],b.sst_transport->reference[slot],"reference");
 if(std::memcmp(a.force.data(),b.force.data(),sizeof(a.force))||std::memcmp(a.pressure_force.data(),b.pressure_force.data(),sizeof(a.pressure_force))||std::memcmp(a.viscous_force.data(),b.viscous_force.data(),sizeof(a.viscous_force)))throw std::runtime_error("stage force bytes differ");
 return bytes;
}
std::pair<uint64_t,uint64_t>partial(const m6::Solver&a,const m6::Solver&b){
 uint64_t bytes=common(a,b),different_k_cells=0;
 bytes+=same(a.sst_transport->reference[1],b.sst_transport->reference[1],"reference_omega");
 if(a.lo.size()!=b.lo.size()||a.grad.size()!=b.grad.size()||a.scalar_gradient.size()!=b.scalar_gradient.size())throw std::runtime_error("partial sizes differ");
 for(size_t i=0;i<a.lo.size();++i){
  for(int k=0;k<7;++k)if(k!=5){if(std::memcmp(&a.lo[i][k],&b.lo[i][k],sizeof(double)))throw std::runtime_error("live partial limiter differs");bytes+=sizeof(double);}
  for(int k=0;k<21;++k)if(k<15||k>=18){if(std::memcmp(&a.grad[i][k],&b.grad[i][k],sizeof(double)))throw std::runtime_error("live partial gradient differs");bytes+=sizeof(double);}
  if(std::memcmp(a.scalar_gradient[i].data()+3,b.scalar_gradient[i].data()+3,3*sizeof(double)))throw std::runtime_error("omega physical gradient differs");
  bytes+=3*sizeof(double);
  if(std::memcmp(&a.lo[i][5],&b.lo[i][5],sizeof(double))||std::memcmp(a.grad[i].data()+15,b.grad[i].data()+15,3*sizeof(double))||
     std::memcmp(a.scalar_gradient[i].data(),b.scalar_gradient[i].data(),3*sizeof(double)))++different_k_cells;
 }
 return {bytes,different_k_cells};
}
m6::Options options(const std::string&input,const std::string&restart,bool synthetic){
 m6::Options o;o.input=input;o.restart=restart;o.iterations=0;o.checkpoint=0;
 o.coupled_ausm=true;o.turbulence_gauss=true;o.harmonic_scalars=true;o.independent_face_thermo=true;o.boundary_face_ray=true;
 o.exact_flux_jacobian=true;o.linear_steps=24;o.linear_tolerance=.1;o.linear_audit_period=1;o.flow_relax=.75;o.turb_relax=.8;
 o.batched_orthogonalization=true;o.steady_turbulence=true;o.turbulence_k_first=true;o.characteristic_farfield=true;
 o.wall_slip_gradients=true;o.coupled_wall_reconstruction=true;o.symmetry_normal_lsq=true;o.enhanced_thermal_wall=true;
 o.symmetry_face_gradients=true;o.audit_turbulence_linear=true;o.coupled_farfield_reconstruction=true;
 o.symmetry_viscous_flux=true;o.farfield_zero_viscous_flux=true;o.normal_sst_transport=true;o.lagged_sst_mass_flux=true;
 o.normal_coupled_timestep=true;o.normal_scalar_step=true;o.local_preconditioning=true;o.case_material_polynomials=!synthetic;
 o.normal_approximate_jacobian=true;o.normal_farfield_jacobian=true;o.primitive_flow_update=true;
 o.skip_unused_face_work=true;o.case_flow_update_limits=!synthetic;o.normal_wall_jacobian=true;o.normal_symmetry_jacobian=true;
 return o;
}
}

int main(int argc,char**argv){m6::MPIWorld world(argc,argv);try{
 std::string input,output,restart,cpus;bool synthetic=false;
 for(int i=1;i<argc;++i){std::string arg=argv[i];auto value=[&](){if(++i>=argc)throw std::runtime_error("missing argument");return std::string(argv[i]);};
  if(arg=="--input")input=value();else if(arg=="--output")output=value();else if(arg=="--restart")restart=value();
  else if(arg=="--rank-cpus")cpus=value();else if(arg=="--synthetic"){synthetic=true;output=value();}
  else throw std::runtime_error("unknown post-flow diagnostic argument "+arg);
 }
 if(output.empty())throw std::runtime_error("dedicated probe output required");
 if(synthetic){if(world.size!=1||!input.empty()||!restart.empty())throw std::runtime_error("synthetic mode is one rank with no external state");
  input=(std::filesystem::path(output)/"synthetic-input").string();fixture(input);}
 if(input.empty())throw std::runtime_error("original input partition required");
 m6::apply_cpu_affinity(world,cpus);std::filesystem::create_directories(output);
 auto opts=options(input,restart,synthetic);m6::Solver legacy(world,opts),candidate(world,opts);legacy.validate();candidate.validate();
 // Match run()'s cold initialization before entering a post-flow caller.
 // In particular properties() populates reference_velocity; the retained
 // properties(true,false) later deliberately does not initialize that cache.
 legacy.gradients();legacy.properties();candidate.gradients();candidate.properties();
 int rejected=0;const auto frozen=candidate.lo;
 try{candidate.gradients(&frozen,true);}catch(const std::invalid_argument&){++rejected;}
 candidate.options.normal_sst_transport=false;
 try{candidate.gradients(nullptr,true);}catch(const std::invalid_argument&){++rejected;}
 candidate.options.normal_sst_transport=true;
 if(rejected!=2)throw std::runtime_error("invalid caller contexts were not rejected before mutation");
 uint64_t bytes=full(legacy,candidate),different=0;int cases=synthetic?4:2;
 for(int mode=0;mode<cases;++mode){
  if(mode>=2){
   // Pure synthetic robustness, not an actual production integration state.
   for(int i=0;i<legacy.mesh.n;++i)for(int k=0;k<7;++k){double factor=1.+.003*(k+1)*std::sin(.13*(i+1)+.2*k+mode);
    legacy.w[i][k]*=factor;candidate.w[i][k]*=factor;}
   for(int i=0;i<legacy.mesh.n;++i){legacy.q[i]=legacy.gas.conserved(legacy.w[i]);candidate.q[i]=candidate.gas.conserved(candidate.w[i]);}
   if(mode==3)for(int slot=0;slot<2;++slot)for(size_t fi:legacy.sst_transport->faces){
    double value=(slot?.8:.9)*legacy.sst_transport->boundary[slot][fi];unsigned char mask=static_cast<unsigned char>((fi+slot)%2);
    legacy.sst_transport->boundary[slot][fi]=candidate.sst_transport->boundary[slot][fi]=value;
    legacy.sst_transport->allocated[slot][fi]=candidate.sst_transport->allocated[slot][fi]=mask;
   }
  }
  legacy.update_wall_temperature();candidate.update_wall_temperature();
  legacy.gradients();candidate.gradients(nullptr,true);
  auto checked=partial(legacy,candidate);bytes+=checked.first;different+=checked.second;
  legacy.sst_transport->refresh_farfield(0,legacy.w,legacy.freestream());legacy.reconstruct_normal_scalar(0);
  candidate.sst_transport->refresh_farfield(0,candidate.w,candidate.freestream());candidate.reconstruct_normal_scalar(0);
  bytes+=full(legacy,candidate);
  legacy.properties(true,false);candidate.properties(true,false);bytes+=full(legacy,candidate);
  legacy.assemble(false);candidate.assemble(false);bytes+=full(legacy,candidate);
 }
 auto allbytes=static_cast<uint64_t>(world.sum(static_cast<double>(bytes)));
 auto alldifferent=static_cast<uint64_t>(world.sum(static_cast<double>(different)));
 auto allcells=static_cast<uint64_t>(world.sum(static_cast<double>(legacy.mesh.nt)));
 if(synthetic&&alldifferent==0)throw std::runtime_error("synthetic fixture did not exercise distinct deferred intermediates");
 std::ostringstream filename;filename<<"rank-"<<std::setfill('0')<<std::setw(3)<<world.rank<<".json";
 std::ofstream perrank(std::filesystem::path(output)/filename.str());
 perrank<<"{\"rank\":"<<world.rank<<",\"ranks\":"<<world.size<<",\"cases\":"<<cases<<",\"compared_bytes\":"<<bytes
  <<",\"intermediate_k_different_cell_visits\":"<<different<<",\"owned_cells\":"<<legacy.mesh.n<<",\"owned_and_ghost_cells\":"<<legacy.mesh.nt
  <<",\"caller_context_negatives_rejected\":"<<rejected<<",\"nonlinear_updates\":0}\n";
 perrank.close();if(!perrank)throw std::runtime_error("per-rank result write failed");
 if(world.rank==0){std::ofstream result(std::filesystem::path(output)/"result.json");
  result<<"{\"status\":\"postflow_k_live_stage_bits_passed\",\"ranks\":"<<world.size<<",\"cases\":"<<cases<<",\"compared_bytes\":"<<allbytes
   <<",\"owned_and_ghost_cells\":"<<allcells<<",\"intermediate_k_different_cell_visits\":"<<alldifferent
   <<",\"negative_controls_per_rank\":2,\"synthetic\":"<<(synthetic?"true":"false")<<",\"nonlinear_updates\":0,\"reference_fields_used\":false}\n";
  result.close();if(!result)throw std::runtime_error("global result write failed");
  std::cout<<"M6_POSTFLOW_K_DIRECT cases="<<cases<<" bytes="<<allbytes<<" nonlinear_updates=0 passed\n";
 }
 return 0;
 }catch(const std::exception&e){std::cerr<<"M6_POSTFLOW_K_DIRECT_ERROR rank="<<world.rank<<' '<<e.what()<<std::endl;MPI_Abort(MPI_COMM_WORLD,2);return 2;}}
