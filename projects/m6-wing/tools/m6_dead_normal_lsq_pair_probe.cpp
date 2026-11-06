// Test-only two-Solver comparison. Both variants use the actual same source
// body, common MPI/physics/mesh types, and no changes to the production API.
// Large real 32-rank arrays are compared directly in memory, not serialized.
#ifdef M6_SKIP_UNUSED_NORMAL_LSQ
#error The dual probe must compile its first Solver without the candidate flag
#endif
#define main m6_fixture_probe_unused_main
#include "m6_dead_normal_lsq_probe.cpp"
#undef main

namespace m6_lsq_candidate {using namespace m6;}
#define m6 m6_lsq_candidate
#define M6_SKIP_UNUSED_NORMAL_LSQ 1
#include "../src/m6.cpp"
#undef M6_SKIP_UNUSED_NORMAL_LSQ
#undef m6
#include <cstring>

namespace {
template<class T>uint64_t same(const std::vector<T>&a,const std::vector<T>&b,const std::string&name){
 if(a.size()!=b.size())throw std::runtime_error(name+" vector lengths differ");
 size_t bytes=a.size()*sizeof(T);
 if(bytes&&std::memcmp(a.data(),b.data(),bytes)!=0){
  auto left=reinterpret_cast<const unsigned char*>(a.data()),right=reinterpret_cast<const unsigned char*>(b.data());
  size_t first=0;while(first<bytes&&left[first]==right[first])++first;
  throw std::runtime_error(name+" bit mismatch at byte "+std::to_string(first));
 }
 return bytes;
}
struct Check{uint64_t bytes=0,scratch_cells_different=0;};
Check same_live(const m6::Solver&a,const m6_lsq_candidate::Solver&b,bool compare_dead){
 Check r;
 r.bytes+=same(a.w,b.w,"w");r.bytes+=same(a.q,b.q,"q");r.bytes+=same(a.lo,b.lo,"lo");r.bytes+=same(a.grad,b.grad,"grad");
 if(a.hi.size()!=b.hi.size())throw std::runtime_error("hi sizes differ");
 for(size_t i=0;i<a.hi.size();++i){
  if(std::memcmp(a.hi[i].data(),b.hi[i].data(),5*sizeof(double))!=0)throw std::runtime_error("live hi meanflow columns differ");
  r.bytes+=5*sizeof(double);
  if(std::memcmp(a.hi[i].data()+5,b.hi[i].data()+5,2*sizeof(double))!=0){
   if(compare_dead)throw std::runtime_error("frozen fallback hi scalar columns differ");
   ++r.scratch_cells_different;
  }
  if(compare_dead)r.bytes+=2*sizeof(double);
 }
 r.bytes+=same(a.physical_gradient,b.physical_gradient,"physical_gradient");r.bytes+=same(a.scalar_gradient,b.scalar_gradient,"scalar_gradient");
 r.bytes+=same(a.thermodynamics,b.thermodynamics,"thermodynamics");r.bytes+=same(a.thermo_gradient,b.thermo_gradient,"thermo_gradient");
 r.bytes+=same(a.thermo_limiter,b.thermo_limiter,"thermo_limiter");r.bytes+=same(a.total_enthalpy,b.total_enthalpy,"total_enthalpy");
 r.bytes+=same(a.enthalpy_gradient,b.enthalpy_gradient,"enthalpy_gradient");r.bytes+=same(a.wall_virtual_jacobian,b.wall_virtual_jacobian,"wall_virtual_jacobian");
 r.bytes+=same(a.wall_temperature,b.wall_temperature,"wall_temperature");r.bytes+=same(a.wall_temperature_jacobian,b.wall_temperature_jacobian,"wall_temperature_jacobian");
 if(!a.sst_transport||!b.sst_transport)throw std::runtime_error("dual probe requires both normal SST caches");
 for(int slot=0;slot<2;++slot){r.bytes+=same(a.sst_transport->reference[slot],b.sst_transport->reference[slot],"sst.reference");
  r.bytes+=same(a.sst_transport->boundary[slot],b.sst_transport->boundary[slot],"sst.boundary");
  r.bytes+=same(a.sst_transport->allocated[slot],b.sst_transport->allocated[slot],"sst.allocated");}
 return r;
}
template<class O>O actual_options(const std::string&input,const std::string&restart,bool synthetic){
 // Exact non-timing 970302 production option values. No reference states.
 O o;o.input=input;o.restart=restart;o.iterations=0;o.checkpoint=0;
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
 for(int i=1;i<argc;++i){std::string argument=argv[i];
  auto value=[&](){if(++i>=argc)throw std::runtime_error("missing probe argument "+argument);return std::string(argv[i]);};
  if(argument=="--input")input=value();else if(argument=="--output")output=value();else if(argument=="--restart")restart=value();
  else if(argument=="--rank-cpus")cpus=value();else if(argument=="--synthetic"){synthetic=true;output=value();}
  else throw std::runtime_error("unknown dual gradient probe argument "+argument);
 }
 if(output.empty())throw std::runtime_error("dedicated output directory required");
 if(synthetic){if(world.size!=1||!input.empty()||!restart.empty())throw std::runtime_error("synthetic probe is one rank with no external state");
  input=(std::filesystem::path(output)/"synthetic-input").string();fixture(input);}
 if(input.empty())throw std::runtime_error("original partition required");
 m6::apply_cpu_affinity(world,cpus);std::filesystem::create_directories(output);
 auto legacy_options=actual_options<m6::Options>(input,restart,synthetic);
 auto candidate_options=actual_options<m6_lsq_candidate::Options>(input,restart,synthetic);
 m6::Solver legacy(world,legacy_options);m6_lsq_candidate::Solver candidate(world,candidate_options);
 legacy.validate();candidate.validate();
 uint64_t compared=0,dead_difference=0;
 for(int mode=0;mode<4;++mode){
  auto legacy_frozen=legacy.lo,candidate_frozen=candidate.lo;
  if(mode<2){legacy.gradients();candidate.gradients();}
  else if(mode==2){legacy.gradients(&legacy_frozen);candidate.gradients(&candidate_frozen);}
  else{legacy.gradients(&legacy.lo);candidate.gradients(&candidate.lo);}
  Check checked=same_live(legacy,candidate,mode>=2);compared+=checked.bytes;dead_difference+=checked.scratch_cells_different;
 }
 // The comparisons include all ghost entries and the two frozen fallbacks.
 std::ostringstream name;name<<"rank-"<<std::setfill('0')<<std::setw(3)<<world.rank<<".json";
 std::ofstream rank(std::filesystem::path(output)/name.str());rank<<"{\"rank\":"<<world.rank<<",\"ranks\":"<<world.size
  <<",\"owned_cells\":"<<legacy.mesh.n<<",\"owned_and_ghost_cells\":"<<legacy.mesh.nt<<",\"compared_bytes\":"<<compared
  <<",\"dead_scratch_cell_disagreements\":"<<dead_difference<<",\"nonnull_frozen_hi_also_identical\":true,\"comparisons_passed\":4}\n";
 rank.close();if(!rank)throw std::runtime_error("rank result write failed");
 auto total_bytes=static_cast<uint64_t>(world.sum(static_cast<double>(compared)));
 auto total_dead=static_cast<uint64_t>(world.sum(static_cast<double>(dead_difference)));
 auto total_cells=static_cast<uint64_t>(world.sum(static_cast<double>(legacy.mesh.nt)));
 if(world.rank==0){std::ofstream result(std::filesystem::path(output)/"result.json");
  result<<"{\"status\":\"live_arrays_direct_memcmp_passed\",\"ranks\":"<<world.size<<",\"owned_and_ghost_cells\":"<<total_cells
   <<",\"comparisons_passed\":4,\"compared_bytes\":"<<total_bytes<<",\"dead_scratch_cell_disagreements\":"<<total_dead
   <<",\"synthetic_fixture\":"<<(synthetic?"true":"false")<<",\"nonlinear_updates\":0,\"reference_state_used\":false}\n";
  result.close();if(!result)throw std::runtime_error("global probe result write failed");
  std::cout<<"M6_DEAD_NORMAL_LSQ_DIRECT ranks="<<world.size<<" comparisons=4 bytes="<<total_bytes<<" nonlinear_updates=0 passed\n";
 }
 return 0;
 }catch(const std::exception&e){std::cerr<<"M6_DEAD_NORMAL_LSQ_DIRECT_ERROR rank="<<world.rank<<' '<<e.what()<<std::endl;MPI_Abort(MPI_COMM_WORLD,2);return 2;}}
