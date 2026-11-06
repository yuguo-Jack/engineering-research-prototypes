// Independent small-mesh calls to the real Solver bodies, not a surrogate
// scalar/flow model. Both definitions share identical physics/MPI types.
#ifdef M6_SKIP_DISCARDED_POSTFLOW_ASSEMBLY
#error Build the paired probe without the candidate global definition
#endif
#ifdef M6_SKIP_UNUSED_NORMAL_LSQ
#error This probe isolates post-flow assembly from the LSQ experiment
#endif
#define main m6_fixture_probe_unused_main
#include "m6_dead_normal_lsq_probe.cpp"
#undef main
namespace m6_postflow_candidate {using namespace m6;}
#define m6 m6_postflow_candidate
#define M6_SKIP_DISCARDED_POSTFLOW_ASSEMBLY 1
#include "../src/m6.cpp"
#undef M6_SKIP_DISCARDED_POSTFLOW_ASSEMBLY
#undef m6
#include <cstring>

namespace {
uint64_t compared_bytes=0;int checks=0,dead_cells_different=0;
void check(bool ok,const std::string&name){if(!ok)throw std::runtime_error(name);++checks;}
template<class T>void bytes(const T&a,const T&b,const std::string&name){
 static_assert(std::is_trivially_copyable<T>::value,"value comparison only");
 check(std::memcmp(&a,&b,sizeof(T))==0,name+" bit mismatch");compared_bytes+=sizeof(T);
}
template<class T>void same(const std::vector<T>&a,const std::vector<T>&b,const std::string&name){
 check(a.size()==b.size(),name+" size mismatch");size_t n=a.size()*sizeof(T);
 if(n&&std::memcmp(a.data(),b.data(),n)!=0){auto l=reinterpret_cast<const unsigned char*>(a.data()),r=reinterpret_cast<const unsigned char*>(b.data());
  size_t p=0;while(p<n&&l[p]==r[p])++p;throw std::runtime_error(name+" bit mismatch at byte "+std::to_string(p));}
 ++checks;compared_bytes+=n;
}
template<class A,class B>void sparse(const A&a,const B&b,const std::string&name){
 same(a.row,b.row,name+".row");same(a.col,b.col,name+".col");same(a.diag,b.diag,name+".diag");
 same(a.a,b.a,name+".a");same(a.lu,b.lu,name+".lu");same(a.dinv,b.dinv,name+".dinv");
}
void live(const m6::Solver&a,const m6_postflow_candidate::Solver&b,bool report,bool dead){
#define V(member) same(a.member,b.member,#member)
#define B(member) bytes(a.member,b.member,#member)
 V(w);V(q);V(lo);V(hi);V(grad);V(inverse_lsq);V(scalar_dt);V(prop);V(massflux);V(wall_weight);V(wall_pk);
 V(physical_gradient);V(scalar_gradient);V(reference_velocity);V(total_enthalpy);V(enthalpy_gradient);
 V(thermodynamics);V(thermo_limiter);V(thermo_gradient);V(wall_virtual_jacobian);V(wall_temperature);V(wall_temperature_jacobian);
 V(flow_diffusion_cell_cache);V(flow_diffusion_face_cache);V(recycled_flow_direction);
 B(flow_diffusion_cache_active);B(flow_diffusion_cache_builds);B(scalar_negative_repairs);B(scalar_negative_source_rows);
 B(norms);B(linear_estimate);B(linear_true_residual);B(maximum_linear_true_residual);B(linear_audits);B(linear_target_missed);
 B(turbulence_linear_estimate);B(turbulence_linear_true);B(turbulence_linear_maximum);B(turbulence_linear_steps);
 B(turbulence_linear_audits);B(turbulence_linear_misses);B(scalar_merit);B(limited_updates);B(aa_accepted);B(aa_rejected);
 B(recycle_flow_direction);B(recycled_attempts);B(recycled_accepted);B(recycled_initial_ratio_sum);
 B(linear_start.relative_residual);B(linear_start.coefficient);B(linear_start.trial_relative_residual);B(linear_start.recycled);
 sparse(a.flow,b.flow,"flow");sparse(a.scalar,b.scalar,"scalar");
 check(a.residual.size()==b.residual.size(),"residual sizes");
 for(size_t i=0;i<a.residual.size();++i){bytes(a.residual[i][0],b.residual[i][0],"retained continuity residual");
  if(std::memcmp(a.residual[i].data()+1,b.residual[i].data()+1,4*sizeof(double)))++dead_cells_different;}
 if(dead){V(residual);V(pseudo);}
 if(report){B(force);B(pressure_force);B(viscous_force);}
 check(bool(a.sst_transport)==bool(b.sst_transport),"SST ownership differs");
 if(a.sst_transport)for(int k=0;k<2;++k){same(a.sst_transport->reference[k],b.sst_transport->reference[k],"SST.reference");
  same(a.sst_transport->boundary[k],b.sst_transport->boundary[k],"SST.boundary");same(a.sst_transport->allocated[k],b.sst_transport->allocated[k],"SST.allocated");}
#undef V
#undef B
}
template<class O>O options(const std::filesystem::path&input){
 O o;o.input=input.string();o.iterations=0;o.checkpoint=0;o.linear_steps=24;o.linear_tolerance=.1;o.linear_audit_period=1;
 o.flow_relax=.75;o.turb_relax=.8;o.batched_orthogonalization=true;o.steady_turbulence=true;o.turbulence_k_first=true;
 o.coupled_ausm=true;o.turbulence_gauss=true;o.harmonic_scalars=true;o.independent_face_thermo=true;o.boundary_face_ray=true;
 o.exact_flux_jacobian=true;o.characteristic_farfield=true;o.wall_slip_gradients=true;o.coupled_wall_reconstruction=true;
 o.symmetry_normal_lsq=true;o.enhanced_thermal_wall=true;o.symmetry_face_gradients=true;o.audit_turbulence_linear=true;
 o.coupled_farfield_reconstruction=true;o.symmetry_viscous_flux=true;o.farfield_zero_viscous_flux=true;
 o.normal_sst_transport=true;o.lagged_sst_mass_flux=true;o.normal_coupled_timestep=true;o.normal_scalar_step=true;
 o.local_preconditioning=true;o.normal_approximate_jacobian=true;o.normal_farfield_jacobian=true;o.primitive_flow_update=true;
 o.skip_unused_face_work=true;o.normal_wall_jacobian=true;o.normal_symmetry_jacobian=true;
 return o;
}
void small_input(const std::filesystem::path&input){
 fixture(input);constexpr int n=27;
 // Entirely synthetic nonuniform flow, written in the shard's PHYSICAL
 // rho/u/v/w/T/k/omega convention; Mesh then nondimensionalizes it.
 const std::streamoff start=8+4*4+8+48*8+n*(4+8+3*8+8);
 std::fstream mesh(input/"rank-000.bin",std::ios::binary|std::ios::in|std::ios::out);mesh.seekp(start);
 for(int z=0;z<3;++z)for(int y=0;y<3;++y)for(int x=0;x<3;++x){m6::W7 v={.715+.002*x+.001*y-.001*z,.55+.012*x-.005*z,.01+.002*y-.001*x,.006-.001*z+.0005*y,
   1.02+.004*x-.002*y+.003*z,.0002+.00001*x+.00002*y+.000015*z,.8+.03*x-.02*y+.01*z};
  constexpr double rho0=1.1361339060496805,a0=352.517926811456,T0=298.15,R=299.124701264874;
  v[0]*=rho0*a0*a0/(R*T0*v[4]);for(int j=1;j<=3;++j)v[j]*=a0;v[4]*=T0;v[5]*=a0*a0;v[6]*=a0;
  mesh.write(reinterpret_cast<const char*>(v.data()),sizeof(v));}
 check(bool(mesh),"synthetic initial-state write failed");mesh.close();
 std::fstream sidecar(input/"scalar-boundary-init-rank-000.bin",std::ios::binary|std::ios::in|std::ios::out);
 char magic[8];int32_t h[6];sidecar.read(magic,8);sidecar.read(reinterpret_cast<char*>(h),sizeof(h));sidecar.seekp(8+sizeof(h)+h[1]*2*sizeof(int32_t));
 for(int i=0;i<h[1];++i){std::array<double,2>v={(.0002+.00001*(i%7))*352.517926811456*352.517926811456,(.8+.03*(i%11))*352.517926811456};sidecar.write(reinterpret_cast<const char*>(v.data()),sizeof(v));}
 check(bool(sidecar),"synthetic boundary write failed");
}
template<class S>void initial(S&s){s.validate();s.update_wall_temperature();s.gradients();s.properties();s.assemble(true);}
template<class S>void finish(S&s){s.gradients();s.properties();s.forces_only();for(int i=0;i<s.mesh.n;++i)s.q[i]=s.gas.conserved(s.w[i]);}
std::string read_file(const std::filesystem::path&path){std::ifstream f(path,std::ios::binary);check(bool(f),"missing output "+path.string());return {std::istreambuf_iterator<char>(f),{}};}
void files(const std::filesystem::path&a,const std::filesystem::path&b){
 for(const char*name:{"state-rank-000.bin","wall-thermal-rank-000.bin","scalar-boundary-rank-000.bin","wall-rank-000.csv","linear-audit.csv","turbulence-linear-audit.csv"})
  check(read_file(a/name)==read_file(b/name),std::string("run output differs ")+name);
}
template<class F>std::string error(F&&f){try{f();}catch(const std::exception&e){return e.what();}return {};}
}

int main(int argc,char**argv){m6::MPIWorld world(argc,argv);try{
 if(argc!=2||world.size!=1)throw std::runtime_error("probe requires one rank and a new dedicated output directory");
 std::filesystem::path out=argv[1];if(std::filesystem::exists(out))throw std::runtime_error("probe output must be new; retain earlier failures");
 std::filesystem::create_directories(out);small_input(out/"input");auto off=options<m6::Options>(out/"input");auto on=options<m6_postflow_candidate::Options>(out/"input");
 m6::Solver a(world,off);m6_postflow_candidate::Solver b(world,on);initial(a);initial(b);live(a,b,true,true);
 check(!a.may_skip_discarded_postflow_assembly(true),"OFF flag activated");check(b.may_skip_discarded_postflow_assembly(true),"ON production guard missed");
 for(int step=1;step<=1;++step){
  std::cout<<"PROBE_STAGE paired-step="<<step<<" meanflow"<<std::endl;
  a.meanflow(true,true);b.meanflow(true,true);live(a,b,true,true);
  std::cout<<"PROBE_STAGE paired-step="<<step<<" advance-turbulence"<<std::endl;
  a.advance_turbulence(true,true);b.advance_turbulence(true,true);live(a,b,false,false);
  check(b.skipped_postflow_assemblies==static_cast<uint64_t>(step),"production skip count wrong");
  std::cout<<"PROBE_STAGE paired-step="<<step<<" finish"<<std::endl;
  finish(a);finish(b);live(a,b,true,false);
  // Diagnostic RHS queries are separate from the actual production updates.
  for(int variable=5;variable<=6;++variable){std::cout<<"PROBE_STAGE paired-step="<<step<<" diagnostic-rhs="<<variable<<std::endl;m6::Field<1>ar,br;a.turbulence(variable,false,false,&ar);b.turbulence(variable,false,false,&br);same(ar,br,"post-update scalar RHS");live(a,b,true,false);}
  a.update_wall_temperature();b.update_wall_temperature();a.assemble(true);b.assemble(true);live(a,b,true,true);
 }
 check(dead_cells_different>0,"probe failed to distinguish discarded residual slots");
 // Direct-call and changed Jacobian lifecycle fallbacks preserve ALL outputs.
 std::ofstream fallbacks(out/"fallback-checks.json");fallbacks<<"[";
 for(int mode=0;mode<8;++mode){auto oa=off;auto ob=on;if(mode==1){oa.jacobian_period=2;ob.jacobian_period=2;}
  if(mode==2){oa.normal_scalar_step=false;ob.normal_scalar_step=false;}
  if(mode>=3){oa.normal_scalar_step=false;ob.normal_scalar_step=false;oa.lagged_sst_mass_flux=false;ob.lagged_sst_mass_flux=false;}
  if(mode>=4){oa.normal_coupled_timestep=false;ob.normal_coupled_timestep=false;}
  if(mode>=5){oa.normal_sst_transport=false;ob.normal_sst_transport=false;}
  if(mode==6){oa.turbulence_k_first=false;ob.turbulence_k_first=false;}
  if(mode==7){oa.enhanced_thermal_wall=false;ob.enhanced_thermal_wall=false;oa.anderson_depth=1;ob.anderson_depth=1;}
  m6::Solver x(world,oa);m6_postflow_candidate::Solver y(world,ob);initial(x);initial(y);x.meanflow(true,true);y.meanflow(true,true);
  auto ex=error([&]{x.advance_turbulence(true,mode!=0);});auto ey=error([&]{y.advance_turbulence(true,mode!=0);});
  check(ex==ey,"fallback first exceptions differ "+std::to_string(mode));live(x,y,true,true);check(y.skipped_postflow_assemblies==0,"fallback skipped assembly");
  if(mode)fallbacks<<',';fallbacks<<"{\"case\":"<<mode<<",\"skipped\":0,\"all_bytes_same\":true,\"both_returned\":"<<(ex.empty()?"true":"false")<<",\"same_error\":\""<<ex<<"\"}";
 }fallbacks<<"]\n";fallbacks.close();check(bool(fallbacks),"fallback evidence write failed");
 // Every boundary guard is checked explicitly; invalid constructor option
 // combinations are never integrated to manufacture a successful fallback.
 auto original=b.options;
 auto rejects=[&](const char*name){check(!b.may_skip_discarded_postflow_assembly(true),std::string("guard accepted ")+name);b.options=original;};
 b.options.normal_sst_transport=false;rejects("non-normal");b.options.normal_scalar_step=false;rejects("non-normal scalar");
 b.options.lagged_sst_mass_flux=false;rejects("unretained mass");b.options.normal_coupled_timestep=false;rejects("unretained timestep");
 b.options.anderson_depth=1;rejects("AA");b.options.jacobian_period=2;rejects("reused Jacobian");b.options.validate=true;rejects("validate");
 b.options.check_viscous=true;rejects("viscous diagnostic");b.options.check_linear=true;rejects("linear diagnostic");b.options.export_linear_operator=true;rejects("operator diagnostic");
 check(!b.may_skip_discarded_postflow_assembly(false),"direct diagnostic accepted");
 // Same real run() callers, all checkpoint/cache/wall/audit output bytes.
 off.iterations=1;on.iterations=1;off.output=(out/"run-off").string();on.output=(out/"run-on").string();
 m6::Solver run_a(world,off);m6_postflow_candidate::Solver run_b(world,on);run_a.run();run_b.run();live(run_a,run_b,true,false);files(off.output,on.output);
 check(run_b.skipped_postflow_assemblies==1,"real run did not activate exactly once");
 for(const auto&v:run_b.w){for(double x:v)check(std::isfinite(x),"nonfinite completed state");check(v[0]>0&&v[4]>0&&v[5]>0&&v[6]>0,"nonpositive completed state");}
 // Exercise the actual run() guard and the inherited019 report-stage caller,
 // not only a direct advance_turbulence surrogate. Preserve the original
 // default-report run and all later fallback/error cases unchanged.
 auto stage_off=off;auto stage_on=on;
 stage_off.coupled_force_report_stage=true;stage_on.coupled_force_report_stage=true;
 stage_off.output=(out/"stage-true-run-off").string();stage_on.output=(out/"stage-true-run-on").string();
 m6::Solver stage_a(world,stage_off);m6_postflow_candidate::Solver stage_b(world,stage_on);
 stage_a.run();stage_b.run();live(stage_a,stage_b,true,false);files(stage_off.output,stage_on.output);
 check(stage_a.skipped_postflow_assemblies==0,"stageTRUE OFF run skipped assembly");
 check(stage_b.skipped_postflow_assemblies==1,"stageTRUE ON real run did not activate exactly once");
 for(const auto&v:stage_b.w){for(double x:v)check(std::isfinite(x),"nonfinite stageTRUE completed state");check(v[0]>0&&v[4]>0&&v[5]>0&&v[6]>0,"nonpositive stageTRUE completed state");}
 std::ofstream stage(out/"report-stage-coverage.json");
 stage<<"{\"default_report_paired_run_updates\":1,\"coupled_stage_true_paired_run_updates\":1,\"actual_run_called_both_modes\":true,\"stage_true_OFF_skips\":0,\"stage_true_ON_skips\":1,\"all_live_and_payload_bytes_equal\":true,\"reference_fields_used\":false}\n";
 stage.close();check(bool(stage),"report-stage evidence write failed");
 // This deliberately retained artificial boundary configuration is rejected
 // at the second scalar update by the original code. Run BOTH arms fully to
 // that same error instead of tuning the physics/limiter until it passes.
 off.iterations=2;on.iterations=2;off.output=(out/"invalid-two-step-off").string();on.output=(out/"invalid-two-step-on").string();
 m6::Solver invalid_a(world,off);m6_postflow_candidate::Solver invalid_b(world,on);
 auto two_off=error([&]{invalid_a.run();});auto two_on=error([&]{invalid_b.run();});
 check(!two_off.empty()&&two_off==two_on,"second-step original and candidate errors differ OFF="+two_off+" ON="+two_on);
 check(two_off.find("nonpositive scalar face state")!=std::string::npos,"unexpected second-step rejection stage");
 live(invalid_a,invalid_b,false,false);
 std::ofstream two(out/"second-step-original-domain-rejection.json");two<<"{\"both_arms_executed\":true,\"completed_steps_each\":1,\"rejected_step\":2,\"stage\":\"advance_turbulence_scalar\",\"same_error\":true,\"off_error\":\""<<two_off<<"\",\"on_error\":\""<<two_on<<"\"}\n";two.close();check(bool(two),"second-step evidence write failed");
 // Original domain errors must remain observable with identical messages.
 std::ofstream failures(out/"domain-checks.json");failures<<"[";bool comma=false;
 for(int mode=0;mode<16;++mode){auto oa=options<m6::Options>(out/"input");auto ob=options<m6_postflow_candidate::Options>(out/"input");
  auto enthalpy_options=[](auto&o){o.reconstruct_enthalpy=true;o.independent_face_thermo=false;o.exact_flux_jacobian=false;
   o.normal_approximate_jacobian=false;o.normal_farfield_jacobian=false;o.normal_wall_jacobian=false;o.normal_symmetry_jacobian=false;};
  if(mode==15){enthalpy_options(oa);enthalpy_options(ob);}
  m6::Solver x(world,oa);m6_postflow_candidate::Solver y(world,ob);initial(x);initial(y);
  auto corrupt=[&](auto&s){switch(mode){
   case 0:s.w[0][0]=-1;break;case 1:s.prop[0][0]=-1;break;case 2:s.gas.mu=0;break;case 3:s.mesh.meta[26]=0;break;
   case 4:s.lo[0][0]=1;s.grad[0][0]=-1e8;break;case 5:s.lo[0][4]=1;s.grad[0][12]=-1e8;break;
   case 6:s.thermodynamics[0][0]=-1;s.thermo_limiter[0]={0,0};break;
   case 7:s.thermodynamics[0][1]=std::numeric_limits<double>::quiet_NaN();s.thermo_limiter[0]={0,0};break;
   case 8:s.reference_velocity[0][0]=0;break;
   case 9:s.mesh.faces[0].area=0;break;
   case 10:{auto f=s.mesh.faces[0];s.mesh.center[f.r]=s.mesh.center[f.l];break;}
   case 11:for(auto&f:s.mesh.faces)if(f.type==7){f.dn=0;break;}break;
   case 12:s.mesh.faces[0].r=-1;s.mesh.faces[0].type=9;s.w[0][1]=std::numeric_limits<double>::quiet_NaN();break;
   case 13:s.mesh.faces[0].r=-1;s.mesh.faces[0].type=9;s.w[0][1]=std::numeric_limits<double>::infinity();break;
   case 14:s.w[0][4]=.5;s.gas.cp0=100;break;
   case 15:s.total_enthalpy[0]=std::numeric_limits<double>::quiet_NaN();break;
  }};corrupt(x);corrupt(y);auto ex=error([&]{x.assemble(false);}),ey=error([&]{y.validate_discarded_postflow_assembly();});
  check(!ex.empty(),"invalid case did not throw "+std::to_string(mode));check(ex==ey,"domain exception differs "+std::to_string(mode)+" OFF="+ex+" ON="+ey);
  if(comma)failures<<',';comma=true;failures<<"{\"case\":"<<mode<<",\"same_error\":true,\"error\":\""<<ex<<"\"}";
 }failures<<"]\n";failures.close();check(bool(failures),"domain evidence write failed");
 std::ofstream result(out/"result.json");result<<"{\"status\":\"small_mesh_live_bits_and_domain_checks_passed\",\"ranks\":1,\"synthetic_cells\":27,\"paired_stage_updates\":1,\"paired_run_updates\":1,\"coupled_stage_true_paired_run_updates\":1,\"total_paired_real_run_updates\":2,\"checks\":"<<checks<<",\"compared_bytes\":"<<compared_bytes<<",\"discarded_residual_cell_differences\":"<<dead_cells_different<<",\"domain_cases\":16,\"same_second_step_domain_rejection\":true,\"external_case_inputs\":false,\"real32_validated\":false,\"performance_claim\":false}\n";
 result.close();check(bool(result),"result write failed");std::cout<<"M6_SKIP_POSTFLOW_ASSEMBLY_PROBE checks="<<checks<<" compared_bytes="<<compared_bytes<<" domain_cases=16 complete\n";return 0;
 }catch(const std::exception&e){std::cerr<<"M6_SKIP_POSTFLOW_ASSEMBLY_PROBE_ERROR "<<e.what()<<std::endl;MPI_Abort(MPI_COMM_WORLD,2);return 2;}}
