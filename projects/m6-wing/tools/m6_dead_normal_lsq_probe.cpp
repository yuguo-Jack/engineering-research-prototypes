// Small original synthetic meshes exercise the actual Solver::gradients body.
// No nonlinear updates, native/reference states or vendor libraries are used.
#define M6_NO_MAIN
#include "../src/m6.cpp"
#include <type_traits>

namespace {
template<class T>void write_items(std::ofstream&out,const std::vector<T>&v){
 static_assert(std::is_trivially_copyable<T>::value,"binary probe needs value objects");
 uint64_t n=v.size();out.write(reinterpret_cast<const char*>(&n),sizeof(n));
 out.write(reinterpret_cast<const char*>(v.data()),static_cast<std::streamsize>(v.size()*sizeof(T)));
 if(!out)throw std::runtime_error("probe snapshot write failed");
}
template<size_t N>void field(std::ofstream&out,const std::vector<std::array<double,N>>&v){
 for(const auto&a:v)for(double x:a)if(!std::isfinite(x))throw std::runtime_error("nonfinite probe field");
 write_items(out,v);
}
void field(std::ofstream&out,const std::vector<double>&v){
 for(double x:v)if(!std::isfinite(x))throw std::runtime_error("nonfinite scalar probe field");
 write_items(out,v);
}

void fixture(const std::filesystem::path&path){
 std::filesystem::create_directories(path);
 constexpr int side=3,n=side*side*side;
 std::array<double,48>meta{};
 meta[0]=1.1361339060496805;meta[1]=352.517926811456;meta[2]=298.15;meta[3]=299.124701264874;
 meta[4]=101325;meta[5]=1.72e-5;meta[6]=.0454;meta[7]=1059.4853515625;
 std::array<double,5>low={985.1700244,.4572524926,-.0010998299666,1.5381987004e-6,-6.675425469e-10};
 std::array<double,5>high={901.5941633,.4449195144,-.00015867025232,2.7184049518e-8,-1.7782013004e-12};
 for(int i=0;i<5;++i){meta[8+i]=low[i];meta[13+i]=high[i];}
 meta[26]=.85;meta[27]=.075;meta[28]=.0828;meta[29]=1/.85;meta[30]=1.;meta[31]=2.;meta[32]=1/.856;
 meta[33]=.31;meta[34]=.09;meta[35]=.41;meta[36]=10;
 meta[39]=.83;meta[40]=1;meta[41]=.04;meta[42]=.1;meta[43]=298.15;meta[44]=101325;meta[45]=.04;meta[46]=80;meta[47]=1;
 std::vector<int>ids;std::vector<double>volumes,distances;
 std::vector<m6::D3>centers;std::vector<m6::W7>initial;
 const m6::D3 widths={.04,.03,.02};
 auto index=[=](int x,int y,int z){return x+side*(y+side*z);};
 for(int z=0;z<side;++z)for(int y=0;y<side;++y)for(int x=0;x<side;++x){
  ids.push_back(index(x,y,z));volumes.push_back(widths[0]*widths[1]*widths[2]);distances.push_back(.005+.002*y);
  centers.push_back({(x+.5)*widths[0],(y+.5)*widths[1],(z+.5)*widths[2]});
  initial.push_back({1.05+.018*x+.009*y-.006*z,150.+9*x-4*z,4.+2*y-1.5*x,7.-2*z+.75*y,
    304.+5*x-2*y+3*z,.025+.019*x+.007*y+.011*z,50.+12*x-7*y+4*z});
 }
 std::vector<int>left,right,type,face_id;
 std::vector<m6::D3>face_center,area_vector;
 auto face=[&](int cell,int neighbor,int axis,int sign,int kind){
  auto center=centers[cell];center[axis]+=.5*sign*widths[axis];
  m6::D3 area{};area[axis]=sign*volumes[cell]/widths[axis];
  left.push_back(cell);right.push_back(neighbor);type.push_back(kind);face_id.push_back(static_cast<int>(face_id.size()));
  face_center.push_back(center);area_vector.push_back(area);
 };
 for(int z=0;z<side;++z)for(int y=0;y<side;++y)for(int x=0;x<side;++x){
  int cell=index(x,y,z);int coordinate[3]={x,y,z};
  for(int axis=0;axis<3;++axis){
   if(coordinate[axis]<side-1){int next[3]={x,y,z};++next[axis];face(cell,index(next[0],next[1],next[2]),axis,1,2);}
   if(coordinate[axis]==0)face(cell,-1,axis,-1,axis==1?3:axis==2?7:9);
   if(coordinate[axis]==side-1)face(cell,-1,axis,1,axis==1?3:9);
  }
 }
 std::ofstream mesh(path/"rank-000.bin",std::ios::binary);mesh.write("CFDM6A01",8);
 uint32_t counts[4]={n,n,static_cast<uint32_t>(left.size()),0};uint64_t global=n;
 mesh.write(reinterpret_cast<const char*>(counts),sizeof(counts));mesh.write(reinterpret_cast<const char*>(&global),sizeof(global));
 mesh.write(reinterpret_cast<const char*>(meta.data()),sizeof(meta));
 // Mesh::read format deliberately has no per-vector length prefixes.
 auto raw=[&](const auto&v){mesh.write(reinterpret_cast<const char*>(v.data()),static_cast<std::streamsize>(v.size()*sizeof(v[0])));};
 raw(ids);raw(volumes);raw(centers);raw(distances);raw(initial);
 raw(left);raw(right);raw(type);raw(face_id);raw(face_center);raw(area_vector);
 mesh.close();if(!mesh)throw std::runtime_error("synthetic mesh write failed");
 std::vector<int32_t>boundary_ids,masks;std::vector<std::array<double,2>>values;
 for(size_t i=0;i<left.size();++i)if(right[i]<0&&(type[i]==3||type[i]==9)){
  boundary_ids.push_back(face_id[i]);masks.push_back(static_cast<int32_t>(i%4));
  values.push_back({.035+.004*(i%7),70.+3*(i%11)});
 }
 std::ofstream sidecar(path/"scalar-boundary-init-rank-000.bin",std::ios::binary);sidecar.write("M6SBI001",8);
 int32_t header[6]={1,static_cast<int32_t>(values.size()),0,1,n,static_cast<int32_t>(values.size())};
 sidecar.write(reinterpret_cast<const char*>(header),sizeof(header));
 sidecar.write(reinterpret_cast<const char*>(boundary_ids.data()),boundary_ids.size()*sizeof(boundary_ids[0]));
 sidecar.write(reinterpret_cast<const char*>(masks.data()),masks.size()*sizeof(masks[0]));
 sidecar.write(reinterpret_cast<const char*>(values.data()),values.size()*sizeof(values[0]));
 sidecar.close();if(!sidecar)throw std::runtime_error("synthetic scalar input write failed");
}

struct Mode{const char*name;bool normal,coupled,enthalpy,scalar_gauss;int order,frozen;};
const std::array<Mode,11>modes={{{"normal_basic",true,false,false,true,2,0},
 {"normal_order1",true,false,false,true,1,0},{"normal_thermal",true,true,false,true,2,0},
 {"normal_enthalpy",true,true,true,true,2,0},{"normal_frozen_copy",true,true,false,true,2,1},
 {"normal_frozen_lo_alias",true,true,false,true,2,2},{"legacy_gauss",false,false,false,true,2,0},
 {"legacy_lsq",false,false,false,false,2,0},{"legacy_thermal",false,true,false,true,2,0},
 {"legacy_frozen_copy",false,false,false,true,2,1},{"legacy_frozen_lo_alias",false,false,false,false,1,2}}};

void snapshot(const m6::Solver&s,const std::filesystem::path&path,const std::filesystem::path&scratch){
 std::ofstream out(path,std::ios::binary);out.write("M6GLIVE1",8);
 field(out,s.w);field(out,s.q);field(out,s.lo);field(out,s.grad);
 std::vector<m6::Q5>live_hi(s.hi.size());std::vector<std::array<double,2>>dead_hi(s.hi.size());
 for(size_t i=0;i<s.hi.size();++i){for(int j=0;j<5;++j)live_hi[i][j]=s.hi[i][j];dead_hi[i]={s.hi[i][5],s.hi[i][6]};}
 field(out,live_hi);field(out,s.physical_gradient);field(out,s.scalar_gradient);
 field(out,s.thermodynamics);field(out,s.thermo_gradient);field(out,s.thermo_limiter);
 field(out,s.total_enthalpy);field(out,s.enthalpy_gradient);field(out,s.wall_virtual_jacobian);
 field(out,s.wall_temperature);field(out,s.wall_temperature_jacobian);
 if(s.sst_transport)for(int slot=0;slot<2;++slot){field(out,s.sst_transport->reference[slot]);
  field(out,s.sst_transport->boundary[slot]);write_items(out,s.sst_transport->allocated[slot]);}
 out.close();if(!out)throw std::runtime_error("live snapshot failed");
 std::ofstream dead(scratch,std::ios::binary);field(dead,dead_hi);dead.close();if(!dead)throw std::runtime_error("scratch snapshot failed");
}
}

int main(int argc,char**argv){m6::MPIWorld world(argc,argv);try{
 if(argc!=2||world.size!=1)throw std::runtime_error("probe requires one rank and a dedicated output directory");
 std::filesystem::path root=argv[1];fixture(root/"synthetic-input");
 std::ofstream metadata(root/"manifest.json");metadata<<"{\"fixture\":\"synthetic-cartesian-27\",\"nonlinear_updates\":0,\"ranks\":1,\"cases\":[";
 bool comma=false;
 for(const auto&mode:modes){
  m6::Options o;o.input=(root/"synthetic-input").string();o.iterations=0;o.order=mode.order;
  o.normal_sst_transport=mode.normal;o.turbulence_k_first=true;o.steady_turbulence=true;
  o.coupled_ausm=true;o.harmonic_scalars=true;o.characteristic_farfield=true;o.turbulence_gauss=mode.scalar_gauss;
  o.wall_slip_gradients=mode.coupled;o.coupled_wall_reconstruction=mode.coupled;
  o.boundary_face_ray=mode.coupled;o.enhanced_thermal_wall=mode.coupled;
  o.symmetry_normal_lsq=mode.coupled;o.symmetry_face_gradients=mode.coupled;
  o.coupled_farfield_reconstruction=mode.coupled;o.independent_face_thermo=!mode.enthalpy;
  o.reconstruct_enthalpy=mode.enthalpy;
  m6::Solver solver(world,o);solver.validate();
  for(int pass=0;pass<2;++pass){
   // Change the scalar/meanflow patterns between calls to exercise scratch reset.
   for(int i=0;i<solver.mesh.n;++i){for(int j=0;j<7;++j)solver.w[i][j]*=1.+.009*(j+1)*std::sin(.4*(i+1)+.3*j+.2*pass);
    solver.q[i]=solver.gas.conserved(solver.w[i]);}
   solver.update_wall_temperature();
   // Use an admissible limiter snapshot, as the real directional checks do.
   // Arbitrary large frozen values can make the original scalar face negative.
   if(mode.frozen)solver.gradients();
   std::vector<m6::W7>frozen=solver.lo;
   if(mode.frozen==2)solver.lo=frozen;
   solver.update_wall_temperature();solver.gradients(mode.frozen==1?&frozen:mode.frozen==2?&solver.lo:nullptr);
   std::string name=std::string(mode.name)+"-"+std::to_string(pass);
   snapshot(solver,root/(name+".live.bin"),root/(name+".scratch.bin"));
   if(comma)metadata<<',';comma=true;
   metadata<<"{\"name\":\""<<name<<"\",\"optimized_eligible\":"<<(mode.normal&&mode.frozen==0?"true":"false")<<'}';
  }
 }
 metadata<<"]}\n";metadata.close();if(!metadata)throw std::runtime_error("probe manifest failed");
 std::cout<<"M6_DEAD_NORMAL_LSQ_PROBE modes="<<modes.size()<<" passes=2 cells=27 nonlinear_updates=0 complete\n";
 return 0;
 }catch(const std::exception&e){std::cerr<<"M6_DEAD_NORMAL_LSQ_PROBE_ERROR "<<e.what()<<std::endl;MPI_Abort(MPI_COMM_WORLD,2);return 2;}}
