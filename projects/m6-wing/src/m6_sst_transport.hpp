#pragma once
#include "m6_scalar_reconstruction.hpp"

namespace m6 {
// Boundary scalar values are part of the evolving solver state. Initial values
// come solely from the original input data and case-derived initialization;
// neither reference diagnostics nor Fluent libraries are admitted here.
struct SSTTransportCache {
 const Mesh&mesh;const Gas&gas;const MPIWorld&world;
 std::array<std::vector<double>,2>boundary,reference;
 std::array<std::vector<unsigned char>,2>allocated;
 std::vector<size_t>faces;uint64_t global_faces=0;
 SSTTransportCache(const Mesh&m,const Gas&g,const MPIWorld&w,const std::string&directory):mesh(m),gas(g),world(w){
  for(int k=0;k<2;++k){boundary[k].resize(mesh.faces.size());allocated[k].resize(mesh.faces.size());reference[k].resize(mesh.nt);}
  for(size_t fi=0;fi<mesh.faces.size();++fi)if(mesh.faces[fi].r<0&&(mesh.faces[fi].type==3||mesh.faces[fi].type==9))faces.push_back(fi);
  global_faces=static_cast<uint64_t>(world.sum(static_cast<double>(faces.size())));
  std::ostringstream name;name<<"scalar-boundary-init-rank-"<<std::setfill('0')<<std::setw(3)<<world.rank<<".bin";
  std::ifstream f(std::filesystem::path(directory)/name.str(),std::ios::binary);if(!f)throw std::runtime_error("normal SST requires the original-input scalar boundary sidecar");char magic[8];int32_t h[6];Mesh::read(f,magic,8);Mesh::read(f,h,6);
  if(std::string(magic,8)!="M6SBI001"||h[0]!=1||h[1]!=faces.size()||h[2]!=world.rank||h[3]!=world.size||h[4]!=mesh.global_cells||h[5]!=global_faces)throw std::runtime_error("original scalar boundary input identity mismatch");
  std::vector<int32_t>ids(faces.size()),masks(faces.size());std::vector<std::array<double,2>>v(faces.size());Mesh::read(f,ids.data(),ids.size());Mesh::read(f,masks.data(),masks.size());Mesh::read(f,v.data(),v.size());if(f.peek()!=std::char_traits<char>::eof())throw std::runtime_error("trailing original scalar boundary input bytes");
  for(size_t j=0;j<faces.size();++j){size_t fi=faces[j];if(ids[j]!=mesh.faces[fi].id||masks[j]<0||masks[j]>3)throw std::runtime_error("original scalar boundary face map mismatch");for(int k=0;k<2;++k){double value=v[j][k]/(k==0?gas.a0*gas.a0:gas.a0);if(!(value>=0&&std::isfinite(value)))throw std::runtime_error("invalid original scalar boundary value");boundary[k][fi]=value;allocated[k][fi]=(masks[j]&(1<<k))?1:0;}}
 }
 ScalarGaussReconstruction reconstruct(int slot,const std::vector<W7>&w,const std::vector<W7>*frozen=nullptr){
  for(int i=0;i<mesh.nt;++i)reference[slot][i]=w[i][5+slot];std::vector<std::array<double,1>>limiter;
  if(frozen){limiter.resize(mesh.nt);for(int i=0;i<mesh.nt;++i)limiter[i][0]=(*frozen)[i][5+slot];}
  return scalar_gauss_reconstruction(mesh,reference[slot],boundary[slot],allocated[slot],frozen?&limiter:nullptr);
 }
 void refresh_farfield(int slot,const std::vector<W7>&w,const W7&freestream){for(size_t fi:faces){auto&f=mesh.faces[fi];if(f.type!=9)continue;double un=0;for(int j=0;j<3;++j)un+=w[f.l][j+1]*f.normal[j];boundary[slot][fi]=un>=0?w[f.l][5+slot]:freestream[5+slot];allocated[slot][fi]=1;}}
 void refresh_wall_k(const std::vector<W7>&w){for(size_t fi:faces){auto&f=mesh.faces[fi];if(f.type==3){boundary[0][fi]=w[f.l][5];allocated[0][fi]=1;}}}
 void refresh_wall_omega(const std::vector<W7>&w){for(size_t fi:faces){auto&f=mesh.faces[fi];if(f.type!=3)continue;const auto&v=w[f.l];D3 velocity={v[1],v[2],v[3]};double un=dot3(velocity,f.normal);for(int j=0;j<3;++j)velocity[j]-=un*f.normal[j];double rho=gas.rho(v);auto model=wall_model(rho,gas.mu,v[5],v[6],length3(velocity),f.dn);boundary[1][fi]=smooth_wall_omega_face(rho,gas.mu,model.ustar,f.dn,v[6]);allocated[1][fi]=1;}}
 std::vector<std::array<double,3>>snapshot()const{std::vector<std::array<double,3>>values(faces.size());for(size_t j=0;j<faces.size();++j){size_t fi=faces[j];values[j]={boundary[0][fi],boundary[1][fi],static_cast<double>((allocated[0][fi]?1:0)|(allocated[1][fi]?2:0))};}return values;}
 void restore(const std::vector<std::array<double,3>>&values){if(values.size()!=faces.size())throw std::runtime_error("incomplete saved scalar boundary state");for(size_t j=0;j<faces.size();++j){size_t fi=faces[j];int mask=static_cast<int>(values[j][2]);if(mask<0||mask>3||values[j][2]!=mask)throw std::runtime_error("invalid scalar boundary allocation mask");for(int k=0;k<2;++k){if(!(values[j][k]>=0&&std::isfinite(values[j][k])))throw std::runtime_error("invalid saved scalar boundary value");boundary[k][fi]=values[j][k];allocated[k][fi]=(mask&(1<<k))?1:0;}}}
 void save(const std::string&directory,int iteration,const std::vector<std::array<double,3>>&values)const{
  if(values.size()!=faces.size())throw std::runtime_error("incomplete scalar boundary checkpoint");std::ostringstream name;name<<"scalar-boundary-rank-"<<std::setfill('0')<<std::setw(3)<<world.rank<<".bin";std::ofstream out(std::filesystem::path(directory)/name.str(),std::ios::binary);out.write("M6SBP001",8);int32_t h[4]={static_cast<int32_t>(faces.size()),iteration,world.rank,world.size};out.write(reinterpret_cast<const char*>(h),sizeof(h));std::vector<int32_t>ids(faces.size()),masks(faces.size());std::vector<std::array<double,2>>v(faces.size());for(size_t j=0;j<faces.size();++j){ids[j]=mesh.faces[faces[j]].id;masks[j]=static_cast<int32_t>(values[j][2]);v[j]={values[j][0],values[j][1]};}out.write(reinterpret_cast<const char*>(ids.data()),ids.size()*sizeof(int32_t));out.write(reinterpret_cast<const char*>(masks.data()),masks.size()*sizeof(int32_t));out.write(reinterpret_cast<const char*>(v.data()),v.size()*sizeof(v[0]));out.flush();if(!out)throw std::runtime_error("scalar boundary checkpoint write failed");
 }
 void load(const std::string&directory,int iteration){std::ostringstream name;name<<"scalar-boundary-rank-"<<std::setfill('0')<<std::setw(3)<<world.rank<<".bin";std::ifstream input(std::filesystem::path(directory)/name.str(),std::ios::binary);if(!input)throw std::runtime_error("normal SST restart requires its saved scalar boundary cache");char magic[8];int32_t h[4];Mesh::read(input,magic,8);Mesh::read(input,h,4);if(std::string(magic,8)!="M6SBP001"||h[0]!=faces.size()||h[1]!=iteration||h[2]!=world.rank||h[3]!=world.size)throw std::runtime_error("scalar boundary restart identity mismatch");std::vector<int32_t>ids(faces.size()),masks(faces.size());std::vector<std::array<double,2>>v(faces.size());Mesh::read(input,ids.data(),ids.size());Mesh::read(input,masks.data(),masks.size());Mesh::read(input,v.data(),v.size());if(input.peek()!=std::char_traits<char>::eof())throw std::runtime_error("trailing scalar boundary restart bytes");std::vector<std::array<double,3>>values(faces.size());for(size_t j=0;j<faces.size();++j){if(ids[j]!=mesh.faces[faces[j]].id)throw std::runtime_error("scalar boundary restart face ordering mismatch");values[j]={v[j][0],v[j][1],static_cast<double>(masks[j])};}restore(values);}
};
}
