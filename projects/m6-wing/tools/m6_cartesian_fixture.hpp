#pragma once
// Original27-cell synthetic input, reused without formula/data changes.
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
