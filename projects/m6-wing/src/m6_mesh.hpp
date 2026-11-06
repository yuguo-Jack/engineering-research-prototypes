#pragma once
#include "m6_physics.hpp"
#include <mpi.h>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace m6 {
struct MPIWorld {
 int rank=0,size=1;
 MPIWorld(int&argc,char**&argv){int provided=0;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);}
 ~MPIWorld(){MPI_Finalize();}
 double sum(double x)const{double y;MPI_Allreduce(&x,&y,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);return y;}
 double max(double x)const{double y;MPI_Allreduce(&x,&y,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);return y;}
 double min(double x)const{double y;MPI_Allreduce(&x,&y,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);return y;}
 template<size_t N>std::array<double,N>sum(const std::array<double,N>&x)const{std::array<double,N>y;MPI_Allreduce(x.data(),y.data(),static_cast<int>(N),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);return y;}
};
struct Peer{int rank;std::vector<int>send,recv;mutable std::vector<double>sb,rb;};
struct Face {int l,r,type,id,lr=-1,rl=-1;D3 center,normal;double area,dn;};
struct Mesh {
 int n=0,nt=0;uint64_t global_cells=0;std::array<double,48>meta{};
 std::vector<int>global;std::vector<double>volume,distance;std::vector<D3>center;
 std::vector<W7>initial;std::vector<Face>faces;std::vector<Peer>peers;
 template<class T>static void read(std::ifstream&f,T*ptr,size_t n){f.read(reinterpret_cast<char*>(ptr),static_cast<std::streamsize>(sizeof(T)*n));if(!f)throw std::runtime_error("truncated M6 shard");}
 Mesh(const std::string&directory,const MPIWorld&world){
  std::ostringstream filename;filename<<"rank-"<<std::setfill('0')<<std::setw(3)<<world.rank<<".bin";
  std::ifstream f(std::filesystem::path(directory)/filename.str(),std::ios::binary);if(!f)throw std::runtime_error("cannot open M6 shard "+filename.str());
  char magic[8];read(f,magic,8);if(std::string(magic,8)!="CFDM6A01")throw std::runtime_error("unsupported M6 shard version");
  uint32_t counts[4];read(f,counts,4);read(f,&global_cells,1);n=static_cast<int>(counts[0]);nt=static_cast<int>(counts[1]);int nf=static_cast<int>(counts[2]);
  read(f,meta.data(),meta.size());global.resize(nt);volume.resize(nt);distance.resize(nt);center.resize(nt);initial.resize(nt);
  read(f,global.data(),nt);read(f,volume.data(),nt);read(f,center.data(),nt);read(f,distance.data(),nt);read(f,initial.data(),nt);
  std::vector<int>l(nf),r(nf),t(nf),id(nf);std::vector<D3>fc(nf),av(nf);
  read(f,l.data(),nf);read(f,r.data(),nf);read(f,t.data(),nf);read(f,id.data(),nf);read(f,fc.data(),nf);read(f,av.data(),nf);
  faces.resize(nf);for(int j=0;j<nf;++j){double ar=length3(av[j]);D3 normal={av[j][0]/ar,av[j][1]/ar,av[j][2]/ar};D3 dx=r[j]>=0?sub3(center[r[j]],center[l[j]]):sub3(fc[j],center[l[j]]);double dn=dot3(dx,normal);
   if(l[j]<0||l[j]>=n||r[j]>=nt||!(ar>0)||!(dn>0))throw std::runtime_error("invalid original face geometry/ownership");
   faces[j]={l[j],r[j],t[j],id[j],-1,-1,fc[j],normal,ar,dn};
  }
  peers.resize(counts[3]);for(auto&p:peers){uint32_t a[3];read(f,a,3);p.rank=static_cast<int>(a[0]);p.send.resize(a[1]);p.recv.resize(a[2]);read(f,p.send.data(),p.send.size());read(f,p.recv.data(),p.recv.size());if(p.rank>=world.size)throw std::runtime_error("MPI rank count differs from prepared partition");}
  if(world.sum(static_cast<double>(n))!=static_cast<double>(global_cells))throw std::runtime_error("partition ownership does not cover original mesh");
  Gas gas(meta);for(auto&w:initial){double rho=w[0];w[0]=rho*meta[3]*w[4]/(meta[0]*meta[1]*meta[1]);for(int k=1;k<=3;++k)w[k]/=meta[1];w[4]/=meta[2];w[5]/=(meta[1]*meta[1]);w[6]/=meta[1];}
 }
 template<size_t N>void exchange(std::vector<std::array<double,N>>&values)const{
  if(peers.empty())return;if(values.size()<static_cast<size_t>(nt))values.resize(nt);
  std::vector<MPI_Request>requests;requests.reserve(2*peers.size());
  for(auto&p:peers){p.sb.resize(p.send.size()*N);p.rb.resize(p.recv.size()*N);for(size_t j=0;j<p.send.size();++j)std::copy(values[p.send[j]].begin(),values[p.send[j]].end(),p.sb.begin()+j*N);
   MPI_Request rr,rs;MPI_Irecv(p.rb.data(),static_cast<int>(p.rb.size()),MPI_DOUBLE,p.rank,110+static_cast<int>(N),MPI_COMM_WORLD,&rr);MPI_Isend(p.sb.data(),static_cast<int>(p.sb.size()),MPI_DOUBLE,p.rank,110+static_cast<int>(N),MPI_COMM_WORLD,&rs);requests.push_back(rr);requests.push_back(rs);
  }
  MPI_Waitall(static_cast<int>(requests.size()),requests.data(),MPI_STATUSES_IGNORE);
  for(auto&p:peers)for(size_t j=0;j<p.recv.size();++j)std::copy(p.rb.begin()+j*N,p.rb.begin()+(j+1)*N,values[p.recv[j]].begin());
 }
};
}
