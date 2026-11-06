#pragma once
#include "m6_mesh.hpp"
namespace m6 {
// Short-history Anderson Type-II, reorthogonalized QR. MPI owns disjoint rows.
// The caller evaluates the actual discretized equations before accepting it.
class Anderson7 {
 const MPIWorld&world;int depth;bool valid=false;
 using F=std::vector<std::array<double,7>>;
 F oldy,oldstep;std::vector<F>df,dg;
 double dot(const F&a,const F&b)const{double v=0;for(size_t i=0;i<a.size();++i)for(int k=0;k<7;++k)v+=a[i][k]*b[i][k];return world.sum(v);}
 static void axpy(F&a,double t,const F&b){for(size_t i=0;i<a.size();++i)for(int k=0;k<7;++k)a[i][k]+=t*b[i][k];}
 public:
 Anderson7(const MPIWorld&w,int m):world(w),depth(m){}
 void clear(){valid=false;oldy.clear();oldstep.clear();df.clear();dg.clear();}
 F propose(const F&y,const F&image){
  F step=image;axpy(step,-1,y);
  if(valid){F f=step,g=image;axpy(f,-1,oldstep);axpy(g,-1,oldy);axpy(g,-1,oldstep);df.push_back(std::move(f));dg.push_back(std::move(g));if(static_cast<int>(df.size())>depth){df.erase(df.begin());dg.erase(dg.begin());}}
  oldy=y;oldstep=step;valid=true;F answer=image;
  int m=static_cast<int>(df.size());std::vector<F>basis;std::vector<int>active;std::vector<std::vector<double>>r(m+1,std::vector<double>(m+1));std::vector<double>b;
  for(int j=0;j<m;++j){F v=df[j];int n=static_cast<int>(basis.size());for(int i=0;i<=n;++i)r[i][n]=0;
   for(int pass=0;pass<2;++pass)for(int i=0;i<n;++i){double h=dot(basis[i],v);r[i][n]+=h;axpy(v,-h,basis[i]);}
   double length=std::sqrt(dot(v,v));if(length<1e-11*std::max(std::sqrt(dot(df[j],df[j])),1e-30))continue;r[n][n]=length;for(auto&item:v)for(double&t:item)t/=length;basis.push_back(std::move(v));active.push_back(j);b.push_back(dot(basis.back(),step));
  }
  int n=static_cast<int>(basis.size());std::vector<double>gamma(n);
  for(int i=n-1;i>=0;--i){double v=b[i];for(int j=i+1;j<n;++j)v-=r[i][j]*gamma[j];gamma[i]=v/r[i][i];if(!std::isfinite(gamma[i])||std::abs(gamma[i])>1e3){df.clear();dg.clear();return image;}}
  for(int i=0;i<n;++i)axpy(answer,-gamma[i],dg[active[i]]);return answer;
 }
};
}
