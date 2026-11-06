#pragma once
#include "m6_mesh.hpp"

namespace m6 {
struct ScalarGaussReconstruction {
 std::vector<D3>raw,limited,physical;
 std::vector<std::array<double,1>>limiter;
};

// Positive turbulent scalars: harmonic interior GG reconstruction, type-aware
// boundary extrema, and arithmetic reconstructed-face physical derivatives.
// Boundary values and allocation are explicit inputs; their lifecycle belongs
// to the caller. This routine never changes any cell or boundary value.
inline ScalarGaussReconstruction scalar_gauss_reconstruction(
 const Mesh&mesh,const std::vector<double>&phi,
 const std::vector<double>&boundary,const std::vector<unsigned char>&allocated,
 const std::vector<std::array<double,1>>*frozen_limiter=nullptr){
 if(phi.size()!=static_cast<size_t>(mesh.nt)||boundary.size()!=mesh.faces.size()||allocated.size()!=mesh.faces.size())throw std::runtime_error("scalar reconstruction input dimensions differ from mesh");
 ScalarGaussReconstruction result;result.raw.resize(mesh.nt);result.limited.resize(mesh.nt);result.physical.resize(mesh.nt);result.limiter.resize(mesh.nt);
 std::vector<double>low=phi,high=phi;
 for(double v:phi)if(!(v>0&&std::isfinite(v)))throw std::runtime_error("scalar reconstruction requires finite positive cells");
 for(size_t fi=0;fi<mesh.faces.size();++fi){auto&f=mesh.faces[fi];int l=f.l,r=f.r;double face;
  if(r>=0){face=2*phi[l]*phi[r]/(phi[l]+phi[r]);low[l]=std::min(low[l],phi[r]);high[l]=std::max(high[l],phi[r]);if(r<mesh.n){low[r]=std::min(low[r],phi[l]);high[r]=std::max(high[r],phi[l]);}}
  else{face=allocated[fi]?boundary[fi]:phi[l];if(!(face>0&&std::isfinite(face)))throw std::runtime_error("scalar boundary cache must be finite and positive");double bound=f.type==3?2*face-phi[l]:face;low[l]=std::min(low[l],bound);high[l]=std::max(high[l],bound);}
  for(int j=0;j<3;++j){double a=f.area*f.normal[j];result.raw[l][j]+=(face-phi[l])*a;if(r>=0&&r<mesh.n)result.raw[r][j]-=(face-phi[r])*a;}
 }
 for(int i=0;i<mesh.n;++i){result.limiter[i][0]=1;for(double&x:result.raw[i])x/=mesh.volume[i];}
 for(auto&f:mesh.faces)for(int side=0;side<2;++side){int i=side?f.r:f.l;if(i<0||i>=mesh.n)continue;double delta=dot3(result.raw[i],sub3(f.center,mesh.center[i])),d=std::abs(delta),bound=delta>0?high[i]-phi[i]:phi[i]-low[i];if(bound>=2*d)continue;double t=bound*bound+4*d*d;double factor=bound*t/(bound*bound*bound+d*t+1e-100);result.limiter[i][0]=std::min(result.limiter[i][0],factor);}
 if(frozen_limiter){if(frozen_limiter->size()!=static_cast<size_t>(mesh.nt))throw std::runtime_error("frozen scalar limiter size mismatch");for(int i=0;i<mesh.n;++i)result.limiter[i]=(*frozen_limiter)[i];}
 for(int i=0;i<mesh.n;++i)for(int j=0;j<3;++j)result.limited[i][j]=result.limiter[i][0]*result.raw[i][j];mesh.exchange(result.limited);mesh.exchange(result.raw);mesh.exchange(result.limiter);
 for(size_t fi=0;fi<mesh.faces.size();++fi){auto&f=mesh.faces[fi];int l=f.l,r=f.r;double left=phi[l]+dot3(result.limited[l],sub3(f.center,mesh.center[l]));double face=r>=0?.5*(left+phi[r]+dot3(result.limited[r],sub3(f.center,mesh.center[r]))):allocated[fi]?boundary[fi]:left;
  for(int j=0;j<3;++j){double a=f.area*f.normal[j];result.physical[l][j]+=(face-phi[l])*a;if(r>=0&&r<mesh.n)result.physical[r][j]-=(face-phi[r])*a;}
 }
 for(int i=0;i<mesh.n;++i)for(double&x:result.physical[i])x/=mesh.volume[i];mesh.exchange(result.physical);
 return result;
}
}
