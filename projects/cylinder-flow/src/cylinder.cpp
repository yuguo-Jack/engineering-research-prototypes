// Independent finite-volume Euler prototype. No Fluent binary/library is linked.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#include "coupled_flux.hpp"
#ifndef CFD_CYLINDER_SKIP_UNUSED_RADIUS_WORK
#define CFD_CYLINDER_SKIP_UNUSED_RADIUS_WORK 0
#endif
#ifndef CFD_CYLINDER_INCREMENTAL_ANDERSON_QR
#define CFD_CYLINDER_INCREMENTAL_ANDERSON_QR 0
#endif

using V=std::array<double,4>;
using M=std::array<double,16>;
using VV=std::vector<V>;
V operator+(V a,const V&b){for(int k=0;k<4;++k)a[k]+=b[k];return a;}
V operator-(V a,const V&b){for(int k=0;k<4;++k)a[k]-=b[k];return a;}
V operator*(double s,V a){for(double&v:a)v*=s;return a;}
V& operator+=(V&a,const V&b){for(int k=0;k<4;++k)a[k]+=b[k];return a;}
V& operator-=(V&a,const V&b){for(int k=0;k<4;++k)a[k]-=b[k];return a;}
V mv(const M&a,const V&b){V r{};for(int i=0;i<4;++i)for(int j=0;j<4;++j)r[i]+=a[4*i+j]*b[j];return r;}
M mm(const M&a,const M&b){M r{};for(int i=0;i<4;++i)for(int k=0;k<4;++k)for(int j=0;j<4;++j)r[4*i+j]+=a[4*i+k]*b[4*k+j];return r;}
M inv(M a){
 M b{};for(int i=0;i<4;++i)b[5*i]=1;
 for(int k=0;k<4;++k){int p=k;for(int i=k+1;i<4;++i)if(std::abs(a[4*i+k])>std::abs(a[4*p+k]))p=i;
  if(std::abs(a[4*p+k])<1e-20)throw std::runtime_error("singular 4x4 block");
  for(int j=0;j<4;++j){std::swap(a[4*k+j],a[4*p+j]);std::swap(b[4*k+j],b[4*p+j]);}
  const double d=a[5*k];for(int j=0;j<4;++j){a[4*k+j]/=d;b[4*k+j]/=d;}
  for(int i=0;i<4;++i)if(i!=k){const double f=a[4*i+k];for(int j=0;j<4;++j){a[4*i+j]-=f*a[4*k+j];b[4*i+j]-=f*b[4*k+j];}}
 }return b;
}
double dot(const VV&a,const VV&b){double s=0;for(size_t i=0;i<a.size();++i)for(int k=0;k<4;++k)s+=a[i][k]*b[i][k];return s;}
double norm(const VV&a){return std::sqrt(dot(a,a));}
void axpy(VV&a,double s,const VV&b){for(size_t i=0;i<a.size();++i)for(int k=0;k<4;++k)a[i][k]+=s*b[i][k];}

#if CFD_CYLINDER_INCREMENTAL_ANDERSON_QR
struct Anderson {
 int depth;bool valid=false;VV oldq,oldstep;std::vector<VV>df,dg,basis;
 std::vector<double>dfnorms;std::vector<std::vector<double>>r;std::vector<int>active;
 bool qr_valid=true;int qr_age=0;
 int full_rebuilds=0,incremental_columns=0,downdates=0,rank_fallbacks=0;
 explicit Anderson(int m):depth(m){}
 void rebuild(){
  int m=static_cast<int>(df.size());basis.clear();active.clear();r.assign(m,std::vector<double>(m));
  for(int j=0;j<m;++j){VV v=df[j];int n=static_cast<int>(basis.size());for(int i=0;i<=n;++i)r[i][n]=0;
   for(int pass=0;pass<2;++pass)for(int i=0;i<n;++i){double h=dot(basis[i],v);r[i][n]+=h;axpy(v,-h,basis[i]);}
   double length=norm(v);if(length<1e-10*std::max(dfnorms[j],1e-30))continue;
   r[n][n]=length;for(auto&item:v)for(double&x:item)x/=length;
   basis.push_back(std::move(v));active.push_back(j);
  }
  qr_valid=basis.size()==df.size();qr_age=0;++full_rebuilds;
 }
 bool drop_first(){
  int n=static_cast<int>(basis.size());if(n==0)return false;
  std::vector<std::vector<double>>next(n,std::vector<double>(n-1));
  for(int i=0;i<n;++i)for(int j=0;j<n-1;++j)next[i][j]=r[i][j+1];
  for(int i=0;i<n-1;++i){double a=next[i][i],b=next[i+1][i],length=std::hypot(a,b);
   if(!(length>0&&std::isfinite(length)))return false;
   double c=a/length,s=b/length;
   for(int j=i;j<n-1;++j){double top=next[i][j],bottom=next[i+1][j];
    next[i][j]=c*top+s*bottom;next[i+1][j]=-s*top+c*bottom;
   }next[i+1][i]=0;
   for(size_t cell=0;cell<basis[i].size();++cell)for(int k=0;k<4;++k){
    double left=basis[i][cell][k],right=basis[i+1][cell][k];
    basis[i][cell][k]=c*left+s*right;basis[i+1][cell][k]=-s*left+c*right;
   }
  }
  basis.pop_back();r.assign(n-1,std::vector<double>(n-1));
  for(int i=0;i<n-1;++i)for(int j=i;j<n-1;++j)r[i][j]=next[i][j];
  for(int j=0;j<n-1;++j)if(std::abs(r[j][j])<1e-10*std::max(dfnorms[j],1e-30))return false;
  ++downdates;return true;
 }
 bool append_latest(){
  int n=static_cast<int>(basis.size());VV v=df.back();
  for(auto&row:r)row.resize(n+1,0);r.resize(n+1,std::vector<double>(n+1));
  for(int pass=0;pass<2;++pass)for(int i=0;i<n;++i){double h=dot(basis[i],v);r[i][n]+=h;axpy(v,-h,basis[i]);}
  double length=norm(v);if(length<1e-10*std::max(dfnorms.back(),1e-30)||!std::isfinite(length))return false;
  r[n][n]=length;for(auto&item:v)for(double&x:item)x/=length;basis.push_back(std::move(v));
  active.resize(basis.size());std::iota(active.begin(),active.end(),0);++incremental_columns;++qr_age;return true;
 }
 VV propose(const VV&q,const VV&step){
  const size_t old_columns=df.size();bool dropped=false;
  if(valid){VV f=step,g=q;axpy(f,-1,oldstep);axpy(g,1,step);axpy(g,-1,oldq);axpy(g,-1,oldstep);
   dfnorms.push_back(norm(f));df.push_back(std::move(f));dg.push_back(std::move(g));
   if(static_cast<int>(df.size())>depth){df.erase(df.begin());dg.erase(dg.begin());dfnorms.erase(dfnorms.begin());dropped=true;}
  }
  oldq=q;oldstep=step;valid=true;VV answer=step;
  bool incremental=qr_valid&&basis.size()==old_columns&&qr_age<8&&!df.empty();
  if(incremental&&dropped)incremental=drop_first();
  if(incremental)incremental=append_latest();
  if(!incremental){if(!df.empty())++rank_fallbacks;rebuild();}
  int n=static_cast<int>(basis.size());std::vector<double>b(n),gamma(n);
  for(int i=0;i<n;++i)b[i]=dot(basis[i],step);
  for(int i=n-1;i>=0;--i){double value=b[i];for(int j=i+1;j<n;++j)value-=r[i][j]*gamma[j];gamma[i]=value/r[i][i];
   if(!std::isfinite(gamma[i])||std::abs(gamma[i])>1e4){restart();return step;}
  }
  for(int i=0;i<n;++i)axpy(answer,-gamma[i],dg[active[i]]);return answer;
 }
 void restart(){df.clear();dg.clear();dfnorms.clear();basis.clear();r.clear();active.clear();qr_valid=true;qr_age=0;}
};
#else
struct Anderson {
 int depth;bool valid=false;VV oldq,oldstep;std::vector<VV>df,dg;
 explicit Anderson(int m):depth(m){}
 VV propose(const VV&q,const VV&step){
  if(valid){VV f=step,g=q;axpy(f,-1,oldstep);axpy(g,1,step);axpy(g,-1,oldq);axpy(g,-1,oldstep);df.push_back(std::move(f));dg.push_back(std::move(g));if(static_cast<int>(df.size())>depth){df.erase(df.begin());dg.erase(dg.begin());}}
  oldq=q;oldstep=step;valid=true;VV answer=step;
  // Reorthogonalized thin QR; dependent columns are omitted instead of solving
  // normal equations. The small history is rebuilt for transparent diagnostics.
  int m=static_cast<int>(df.size());std::vector<VV>basis;std::vector<int>active;std::vector<std::vector<double>>r(m,std::vector<double>(m));std::vector<double>b;
  for(int j=0;j<m;++j){VV v=df[j];int n=static_cast<int>(basis.size());for(int i=0;i<=n;++i)r[i][n]=0;for(int pass=0;pass<2;++pass)for(int i=0;i<n;++i){double h=dot(basis[i],v);r[i][n]+=h;axpy(v,-h,basis[i]);}double length=norm(v);if(length<1e-10*std::max(norm(df[j]),1e-30))continue;r[n][n]=length;for(auto&item:v)for(double&x:item)x/=length;basis.push_back(std::move(v));active.push_back(j);b.push_back(dot(basis.back(),step));}
  int n=static_cast<int>(basis.size());std::vector<double>gamma(n);for(int i=n-1;i>=0;--i){double v=b[i];for(int j=i+1;j<n;++j)v-=r[i][j]*gamma[j];gamma[i]=v/r[i][i];if(!std::isfinite(gamma[i])||std::abs(gamma[i])>1e4){df.clear();dg.clear();return step;}}
  for(int i=0;i<n;++i)axpy(answer,-gamma[i],dg[active[i]]);return answer;
 }
 void restart(){df.clear();dg.clear();}
};
#endif

struct Gas {
 double g=1.4, mach=.1;bool precondition=true;
 V conserved(const V&w)const{return {w[0],w[0]*w[1],w[0]*w[2],w[3]/(g-1)+.5*w[0]*(w[1]*w[1]+w[2]*w[2])};}
 V primitive(const V&q)const{const double u=q[1]/q[0],v=q[2]/q[0];return {q[0],u,v,(g-1)*(q[3]-.5*q[0]*(u*u+v*v))};}
 V physical(const V&w,double nx,double ny)const{double un=w[1]*nx+w[2]*ny;V q=conserved(w);return {q[0]*un,q[1]*un+w[3]*nx,q[2]*un+w[3]*ny,(q[3]+w[3])*un};}
 struct LinearState{
  double rho,u,v,p,T,H,cp,theta,un,beta2,sound;
  V ci(const V&z)const{
   double du=(z[1]-u*z[0])/rho,dv=(z[2]-v*z[0])/rho;
   double gamma=cp/(cp-1/(T*rho/p));
   double dp=(gamma-1)*(z[3]-.5*(u*u+v*v)*z[0]-rho*(u*du+v*dv));
   return {dp,du,dv,T*(dp/p-z[0]/rho)};
  }
  V gamma_apply(const V&z)const{
   double dm=theta*z[0]-rho/T*z[3];
   return {dm,u*dm+rho*z[1],v*dm+rho*z[2],H*dm-z[0]+rho*(u*z[1]+v*z[2]+cp*z[3])};
  }
  V gamma_inverse(const V&z)const{
   double du=(z[1]-u*z[0])/rho,dv=(z[2]-v*z[0])/rho;
   double dp=(z[3]-.5*(u*u+v*v)*z[0]-rho*(u*du+v*dv))/(cp*T*theta-1);
   return {dp,du,dv,T/rho*(theta*dp-z[0])};
  }
  V a_apply(const V&z,double nx,double ny)const{
   double dr=rho/p*z[0]-rho/T*z[3],dn=nx*z[1]+ny*z[2];
   double dm=un*dr+rho*dn;
   return {dm,u*dm+rho*un*z[1]+nx*z[0],v*dm+rho*un*z[2]+ny*z[0],H*dm+rho*un*(cp*z[3]+u*z[1]+v*z[2])};
  }
 };
 LinearState state(double rho,double u,double v,double H,double nx,double ny,double explicit_uref2=-1)const{
  double cp=1/(g-1),T=(H-.5*(u*u+v*v))/cp;
  double a2=T,beta2=precondition?(explicit_uref2>=0?std::min(1.,std::max(1e-20,explicit_uref2/a2)):std::min(1.,std::max(mach*mach,(u*u+v*v)/a2))):1.;
  double theta=1/(beta2*a2)+1/(cp*T);
  return {rho,u,v,rho*T/g,T,H,cp,theta,u*nx+v*ny,beta2,std::sqrt(a2)};
 }
 LinearState state(const V&w,double nx,double ny,double explicit_uref2=-1)const{
  double H=g/(g-1)*w[3]/w[0]+.5*(w[1]*w[1]+w[2]*w[2]);
  return state(w[0],w[1],w[2],H,nx,ny,explicit_uref2);
 }
 double radius(const V&w,double nx,double ny,double explicit_uref2=-1)const{
  auto s=state(w,nx,ny,explicit_uref2);double up=.5*(1+s.beta2)*s.un;
  double c=std::sqrt(.25*std::pow((1-s.beta2)*s.un,2)+s.beta2*s.sound*s.sound);
  return std::abs(up)+c;
 }
 V roe(const V&l,const V&r,double nx,double ny,double explicit_uref2=-1)const{
  double sl=std::sqrt(l[0]),sr=std::sqrt(r[0]),den=sl+sr;
  double hl=g/(g-1)*l[3]/l[0]+.5*(l[1]*l[1]+l[2]*l[2]);
  double hr=g/(g-1)*r[3]/r[0]+.5*(r[1]*r[1]+r[2]*r[2]);
  auto s=state(sl*sr,(sl*l[1]+sr*r[1])/den,(sl*l[2]+sr*r[2])/den,(sl*hl+sr*hr)/den,nx,ny,explicit_uref2);
  V dw=s.ci(conserved(r)-conserved(l));
  auto k_apply=[&](const V&z){return s.gamma_inverse(s.a_apply(z,nx,ny))-s.un*z;};
  V k1=k_apply(dw),k2=k_apply(k1);
  double up=.5*(1+s.beta2)*s.un;
  double c=std::sqrt(.25*std::pow((1-s.beta2)*s.un,2)+s.beta2*s.sound*s.sound);
  double lp=up+c,lm=up-c,dp=lp-s.un,dm=lm-s.un,ab=std::abs(s.un);
  auto entropy=[&](double x){double d=.05*s.sound*std::sqrt(s.beta2);x=std::abs(x);return x<d?.5*(x*x/d+d):x;};
  double fp=entropy(lp)-ab,fm=entropy(lm)-ab;
  double a2=(fp/dp-fm/dm)/(dp-dm),a1=fp/dp-a2*dp;
  V diss=s.gamma_apply(ab*dw+a1*k1+a2*k2);
  return .5*(physical(l,nx,ny)+physical(r,nx,ny)-diss);
 }
 V coupled_roe(const V&l,const V&r,double nx,double ny,double uref_l,double uref_r,double cell_a2_l,double cell_a2_r)const{
  return coupled_roe_split(l,r,nx,ny,g,uref_l,uref_r,cell_a2_l,cell_a2_r);
 }
 V farfield(const V&w,const V&inf,double nx,double ny)const{
  double un=w[1]*nx+w[2]*ny,ui=inf[1]*nx+inf[2]*ny;
  double a=std::sqrt(g*w[3]/w[0]),ai=std::sqrt(g*inf[3]/inf[0]);
  if(un>=a)return w;if(un<=-a)return inf;
  double jplus=un+2*a/(g-1),jminus=ui-2*ai/(g-1);
  double ub=.5*(jplus+jminus),ab=.25*(g-1)*(jplus-jminus);
  const V&up=ub>=0?w:inf;double entropy=up[3]/std::pow(up[0],g);
  double rb=std::pow(ab*ab/(g*entropy),1/(g-1));
  double normal=up[1]*nx+up[2]*ny;
  return {rb,up[1]+(ub-normal)*nx,up[2]+(ub-normal)*ny,rb*ab*ab/g};
 }
};

struct Face{int l,r,type;double x,y,nx,ny,area;int ll,lr,rl,rr;};
struct Mesh{
 int n=0;double R=0,cp=0,rho0=0,a0=0,pinf=0,Tinf=0;Gas gas;
 std::vector<double>volume;std::vector<std::array<double,2>>center;std::vector<Face>faces;VV initial;
 template<class T>static void read(std::ifstream&f,T*ptr,size_t n){f.read(reinterpret_cast<char*>(ptr),static_cast<std::streamsize>(n*sizeof(T)));if(!f)throw std::runtime_error("truncated input");}
 explicit Mesh(const std::string&path){
  std::ifstream f(path,std::ios::binary);if(!f)throw std::runtime_error("cannot open mesh input");
  char magic[8];read(f,magic,8);if(std::string(magic,8)!="CFDCYL01")throw std::runtime_error("wrong input magic");
  uint32_t nn,nf;read(f,&nn,1);read(f,&nf,1);n=static_cast<int>(nn);
  double meta[7];read(f,meta,7);gas.g=meta[0];R=meta[1];cp=meta[2];rho0=meta[3];a0=meta[4];pinf=meta[5];Tinf=meta[6];
  volume.resize(n);center.resize(n);read(f,volume.data(),n);read(f,center.data(),n);
  std::vector<int32_t>l(nf),r(nf),t(nf);std::vector<std::array<double,2>>fc(nf),av(nf);
  read(f,l.data(),nf);read(f,r.data(),nf);read(f,t.data(),nf);read(f,fc.data(),nf);read(f,av.data(),nf);
  faces.resize(nf);for(size_t i=0;i<nf;++i){double ar=std::hypot(av[i][0],av[i][1]);faces[i]={l[i],r[i],t[i],fc[i][0],fc[i][1],av[i][0]/ar,av[i][1]/ar,ar,-1,-1,-1,-1};}
  initial.resize(n);read(f,initial.data(),n);gas.mach=initial[0][1];
 }
};

struct Sparse{
 int n;std::vector<int>row,col,diag;std::vector<M>a,lu,dinv;
 explicit Sparse(Mesh&m,int fill_level=0):n(m.n),row(n+1),diag(n){
  std::vector<std::vector<int>>adj(n);for(int i=0;i<n;++i)adj[i].push_back(i);
  for(auto&f:m.faces)if(f.r>=0){adj[f.l].push_back(f.r);adj[f.r].push_back(f.l);}
  if(fill_level>0){std::vector<std::map<int,int>>level(n);for(int i=0;i<n;++i)for(int j:adj[i])level[i][j]=0;for(int i=0;i<n;++i){for(auto it=level[i].begin();it!=level[i].end()&&it->first<i;++it){int j=it->first;for(auto jt=level[j].upper_bound(j);jt!=level[j].end();++jt){int lev=it->second+jt->second+1;if(lev>fill_level)continue;auto old=level[i].find(jt->first);if(old==level[i].end())level[i][jt->first]=lev;else old->second=std::min(old->second,lev);}}adj[i].clear();for(const auto&entry:level[i])adj[i].push_back(entry.first);}}
  for(int i=0;i<n;++i){auto&v=adj[i];std::sort(v.begin(),v.end());v.erase(std::unique(v.begin(),v.end()),v.end());row[i]=static_cast<int>(col.size());for(int j:v)col.push_back(j);}
  row[n]=static_cast<int>(col.size());a.resize(col.size());lu.resize(col.size());dinv.resize(n);
  for(int i=0;i<n;++i)diag[i]=index(i,i);
  for(auto&f:m.faces){f.ll=diag[f.l];if(f.r>=0){f.lr=index(f.l,f.r);f.rl=index(f.r,f.l);f.rr=diag[f.r];}}
 }
 int index(int i,int j)const{auto it=std::lower_bound(col.begin()+row[i],col.begin()+row[i+1],j);return it!=col.begin()+row[i+1]&&*it==j?static_cast<int>(it-col.begin()):-1;}
 void clear(){std::fill(a.begin(),a.end(),M{});}
 void factor(){
  lu=a;for(int i=0;i<n;++i){for(int p=row[i];p<diag[i];++p){int j=col[p];M l=mm(lu[p],dinv[j]);lu[p]=l;
    for(int q=diag[j]+1;q<row[j+1];++q){int pos=index(i,col[q]);if(pos>=0){M prod=mm(l,lu[q]);for(int k=0;k<16;++k)lu[pos][k]-=prod[k];}}
   }dinv[i]=inv(lu[diag[i]]);
  }
 }
 VV multiply(const VV&x)const{VV y(n);for(int i=0;i<n;++i)for(int p=row[i];p<row[i+1];++p)y[i]+=mv(a[p],x[col[p]]);return y;}
 VV solve(const VV&b)const{VV x=b;for(int i=0;i<n;++i)for(int p=row[i];p<diag[i];++p)x[i]-=mv(lu[p],x[col[p]]);for(int i=n-1;i>=0;--i){for(int p=diag[i]+1;p<row[i+1];++p)x[i]-=mv(lu[p],x[col[p]]);x[i]=mv(dinv[i],x[i]);}return x;}
 int gmres(const VV&b,VV&x,int restart,int maxsteps,double tolerance)const{
  x.assign(n,V{});double bn=norm(b);if(bn<1e-30)return 0;int steps=0;
  while(steps<maxsteps){VV res=b;VV ax=multiply(x);axpy(res,-1,ax);double beta=norm(res);if(beta<=tolerance*bn)return steps;
   std::vector<VV>v,z;v.push_back(res);for(auto&item:v[0])for(double&d:item)d/=beta;
   std::vector<std::vector<double>>h(restart+1,std::vector<double>(restart));std::vector<double>cs(restart),sn(restart),s(restart+1);s[0]=beta;int used=0;
   for(int j=0;j<restart&&steps<maxsteps;++j,++steps){z.push_back(solve(v[j]));VV w=multiply(z[j]);
    for(int i=0;i<=j;++i){h[i][j]=dot(w,v[i]);axpy(w,-h[i][j],v[i]);}
    // Reorthogonalization protects the low-Mach block system from loss of basis orthogonality.
    for(int i=0;i<=j;++i){double corr=dot(w,v[i]);h[i][j]+=corr;axpy(w,-corr,v[i]);}
    h[j+1][j]=norm(w);if(h[j+1][j]>1e-25)for(auto&item:w)for(double&d:item)d/=h[j+1][j];v.push_back(std::move(w));
    for(int i=0;i<j;++i){double t=cs[i]*h[i][j]+sn[i]*h[i+1][j];h[i+1][j]=-sn[i]*h[i][j]+cs[i]*h[i+1][j];h[i][j]=t;}
    double den=std::hypot(h[j][j],h[j+1][j]);cs[j]=den>0?h[j][j]/den:1;sn[j]=den>0?h[j+1][j]/den:0;
    h[j][j]=cs[j]*h[j][j]+sn[j]*h[j+1][j];h[j+1][j]=0;s[j+1]=-sn[j]*s[j];s[j]*=cs[j];used=j+1;
    if(std::abs(s[j+1])<tolerance*bn){++steps;break;}
   }
   std::vector<double>y(used);for(int i=used-1;i>=0;--i){double sum=s[i];for(int j=i+1;j<used;++j)sum-=h[i][j]*y[j];y[i]=sum/h[i][i];}
   for(int i=0;i<used;++i)axpy(x,y[i],z[i]);
  }return steps;
 }
};

struct Solver{
 Mesh mesh;Sparse matrix;VV q,w,gx,gy;std::vector<std::array<double,3>>lsq,scalarlsq;std::vector<double>uref2;V inf;int order=2;std::string limiter="barth";bool wall_extrapolate=false,native_numerics=false,coupled_flux=false,analytic_jacobian=false;
 explicit Solver(const std::string&input,int fill_level=0):mesh(input),matrix(mesh,fill_level),q(mesh.n),w(mesh.n),gx(mesh.n),gy(mesh.n),lsq(mesh.n),scalarlsq(mesh.n){
  inf=mesh.initial[0];for(int i=0;i<mesh.n;++i)q[i]=mesh.gas.conserved(mesh.initial[i]);
  for(auto&f:mesh.faces){auto add=[&](int i,int j){double dx=j>=0?mesh.center[j][0]-mesh.center[i][0]:2*(f.x-mesh.center[i][0]);double dy=j>=0?mesh.center[j][1]-mesh.center[i][1]:2*(f.y-mesh.center[i][1]);double z=1/(dx*dx+dy*dy);lsq[i][0]+=dx*dx*z;lsq[i][1]+=dx*dy*z;lsq[i][2]+=dy*dy*z;if(j>=0||f.type!=3){scalarlsq[i][0]+=dx*dx*z;scalarlsq[i][1]+=dx*dy*z;scalarlsq[i][2]+=dy*dy*z;}};add(f.l,f.r);if(f.r>=0)add(f.r,f.l);}
  for(auto&v:lsq){double d=v[0]*v[2]-v[1]*v[1];v={v[2]/d,-v[1]/d,v[0]/d};}
  for(auto&v:scalarlsq){double d=v[0]*v[2]-v[1]*v[1];v={v[2]/d,-v[1]/d,v[0]/d};}
  uref2.resize(mesh.n);
 }
 void update_reference_velocity(){
  for(int i=0;i<mesh.n;++i)uref2[i]=w[i][1]*w[i][1]+w[i][2]*w[i][2];
  for(const auto&f:mesh.faces)if(f.r>=0){double jump=4*std::abs(w[f.l][3]-w[f.r][3])/(w[f.l][0]+w[f.r][0]);uref2[f.l]=std::max(uref2[f.l],jump);uref2[f.r]=std::max(uref2[f.r],jump);}
  for(int i=0;i<mesh.n;++i){double a2=mesh.gas.g*w[i][3]/w[i][0];double u2=uref2[i];uref2[i]=std::min(a2,std::max(1e-20*a2,u2/(.5-std::min(u2/a2,.25))));}
 }
 V boundary(const V&wl,const Face&f)const{if(f.type==3){V v=wl;double un=wl[1]*f.nx+wl[2]*f.ny;v[1]-=2*un*f.nx;v[2]-=2*un*f.ny;return v;}return mesh.gas.farfield(wl,inf,f.nx,f.ny);}
 void gradients(){
  for(int i=0;i<mesh.n;++i)w[i]=mesh.gas.primitive(q[i]);std::fill(gx.begin(),gx.end(),V{});std::fill(gy.begin(),gy.end(),V{});if(order==1)return;
  VV low=w,high=w;
  auto accumulate=[&](int i,int j,const Face&f){
   double dx=j>=0?mesh.center[j][0]-mesh.center[i][0]:2*(f.x-mesh.center[i][0]);double dy=j>=0?mesh.center[j][1]-mesh.center[i][1]:2*(f.y-mesh.center[i][1]);
   V wr=j>=0?w[j]:boundary(w[i],f);if(j<0&&f.type==9)wr=2*wr-w[i];
   double d=1/(dx*dx+dy*dy);for(int k=0;k<4;++k){double delta=wr[k]-w[i][k];gx[i][k]+=d*dx*delta;gy[i][k]+=d*dy*delta;if(j>=0||limiter!="standard"){low[i][k]=std::min(low[i][k],wr[k]);high[i][k]=std::max(high[i][k],wr[k]);}}
  };
  for(auto&f:mesh.faces){accumulate(f.l,f.r,f);if(f.r>=0)accumulate(f.r,f.l,f);}
  for(int i=0;i<mesh.n;++i)for(int k=0;k<4;++k){const auto&mat=wall_extrapolate&&(k==0||k==3)?scalarlsq[i]:lsq[i];double x=mat[0]*gx[i][k]+mat[1]*gy[i][k];double y=mat[1]*gx[i][k]+mat[2]*gy[i][k];gx[i][k]=x;gy[i][k]=y;}
  if(limiter=="none")return;
  if(limiter=="standard")for(const auto&f:mesh.faces)if(f.r<0){int i=f.l;V bc=boundary(w[i],f);if(f.type==3){double dn=(f.x-mesh.center[i][0])*f.nx+(f.y-mesh.center[i][1])*f.ny;bc=w[i];for(int k=0;k<4;++k)bc[k]+=2*dn*(gx[i][k]*f.nx+gy[i][k]*f.ny);double un=bc[1]*f.nx+bc[2]*f.ny;bc[1]-=un*f.nx;bc[2]-=un*f.ny;}for(int k=0;k<4;++k){low[i][k]=std::min(low[i][k],bc[k]);high[i][k]=std::max(high[i][k],bc[k]);}}
  VV phi(mesh.n,V{1,1,1,1});V global_low=w[0],global_high=w[0];for(const auto&v:w)for(int k=0;k<4;++k){global_low[k]=std::min(global_low[k],v[k]);global_high[k]=std::max(global_high[k],v[k]);}
  auto limit=[&](int i,const Face&f){double dx=f.x-mesh.center[i][0],dy=f.y-mesh.center[i][1];for(int k=0;k<4;++k){double dq=dx*gx[i][k]+dy*gy[i][k];if(std::abs(dq)<1e-16)continue;double bound=dq>0?high[i][k]-w[i][k]:low[i][k]-w[i][k];double ratio=bound/dq;
    if(limiter=="venkat"){double eps=std::pow(.1*std::sqrt(mesh.volume[i]),3);double d2=bound*bound;ratio=(d2+2*bound*dq+eps)/(d2+bound*dq+2*dq*dq+eps);}
    if(limiter=="differentiable"){double eps=std::pow(.05*(global_high[k]-global_low[k]),2);double d2=bound*bound;ratio=(d2+2*bound*dq+eps)/(d2+bound*dq+2*dq*dq+eps);}
    if(limiter=="standard"){double d=std::abs(dq),b=std::abs(bound),poly=b*b+4*d*d;ratio=b>=2*d?1:(b*poly+1e-32)/(b*b*b+d*poly+1e-32);}
    phi[i][k]=std::max(0.,std::min(phi[i][k],ratio));}};
  for(auto&f:mesh.faces){limit(f.l,f);if(f.r>=0)limit(f.r,f);}
  for(int i=0;i<mesh.n;++i)for(int k=0;k<4;++k){gx[i][k]*=phi[i][k];gy[i][k]*=phi[i][k];}
 }
 V reconstruct(int i,const Face&f,const V&u)const{V value=mesh.gas.primitive(u);double dx=f.x-mesh.center[i][0],dy=f.y-mesh.center[i][1];for(int k=0;k<4;++k)value[k]+=dx*gx[i][k]+dy*gy[i][k];return value;}
 V faceflux(const Face&f,const V&ul,const V&ur)const{
  V l=reconstruct(f.l,f,ul);if(l[0]<=0||l[3]<=0)throw std::runtime_error("nonpositive reconstructed left state");
  if(f.r<0){V right=boundary(l,f);if(f.type==3){if(native_numerics){double dp=l[0]*std::sqrt(uref2[f.l])*(l[1]*f.nx+l[2]*f.ny);return {0,(l[3]+dp)*f.nx,(l[3]+dp)*f.ny,0};}return mesh.gas.roe(l,right,f.nx,f.ny);}return mesh.gas.physical(right,f.nx,f.ny);}
  V r=reconstruct(f.r,f,ur);if(r[0]<=0||r[3]<=0)throw std::runtime_error("nonpositive reconstructed right state");if(coupled_flux)return mesh.gas.coupled_roe(l,r,f.nx,f.ny,uref2[f.l],uref2[f.r],mesh.gas.g*w[f.l][3]/w[f.l][0],mesh.gas.g*w[f.r][3]/w[f.r][0]);return mesh.gas.roe(l,r,f.nx,f.ny,native_numerics?.5*(uref2[f.l]+uref2[f.r]):-1);
 }
 VV residual(bool assemble,double cfl){
  gradients();if(native_numerics)update_reference_velocity();VV res(mesh.n);
#if CFD_CYLINDER_SKIP_UNUSED_RADIUS_WORK
  std::vector<double>sum;if(assemble)sum.assign(mesh.n,0);
#else
  std::vector<double>sum(mesh.n,0);
#endif
  if(assemble)matrix.clear();
  for(auto&f:mesh.faces){V right=f.r>=0?q[f.r]:q[f.l],base=faceflux(f,q[f.l],right);V value=f.area*base;res[f.l]+=value;if(f.r>=0)res[f.r]-=value;
#if CFD_CYLINDER_SKIP_UNUSED_RADIUS_WORK
   if(assemble){
#endif
   double lam=f.area*mesh.gas.radius(w[f.l],f.nx,f.ny,native_numerics?uref2[f.l]:-1);sum[f.l]+=lam;if(f.r>=0)sum[f.r]+=f.area*mesh.gas.radius(w[f.r],f.nx,f.ny,native_numerics?uref2[f.r]:-1);
#if CFD_CYLINDER_SKIP_UNUSED_RADIUS_WORK
   }
#endif
   if(assemble&&analytic_jacobian&&coupled_flux&&f.r>=0){
    V left=reconstruct(f.l,f,q[f.l]),rightp=reconstruct(f.r,f,right);std::array<Dual8,4>ld,rd;
    auto seed=[&](std::array<Dual8,4>&v,const V&rec,int cell,int offset){for(int k=0;k<4;++k)v[k].value=rec[k];const V&c=w[cell];v[0].d[offset]=1;v[1].d[offset]=-c[1]/c[0];v[1].d[offset+1]=1/c[0];v[2].d[offset]=-c[2]/c[0];v[2].d[offset+2]=1/c[0];v[3].d[offset]=.5*(mesh.gas.g-1)*(c[1]*c[1]+c[2]*c[2]);v[3].d[offset+1]=-(mesh.gas.g-1)*c[1];v[3].d[offset+2]=-(mesh.gas.g-1)*c[2];v[3].d[offset+3]=mesh.gas.g-1;};
    seed(ld,left,f.l,0);seed(rd,rightp,f.r,4);
    auto flux=coupled_roe_split(ld,rd,f.nx,f.ny,mesh.gas.g,uref2[f.l],uref2[f.r],mesh.gas.g*w[f.l][3]/w[f.l][0],mesh.gas.g*w[f.r][3]/w[f.r][0]);
    for(int k=0;k<4;++k)for(int j=0;j<4;++j){double dl=f.area*flux[k].d[j],dr=f.area*flux[k].d[j+4];matrix.a[f.ll][4*k+j]+=dl;matrix.a[f.rl][4*k+j]-=dl;matrix.a[f.lr][4*k+j]+=dr;matrix.a[f.rr][4*k+j]-=dr;}
   }else if(assemble)for(int j=0;j<4;++j){V pert=q[f.l];double eps=1e-7*std::max(1.,std::abs(pert[j]));pert[j]+=eps;V df=(f.area/eps)*(faceflux(f,pert,right)-base);
     for(int k=0;k<4;++k){matrix.a[f.ll][4*k+j]+=df[k];if(f.r>=0)matrix.a[f.rl][4*k+j]-=df[k];}
     if(f.r>=0){pert=right;eps=1e-7*std::max(1.,std::abs(pert[j]));pert[j]+=eps;df=(f.area/eps)*(faceflux(f,q[f.l],pert)-base);for(int k=0;k<4;++k){matrix.a[f.lr][4*k+j]+=df[k];matrix.a[f.rr][4*k+j]-=df[k];}}
   }
  }
  if(assemble)for(int i=0;i<mesh.n;++i){auto state=mesh.gas.state(w[i],1,0,native_numerics?uref2[i]:-1);double scale=sum[i]/(2*cfl);for(int j=0;j<4;++j){V unit{};unit[j]=1;V column=state.gamma_apply(state.ci(unit));for(int k=0;k<4;++k)matrix.a[matrix.diag[i]][4*k+j]+=scale*column[k];}}
  return res;
 }
 std::array<double,2>forces(std::ofstream*wall=nullptr){
  gradients();if(native_numerics)update_reference_velocity();std::array<double,2>force{};if(wall)*wall<<"x,y,numerical_wall_pressure,reconstructed_pressure,cell_pressure,cell\n"<<std::setprecision(16);
  for(auto&f:mesh.faces)if(f.type==3){V flux=faceflux(f,q[f.l],q[f.l]);force[0]+=flux[1]*f.area;force[1]+=flux[2]*f.area;
   if(wall)*wall<<f.x<<','<<f.y<<','<<(flux[1]*f.nx+flux[2]*f.ny)*mesh.rho0*mesh.a0*mesh.a0<<','<<reconstruct(f.l,f,q[f.l])[3]*mesh.rho0*mesh.a0*mesh.a0<<','<<w[f.l][3]*mesh.rho0*mesh.a0*mesh.a0<<','<<f.l+1<<'\n';
  }for(double&f:force)f*=mesh.rho0*mesh.a0*mesh.a0;return force;
 }
};

void tests(){
 Gas gas;V l{1,.1,.02,1/gas.g},r{.9,-.03,.07,.95/gas.g};
 for(bool pre:{false,true}){gas.precondition=pre;V physical=gas.physical(l,.6,.8),same=gas.roe(l,l,.6,.8);for(int k=0;k<4;++k)if(std::abs(physical[k]-same[k])>1e-12)throw std::runtime_error("constant-state consistency failed");
  V a=gas.roe(l,r,.6,.8),b=gas.roe(r,l,-.6,-.8);for(int k=0;k<4;++k)if(std::abs(a[k]+b[k])>1e-11)throw std::runtime_error("interior face conservation failed");
  V c=gas.coupled_roe(l,l,.6,.8,.02,.02,gas.g*l[3]/l[0],gas.g*l[3]/l[0]);for(int k=0;k<4;++k)if(std::abs(c[k]-physical[k])>1e-12)throw std::runtime_error("coupled flux constant-state consistency failed");
  V ca=gas.coupled_roe(l,r,.6,.8,.02,.03,gas.g*l[3]/l[0],gas.g*r[3]/r[0]),cb=gas.coupled_roe(r,l,-.6,-.8,.03,.02,gas.g*r[3]/r[0],gas.g*l[3]/l[0]);for(int k=0;k<4;++k)if(std::abs(ca[k]+cb[k])>1e-11)throw std::runtime_error("coupled flux conservation failed");
  auto state=gas.state(l,.6,.8);V z{.12,.07,-.11,.03},back=state.gamma_inverse(state.gamma_apply(z));for(int k=0;k<4;++k)if(std::abs(z[k]-back[k])>1e-11)throw std::runtime_error("preconditioner inverse failed");
 }
 M m{4,1,2,0,1,5,0,1,2,0,6,1,0,1,1,3};M p=mm(m,inv(m));for(int i=0;i<4;++i)for(int j=0;j<4;++j)if(std::abs(p[4*i+j]-(i==j?1:0))>1e-12)throw std::runtime_error("block inverse failed");
 std::array<Dual8,4>ld,rd;for(int j=0;j<4;++j){ld[j].value=l[j];ld[j].d[j]=1;rd[j].value=r[j];rd[j].d[j+4]=1;}
 auto fd=coupled_roe_split(ld,rd,.6,.8,gas.g,.02,.03,1.,1.1);
 for(int j=0;j<8;++j){V lp=l,lm=l,rp=r,rm=r;double eps=1e-6;if(j<4){lp[j]+=eps;lm[j]-=eps;}else{rp[j-4]+=eps;rm[j-4]-=eps;}V fp=gas.coupled_roe(lp,rp,.6,.8,.02,.03,1.,1.1),fm=gas.coupled_roe(lm,rm,.6,.8,.02,.03,1.,1.1);for(int k=0;k<4;++k)if(std::abs((fp[k]-fm[k])/(2*eps)-fd[k].d[j])>1e-7)throw std::runtime_error("analytic face Jacobian failed central-difference verification");}
 Anderson aa(3);VV x0(1,V{}),s0(1,V{1,2,3,4});aa.propose(x0,s0);VV x1=s0,s1=s0;for(auto&v:s1)for(double&z:v)z*=.9;VV update=aa.propose(x1,s1);for(int j=0;j<4;++j)if(std::abs(x1[0][j]+update[0][j]-10*s0[0][j])>1e-10)throw std::runtime_error("Anderson linear fixed-point check failed");
 std::cout<<"constant-state, face conservation, preconditioner inverse and block inverse: passed\n";
}

int cylinder_run_once(int argc,char**argv){try{
 std::string input="data/cylinder/cylinder.bin",out="runs/manual",limiter="barth",linear_mode="gmres",jacobian="finite-difference";int maxiter=500,order=2,linearmax=60,jacperiod=1,precperiod=1,fill_level=0,warmup=0,aa_depth=0,aa_start=5;double cfl=50,tol=1e-7,ltol=.03;bool pre=true,wall_extra=false,native=false,coupled=false;
 for(int i=1;i<argc;++i){std::string s=argv[i];auto val=[&](){if(++i>=argc)throw std::runtime_error("missing argument");return std::string(argv[i]);};
  if(s=="--self-test"){tests();return 0;}else if(s=="--input")input=val();else if(s=="--output")out=val();else if(s=="--iterations")maxiter=std::stoi(val());else if(s=="--order")order=std::stoi(val());else if(s=="--limiter")limiter=val();else if(s=="--tol")tol=std::stod(val());else if(s=="--linear-tol")ltol=std::stod(val());else if(s=="--linear-max")linearmax=std::stoi(val());else if(s=="--jacobian-period")jacperiod=std::stoi(val());else if(s=="--preconditioner-period")precperiod=std::stoi(val());else if(s=="--ilu-fill")fill_level=std::stoi(val());else if(s=="--jacobian-warmup")warmup=std::stoi(val());else if(s=="--anderson-depth")aa_depth=std::stoi(val());else if(s=="--anderson-start")aa_start=std::stoi(val());else if(s=="--jacobian")jacobian=val();else if(s=="--linear-mode")linear_mode=val();else if(s=="--no-precondition")pre=false;else if(s=="--native-numerics")native=true;else if(s=="--coupled-roe"){coupled=true;native=true;}else if(s=="--wall-extrapolate")wall_extra=true;else if(s=="--cfl")cfl=std::stod(val());else throw std::runtime_error("unknown option: "+s);
 }
 if(jacperiod<1||linearmax<1||maxiter<1||cfl<=0||tol<=0||ltol<=0)throw std::runtime_error("iteration counts and tolerances must be positive");
 if(precperiod<1||fill_level<0||fill_level>3||warmup<0)throw std::runtime_error("invalid preconditioner, fill or warmup setting");
 if(aa_depth<0||aa_depth>20||aa_start<0)throw std::runtime_error("invalid Anderson history settings");
 if(linear_mode!="gmres"&&linear_mode!="ilu")throw std::runtime_error("linear mode must be gmres or ilu");
 if(jacobian!="finite-difference"&&jacobian!="analytic")throw std::runtime_error("jacobian must be finite-difference or analytic");
 if(limiter!="barth"&&limiter!="venkat"&&limiter!="differentiable"&&limiter!="standard"&&limiter!="none")throw std::runtime_error("unknown limiter");
 auto start=std::chrono::steady_clock::now();Solver solver(input,fill_level);solver.order=order;solver.limiter=limiter;solver.mesh.gas.precondition=pre;solver.wall_extrapolate=wall_extra;solver.native_numerics=native;solver.coupled_flux=coupled;solver.analytic_jacobian=jacobian=="analytic";
 std::filesystem::create_directories(out);std::ofstream history(out+"/history.csv");history<<"iteration,res_mass,res_xmom,res_ymom,res_energy,scaled_max,linear_iterations,alpha,wall_seconds,jacobian_refreshed\n"<<std::setprecision(16);
 V scale{};double rel=1;int done=0,alllinear=0,jacobians=0,aa_accepted=0,aa_rejected=0,cached_evaluations=0;Anderson anderson(aa_depth);VV cached_residual;bool cache_valid=false;double assembly=0,linear=0,acceleration=0;auto iterationstart=std::chrono::steady_clock::now();
 for(int it=0;it<maxiter;++it){bool refresh=it<warmup||it%jacperiod==0;auto t0=std::chrono::steady_clock::now();VV residual;if(!refresh&&cache_valid){residual=std::move(cached_residual);++cached_evaluations;}else residual=solver.residual(refresh,cfl);cache_valid=false;auto t1=std::chrono::steady_clock::now();assembly+=std::chrono::duration<double>(t1-t0).count();V norms{};for(auto&r:residual)for(int k=0;k<4;++k)norms[k]+=std::abs(r[k]);
  if(it<5)for(int k=0;k<4;++k)scale[k]=std::max(scale[k],norms[k]);rel=0;for(int k=0;k<4;++k)rel=std::max(rel,norms[k]/std::max(scale[k],1e-30));
  if(rel<tol){done=it;break;}
  if(refresh){++jacobians;if(it<warmup||it%precperiod==0)solver.matrix.factor();}for(auto&r:residual)for(double&v:r)v=-v;VV delta;int nl=1;if(linear_mode=="ilu")delta=solver.matrix.solve(residual);else nl=solver.matrix.gmres(residual,delta,25,linearmax,ltol);alllinear+=nl;auto t2=std::chrono::steady_clock::now();linear+=std::chrono::duration<double>(t2-t1).count();
  if(aa_depth>0&&it>=aa_start){VV candidate=anderson.propose(solver.q,delta);if(it>aa_start){VV saved=solver.q,rr;bool okay=true;for(int i=0;i<solver.mesh.n;++i){solver.q[i]+=candidate[i];V state=solver.mesh.gas.primitive(solver.q[i]);if(!(state[0]>1e-6&&state[3]>1e-6&&std::isfinite(state[3])))okay=false;}double trial=std::numeric_limits<double>::infinity();if(okay)try{rr=solver.residual(false,cfl);V nn{};for(const auto&v:rr)for(int k=0;k<4;++k)nn[k]+=std::abs(v[k]);trial=0;for(int k=0;k<4;++k)trial=std::max(trial,nn[k]/std::max(scale[k],1e-30));}catch(const std::runtime_error&){okay=false;}solver.q=std::move(saved);if(okay&&trial<rel*.99){delta=std::move(candidate);cached_residual=std::move(rr);cache_valid=true;++aa_accepted;}else{anderson.restart();++aa_rejected;}}acceleration+=std::chrono::duration<double>(std::chrono::steady_clock::now()-t2).count();}
  double alpha=1;bool admissible=false;
  for(int trial=0;trial<20;++trial){admissible=true;for(int i=0;i<solver.mesh.n;++i){V w=solver.mesh.gas.primitive(solver.q[i]+alpha*delta[i]);if(!(w[0]>1e-6&&w[3]>1e-6&&std::isfinite(w[3]))){admissible=false;break;}}if(admissible)break;alpha*=.5;}
  if(!admissible)throw std::runtime_error("positivity line search exhausted");
  if(alpha!=1)cache_valid=false;
  for(int i=0;i<solver.mesh.n;++i)solver.q[i]+=alpha*delta[i];done=it+1;double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-iterationstart).count();
  history<<done;for(double v:norms)history<<','<<v;history<<','<<rel<<','<<nl<<','<<alpha<<','<<elapsed<<','<<refresh<<'\n';
  if(it<5||done%20==0)std::cout<<done<<" residual="<<rel<<" linear="<<nl<<" alpha="<<alpha<<" time="<<elapsed<<'\n';
  if(!std::isfinite(rel)||rel>1e9)throw std::runtime_error("nonlinear divergence");
 }
 double itertime=std::chrono::duration<double>(std::chrono::steady_clock::now()-iterationstart).count();std::ofstream wall(out+"/wall_pressure.csv");auto force=solver.forces(&wall);
 std::ofstream field(out+"/field.csv");field<<"cell,x,y,rho,u,v,p\n"<<std::setprecision(16);for(int i=0;i<solver.mesh.n;++i){V w=solver.mesh.gas.primitive(solver.q[i]);field<<i+1<<','<<solver.mesh.center[i][0]<<','<<solver.mesh.center[i][1]<<','<<w[0]*solver.mesh.rho0<<','<<w[1]*solver.mesh.a0<<','<<w[2]*solver.mesh.a0<<','<<w[3]*solver.mesh.rho0*solver.mesh.a0*solver.mesh.a0<<'\n';}
 double dynamic=.5*solver.mesh.rho0*std::pow(solver.inf[1]*solver.mesh.a0,2);
 std::ofstream summary(out+"/summary.json");summary<<std::setprecision(16)<<"{\n\"cached_residual_evaluations\":"<<cached_evaluations<<",\n\"jacobian\":\""<<jacobian<<"\",\n\"ilu_fill\":"<<fill_level<<",\n\"preconditioner_period\":"<<precperiod<<",\n\"jacobian_warmup\":"<<warmup<<",\n\"anderson_depth\":"<<aa_depth<<",\n\"anderson_start\":"<<aa_start<<",\n\"anderson_accepted\":"<<aa_accepted<<",\n\"anderson_rejected\":"<<aa_rejected<<",\n\"acceleration_seconds\":"<<acceleration<<",\n\"jacobian_period\":"<<jacperiod<<",\n\"jacobian_assemblies\":"<<jacobians<<",\n\"linear_mode\":\""<<linear_mode<<"\",\n\"coupled_roe\":"<<(coupled?"true":"false")<<",\n\"local_preconditioning_and_wall_dissipation\":"<<(native?"true":"false")<<",\n\"residual_tolerance\":"<<tol<<",\n\"linear_tolerance\":"<<ltol<<",\n\"iterations\":"<<done<<",\n\"cfl\":"<<cfl<<",\n\"order\":"<<order<<",\n\"limiter\":\""<<limiter<<"\",\n\"precondition\":"<<(pre?"true":"false")<<",\n\"residual_scaled_max\":"<<rel<<",\n\"converged\":"<<(rel<tol?"true":"false")<<",\n\"iteration_seconds\":"<<itertime<<",\n\"assembly_seconds\":"<<assembly<<",\n\"linear_seconds\":"<<linear<<",\n\"linear_iterations_total\":"<<alllinear<<",\n\"fx\":"<<force[0]<<",\n\"fy\":"<<force[1]<<",\n\"cd\":"<<force[0]/dynamic<<",\n\"cl\":"<<force[1]/dynamic<<",\n\"total_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\n}\n";
#if CFD_CYLINDER_INCREMENTAL_ANDERSON_QR
 std::cout<<"Anderson QR incremental="<<anderson.incremental_columns<<" downdates="<<anderson.downdates<<" rebuilds="<<anderson.full_rebuilds<<'\n';
#endif
 std::cout<<"RESULT iterations="<<done<<" residual="<<rel<<" seconds="<<itertime<<" Cd="<<force[0]/dynamic<<" Cl="<<force[1]/dynamic<<'\n';return rel<tol?0:2;
 }catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}}

int main(int argc,char**argv){try{
 int repeats=1;std::vector<char*>args;
 for(int i=0;i<argc;++i){if(std::string(argv[i])=="--profile-repeats"){if(++i>=argc)throw std::runtime_error("missing profile repetition count");repeats=std::stoi(argv[i]);}else args.push_back(argv[i]);}
 if(repeats<1||repeats>1000)throw std::runtime_error("profile repetitions must be 1..1000");
 for(int i=0;i<repeats;++i){int result=cylinder_run_once(static_cast<int>(args.size()),args.data());if(result)return result;}
 return 0;
 }catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}}
