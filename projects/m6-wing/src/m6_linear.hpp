#pragma once
#include "m6_mesh.hpp"
#include <limits>
#include <functional>

namespace m6 {
template<int N>using Vec=std::array<double,N>;
template<int N>using Block=std::array<double,N*N>;
template<int N>using Field=std::vector<Vec<N>>;
struct KrylovStart{double relative_residual=1,coefficient=0,trial_relative_residual=1;bool recycled=false;};
// This tiny kernel is called once per sparse block, including triangular solves.
// Keep its arithmetic order while allowing the caller to eliminate the returned
// temporary and expose fixed-size loops to the optimizer.
#if defined(_MSC_VER)
#define M6_BLOCK_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define M6_BLOCK_INLINE inline __attribute__((always_inline))
#else
#define M6_BLOCK_INLINE inline
#endif
template<int N>M6_BLOCK_INLINE Vec<N>bmultiply(const Block<N>&a,const Vec<N>&x){Vec<N>y{};for(int i=0;i<N;++i)for(int j=0;j<N;++j)y[i]+=a[i*N+j]*x[j];return y;}
#undef M6_BLOCK_INLINE
template<int N>Block<N>bproduct(const Block<N>&a,const Block<N>&b){Block<N>c{};for(int i=0;i<N;++i)for(int k=0;k<N;++k)for(int j=0;j<N;++j)c[i*N+j]+=a[i*N+k]*b[k*N+j];return c;}
template<int N>Block<N>binverse(Block<N>a){
 Block<N>b{};for(int i=0;i<N;++i)b[i*(N+1)]=1;
 double scale=0;for(double v:a)scale=std::max(scale,std::abs(v));
 for(int k=0;k<N;++k){int p=k;for(int i=k+1;i<N;++i)if(std::abs(a[i*N+k])>std::abs(a[p*N+k]))p=i;
  if(std::abs(a[p*N+k])<1e-17*scale||!std::isfinite(a[p*N+k]))throw std::runtime_error("singular/nonfinite implicit block");
  for(int j=0;j<N;++j){std::swap(a[k*N+j],a[p*N+j]);std::swap(b[k*N+j],b[p*N+j]);}double d=a[k*N+k];for(int j=0;j<N;++j){a[k*N+j]/=d;b[k*N+j]/=d;}
  for(int i=0;i<N;++i)if(i!=k){double t=a[i*N+k];for(int j=0;j<N;++j){a[i*N+j]-=t*a[k*N+j];b[i*N+j]-=t*b[k*N+j];}}
 }return b;
}
template<int N>class Sparse {
 const Mesh&mesh;const MPIWorld&world;
 public:
 int n;std::vector<int>row,col,diag;std::vector<Block<N>>a,lu,dinv;
 Sparse(const Mesh&m,const MPIWorld&w):mesh(m),world(w),n(m.n),row(n+1),diag(n){
  std::vector<std::vector<int>>adj(n);for(int i=0;i<n;++i)adj[i].push_back(i);
  for(auto&f:m.faces)if(f.r>=0){adj[f.l].push_back(f.r);if(f.r<n)adj[f.r].push_back(f.l);}
  for(int i=0;i<n;++i){auto&v=adj[i];std::sort(v.begin(),v.end());v.erase(std::unique(v.begin(),v.end()),v.end());row[i]=static_cast<int>(col.size());col.insert(col.end(),v.begin(),v.end());}row[n]=static_cast<int>(col.size());
  a.resize(col.size());lu.resize(col.size());dinv.resize(n);for(int i=0;i<n;++i)diag[i]=index(i,i);
 }
 int index(int i,int j)const{auto it=std::lower_bound(col.begin()+row[i],col.begin()+row[i+1],j);return it!=col.begin()+row[i+1]&&*it==j?static_cast<int>(it-col.begin()):-1;}
 void clear(){std::fill(a.begin(),a.end(),Block<N>{});}
 void factor(bool full=true){
  lu=a;for(int i=0;i<n;++i){if(full)for(int p=row[i];p<diag[i];++p){int j=col[p];auto l=bproduct<N>(lu[p],dinv[j]);lu[p]=l;
    for(int q=diag[j]+1;q<row[j+1]&&col[q]<n;++q){int where=index(i,col[q]);if(where>=0){auto v=bproduct<N>(l,lu[q]);for(int k=0;k<N*N;++k)lu[where][k]-=v[k];}}
   }dinv[i]=binverse<N>(lu[diag[i]]);
  }
 }
 void multiply(Field<N>&x,Field<N>&y)const{
  mesh.exchange(x);y.assign(n,Vec<N>{});for(int i=0;i<n;++i)for(int p=row[i];p<row[i+1];++p){auto v=bmultiply<N>(a[p],x[col[p]]);for(int k=0;k<N;++k)y[i][k]+=v[k];}
 }
 void solve(const Field<N>&b,Field<N>&x)const{
  x=b;for(int i=0;i<n;++i)for(int p=row[i];p<diag[i];++p){auto v=bmultiply<N>(lu[p],x[col[p]]);for(int k=0;k<N;++k)x[i][k]-=v[k];}
  for(int i=n-1;i>=0;--i){for(int p=diag[i]+1;p<row[i+1]&&col[p]<n;++p){auto v=bmultiply<N>(lu[p],x[col[p]]);for(int k=0;k<N;++k)x[i][k]-=v[k];}x[i]=bmultiply<N>(dinv[i],x[i]);}
 }
 double dot(const Field<N>&x,const Field<N>&y)const{double t=0;for(int i=0;i<n;++i)for(int k=0;k<N;++k)t+=x[i][k]*y[i][k];return world.sum(t);}
 static void axpy(Field<N>&x,double s,const Field<N>&y,int count){for(int i=0;i<count;++i)for(int k=0;k<N;++k)x[i][k]+=s*y[i][k];}
 int gmres(const Field<N>&b,Field<N>&x,int maxsteps=8,double tolerance=.05,const std::function<void(const Field<N>&,Field<N>&)>&extra={},bool batched_orthogonalization=false,double*relative_residual=nullptr,const Field<N>*recycled_direction=nullptr,KrylovStart*start_info=nullptr,const std::function<void(const Field<N>&,Field<N>&)>&preconditioner={})const{
  x.assign(n,Vec<N>{});if(start_info)*start_info=KrylovStart{};double bn=std::sqrt(dot(b,b));if(relative_residual)*relative_residual=0;if(bn<1e-35){if(start_info)start_info->relative_residual=0;return 0;}
  Field<N>initial_residual=b;double beta=bn;
  // A one-dimensional minimum-residual prediction in the CURRENT operator.
  // The convergence target stays relative to the original right-hand side.
  if(recycled_direction&&!recycled_direction->empty()){
   if(recycled_direction->size()<static_cast<size_t>(n))throw std::runtime_error("short recycled Krylov direction");
   Field<N>direction(recycled_direction->begin(),recycled_direction->begin()+n),applied;multiply(direction,applied);if(extra)extra(direction,applied);
   std::array<double,2>products{};for(int i=0;i<n;++i)for(int j=0;j<N;++j){products[0]+=b[i][j]*applied[i][j];products[1]+=applied[i][j]*applied[i][j];}products=world.sum(products);
   if(products[1]>1e-60&&std::isfinite(products[0])&&std::isfinite(products[1])){double coefficient=products[0]/products[1];
    for(int i=0;i<n;++i)for(int j=0;j<N;++j)initial_residual[i][j]=b[i][j]-coefficient*applied[i][j];double trial=std::sqrt(dot(initial_residual,initial_residual));
    if(start_info)start_info->trial_relative_residual=trial/bn;
    if(std::isfinite(trial)&&trial<=.95*bn){beta=trial;for(int i=0;i<n;++i)for(int j=0;j<N;++j)x[i][j]=coefficient*direction[i][j];if(start_info){start_info->relative_residual=beta/bn;start_info->coefficient=coefficient;start_info->recycled=true;}}
    else initial_residual=b;
   }
  }
  if(beta<=tolerance*bn){if(relative_residual)*relative_residual=beta/bn;return 0;}
  std::vector<Field<N>>v(maxsteps+1),z(maxsteps);v[0]=std::move(initial_residual);for(auto&item:v[0])for(double&t:item)t/=beta;
  std::vector<std::vector<double>>h(maxsteps+1,std::vector<double>(maxsteps));std::vector<double>cs(maxsteps),sn(maxsteps),s(maxsteps+1);s[0]=beta;int used=0;
  for(int j=0;j<maxsteps;++j){if(preconditioner)preconditioner(v[j],z[j]);else solve(v[j],z[j]);multiply(z[j],v[j+1]);if(extra)extra(z[j],v[j+1]);
   // Both paths use two passes. CGS2 batches the dot products into one
   // collective per pass, reducing MPI synchronization as the basis grows.
   for(int pass=0;pass<2;++pass){
    if(batched_orthogonalization){std::vector<double>coefficients(j+1);
     for(int i=0;i<n;++i)for(int component=0;component<N;++component){double value=v[j+1][i][component];for(int k=0;k<=j;++k)coefficients[k]+=v[k][i][component]*value;}
     MPI_Allreduce(MPI_IN_PLACE,coefficients.data(),j+1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
     for(int k=0;k<=j;++k)h[k][j]+=coefficients[k];
     for(int i=0;i<n;++i)for(int component=0;component<N;++component){double correction=0;for(int k=0;k<=j;++k)correction+=coefficients[k]*v[k][i][component];v[j+1][i][component]-=correction;}
    }else for(int k=0;k<=j;++k){double t=dot(v[k],v[j+1]);h[k][j]+=t;axpy(v[j+1],-t,v[k],n);}
   }
   h[j+1][j]=std::sqrt(dot(v[j+1],v[j+1]));if(h[j+1][j]>1e-35)for(int i=0;i<n;++i)for(double&t:v[j+1][i])t/=h[j+1][j];
   for(int k=0;k<j;++k){double t=cs[k]*h[k][j]+sn[k]*h[k+1][j];h[k+1][j]=-sn[k]*h[k][j]+cs[k]*h[k+1][j];h[k][j]=t;}
   double d=std::hypot(h[j][j],h[j+1][j]);cs[j]=d>0?h[j][j]/d:1;sn[j]=d>0?h[j+1][j]/d:0;h[j][j]=d;h[j+1][j]=0;
   s[j+1]=-sn[j]*s[j];s[j]*=cs[j];used=j+1;if(std::abs(s[j+1])<tolerance*bn)break;
  }
  std::vector<double>y(used);for(int i=used-1;i>=0;--i){double t=s[i];for(int j=i+1;j<used;++j)t-=h[i][j]*y[j];y[i]=t/h[i][i];}
  for(int j=0;j<used;++j)axpy(x,y[j],z[j],n);if(relative_residual)*relative_residual=std::abs(s[used])/bn;return used;
 }
};
}
