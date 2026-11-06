#pragma once
#include <array>
#include <cmath>

// Forward differentiation in the eight conservative inputs of an interior face.
// The reconstruction increments and preconditioning coefficients are frozen
// during a linearization, exactly as in the finite-difference baseline.
struct Dual8 {
 double value=0;std::array<double,8> d{};
 Dual8()=default;Dual8(double v):value(v){}
};
inline Dual8 operator+(const Dual8&a,const Dual8&b){Dual8 c(a.value+b.value);for(int i=0;i<8;++i)c.d[i]=a.d[i]+b.d[i];return c;}
inline Dual8 operator-(const Dual8&a,const Dual8&b){Dual8 c(a.value-b.value);for(int i=0;i<8;++i)c.d[i]=a.d[i]-b.d[i];return c;}
inline Dual8 operator*(const Dual8&a,const Dual8&b){Dual8 c(a.value*b.value);for(int i=0;i<8;++i)c.d[i]=a.d[i]*b.value+a.value*b.d[i];return c;}
inline Dual8 operator/(const Dual8&a,const Dual8&b){Dual8 c(a.value/b.value);for(int i=0;i<8;++i)c.d[i]=(a.d[i]-c.value*b.d[i])/b.value;return c;}
inline Dual8 sqrt(const Dual8&a){Dual8 c(std::sqrt(a.value));for(int i=0;i<8;++i)c.d[i]=a.d[i]/(2*c.value);return c;}
inline Dual8 abs(const Dual8&a){Dual8 c(std::abs(a.value));double s=a.value>0?1:a.value<0?-1:0;for(int i=0;i<8;++i)c.d[i]=s*a.d[i];return c;}
inline double scalar_value(double v){return v;}
inline double scalar_value(const Dual8&v){return v.value;}

template<class T>
std::array<T,4> coupled_roe_split(const std::array<T,4>&l,const std::array<T,4>&r,
 double nx,double ny,double g,double uref_l,double uref_r,double cell_a2_l,double cell_a2_r){
 using std::sqrt;using std::abs;
 T un=.5*((l[1]+r[1])*nx+(l[2]+r[2])*ny);
 const auto&up=scalar_value(un)>0?l:r;
 auto enthalpy=[&](const std::array<T,4>&v){return g/(g-1)*v[3]/v[0]+.5*(v[1]*v[1]+v[2]*v[2]);};
 T mass=up[0]*un,pbar=.5*(l[3]+r[3]);
 std::array<T,4>f{mass,mass*up[1]+pbar*nx,mass*up[2]+pbar*ny,mass*enthalpy(up)};
 double alpha=.5*(1-.5*(uref_l/cell_a2_l+uref_r/cell_a2_r)),u2=.5*(uref_l+uref_r);
 T c=sqrt(alpha*alpha*un*un+u2),lp=abs((1-alpha)*un+c),lm=abs((1-alpha)*un-c);
 T dl=(lp-lm)/(2*c),spread=.5*(lp+lm)-abs(un),du=(l[1]-r[1])*nx+(l[2]-r[2])*ny,dp=l[3]-r[3];
 T dm=.5*(dl*up[0]*du+dp/u2*(2*alpha*abs(un)-alpha*dl*un+spread));
 T pm=.5*(up[0]*(spread+alpha*un*dl)*du+dp*dl);
 const auto&diss_up=scalar_value(dm)>0?l:r;
 f[0]=f[0]+dm;f[1]=f[1]+dm*diss_up[1]+pm*nx;f[2]=f[2]+dm*diss_up[2]+pm*ny;f[3]=f[3]+dm*enthalpy(diss_up);
 return f;
}
