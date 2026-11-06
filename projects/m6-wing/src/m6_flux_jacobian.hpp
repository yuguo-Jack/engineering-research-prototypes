#pragma once
#include "m6_physics.hpp"

namespace m6 {
// Differentiate both five-variable face states in a single forward pass.
// Primitive reconstruction increments, cell acoustics and Uref are frozen.
// This is a local flux Jacobian, not the complete second-order stencil.
struct FaceDual {
 double value=0;std::array<double,10>d{};
 FaceDual()=default;FaceDual(double v):value(v){}
};
inline FaceDual operator+(const FaceDual&a,const FaceDual&b){FaceDual c(a.value+b.value);for(int j=0;j<10;++j)c.d[j]=a.d[j]+b.d[j];return c;}
inline FaceDual operator-(const FaceDual&a,const FaceDual&b){FaceDual c(a.value-b.value);for(int j=0;j<10;++j)c.d[j]=a.d[j]-b.d[j];return c;}
inline FaceDual operator*(const FaceDual&a,const FaceDual&b){FaceDual c(a.value*b.value);for(int j=0;j<10;++j)c.d[j]=a.d[j]*b.value+a.value*b.d[j];return c;}
inline FaceDual operator/(const FaceDual&a,const FaceDual&b){FaceDual c(a.value/b.value);for(int j=0;j<10;++j)c.d[j]=(a.d[j]-c.value*b.d[j])/b.value;return c;}
inline FaceDual face_sqrt(const FaceDual&a){FaceDual c(std::sqrt(a.value));for(int j=0;j<10;++j)c.d[j]=a.d[j]/(2*c.value);return c;}
inline FaceDual face_abs(const FaceDual&a){FaceDual c(std::abs(a.value));double s=a.value>0?1.:a.value<0?-1.:0.;for(int j=0;j<10;++j)c.d[j]=s*a.d[j];return c;}
struct FluxJacobians {B5 left{},right{};};
inline FluxJacobians coupled_flux_jacobian(const Gas&gas,const W7&cell_l,const W7&cell_r,const W7&face_l,const W7&face_r,const D3&n,double a2_l,double a2_r,double ref_l,double ref_r,const std::array<double,4>*face_thermo=nullptr){
 using D=FaceDual;std::array<D,5>l{},r{};
 auto seed=[&](const W7&cell,const W7&face,int offset,std::array<D,5>&v){auto dp=gas.pressure_jacobian(cell);auto dw=gas.primitive_jacobian(cell);for(int k=0;k<5;++k){v[k].value=face[k];for(int j=0;j<5;++j)v[k].d[offset+j]=k==0?dp[j]:dw[5*k+j];}};
 seed(cell_l,face_l,0,l);seed(cell_r,face_r,5,r);
 double a2=.5*(a2_l+a2_r),a=std::sqrt(a2),beta=std::clamp(2*std::max(ref_l,ref_r)/a2,1e-20,1.);
 D ml=(l[1]*n[0]+l[2]*n[1]+l[3]*n[2])/a,mr=(r[1]*n[0]+r[2]*n[1]+r[3]*n[2])/a,average=.5*(face_abs(ml)+face_abs(mr));
 D factor=face_sqrt(4*beta+(1-beta)*(1-beta)*average*average)/(1+beta),c=a*factor;
 D left=((1+beta)*ml+(1-beta)*mr)/(2*factor),right=((1-beta)*ml+(1+beta)*mr)/(2*factor);
 auto mp=[](const D&m)->D{return std::abs(m.value)>=1?.5*(m+face_abs(m)):.25*(m+1)*(m+1)*(1+.5*(m-1)*(m-1));};
 auto mm=[](const D&m)->D{return std::abs(m.value)>=1?.5*(m-face_abs(m)):(0.-.25)*(m-1)*(m-1)*(1+.5*(m+1)*(m+1));};
 auto pp=[](const D&m)->D{return std::abs(m.value)>=1?D(m.value>=0?1.:0.):.25*(m+1)*(m+1)*(2-m)+.1875*m*(m*m-1)*(m*m-1);};
 auto pm=[](const D&m)->D{return std::abs(m.value)>=1?D(m.value<=0?1.:0.):.25*(m-1)*(m-1)*(2+m)-.1875*m*(m*m-1)*(m*m-1);};
 D mplus=mp(left),mminus=mm(right),mface=mplus+mminus,rl=l[0]/(gas.rg*l[4]),rr=r[0]/(gas.rg*r[4]);
 if(face_thermo){rl=D((*face_thermo)[0]);rr=D((*face_thermo)[1]);rl.d[0]=1;rr.d[5]=1;}
 D mass=c*mface*(mface.value>=0?rl:rr)+(1/beta-1)*c*(l[0]-r[0])*((0.-.5)*(left+face_abs(left))+.5*(right-face_abs(right))+mplus-mminus)/(l[0]/rl+r[0]/rr);
 D pplus=pp(left),pminus=pm(right),pressure=l[0]*pplus+r[0]*pminus+.5*(rl+rr)*c*c*(left-right)*pplus*pminus;
 const auto&up=mass.value>=0?l:r;const auto&face_up=mass.value>=0?face_l:face_r;
 D enthalpy(gas.enthalpy(face_up));for(int j=0;j<10;++j){enthalpy.d[j]=gas.cp(face_up[4])*up[4].d[j];for(int k=1;k<=3;++k)enthalpy.d[j]+=face_up[k]*up[k].d[j];}
 if(face_thermo){bool left_up=mass.value>=0;const auto&cell=left_up?cell_l:cell_r;auto dp=gas.pressure_jacobian(cell);double rho=gas.rho(cell),H=gas.enthalpy(cell);enthalpy=D((*face_thermo)[left_up?2:3]);for(int j=0;j<5;++j)enthalpy.d[(left_up?0:5)+j]=((j==4?1.:0.)+dp[j]-H*(j==0?1.:0.))/rho;}
 std::array<D,5>f={mass,mass*up[1]+pressure*n[0],mass*up[2]+pressure*n[1],mass*up[3]+pressure*n[2],mass*enthalpy};
 FluxJacobians result;for(int row=0;row<5;++row)for(int col=0;col<5;++col){result.left[5*row+col]=f[row].d[col];result.right[5*row+col]=f[row].d[5+col];}return result;
}
}
