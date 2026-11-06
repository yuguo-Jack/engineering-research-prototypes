#pragma once
// Independent calorically-imperfect mixture, AUSM and SST wall operators.
// Public model equations and supplied material coefficients; no vendor linkage.
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <memory>
#include "m6_material.hpp"

namespace m6 {
using D3=std::array<double,3>;
using W7=std::array<double,7>; // p, u, v, w, T, k, omega, nondimensional
using Q5=std::array<double,5>;
using B5=std::array<double,25>;
inline double dot3(const D3&a,const D3&b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
inline D3 sub3(const D3&a,const D3&b){return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};}
inline double length3(const D3&a){return std::sqrt(dot3(a,a));}
// Both projections reproduce a linear field. The edge projection remains
// bounded as the connector becomes tangent to the face (vanishing n dot dx).
inline D3 correct_face_gradient(D3 g,const D3&dx,const D3&normal,double delta,bool edge,double limit){
 double dn=dot3(dx,normal),distance2=dot3(dx,dx),error=delta-dot3(g,dx);
 if(edge&&dn/std::sqrt(distance2)<.2){for(int j=0;j<3;++j)g[j]+=error*dx[j]/distance2;return g;}
 double correction=error/dn;
 if(limit>=0&&dn/std::sqrt(distance2)<.2){double orthogonal=delta/dn,gn=dot3(g,normal),skew=gn-dot3(g,dx)/dn;correction=orthogonal+std::clamp(skew,-limit*std::abs(orthogonal),limit*std::abs(orthogonal))-gn;}
 for(int j=0;j<3;++j)g[j]+=correction*normal[j];return g;
}
inline double poly(const std::array<double,5>&c,double T){return ((((c[4]*T+c[3])*T+c[2])*T+c[1])*T+c[0]);}
inline double integral(const std::array<double,5>&c,double a,double b){
 double result=0,ap=a,bp=b;for(int j=0;j<5;++j){result+=c[j]*(bp-ap)/(j+1);ap*=a;bp*=b;}return result;
}
struct Gas {
 double rho0,a0,T0,R,pinf,mu,conductivity,cp0,rg;
 std::array<double,5>low,high;
 std::shared_ptr<const MixtureMaterial>material;
 explicit Gas(const std::array<double,48>&m):rho0(m[0]),a0(m[1]),T0(m[2]),R(m[3]),pinf(m[4]),mu(m[5]/(rho0*a0)),conductivity(m[6]*T0/(rho0*a0*a0*a0)),cp0(m[7]),rg(R*T0/(a0*a0)){
  for(int j=0;j<5;++j){low[j]=m[8+j];high[j]=m[13+j];}
 }
 void load_material(const std::string&path){material=std::make_shared<MixtureMaterial>(path);}
 double cp_phys(double T)const{if(material)return material->cp(T);return T<300?cp0:(T<1000?poly(low,T):poly(high,std::min(T,5000.)));}
 double h_phys(double T)const{
  if(material)return material->h(T);
  if(T<=300)return cp0*T;
  double h=cp0*300+integral(low,300,std::min(T,1000.));
  if(T>1000)h+=integral(high,1000,std::min(T,5000.));
  if(T>5000)h+=(T-5000)*poly(high,5000);return h;
 }
 double cp(double t)const{return cp_phys(t*T0)*T0/(a0*a0);}
 double cp_derivative(double t)const{double T=t*T0;if(material)return material->cp_derivative(T)*T0*T0/(a0*a0);if(T<300||T>5000)return 0;const auto&c=T<1000?low:high;return (((4*c[4]*T+3*c[3])*T+2*c[2])*T+c[1])*T0*T0/(a0*a0);}
 double h(double t)const{return h_phys(t*T0)/(a0*a0);}
 double temperature_from_h(double target,double guess)const{double lower_cp=material?material->cp(300):cp0;if(target<=lower_cp*300/(a0*a0))return target*a0*a0/(lower_cp*T0);double t=guess;for(int j=0;j<10;++j){double step=(h(t)-target)/cp(t);t-=step;if(!(t>0)||!std::isfinite(t))throw std::runtime_error("nonpositive reconstructed enthalpy temperature");if(std::abs(step)<1e-13*std::max(1.,t))break;}return t;}
 double rho(const W7&w)const{return w[0]/(rg*w[4]);}
 double sound(const W7&w)const{double c=cp(w[4]);return std::sqrt(c/(c-rg)*rg*w[4]);}
 double enthalpy(const W7&w)const{return h(w[4])+.5*(w[1]*w[1]+w[2]*w[2]+w[3]*w[3]);}
 Q5 conserved(const W7&w)const{double r=rho(w);return {r,r*w[1],r*w[2],r*w[3],r*(enthalpy(w)-rg*w[4])};}
 W7 primitive(const Q5&q,double tguess,double k,double omega)const{
  if(!(q[0]>0))throw std::runtime_error("nonpositive density in mixture inversion");
  double u=q[1]/q[0],v=q[2]/q[0],w=q[3]/q[0],e=q[4]/q[0]-.5*(u*u+v*v+w*w),t=tguess;
  for(int it=0;it<12;++it){double delta=(h(t)-rg*t-e)/(cp(t)-rg);t-=delta;if(!(t>0)||!std::isfinite(t))throw std::runtime_error("nonpositive temperature in mixture inversion");if(std::abs(delta)<1e-13*std::max(t,1.))break;}
  return {q[0]*rg*t,u,v,w,t,k,omega};
 }
 Q5 physical(const W7&w,const D3&n)const{auto q=conserved(w);double un=w[1]*n[0]+w[2]*n[1]+w[3]*n[2];return {q[0]*un,q[1]*un+w[0]*n[0],q[2]*un+w[0]*n[1],q[3]*un+w[0]*n[2],(q[4]+w[0])*un};}
 B5 jacobian(const W7&w,const D3&n)const{
  B5 a{};double r=rho(w),un=w[1]*n[0]+w[2]*n[1]+w[3]*n[2],cv=cp(w[4])-rg;auto q=conserved(w);
  for(int j=0;j<5;++j){double dr=j==0?1.:0.,dE=j==4?1.:0.;D3 du{};
   for(int k=0;k<3;++k)du[k]=((j==k+1?1.:0.)-w[k+1]*dr)/r;
   double de=(dE-q[4]/r*dr)/r;for(int k=0;k<3;++k)de-=w[k+1]*du[k];double dt=de/cv,dp=rg*(w[4]*dr+r*dt),dun=dot3(du,n),dm=dr*un+r*dun;
   a[j]=dm;for(int k=0;k<3;++k)a[5*(k+1)+j]=w[k+1]*dm+r*un*du[k]+dp*n[k];a[20+j]=(dE+dp)*un+(q[4]+w[0])*dun;
  }return a;
 }
 // Derivatives of velocities and temperature with respect to Q.
 B5 primitive_jacobian(const W7&w)const{
  B5 d{};double r=rho(w),cv=cp(w[4])-rg;auto q=conserved(w);
  for(int j=0;j<5;++j){double dr=j==0?1.:0.;D3 du{};for(int k=0;k<3;++k){du[k]=((j==k+1?1.:0.)-w[k+1]*dr)/r;d[5*(k+1)+j]=du[k];}
   double de=((j==4?1.:0.)-q[4]/r*dr)/r;for(int k=0;k<3;++k)de-=w[k+1]*du[k];d[20+j]=de/cv;d[j]=j==0?1.:0.;
  }return d;
 }
 Q5 pressure_jacobian(const W7&w)const{auto d=primitive_jacobian(w);Q5 p{};double r=rho(w);for(int j=0;j<5;++j)p[j]=rg*(w[4]*(j==0?1.:0.)+r*d[20+j]);return p;}
 B5 preconditioning_matrix(const W7&w,double uref2)const{
  B5 p{};for(int j=0;j<5;++j)p[6*j]=1.;double alpha=std::max(0.,1/uref2-1/std::pow(sound(w),2));auto dp=pressure_jacobian(w);Q5 v={1,w[1],w[2],w[3],enthalpy(w)};
  for(int row=0;row<5;++row)for(int col=0;col<5;++col)p[5*row+col]+=alpha*v[row]*dp[col];return p;
 }
 double preconditioned_radius(const W7&w,const D3&n,double uref2)const{double a2=std::pow(sound(w),2),beta=std::clamp(uref2/a2,1e-20,1.),un=w[1]*n[0]+w[2]*n[1]+w[3]*n[2],u=.5*(1+beta)*un;return std::abs(u)+std::sqrt(std::pow(.5*(1-beta)*un,2)+uref2);}
 // Static Roe mass prototype. Its registered callback was not invoked in
 // the observed normal step; it is not the production scalar mass path.
 double coupled_roe_mass(const W7&l,const W7&r,const D3&n,double a2_l,double a2_r,double ref_l,double ref_r)const{
  double un=.5*((l[1]+r[1])*n[0]+(l[2]+r[2])*n[1]+(l[3]+r[3])*n[2]),rho_up=rho(un>=0?l:r);
  double alpha=.5*(1-.5*(ref_l/a2_l+ref_r/a2_r)),u2=std::max(.5*(ref_l+ref_r),1e-30),c=std::sqrt(alpha*alpha*un*un+u2),lp=std::abs((1-alpha)*un+c),lm=std::abs((1-alpha)*un-c),dl=(lp-lm)/(2*c),spread=.5*(lp+lm)-std::abs(un);
  double du=(l[1]-r[1])*n[0]+(l[2]-r[2])*n[1]+(l[3]-r[3])*n[2],dp=l[0]-r[0];
  return rho_up*un+.5*(dl*rho_up*du+dp/u2*(2*alpha*std::abs(un)-alpha*dl*un+spread));
 }
 Q5 ausm(const W7&l,const W7&r,const D3&n)const{
  const double a=.5*(sound(l)+sound(r)),ul=l[1]*n[0]+l[2]*n[1]+l[3]*n[2],ur=r[1]*n[0]+r[2]*n[1]+r[3]*n[2],ml=ul/a,mr=ur/a;
  auto mp=[](double m){return std::abs(m)>=1?.5*(m+std::abs(m)):.25*(m+1)*(m+1);};
  auto mm=[](double m){return std::abs(m)>=1?.5*(m-std::abs(m)):-.25*(m-1)*(m-1);};
  auto pp=[](double m){return std::abs(m)>=1?(m>=0?1.:0.):.25*(m+1)*(m+1)*(2-m);};
  auto pm=[](double m){return std::abs(m)>=1?(m<=0?1.:0.):.25*(m-1)*(m-1)*(2+m);};
  double speed=a*(mp(ml)+mm(mr)),pressure=pp(ml)*l[0]+pm(mr)*r[0];const W7&up=speed>=0?l:r;double mass=rho(up)*speed;
  return {mass,mass*up[1]+pressure*n[0],mass*up[2]+pressure*n[1],mass*up[3]+pressure*n[2],mass*enthalpy(up)};
 }
 // Independent implementation of the observed M4/P5 all-speed AUSM algebra.
 // Cached cell acoustics and local reference velocities remain distinct from
 // reconstructed face thermodynamics. Inputs and outputs are nondimensional.
 Q5 coupled_ausm(const W7&l,const W7&r,const D3&n,double cell_a2_l,double cell_a2_r,double ref_l,double ref_r,const std::array<double,4>*face_thermo=nullptr)const{
  double a2=.5*(cell_a2_l+cell_a2_r),a=std::sqrt(a2),ul=l[1]*n[0]+l[2]*n[1]+l[3]*n[2],ur=r[1]*n[0]+r[2]*n[1]+r[3]*n[2];
  double beta=std::clamp(2*std::max(ref_l,ref_r)/a2,1e-20,1.),ml=ul/a,mr=ur/a,average=.5*(std::abs(ml)+std::abs(mr));
  double factor=std::sqrt(4*beta+(1-beta)*(1-beta)*average*average)/(1+beta),c=a*factor;
  double left=((1+beta)*ml+(1-beta)*mr)/(2*factor),right=((1-beta)*ml+(1+beta)*mr)/(2*factor);
  auto mp=[](double m){return std::abs(m)>=1?.5*(m+std::abs(m)):.25*(m+1)*(m+1)*(1+.5*(m-1)*(m-1));};
  auto mm=[](double m){return std::abs(m)>=1?.5*(m-std::abs(m)):-.25*(m-1)*(m-1)*(1+.5*(m+1)*(m+1));};
  auto pp=[](double m){return std::abs(m)>=1?(m>=0?1.:0.):.25*(m+1)*(m+1)*(2-m)+.1875*m*(m*m-1)*(m*m-1);};
  auto pm=[](double m){return std::abs(m)>=1?(m<=0?1.:0.):.25*(m-1)*(m-1)*(2+m)-.1875*m*(m*m-1)*(m*m-1);};
  double mplus=mp(left),mminus=mm(right),mface=mplus+mminus,rl=face_thermo?(*face_thermo)[0]:rho(l),rr=face_thermo?(*face_thermo)[1]:rho(r);
  double mass=c*mface*(mface>=0?rl:rr)+(1/beta-1)*c*(l[0]-r[0])*(-.5*(left+std::abs(left))+.5*(right-std::abs(right))+mplus-mminus)/(l[0]/rl+r[0]/rr);
  double pplus=pp(left),pminus=pm(right),pressure=l[0]*pplus+r[0]*pminus+.5*(rl+rr)*c*c*(left-right)*pplus*pminus;
  const W7&up=mass>=0?l:r;double H=face_thermo?(*face_thermo)[mass>=0?2:3]:enthalpy(up);return {mass,mass*up[1]+pressure*n[0],mass*up[2]+pressure*n[1],mass*up[3]+pressure*n[2],mass*H};
 }
};

inline double coupled_wall_pressure(const W7&face,const D3&normal,double reference2,double density){
 double un=0;for(int j=0;j<3;++j)un+=face[j+1]*normal[j];return face[0]+density*std::sqrt(reference2)*un;
}
inline Q5 coupled_wall_pressure_jacobian(const Gas&gas,const W7&cell,const W7&face,const D3&normal,double reference2,double density,bool independent_density){
 auto primitive=gas.primitive_jacobian(cell);auto dp=gas.pressure_jacobian(cell);Q5 result{};double un=0;for(int k=0;k<3;++k)un+=face[k+1]*normal[k];
 for(int j=0;j<5;++j){double dun=0;for(int k=0;k<3;++k)dun+=normal[k]*primitive[5*(k+1)+j];double drho=independent_density?(j==0?1.:0.):dp[j]/(gas.rg*face[4])-density/face[4]*primitive[20+j];result[j]=dp[j]+std::sqrt(reference2)*(drho*un+density*dun);}return result;
}
struct Wall {
 double shear,omega,production,velocity_coefficient,yplus,omega_outward_derivative,ustar;
};
struct OmegaSource {double value,negative_derivative,cross;};
// The supplied Case uses coeff-kw-mu-dok=.001. Runtime source observations
// additionally select |Dw| <= 1e6*Yw and linearize both signs of Dw.
// The formula is homogeneous: it accepts physical or consistent solver units.
inline OmegaSource sst_omega_source(double rho,double mut,double mu,double pk,double omega,double beta,double alpha,double raw_cross,double laminar_factor=.001,double cross_limit=1e6,bool linearize_negative=true){
 double destruction=beta*rho*omega*omega,cross=raw_cross;
 if(cross_limit>0)cross=std::clamp(cross,-cross_limit*destruction,cross_limit*destruction);
 double production=alpha*rho/std::max(mut+laminar_factor*mu,1e-30)*pk;
 double derivative=2*beta*rho*omega+(linearize_negative?std::abs(cross):std::max(cross,0.))/std::max(omega,1e-30);
 return {production-destruction+cross,derivative,cross};
}
// Fluent's documented smooth correlation treatment. k has zero normal flux.
inline Wall wall_model(double rho,double mu,double k,double omega,double speed,double d){
 d=std::max(d,1e-14);const double nu=mu/rho,kappa=.4187,beta=.072,betaStar=.09,E=9.793;
 double ustar=std::sqrt(std::max(0.,nu*std::abs(speed)/d+std::sqrt(betaStar)*k));
 double yplus=std::max(d*ustar/nu,1e-15),logp=std::log(E*std::max(yplus,.2))/kappa;
 double inverse_uplus=std::pow(std::pow(yplus,-4)+std::pow(logp,-4),.25);
 double shear=rho*speed*inverse_uplus*ustar,coefficient=rho*inverse_uplus*ustar;
 double wpLog=1/(std::sqrt(betaStar)*kappa*std::max(yplus,.2));
 double lamOmega=2*nu/(beta*d*d),logOmega=ustar*ustar/nu*wpLog;
 double wallOmega=std::pow(std::pow(lamOmega,1.3)+std::pow(logOmega,1.3),1/1.3);
 double pl=rho*k/std::max(omega,1e-20)*std::pow(shear/mu,2);
 double pt=yplus>.2?shear*shear/(mu*kappa*yplus):1e100;
 double pk=pl>0&&pt>0?pl*pt/(pl+pt):0;
 double ratio=(2*std::pow(lamOmega,1.3)+std::pow(logOmega,1.3))/(std::pow(lamOmega,1.3)+std::pow(logOmega,1.3));
 return {shear,wallOmega,pk,coefficient,yplus,wallOmega*ratio/d,ustar};
}
// Smooth-wall face omega is an extrapolation using the blended gradient
// magnitudes. It is distinct from differentiating the blended cell constraint.
inline double smooth_wall_omega_face(double rho,double mu,double ustar,double d,double cell_omega){
 if(!(rho>0&&mu>0&&d>0&&ustar>=0&&cell_omega>0))throw std::runtime_error("invalid smooth wall omega face inputs");
 double y=rho*d*ustar/mu,laminar=-4/(.072*std::pow(y+1e-12,3)),logarithmic=-1/(std::sqrt(.09)*.4187*std::pow(std::max(y,.2),2));
 return cell_omega-laminar*std::pow(1+std::pow(logarithmic/laminar,1.3),1/1.3)*y*rho*ustar*ustar/mu;
}
}
