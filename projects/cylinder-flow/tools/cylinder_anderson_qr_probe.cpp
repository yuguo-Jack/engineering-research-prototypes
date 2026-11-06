// Independent factorization checks and optional small NumPy-oracle snapshots.
#define main cylinder_driver_main
#include "../src/cylinder.cpp"
#undef main

int main(int argc,char**argv){try{
 std::ofstream samples;if(argc==2){samples.open(argv[1]);if(!samples)throw std::runtime_error("cannot create QR samples");
  samples<<"sample,dof,columns,step,answer";for(int j=0;j<8;++j)samples<<",df"<<j;for(int j=0;j<8;++j)samples<<",dg"<<j;samples<<'\n'<<std::setprecision(17);
 }else if(argc!=1)throw std::runtime_error("usage: cylinder_anderson_qr_probe [samples.csv]");
 double maximum_factorization=0,maximum_orthogonality=0;int proposals=0,snapshot=0,total_downdates=0,total_incremental=0,total_rebuilds=0;
 for(int depth:{2,3,6,8})for(int dependent=0;dependent<2;++dependent){Anderson aa(depth);
  for(int iteration=0;iteration<64;++iteration){VV q(4),step(4);
   for(int i=0;i<4;++i)for(int k=0;k<4;++k){int dof=4*i+k+1;
    q[i][k]=std::sin((iteration+1)*dof*.071)+.01*iteration*std::cos(dof*.73);
    step[i][k]=dependent?std::pow(.9,iteration)*std::sin(dof*.37)+1e-13*std::sin(iteration*dof*.51):
        std::sin((iteration+1)*dof*.173)+.2*std::cos((iteration+2)*dof*.41);
   }
   VV answer=aa.propose(q,step);++proposals;int n=static_cast<int>(aa.basis.size());
   for(int i=0;i<n;++i)for(int j=0;j<n;++j){double error=std::abs(dot(aa.basis[i],aa.basis[j])-(i==j?1.:0.));maximum_orthogonality=std::max(maximum_orthogonality,error);}
   for(int j=0;j<n;++j){VV reconstructed(4);for(int i=0;i<=j;++i)axpy(reconstructed,aa.r[i][j],aa.basis[i]);
    axpy(reconstructed,-1,aa.df[aa.active[j]]);double error=norm(reconstructed)/std::max(norm(aa.df[aa.active[j]]),1e-30);maximum_factorization=std::max(maximum_factorization,error);
   }
   if(!dependent&&samples&&aa.qr_valid&&!aa.df.empty()&&iteration%8==7){int m=static_cast<int>(aa.df.size());
    for(int i=0;i<4;++i)for(int k=0;k<4;++k){samples<<snapshot<<','<<4*i+k<<','<<m<<','<<step[i][k]<<','<<answer[i][k];
     for(int j=0;j<8;++j)samples<<','<<(j<m?aa.df[j][i][k]:0.);for(int j=0;j<8;++j)samples<<','<<(j<m?aa.dg[j][i][k]:0.);samples<<'\n';
    }++snapshot;
   }
   if(iteration==23||iteration==47)aa.restart();
  }
  total_downdates+=aa.downdates;total_incremental+=aa.incremental_columns;total_rebuilds+=aa.full_rebuilds;
 }
 if(maximum_factorization>5e-10||maximum_orthogonality>5e-10||total_downdates==0||total_incremental==0)throw std::runtime_error("QR factorization/orthogonality or actual incremental coverage failed");
 std::cout<<std::setprecision(17)<<"{\"status\":\"sliding_QR_factorization_checks_passed\",\"proposals\":"<<proposals
  <<",\"max_factorization_relative\":"<<maximum_factorization<<",\"max_orthogonality\":"<<maximum_orthogonality
  <<",\"downdates\":"<<total_downdates<<",\"incremental_columns\":"<<total_incremental<<",\"full_rebuilds\":"<<total_rebuilds
  <<",\"numpy_snapshots\":"<<snapshot<<",\"integration_steps\":0}\n";
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
