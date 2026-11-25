#include <Rcpp.h>
#include <algorithm>
#include <utility>
#include <vector>
using namespace Rcpp;

// Reapply genome-wide BH at the original threshold, then intersect with RBP.
// Every gene is evaluated jointly across all three comparisons in each draw.
// [[Rcpp::export]]
NumericMatrix bootstrap_counts(List delta, NumericMatrix beta, NumericMatrix se,
                               LogicalMatrix valid, LogicalVector rbp,
                               double alpha, double lfc_cut) {
  NumericMatrix dn=delta[0];
  const int G=beta.nrow(),B=dn.ncol(),K=beta.ncol();
  NumericMatrix out(B,3*K);
  const double zcut=-R::qnorm(alpha/2.0,0.0,1.0,true,false);
  for(int k=0;k<K;++k) {
    NumericMatrix dt=delta[k+1];
    int m=0;for(int g=0;g<G;++g) if(valid(g,k)==TRUE) ++m;
    std::vector<double> effect(G);
    std::vector<std::pair<double,int> > pv;
    pv.reserve(G);
    for(int b=0;b<B;++b) {
      pv.clear();
      for(int g=0;g<G;++g) {
        if(valid(g,k)!=TRUE) continue;
        const double e=beta(g,k)+dt(g,b)-dn(g,b);
        effect[g]=e;
        const double z=std::abs(e/se(g,k));
        if(z>=zcut) {
          const double p=2.0*R::pnorm(-z,0.0,1.0,true,false);
          if(p<=alpha) pv.emplace_back(p,g);
        }
      }
      std::sort(pv.begin(),pv.end());
      int reject=0;
      for(int i=0;i<(int)pv.size();++i)
        if(pv[i].first<alpha*(i+1.0)/m) reject=i+1;
      int up=0,down=0;
      for(int i=0;i<reject;++i) {
        const int g=pv[i].second;
        if(rbp[g]!=TRUE) continue;
        if(effect[g]>=lfc_cut) ++up;
        if(effect[g]<=-lfc_cut) ++down;
      }
      out(b,3*k)=up;out(b,3*k+1)=down;out(b,3*k+2)=reject;
    }
  }
  return out;
}
