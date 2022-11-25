#ifndef __sliceSampling__
#define __sliceSampling__

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// truncated samplings
double rtbeta_cpp(const double &alpha, const double &beta, const double &a, const double &b);
double rtgamma_cpp(const double &shape, const double &scale, const double &a, const double &b);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Confluent Hypergeometric, v~CH(a, b,s)
// p(v) \propto v^(a-1) (1-v)^(b-1) exp(-sv)
void rCH_void(double& v_, double& l_, 
              const double& aa, const double& bb, const double& ss);
arma::vec rCH_cpp(int nsamp, int burnin,
                  const double& aa, const double& bb, const double& ss);

// Gaussian Hypergeometric, v~GH(a, b, x, z)
// p(v) \propto v^(a-1) (1-v)^(b-1) (1 + xv)^(-z)
void rGH_void(double& v_, double& t_, double& l_,
              const double& aa, const double& bb, const double& xx, const double& zz);
arma::vec rGH_cpp(int nsamp, int burnin,
                  const double& aa, const double& bb, const double& xx, const double& zz);

// truncated Compound Confluent Hypergeometric, u~tCCH(a, b, z, s, nu, theta)
// p(u) \propto u^(a-1) (1-nu u)^(b-1) (theta + (1-theta)nu u)^(-z) exp(-s u)
// Compound Confluent Hypergeometric, v~tCCH(a, b, z, tilde_s, x) (tilde_s = s/nu, x = 1/theta-1)
// p(v) \propto v^(a-1) (1-v)^(b-1) (1 + x v)^(-z) exp(-tilde_s u)
void rCCH_void(double& v_, double& t_, double& l_,double& m_,
               const double& aa, const double& bb, const double& zz, const double& ss, const double& xx);
arma::vec rtCCH_cpp(int nsamp, int burnin,
                    const double& aa, const double& bb, const double& zz, const double& ss, const double &nu, const double& theta);

// Appell distribution, u~APL(a,b,z,w,x,y)
// p(u) \propto u^(a-1) (1-u)^(b-1) (1+xu)^(-z) (1+yu)^(-w)
void rAPL_void(double& v_, double& t_, double& q_, double& l_,double& m_,
               const double& aa, const double& bb, const double& zz, const double& ww, const double& xx, const double& yy);
arma::vec rAPL_cpp(int nsamp, int burnin,
                   const double& aa, const double& bb, const double& zz, const double& ww, const double& xx, const double& yy);

#endif // __sliceSampling__