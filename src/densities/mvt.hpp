/*
 Multivariate t distribution with user supplied scale matrix
Class to evaluate the negative log density of a mean zero
multivariate t distributed variable with general covariance matrix Sigma.
Intended for small dense covariance matrices.
*/
using namespace density;

template <class Type>
class MVT_t: public MVNORM_t<Type>
{
  Type df;

public:
  MVT_t(Type df_)
    : MVNORM_t<Type>()
    {
      df = df_;
    }
  MVT_t(matrix<Type> Sigma_, Type df_)
    : MVNORM_t<Type>(Sigma_)
    {
      df = df_;
    }

  Type getdf(){
    return df;
  }
  void setdf(Type df_){
    df = df_;
  }

  /** \brief Evaluate the negative log density */
  Type operator()(vector<Type> x){
    Type p = x.size();
    //Lange et al. 1989 http://www.jstor.org/stable/2290063
    return -lgamma(Type(0.5)*(df+p))+lgamma(Type(0.5)*df)+p*Type(0.5)*log(df)+p*lgamma(Type(0.5))-Type(0.5)*this->logdetQ + Type(0.5)*(df+p)*log(Type(1.0)+this->Quadform(x)/df);

  }
};
