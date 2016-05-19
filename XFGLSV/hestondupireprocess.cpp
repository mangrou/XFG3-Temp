// For XFG3 LSV implementation illustration

#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/math/distributions/chisquaredistribution.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/processes/eulerdiscretization.hpp>
#include <ql/experimental/finitedifferences/hestondupireprocess.hpp>

#include <boost/shared_array.hpp>
#include <ql/experimental/finitedifferences/leveragesurface.hpp>

namespace QuantLib {

    HestonDupireProcess::HestonDupireProcess(
                              const Handle<YieldTermStructure>& riskFreeRate,
                              const Handle<YieldTermStructure>& dividendYield,
                              const Handle<Quote>& s0,
                              Real v0, Real kappa,
                              Real theta, Real sigma, Real rho,
							  const boost::shared_ptr<LeverageSurface> leverSurf,
                              Discretization d)
    : StochasticProcess(boost::shared_ptr<discretization>(new EulerDiscretization)),
      riskFreeRate_(riskFreeRate), dividendYield_(dividendYield), s0_(s0),
      v0_(v0), kappa_(kappa), theta_(theta), sigma_(sigma), rho_(rho),
	  leverSurf_(leverSurf),
      discretization_(d) 
	{
        registerWith(riskFreeRate_);
        registerWith(dividendYield_);
        registerWith(s0_);
    }

    Size HestonDupireProcess::size() const 
	{
        return 2;
    }

    Disposable<Array> HestonDupireProcess::evolve(Time t0, const Array& x0,
                                            Time dt, const Array& dw) const 
	{
        Array retVal(2);
        Real mu;

        const Real sdt = std::sqrt(dt);
        const Real sqrhov = std::sqrt(1.0 - rho_*rho_);

		Real leverageVal = leverSurf_->localVol(t0, x0[0]);;

        switch (discretization_) 
		{
          case QuadraticExponential:
          {
            // for details of the quadratic exponential discretization scheme
            // see Leif Andersen,
            // Efficient Simulation of the Heston Stochastic Volatility Model
		    // For local freezing of lsv leverage
		    // see Anthonie W Van der Stoep
		    // The heston stochastic local volatility model: efficient monte carlo simulation
            const Real ex = std::exp(-kappa_*dt);

            const Real m  =  theta_+(x0[1]-theta_)*ex;
            const Real s2 =  x0[1]*sigma_*sigma_*ex/kappa_*(1-ex)
                           + theta_*sigma_*sigma_/(2*kappa_)*(1-ex)*(1-ex);
            const Real psi = s2/(m*m);

            const Real g1 =  0.5;
            const Real g2 =  0.5;
                  Real k0 = -rho_*kappa_*theta_*dt/sigma_;
            const Real k1 =  g1*dt*(kappa_*rho_/sigma_-0.5)-rho_/sigma_;
            const Real k2 =  g2*dt*(kappa_*rho_/sigma_-0.5)+rho_/sigma_;
            const Real k3 =  g1*dt*(1-rho_*rho_);
            const Real k4 =  g2*dt*(1-rho_*rho_);
            const Real A  =  k2+0.5*k4;

            if (psi < 1.5) 
			{
                const Real b2 = 2/psi-1+std::sqrt(2/psi*(2/psi-1));
                const Real b  = std::sqrt(b2);
                const Real a  = m/(1+b2);

                retVal[1] = a*(b+dw[1])*(b+dw[1]);
            }
            else 
			{
                const Real p = (psi-1)/(psi+1);
                const Real beta = (1-p)/m;

                const Real u = CumulativeNormalDistribution()(dw[1]);

                retVal[1] = ((u <= p) ? 0.0 : std::log((1-p)/(1-u))/beta);
            }

            mu =   riskFreeRate_->forwardRate(t0, t0+dt, Continuous)
                 - dividendYield_->forwardRate(t0, t0+dt, Continuous);

			// local freezing part
			double c0 = leverageVal*leverageVal*x0[1]*dt;
			double k5 = rho_*(retVal[1]-x0[1]+kappa_*dt*x0[1]-kappa_*dt*theta_)/sigma_;

            retVal[0] = x0[0]*std::exp(mu*dt - 0.5*c0 + k5 + sqrhov*sqrt(c0)*dw[0]);
          }
          break;
          default:
            QL_FAIL("unknown discretization schema");
        }

        return retVal;
    }

    const Handle<Quote>& HestonDupireProcess::s0() const 
	{
        return s0_;
    }

    const Handle<YieldTermStructure>& HestonDupireProcess::dividendYield() const 
	{
        return dividendYield_;
    }

    const Handle<YieldTermStructure>& HestonDupireProcess::riskFreeRate() const 
	{
        return riskFreeRate_;
    }

}
