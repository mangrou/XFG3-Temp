// For XFG3 LSV implementation illustration

#ifndef quantlib_heston_dupire_process_hpp
#define quantlib_heston_dupire_process_hpp

#include <ql/stochasticprocess.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/quote.hpp>


namespace QuantLib {

	class LeverageSurface;


    class HestonDupireProcess : public StochasticProcess 
	{
      public:
        enum Discretization { QuadraticExponential};

        HestonDupireProcess(const Handle<YieldTermStructure>& riskFreeRate,
                      const Handle<YieldTermStructure>& dividendYield,
                      const Handle<Quote>& s0,
                      Real v0, Real kappa,
                      Real theta, Real sigma, Real rho,
					  const boost::shared_ptr<LeverageSurface> leverSurf,
                      Discretization d = QuadraticExponential);

        Size size() const;
        Disposable<Array> evolve(Time t0, const Array& x0,
                                 Time dt, const Array& dw) const;

        Real v0()    const { return v0_; }
        Real rho()   const { return rho_; }
        Real kappa() const { return kappa_; }
        Real theta() const { return theta_; }
        Real sigma() const { return sigma_; }

        const Handle<Quote>& s0() const;
        const Handle<YieldTermStructure>& dividendYield() const;
        const Handle<YieldTermStructure>& riskFreeRate() const;

		boost::shared_ptr<LeverageSurface> getLeverSurf() const { return leverSurf_; }

      private:
        Handle<YieldTermStructure> riskFreeRate_, dividendYield_;
        Handle<Quote> s0_;
        Real v0_, kappa_, theta_, sigma_, rho_;
        Discretization discretization_;

		const boost::shared_ptr<LeverageSurface> leverSurf_;
    };

}


#endif
