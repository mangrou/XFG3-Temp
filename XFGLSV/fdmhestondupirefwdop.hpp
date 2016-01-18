// for XFG3 LSV implementation illustration

#ifndef quantlib_fdm_hestondupire_fwd_op_hpp
#define quantlib_fdm_hestondupire_fwd_op_hpp

#include <ql/methods/finitedifferences/operators/fdmlinearopcomposite.hpp>

namespace QuantLib {
    class FdmMesher;
	class HestonDupireProcess;
	class YieldTermStructure;
    class ModTripleBandLinearOp;
	class ModNinePointLinearOp;
	class LeverageSurface;


	class FdmHestonDupireVarianceFwdOp
	{
	  public:
		FdmHestonDupireVarianceFwdOp( const boost::shared_ptr<FdmMesher>& mesher, 
			Real kappa, Real theta, Real sigma, Real v0 );
		void setTime(Time t1, Time t2) const;
		boost::shared_ptr<ModTripleBandLinearOp> getMap() const { return mapT_; }

	  private:
        const Real v0_, kappa_, theta_, sigma_, a_;

		Array varianceValues_;

		const boost::shared_ptr<FdmMesher> mesher_;

		const boost::shared_ptr<ModTripleBandLinearOp> dvMap_;
		const boost::shared_ptr<ModTripleBandLinearOp> dvvMap_;

		boost::shared_ptr<ModTripleBandLinearOp> mapT_;
	};

	class FdmHestonDupireEquityFwdOp
	{
	  public:
		FdmHestonDupireEquityFwdOp( const boost::shared_ptr<FdmMesher>& mesher, 
			const boost::shared_ptr<HestonDupireProcess>& process );

        void setTime(double t1, double t2, bool isExplicit=true) const;

		boost::shared_ptr<ModTripleBandLinearOp> getMap() const { return mapT_; }

	  protected:

		const Real v0_;

		Array varianceValues_;

		const boost::shared_ptr<FdmMesher> mesher_;

		const boost::shared_ptr<LeverageSurface> leverSurf_;

		const boost::shared_ptr<YieldTermStructure> rTS_;
        const boost::shared_ptr<YieldTermStructure> qTS_;

		const boost::shared_ptr<ModTripleBandLinearOp> dxMap_;
		const boost::shared_ptr<ModTripleBandLinearOp> dxxMap_;

		boost::shared_ptr<ModTripleBandLinearOp> mapT_;
	};

	class FdmHestonDupireCorrelationFwdOp
	{
	  public:
		FdmHestonDupireCorrelationFwdOp( const boost::shared_ptr<FdmMesher>& mesher, 
			const boost::shared_ptr<HestonDupireProcess>& process );

		void setTime(double t1, double t2) const;

		boost::shared_ptr<ModNinePointLinearOp> getMap() const { return mapT_; }

	  protected:

		const boost::shared_ptr<FdmMesher> mesher_;
		const boost::shared_ptr<LeverageSurface> leverSurf_;

		const boost::shared_ptr<ModNinePointLinearOp> dxyMap_;
		boost::shared_ptr<ModNinePointLinearOp> mapT_;
	};

	class FdmHestonDupireFwdOp : public FdmLinearOpComposite {
      public:
        FdmHestonDupireFwdOp(
            const boost::shared_ptr<FdmMesher>& mesher,
            const boost::shared_ptr<HestonDupireProcess>& process);

		Size size() const;
        void setTime(Time t1, Time t2) const;

        Disposable<Array> apply(const Array& r) const;
        Disposable<Array> apply_mixed(const Array& r) const;

        Disposable<Array> apply_direction(Size direction,
                                          const Array& r) const;
        Disposable<Array> solve_splitting(Size direction,
                                          const Array& r, Real s) const;

	  private:
        const Real kappa_, theta_, sigma_, rho_, v0_;
		mutable Time t1_, t2_;

		const boost::shared_ptr<FdmMesher> mesher_;
		boost::shared_ptr<FdmMesher> logSpotMesher_;
		const boost::shared_ptr<LeverageSurface> leverSurf_;

		const boost::shared_ptr<YieldTermStructure> rTS_;
        const boost::shared_ptr<YieldTermStructure> qTS_;

		boost::shared_ptr<FdmHestonDupireCorrelationFwdOp> correlationMap_;
		boost::shared_ptr<FdmHestonDupireEquityFwdOp> xMap_;
		boost::shared_ptr<FdmHestonDupireVarianceFwdOp> yMap_;
	};

}

#endif
