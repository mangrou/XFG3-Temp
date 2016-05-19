// for XFG3 LSV implementation illustration

#ifndef quantlib_fdm_hestondupire_op_hpp
#define quantlib_fdm_hestondupire_op_hpp

#include <ql/methods/finitedifferences/operators/fdmlinearopcomposite.hpp>

namespace QuantLib {
    class FdmMesher;
	class HestonDupireProcess;
	class YieldTermStructure;
	class LeverageSurface;
	class ModNinePointLinearOp;



	class FdmHestonDupireVarianceOp
	{
	  public:
		FdmHestonDupireVarianceOp( const boost::shared_ptr<FdmMesher>& mesher, 
			const boost::shared_ptr<HestonDupireProcess>& process );

		void setTime(Time t1, Time t2) const;
		const boost::shared_ptr<TripleBandLinearOp> getMap() const { return mapT_; }

	  protected:
		const boost::shared_ptr<TripleBandLinearOp>  dyMap_;
		boost::shared_ptr<TripleBandLinearOp> mapT_;

		const boost::shared_ptr<YieldTermStructure> rTS_;
	};



	class FdmHestonDupireEquityOp
	{
	  public:
		FdmHestonDupireEquityOp( const boost::shared_ptr<FdmMesher>& mesher, 
			const boost::shared_ptr<HestonDupireProcess>& process);

        void setTime(Time t1, Time t2) const;
		const boost::shared_ptr<TripleBandLinearOp> getMap() const { return mapT_; }

	  protected:

		Array equityArray_;
        Array varianceValues_;

		const boost::shared_ptr<TripleBandLinearOp>  dxMap_;
		const boost::shared_ptr<TripleBandLinearOp> dxxMap_;
		boost::shared_ptr<TripleBandLinearOp> mapT_;

        const boost::shared_ptr<FdmMesher> mesher_;
		const boost::shared_ptr<YieldTermStructure> rTS_;
        const boost::shared_ptr<YieldTermStructure> qTS_;
		const boost::shared_ptr<LeverageSurface> leverSurf_;
	};

	class FdmHestonDupireCorrelationOp
	{
	  public:
		FdmHestonDupireCorrelationOp( const boost::shared_ptr<FdmMesher>& mesher, 
			const boost::shared_ptr<HestonDupireProcess>& process);
		void setTime(Time t1, Time t2) const;
		const boost::shared_ptr<ModNinePointLinearOp> getMap() const { return mapT_; }

	  protected:
		Array equityArray_;

		const boost::shared_ptr<FdmMesher> mesher_;
		const boost::shared_ptr<LeverageSurface> leverSurf_;

		const boost::shared_ptr<ModNinePointLinearOp> dxyMap_;
		boost::shared_ptr<ModNinePointLinearOp> mapT_;
	};

	class FdmHestonDupireOp : public FdmLinearOpComposite 
	{
	  public:
		FdmHestonDupireOp( const boost::shared_ptr<FdmMesher>& mesher,
            const boost::shared_ptr<HestonDupireProcess>& process);

		Size size() const;
        void setTime(Time t1, Time t2) const;

        Disposable<Array> apply(const Array& r) const;
        Disposable<Array> apply_mixed(const Array& r) const;

        Disposable<Array> apply_direction(Size direction, const Array& r) const;
        Disposable<Array> solve_splitting(Size direction, const Array& r, Real a) const;

	  private:
		boost::shared_ptr<FdmHestonDupireCorrelationOp> correlationMap_;
		boost::shared_ptr<FdmHestonDupireEquityOp> xMap_;
		boost::shared_ptr<FdmHestonDupireVarianceOp> yMap_;
	};

}

#endif
