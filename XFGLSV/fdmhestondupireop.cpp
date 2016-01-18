// for XFG3 LSV implementation illustration

#include <ql/math/functional.hpp>
#include <ql/methods/finitedifferences/meshers/fdmmesher.hpp>
#include <ql/methods/finitedifferences/meshers/fdmmeshercomposite.hpp>
#include <ql/methods/finitedifferences/meshers/Predefined1dMesher.hpp>
#include <ql/methods/finitedifferences/operators/fdmlinearoplayout.hpp>
#include <ql/methods/finitedifferences/operators/firstderivativeop.hpp>
#include <ql/methods/finitedifferences/operators/secondderivativeop.hpp>
#include <ql/methods/finitedifferences/operators/secondordermixedderivativeop.hpp>

#include <ql/experimental/finitedifferences/hestondupireprocess.hpp>
#include <ql/experimental/finitedifferences/fdmhestondupireop.hpp>
#include <ql/experimental/finitedifferences/modninepointlinearop.hpp>
#include <ql/experimental/finitedifferences/leveragesurface.hpp>

namespace QuantLib {

	


	// ---------------------------   FdmHestonDupireVarianceOp    -----------------

	FdmHestonDupireVarianceOp::FdmHestonDupireVarianceOp( const boost::shared_ptr<FdmMesher>& mesher, 
			const boost::shared_ptr<HestonDupireProcess>& process )
	: dyMap_( new TripleBandLinearOp(SecondDerivativeOp(1, mesher)
	.mult(0.5*process->sigma()*process->sigma()*mesher->locations(1))
             .add(FirstDerivativeOp(1, mesher)
			 .mult(process->kappa()*(process->theta() - mesher->locations(1))))) ),
	mapT_( new TripleBandLinearOp(1, mesher) ),
	rTS_(process->riskFreeRate().currentLink())
	{
		;
	}

	void FdmHestonDupireVarianceOp::setTime(Time t1, Time t2) const
	{
		const Rate r = rTS_->forwardRate(t1, t2, Continuous).rate();

		mapT_->axpyb( Array(), *dyMap_, *dyMap_, Array(1,-0.5*r) );
	}




	// ---------------------------   FdmHestonDupireEquityOp    -------------------

	FdmHestonDupireEquityOp::FdmHestonDupireEquityOp( const boost::shared_ptr<FdmMesher>& mesher, 
			const boost::shared_ptr<HestonDupireProcess>& process)
	: varianceValues_( mesher->locations(1) * 0.5 ),
	  equityArray_( mesher->locations(0) ),
	  dxMap_( new FirstDerivativeOp(0, mesher)),
	  dxxMap_( new TripleBandLinearOp(SecondDerivativeOp(0, mesher).mult( mesher->locations(1)*0.5 )) ),
	  mapT_( new TripleBandLinearOp(0, mesher) ),
	  mesher_(mesher),
	  leverSurf_(process->getLeverSurf()),
	  rTS_(process->riskFreeRate().currentLink()),
      qTS_(process->dividendYield().currentLink())
	{
		// on the boundary s_min and s_max the second derivative
		// d^2V/dS^2 is zero and due to Ito's Lemma the variance term
		// in the drift should vanish.
		const boost::shared_ptr<FdmLinearOpLayout> layout = mesher_->layout();
		const FdmLinearOpIterator endIter = layout->end();
		for (FdmLinearOpIterator iter = layout->begin(); iter!=endIter; ++iter) 
		{
			if ( iter.coordinates()[0] == 0 || iter.coordinates()[0] == layout->dim()[0]-1) 
			{
				varianceValues_[iter.index()] = 0.0;
			}
		}

		for (Size i = 0; i < equityArray_.size(); i++)
		{
			equityArray_[i] =  exp(equityArray_[i]);
		}
	}

	void FdmHestonDupireEquityOp::setTime(Time t1, Time t2) const
	{
		const Rate r = rTS_->forwardRate(t1, t2, Continuous).rate();
        const Rate q = qTS_->forwardRate(t1, t2, Continuous).rate();

		Real t = (t1+t2) / 2;
		Array leverageArray(equityArray_);
		for (Size i = 0; i < equityArray_.size(); i++)
		{
			leverageArray[i] = pow(leverSurf_->localVol(t, equityArray_[i]), 2);
		}

		mapT_->axpyb( r-q-leverageArray*varianceValues_, *dxMap_, dxxMap_->mult(leverageArray), Array(1, -0.5*r) );
	}





	// --------------------------- FdmHestonDupireCorrelationOp -------------------

	FdmHestonDupireCorrelationOp::FdmHestonDupireCorrelationOp( const boost::shared_ptr<FdmMesher>& mesher, 
			const boost::shared_ptr<HestonDupireProcess>& process )
	: dxyMap_( new ModNinePointLinearOp(SecondOrderMixedDerivativeOp(0, 1, mesher).mult(mesher->locations(1)*process->rho()*process->sigma())) ),
	mapT_( new ModNinePointLinearOp(0, 1, mesher) ),
	mesher_(mesher),
	leverSurf_(process->getLeverSurf()),
	equityArray_(mesher->locations(0))
	{
		for (Size i = 0; i < equityArray_.size(); i++)
		{
			equityArray_[i] = exp(equityArray_[i]);
		}
	}

	void FdmHestonDupireCorrelationOp::setTime(double t1, double t2) const
	{
		double t = (t1 + t2) / 2;
		Array leverageArray(equityArray_);
		for (Size i = 0; i < equityArray_.size(); i++)
		{
			leverageArray[i] = leverSurf_->localVol(t,equityArray_[i]);
		}
		mapT_->amx(leverageArray, *dxyMap_);
	}


	// ---------------------------        FdmHestonDupireOp     -------------------

	FdmHestonDupireOp::FdmHestonDupireOp( const boost::shared_ptr<FdmMesher>& mesher,
            const boost::shared_ptr<HestonDupireProcess>& process )
	: correlationMap_( new FdmHestonDupireCorrelationOp(mesher, process) ),
	yMap_( new FdmHestonDupireVarianceOp(mesher, process) ),
	xMap_( new FdmHestonDupireEquityOp(mesher, process) )
	{
	}

	Disposable<Array> FdmHestonDupireOp::apply(const Array& u) const 
	{
		return xMap_->getMap()->apply(u) + yMap_->getMap()->apply(u) + correlationMap_->getMap()->apply(u);
	}

	Size FdmHestonDupireOp::size() const 
	{
		return 2;
	}

	void FdmHestonDupireOp::setTime(Time t1, Time t2) const
	{
		correlationMap_->setTime(t1, t2);
		xMap_->setTime(t1, t2);
		yMap_->setTime(t1, t2);
	}

	Disposable<Array> FdmHestonDupireOp::apply_mixed(const Array& r) const 
	{
		return correlationMap_->getMap()->apply(r);
	}

	Disposable<Array> FdmHestonDupireOp::apply_direction(Size direction, const Array& r) const 
	{
		if (direction == 0)
			return xMap_->getMap()->apply(r);
		else if (direction == 1)
			return yMap_->getMap()->apply(r);
		else
			QL_FAIL("direction too large");
	}

	Disposable<Array> FdmHestonDupireOp::solve_splitting(Size direction, const Array& r, double a) const 
	{
		if (direction == 0) 
		{
			return xMap_->getMap()->solve_splitting(r, a, 1.0);
		}
		else if (direction == 1) 
		{
			return yMap_->getMap()->solve_splitting(r, a, 1.0);
		}
		else
			QL_FAIL("direction too large");
	}


}
