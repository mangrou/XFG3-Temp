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
#include <ql/experimental/finitedifferences/fdmhestondupirefwdop.hpp>
#include <ql/experimental/finitedifferences/modtriplebandlinearop.hpp>
#include <ql/experimental/finitedifferences/modninepointlinearop.hpp>
#include <ql/experimental/finitedifferences/leveragesurface.hpp>

namespace QuantLib {


	// --------------------------- FdmHestonDupireVarianceFwdOp -------------------

	FdmHestonDupireVarianceFwdOp::FdmHestonDupireVarianceFwdOp(
             const boost::shared_ptr<FdmMesher>& mesher, 
			Real kappa, Real theta, Real sigma, Real v0)
    : kappa_(kappa),
      theta_(theta),
      sigma_(sigma),
      v0_  (v0),
	  a_(kappa*theta-0.5*sigma_*sigma_),
	  dvMap_( new ModTripleBandLinearOp(FirstDerivativeOp(1, mesher)) ),
	  dvvMap_( new ModTripleBandLinearOp(SecondDerivativeOp(1, mesher).mult(Array(mesher->layout()->size(),0.5*sigma*sigma))) ),
	  mapT_( new ModTripleBandLinearOp(1, mesher) ),
	  mesher_(mesher),
	  varianceValues_(mesher->locations(1))
	{
		for (Size i = 0; i < varianceValues_.size(); i++)
		{
			varianceValues_[i] =  1/ (exp(varianceValues_[i]) * v0_);
		}
	}

	void FdmHestonDupireVarianceFwdOp::setTime(Real t1, Real t2) const
	{
		boost::shared_array<Size> i0 = mapT_->i0(), i2 = mapT_->i2();

		const boost::shared_ptr<FdmLinearOpLayout> layout = mesher_->layout();
		const FdmLinearOpIterator endIter = layout->end();
		Real varLower, varDiag, varUpper;
		Size lowerIdx, upperIdx;
        for (FdmLinearOpIterator iter = layout->begin(); iter != endIter; ++iter) 
		{
			const Size idx = iter.index();

			lowerIdx = i0[idx];
			upperIdx = i2[idx];

			varLower = varianceValues_[lowerIdx];
			varDiag = varianceValues_[idx];
			varUpper = varianceValues_[upperIdx];

			mapT_->lower()[idx] = -dvMap_->lower()[idx] * (a_ * varLower - kappa_) + dvvMap_->lower()[idx] * varLower;
			mapT_->diag()[idx] = -dvMap_->diag()[idx] * (a_ * varDiag - kappa_) + dvvMap_->diag()[idx] * varDiag;
			mapT_->upper()[idx] = -dvMap_->upper()[idx] * (a_ * varUpper - kappa_) + dvvMap_->upper()[idx] * varUpper;
		}
	}

	// --------------------------- FdmHestonDupireEquityFwdOp -------------------

	FdmHestonDupireEquityFwdOp::FdmHestonDupireEquityFwdOp( const boost::shared_ptr<FdmMesher>& mesher, 
			const boost::shared_ptr<HestonDupireProcess>& process )
	: v0_(process->v0()),
	  rTS_(process->riskFreeRate().currentLink()),
      qTS_(process->dividendYield().currentLink()),
	  leverSurf_(process->getLeverSurf()),
	  dxxMap_( new ModTripleBandLinearOp(SecondDerivativeOp(0, mesher)) ),
	  dxMap_( new ModTripleBandLinearOp(FirstDerivativeOp(0, mesher)) ),
	  mapT_( new ModTripleBandLinearOp(0, mesher) ),
      mesher_(mesher),
	  varianceValues_(mesher->locations(1))
	{
		for (Size i = 0; i < varianceValues_.size(); i++)
		{
			varianceValues_[i] =  0.5 * v0_ * exp(varianceValues_[i]);
		}
	}

	void FdmHestonDupireEquityFwdOp::setTime(Real t1, Real t2, bool isExplicit) const
	{
		const Rate r = rTS_->forwardRate(t1, t2, Continuous).rate();
        const Rate q = qTS_->forwardRate(t1, t2, Continuous).rate();

		boost::shared_array<Size> i0 = mapT_->i0(), i2 = mapT_->i2();

		boost::shared_array<Real> leverageT1 = leverSurf_->getLeverageT1();      // leverage values at time t1
		boost::shared_array<Real> leverageT2 = leverSurf_->getLeverageT2();      // leverage values at time t1

		const boost::shared_ptr<FdmLinearOpLayout> layout = mesher_->layout();
		const FdmLinearOpIterator endIter = layout->end();
		Real varLower, varDiag, varUpper;
		Size lowerIdx, upperIdx;
		for (FdmLinearOpIterator iter = layout->begin(); iter!=endIter; ++iter) 
		{
			const Size idx = iter.index();

			lowerIdx = i0[idx];
			upperIdx = i2[idx];

			if (isExplicit == true)
			{
				varLower = leverageT1[lowerIdx];
				varDiag  = leverageT1[idx];
				varUpper = leverageT1[upperIdx];
			}
			else
			{
				varLower = leverageT2[lowerIdx];
				varDiag  = leverageT2[idx];
				varUpper = leverageT2[upperIdx];
			}

			varLower = varianceValues_[lowerIdx] * varLower * varLower;
			varDiag = varianceValues_[idx] * varDiag * varDiag;
			varUpper = varianceValues_[upperIdx] * varUpper * varUpper;

			mapT_->lower()[idx] = dxMap_->lower()[idx] * (varLower - r + q) + dxxMap_->lower()[idx] * varLower;
			mapT_->diag()[idx] = dxMap_->diag()[idx] * (varDiag - r + q) + dxxMap_->diag()[idx] * varDiag;
			mapT_->upper()[idx] = dxMap_->upper()[idx] * (varUpper - r + q) + dxxMap_->upper()[idx] * varUpper;
		}
	}


	// ----------------------    FdmHestonDupireCorrelationFwdOp    -----------------

	FdmHestonDupireCorrelationFwdOp::FdmHestonDupireCorrelationFwdOp( const boost::shared_ptr<FdmMesher>& mesher, 
			const boost::shared_ptr<HestonDupireProcess>& process )
	:
	dxyMap_( new ModNinePointLinearOp( SecondOrderMixedDerivativeOp(0, 1, mesher)
		.mult(Array(mesher->layout()->size(), process->rho()*process->sigma())) ) ),
	mapT_( new ModNinePointLinearOp(0, 1, mesher) ),
	leverSurf_(process->getLeverSurf()),
    mesher_(mesher)
	{
	}

	void FdmHestonDupireCorrelationFwdOp::setTime(Real t1, Real t2) const
	{
		boost::shared_array<Size> i00 = dxyMap_->i00(), i10 = dxyMap_->i10(), i20 = dxyMap_->i20();
		boost::shared_array<Size> i01 = dxyMap_->i01(), i21 = dxyMap_->i21();
		boost::shared_array<Size> i02 = dxyMap_->i02(), i12 = dxyMap_->i12(), i22 = dxyMap_->i22();

		boost::shared_array<Real> leverageT1 = leverSurf_->getLeverageT1();      // leverage values at time t1

		const boost::shared_ptr<FdmLinearOpLayout> layout = mesher_->layout();
		const FdmLinearOpIterator endIter = layout->end();
		for (FdmLinearOpIterator iter = layout->begin(); iter!=endIter; ++iter) 
		{
			const Size idx = iter.index();

			mapT_->a00()[idx] = dxyMap_->a00()[idx] * leverageT1[i00[idx]];
			mapT_->a10()[idx] = dxyMap_->a10()[idx] * leverageT1[i10[idx]];
			mapT_->a20()[idx] = dxyMap_->a20()[idx] * leverageT1[i20[idx]];

			mapT_->a01()[idx] = dxyMap_->a01()[idx] * leverageT1[i01[idx]];
			mapT_->a11()[idx] = dxyMap_->a11()[idx] * leverageT1[idx];
			mapT_->a21()[idx] = dxyMap_->a21()[idx] * leverageT1[i21[idx]];

			mapT_->a02()[idx] = dxyMap_->a02()[idx] * leverageT1[i02[idx]];
			mapT_->a12()[idx] = dxyMap_->a12()[idx] * leverageT1[i12[idx]];
			mapT_->a22()[idx] = dxyMap_->a22()[idx] * leverageT1[i22[idx]];
		}
	}



	// ---------------------------    FdmHestonDupireFwdOp    -----------------------

    FdmHestonDupireFwdOp::FdmHestonDupireFwdOp(
            const boost::shared_ptr<FdmMesher>& mesher,
            const boost::shared_ptr<HestonDupireProcess>& process )
    : kappa_(process->kappa()),
      theta_(process->theta()),
      sigma_(process->sigma()),
      rho_  (process->rho()),
      v0_   (process->v0()),
	  t1_(0),
	  t2_(0),
	  leverSurf_(process->getLeverSurf()),
	  mesher_(mesher),
	  xMap_( new FdmHestonDupireEquityFwdOp(mesher, process) ),
	  yMap_( new FdmHestonDupireVarianceFwdOp(mesher, process->kappa(), process->theta(), process->sigma(), process->v0()) ),
	  correlationMap_( new FdmHestonDupireCorrelationFwdOp(mesher, process) ),
      rTS_  (process->riskFreeRate().currentLink()),
      qTS_  (process->dividendYield().currentLink())
	{
		Real s0 = process->s0()->value();
		Real logS0 = log(s0);

		const boost::shared_ptr<FdmLinearOpLayout> layout = mesher_->layout();
		const FdmLinearOpIterator endIter = layout->end();
		Size equitySize = layout->dim()[0], varSize = layout->dim()[1];

		std::vector<Real> equityVec, varVec;
		equityVec.reserve(equitySize);
        varVec.reserve(varSize);

		for (FdmLinearOpIterator iter = layout->begin(); iter!=endIter; ++iter) 
		{
			if (!iter.coordinates()[1]) 
                equityVec.push_back(mesher->location(iter, 0));

            if (!iter.coordinates()[0]) 
                varVec.push_back(mesher->location(iter, 1));
		}

		for ( Size i = 0; i < equitySize; i++ )
			equityVec[i] = equityVec[i] + logS0;

		const boost::shared_ptr<Fdm1dMesher> equityMesher(new Predefined1dMesher(equityVec));
		const boost::shared_ptr<Fdm1dMesher> varianceMesher(new Predefined1dMesher(varVec));
		logSpotMesher_ = boost::shared_ptr<FdmMesher>( new FdmMesherComposite(equityMesher, varianceMesher) );
	}

	Size FdmHestonDupireFwdOp::size() const 
	{
		return 2;
	}

	Disposable<Array> FdmHestonDupireFwdOp::apply(const Array& u) const 
	{
		return xMap_->getMap()->apply(u)
                + yMap_->getMap()->apply(u)
                + correlationMap_->getMap()->apply(u);
    }

	Disposable<Array> FdmHestonDupireFwdOp::apply_mixed(const Array& u) const
	{
		return correlationMap_->getMap()->apply(u);
	}

	Disposable<Array> FdmHestonDupireFwdOp::apply_direction(Size direction, const Array& u) const
	{
		if (direction == 0)
			return xMap_->getMap()->apply(u);
		else if (direction == 1)
			return yMap_->getMap()->apply(u);
		else
			QL_FAIL("direction too large");
	}

	Disposable<Array> FdmHestonDupireFwdOp::solve_splitting(Size direction, const Array& u, Real s) const
	{
		if (direction == 0)
		{
			xMap_->setTime(t1_, t2_, false);
			Array arr = xMap_->getMap()->solve_splitting(u, s, 1.0);
			xMap_->setTime(t1_, t2_, true);
			return arr;
		}
		else if (direction == 1)
			return yMap_->getMap()->solve_splitting(u, s, 1.0);
		else
			QL_FAIL("direction too large");
	}

	void FdmHestonDupireFwdOp::setTime(Time t1, Time t2) const
	{
		t1_ = t1;
		t2_ = t2;
		leverSurf_->setTime( t1, t2, logSpotMesher_ );

		correlationMap_->setTime(t1, t2);
		xMap_->setTime(t1, t2);
		yMap_->setTime(t1, t2);
	}


}
