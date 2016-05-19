// for XFG3 LSV implementation illustration

#ifndef quantlib_leverage_surface_hpp
#define quantlib_leverage_surface_hpp

#include <ql/methods/finitedifferences/operators/fdmlinearopcomposite.hpp>
#include <ql/methods/finitedifferences/meshers/fdmmesher.hpp>

namespace QuantLib {
    

	class LeverageSurface
	{
	  protected:
		boost::shared_array<Real> leverageT1_;
		boost::shared_array<Real> leverageT2_;

	  public:

		void setTime(Real t1, Real t2, const boost::shared_ptr<FdmMesher>& mesher);

		boost::shared_array<Real> getLeverageT1() const { return leverageT1_; }
		boost::shared_array<Real> getLeverageT2() const { return leverageT2_; }

		virtual Real localVol(Real time, Real spot) const = 0;
	};

}

#endif
