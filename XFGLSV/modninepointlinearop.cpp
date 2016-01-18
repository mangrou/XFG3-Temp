// modified for XFG3 LSV implementation illustration

#include <ql/experimental/finitedifferences/modninepointlinearop.hpp>
#include <ql/methods/finitedifferences/operators/fdmlinearoplayout.hpp>

namespace QuantLib {

    void ModNinePointLinearOp::amx(const Array& a, const ModNinePointLinearOp& x) const
	{
		Size size = mesher_->layout()->size();

		for (Size i=0; i < size; ++i) 
		{
			const double s = a[i];
			a11_[i] = x.a11_[i]*s;		a00_[i] = x.a00_[i] * s;
			a01_[i] = x.a01_[i]*s;		a02_[i] = x.a02_[i] * s;
			a10_[i] = x.a10_[i]*s;		a20_[i] = x.a20_[i] * s;
			a21_[i] = x.a21_[i]*s;		a12_[i] = x.a12_[i] * s;
			a22_[i] = x.a22_[i]*s;
		}
	}


}

