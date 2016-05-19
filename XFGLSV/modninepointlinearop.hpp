// modified for XFG3 LSV implementation illustration

#ifndef quantlib_mod_nine_point_linear_op_hpp
#define quantlib_mod_nine_point_linear_op_hpp

#include <ql/methods/finitedifferences/meshers/fdmmesher.hpp>
#include <ql/methods/finitedifferences/operators/ninepointlinearop.hpp>

namespace QuantLib {

    class ModNinePointLinearOp : public NinePointLinearOp {
      public:
        ModNinePointLinearOp(Size d0, Size d1,
                              const boost::shared_ptr<FdmMesher>& mesher)
        : NinePointLinearOp(d0, d1, mesher) { }

        ModNinePointLinearOp(const NinePointLinearOp& m)
        : NinePointLinearOp(m) { }

		void amx(const Array & u, const ModNinePointLinearOp& m) const;

        boost::shared_array<Real>& a00() { return a00_; }
        boost::shared_array<Real>& a10() { return a10_; }
        boost::shared_array<Real>& a20() { return a20_; }
		boost::shared_array<Real>& a01() { return a01_; }
        boost::shared_array<Real>& a11() { return a11_; }
        boost::shared_array<Real>& a21() { return a21_; }
		boost::shared_array<Real>& a02() { return a02_; }
        boost::shared_array<Real>& a12() { return a12_; }
        boost::shared_array<Real>& a22() { return a22_; }

		boost::shared_array<Size>& i00() { return i00_; }
		boost::shared_array<Size>& i10() { return i10_; }
		boost::shared_array<Size>& i20() { return i20_; }
		boost::shared_array<Size>& i01() { return i01_; }
		boost::shared_array<Size>& i21() { return i21_; }
		boost::shared_array<Size>& i02() { return i02_; }
		boost::shared_array<Size>& i12() { return i12_; }
		boost::shared_array<Size>& i22() { return i22_; }
    };
}

#endif
