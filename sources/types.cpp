// This software is MIT licensed (see LICENSE)

#define BUILD_DLL

#include <math.h>

#include <simd/core.h>

namespace simd
{
	template<>	FORCE_INLINE	uint8_t		abs<uint8_t>(const uint8_t& v) {return v;}
	template<>	FORCE_INLINE	uint16_t	abs<uint16_t>(const uint16_t& v) {return v;}
	template<>	FORCE_INLINE	uint32_t	abs<uint32_t>(const uint32_t& v) {return v;}
	template<>	FORCE_INLINE	uint64_t	abs<uint64_t>(const uint64_t& v) {return v;}
	template<>	FORCE_INLINE	int8_t		abs<int8_t>(const int8_t& v) {return v > 0 ? v : -v;}
	template<>	FORCE_INLINE	int16_t		abs<int16_t>(const int16_t& v) {return v > 0 ? v : -v;}
	template<>	FORCE_INLINE	int32_t		abs<int32_t>(const int32_t& v) {return v > 0 ? v : -v;}
	template<>	FORCE_INLINE	int64_t		abs<int64_t>(const int64_t& v) {return v > 0 ? v : -v;}
	template<>	DLL_FUNCTION(float)		abs<float>(const float& v) {return fabs(v);}
	template<>	DLL_FUNCTION(double)		abs<double>(const double& v) {return fabs(v);}

	
	template<>	DLL_FUNCTION(float)	radians(const float& degrees)	{return degrees * static_cast<float>(0.01745329251994329576923690768489);}
	template<>	DLL_FUNCTION(double)	radians(const double& degrees)	{return degrees * static_cast<double>(0.01745329251994329576923690768489);}
}

namespace simd
{
#define SIMD_ITEM_FUNCTION(_T)	\
template<>	DLL_FUNCTION(_T)						max(const _T& a, const _T& b)	{return a >= b ? a : b;} \
template<>	DLL_FUNCTION(_T)						min(const _T& a, const _T& b)	{return a <= b ? a : b;} \
template<>	DLL_FUNCTION(_T)						sqrt(const _T& v)				{return _T(::sqrt(v));} \
template<>	DLL_FUNCTION(_T)						rcp(const _T& v)				{return _T(1.0 / v);} \
template<>	DLL_FUNCTION(_T)						rsq(const _T& v)				{return _T(1.0 / ::sqrt(v));} \
template<>	DLL_FUNCTION(_T)						sin<_T>(const _T& v)			{return _T(::sin(v));} \
template<>	DLL_FUNCTION(_T)						cos<_T>(const _T& v)			{return _T(::cos(v));} \
template<>	DLL_FUNCTION(_T)						tan<_T>(const _T& v)			{return _T(::tan(v));} \
template<>	DLL_FUNCTION(lengthi_t<_T>::type)		length<_T>(const _T& v)			{return simd::abs(v);} \


	
	SIMD_ITEM_FUNCTION(int8_t)
	SIMD_ITEM_FUNCTION(uint8_t)
	SIMD_ITEM_FUNCTION(int16_t)
	SIMD_ITEM_FUNCTION(uint16_t)
	SIMD_ITEM_FUNCTION(int32_t)
	SIMD_ITEM_FUNCTION(uint32_t)
	SIMD_ITEM_FUNCTION(int64_t)
	SIMD_ITEM_FUNCTION(uint64_t)
	SIMD_ITEM_FUNCTION(float)
	SIMD_ITEM_FUNCTION(double)
}






















//namespace sml
//{
//	namespace detail
//	{
//		template<typename T>
//		DLL_FUNCTION(T) taylorCos(T const& x)
//		{
//			return static_cast<T>(1)
//			- (x * x) * (1.f / 2.f)
//			+ ((x * x) * (x * x)) * (1.f / 24.f)
//			- (((x * x) * (x * x)) * (x * x)) * (1.f / 720.f)
//			+ (((x * x) * (x * x)) * ((x * x) * (x * x))) * (1.f / 40320.f);
//		}
//
//		template<typename T>
//		DLL_FUNCTION(T) cos_52s(T x)
//		{
//			T const xx(x * x);
//			return (T(0.9999932946) + xx * (T(-0.4999124376) + xx * (T(0.0414877472) + xx * T(-0.0012712095))));
//		}
//
////		template<length_t L, typename T, qualifier Q>
////		SIMD_FUNC_QUALIFIER vec<L, T, Q> cos_52s(vec<L, T, Q> const& x)
////		{
////			return detail::functor1<L, T, T, Q>::call(cos_52s, x);
////		}
//	}//namespace detail
//
//	// wrapAngle
//	template<typename T>
//	DLL_FUNCTION(T) wrapAngle(T angle)
//	{
//		return abs<T>(mod<T>(angle, two_pi<T>()));
//	}
//
////	template<length_t L, typename T, qualifier Q>
////	SIMD_FUNC_QUALIFIER vec<L, T, Q> wrapAngle(vec<L, T, Q> const& x)
////	{
////		return detail::functor1<L, T, T, Q>::call(wrapAngle, x);
////	}
////
//	// cos
//	template<typename T>
//	DLL_FUNCTION(T) fastCos(T x)
//	{
//		T const angle(wrapAngle<T>(x));
//
//		if(angle < half_pi<T>())
//			return detail::cos_52s(angle);
//		if(angle < pi<T>())
//			return -detail::cos_52s(pi<T>() - angle);
//		if(angle < (T(3) * half_pi<T>()))
//			return -detail::cos_52s(angle - pi<T>());
//
//		return detail::cos_52s(two_pi<T>() - angle);
//	}
//
//	template<length_t L, typename T, qualifier Q>
//	SIMD_FUNC_QUALIFIER vec<L, T, Q> fastCos(vec<L, T, Q> const& x)
//	{
//		return detail::functor1<L, T, T, Q>::call(fastCos, x);
//	}
//
//	// sin
//	template<typename T>
//	DLL_FUNCTION(T) fastSin(T x)
//	{
//		return fastCos<T>(half_pi<T>() - x);
//	}
//
//	template<length_t L, typename T, qualifier Q>
//	SIMD_FUNC_QUALIFIER vec<L, T, Q> fastSin(vec<L, T, Q> const& x)
//	{
//		return detail::functor1<L, T, T, Q>::call(fastSin, x);
//	}
//
//	// tan
//	template<typename T>
//	DLL_FUNCTION(T) fastTan(T x)
//	{
//		return x + (x * x * x * T(0.3333333333)) + (x * x * x * x * x * T(0.1333333333333)) + (x * x * x * x * x * x * x * T(0.0539682539));
//	}
//
//	template<length_t L, typename T, qualifier Q>
//	SIMD_FUNC_QUALIFIER vec<L, T, Q> fastTan(vec<L, T, Q> const& x)
//	{
//		return detail::functor1<L, T, T, Q>::call(fastTan, x);
//	}
//
//	// asin
//	template<typename T>
//	DLL_FUNCTION(T) fastAsin(T x)
//	{
//		return x + (x * x * x * T(0.166666667)) + (x * x * x * x * x * T(0.075)) + (x * x * x * x * x * x * x * T(0.0446428571)) + (x * x * x * x * x * x * x * x * x * T(0.0303819444));// + (x * x * x * x * x * x * x * x * x * x * x * T(0.022372159));
//	}
//
//	template<length_t L, typename T, qualifier Q>
//	SIMD_FUNC_QUALIFIER vec<L, T, Q> fastAsin(vec<L, T, Q> const& x)
//	{
//		return detail::functor1<L, T, T, Q>::call(fastAsin, x);
//	}
//
//	// acos
//	template<typename T>
//	DLL_FUNCTION(T) fastAcos(T x)
//	{
//		return T(1.5707963267948966192313216916398) - fastAsin(x); //(PI / 2)
//	}
//
//	template<length_t L, typename T, qualifier Q>
//	SIMD_FUNC_QUALIFIER vec<L, T, Q> fastAcos(vec<L, T, Q> const& x)
//	{
//		return detail::functor1<L, T, T, Q>::call(fastAcos, x);
//	}
//
//	// atan
//	template<typename T>
//	DLL_FUNCTION(T) fastAtan(T y, T x)
//	{
//		T sgn = sign(y) * sign(x);
//		return abs(fastAtan(y / x)) * sgn;
//	}
//
//	template<length_t L, typename T, qualifier Q>
//	SIMD_FUNC_QUALIFIER vec<L, T, Q> fastAtan(vec<L, T, Q> const& y, vec<L, T, Q> const& x)
//	{
//		return detail::functor2<L, T, Q>::call(fastAtan, y, x);
//	}
//
//	template<typename T>
//	DLL_FUNCTION(T) fastAtan(T x)
//	{
//		return x - (x * x * x * T(0.333333333333)) + (x * x * x * x * x * T(0.2)) - (x * x * x * x * x * x * x * T(0.1428571429)) + (x * x * x * x * x * x * x * x * x * T(0.111111111111)) - (x * x * x * x * x * x * x * x * x * x * x * T(0.0909090909));
//	}
//
//	template<length_t L, typename T, qualifier Q>
//	SIMD_FUNC_QUALIFIER vec<L, T, Q> fastAtan(vec<L, T, Q> const& x)
//	{
//		return detail::functor1<L, T, T, Q>::call(fastAtan, x);
//	}
//}//namespace glm



namespace glm
{
////	template<length_t L, typename T, qualifier Q>
////	SIMD_FUNC_QUALIFIER SIMD_CONSTEXPR vec<L, T, Q> radians(vec<L, T, Q> const& v)
////	{
////		return detail::functor1<L, T, T, Q>::call(radians, v);
////	}
//
//	// degrees
//	template<typename T>
//	SIMD_FUNC_QUALIFIER SIMD_CONSTEXPR T degrees(T radians)
//	{
//		SIMD_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'degrees' only accept floating-point input");
//
//		return radians * static_cast<T>(57.295779513082320876798154814105);
//	}

//	template<length_t L, typename T, qualifier Q>
//	SIMD_FUNC_QUALIFIER SIMD_CONSTEXPR vec<L, T, Q> degrees(vec<L, T, Q> const& v)
//	{
//		return detail::functor1<L, T, T, Q>::call(degrees, v);
//	}
//
//	// sin
//	using ::std::sin;
//
//	template<length_t L, typename T, qualifier Q>
//	SIMD_FUNC_QUALIFIER vec<L, T, Q> sin(vec<L, T, Q> const& v)
//	{
//		return detail::functor1<L, T, T, Q>::call(sin, v);
//	}
//
//	// cos
//	using std::cos;
//
//	template<length_t L, typename T, qualifier Q>
//	SIMD_FUNC_QUALIFIER vec<L, T, Q> cos(vec<L, T, Q> const& v)
//	{
//		return detail::functor1<L, T, T, Q>::call(cos, v);
//	}
//
//	// tan
//	using std::tan;
//
//	template<length_t L, typename T, qualifier Q>
//	SIMD_FUNC_QUALIFIER vec<L, T, Q> tan(vec<L, T, Q> const& v)
//	{
//		return detail::functor1<L, T, T, Q>::call(tan, v);
//	}
//
//	// asin
//	using std::asin;
//
//	template<length_t L, typename T, qualifier Q>
//	SIMD_FUNC_QUALIFIER vec<L, T, Q> asin(vec<L, T, Q> const& v)
//	{
//		return detail::functor1<L, T, T, Q>::call(asin, v);
//	}
//
//	// acos
//	using std::acos;
//
//	template<length_t L, typename T, qualifier Q>
//	SIMD_FUNC_QUALIFIER vec<L, T, Q> acos(vec<L, T, Q> const& v)
//	{
//		return detail::functor1<L, T, T, Q>::call(acos, v);
//	}
//
//	// atan
//	template<typename T>
//	DLL_FUNCTION(T) atan(T y, T x)
//	{
//		SIMD_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'atan' only accept floating-point input");
//
//		return ::std::atan2(y, x);
//	}
//
//	template<length_t L, typename T, qualifier Q>
//	SIMD_FUNC_QUALIFIER vec<L, T, Q> atan(vec<L, T, Q> const& a, vec<L, T, Q> const& b)
//	{
//		return detail::functor2<L, T, Q>::call(::std::atan2, a, b);
//	}
//
//	using std::atan;
//
//	template<length_t L, typename T, qualifier Q>
//	SIMD_FUNC_QUALIFIER vec<L, T, Q> atan(vec<L, T, Q> const& v)
//	{
//		return detail::functor1<L, T, T, Q>::call(atan, v);
//	}
//
//	// sinh
//	using std::sinh;
//
//	template<length_t L, typename T, qualifier Q>
//	SIMD_FUNC_QUALIFIER vec<L, T, Q> sinh(vec<L, T, Q> const& v)
//	{
//		return detail::functor1<L, T, T, Q>::call(sinh, v);
//	}
//
//	// cosh
//	using std::cosh;
//
//	template<length_t L, typename T, qualifier Q>
//	SIMD_FUNC_QUALIFIER vec<L, T, Q> cosh(vec<L, T, Q> const& v)
//	{
//		return detail::functor1<L, T, T, Q>::call(cosh, v);
//	}
//
//	// tanh
//	using std::tanh;
//
//	template<length_t L, typename T, qualifier Q>
//	SIMD_FUNC_QUALIFIER vec<L, T, Q> tanh(vec<L, T, Q> const& v)
//	{
//		return detail::functor1<L, T, T, Q>::call(tanh, v);
//	}
//
//	// asinh
//#	if SIMD_HAS_CXX11_STL
//	using std::asinh;
//#	else
//	template<typename T>
//	DLL_FUNCTION(T) asinh(T x)
//	{
//		SIMD_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'asinh' only accept floating-point input");
//
//		return (x < static_cast<T>(0) ? static_cast<T>(-1) : (x > static_cast<T>(0) ? static_cast<T>(1) : static_cast<T>(0))) * log(std::abs(x) + sqrt(static_cast<T>(1) + x * x));
//	}
//#	endif
//
//	template<length_t L, typename T, qualifier Q>
//	SIMD_FUNC_QUALIFIER vec<L, T, Q> asinh(vec<L, T, Q> const& v)
//	{
//		return detail::functor1<L, T, T, Q>::call(asinh, v);
//	}
//
//	// acosh
//#	if SIMD_HAS_CXX11_STL
//	using std::acosh;
//#	else
//	template<typename T>
//	DLL_FUNCTION(T) acosh(T x)
//	{
//		SIMD_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'acosh' only accept floating-point input");
//
//		if(x < static_cast<T>(1))
//			return static_cast<T>(0);
//		return log(x + sqrt(x * x - static_cast<T>(1)));
//	}
//#	endif
//
//	template<length_t L, typename T, qualifier Q>
//	SIMD_FUNC_QUALIFIER vec<L, T, Q> acosh(vec<L, T, Q> const& v)
//	{
//		return detail::functor1<L, T, T, Q>::call(acosh, v);
//	}
//
//	// atanh
//#	if SIMD_HAS_CXX11_STL
//	using std::atanh;
//#	else
//	template<typename T>
//	DLL_FUNCTION(T) atanh(T x)
//	{
//		SIMD_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'atanh' only accept floating-point input");
//
//		if(std::abs(x) >= static_cast<T>(1))
//			return 0;
//		return static_cast<T>(0.5) * log((static_cast<T>(1) + x) / (static_cast<T>(1) - x));
//	}
//#	endif
//
//	template<length_t L, typename T, qualifier Q>
//	SIMD_FUNC_QUALIFIER vec<L, T, Q> atanh(vec<L, T, Q> const& v)
//	{
//		return detail::functor1<L, T, T, Q>::call(atanh, v);
//	}
}

//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//namespace glm{
//	namespace detail
//	{
//		/// Make a linear combination of two vectors and return the result.
//		// result = (a * ascl) + (b * bscl)
//		template<typename T, qualifier Q>
//		SIMD_FUNC_QUALIFIER vec<3, T, Q> combine(
//												vec<3, T, Q> const& a,
//												vec<3, T, Q> const& b,
//												T ascl, T bscl)
//		{
//			return (a * ascl) + (b * bscl);
//		}
//
//		template<typename T, qualifier Q>
//		SIMD_FUNC_QUALIFIER vec<3, T, Q> scale(vec<3, T, Q> const& v, T desiredLength)
//		{
//			return v * desiredLength / length(v);
//		}
//	}//namespace detail
//
//	// Matrix decompose
//	// http://www.opensource.apple.com/source/WebCore/WebCore-514/platform/graphics/transforms/TransformationMatrix.cpp
//	// Decomposes the mode matrix to translations,rotation scale components
//
//	template<typename T, qualifier Q>
//	SIMD_FUNC_QUALIFIER bool decompose(mat<4, 4, T, Q> const& ModelMatrix, vec<3, T, Q> & Scale, tquat<T, Q> & Orientation, vec<3, T, Q> & Translation, vec<3, T, Q> & Skew, vec<4, T, Q> & Perspective)
//	{
//		mat<4, 4, T, Q> LocalMatrix(ModelMatrix);
//
//		// Normalize the matrix.
//		if(epsilonEqual(LocalMatrix[3][3], static_cast<T>(0), epsilon<T>()))
//			return false;
//
//		for(length_t i = 0; i < 4; ++i)
//			for(length_t j = 0; j < 4; ++j)
//				LocalMatrix[i][j] /= LocalMatrix[3][3];
//
//		// perspectiveMatrix is used to solve for perspective, but it also provides
//		// an easy way to test for singularity of the upper 3x3 component.
//		mat<4, 4, T, Q> PerspectiveMatrix(LocalMatrix);
//
//		for(length_t i = 0; i < 3; i++)
//			PerspectiveMatrix[i][3] = static_cast<T>(0);
//		PerspectiveMatrix[3][3] = static_cast<T>(1);
//
//		/// TODO: Fixme!
//		if(epsilonEqual(determinant(PerspectiveMatrix), static_cast<T>(0), epsilon<T>()))
//			return false;
//
//		// First, isolate perspective.  This is the messiest.
//		if(
//		   epsilonNotEqual(LocalMatrix[0][3], static_cast<T>(0), epsilon<T>()) ||
//		   epsilonNotEqual(LocalMatrix[1][3], static_cast<T>(0), epsilon<T>()) ||
//		   epsilonNotEqual(LocalMatrix[2][3], static_cast<T>(0), epsilon<T>()))
//		{
//			// rightHandSide is the right hand side of the equation.
//			vec<4, T, Q> RightHandSide;
//			RightHandSide[0] = LocalMatrix[0][3];
//			RightHandSide[1] = LocalMatrix[1][3];
//			RightHandSide[2] = LocalMatrix[2][3];
//			RightHandSide[3] = LocalMatrix[3][3];
//
//			// Solve the equation by inverting PerspectiveMatrix and multiplying
//			// rightHandSide by the inverse.  (This is the easiest way, not
//			// necessarily the best.)
//			mat<4, 4, T, Q> InversePerspectiveMatrix = glm::inverse(PerspectiveMatrix);//   inverse(PerspectiveMatrix, inversePerspectiveMatrix);
//			mat<4, 4, T, Q> TransposedInversePerspectiveMatrix = glm::transpose(InversePerspectiveMatrix);//   transposeMatrix4(inversePerspectiveMatrix, transposedInversePerspectiveMatrix);
//
//			Perspective = TransposedInversePerspectiveMatrix * RightHandSide;
//			//  v4MulPointByMatrix(rightHandSide, transposedInversePerspectiveMatrix, perspectivePoint);
//
//			// Clear the perspective partition
//			LocalMatrix[0][3] = LocalMatrix[1][3] = LocalMatrix[2][3] = static_cast<T>(0);
//			LocalMatrix[3][3] = static_cast<T>(1);
//		}
//		else
//		{
//			// No perspective.
//			Perspective = vec<4, T, Q>(0, 0, 0, 1);
//		}
//
//		// Next take care of translation (easy).
//		Translation = vec<3, T, Q>(LocalMatrix[3]);
//		LocalMatrix[3] = vec<4, T, Q>(0, 0, 0, LocalMatrix[3].w);
//
//		vec<3, T, Q> Row[3], Pdum3;
//
//		// Now get scale and shear.
//		for(length_t i = 0; i < 3; ++i)
//			for(length_t j = 0; j < 3; ++j)
//				Row[i][j] = LocalMatrix[i][j];
//
//		// Compute X scale factor and normalize first row.
//		Scale.x = length(Row[0]);// v3Length(Row[0]);
//
//		Row[0] = detail::scale(Row[0], static_cast<T>(1));
//
//		// Compute XY shear factor and make 2nd row orthogonal to 1st.
//		Skew.z = dot(Row[0], Row[1]);
//		Row[1] = detail::combine(Row[1], Row[0], static_cast<T>(1), -Skew.z);
//
//		// Now, compute Y scale and normalize 2nd row.
//		Scale.y = length(Row[1]);
//		Row[1] = detail::scale(Row[1], static_cast<T>(1));
//		Skew.z /= Scale.y;
//
//		// Compute XZ and YZ shears, orthogonalize 3rd row.
//		Skew.y = glm::dot(Row[0], Row[2]);
//		Row[2] = detail::combine(Row[2], Row[0], static_cast<T>(1), -Skew.y);
//		Skew.x = glm::dot(Row[1], Row[2]);
//		Row[2] = detail::combine(Row[2], Row[1], static_cast<T>(1), -Skew.x);
//
//		// Next, get Z scale and normalize 3rd row.
//		Scale.z = length(Row[2]);
//		Row[2] = detail::scale(Row[2], static_cast<T>(1));
//		Skew.y /= Scale.z;
//		Skew.x /= Scale.z;
//
//		// At this point, the matrix (in rows[]) is orthonormal.
//		// Check for a coordinate system flip.  If the determinant
//		// is -1, then negate the matrix and the scaling factors.
//		Pdum3 = cross(Row[1], Row[2]); // v3Cross(row[1], row[2], Pdum3);
//		if(dot(Row[0], Pdum3) < 0)
//		{
//			for(length_t i = 0; i < 3; i++)
//			{
//				Scale[i] *= static_cast<T>(-1);
//				Row[i] *= static_cast<T>(-1);
//			}
//		}
//
//		// Now, get the rotations out, as described in the gem.
//
//		// FIXME - Add the ability to return either quaternions (which are
//		// easier to recompose with) or Euler angles (rx, ry, rz), which
//		// are easier for authors to deal with. The latter will only be useful
//		// when we fix https://bugs.webkit.org/show_bug.cgi?id=23799, so I
//		// will leave the Euler angle code here for now.
//
//		// ret.rotateY = asin(-Row[0][2]);
//		// if (cos(ret.rotateY) != 0) {
//		//     ret.rotateX = atan2(Row[1][2], Row[2][2]);
//		//     ret.rotateZ = atan2(Row[0][1], Row[0][0]);
//		// } else {
//		//     ret.rotateX = atan2(-Row[2][0], Row[1][1]);
//		//     ret.rotateZ = 0;
//		// }
//
//		int i, j, k = 0;
//		float root, trace = Row[0].x + Row[1].y + Row[2].z;
//		if(trace > static_cast<T>(0))
//		{
//			root = sqrt(trace + static_cast<T>(1.0));
//			Orientation.w = static_cast<T>(0.5) * root;
//			root = static_cast<T>(0.5) / root;
//			Orientation.x = root * (Row[1].z - Row[2].y);
//			Orientation.y = root * (Row[2].x - Row[0].z);
//			Orientation.z = root * (Row[0].y - Row[1].x);
//		} // End if > 0
//		else
//		{
//			static int Next[3] = {1, 2, 0};
//			i = 0;
//			if(Row[1].y > Row[0].x) i = 1;
//			if(Row[2].z > Row[i][i]) i = 2;
//			j = Next[i];
//			k = Next[j];
//
//			root = sqrt(Row[i][i] - Row[j][j] - Row[k][k] + static_cast<T>(1.0));
//
//			Orientation[i] = static_cast<T>(0.5) * root;
//			root = static_cast<T>(0.5) / root;
//			Orientation[j] = root * (Row[i][j] + Row[j][i]);
//			Orientation[k] = root * (Row[i][k] + Row[k][i]);
//			Orientation.w = root * (Row[j][k] - Row[k][j]);
//		} // End if <= 0
//
//		return true;
//	}
//}//namespace glm



