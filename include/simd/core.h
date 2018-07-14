// This software is MIT licensed (see LICENSE)

#pragma once

#include <pre.h>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Preprocessor directives
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define SIMD_VERSION_MAJOR			0
#define SIMD_VERSION_MINOR			9
#define SIMD_VERSION_PATCH			9
#define SIMD_VERSION_REVISION		0
#define SIMD_VERSION					990

#define SIMD_SETUP_INCLUDED SIMD_VERSION

///////////////////////////////////////////////////////////////////////////////////
// Qualifiers

#if COMPILER_PLATFORM & COMPILER_PLATFORM_CUDA
#	define SIMD_CUDA_FUNC_DEF __device__ __host__
#	define SIMD_CUDA_FUNC_DECL __device__ __host__
#else
#	define SIMD_CUDA_FUNC_DEF
#	define SIMD_CUDA_FUNC_DECL
#endif

#if COMPILER_PLATFORM & COMPILER_PLATFORM_GCC
#	define SIMD_VAR_USED __attribute__ ((unused))
#else
#	define SIMD_VAR_USED
#endif


#define SIMD_STATIC_INLINE	static	FORCE_INLINE

#define SIMD_FUNC_DECL SIMD_CUDA_FUNC_DECL
#define SIMD_FUNC_QUALIFIER SIMD_CUDA_FUNC_DEF FORCE_INLINE

///////////////////////////////////////////////////////////////////////////////////
// Swizzle operators

// User defines: PRE_FORCE_SWIZZLE

#define SIMD_SWIZZLE_ENABLED 1
#define SIMD_SWIZZLE_DISABLE 0

#if defined(PRE_FORCE_SWIZZLE)
#	define SIMD_SWIZZLE SIMD_SWIZZLE_ENABLED
#else
#	define SIMD_SWIZZLE SIMD_SWIZZLE_DISABLE
#endif

#if PRE_MESSAGES == PRE_MESSAGES_ENABLED && !defined(PRE_MESSAGE_SWIZZLE_DISPLAYED)
#	define PRE_MESSAGE_SWIZZLE_DISPLAYED
#	if SIMD_SWIZZLE == SIMD_SWIZZLE_ENABLED
#		pragma message("PRE: Swizzling operators enabled")
#	else
#		pragma message("PRE: Swizzling operators disabled, #define SIMD_SWIZZLE to enable swizzle operators")
#	endif
#endif//PRE_MESSAGES

///////////////////////////////////////////////////////////////////////////////////
// Clip control

#define SIMD_DEPTH_ZERO_TO_ONE				0x00000001
#define SIMD_DEPTH_NEGATIVE_ONE_TO_ONE		0x00000002

#ifdef PRE_FORCE_DEPTH_ZERO_TO_ONE
#	define SIMD_DEPTH_CLIP_SPACE SIMD_DEPTH_ZERO_TO_ONE
#else
#	define SIMD_DEPTH_CLIP_SPACE SIMD_DEPTH_NEGATIVE_ONE_TO_ONE
#endif

#if PRE_MESSAGES == PRE_MESSAGES_ENABLED && !defined(PRE_MESSAGE_DEPTH_DISPLAYED)
#	define PRE_MESSAGE_DEPTH_DISPLAYED
#	if SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_ZERO_TO_ONE
#		pragma message("PRE: Depth clip space: Zero to one")
#	else
#		pragma message("PRE: Depth clip space: negative one to one")
#	endif
#endif//PRE_MESSAGES

///////////////////////////////////////////////////////////////////////////////////
// Coordinate system, define PRE_FORCE_LEFT_HANDED before including GLM
// to use left handed coordinate system by default.

#define SIMD_LEFT_HANDED				0x00000001	// For DirectX, Metal, Vulkan
#define SIMD_RIGHT_HANDED			0x00000002	// For OpenGL, default in GLM

#ifdef PRE_FORCE_LEFT_HANDED
#	define SIMD_COORDINATE_SYSTEM SIMD_LEFT_HANDED
#else
#	define SIMD_COORDINATE_SYSTEM SIMD_RIGHT_HANDED
#endif

#if PRE_MESSAGES == PRE_MESSAGES_ENABLED && !defined(PRE_MESSAGE_HANDED_DISPLAYED)
#	define PRE_MESSAGE_HANDED_DISPLAYED
#	if SIMD_COORDINATE_SYSTEM == SIMD_LEFT_HANDED
#		pragma message("PRE: Coordinate system: left handed")
#	else
#		pragma message("PRE: Coordinate system: right handed")
#	endif
#endif//PRE_MESSAGES

// With MinGW-W64, including intrinsic headers before intrin.h will produce some errors. The problem is
// that windows.h (and maybe other headers) will silently include intrin.h, which of course causes problems.
// To fix, we just explicitly include intrin.h here.
#if defined(__MINGW64__) && (CPU_ARCH != CPU_ARCH_PURE)
#	include <intrin.h>
#endif

#if CPU_ARCH & CPU_ARCH_AVX2_BIT
#	include <immintrin.h>
#elif CPU_ARCH & CPU_ARCH_AVX_BIT
#	include <immintrin.h>
#elif CPU_ARCH & CPU_ARCH_SSE42_BIT
#	if COMPILER_PLATFORM & COMPILER_PLATFORM_CLANG
#		include <popcntintrin.h>
#	endif
#	include <nmmintrin.h>
#elif CPU_ARCH & CPU_ARCH_SSE41_BIT
#	include <smmintrin.h>
#elif CPU_ARCH & CPU_ARCH_SSSE3_BIT
#	include <tmmintrin.h>
#elif CPU_ARCH & CPU_ARCH_SSE3_BIT
#	include <pmmintrin.h>
#elif CPU_ARCH & CPU_ARCH_SSE2_BIT
#	include <mmintrin.h>
#	include <xmmintrin.h>
#	include <emmintrin.h>
#elif CPU_ARCH & CPU_ARCH_NEON_BIT
#	include <arm_neon.h>
#endif//CPU_ARCH

#include <stdint.h>
#include <assert.h>
#include <cstddef>

///////////////////////////////////////////////////////////////////////////////////
// Length type: all length functions returns a length_t type.
// When PRE_FORCE_SIZE_T_LENGTH is defined, length_t is a typedef of size_t otherwise
// length_t is a typedef of int like GLSL defines it.

// User define: PRE_FORCE_SIZE_T_LENGTH

namespace glm
{
	using std::size_t;
#	if defined(PRE_FORCE_SIZE_T_LENGTH)
	typedef size_t length_t;
#	else
	typedef int length_t;
#	endif
}//namespace glm

#if PRE_MESSAGES == PRE_MESSAGES_ENABLED && !defined(PRE_MESSAGE_FORCE_SIZE_T_LENGTH)
#	define PRE_MESSAGE_FORCE_SIZE_T_LENGTH
#	if defined PRE_FORCE_SIZE_T_LENGTH
#		pragma message("GLM: .length() returns glm::length_t, a typedef of std::size_t")
#	else
#		pragma message("GLM: .length() returns glm::length_t, a typedef of int following the GLSL specification")
#	endif
#endif//PRE_MESSAGES

///////////////////////////////////////////////////////////////////////////////////
// countof

#if SIMD_HAS_CONSTEXPR_PARTIAL
namespace sml
{
	template<typename T, std::size_t N>
	constexpr std::size_t countof(T const (&)[N])
	{
		return N;
	}
}//namespace glm
#	define SIMD_COUNTOF(arr) glm::countof(arr)
#elif defined(_MSC_VER)
#	define SIMD_COUNTOF(arr) _countof(arr)
#else
#	define SIMD_COUNTOF(arr) sizeof(arr) / sizeof(arr[0])
#endif


///////////////////////////////////////////////////////////////////////////////////
// Function prototypes

#define SIMD_DLLFNC_D(_Ret)		template <int D>	DLL_FUNCTION(_Ret)
#define SIMD_DLLFNC_DT(_Ret)	template <int D, typename T>	DLL_FUNCTION(_Ret)
#define SIMD_DLLFNC_T(_Ret)		template <typename T>	DLL_FUNCTION(_Ret)
#define SIMD_DLLFNC_T2(_Ret)	template <typename T1, typename T2>	DLL_FUNCTION(_Ret)
#define SIMD_DLLFNC_T3(_Ret)	template <typename T1, typename T2, typename T3>	DLL_FUNCTION(_Ret)
#define SIMD_INLFNC_T(_Ret)		template <typename T>	FORCE_INLINE	_Ret
#define SIMD_INLFNC_T2(_Ret)	template <typename T1, typename T2>	FORCE_INLINE	_Ret
#define SIMD_INLFNC_T3(_Ret)	template <typename T1, typename T2, typename T3>	FORCE_INLINE	_Ret
#define DEFINE_VECTOR_OP_EQ(_dim) \
	FORCE_INLINE	auto	operator [] (int Index) -> T& 				{assert(0 <= Index && Index < component_count()); return ((T*)(&m))[Index];}; \
	FORCE_INLINE	auto	operator [] (int Index) const -> const T& 	{assert(0 <= Index && Index < component_count()); return ((const T*)(&m))[Index];}; \
	template <typename T2> \
	FORCE_INLINE	void	operator = (const pack_t<_dim, T2>& v) {for (int i = 0; i < _dim; i++)channel(i) = T(v[i]);} \
	DLL_FUNCTION(void)	operator += (const pack_t<_dim, T>&); \
	DLL_FUNCTION(void)	operator += (const T&); \
	DLL_FUNCTION(void)	operator -= (const pack_t<_dim, T>&); \
	DLL_FUNCTION(void)	operator -= (const T&); \
	DLL_FUNCTION(void)	operator *= (const pack_t<_dim, T>&); \
	DLL_FUNCTION(void)	operator *= (const T&); \
	DLL_FUNCTION(void)	operator /= (const pack_t<_dim, T>&); \
	DLL_FUNCTION(void)	operator /= (const T&);
#define SIMD_DECLARE_EXTERN_PACK(_T) \
	extern template struct pack_t<1, _T>; \
	extern template struct pack_t<2, _T>; \
	extern template struct pack_t<3, _T>; \
	extern template struct pack_t<4, _T>; \
	extern template struct pack_t<8, _T>; \
	extern template struct pack_t<16, _T>;
#define SIMD_DEFINE_PACK_T_P(_T, _P) \
	typedef pack_t<16, _T >		pack16_##_P; \
	typedef pack_t<8, _T >		pack8_##_P; \
	typedef pack_t<4, _T >		pack4_##_P; \
	typedef pack_t<3, _T >		pack3_##_P; \
	typedef pack_t<2, _T >		pack2_##_P; \
	typedef pack_t<1, _T >		pack1_##_P; \
	typedef pack_t<16, _T >		_P##pack16; \
	typedef pack_t<8, _T >		_P##pack8; \
	typedef pack_t<4, _T >		_P##pack4; \
	typedef pack_t<3, _T >		_P##pack3; \
	typedef pack_t<2, _T >		_P##pack2; \
	typedef pack_t<1, _T >		_P##pack1; \
	typedef pack_t<16, _T >		vec16d_##_P; \
	typedef pack_t<8, _T >		vec8d_##_P; \
	typedef pack_t<4, _T >		vec4d_##_P; \
	typedef pack_t<3, _T >		vec3d_##_P; \
	typedef pack_t<2, _T >		vec2d_##_P; \
	typedef pack_t<1, _T >		vec1d_##_P;
#define SIMD_MAT_XYZW_SO_ON template< \
	typename X1, typename Y1, typename Z1, typename W1, \
	typename X2, typename Y2, typename Z2, typename W2, \
	typename X3, typename Y3, typename Z3, typename W3, \
	typename X4, typename Y4, typename Z4, typename W4>
#define SIMD_MAT_V1V2V3V4 template<typename V1, typename V2, typename V3, typename V4>
#define SIMD_TMAT_EXPL_CRU template<typename U>
#define SIMD_DECLARE_EXTERN_MAT(_T) extern template struct mat_t<4, 4, _T>;
#define SIMD_INL_U 	template<typename U>	FORCE_INLINE
#define SIMD_INL_UV	template<typename U, typename V>	FORCE_INLINE






















#if CPU_ARCH & CPU_ARCH_SSE2_BIT
typedef __m64		SIMD_vec2;
typedef __m128		SIMD_vec4;
typedef __m128i		SIMD_ivec4;
typedef __m128i		SIMD_uvec4;
#endif

#if CPU_ARCH & CPU_ARCH_NEON_BIT
typedef float32x2_t	SIMD_vec2;
typedef float32x4_t	SIMD_vec4;
typedef int32x4_t	SIMD_ivec4;
typedef uint32x4_t	SIMD_uvec4;
#endif

#if CPU_ARCH & CPU_ARCH_AVX_BIT
typedef __m256d		SIMD_dvec4;
#endif

#if CPU_ARCH & CPU_ARCH_AVX2_BIT
typedef __m256i		SIMD_i64vec4;
typedef __m256i		SIMD_u64vec4;
#endif

namespace simd
{
	template <typename T>	struct length_t				{typedef	double	type;};
	template <>				struct length_t<float>		{typedef	float	type;};
	template <typename T>	struct lengthi_t			{typedef	T	type;};
	template <>				struct lengthi_t<float>		{typedef	float	type;};
	template <>				struct lengthi_t<double>	{typedef	double	type;};
	
	template <typename T1, typename T2>	struct high_precission_t					{typedef	T1	type;};
	template <typename T>				struct high_precission_t<T, T>				{typedef	T	type;};
	template <>							struct high_precission_t<double, double>	{typedef	double	type;};
	template <typename T>				struct high_precission_t<T, double>			{typedef	double	type;};
	template <typename T>				struct high_precission_t<double, T>			{typedef	double	type;};
	
	template<typename T, int D>	struct reg	{typedef T type_t[D];};
	
#if CPU_ARCH & CPU_ARCH_NEON_BIT
	template<>	struct reg<float, 2>		{typedef SIMD_vec2 type_t;};
	template<>	struct reg<float, 4>		{typedef SIMD_vec4 type_t;};
	template<>	struct reg<int, 4>			{typedef SIMD_ivec4 type_t;};
	template<>	struct reg<double, 4>		{typedef float64x1x4_t type_t;};
#elif CPU_ARCH & CPU_ARCH_SSE2_BIT
	template<>	struct reg<float, 2>		{typedef SIMD_vec2 type_t;};
	template<>	struct reg<float, 4>		{typedef SIMD_vec4 type_t;};
	template<>	struct reg<int, 4>			{typedef __m128i type_t;};
	template<>	struct reg<double, 4>		{typedef __m128d type_t;};
#elif defined(_MSC_VER)
	template<>	struct reg<float, 2>		{typedef float	type_t[2];};
	template<>	struct reg<float, 4>		{typedef float	type_t	__attribute__	((__vector_size__ (16), __may_alias__));};
	template<>	struct reg<int, 4>			{typedef int 	type_t	__attribute__	((__vector_size__ (16), __may_alias__));};
	template<>	struct reg<double, 4>		{typedef double type_t	__attribute__	((__vector_size__ (16), __may_alias__));};
#else
	template<>	struct reg<float, 2>		{typedef float	type_t	__attribute__	((__vector_size__ (8), __may_alias__));};
	template<>	struct reg<float, 4>		{typedef float	type_t	__attribute__ 	((__vector_size__ (16), __may_alias__));};
	template<>	struct reg<int, 4>			{typedef int	type_t	__attribute__ 	((__vector_size__ (16), __may_alias__));};
	template<>	struct reg<double, 4>		{typedef double type_t	__attribute__ 	((__vector_size__ (16), __may_alias__));};
#endif
	
	
	
	
#if defined(_MSC_VER)
#else
	template<>	struct reg<int8_t,2>{typedef int8_t type_t __attribute__		((__vector_size__ (2), __may_alias__));};
	//	template<>	struct reg<int8_t,3>{typedef int8_t type_t __attribute__		((__vector_size__ (4), __may_alias__));};
	template<>	struct reg<int8_t,4>{typedef int8_t type_t __attribute__		((__vector_size__ (4), __may_alias__));};
	template<>	struct reg<int8_t,8>{typedef int8_t type_t __attribute__		((__vector_size__ (8), __may_alias__));};
	template<>	struct reg<int8_t,16>{typedef int8_t type_t __attribute__		((__vector_size__ (16), __may_alias__));};
	template<>	struct reg<uint8_t,2>{typedef uint8_t type_t __attribute__		((__vector_size__ (2), __may_alias__));};
	//	template<>	struct reg<uint8_t,3>{typedef uint8_t type_t __attribute__		((__vector_size__ (4), __may_alias__));};
	template<>	struct reg<uint8_t,4>{typedef uint8_t type_t __attribute__		((__vector_size__ (4), __may_alias__));};
	template<>	struct reg<uint8_t,8>{typedef uint8_t type_t __attribute__		((__vector_size__ (8), __may_alias__));};
	template<>	struct reg<uint8_t,16>{typedef uint8_t type_t __attribute__		((__vector_size__ (16), __may_alias__));};
	
	template<>	struct reg<int16_t,2>{typedef int16_t type_t __attribute__		((__vector_size__ (4), __may_alias__));};
	//	template<>	struct reg<int16_t,3>{typedef int16_t type_t __attribute__		((__vector_size__ (8), __may_alias__));};
	template<>	struct reg<int16_t,4>{typedef int16_t type_t __attribute__		((__vector_size__ (8), __may_alias__));};
	template<>	struct reg<int16_t,8>{typedef int16_t type_t __attribute__		((__vector_size__ (16), __may_alias__));};
	template<>	struct reg<int16_t,16>{typedef int16_t type_t __attribute__		((__vector_size__ (16), __may_alias__));};
	template<>	struct reg<uint16_t,2>{typedef uint16_t type_t __attribute__	((__vector_size__ (4), __may_alias__));};
	//	template<>	struct reg<uint16_t,3>{typedef uint16_t type_t __attribute__	((__vector_size__ (8), __may_alias__));};
	template<>	struct reg<uint16_t,4>{typedef uint16_t type_t __attribute__	((__vector_size__ (8), __may_alias__));};
	template<>	struct reg<uint16_t,8>{typedef uint16_t type_t __attribute__	((__vector_size__ (16), __may_alias__));};
	template<>	struct reg<uint16_t,16>{typedef uint16_t type_t __attribute__	((__vector_size__ (16), __may_alias__));};
	
	template<>	struct reg<int32_t,2>{typedef int32_t type_t __attribute__		((__vector_size__ (8), __may_alias__));};
	//	template<>	struct reg<int32_t,3>{typedef int32_t type_t __attribute__		((__vector_size__ (16), __may_alias__));};
	//	template<>	struct reg<int32_t,4>{typedef int32_t type_t __attribute__		((__vector_size__ (16), __may_alias__));};
	template<>	struct reg<int32_t,8>{typedef int32_t type_t __attribute__		((__vector_size__ (16), __may_alias__));};
	template<>	struct reg<int32_t,16>{typedef int32_t type_t __attribute__		((__vector_size__ (16), __may_alias__));};
	template<>	struct reg<uint32_t,2>{typedef uint32_t type_t __attribute__	((__vector_size__ (8), __may_alias__));};
	//	template<>	struct reg<uint32_t,3>{typedef uint32_t type_t __attribute__	((__vector_size__ (16), __may_alias__));};
	template<>	struct reg<uint32_t,4>{typedef uint32_t type_t __attribute__	((__vector_size__ (16), __may_alias__));};
	template<>	struct reg<uint32_t,8>{typedef uint32_t type_t __attribute__	((__vector_size__ (16), __may_alias__));};
	template<>	struct reg<uint32_t,16>{typedef uint32_t type_t __attribute__	((__vector_size__ (16), __may_alias__));};
	
	template<>	struct reg<int64_t,2>{typedef int64_t type_t __attribute__		((__vector_size__ (16), __may_alias__));};
	//	template<>	struct reg<int64_t,3>{typedef int64_t type_t __attribute__		((__vector_size__ (16), __may_alias__));};
	template<>	struct reg<int64_t,4>{typedef int64_t type_t __attribute__		((__vector_size__ (16), __may_alias__));};
	template<>	struct reg<int64_t,8>{typedef int64_t type_t __attribute__		((__vector_size__ (16), __may_alias__));};
	template<>	struct reg<int64_t,16>{typedef int64_t type_t __attribute__		((__vector_size__ (16), __may_alias__));};
	template<>	struct reg<uint64_t,2>{typedef uint64_t type_t __attribute__	((__vector_size__ (16), __may_alias__));};
	//	template<>	struct reg<uint64_t,3>{typedef uint64_t type_t __attribute__	((__vector_size__ (16), __may_alias__));};
	template<>	struct reg<uint64_t,4>{typedef uint64_t type_t __attribute__	((__vector_size__ (16), __may_alias__));};
	template<>	struct reg<uint64_t,8>{typedef uint64_t type_t __attribute__	((__vector_size__ (16), __may_alias__));};
	template<>	struct reg<uint64_t,16>{typedef uint64_t type_t __attribute__	((__vector_size__ (16), __may_alias__));};
	
	//	template<>	struct reg<float,2>{typedef float type_t __attribute__			((__vector_size__ (8), __may_alias__));};
	//	template<>	struct reg<float,3>{typedef float type_t __attribute__			((__vector_size__ (16), __may_alias__));};
	//	template<>	struct reg<float,4>{typedef float type_t __attribute__			((__vector_size__ (16), __may_alias__));};
	template<>	struct reg<float,8>{typedef float type_t __attribute__			((__vector_size__ (16), __may_alias__));};
	template<>	struct reg<float,16>{typedef float type_t __attribute__			((__vector_size__ (16), __may_alias__));};
	
	template<>	struct reg<double,2>{typedef double type_t __attribute__		((__vector_size__ (16), __may_alias__));};
	//	template<>	struct reg<double,3>{typedef double type_t __attribute__		((__vector_size__ (16), __may_alias__));};
	//	template<>	struct reg<double,4>{typedef double type_t __attribute__		((__vector_size__ (16), __may_alias__));};
	template<>	struct reg<double,8>{typedef double type_t __attribute__		((__vector_size__ (16), __may_alias__));};
	template<>	struct reg<double,16>{typedef double type_t __attribute__		((__vector_size__ (16), __may_alias__));};
#endif
	
}




























////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Core classes
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace simd
{
	template <typename T, bool S>
	struct normalized_t
	{
		T	m;
	};
	
	typedef normalized_t<int8_t, true>		nint8;
	typedef normalized_t<int16_t, true> 	nint16;
	typedef normalized_t<int32_t, true> 	nint32;
	typedef normalized_t<int64_t, true> 	nint64;
	typedef normalized_t<uint8_t, false>	unint8;
	typedef normalized_t<uint16_t, false>	unint16;
	typedef normalized_t<uint32_t, false>	unint32;
	typedef normalized_t<uint64_t, false>	unint64;
}

namespace simd
{
	SIMD_DLLFNC_T2(T1)		add(const T1&, const T2&);
	SIMD_DLLFNC_T2(T1)		sub(const T1&, const T2&);
	SIMD_DLLFNC_T2(T1)		mul(const T1&, const T2&);
	SIMD_DLLFNC_T2(T1)		div(const T1&, const T2&);
	SIMD_DLLFNC_T(T)		min(const T&, const T&);
	SIMD_DLLFNC_T(T)		max(const T&, const T&);
	SIMD_DLLFNC_T(T)		sqrt(const T& v);
	SIMD_DLLFNC_T(T)		rsq(const T& v);	// rcp (reciprocal square root) = 1 / sqrt(v)
	SIMD_DLLFNC_T(T)		rcp(const T& v);	// rcp (reciprocal) = 1 / v

	SIMD_DLLFNC_T2(T1)		mulComponents(const T1&, const T2&);
	SIMD_DLLFNC_T2(T1)		mod(const T1&, const T2&);
	SIMD_DLLFNC_T2(T1)		modf(const T1&, const T2&);
	SIMD_DLLFNC_T(T)		abs(const T&);
	SIMD_DLLFNC_T(T)		sign(const T&);
	SIMD_DLLFNC_T(T)		sign(const T&);
	SIMD_INLFNC_T3(T1)		mad(const T1& a, const T2& b, const T3& c)	{return add(mul(a, b), c);}
	SIMD_INLFNC_T3(T1)		fma(const T1& a, const T2& b, const T3& c)	{return mad(a, b, c);}

	SIMD_DLLFNC_T(T)		floor(const T&);
	SIMD_DLLFNC_T(T)		trunc(const T&);
	SIMD_DLLFNC_T(T)		round(const T&);
	SIMD_DLLFNC_T(T)		roundEven(const T&);
	SIMD_DLLFNC_T(T)		ceil(const T&);
	SIMD_DLLFNC_T(T)		fract(const T&);
	SIMD_DLLFNC_T(T)		ceil(const T&);
	SIMD_DLLFNC_T(T)		ceil(const T&);

	SIMD_DLLFNC_T(bool)		isnan(const T&);
	SIMD_DLLFNC_T(bool)		isinf(const T&);
	SIMD_DLLFNC_T(bool)		isfinite(const T&);

	SIMD_DLLFNC_T2(T1)		clamp(const T1&, const T2&);
	SIMD_DLLFNC_T(T)		saturate(const T&);
	SIMD_DLLFNC_T2(T1)		step(const T1& x, const T1& y, const T2& a);
	SIMD_DLLFNC_T2(T1)		smoothstep(const T1& x, const T1& y, const T2& a);
	SIMD_DLLFNC_T2(T1)		mix(const T1& x, const T1& y, const T2& a);
	SIMD_INLFNC_T2(T1)		lerp(const T1& x, const T1& y, const T2& a){return mix(x, y, a);}

	SIMD_DLLFNC_T2(T1)		frexp(const T1&, const T2&);
	SIMD_DLLFNC_T2(T1)		ldexp(const T1&, const T2&);

	DLL_FUNCTION(int)		floatBitsToInt(const float&);
	DLL_FUNCTION(auto)		floatBitsToUint(const float&) -> unsigned int;
	DLL_FUNCTION(auto)		intBitsToFloat(const int& v) -> float;
	DLL_FUNCTION(auto)		uintBitsToFloat(const unsigned int& v) -> float;

	SIMD_DLLFNC_T(auto)		length(const T&) -> typename lengthi_t<T>::type;
	SIMD_INLFNC_T(auto)		distance(const T& a, const T& b) -> typename lengthi_t<T>::type {return length(b - a);}
	SIMD_DLLFNC_T(auto)		dot(const T&, const T&) -> typename length_t<T>::type;
	SIMD_DLLFNC_T(auto)		dot3(const T&, const T&) -> typename length_t<T>::type;
	SIMD_DLLFNC_T(T)		normalize(const T&);
	SIMD_DLLFNC_T(T)		faceforward(const T& N, const T& I, const T& Nref);
	SIMD_DLLFNC_T(T)		reflect(const T& I, const T& N);
	SIMD_DLLFNC_T(T)		refract(const T& I, const T& N, T eta);

	SIMD_DLLFNC_T(T)		radians(const T& degrees);
	SIMD_DLLFNC_T(T)		degrees(const T& degrees);
	SIMD_DLLFNC_T(T)		cos(const T& angle);
	SIMD_DLLFNC_T(T)		sin(const T& angle);
	SIMD_DLLFNC_T(T)		tan(const T& angle);
	SIMD_DLLFNC_T(T)		cosh(const T& angle);
	SIMD_DLLFNC_T(T)		sinh(const T& angle);
	SIMD_DLLFNC_T(T)		tanh(const T& angle);
	SIMD_DLLFNC_T(T)		acos(const T& angle);
	SIMD_DLLFNC_T(T)		asin(const T& angle);
	SIMD_DLLFNC_T(T)		atan(const T& angle);
	SIMD_DLLFNC_T(T)		atan2(const T& x, const T& y);
	SIMD_DLLFNC_T(T)		acosh(const T& angle);
	SIMD_DLLFNC_T(T)		asinh(const T& angle);
	SIMD_DLLFNC_T(T)		atanh(const T& angle);
	SIMD_DLLFNC_T(T)		wrapAngle(T angle);
	SIMD_DLLFNC_T(T)		fastSin(T angle);
	SIMD_DLLFNC_T(T)		fastCos(T angle);
	SIMD_DLLFNC_T(T)		fastTan(T angle);
	SIMD_DLLFNC_T(T)		fastAsin(T angle);
	SIMD_DLLFNC_T(T)		fastAcos(T angle);
	SIMD_DLLFNC_T(T)		fastAtan(T y, T x);
	SIMD_DLLFNC_T(T)		fastAtan(T angle);
}

namespace simd
{
	template <typename R, typename T>
	FORCE_INLINE	void	create_register(R& o, T x, T y, T z, T w)
	{
		((T*)&o)[0] = x;
		((T*)&o)[1] = y;
		((T*)&o)[2] = z;
		((T*)&o)[3] = w;
	}

#if CPU_ARCH & (CPU_ARCH_SSE2_BIT | CPU_ARCH_NEON_BIT)
	FORCE_INLINE void SIMD_vec4_create(SIMD_vec4& o, float x, float y, float z, float w)
	{
		o = { x, y, z, w };
	}
#endif

	
	template <int D, typename T>
	struct pack_t
	{
		FORCE_INLINE	int	component_count() const {return D;}
	};

	template <typename T>
	struct pack_t<2, T>
	{
		union
		{
			struct
			{
				T x, y;
			};
			typename reg<T, 2>::type_t	m;
		};
		FORCE_INLINE	int	component_count() const {return 2;}
		FORCE_INLINE	pack_t(){}
		FORCE_INLINE	pack_t(T v):x(v), y(v){}
		FORCE_INLINE	pack_t(T ix, T iy):x(ix), y(iy){}
		FORCE_INLINE	pack_t(typename reg<T, 2>::type_t r): m(r){}
		template <typename T2>
		FORCE_INLINE	explicit pack_t(const pack_t<2, T2>& v)		{for (int i = 0; i < component_count(); i++)channel(i) = T(v[i]);}

		FORCE_INLINE	auto	channel(int i) -> T& 				{assert(0 <= i && i < component_count()); return ((T*)this)[i];}
		FORCE_INLINE	auto	channel(int i) const -> const T&	{assert(0 <= i && i < component_count()); return ((T*)this)[i];}

		DEFINE_VECTOR_OP_EQ(2)
	};
	
	template <typename T>
	struct pack_t<3, T>
	{
		union
		{
			struct
			{
				T x, y, z;
			};
			typename reg<T, 3>::type_t	m;
		};
		FORCE_INLINE	int	component_count() const {return 3;}
		FORCE_INLINE	pack_t(){}
		FORCE_INLINE	pack_t(T v):x(v), y(v), z(v){}
		FORCE_INLINE	pack_t(T ix, T iy, T iz):x(ix), y(iy), z(iz){}
		FORCE_INLINE	pack_t(const pack_t<2, T>& v, T iz) {create_register(m, v.x, v.y, iz);}
		FORCE_INLINE	pack_t(typename reg<T, 3>::type_t r): m(r){}
		template <typename T2>
		FORCE_INLINE	explicit pack_t(const pack_t<3, T2>& v)		{for (int i = 0; i < component_count(); i++)channel(i) = T(v[i]);}

		FORCE_INLINE	auto	channel(int i) -> T& 				{assert(0 <= i && i < component_count()); return ((T*)this)[i];}
		FORCE_INLINE	auto	channel(int i) const -> const T&	{assert(0 <= i && i < component_count()); return ((T*)this)[i];}

		DEFINE_VECTOR_OP_EQ(3)
	};
	
	template <typename T>
	struct pack_t<4, T>
	{
		union
		{
			struct
			{
				T x, y, z, w;
			};
			typename reg<T, 4>::type_t	m;
		};
		FORCE_INLINE	int	component_count() const {return 4;}
		FORCE_INLINE	pack_t(){}
		FORCE_INLINE	pack_t(T v) {create_register(m, v, v, v, v);}
		FORCE_INLINE	pack_t(T ix, T iy, T iz, T iw) {create_register(m, ix, iy, iz, iw);}
		FORCE_INLINE	pack_t(const pack_t<3, T>& v, T iw) {create_register(m, v.x, v.y, v.z, iw);}
		SIMD_INL_UV		pack_t(const pack_t<3, U>& v, V iw) {channel(0) = T(v.x); channel(1) = T(v.y); channel(2) = T(v.z); channel(3) = T(iw);}
		FORCE_INLINE	pack_t(typename reg<T, 4>::type_t r): m(r){}
		template <typename T2>
		FORCE_INLINE	explicit pack_t(const pack_t<4, T2>& v)		{for (int i = 0; i < component_count(); i++)channel(i) = T(v[i]);}

		FORCE_INLINE	auto	channel(int i) -> T& 				{assert(0 <= i && i < component_count()); return ((T*)this)[i];}
		FORCE_INLINE	auto	channel(int i) const -> const T&	{assert(0 <= i && i < component_count()); return ((T*)this)[i];}

		DEFINE_VECTOR_OP_EQ(4)
	};
	
	template <int D, typename T>	struct length_t<pack_t<D, T>>		{typedef	double	type;};
	template <int D>				struct length_t<pack_t<D, float>>	{typedef	float	type;};
	
	template <int D, typename T>	struct lengthi_t<pack_t<D, T>>		{typedef	T	type;};
	template <int D>				struct lengthi_t<pack_t<D, float>>	{typedef	float	type;};
	template <int D>				struct lengthi_t<pack_t<D, double>>	{typedef	double	type;};
	
	SIMD_DECLARE_EXTERN_PACK(int8_t)
	SIMD_DECLARE_EXTERN_PACK(uint8_t)
	SIMD_DECLARE_EXTERN_PACK(int16_t)
	SIMD_DECLARE_EXTERN_PACK(uint16_t)
	SIMD_DECLARE_EXTERN_PACK(int32_t)
	SIMD_DECLARE_EXTERN_PACK(uint32_t)
	SIMD_DECLARE_EXTERN_PACK(int64_t)
	SIMD_DECLARE_EXTERN_PACK(uint64_t)
	SIMD_DECLARE_EXTERN_PACK(nint8)
	SIMD_DECLARE_EXTERN_PACK(unint8)
	SIMD_DECLARE_EXTERN_PACK(nint16)
	SIMD_DECLARE_EXTERN_PACK(unint16)
	SIMD_DECLARE_EXTERN_PACK(nint32)
	SIMD_DECLARE_EXTERN_PACK(unint32)
	SIMD_DECLARE_EXTERN_PACK(nint64)
	SIMD_DECLARE_EXTERN_PACK(unint64)
	SIMD_DECLARE_EXTERN_PACK(float)
	SIMD_DECLARE_EXTERN_PACK(double)
	
	SIMD_DEFINE_PACK_T_P(uint8_t, ui8)
	SIMD_DEFINE_PACK_T_P(int8_t, i8)
	SIMD_DEFINE_PACK_T_P(uint16_t, ui16)
	SIMD_DEFINE_PACK_T_P(int16_t, i16)
	SIMD_DEFINE_PACK_T_P(uint32_t, ui32)
	SIMD_DEFINE_PACK_T_P(int32_t, i32)
	SIMD_DEFINE_PACK_T_P(uint64_t, ui64)
	SIMD_DEFINE_PACK_T_P(int64_t, i64)
	SIMD_DEFINE_PACK_T_P(float, f32)
	SIMD_DEFINE_PACK_T_P(double, f64)

}

#include "inline/core_vec_op.h"


namespace simd
{
	template <int C, int R, typename T>
	struct mat_t
	{
	public:
		typedef pack_t<R, T> col_type;
		typedef pack_t<C, T> row_type;
	public:
		col_type	internal_vector[C];
	public:
		static	FORCE_INLINE	int	row_count() 	{return R;}
		static	FORCE_INLINE	int	column_count()	{return C;}
		FORCE_INLINE	auto	operator [] (int Index) -> col_type& 					{assert(0 <= Index && Index < C); return internal_vector[Index];};
		FORCE_INLINE	auto	operator [] (int Index) const -> const col_type& 		{assert(0 <= Index && Index < C); return internal_vector[Index];};
		FORCE_INLINE	auto	operator () (int AColumn, int ARow) -> T& 				{assert(0 <= AColumn && AColumn < C && 0 <= ARow && ARow < R); return internal_vector[AColumn][ARow];};
		FORCE_INLINE	auto	operator () (int AColumn, int ARow) const -> const T& 	{assert(0 <= AColumn && AColumn < C && 0 <= ARow && ARow < R); return internal_vector[AColumn][ARow];};

		FORCE_INLINE 			mat_t() {}
								mat_t(mat_t<C, R, T> const& m)							{for (int c = 0; c < C; c++) internal_vector[c] = m[c];}
		SIMD_TMAT_EXPL_CRU		mat_t(mat_t<C, R, U> const& m)							{for (int c = 0; c < C; c++) internal_vector[c] = m[c];}
		explicit				mat_t(const T& x);
								mat_t<C, R, T>&	operator = (mat_t<C, R, T> const& m)	{for (int c = 0; c < C; c++) internal_vector[c] = m[c];return *this;}
		template<typename U>	mat_t<C, R, T>&	operator = (mat_t<C, R, U> const& m)	{for (int c = 0; c < C; c++) internal_vector[c] = m[c];return *this;}
		template<typename U>	mat_t<C, R, T>&	operator += (U s);
		template<typename U>	mat_t<C, R, T>&	operator += (mat_t<C, R, U> const& m);
		template<typename U>	mat_t<C, R, T>&	operator -= (U s);
		template<typename U>	mat_t<C, R, T>&	operator -= (mat_t<C, R, U> const& m);
		template<typename U>	mat_t<C, R, T>&	operator *= (U s);
		template<typename U>	mat_t<C, R, T>&	operator *= (mat_t<C, R, U> const& m);
		template<typename U>	mat_t<C, R, T>&	operator /= (U s);
		template<typename U>	mat_t<C, R, T>&	operator /= (mat_t<C, R, U> const& m);
		
		mat_t<C, R, T>&	operator++();
		mat_t<C, R, T>&	operator--();
		mat_t<C, R, T>	operator++(int);
		mat_t<C, R, T>	operator--(int);
		
		auto	row(int Index) const -> row_type;
		FORCE_INLINE	auto	column(int Index) const -> const col_type& {return internal_vector[Index];}
	};

	template<int C, int R, typename T>	auto	operator + (mat_t<C, R, T> const& m) -> mat_t<C, R, T>;
	template<int C, int R, typename T>	auto	operator - (mat_t<C, R, T> const& m) -> mat_t<C, R, T>;
	
	template<int C, int R, typename T>	auto	operator + (mat_t<C, R, T> const& m,	T const& s) -> mat_t<C, R, T>;
	template<int C, int R, typename T>	auto	operator - (mat_t<C, R, T> const& m, T const& s) -> mat_t<C, R, T>;
	template<int C, int R, typename T>	auto	operator * (mat_t<C, R, T> const& m, T const& s) -> mat_t<C, R, T>;
	template<int C, int R, typename T>	auto	operator / (mat_t<C, R, T> const& m, T const& s) -> mat_t<C, R, T>;

	template<int C, int R, typename T>	auto	operator + (T const& s, mat_t<C, R, T> const& m) -> mat_t<C, R, T>;
	template<int C, int R, typename T>	auto	operator - (T const& s, mat_t<C, R, T> const& m) -> mat_t<C, R, T>;
	template<int C, int R, typename T>	auto	operator * (T const& s, mat_t<C, R, T> const& m) -> mat_t<C, R, T>;
	template<int C, int R, typename T>	auto	operator / (T const& s, mat_t<C, R, T> const& m) -> mat_t<C, R, T>;

	template<int C, int R, typename T>	auto	operator + (mat_t<C, R, T> const& m1, mat_t<C, R, T> const& m2) -> mat_t<C, R, T>;
	template<int C, int R, typename T>	auto	operator - (mat_t<C, R, T> const& m1, mat_t<C, R, T> const& m2) -> mat_t<C, R, T>;
	template<int C, int R, typename T>	auto	operator * (mat_t<C, R, T> const& m1, mat_t<C, R, T> const& m2) -> mat_t<C, R, T>;
	template<int C, int R, typename T>	auto	operator / (mat_t<C, R, T> const& m1, mat_t<C, R, T> const& m2) -> mat_t<C, R, T>;

	template<int C, int R, typename T>	auto	operator * (const mat_t<C, R, T>& m, const pack_t<C, T>& v) -> pack_t<R, T>;
	template<int C, int R, typename T>	auto	operator / (mat_t<C, R, T> const& m, typename mat_t<C, R, T>::row_type const& v) -> typename mat_t<C, R, T>::col_type;

	template<int C, int R, typename T>	auto	operator * (typename mat_t<C, R, T>::col_type const& v, mat_t<C, R, T> const& m) -> typename mat_t<C, R, T>::row_type;
	template<int C, int R, typename T>	auto	operator / (typename mat_t<C, R, T>::col_type const& v, mat_t<C, R, T> const& m) -> typename mat_t<C, R, T>::row_type;
	
	template<int C, int R, typename T>	bool	operator == (mat_t<C, R, T> const& m1, mat_t<C, R, T> const& m2);
	template<int C, int R, typename T>	bool	operator != (mat_t<C, R, T> const& m1, mat_t<C, R, T> const& m2);

	template<typename T>	auto	operator * (mat_t<4, 4, T> const& m1, mat_t<2, 4, T> const& m2) -> mat_t<2, 4, T>;
	template<typename T>	auto	operator * (mat_t<4, 4, T> const& m1, mat_t<3, 4, T> const& m2) -> mat_t<3, 4, T>;

	template<int C, int R, typename T>	mat_t<R, C, T>	transpose(mat_t<C, R, T> const& m);
	template<int C, int R, typename T>	mat_t<R, C, T>	inverse(mat_t<C, R, T> const& m);
	template<int C, int R, typename T>	mat_t<R, C, T>	inverse_fast(mat_t<C, R, T> const& m);
	template<int C, int R, typename T>	mat_t<R, C, T>	inverse_homogeneous(mat_t<C, R, T> const& m);
	
	SIMD_DECLARE_EXTERN_MAT(float)
	SIMD_DECLARE_EXTERN_MAT(double)
}



















namespace cl
{
#if defined(__CAN_USE_QUAD_FLOATING_PRECISSION)
	//typedef __float128 quad;
#endif
	typedef uint64_t	ulong;
	typedef int64_t		slong;
	typedef uint32_t	uint;
	typedef int32_t		sint;
	typedef uint16_t	ushort;
	typedef int16_t		sshort;
	typedef uint8_t		uchar;
	typedef int8_t		schar;

#if defined(__CAN_USE_QUAD_FLOATING_PRECISSION)
	typedef simd::pack_t<quad,16, false>	quad16;
	typedef simd::pack_t<quad,8, false>	quad8;
	typedef simd::pack_t<quad,4, false>	quad4;
	typedef simd::pack_t<quad,3, false>	quad3;
	typedef simd::pack_t<quad,2, false>	quad2;
#endif
	
#define TYPEDEF_CL_T_P(_T, _P) \
typedef simd::pack_t<16, _T >		_P##16; \
typedef simd::pack_t<8, _T >		_P##8; \
typedef simd::pack_t<4, _T >		_P##4; \
typedef simd::pack_t<3, _T >		_P##3; \
typedef simd::pack_t<2, _T >		_P##2;

	TYPEDEF_CL_T_P(double, double)
	TYPEDEF_CL_T_P(float, float)
	TYPEDEF_CL_T_P(ulong, ulong)
	TYPEDEF_CL_T_P(slong, long)
	TYPEDEF_CL_T_P(uint, uint)
	TYPEDEF_CL_T_P(sint, int)
	TYPEDEF_CL_T_P(ushort, ushort)
	TYPEDEF_CL_T_P(sshort, short)
	TYPEDEF_CL_T_P(uchar, uchar)
	TYPEDEF_CL_T_P(schar, char)
}


namespace glsl
{
#define DEFINE_TYPEDEF_VEC(_T, _pre) \
	typedef	simd::pack_t<2, _T>	_pre##vec2; \
	typedef	simd::pack_t<3, _T>	_pre##vec3; \
	typedef	simd::pack_t<4, _T>	_pre##vec4;
	
	DEFINE_TYPEDEF_VEC(int32_t, i)
	DEFINE_TYPEDEF_VEC(uint32_t, ui)
	DEFINE_TYPEDEF_VEC(float, )
	DEFINE_TYPEDEF_VEC(float, f32)
	DEFINE_TYPEDEF_VEC(double, f64)
	
#define DEFINE_TYPEDEF_MAT(_T, _pre) \
	typedef	simd::mat_t<4, 4, _T>	_pre##mat4;
	
	DEFINE_TYPEDEF_MAT(float, )
	DEFINE_TYPEDEF_MAT(float, f32)
	DEFINE_TYPEDEF_MAT(double, f64)

	
	inline	simd::mat_t<4, 4, float>	mat4_make(simd::pack_t<4, float> const& v1, simd::pack_t<4, float> const& v2, simd::pack_t<4, float> const& v3, simd::pack_t<4, float> const& v4)
	{
		simd::mat_t<4, 4, float> ret;
		ret.internal_vector[0] = v1;
		ret.internal_vector[1] = v2;
		ret.internal_vector[2] = v3;
		ret.internal_vector[3] = v4;
		return ret;
	}
	
	inline	simd::mat_t<4, 4, float>	mat4_make(	float const& x0, float const& y0, float const& z0, float const& w0,
												  	float const& x1, float const& y1, float const& z1, float const& w1,
													float const& x2, float const& y2, float const& z2, float const& w2,
													float const& x3, float const& y3, float const& z3, float const& w3)
	{
		return mat4_make(simd::pack4_f32(x0, y0, z0, w0), simd::pack4_f32(x1, y1, z1, w1), simd::pack4_f32(x2, y2, z2, w2), simd::pack4_f32(x3, y3, z3, w3));
	}
}

namespace hlsl
{
	
}

namespace cl
{
	
}


