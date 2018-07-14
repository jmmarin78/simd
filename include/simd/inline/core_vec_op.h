// This software is MIT licensed (see LICENSE)

namespace simd
{
#if CPU_ARCH & CPU_ARCH_SSE2_BIT
	
	FORCE_INLINE SIMD_vec4 SIMD_vec1_add(SIMD_vec4 a, SIMD_vec4 b)
	{
		return _mm_add_ss(a, b);
	}
	
	
	FORCE_INLINE SIMD_vec4 SIMD_vec1_sub(SIMD_vec4 a, SIMD_vec4 b)
	{
		return _mm_sub_ss(a, b);
	}
	
	FORCE_INLINE SIMD_vec4 SIMD_vec1_mul(SIMD_vec4 a, SIMD_vec4 b)
	{
		return _mm_mul_ss(a, b);
	}
	
	FORCE_INLINE SIMD_vec4 SIMD_vec1_div(SIMD_vec4 a, SIMD_vec4 b)
	{
		return _mm_div_ss(a, b);
	}
	
	FORCE_INLINE SIMD_vec4 SIMD_vec4_swizzle_xyzw(SIMD_vec4 a)
	{
#	if CPU_ARCH & CPU_ARCH_AVX2_BIT
		return _mm_permute_ps(a, _MM_SHUFFLE(3, 2, 1, 0));
#	else
		return _mm_shuffle_ps(a, a, _MM_SHUFFLE(3, 2, 1, 0));
#	endif
	}
	
	FORCE_INLINE SIMD_vec4 SIMD_vec1_fma(SIMD_vec4 a, SIMD_vec4 b, SIMD_vec4 c)
	{
#	if (CPU_ARCH & CPU_ARCH_AVX2_BIT) && !(SIMD_COMPILER & SIMD_COMPILER_CLANG)
		return _mm_fmadd_ss(a, b, c);
#	else
		return _mm_add_ss(_mm_mul_ss(a, b), c);
#	endif
	}
	
//	FORCE_INLINE SIMD_vec4 SIMD_vec4_abs(SIMD_vec4 x)
//	{
//		return _mm_and_ps(x, _mm_castsi128_ps(_mm_set1_epi32(0x7FFFFFFF)));
//	}
	extern __m128 SIGNMASK_VEC4;
	FORCE_INLINE	SIMD_ivec4 SIMD_vec4_abs(SIMD_ivec4 x)
	{
#	if CPU_ARCH & CPU_ARCH_SSSE3_BIT
		return _mm_andnot_ps(SIGNMASK_VEC4, x);
#	endif
	}
	
	FORCE_INLINE	SIMD_ivec4 SIMD_vec4_neg(SIMD_ivec4 x)
	{
#	if CPU_ARCH & CPU_ARCH_SSSE3_BIT
		return _mm_xor_ps(x, SIGNMASK_VEC4);
#	endif
	}
	
	DLL_FUNCTION(SIMD_vec4) SIMD_vec4_sign(SIMD_vec4 x);
	DLL_FUNCTION(SIMD_vec4) SIMD_vec4_dot(SIMD_vec4 v1, SIMD_vec4 v2);
	
	FORCE_INLINE SIMD_vec4 SIMD_vec1_sqrt_lowp(SIMD_vec4 x)
	{
		return _mm_mul_ss(_mm_rsqrt_ss(x), x);
	}
	
	FORCE_INLINE SIMD_vec4 SIMD_vec4_sqrt_lowp(SIMD_vec4 x)
	{
		return _mm_mul_ps(_mm_rsqrt_ps(x), x);
	}
#endif//CPU_ARCH & CPU_ARCH_SSE2_BIT

}

#define SIMD_DECLARE_PACK_OPERATOR(_op) \
	template <int D, typename T>				DLL_APICALL pack_t<D, T>	DLL_CALLING_CONVENTION	operator _op (const pack_t<D, T>&, const pack_t<D, T>&);\
	template <int D, typename T, typename T2>	DLL_APICALL pack_t<D, T>	DLL_CALLING_CONVENTION 	operator _op (const T2&, const pack_t<D, T>&);\
	template <int D, typename T, typename T2>	DLL_APICALL pack_t<D, T>	DLL_CALLING_CONVENTION 	operator _op (const pack_t<D, T>&, const T2&);
#define SIMD_DEFINE_PACK4_OPERATOR(_T, _O, _F) \
	template <bool N> \
	FORCE_INLINE	pack_t<4, _T>	operator _O (const pack_t<4, _T>& a, const pack_t<4, _T>& b) \
	{ \
		return _F(a.m, b.m); \
	}

namespace simd
{
	SIMD_DECLARE_PACK_OPERATOR(*)
	SIMD_DECLARE_PACK_OPERATOR(+)
	SIMD_DECLARE_PACK_OPERATOR(-)
	SIMD_DECLARE_PACK_OPERATOR(/)
	
	SIMD_DLLFNC_DT(auto)	isnan(const pack_t<D, T>&) -> pack_t<D, bool>;
	SIMD_DLLFNC_DT(auto)	isinf(const pack_t<D, T>&) -> pack_t<D, bool>;
	SIMD_DLLFNC_D(auto)		floatBitsToInt(const pack_t<D, float>&) -> pack_t<D, int>;
	SIMD_DLLFNC_D(auto)		floatBitsToUint(const pack_t<D, float>&) -> pack_t<D, uint32_t>;
	SIMD_DLLFNC_D(auto)		intBitsToFloat(const pack_t<D, int>&) -> pack_t<D, float>;
	SIMD_DLLFNC_D(auto)		uintBitsToFloat(const pack_t<D, unsigned int>&) -> pack_t<D, float>;
	SIMD_DLLFNC_T2(auto)	cross(const pack_t<3, T1>&, const pack_t<3, T2>&) -> pack_t<3, typename high_precission_t<T1, T2>::type>;
	SIMD_DLLFNC_T2(auto)	cross(const pack_t<2, T1>&, const pack_t<2, T2>&) -> typename high_precission_t<T1, T2>::type;
}

namespace simd
{
#if CPU_ARCH & CPU_ARCH_SSE2_BIT
	FORCE_INLINE	pack_t<4, float> add (const pack_t<4, float>& a, const pack_t<4, float>& b)			{return _mm_add_ps(a.m, b.m);}
	FORCE_INLINE	pack_t<4, float> sub (const pack_t<4, float>& a, const pack_t<4, float>& b)			{return _mm_sub_ps(a.m, b.m);}
	FORCE_INLINE	pack_t<4, float> mul (const pack_t<4, float>& a, const pack_t<4, float>& b)			{return _mm_mul_ps(a.m, b.m);}
	FORCE_INLINE	pack_t<4, float> div (const pack_t<4, float>& a, const pack_t<4, float>& b)			{return _mm_div_ps(a.m, b.m);}
	FORCE_INLINE	pack_t<4, float> sqrt(const pack_t<4, float>& a)									{return _mm_sqrt_ps(a.m);}
	FORCE_INLINE	pack_t<4, float> rcp(const pack_t<4, float>& a)										{return _mm_rcp_ps(a.m);}
	FORCE_INLINE	pack_t<4, float> rsq(const pack_t<4, float>& a)										{return _mm_rsqrt_ps(a.m);}
	FORCE_INLINE	pack_t<4, float> min(const pack_t<4, float>& a, const pack_t<4, float>& b)			{return _mm_min_ps(a.m, b.m);}
	FORCE_INLINE	pack_t<4, float> max(const pack_t<4, float>& a, const pack_t<4, float>& b)			{return _mm_max_ps(a.m, b.m);}

	FORCE_INLINE	pack_t<4, float> operator + (const pack_t<4, float>& a, const pack_t<4, float>& b)	{return _mm_add_ps(a.m, b.m);}
	FORCE_INLINE	pack_t<4, float> operator - (const pack_t<4, float>& a, const pack_t<4, float>& b)	{return _mm_sub_ps(a.m, b.m);}
	FORCE_INLINE	pack_t<4, float> operator * (const pack_t<4, float>& a, const pack_t<4, float>& b)	{return _mm_mul_ps(a.m, b.m);}
	FORCE_INLINE	pack_t<4, float> operator / (const pack_t<4, float>& a, const pack_t<4, float>& b)	{return _mm_div_ps(a.m, b.m);}
	
	FORCE_INLINE	pack_t<4, float> div_fast(const pack_t<4, float>& a, const pack_t<4, float>& b)		{return _mm_mul_ps(a.m, _mm_rcp_ps(b.m));}
	FORCE_INLINE	pack_t<4, float> abs(const pack_t<4, float>& a)										{return SIMD_vec4_abs(a.m);}
	FORCE_INLINE	pack_t<4, float> sign(const pack_t<4, float>& a)									{return SIMD_vec4_sign(a.m);}
	FORCE_INLINE	pack_t<4, float> sign(const pack_t<4, float>& a, const pack_t<4, float>& b)			{return SIMD_vec4_dot(a.m, b.m);}

#endif
	
#if (CPU_ARCH & CPU_ARCH_AVX2_BIT) && !(SIMD_COMPILER & SIMD_COMPILER_CLANG)
	FORCE_INLINE pack_t<4, float> mad(const pack_t<4, float>& a, const pack_t<4, float>& b, const pack_t<4, float>& c)
	{
		return _mm_fmadd_ps(a.m, b.m, c.m);
	}
#endif

}





