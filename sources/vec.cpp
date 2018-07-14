// This software is MIT licensed (see LICENSE)

#define BUILD_DLL

#include <math.h>

#include <pre/number.h>
#include <simd/core.h>

#include "defines.h"

namespace simd
{
	IMPLEMENT_VECTOR_OP_EQ(2)
	IMPLEMENT_VECTOR_OP_EQ(3)
	IMPLEMENT_VECTOR_OP_EQ(4)
}



namespace simd
{
#if CPU_ARCH & CPU_ARCH_SSE2_BIT
	
	SIMD_uvec4 SIMD_i128_interleave(SIMD_uvec4 x)
	{
		SIMD_uvec4 const Mask4 = _mm_set1_epi32(0x0000FFFF);
		SIMD_uvec4 const Mask3 = _mm_set1_epi32(0x00FF00FF);
		SIMD_uvec4 const Mask2 = _mm_set1_epi32(0x0F0F0F0F);
		SIMD_uvec4 const Mask1 = _mm_set1_epi32(0x33333333);
		SIMD_uvec4 const Mask0 = _mm_set1_epi32(0x55555555);
		
		SIMD_uvec4 Reg1;
		SIMD_uvec4 Reg2;
		
		// REG1 = x;
		// REG2 = y;
		//Reg1 = _mm_unpacklo_epi64(x, y);
		Reg1 = x;
		
		//REG1 = ((REG1 << 16) | REG1) & SML::uint64(0x0000FFFF0000FFFF);
		//REG2 = ((REG2 << 16) | REG2) & SML::uint64(0x0000FFFF0000FFFF);
		Reg2 = _mm_slli_si128(Reg1, 2);
		Reg1 = _mm_or_si128(Reg2, Reg1);
		Reg1 = _mm_and_si128(Reg1, Mask4);
		
		//REG1 = ((REG1 <<  8) | REG1) & SML::uint64(0x00FF00FF00FF00FF);
		//REG2 = ((REG2 <<  8) | REG2) & SML::uint64(0x00FF00FF00FF00FF);
		Reg2 = _mm_slli_si128(Reg1, 1);
		Reg1 = _mm_or_si128(Reg2, Reg1);
		Reg1 = _mm_and_si128(Reg1, Mask3);
		
		//REG1 = ((REG1 <<  4) | REG1) & SML::uint64(0x0F0F0F0F0F0F0F0F);
		//REG2 = ((REG2 <<  4) | REG2) & SML::uint64(0x0F0F0F0F0F0F0F0F);
		Reg2 = _mm_slli_epi32(Reg1, 4);
		Reg1 = _mm_or_si128(Reg2, Reg1);
		Reg1 = _mm_and_si128(Reg1, Mask2);
		
		//REG1 = ((REG1 <<  2) | REG1) & SML::uint64(0x3333333333333333);
		//REG2 = ((REG2 <<  2) | REG2) & SML::uint64(0x3333333333333333);
		Reg2 = _mm_slli_epi32(Reg1, 2);
		Reg1 = _mm_or_si128(Reg2, Reg1);
		Reg1 = _mm_and_si128(Reg1, Mask1);
		
		//REG1 = ((REG1 <<  1) | REG1) & SML::uint64(0x5555555555555555);
		//REG2 = ((REG2 <<  1) | REG2) & SML::uint64(0x5555555555555555);
		Reg2 = _mm_slli_epi32(Reg1, 1);
		Reg1 = _mm_or_si128(Reg2, Reg1);
		Reg1 = _mm_and_si128(Reg1, Mask0);
		
		//return REG1 | (REG2 << 1);
		Reg2 = _mm_slli_epi32(Reg1, 1);
		Reg2 = _mm_srli_si128(Reg2, 8);
		Reg1 = _mm_or_si128(Reg1, Reg2);
		
		return Reg1;
	}
	
	SIMD_uvec4 SIMD_i128_interleave2(SIMD_uvec4 x, SIMD_uvec4 y)
	{
		SIMD_uvec4 const Mask4 = _mm_set1_epi32(0x0000FFFF);
		SIMD_uvec4 const Mask3 = _mm_set1_epi32(0x00FF00FF);
		SIMD_uvec4 const Mask2 = _mm_set1_epi32(0x0F0F0F0F);
		SIMD_uvec4 const Mask1 = _mm_set1_epi32(0x33333333);
		SIMD_uvec4 const Mask0 = _mm_set1_epi32(0x55555555);
		
		SIMD_uvec4 Reg1;
		SIMD_uvec4 Reg2;
		
		// REG1 = x;
		// REG2 = y;
		Reg1 = _mm_unpacklo_epi64(x, y);
		
		//REG1 = ((REG1 << 16) | REG1) & SML::uint64(0x0000FFFF0000FFFF);
		//REG2 = ((REG2 << 16) | REG2) & SML::uint64(0x0000FFFF0000FFFF);
		Reg2 = _mm_slli_si128(Reg1, 2);
		Reg1 = _mm_or_si128(Reg2, Reg1);
		Reg1 = _mm_and_si128(Reg1, Mask4);
		
		//REG1 = ((REG1 <<  8) | REG1) & SML::uint64(0x00FF00FF00FF00FF);
		//REG2 = ((REG2 <<  8) | REG2) & SML::uint64(0x00FF00FF00FF00FF);
		Reg2 = _mm_slli_si128(Reg1, 1);
		Reg1 = _mm_or_si128(Reg2, Reg1);
		Reg1 = _mm_and_si128(Reg1, Mask3);
		
		//REG1 = ((REG1 <<  4) | REG1) & SML::uint64(0x0F0F0F0F0F0F0F0F);
		//REG2 = ((REG2 <<  4) | REG2) & SML::uint64(0x0F0F0F0F0F0F0F0F);
		Reg2 = _mm_slli_epi32(Reg1, 4);
		Reg1 = _mm_or_si128(Reg2, Reg1);
		Reg1 = _mm_and_si128(Reg1, Mask2);
		
		//REG1 = ((REG1 <<  2) | REG1) & SML::uint64(0x3333333333333333);
		//REG2 = ((REG2 <<  2) | REG2) & SML::uint64(0x3333333333333333);
		Reg2 = _mm_slli_epi32(Reg1, 2);
		Reg1 = _mm_or_si128(Reg2, Reg1);
		Reg1 = _mm_and_si128(Reg1, Mask1);
		
		//REG1 = ((REG1 <<  1) | REG1) & SML::uint64(0x5555555555555555);
		//REG2 = ((REG2 <<  1) | REG2) & SML::uint64(0x5555555555555555);
		Reg2 = _mm_slli_epi32(Reg1, 1);
		Reg1 = _mm_or_si128(Reg2, Reg1);
		Reg1 = _mm_and_si128(Reg1, Mask0);
		
		//return REG1 | (REG2 << 1);
		Reg2 = _mm_slli_epi32(Reg1, 1);
		Reg2 = _mm_srli_si128(Reg2, 8);
		Reg1 = _mm_or_si128(Reg1, Reg2);
		
		return Reg1;
	}
	
#endif//CPU_ARCH & CPU_ARCH_SSE2_BIT
	
#if CPU_ARCH & CPU_ARCH_SSE2_BIT
	
	FORCE_INLINE SIMD_vec4 SIMD_vec4_fma(SIMD_vec4 a, SIMD_vec4 b, SIMD_vec4 c)
	{
		return _mm_add_ps(_mm_mul_ps(a, b), c);
	}

	
	DLL_FUNCTION(SIMD_vec4) SIMD_vec4_dot(SIMD_vec4 v1, SIMD_vec4 v2)
	{
#	if CPU_ARCH & CPU_ARCH_AVX_BIT
		return _mm_dp_ps(v1, v2, 0xff);
#	elif CPU_ARCH & CPU_ARCH_SSE3_BIT
		SIMD_vec4 const mul0 = _mm_mul_ps(v1, v2);
		SIMD_vec4 const hadd0 = _mm_hadd_ps(mul0, mul0);
		SIMD_vec4 const hadd1 = _mm_hadd_ps(hadd0, hadd0);
		return hadd1;
#	else
		SIMD_vec4 const mul0 = _mm_mul_ps(v1, v2);
		SIMD_vec4 const swp0 = _mm_shuffle_ps(mul0, mul0, _MM_SHUFFLE(2, 3, 0, 1));
		SIMD_vec4 const add0 = _mm_add_ps(mul0, swp0);
		SIMD_vec4 const swp1 = _mm_shuffle_ps(add0, add0, _MM_SHUFFLE(0, 1, 2, 3));
		SIMD_vec4 const add1 = _mm_add_ps(add0, swp1);
		return add1;
#	endif
	}
	
	SIMD_vec4 SIMD_vec4_length(SIMD_vec4 x)
	{
		SIMD_vec4 const dot0 = SIMD_vec4_dot(x, x);
		SIMD_vec4 const sqt0 = _mm_sqrt_ps(dot0);
		return sqt0;
	}
	
	SIMD_vec4 SIMD_vec4_distance(SIMD_vec4 p0, SIMD_vec4 p1)
	{
		SIMD_vec4 const sub0 = _mm_sub_ps(p0, p1);
		SIMD_vec4 const len0 = SIMD_vec4_length(sub0);
		return len0;
	}
	
	SIMD_vec4 SIMD_vec1_dot(SIMD_vec4 v1, SIMD_vec4 v2)
	{
#	if CPU_ARCH & CPU_ARCH_AVX_BIT
		return _mm_dp_ps(v1, v2, 0xff);
#	elif CPU_ARCH & CPU_ARCH_SSE3_BIT
		SIMD_vec4 const mul0 = _mm_mul_ps(v1, v2);
		SIMD_vec4 const had0 = _mm_hadd_ps(mul0, mul0);
		SIMD_vec4 const had1 = _mm_hadd_ps(had0, had0);
		return had1;
#	else
		SIMD_vec4 const mul0 = _mm_mul_ps(v1, v2);
		SIMD_vec4 const mov0 = _mm_movehl_ps(mul0, mul0);
		SIMD_vec4 const add0 = _mm_add_ps(mov0, mul0);
		SIMD_vec4 const swp1 = _mm_shuffle_ps(add0, add0, 1);
		SIMD_vec4 const add1 = _mm_add_ss(add0, swp1);
		return add1;
#	endif
	}
		
	SIMD_vec4 SIMD_vec4_normalize(SIMD_vec4 v)
	{
		SIMD_vec4 const dot0 = SIMD_vec4_dot(v, v);
		SIMD_vec4 const isr0 = _mm_rsqrt_ps(dot0);
		SIMD_vec4 const mul0 = _mm_mul_ps(v, isr0);
		return mul0;
	}
	
	SIMD_vec4 SIMD_vec4_faceforward(SIMD_vec4 N, SIMD_vec4 I, SIMD_vec4 Nref)
	{
		SIMD_vec4 const dot0 = SIMD_vec4_dot(Nref, I);
		SIMD_vec4 const sgn0 = SIMD_vec4_sign(dot0);
		SIMD_vec4 const mul0 = _mm_mul_ps(sgn0, _mm_set1_ps(-1.0f));
		SIMD_vec4 const mul1 = _mm_mul_ps(N, mul0);
		return mul1;
	}
	
	SIMD_vec4 SIMD_vec4_reflect(SIMD_vec4 I, SIMD_vec4 N)
	{
		SIMD_vec4 const dot0 = SIMD_vec4_dot(N, I);
		SIMD_vec4 const mul0 = _mm_mul_ps(N, dot0);
		SIMD_vec4 const mul1 = _mm_mul_ps(mul0, _mm_set1_ps(2.0f));
		SIMD_vec4 const sub0 = _mm_sub_ps(I, mul1);
		return sub0;
	}
	
	__m128 SIMD_vec4_refract(SIMD_vec4 I, SIMD_vec4 N, SIMD_vec4 eta)
	{
		SIMD_vec4 const dot0 = SIMD_vec4_dot(N, I);
		SIMD_vec4 const mul0 = _mm_mul_ps(eta, eta);
		SIMD_vec4 const mul1 = _mm_mul_ps(dot0, dot0);
		SIMD_vec4 const sub0 = _mm_sub_ps(_mm_set1_ps(1.0f), mul0);
		SIMD_vec4 const sub1 = _mm_sub_ps(_mm_set1_ps(1.0f), mul1);
		SIMD_vec4 const mul2 = _mm_mul_ps(sub0, sub1);
		
		if(_mm_movemask_ps(_mm_cmplt_ss(mul2, _mm_set1_ps(0.0f))) == 0)
			return _mm_set1_ps(0.0f);
		
		SIMD_vec4 const sqt0 = _mm_sqrt_ps(mul2);
		SIMD_vec4 const mad0 = SIMD_vec4_fma(eta, dot0, sqt0);
		SIMD_vec4 const mul4 = _mm_mul_ps(mad0, N);
		SIMD_vec4 const mul5 = _mm_mul_ps(eta, I);
		SIMD_vec4 const sub2 = _mm_sub_ps(mul5, mul4);
		
		return sub2;
	}
	
#endif//CPU_ARCH & CPU_ARCH_SSE2_BIT
	
	
	
	
	
	
	
	
	
	
	
	
	
	
#if CPU_ARCH & CPU_ARCH_SSE2_BIT
	
	SIMD_ivec4 SIMD_ivec4_abs(SIMD_ivec4 x)
	{
#	if CPU_ARCH & CPU_ARCH_SSSE3_BIT
		return _mm_sign_epi32(x, x);
#	else
		SIMD_ivec4 const sgn0 = _mm_srai_epi32(x, 31);
		SIMD_ivec4 const inv0 = _mm_xor_si128(x, sgn0);
		SIMD_ivec4 const sub0 = _mm_sub_epi32(inv0, sgn0);
		return sub0;
#	endif
	}
	
	__m128 SIGNMASK_VEC4 = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
	
	
	DLL_FUNCTION(SIMD_vec4) SIMD_vec4_sign(SIMD_vec4 x)
	{
		SIMD_vec4 const zro0 = _mm_setzero_ps();
		SIMD_vec4 const cmp0 = _mm_cmplt_ps(x, zro0);
		SIMD_vec4 const cmp1 = _mm_cmpgt_ps(x, zro0);
		SIMD_vec4 const and0 = _mm_and_ps(cmp0, _mm_set1_ps(-1.0f));
		SIMD_vec4 const and1 = _mm_and_ps(cmp1, _mm_set1_ps(1.0f));
		SIMD_vec4 const or0 = _mm_or_ps(and0, and1);;
		return or0;
	}
	
	SIMD_vec4 SIMD_vec4_round(SIMD_vec4 x)
	{
#	if CPU_ARCH & CPU_ARCH_SSE41_BIT
		return _mm_round_ps(x, _MM_FROUND_TO_NEAREST_INT);
#	else
		SIMD_vec4 const sgn0 = _mm_castsi128_ps(_mm_set1_epi32(int(0x80000000)));
		SIMD_vec4 const and0 = _mm_and_ps(sgn0, x);
		SIMD_vec4 const or0 = _mm_or_ps(and0, _mm_set_ps1(8388608.0f));
		SIMD_vec4 const add0 = _mm_add_ps(x, or0);
		SIMD_vec4 const sub0 = _mm_sub_ps(add0, or0);
		return sub0;
#	endif
	}
	
	SIMD_vec4 SIMD_vec4_floor(SIMD_vec4 x)
	{
#	if CPU_ARCH & CPU_ARCH_SSE41_BIT
		return _mm_floor_ps(x);
#	else
		SIMD_vec4 const rnd0 = SIMD_vec4_round(x);
		SIMD_vec4 const cmp0 = _mm_cmplt_ps(x, rnd0);
		SIMD_vec4 const and0 = _mm_and_ps(cmp0, _mm_set1_ps(1.0f));
		SIMD_vec4 const sub0 = _mm_sub_ps(rnd0, and0);
		return sub0;
#	endif
	}
	
	/* trunc TODO
	 SIMD_vec4 SIMD_vec4_trunc(SIMD_vec4 x)
	 {
	 return SIMD_vec4();
	 }
	 */
	
	//roundEven
	SIMD_vec4 SIMD_vec4_roundEven(SIMD_vec4 x)
	{
		SIMD_vec4 const sgn0 = _mm_castsi128_ps(_mm_set1_epi32(int(0x80000000)));
		SIMD_vec4 const and0 = _mm_and_ps(sgn0, x);
		SIMD_vec4 const or0 = _mm_or_ps(and0, _mm_set_ps1(8388608.0f));
		SIMD_vec4 const add0 = _mm_add_ps(x, or0);
		SIMD_vec4 const sub0 = _mm_sub_ps(add0, or0);
		return sub0;
	}
	
	SIMD_vec4 SIMD_vec4_ceil(SIMD_vec4 x)
	{
#	if CPU_ARCH & CPU_ARCH_SSE41_BIT
		return _mm_ceil_ps(x);
#	else
		SIMD_vec4 const rnd0 = SIMD_vec4_round(x);
		SIMD_vec4 const cmp0 = _mm_cmpgt_ps(x, rnd0);
		SIMD_vec4 const and0 = _mm_and_ps(cmp0, _mm_set1_ps(1.0f));
		SIMD_vec4 const add0 = _mm_add_ps(rnd0, and0);
		return add0;
#	endif
	}
	
	SIMD_vec4 SIMD_vec4_fract(SIMD_vec4 x)
	{
		SIMD_vec4 const flr0 = SIMD_vec4_floor(x);
		SIMD_vec4 const sub0 = _mm_sub_ps(x, flr0);
		return sub0;
	}
	
	SIMD_vec4 SIMD_vec4_mod(SIMD_vec4 x, SIMD_vec4 y)
	{
		SIMD_vec4 const div0 = _mm_div_ps(x, y);
		SIMD_vec4 const flr0 = SIMD_vec4_floor(div0);
		SIMD_vec4 const mul0 = _mm_mul_ps(y, flr0);
		SIMD_vec4 const sub0 = _mm_sub_ps(x, mul0);
		return sub0;
	}
	
	SIMD_vec4 SIMD_vec4_clamp(SIMD_vec4 v, SIMD_vec4 minVal, SIMD_vec4 maxVal)
	{
		SIMD_vec4 const min0 = _mm_min_ps(v, maxVal);
		SIMD_vec4 const max0 = _mm_max_ps(min0, minVal);
		return max0;
	}
	
	SIMD_vec4 SIMD_vec4_mix(SIMD_vec4 v1, SIMD_vec4 v2, SIMD_vec4 a)
	{
		SIMD_vec4 const sub0 = _mm_sub_ps(_mm_set1_ps(1.0f), a);
		SIMD_vec4 const mul0 = _mm_mul_ps(v1, sub0);
		SIMD_vec4 const mad0 = SIMD_vec4_fma(v2, a, mul0);
		return mad0;
	}
	
	SIMD_vec4 SIMD_vec4_step(SIMD_vec4 edge, SIMD_vec4 x)
	{
		SIMD_vec4 const cmp = _mm_cmple_ps(x, edge);
		return _mm_movemask_ps(cmp) == 0 ? _mm_set1_ps(1.0f) : _mm_setzero_ps();
	}
	
	SIMD_vec4 SIMD_vec4_smoothstep(SIMD_vec4 edge0, SIMD_vec4 edge1, SIMD_vec4 x)
	{
		SIMD_vec4 const sub0 = _mm_sub_ps(x, edge0);
		SIMD_vec4 const sub1 = _mm_sub_ps(edge1, edge0);
		SIMD_vec4 const div0 = _mm_sub_ps(sub0, sub1);
		SIMD_vec4 const clp0 = SIMD_vec4_clamp(div0, _mm_setzero_ps(), _mm_set1_ps(1.0f));
		SIMD_vec4 const mul0 = _mm_mul_ps(_mm_set1_ps(2.0f), clp0);
		SIMD_vec4 const sub2 = _mm_sub_ps(_mm_set1_ps(3.0f), mul0);
		SIMD_vec4 const mul1 = _mm_mul_ps(clp0, clp0);
		SIMD_vec4 const mul2 = _mm_mul_ps(mul1, sub2);
		return mul2;
	}
	
	// Agner Fog method
	SIMD_vec4 SIMD_vec4_nan(SIMD_vec4 x)
	{
		SIMD_ivec4 const t1 = _mm_castps_si128(x);						// reinterpret as 32-bit integer
		SIMD_ivec4 const t2 = _mm_sll_epi32(t1, _mm_cvtsi32_si128(1));	// shift out sign bit
		SIMD_ivec4 const t3 = _mm_set1_epi32(int(0xFF000000));				// exponent mask
		SIMD_ivec4 const t4 = _mm_and_si128(t2, t3);						// exponent
		SIMD_ivec4 const t5 = _mm_andnot_si128(t3, t2);					// fraction
		SIMD_ivec4 const Equal = _mm_cmpeq_epi32(t3, t4);
		SIMD_ivec4 const Nequal = _mm_cmpeq_epi32(t5, _mm_setzero_si128());
		SIMD_ivec4 const And = _mm_and_si128(Equal, Nequal);
		return _mm_castsi128_ps(And);									// exponent = all 1s and fraction != 0
	}
	
	// Agner Fog method
	SIMD_vec4 SIMD_vec4_inf(SIMD_vec4 x)
	{
		SIMD_ivec4 const t1 = _mm_castps_si128(x);										// reinterpret as 32-bit integer
		SIMD_ivec4 const t2 = _mm_sll_epi32(t1, _mm_cvtsi32_si128(1));					// shift out sign bit
		return _mm_castsi128_ps(_mm_cmpeq_epi32(t2, _mm_set1_epi32(int(0xFF000000))));		// exponent is all 1s, fraction is 0
	}
	
#endif//CPU_ARCH & CPU_ARCH_SSE2_BIT
}



















namespace simd
{
	template<typename T>
	DLL_APICALL pack_t<3, T>	DLL_CALLING_CONVENTION polar
	(
	 pack_t<3, T> const& euclidean
	 )
	{
		T const Length(length(euclidean));
		pack_t<3, T> const tmp(euclidean / Length);
		T const xz_dist(sqrt(tmp.x * tmp.x + tmp.z * tmp.z));
		
		return pack_t<3, T>(
							asin(tmp.y),	// latitude
							atan(tmp.x, tmp.z),		// longitude
							xz_dist);				// xz distance
	}
	
	template<typename T>
	DLL_APICALL pack_t<3, T>	DLL_CALLING_CONVENTION euclidean
	(
	 pack_t<2, T> const& polar
	 )
	{
		T const latitude(polar.x);
		T const longitude(polar.y);
		
		return pack_t<3, T>(
							cos(latitude) * sin(longitude),
							sin(latitude),
							cos(latitude) * cos(longitude));
	}
	
}//namespace glm













//
//namespace sml
//{
//	// sec
//	template<typename T>
//	GLM_FUNC_QUALIFIER T sec(T angle)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'sec' only accept floating-point values");
//		return T(1) / glm::cos(angle);
//	}
//	
//	template<length_t L, typename T>
//	GLM_FUNC_QUALIFIER pack_t<L, T> sec(pack_t<L, T> const& x)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'sec' only accept floating-point inputs");
//		return detail::functor1<L, T, T>::call(sec, x);
//	}
//	
//	// csc
//	template<typename T>
//	GLM_FUNC_QUALIFIER T csc(T angle)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'csc' only accept floating-point values");
//		return T(1) / glm::sin(angle);
//	}
//	
//	template<length_t L, typename T>
//	GLM_FUNC_QUALIFIER pack_t<L, T> csc(pack_t<L, T> const& x)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'csc' only accept floating-point inputs");
//		return detail::functor1<L, T, T>::call(csc, x);
//	}
//	
//	// cot
//	template<typename T>
//	GLM_FUNC_QUALIFIER T cot(T angle)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'cot' only accept floating-point values");
//		
//		T const pi_over_2 = T(3.1415926535897932384626433832795 / 2.0);
//		return glm::tan(pi_over_2 - angle);
//	}
//	
//	template<length_t L, typename T>
//	GLM_FUNC_QUALIFIER pack_t<L, T> cot(pack_t<L, T> const& x)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'cot' only accept floating-point inputs");
//		return detail::functor1<L, T, T>::call(cot, x);
//	}
//	
//	// asec
//	template<typename T>
//	GLM_FUNC_QUALIFIER T asec(T x)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'asec' only accept floating-point values");
//		return acos(T(1) / x);
//	}
//	
//	template<length_t L, typename T>
//	GLM_FUNC_QUALIFIER pack_t<L, T> asec(pack_t<L, T> const& x)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'asec' only accept floating-point inputs");
//		return detail::functor1<L, T, T>::call(asec, x);
//	}
//	
//	// acsc
//	template<typename T>
//	GLM_FUNC_QUALIFIER T acsc(T x)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'acsc' only accept floating-point values");
//		return asin(T(1) / x);
//	}
//	
//	template<length_t L, typename T>
//	GLM_FUNC_QUALIFIER pack_t<L, T> acsc(pack_t<L, T> const& x)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'acsc' only accept floating-point inputs");
//		return detail::functor1<L, T, T>::call(acsc, x);
//	}
//	
//	// acot
//	template<typename T>
//	GLM_FUNC_QUALIFIER T acot(T x)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'acot' only accept floating-point values");
//		
//		T const pi_over_2 = T(3.1415926535897932384626433832795 / 2.0);
//		return pi_over_2 - atan(x);
//	}
//	
//	template<length_t L, typename T>
//	GLM_FUNC_QUALIFIER pack_t<L, T> acot(pack_t<L, T> const& x)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'acot' only accept floating-point inputs");
//		return detail::functor1<L, T, T>::call(acot, x);
//	}
//	
//	// sech
//	template<typename T>
//	GLM_FUNC_QUALIFIER T sech(T angle)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'sech' only accept floating-point values");
//		return T(1) / glm::cosh(angle);
//	}
//	
//	template<length_t L, typename T>
//	GLM_FUNC_QUALIFIER pack_t<L, T> sech(pack_t<L, T> const& x)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'sech' only accept floating-point inputs");
//		return detail::functor1<L, T, T>::call(sech, x);
//	}
//	
//	// csch
//	template<typename T>
//	GLM_FUNC_QUALIFIER T csch(T angle)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'csch' only accept floating-point values");
//		return T(1) / glm::sinh(angle);
//	}
//	
//	template<length_t L, typename T>
//	GLM_FUNC_QUALIFIER pack_t<L, T> csch(pack_t<L, T> const& x)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'csch' only accept floating-point inputs");
//		return detail::functor1<L, T, T>::call(csch, x);
//	}
//	
//	// coth
//	template<typename T>
//	GLM_FUNC_QUALIFIER T coth(T angle)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'coth' only accept floating-point values");
//		return glm::cosh(angle) / glm::sinh(angle);
//	}
//	
//	template<length_t L, typename T>
//	GLM_FUNC_QUALIFIER pack_t<L, T> coth(pack_t<L, T> const& x)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'coth' only accept floating-point inputs");
//		return detail::functor1<L, T, T>::call(coth, x);
//	}
//	
//	// asech
//	template<typename T>
//	GLM_FUNC_QUALIFIER T asech(T x)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'asech' only accept floating-point values");
//		return acosh(T(1) / x);
//	}
//	
//	template<length_t L, typename T>
//	GLM_FUNC_QUALIFIER pack_t<L, T> asech(pack_t<L, T> const& x)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'asech' only accept floating-point inputs");
//		return detail::functor1<L, T, T>::call(asech, x);
//	}
//	
//	// acsch
//	template<typename T>
//	GLM_FUNC_QUALIFIER T acsch(T x)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'acsch' only accept floating-point values");
//		return asinh(T(1) / x);
//	}
//	
//	template<length_t L, typename T>
//	GLM_FUNC_QUALIFIER pack_t<L, T> acsch(pack_t<L, T> const& x)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'acsch' only accept floating-point inputs");
//		return detail::functor1<L, T, T>::call(acsch, x);
//	}
//	
//	// acoth
//	template<typename T>
//	GLM_FUNC_QUALIFIER T acoth(T x)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'acoth' only accept floating-point values");
//		return atanh(T(1) / x);
//	}
//	
//	template<length_t L, typename T>
//	GLM_FUNC_QUALIFIER pack_t<L, T> acoth(pack_t<L, T> const& x)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'acoth' only accept floating-point inputs");
//		return detail::functor1<L, T, T>::call(acoth, x);
//	}
//}//namespace glm



#if GLM_ARCH & GLM_ARCH_SSE2_BIT

namespace sml
{
	namespace detail
	{
//		template<bool N>
//		struct compute_length<4, float, N, true>
//		{
//			GLM_FUNC_QUALIFIER static float call(pack_t<4, float> const& v)
//			{
//				return _mm_cvtss_f32(sml_vec4_length(v.data));
//			}
//		};
//
//		template<bool N>
//		struct compute_distance<4, float, N, true>
//		{
//			GLM_FUNC_QUALIFIER static float call(pack_t<4, float> const& p0, pack_t<4, float> const& p1)
//			{
//				return _mm_cvtss_f32(sml_vec4_distance(p0.data, p1.data));
//			}
//		};
//
//		template<bool N>
//		struct compute_dot<pack_t<4, float>, float, true>
//		{
//		};
//
//
//		template<bool N>
//		struct compute_normalize<4, float, N, true>
//		{
//		};
//
//		template<bool N>
//		struct compute_faceforward<4, float, N, true>
//		{
//			GLM_FUNC_QUALIFIER static pack_t<4, float> call(pack_t<4, float> const& N, pack_t<4, float> const& I, pack_t<4, float> const& Nref)
//			{
//				pack_t<4, float> Result;
//				Result.data = sml_vec4_faceforward(N.data, I.data, Nref.data);
//				return Result;
//			}
//		};
//
//		template<bool N>
//		struct compute_reflect<4, float, N, true>
//		{
//			GLM_FUNC_QUALIFIER static pack_t<4, float> call(pack_t<4, float> const& I, pack_t<4, float> const& N)
//			{
//				pack_t<4, float> Result;
//				Result.data = sml_vec4_reflect(I.data, N.data);
//				return Result;
//			}
//		};
//
//		template<bool N>
//		struct compute_refract<4, float, N, true>
//		{
//			GLM_FUNC_QUALIFIER static pack_t<4, float> call(pack_t<4, float> const& I, pack_t<4, float> const& N, float eta)
//			{
//				pack_t<4, float> Result;
//				Result.data = sml_vec4_refract(I.data, N.data, _mm_set1_ps(eta));
//				return Result;
//			}
//		};
	}//namespace detail
}//namespace glm

#endif//GLM_ARCH & GLM_ARCH_SSE2_BIT

namespace sml
{
	namespace detail
	{
//		template<length_t L, typename T, bool N, bool Aligned>
//		struct compute_length
//		{
//			GLM_FUNC_QUALIFIER static T call(pack_t<L, T> const& v)
//			{
//				return sqrt(dot(v, v));
//			}
//		};
//
//		template<length_t L, typename T, bool N, bool Aligned>
//		struct compute_distance
//		{
//			GLM_FUNC_QUALIFIER static T call(pack_t<L, T> const& p0, pack_t<L, T> const& p1)
//			{
//				return length(p1 - p0);
//			}
//		};
//
//		template<typename V, typename T, bool Aligned>
//		struct compute_dot{};
//
//		template<typename T, bool N, bool Aligned>
//		struct compute_dot<pack_t<1, T>, T, Aligned>
//		{
//			GLM_FUNC_QUALIFIER static T call(pack_t<1, T> const& a, pack_t<1, T> const& b)
//			{
//				return a.x * b.x;
//			}
//		};
//
//		template<typename T, bool N, bool Aligned>
//		struct compute_dot<pack_t<2, T>, T, Aligned>
//		{
//			GLM_FUNC_QUALIFIER static T call(pack_t<2, T> const& a, pack_t<2, T> const& b)
//			{
//				pack_t<2, T> tmp(a * b);
//				return tmp.x + tmp.y;
//			}
//		};
//
//		template<typename T, bool N, bool Aligned>
//		struct compute_dot<pack_t<3, T>, T, Aligned>
//		{
//			GLM_FUNC_QUALIFIER static T call(pack_t<3, T> const& a, pack_t<3, T> const& b)
//			{
//				pack_t<3, T> tmp(a * b);
//				return tmp.x + tmp.y + tmp.z;
//			}
//		};
//
//		template<typename T, bool N, bool Aligned>
//		struct compute_dot<pack_t<4, T>, T, Aligned>
//		{
//			GLM_FUNC_QUALIFIER static T call(pack_t<4, T> const& a, pack_t<4, T> const& b)
//			{
//				pack_t<4, T> tmp(a * b);
//				return (tmp.x + tmp.y) + (tmp.z + tmp.w);
//			}
//		};
//
//
//		template<length_t L, typename T, bool N, bool Aligned>
//		struct compute_normalize
//		{
//			GLM_FUNC_QUALIFIER static pack_t<L, T> call(pack_t<L, T> const& v)
//			{
//				GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'normalize' accepts only floating-point inputs");
//			}
//		};
//
//		template<length_t L, typename T, bool N, bool Aligned>
//		struct compute_faceforward
//		{
//			GLM_FUNC_QUALIFIER static pack_t<L, T> call(pack_t<L, T> const& N, pack_t<L, T> const& I, pack_t<L, T> const& Nref)
//			{
//				GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'normalize' accepts only floating-point inputs");
//
//				return dot(Nref, I) < static_cast<T>(0) ? N : -N;
//			}
//		};
//
//		template<length_t L, typename T, bool N, bool Aligned>
//		struct compute_reflect
//		{
//			GLM_FUNC_QUALIFIER static pack_t<L, T> call(pack_t<L, T> const& I, pack_t<L, T> const& N)
//			{
//				return I - N * dot(N, I) * static_cast<T>(2);
//			}
//		};
//
//		template<length_t L, typename T, bool N, bool Aligned>
//		struct compute_refract
//		{
//			GLM_FUNC_QUALIFIER static pack_t<L, T> call(pack_t<L, T> const& I, pack_t<L, T> const& N, T eta)
//			{
//				T const dotValue(dot(N, I));
//				T const k(static_cast<T>(1) - eta * eta * (static_cast<T>(1) - dotValue * dotValue));
//				return (eta * I - (eta * dotValue + std::sqrt(k)) * N) * static_cast<T>(k >= static_cast<T>(0));
//			}
//		};
	}//namespace detail
//
//	// length
//	template<typename genType>
//	GLM_FUNC_QUALIFIER genType length(genType x)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'length' accepts only floating-point inputs");
//
//		return abs(x);
//	}
//
//	template<length_t L, typename T>
//	GLM_FUNC_QUALIFIER T length(pack_t<L, T> const& v)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'length' accepts only floating-point inputs");
//
//		return detail::compute_length<L, T, N, detail::is_aligned<Q>::value>::call(v);
//	}
//
//	// distance
//	template<typename genType>
//	GLM_FUNC_QUALIFIER genType distance(genType const& p0, genType const& p1)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'distance' accepts only floating-point inputs");
//
//		return length(p1 - p0);
//	}
//
//	template<length_t L, typename T>
//	GLM_FUNC_QUALIFIER T distance(pack_t<L, T> const& p0, pack_t<L, T> const& p1)
//	{
//		return detail::compute_distance<L, T, N, detail::is_aligned<Q>::value>::call(p0, p1);
//	}
//
//	// dot
//	template<typename T>
//	GLM_FUNC_QUALIFIER T dot(T x, T y)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'dot' accepts only floating-point inputs");
//		return x * y;
//	}
//
//	template<length_t L, typename T>
//	GLM_FUNC_QUALIFIER T dot(pack_t<L, T> const& x, pack_t<L, T> const& y)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'dot' accepts only floating-point inputs");
//		return detail::compute_dot<pack_t<L, T>, T, detail::is_aligned<Q>::value>::call(x, y);
//	}
//
	

//
//	// normalize
//	template<typename genType>
//	GLM_FUNC_QUALIFIER genType normalize(genType const& x)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'normalize' accepts only floating-point inputs");
//
//		return x < genType(0) ? genType(-1) : genType(1);
//	}
//
//	template<length_t L, typename T>
//	GLM_FUNC_QUALIFIER pack_t<L, T> normalize(pack_t<L, T> const& x)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'normalize' accepts only floating-point inputs");
//
//		return detail::compute_normalize<L, T, N, detail::is_aligned<Q>::value>::call(x);
//	}
//
//	// faceforward
//	template<typename genType>
//	GLM_FUNC_QUALIFIER genType faceforward(genType const& N, genType const& I, genType const& Nref)
//	{
//		return dot(Nref, I) < static_cast<genType>(0) ? N : -N;
//	}
//
//	template<length_t L, typename T>
//	GLM_FUNC_QUALIFIER pack_t<L, T> faceforward(pack_t<L, T> const& N, pack_t<L, T> const& I, pack_t<L, T> const& Nref)
//	{
//		return detail::compute_faceforward<L, T, N, detail::is_aligned<Q>::value>::call(N, I, Nref);
//	}
//
//	// reflect
//	template<typename genType>
//	GLM_FUNC_QUALIFIER genType reflect(genType const& I, genType const& N)
//	{
//		return I - N * dot(N, I) * genType(2);
//	}
//
//	template<length_t L, typename T>
//	GLM_FUNC_QUALIFIER pack_t<L, T> reflect(pack_t<L, T> const& I, pack_t<L, T> const& N)
//	{
//		return detail::compute_reflect<L, T, N, detail::is_aligned<Q>::value>::call(I, N);
//	}
//
//	// refract
//	template<typename genType>
//	GLM_FUNC_QUALIFIER genType refract(genType const& I, genType const& N, genType eta)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'refract' accepts only floating-point inputs");
//		genType const dotValue(dot(N, I));
//		genType const k(static_cast<genType>(1) - eta * eta * (static_cast<genType>(1) - dotValue * dotValue));
//		return (eta * I - (eta * dotValue + sqrt(k)) * N) * static_cast<genType>(k >= static_cast<genType>(0));
//	}
//
//	template<length_t L, typename T>
//	GLM_FUNC_QUALIFIER pack_t<L, T> refract(pack_t<L, T> const& I, pack_t<L, T> const& N, T eta)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'refract' accepts only floating-point inputs");
//		return detail::compute_refract<L, T, N, detail::is_aligned<Q>::value>::call(I, N, eta);
//	}

}//namespace glm








////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// dot
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace simd
{
	template <int D, typename T>
	SIMD_STATIC_INLINE T perform_dot(pack_t<D, T> const& x, pack_t<D, T> const& y)
	{
		T ret = T(0);
		for (int i = 0; i < D; i++)
			ret += x[i] * y[i];
		return ret;
	}

	template <int D, typename T>
	SIMD_STATIC_INLINE T perform_dot3(pack_t<D, T> const& x, pack_t<D, T> const& y)
	{
		T ret = T(0);
		for (int i = 0; i < pre::min<D, 3>(); i++)
			ret += x[i] * y[i];
		return ret;
	}

#if CPU_ARCH & CPU_ARCH_SSE2_BIT
	template <bool N>
	SIMD_STATIC_INLINE float perform_dot(pack_t<4, float> const& x, pack_t<4, float> const& y)
	{
		return _mm_cvtss_f32(SIMD_vec1_dot(x.m, y.m));
	}
#endif

	VECTOR_L_V_V_F(dot)
	VECTOR_L_V_V_F(dot3)
}






////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// length
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace simd
{
#if CPU_ARCH & CPU_ARCH_SSE2_BIT
	template <bool N>
	SIMD_STATIC_INLINE	float perform_length(pack_t<4, float> const& x)
	{
		return _mm_cvtss_f32(SIMD_vec4_length(x.m));
	}
#endif

	template <int D, typename T>
	SIMD_STATIC_INLINE typename lengthi_t<pack_t<D, T>>::type DLL_CALLING_CONVENTION perform_length(pack_t<D, T> const& x)
	{
		return ::sqrt(dot(x, x));
	}

	VECTOR_IL_V_F(length)
}






////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// cross
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace simd
{
#if CPU_ARCH & CPU_ARCH_SSE2_BIT
	
	SIMD_vec4 SIMD_vec4_cross(SIMD_vec4 v1, SIMD_vec4 v2)
	{
		SIMD_vec4 const swp0 = _mm_shuffle_ps(v1, v1, _MM_SHUFFLE(3, 0, 2, 1));
		SIMD_vec4 const swp1 = _mm_shuffle_ps(v1, v1, _MM_SHUFFLE(3, 1, 0, 2));
		SIMD_vec4 const swp2 = _mm_shuffle_ps(v2, v2, _MM_SHUFFLE(3, 0, 2, 1));
		SIMD_vec4 const swp3 = _mm_shuffle_ps(v2, v2, _MM_SHUFFLE(3, 1, 0, 2));
		SIMD_vec4 const mul0 = _mm_mul_ps(swp0, swp3);
		SIMD_vec4 const mul1 = _mm_mul_ps(swp1, swp2);
		SIMD_vec4 const sub0 = _mm_sub_ps(mul0, mul1);
		return sub0;
	}
	
	pack_t<3, float> SIMD_vec4_cross(pack_t<3, float> const& a, pack_t<3, float> const& b)
	{
		__m128 const set0 = _mm_set_ps(0.0f, a.z, a.y, a.x);
		__m128 const set1 = _mm_set_ps(0.0f, b.z, b.y, b.x);
		__m128 const xpd0 = SIMD_vec4_cross(set0, set1);
		
		pack_t<4, float> Result;
		Result.m = xpd0;
		return pack_t<3, float>(Result.x, Result.y, Result.z);
	}
	
#endif//CPU_ARCH & CPU_ARCH_SSE2_BIT
	
	SIMD_DLLFNC_T2(auto) cross(pack_t<2, T1> const& v, pack_t<2, T2> const& u) -> typename high_precission_t<T1, T2>::type
	{
		return v.x * u.y - u.x * v.y;
	}

	SIMD_DLLFNC_T2(auto)	cross(pack_t<3, T1> const& x, pack_t<3, T2> const& y) -> pack_t<3, typename high_precission_t<T1, T2>::type>
	{
		return pack_t<3, typename high_precission_t<T1, T2>::type>(
							  x.y * y.z - y.y * x.z,
							  x.z * y.x - y.z * x.x,
							  x.x * y.y - y.x * x.y);
	}
	
#define VEC_CROSS_T2(_T1, _T2) \
	template	DLL_FUNCTION(auto)	cross<_T1, _T2>(pack_t<2, _T1> const& v, pack_t<2, _T2> const& u) -> high_precission_t<_T1, _T2>::type; \
	template	DLL_FUNCTION(auto)	cross<_T1, _T2>(pack_t<3, _T1> const& v, pack_t<3, _T2> const& u) -> pack_t<3, high_precission_t<_T1, _T2>::type>; \

	VEC_CROSS_T2(double, double);
	VEC_CROSS_T2(double, float);
	VEC_CROSS_T2(float, double);
	VEC_CROSS_T2(float, float);
}









////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// rsq
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace simd
{
	template <int D, typename T>
	SIMD_STATIC_INLINE pack_t<D, T> perform_rsq(const pack_t<D, T>& v)
	{
		pack_t<D, T> ret;
		for (int i = 0; i < D; i++)
			ret[i] = T(1.0) / ::sqrt(ret[i]);
		return ret;
	}
	
	VECTOR_V_V_F(rsq)
	
}











////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// normalize
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace simd
{
	template <int D, typename T>
	SIMD_STATIC_INLINE pack_t<D, T> perform_normalize(const pack_t<D, T>& v)
	{
		return v * T(rsq(dot(v, v)));
	}

#if CPU_ARCH & CPU_ARCH_SSE2_BIT
	template<bool N>
	SIMD_STATIC_INLINE pack_t<4, float> perform_normalize(const pack_t<4, float>& v)
	{
		return SIMD_vec4_normalize(v.m);
	}
#endif
	
	VECTOR_V_V_F(normalize)

}





























////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// radians
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace simd
{
	template <int D, typename T>
	SIMD_STATIC_INLINE pack_t<D, T> DLL_CALLING_CONVENTION perform_radians(const pack_t<D, T>& x)
	{
		return x * pack_t<D, T>(static_cast<T>(0.01745329251994329576923690768489));
	}

	VECTOR_V_V_F(radians)
}















////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// operators
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace simd
{
#define VECTOR_PERFORM_OPERATOR(_O, _F) \
template<int D, typename T> \
FORCE_INLINE	pack_t<D, T>	DLL_CALLING_CONVENTION perform_operator_##_F(const pack_t<D, T>& a, const pack_t<D, T>& b) \
{ \
	pack_t<D, T> ret; \
	for (int i = 0; i < D; i++) \
		ret[i] = a[i] _O b[i]; \
	return ret; \
}\
template<int D, typename T> \
FORCE_INLINE	pack_t<D, T>	DLL_CALLING_CONVENTION perform_operator_##_F(const pack_t<D, T>& a, const T& b) \
{ \
	pack_t<D, T> ret; \
	for (int i = 0; i < D; i++) \
		ret[i] = a[i] _O b; \
	return ret; \
}\
template<int D, typename T> \
FORCE_INLINE	pack_t<D, T>	DLL_CALLING_CONVENTION perform_operator_##_F(const T& a, const pack_t<D, T>& b) \
{ \
	pack_t<D, T> ret; \
	for (int i = 0; i < D; i++) \
		ret[i] = a _O b[i]; \
	return ret; \
}\

	VECTOR_PERFORM_OPERATOR(+, add)
	VECTOR_PERFORM_OPERATOR(-, sub)
	VECTOR_PERFORM_OPERATOR(*, mul)
	VECTOR_PERFORM_OPERATOR(/, div)

	

	VECTOR_V_V_V_O_F(+, add)
	VECTOR_V_V_V_O_F(-, sub)
	VECTOR_V_V_V_O_F(*, mul)
	VECTOR_V_V_V_O_F(/, div)
}















namespace simd
{
#define VEC_OP_ASSIGN_D_T(_D, _T) \
template	DLL_FUNCTION(void)	pack_t<_D, _T>::operator += (const pack_t<_D, _T>&); \
template	DLL_FUNCTION(void)	pack_t<_D, _T>::operator += (const _T&); \
template	DLL_FUNCTION(void)	pack_t<_D, _T>::operator -= (const pack_t<_D, _T>&); \
template	DLL_FUNCTION(void)	pack_t<_D, _T>::operator -= (const _T&); \
template	DLL_FUNCTION(void)	pack_t<_D, _T>::operator *= (const pack_t<_D, _T>&); \
template	DLL_FUNCTION(void)	pack_t<_D, _T>::operator *= (const _T&); \
template	DLL_FUNCTION(void)	pack_t<_D, _T>::operator /= (const pack_t<_D, _T>&); \
template	DLL_FUNCTION(void)	pack_t<_D, _T>::operator /= (const _T&); \

#define VEC_OP_ASSIGN_T(_T) \
VEC_OP_ASSIGN_D_T(2, _T) \
VEC_OP_ASSIGN_D_T(3, _T) \
VEC_OP_ASSIGN_D_T(4, _T)

	VEC_OP_ASSIGN_T(int8_t)
	VEC_OP_ASSIGN_T(uint8_t)
	VEC_OP_ASSIGN_T(int16_t)
	VEC_OP_ASSIGN_T(uint16_t)
	VEC_OP_ASSIGN_T(int32_t)
	VEC_OP_ASSIGN_T(uint32_t)
	VEC_OP_ASSIGN_T(int64_t)
	VEC_OP_ASSIGN_T(uint64_t)
	VEC_OP_ASSIGN_T(float)
	VEC_OP_ASSIGN_T(double)


	
}
