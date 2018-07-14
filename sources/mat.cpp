// This software is MIT licensed (see LICENSE)

#define BUILD_DLL

#include <simd/core.h>
#include <limits>

#include <pre/number.h>

#include "defines.h"

namespace simd
{
#define DEFINE_EXTERN_MAT(_T) \
template struct mat_t<4, 4, _T>;
	
	DEFINE_EXTERN_MAT(float)
	DEFINE_EXTERN_MAT(double)

}

namespace simd
{
	// Matrix functions
#if CPU_ARCH & CPU_ARCH_SSE2_BIT
	
	SIMD_STATIC_INLINE void SIMD_mat4_matrixCompMult(SIMD_vec4 const in1[4], SIMD_vec4 const in2[4], SIMD_vec4 out[4])
	{
		out[0] = _mm_mul_ps(in1[0], in2[0]);
		out[1] = _mm_mul_ps(in1[1], in2[1]);
		out[2] = _mm_mul_ps(in1[2], in2[2]);
		out[3] = _mm_mul_ps(in1[3], in2[3]);
	}
	
	SIMD_STATIC_INLINE void SIMD_mat4_add(SIMD_vec4 const in1[4], SIMD_vec4 const in2[4], SIMD_vec4 out[4])
	{
		out[0] = _mm_add_ps(in1[0], in2[0]);
		out[1] = _mm_add_ps(in1[1], in2[1]);
		out[2] = _mm_add_ps(in1[2], in2[2]);
		out[3] = _mm_add_ps(in1[3], in2[3]);
	}
	
	SIMD_STATIC_INLINE void SIMD_mat4_sub(SIMD_vec4 const in1[4], SIMD_vec4 const in2[4], SIMD_vec4 out[4])
	{
		out[0] = _mm_sub_ps(in1[0], in2[0]);
		out[1] = _mm_sub_ps(in1[1], in2[1]);
		out[2] = _mm_sub_ps(in1[2], in2[2]);
		out[3] = _mm_sub_ps(in1[3], in2[3]);
	}
	
	SIMD_STATIC_INLINE SIMD_vec4 SIMD_mat4_mul_vec4(SIMD_vec4 const m[4], SIMD_vec4 v)
	{
		__m128 v0 = _mm_shuffle_ps(v, v, _MM_SHUFFLE(0, 0, 0, 0));
		__m128 v1 = _mm_shuffle_ps(v, v, _MM_SHUFFLE(1, 1, 1, 1));
		__m128 v2 = _mm_shuffle_ps(v, v, _MM_SHUFFLE(2, 2, 2, 2));
		__m128 v3 = _mm_shuffle_ps(v, v, _MM_SHUFFLE(3, 3, 3, 3));
		
		__m128 m0 = _mm_mul_ps(m[0], v0);
		__m128 m1 = _mm_mul_ps(m[1], v1);
		__m128 m2 = _mm_mul_ps(m[2], v2);
		__m128 m3 = _mm_mul_ps(m[3], v3);
		
		__m128 a0 = _mm_add_ps(m0, m1);
		__m128 a1 = _mm_add_ps(m2, m3);
		__m128 a2 = _mm_add_ps(a0, a1);
		
		return a2;
	}
	
	SIMD_STATIC_INLINE __m128 SIMD_vec4_mul_mat4(SIMD_vec4 v, SIMD_vec4 const m[4])
	{
		__m128 i0 = m[0];
		__m128 i1 = m[1];
		__m128 i2 = m[2];
		__m128 i3 = m[3];
		
		__m128 m0 = _mm_mul_ps(v, i0);
		__m128 m1 = _mm_mul_ps(v, i1);
		__m128 m2 = _mm_mul_ps(v, i2);
		__m128 m3 = _mm_mul_ps(v, i3);
		
		__m128 u0 = _mm_unpacklo_ps(m0, m1);
		__m128 u1 = _mm_unpackhi_ps(m0, m1);
		__m128 a0 = _mm_add_ps(u0, u1);
		
		__m128 u2 = _mm_unpacklo_ps(m2, m3);
		__m128 u3 = _mm_unpackhi_ps(m2, m3);
		__m128 a1 = _mm_add_ps(u2, u3);
		
		__m128 f0 = _mm_movelh_ps(a0, a1);
		__m128 f1 = _mm_movehl_ps(a1, a0);
		__m128 f2 = _mm_add_ps(f0, f1);
		
		return f2;
	}
	
	bool is_aligned(const volatile void *p, std::size_t n)
	{
		return reinterpret_cast<uintptr_t>(p) % n == 0;
	}
	
	SIMD_STATIC_INLINE void SIMD_mat4_mul(const SIMD_vec4 in1[4], const SIMD_vec4 in2[4], SIMD_vec4 out[4])
	{
		{
			__m128 e0 = _mm_shuffle_ps(in2[0], in2[0], _MM_SHUFFLE(0, 0, 0, 0));
			__m128 e1 = _mm_shuffle_ps(in2[0], in2[0], _MM_SHUFFLE(1, 1, 1, 1));
			__m128 e2 = _mm_shuffle_ps(in2[0], in2[0], _MM_SHUFFLE(2, 2, 2, 2));
			__m128 e3 = _mm_shuffle_ps(in2[0], in2[0], _MM_SHUFFLE(3, 3, 3, 3));
			
			__m128 m0 = _mm_mul_ps(in1[0], e0);
			__m128 m1 = _mm_mul_ps(in1[1], e1);
			__m128 m2 = _mm_mul_ps(in1[2], e2);
			__m128 m3 = _mm_mul_ps(in1[3], e3);
			
			__m128 a0 = _mm_add_ps(m0, m1);
			__m128 a1 = _mm_add_ps(m2, m3);
			__m128 a2 = _mm_add_ps(a0, a1);
			
			out[0] = a2;
		}
		
		{
			__m128 e0 = _mm_shuffle_ps(in2[1], in2[1], _MM_SHUFFLE(0, 0, 0, 0));
			__m128 e1 = _mm_shuffle_ps(in2[1], in2[1], _MM_SHUFFLE(1, 1, 1, 1));
			__m128 e2 = _mm_shuffle_ps(in2[1], in2[1], _MM_SHUFFLE(2, 2, 2, 2));
			__m128 e3 = _mm_shuffle_ps(in2[1], in2[1], _MM_SHUFFLE(3, 3, 3, 3));
			
			__m128 m0 = _mm_mul_ps(in1[0], e0);
			__m128 m1 = _mm_mul_ps(in1[1], e1);
			__m128 m2 = _mm_mul_ps(in1[2], e2);
			__m128 m3 = _mm_mul_ps(in1[3], e3);
			
			__m128 a0 = _mm_add_ps(m0, m1);
			__m128 a1 = _mm_add_ps(m2, m3);
			__m128 a2 = _mm_add_ps(a0, a1);
			
			out[1] = a2;
		}
		
		{
			__m128 e0 = _mm_shuffle_ps(in2[2], in2[2], _MM_SHUFFLE(0, 0, 0, 0));
			__m128 e1 = _mm_shuffle_ps(in2[2], in2[2], _MM_SHUFFLE(1, 1, 1, 1));
			__m128 e2 = _mm_shuffle_ps(in2[2], in2[2], _MM_SHUFFLE(2, 2, 2, 2));
			__m128 e3 = _mm_shuffle_ps(in2[2], in2[2], _MM_SHUFFLE(3, 3, 3, 3));
			
			__m128 m0 = _mm_mul_ps(in1[0], e0);
			__m128 m1 = _mm_mul_ps(in1[1], e1);
			__m128 m2 = _mm_mul_ps(in1[2], e2);
			__m128 m3 = _mm_mul_ps(in1[3], e3);
			
			__m128 a0 = _mm_add_ps(m0, m1);
			__m128 a1 = _mm_add_ps(m2, m3);
			__m128 a2 = _mm_add_ps(a0, a1);
			
			out[2] = a2;
		}
		
		{
			//(__m128&)_mm_shuffle_epi32(__m128i&)in2[0], _MM_SHUFFLE(3, 3, 3, 3))
			__m128 e0 = _mm_shuffle_ps(in2[3], in2[3], _MM_SHUFFLE(0, 0, 0, 0));
			__m128 e1 = _mm_shuffle_ps(in2[3], in2[3], _MM_SHUFFLE(1, 1, 1, 1));
			__m128 e2 = _mm_shuffle_ps(in2[3], in2[3], _MM_SHUFFLE(2, 2, 2, 2));
			__m128 e3 = _mm_shuffle_ps(in2[3], in2[3], _MM_SHUFFLE(3, 3, 3, 3));
			
			__m128 m0 = _mm_mul_ps(in1[0], e0);
			__m128 m1 = _mm_mul_ps(in1[1], e1);
			__m128 m2 = _mm_mul_ps(in1[2], e2);
			__m128 m3 = _mm_mul_ps(in1[3], e3);
			
			__m128 a0 = _mm_add_ps(m0, m1);
			__m128 a1 = _mm_add_ps(m2, m3);
			__m128 a2 = _mm_add_ps(a0, a1);
			
			out[3] = a2;
		}
	}
	
	SIMD_STATIC_INLINE void SIMD_mat4_transpose(SIMD_vec4 const in[4], SIMD_vec4 out[4])
	{
		__m128 tmp0 = _mm_shuffle_ps(in[0], in[1], 0x44);
		__m128 tmp2 = _mm_shuffle_ps(in[0], in[1], 0xEE);
		__m128 tmp1 = _mm_shuffle_ps(in[2], in[3], 0x44);
		__m128 tmp3 = _mm_shuffle_ps(in[2], in[3], 0xEE);
		
		out[0] = _mm_shuffle_ps(tmp0, tmp1, 0x88);
		out[1] = _mm_shuffle_ps(tmp0, tmp1, 0xDD);
		out[2] = _mm_shuffle_ps(tmp2, tmp3, 0x88);
		out[3] = _mm_shuffle_ps(tmp2, tmp3, 0xDD);
	}
	
	SIMD_STATIC_INLINE SIMD_vec4 SIMD_mat4_determinant_highp(SIMD_vec4 const in[4])
	{
		__m128 Fac0;
		{
			//	valType SubFactor00 = m[2][2] * m[3][3] - m[3][2] * m[2][3];
			//	valType SubFactor00 = m[2][2] * m[3][3] - m[3][2] * m[2][3];
			//	valType SubFactor06 = m[1][2] * m[3][3] - m[3][2] * m[1][3];
			//	valType SubFactor13 = m[1][2] * m[2][3] - m[2][2] * m[1][3];
			
			__m128 Swp0a = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(3, 3, 3, 3));
			__m128 Swp0b = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(2, 2, 2, 2));
			
			__m128 Swp00 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(2, 2, 2, 2));
			__m128 Swp01 = _mm_shuffle_ps(Swp0a, Swp0a, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp02 = _mm_shuffle_ps(Swp0b, Swp0b, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp03 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(3, 3, 3, 3));
			
			__m128 Mul00 = _mm_mul_ps(Swp00, Swp01);
			__m128 Mul01 = _mm_mul_ps(Swp02, Swp03);
			Fac0 = _mm_sub_ps(Mul00, Mul01);
		}
		
		__m128 Fac1;
		{
			//	valType SubFactor01 = m[2][1] * m[3][3] - m[3][1] * m[2][3];
			//	valType SubFactor01 = m[2][1] * m[3][3] - m[3][1] * m[2][3];
			//	valType SubFactor07 = m[1][1] * m[3][3] - m[3][1] * m[1][3];
			//	valType SubFactor14 = m[1][1] * m[2][3] - m[2][1] * m[1][3];
			
			__m128 Swp0a = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(3, 3, 3, 3));
			__m128 Swp0b = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(1, 1, 1, 1));
			
			__m128 Swp00 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(1, 1, 1, 1));
			__m128 Swp01 = _mm_shuffle_ps(Swp0a, Swp0a, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp02 = _mm_shuffle_ps(Swp0b, Swp0b, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp03 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(3, 3, 3, 3));
			
			__m128 Mul00 = _mm_mul_ps(Swp00, Swp01);
			__m128 Mul01 = _mm_mul_ps(Swp02, Swp03);
			Fac1 = _mm_sub_ps(Mul00, Mul01);
		}
		
		
		__m128 Fac2;
		{
			//	valType SubFactor02 = m[2][1] * m[3][2] - m[3][1] * m[2][2];
			//	valType SubFactor02 = m[2][1] * m[3][2] - m[3][1] * m[2][2];
			//	valType SubFactor08 = m[1][1] * m[3][2] - m[3][1] * m[1][2];
			//	valType SubFactor15 = m[1][1] * m[2][2] - m[2][1] * m[1][2];
			
			__m128 Swp0a = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(2, 2, 2, 2));
			__m128 Swp0b = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(1, 1, 1, 1));
			
			__m128 Swp00 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(1, 1, 1, 1));
			__m128 Swp01 = _mm_shuffle_ps(Swp0a, Swp0a, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp02 = _mm_shuffle_ps(Swp0b, Swp0b, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp03 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(2, 2, 2, 2));
			
			__m128 Mul00 = _mm_mul_ps(Swp00, Swp01);
			__m128 Mul01 = _mm_mul_ps(Swp02, Swp03);
			Fac2 = _mm_sub_ps(Mul00, Mul01);
		}
		
		__m128 Fac3;
		{
			//	valType SubFactor03 = m[2][0] * m[3][3] - m[3][0] * m[2][3];
			//	valType SubFactor03 = m[2][0] * m[3][3] - m[3][0] * m[2][3];
			//	valType SubFactor09 = m[1][0] * m[3][3] - m[3][0] * m[1][3];
			//	valType SubFactor16 = m[1][0] * m[2][3] - m[2][0] * m[1][3];
			
			__m128 Swp0a = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(3, 3, 3, 3));
			__m128 Swp0b = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(0, 0, 0, 0));
			
			__m128 Swp00 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(0, 0, 0, 0));
			__m128 Swp01 = _mm_shuffle_ps(Swp0a, Swp0a, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp02 = _mm_shuffle_ps(Swp0b, Swp0b, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp03 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(3, 3, 3, 3));
			
			__m128 Mul00 = _mm_mul_ps(Swp00, Swp01);
			__m128 Mul01 = _mm_mul_ps(Swp02, Swp03);
			Fac3 = _mm_sub_ps(Mul00, Mul01);
		}
		
		__m128 Fac4;
		{
			//	valType SubFactor04 = m[2][0] * m[3][2] - m[3][0] * m[2][2];
			//	valType SubFactor04 = m[2][0] * m[3][2] - m[3][0] * m[2][2];
			//	valType SubFactor10 = m[1][0] * m[3][2] - m[3][0] * m[1][2];
			//	valType SubFactor17 = m[1][0] * m[2][2] - m[2][0] * m[1][2];
			
			__m128 Swp0a = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(2, 2, 2, 2));
			__m128 Swp0b = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(0, 0, 0, 0));
			
			__m128 Swp00 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(0, 0, 0, 0));
			__m128 Swp01 = _mm_shuffle_ps(Swp0a, Swp0a, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp02 = _mm_shuffle_ps(Swp0b, Swp0b, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp03 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(2, 2, 2, 2));
			
			__m128 Mul00 = _mm_mul_ps(Swp00, Swp01);
			__m128 Mul01 = _mm_mul_ps(Swp02, Swp03);
			Fac4 = _mm_sub_ps(Mul00, Mul01);
		}
		
		__m128 Fac5;
		{
			//	valType SubFactor05 = m[2][0] * m[3][1] - m[3][0] * m[2][1];
			//	valType SubFactor05 = m[2][0] * m[3][1] - m[3][0] * m[2][1];
			//	valType SubFactor12 = m[1][0] * m[3][1] - m[3][0] * m[1][1];
			//	valType SubFactor18 = m[1][0] * m[2][1] - m[2][0] * m[1][1];
			
			__m128 Swp0a = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(1, 1, 1, 1));
			__m128 Swp0b = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(0, 0, 0, 0));
			
			__m128 Swp00 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(0, 0, 0, 0));
			__m128 Swp01 = _mm_shuffle_ps(Swp0a, Swp0a, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp02 = _mm_shuffle_ps(Swp0b, Swp0b, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp03 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(1, 1, 1, 1));
			
			__m128 Mul00 = _mm_mul_ps(Swp00, Swp01);
			__m128 Mul01 = _mm_mul_ps(Swp02, Swp03);
			Fac5 = _mm_sub_ps(Mul00, Mul01);
		}
		
		__m128 SignA = _mm_set_ps( 1.0f,-1.0f, 1.0f,-1.0f);
		__m128 SignB = _mm_set_ps(-1.0f, 1.0f,-1.0f, 1.0f);
		
		// m[1][0]
		// m[0][0]
		// m[0][0]
		// m[0][0]
		__m128 Temp0 = _mm_shuffle_ps(in[1], in[0], _MM_SHUFFLE(0, 0, 0, 0));
		__m128 Vec0 = _mm_shuffle_ps(Temp0, Temp0, _MM_SHUFFLE(2, 2, 2, 0));
		
		// m[1][1]
		// m[0][1]
		// m[0][1]
		// m[0][1]
		__m128 Temp1 = _mm_shuffle_ps(in[1], in[0], _MM_SHUFFLE(1, 1, 1, 1));
		__m128 Vec1 = _mm_shuffle_ps(Temp1, Temp1, _MM_SHUFFLE(2, 2, 2, 0));
		
		// m[1][2]
		// m[0][2]
		// m[0][2]
		// m[0][2]
		__m128 Temp2 = _mm_shuffle_ps(in[1], in[0], _MM_SHUFFLE(2, 2, 2, 2));
		__m128 Vec2 = _mm_shuffle_ps(Temp2, Temp2, _MM_SHUFFLE(2, 2, 2, 0));
		
		// m[1][3]
		// m[0][3]
		// m[0][3]
		// m[0][3]
		__m128 Temp3 = _mm_shuffle_ps(in[1], in[0], _MM_SHUFFLE(3, 3, 3, 3));
		__m128 Vec3 = _mm_shuffle_ps(Temp3, Temp3, _MM_SHUFFLE(2, 2, 2, 0));
		
		// col0
		// + (Vec1[0] * Fac0[0] - Vec2[0] * Fac1[0] + Vec3[0] * Fac2[0]),
		// - (Vec1[1] * Fac0[1] - Vec2[1] * Fac1[1] + Vec3[1] * Fac2[1]),
		// + (Vec1[2] * Fac0[2] - Vec2[2] * Fac1[2] + Vec3[2] * Fac2[2]),
		// - (Vec1[3] * Fac0[3] - Vec2[3] * Fac1[3] + Vec3[3] * Fac2[3]),
		__m128 Mul00 = _mm_mul_ps(Vec1, Fac0);
		__m128 Mul01 = _mm_mul_ps(Vec2, Fac1);
		__m128 Mul02 = _mm_mul_ps(Vec3, Fac2);
		__m128 Sub00 = _mm_sub_ps(Mul00, Mul01);
		__m128 Add00 = _mm_add_ps(Sub00, Mul02);
		__m128 Inv0 = _mm_mul_ps(SignB, Add00);
		
		// col1
		// - (Vec0[0] * Fac0[0] - Vec2[0] * Fac3[0] + Vec3[0] * Fac4[0]),
		// + (Vec0[0] * Fac0[1] - Vec2[1] * Fac3[1] + Vec3[1] * Fac4[1]),
		// - (Vec0[0] * Fac0[2] - Vec2[2] * Fac3[2] + Vec3[2] * Fac4[2]),
		// + (Vec0[0] * Fac0[3] - Vec2[3] * Fac3[3] + Vec3[3] * Fac4[3]),
		__m128 Mul03 = _mm_mul_ps(Vec0, Fac0);
		__m128 Mul04 = _mm_mul_ps(Vec2, Fac3);
		__m128 Mul05 = _mm_mul_ps(Vec3, Fac4);
		__m128 Sub01 = _mm_sub_ps(Mul03, Mul04);
		__m128 Add01 = _mm_add_ps(Sub01, Mul05);
		__m128 Inv1 = _mm_mul_ps(SignA, Add01);
		
		// col2
		// + (Vec0[0] * Fac1[0] - Vec1[0] * Fac3[0] + Vec3[0] * Fac5[0]),
		// - (Vec0[0] * Fac1[1] - Vec1[1] * Fac3[1] + Vec3[1] * Fac5[1]),
		// + (Vec0[0] * Fac1[2] - Vec1[2] * Fac3[2] + Vec3[2] * Fac5[2]),
		// - (Vec0[0] * Fac1[3] - Vec1[3] * Fac3[3] + Vec3[3] * Fac5[3]),
		__m128 Mul06 = _mm_mul_ps(Vec0, Fac1);
		__m128 Mul07 = _mm_mul_ps(Vec1, Fac3);
		__m128 Mul08 = _mm_mul_ps(Vec3, Fac5);
		__m128 Sub02 = _mm_sub_ps(Mul06, Mul07);
		__m128 Add02 = _mm_add_ps(Sub02, Mul08);
		__m128 Inv2 = _mm_mul_ps(SignB, Add02);
		
		// col3
		// - (Vec1[0] * Fac2[0] - Vec1[0] * Fac4[0] + Vec2[0] * Fac5[0]),
		// + (Vec1[0] * Fac2[1] - Vec1[1] * Fac4[1] + Vec2[1] * Fac5[1]),
		// - (Vec1[0] * Fac2[2] - Vec1[2] * Fac4[2] + Vec2[2] * Fac5[2]),
		// + (Vec1[0] * Fac2[3] - Vec1[3] * Fac4[3] + Vec2[3] * Fac5[3]));
		__m128 Mul09 = _mm_mul_ps(Vec0, Fac2);
		__m128 Mul10 = _mm_mul_ps(Vec1, Fac4);
		__m128 Mul11 = _mm_mul_ps(Vec2, Fac5);
		__m128 Sub03 = _mm_sub_ps(Mul09, Mul10);
		__m128 Add03 = _mm_add_ps(Sub03, Mul11);
		__m128 Inv3 = _mm_mul_ps(SignA, Add03);
		
		__m128 Row0 = _mm_shuffle_ps(Inv0, Inv1, _MM_SHUFFLE(0, 0, 0, 0));
		__m128 Row1 = _mm_shuffle_ps(Inv2, Inv3, _MM_SHUFFLE(0, 0, 0, 0));
		__m128 Row2 = _mm_shuffle_ps(Row0, Row1, _MM_SHUFFLE(2, 0, 2, 0));
		
		//	valType Determinant = m[0][0] * Inverse[0][0]
		//						+ m[0][1] * Inverse[1][0]
		//						+ m[0][2] * Inverse[2][0]
		//						+ m[0][3] * Inverse[3][0];
		__m128 Det0 = SIMD_vec4_dot(in[0], Row2);
		return Det0;
	}
	
	FORCE_INLINE SIMD_vec4 SIMD_mat4_determinant_lowp(SIMD_vec4 const m[4])
	{
		// _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(
		
		//T SubFactor00 = m[2][2] * m[3][3] - m[3][2] * m[2][3];
		//T SubFactor01 = m[2][1] * m[3][3] - m[3][1] * m[2][3];
		//T SubFactor02 = m[2][1] * m[3][2] - m[3][1] * m[2][2];
		//T SubFactor03 = m[2][0] * m[3][3] - m[3][0] * m[2][3];
		//T SubFactor04 = m[2][0] * m[3][2] - m[3][0] * m[2][2];
		//T SubFactor05 = m[2][0] * m[3][1] - m[3][0] * m[2][1];
		
		// First 2 columns
		__m128 Swp2A = _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(m[2]), _MM_SHUFFLE(0, 1, 1, 2)));
		__m128 Swp3A = _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(m[3]), _MM_SHUFFLE(3, 2, 3, 3)));
		__m128 MulA = _mm_mul_ps(Swp2A, Swp3A);
		
		// Second 2 columns
		__m128 Swp2B = _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(m[2]), _MM_SHUFFLE(3, 2, 3, 3)));
		__m128 Swp3B = _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(m[3]), _MM_SHUFFLE(0, 1, 1, 2)));
		__m128 MulB = _mm_mul_ps(Swp2B, Swp3B);
		
		// Columns subtraction
		__m128 SubE = _mm_sub_ps(MulA, MulB);
		
		// Last 2 rows
		__m128 Swp2C = _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(m[2]), _MM_SHUFFLE(0, 0, 1, 2)));
		__m128 Swp3C = _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(m[3]), _MM_SHUFFLE(1, 2, 0, 0)));
		__m128 MulC = _mm_mul_ps(Swp2C, Swp3C);
		__m128 SubF = _mm_sub_ps(_mm_movehl_ps(MulC, MulC), MulC);
		
		//tvec4<T, P> DetCof(
		//	+ (m[1][1] * SubFactor00 - m[1][2] * SubFactor01 + m[1][3] * SubFactor02),
		//	- (m[1][0] * SubFactor00 - m[1][2] * SubFactor03 + m[1][3] * SubFactor04),
		//	+ (m[1][0] * SubFactor01 - m[1][1] * SubFactor03 + m[1][3] * SubFactor05),
		//	- (m[1][0] * SubFactor02 - m[1][1] * SubFactor04 + m[1][2] * SubFactor05));
		
		__m128 SubFacA = _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(SubE), _MM_SHUFFLE(2, 1, 0, 0)));
		__m128 SwpFacA = _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(m[1]), _MM_SHUFFLE(0, 0, 0, 1)));
		__m128 MulFacA = _mm_mul_ps(SwpFacA, SubFacA);
		
		__m128 SubTmpB = _mm_shuffle_ps(SubE, SubF, _MM_SHUFFLE(0, 0, 3, 1));
		__m128 SubFacB = _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(SubTmpB), _MM_SHUFFLE(3, 1, 1, 0)));//SubF[0], SubE[3], SubE[3], SubE[1];
		__m128 SwpFacB = _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(m[1]), _MM_SHUFFLE(1, 1, 2, 2)));
		__m128 MulFacB = _mm_mul_ps(SwpFacB, SubFacB);
		
		__m128 SubRes = _mm_sub_ps(MulFacA, MulFacB);
		
		__m128 SubTmpC = _mm_shuffle_ps(SubE, SubF, _MM_SHUFFLE(1, 0, 2, 2));
		__m128 SubFacC = _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(SubTmpC), _MM_SHUFFLE(3, 3, 2, 0)));
		__m128 SwpFacC = _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(m[1]), _MM_SHUFFLE(2, 3, 3, 3)));
		__m128 MulFacC = _mm_mul_ps(SwpFacC, SubFacC);
		
		__m128 AddRes = _mm_add_ps(SubRes, MulFacC);
		__m128 DetCof = _mm_mul_ps(AddRes, _mm_setr_ps( 1.0f,-1.0f, 1.0f,-1.0f));
		
		//return m[0][0] * DetCof[0]
		//	 + m[0][1] * DetCof[1]
		//	 + m[0][2] * DetCof[2]
		//	 + m[0][3] * DetCof[3];
		
		return SIMD_vec4_dot(m[0], DetCof);
	}
	
	SIMD_STATIC_INLINE SIMD_vec4 SIMD_mat4_determinant(SIMD_vec4 const m[4])
	{
		// _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(add)
		
		//T SubFactor00 = m[2][2] * m[3][3] - m[3][2] * m[2][3];
		//T SubFactor01 = m[2][1] * m[3][3] - m[3][1] * m[2][3];
		//T SubFactor02 = m[2][1] * m[3][2] - m[3][1] * m[2][2];
		//T SubFactor03 = m[2][0] * m[3][3] - m[3][0] * m[2][3];
		//T SubFactor04 = m[2][0] * m[3][2] - m[3][0] * m[2][2];
		//T SubFactor05 = m[2][0] * m[3][1] - m[3][0] * m[2][1];
		
		// First 2 columns
		__m128 Swp2A = _mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 1, 1, 2));
		__m128 Swp3A = _mm_shuffle_ps(m[3], m[3], _MM_SHUFFLE(3, 2, 3, 3));
		__m128 MulA = _mm_mul_ps(Swp2A, Swp3A);
		
		// Second 2 columns
		__m128 Swp2B = _mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(3, 2, 3, 3));
		__m128 Swp3B = _mm_shuffle_ps(m[3], m[3], _MM_SHUFFLE(0, 1, 1, 2));
		__m128 MulB = _mm_mul_ps(Swp2B, Swp3B);
		
		// Columns subtraction
		__m128 SubE = _mm_sub_ps(MulA, MulB);
		
		// Last 2 rows
		__m128 Swp2C = _mm_shuffle_ps(m[2], m[2], _MM_SHUFFLE(0, 0, 1, 2));
		__m128 Swp3C = _mm_shuffle_ps(m[3], m[3], _MM_SHUFFLE(1, 2, 0, 0));
		__m128 MulC = _mm_mul_ps(Swp2C, Swp3C);
		__m128 SubF = _mm_sub_ps(_mm_movehl_ps(MulC, MulC), MulC);
		
		//tvec4<T, P> DetCof(
		//	+ (m[1][1] * SubFactor00 - m[1][2] * SubFactor01 + m[1][3] * SubFactor02),
		//	- (m[1][0] * SubFactor00 - m[1][2] * SubFactor03 + m[1][3] * SubFactor04),
		//	+ (m[1][0] * SubFactor01 - m[1][1] * SubFactor03 + m[1][3] * SubFactor05),
		//	- (m[1][0] * SubFactor02 - m[1][1] * SubFactor04 + m[1][2] * SubFactor05));
		
		__m128 SubFacA = _mm_shuffle_ps(SubE, SubE, _MM_SHUFFLE(2, 1, 0, 0));
		__m128 SwpFacA = _mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(0, 0, 0, 1));
		__m128 MulFacA = _mm_mul_ps(SwpFacA, SubFacA);
		
		__m128 SubTmpB = _mm_shuffle_ps(SubE, SubF, _MM_SHUFFLE(0, 0, 3, 1));
		__m128 SubFacB = _mm_shuffle_ps(SubTmpB, SubTmpB, _MM_SHUFFLE(3, 1, 1, 0));//SubF[0], SubE[3], SubE[3], SubE[1];
		__m128 SwpFacB = _mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(1, 1, 2, 2));
		__m128 MulFacB = _mm_mul_ps(SwpFacB, SubFacB);
		
		__m128 SubRes = _mm_sub_ps(MulFacA, MulFacB);
		
		__m128 SubTmpC = _mm_shuffle_ps(SubE, SubF, _MM_SHUFFLE(1, 0, 2, 2));
		__m128 SubFacC = _mm_shuffle_ps(SubTmpC, SubTmpC, _MM_SHUFFLE(3, 3, 2, 0));
		__m128 SwpFacC = _mm_shuffle_ps(m[1], m[1], _MM_SHUFFLE(2, 3, 3, 3));
		__m128 MulFacC = _mm_mul_ps(SwpFacC, SubFacC);
		
		__m128 AddRes = _mm_add_ps(SubRes, MulFacC);
		__m128 DetCof = _mm_mul_ps(AddRes, _mm_setr_ps( 1.0f,-1.0f, 1.0f,-1.0f));
		
		//return m[0][0] * DetCof[0]
		//	 + m[0][1] * DetCof[1]
		//	 + m[0][2] * DetCof[2]
		//	 + m[0][3] * DetCof[3];
		
		return SIMD_vec4_dot(m[0], DetCof);
	}
	
	SIMD_STATIC_INLINE void SIMD_mat4_inverse(SIMD_vec4 const in[4], SIMD_vec4 out[4])
	{
		__m128 Fac0;
		{
			//	valType SubFactor00 = m[2][2] * m[3][3] - m[3][2] * m[2][3];
			//	valType SubFactor00 = m[2][2] * m[3][3] - m[3][2] * m[2][3];
			//	valType SubFactor06 = m[1][2] * m[3][3] - m[3][2] * m[1][3];
			//	valType SubFactor13 = m[1][2] * m[2][3] - m[2][2] * m[1][3];
			
			__m128 Swp0a = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(3, 3, 3, 3));
			__m128 Swp0b = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(2, 2, 2, 2));
			
			__m128 Swp00 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(2, 2, 2, 2));
			__m128 Swp01 = _mm_shuffle_ps(Swp0a, Swp0a, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp02 = _mm_shuffle_ps(Swp0b, Swp0b, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp03 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(3, 3, 3, 3));
			
			__m128 Mul00 = _mm_mul_ps(Swp00, Swp01);
			__m128 Mul01 = _mm_mul_ps(Swp02, Swp03);
			Fac0 = _mm_sub_ps(Mul00, Mul01);
		}
		
		__m128 Fac1;
		{
			//	valType SubFactor01 = m[2][1] * m[3][3] - m[3][1] * m[2][3];
			//	valType SubFactor01 = m[2][1] * m[3][3] - m[3][1] * m[2][3];
			//	valType SubFactor07 = m[1][1] * m[3][3] - m[3][1] * m[1][3];
			//	valType SubFactor14 = m[1][1] * m[2][3] - m[2][1] * m[1][3];
			
			__m128 Swp0a = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(3, 3, 3, 3));
			__m128 Swp0b = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(1, 1, 1, 1));
			
			__m128 Swp00 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(1, 1, 1, 1));
			__m128 Swp01 = _mm_shuffle_ps(Swp0a, Swp0a, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp02 = _mm_shuffle_ps(Swp0b, Swp0b, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp03 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(3, 3, 3, 3));
			
			__m128 Mul00 = _mm_mul_ps(Swp00, Swp01);
			__m128 Mul01 = _mm_mul_ps(Swp02, Swp03);
			Fac1 = _mm_sub_ps(Mul00, Mul01);
		}
		
		
		__m128 Fac2;
		{
			//	valType SubFactor02 = m[2][1] * m[3][2] - m[3][1] * m[2][2];
			//	valType SubFactor02 = m[2][1] * m[3][2] - m[3][1] * m[2][2];
			//	valType SubFactor08 = m[1][1] * m[3][2] - m[3][1] * m[1][2];
			//	valType SubFactor15 = m[1][1] * m[2][2] - m[2][1] * m[1][2];
			
			__m128 Swp0a = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(2, 2, 2, 2));
			__m128 Swp0b = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(1, 1, 1, 1));
			
			__m128 Swp00 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(1, 1, 1, 1));
			__m128 Swp01 = _mm_shuffle_ps(Swp0a, Swp0a, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp02 = _mm_shuffle_ps(Swp0b, Swp0b, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp03 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(2, 2, 2, 2));
			
			__m128 Mul00 = _mm_mul_ps(Swp00, Swp01);
			__m128 Mul01 = _mm_mul_ps(Swp02, Swp03);
			Fac2 = _mm_sub_ps(Mul00, Mul01);
		}
		
		__m128 Fac3;
		{
			//	valType SubFactor03 = m[2][0] * m[3][3] - m[3][0] * m[2][3];
			//	valType SubFactor03 = m[2][0] * m[3][3] - m[3][0] * m[2][3];
			//	valType SubFactor09 = m[1][0] * m[3][3] - m[3][0] * m[1][3];
			//	valType SubFactor16 = m[1][0] * m[2][3] - m[2][0] * m[1][3];
			
			__m128 Swp0a = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(3, 3, 3, 3));
			__m128 Swp0b = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(0, 0, 0, 0));
			
			__m128 Swp00 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(0, 0, 0, 0));
			__m128 Swp01 = _mm_shuffle_ps(Swp0a, Swp0a, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp02 = _mm_shuffle_ps(Swp0b, Swp0b, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp03 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(3, 3, 3, 3));
			
			__m128 Mul00 = _mm_mul_ps(Swp00, Swp01);
			__m128 Mul01 = _mm_mul_ps(Swp02, Swp03);
			Fac3 = _mm_sub_ps(Mul00, Mul01);
		}
		
		__m128 Fac4;
		{
			//	valType SubFactor04 = m[2][0] * m[3][2] - m[3][0] * m[2][2];
			//	valType SubFactor04 = m[2][0] * m[3][2] - m[3][0] * m[2][2];
			//	valType SubFactor10 = m[1][0] * m[3][2] - m[3][0] * m[1][2];
			//	valType SubFactor17 = m[1][0] * m[2][2] - m[2][0] * m[1][2];
			
			__m128 Swp0a = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(2, 2, 2, 2));
			__m128 Swp0b = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(0, 0, 0, 0));
			
			__m128 Swp00 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(0, 0, 0, 0));
			__m128 Swp01 = _mm_shuffle_ps(Swp0a, Swp0a, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp02 = _mm_shuffle_ps(Swp0b, Swp0b, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp03 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(2, 2, 2, 2));
			
			__m128 Mul00 = _mm_mul_ps(Swp00, Swp01);
			__m128 Mul01 = _mm_mul_ps(Swp02, Swp03);
			Fac4 = _mm_sub_ps(Mul00, Mul01);
		}
		
		__m128 Fac5;
		{
			//	valType SubFactor05 = m[2][0] * m[3][1] - m[3][0] * m[2][1];
			//	valType SubFactor05 = m[2][0] * m[3][1] - m[3][0] * m[2][1];
			//	valType SubFactor12 = m[1][0] * m[3][1] - m[3][0] * m[1][1];
			//	valType SubFactor18 = m[1][0] * m[2][1] - m[2][0] * m[1][1];
			
			__m128 Swp0a = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(1, 1, 1, 1));
			__m128 Swp0b = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(0, 0, 0, 0));
			
			__m128 Swp00 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(0, 0, 0, 0));
			__m128 Swp01 = _mm_shuffle_ps(Swp0a, Swp0a, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp02 = _mm_shuffle_ps(Swp0b, Swp0b, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp03 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(1, 1, 1, 1));
			
			__m128 Mul00 = _mm_mul_ps(Swp00, Swp01);
			__m128 Mul01 = _mm_mul_ps(Swp02, Swp03);
			Fac5 = _mm_sub_ps(Mul00, Mul01);
		}
		
		__m128 SignA = _mm_set_ps( 1.0f,-1.0f, 1.0f,-1.0f);
		__m128 SignB = _mm_set_ps(-1.0f, 1.0f,-1.0f, 1.0f);
		
		// m[1][0]
		// m[0][0]
		// m[0][0]
		// m[0][0]
		__m128 Temp0 = _mm_shuffle_ps(in[1], in[0], _MM_SHUFFLE(0, 0, 0, 0));
		__m128 Vec0 = _mm_shuffle_ps(Temp0, Temp0, _MM_SHUFFLE(2, 2, 2, 0));
		
		// m[1][1]
		// m[0][1]
		// m[0][1]
		// m[0][1]
		__m128 Temp1 = _mm_shuffle_ps(in[1], in[0], _MM_SHUFFLE(1, 1, 1, 1));
		__m128 Vec1 = _mm_shuffle_ps(Temp1, Temp1, _MM_SHUFFLE(2, 2, 2, 0));
		
		// m[1][2]
		// m[0][2]
		// m[0][2]
		// m[0][2]
		__m128 Temp2 = _mm_shuffle_ps(in[1], in[0], _MM_SHUFFLE(2, 2, 2, 2));
		__m128 Vec2 = _mm_shuffle_ps(Temp2, Temp2, _MM_SHUFFLE(2, 2, 2, 0));
		
		// m[1][3]
		// m[0][3]
		// m[0][3]
		// m[0][3]
		__m128 Temp3 = _mm_shuffle_ps(in[1], in[0], _MM_SHUFFLE(3, 3, 3, 3));
		__m128 Vec3 = _mm_shuffle_ps(Temp3, Temp3, _MM_SHUFFLE(2, 2, 2, 0));
		
		// col0
		// + (Vec1[0] * Fac0[0] - Vec2[0] * Fac1[0] + Vec3[0] * Fac2[0]),
		// - (Vec1[1] * Fac0[1] - Vec2[1] * Fac1[1] + Vec3[1] * Fac2[1]),
		// + (Vec1[2] * Fac0[2] - Vec2[2] * Fac1[2] + Vec3[2] * Fac2[2]),
		// - (Vec1[3] * Fac0[3] - Vec2[3] * Fac1[3] + Vec3[3] * Fac2[3]),
		__m128 Mul00 = _mm_mul_ps(Vec1, Fac0);
		__m128 Mul01 = _mm_mul_ps(Vec2, Fac1);
		__m128 Mul02 = _mm_mul_ps(Vec3, Fac2);
		__m128 Sub00 = _mm_sub_ps(Mul00, Mul01);
		__m128 Add00 = _mm_add_ps(Sub00, Mul02);
		__m128 Inv0 = _mm_mul_ps(SignB, Add00);
		
		// col1
		// - (Vec0[0] * Fac0[0] - Vec2[0] * Fac3[0] + Vec3[0] * Fac4[0]),
		// + (Vec0[0] * Fac0[1] - Vec2[1] * Fac3[1] + Vec3[1] * Fac4[1]),
		// - (Vec0[0] * Fac0[2] - Vec2[2] * Fac3[2] + Vec3[2] * Fac4[2]),
		// + (Vec0[0] * Fac0[3] - Vec2[3] * Fac3[3] + Vec3[3] * Fac4[3]),
		__m128 Mul03 = _mm_mul_ps(Vec0, Fac0);
		__m128 Mul04 = _mm_mul_ps(Vec2, Fac3);
		__m128 Mul05 = _mm_mul_ps(Vec3, Fac4);
		__m128 Sub01 = _mm_sub_ps(Mul03, Mul04);
		__m128 Add01 = _mm_add_ps(Sub01, Mul05);
		__m128 Inv1 = _mm_mul_ps(SignA, Add01);
		
		// col2
		// + (Vec0[0] * Fac1[0] - Vec1[0] * Fac3[0] + Vec3[0] * Fac5[0]),
		// - (Vec0[0] * Fac1[1] - Vec1[1] * Fac3[1] + Vec3[1] * Fac5[1]),
		// + (Vec0[0] * Fac1[2] - Vec1[2] * Fac3[2] + Vec3[2] * Fac5[2]),
		// - (Vec0[0] * Fac1[3] - Vec1[3] * Fac3[3] + Vec3[3] * Fac5[3]),
		__m128 Mul06 = _mm_mul_ps(Vec0, Fac1);
		__m128 Mul07 = _mm_mul_ps(Vec1, Fac3);
		__m128 Mul08 = _mm_mul_ps(Vec3, Fac5);
		__m128 Sub02 = _mm_sub_ps(Mul06, Mul07);
		__m128 Add02 = _mm_add_ps(Sub02, Mul08);
		__m128 Inv2 = _mm_mul_ps(SignB, Add02);
		
		// col3
		// - (Vec1[0] * Fac2[0] - Vec1[0] * Fac4[0] + Vec2[0] * Fac5[0]),
		// + (Vec1[0] * Fac2[1] - Vec1[1] * Fac4[1] + Vec2[1] * Fac5[1]),
		// - (Vec1[0] * Fac2[2] - Vec1[2] * Fac4[2] + Vec2[2] * Fac5[2]),
		// + (Vec1[0] * Fac2[3] - Vec1[3] * Fac4[3] + Vec2[3] * Fac5[3]));
		__m128 Mul09 = _mm_mul_ps(Vec0, Fac2);
		__m128 Mul10 = _mm_mul_ps(Vec1, Fac4);
		__m128 Mul11 = _mm_mul_ps(Vec2, Fac5);
		__m128 Sub03 = _mm_sub_ps(Mul09, Mul10);
		__m128 Add03 = _mm_add_ps(Sub03, Mul11);
		__m128 Inv3 = _mm_mul_ps(SignA, Add03);
		
		__m128 Row0 = _mm_shuffle_ps(Inv0, Inv1, _MM_SHUFFLE(0, 0, 0, 0));
		__m128 Row1 = _mm_shuffle_ps(Inv2, Inv3, _MM_SHUFFLE(0, 0, 0, 0));
		__m128 Row2 = _mm_shuffle_ps(Row0, Row1, _MM_SHUFFLE(2, 0, 2, 0));
		
		//	valType Determinant = m[0][0] * Inverse[0][0]
		//						+ m[0][1] * Inverse[1][0]
		//						+ m[0][2] * Inverse[2][0]
		//						+ m[0][3] * Inverse[3][0];
		__m128 Det0 = SIMD_vec4_dot(in[0], Row2);
		__m128 Rcp0 = _mm_div_ps(_mm_set1_ps(1.0f), Det0);
		//__m128 Rcp0 = _mm_rcp_ps(Det0);
		
		//	Inverse /= Determinant;
		out[0] = _mm_mul_ps(Inv0, Rcp0);
		out[1] = _mm_mul_ps(Inv1, Rcp0);
		out[2] = _mm_mul_ps(Inv2, Rcp0);
		out[3] = _mm_mul_ps(Inv3, Rcp0);
	}
	
	SIMD_STATIC_INLINE void SIMD_mat4_inverse_lowp(SIMD_vec4 const in[4], SIMD_vec4 out[4])
	{
		__m128 Fac0;
		{
			//	valType SubFactor00 = m[2][2] * m[3][3] - m[3][2] * m[2][3];
			//	valType SubFactor00 = m[2][2] * m[3][3] - m[3][2] * m[2][3];
			//	valType SubFactor06 = m[1][2] * m[3][3] - m[3][2] * m[1][3];
			//	valType SubFactor13 = m[1][2] * m[2][3] - m[2][2] * m[1][3];
			
			__m128 Swp0a = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(3, 3, 3, 3));
			__m128 Swp0b = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(2, 2, 2, 2));
			
			__m128 Swp00 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(2, 2, 2, 2));
			__m128 Swp01 = _mm_shuffle_ps(Swp0a, Swp0a, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp02 = _mm_shuffle_ps(Swp0b, Swp0b, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp03 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(3, 3, 3, 3));
			
			__m128 Mul00 = _mm_mul_ps(Swp00, Swp01);
			__m128 Mul01 = _mm_mul_ps(Swp02, Swp03);
			Fac0 = _mm_sub_ps(Mul00, Mul01);
		}
		
		__m128 Fac1;
		{
			//	valType SubFactor01 = m[2][1] * m[3][3] - m[3][1] * m[2][3];
			//	valType SubFactor01 = m[2][1] * m[3][3] - m[3][1] * m[2][3];
			//	valType SubFactor07 = m[1][1] * m[3][3] - m[3][1] * m[1][3];
			//	valType SubFactor14 = m[1][1] * m[2][3] - m[2][1] * m[1][3];
			
			__m128 Swp0a = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(3, 3, 3, 3));
			__m128 Swp0b = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(1, 1, 1, 1));
			
			__m128 Swp00 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(1, 1, 1, 1));
			__m128 Swp01 = _mm_shuffle_ps(Swp0a, Swp0a, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp02 = _mm_shuffle_ps(Swp0b, Swp0b, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp03 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(3, 3, 3, 3));
			
			__m128 Mul00 = _mm_mul_ps(Swp00, Swp01);
			__m128 Mul01 = _mm_mul_ps(Swp02, Swp03);
			Fac1 = _mm_sub_ps(Mul00, Mul01);
		}
		
		
		__m128 Fac2;
		{
			//	valType SubFactor02 = m[2][1] * m[3][2] - m[3][1] * m[2][2];
			//	valType SubFactor02 = m[2][1] * m[3][2] - m[3][1] * m[2][2];
			//	valType SubFactor08 = m[1][1] * m[3][2] - m[3][1] * m[1][2];
			//	valType SubFactor15 = m[1][1] * m[2][2] - m[2][1] * m[1][2];
			
			__m128 Swp0a = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(2, 2, 2, 2));
			__m128 Swp0b = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(1, 1, 1, 1));
			
			__m128 Swp00 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(1, 1, 1, 1));
			__m128 Swp01 = _mm_shuffle_ps(Swp0a, Swp0a, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp02 = _mm_shuffle_ps(Swp0b, Swp0b, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp03 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(2, 2, 2, 2));
			
			__m128 Mul00 = _mm_mul_ps(Swp00, Swp01);
			__m128 Mul01 = _mm_mul_ps(Swp02, Swp03);
			Fac2 = _mm_sub_ps(Mul00, Mul01);
		}
		
		__m128 Fac3;
		{
			//	valType SubFactor03 = m[2][0] * m[3][3] - m[3][0] * m[2][3];
			//	valType SubFactor03 = m[2][0] * m[3][3] - m[3][0] * m[2][3];
			//	valType SubFactor09 = m[1][0] * m[3][3] - m[3][0] * m[1][3];
			//	valType SubFactor16 = m[1][0] * m[2][3] - m[2][0] * m[1][3];
			
			__m128 Swp0a = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(3, 3, 3, 3));
			__m128 Swp0b = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(0, 0, 0, 0));
			
			__m128 Swp00 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(0, 0, 0, 0));
			__m128 Swp01 = _mm_shuffle_ps(Swp0a, Swp0a, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp02 = _mm_shuffle_ps(Swp0b, Swp0b, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp03 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(3, 3, 3, 3));
			
			__m128 Mul00 = _mm_mul_ps(Swp00, Swp01);
			__m128 Mul01 = _mm_mul_ps(Swp02, Swp03);
			Fac3 = _mm_sub_ps(Mul00, Mul01);
		}
		
		__m128 Fac4;
		{
			//	valType SubFactor04 = m[2][0] * m[3][2] - m[3][0] * m[2][2];
			//	valType SubFactor04 = m[2][0] * m[3][2] - m[3][0] * m[2][2];
			//	valType SubFactor10 = m[1][0] * m[3][2] - m[3][0] * m[1][2];
			//	valType SubFactor17 = m[1][0] * m[2][2] - m[2][0] * m[1][2];
			
			__m128 Swp0a = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(2, 2, 2, 2));
			__m128 Swp0b = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(0, 0, 0, 0));
			
			__m128 Swp00 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(0, 0, 0, 0));
			__m128 Swp01 = _mm_shuffle_ps(Swp0a, Swp0a, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp02 = _mm_shuffle_ps(Swp0b, Swp0b, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp03 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(2, 2, 2, 2));
			
			__m128 Mul00 = _mm_mul_ps(Swp00, Swp01);
			__m128 Mul01 = _mm_mul_ps(Swp02, Swp03);
			Fac4 = _mm_sub_ps(Mul00, Mul01);
		}
		
		__m128 Fac5;
		{
			//	valType SubFactor05 = m[2][0] * m[3][1] - m[3][0] * m[2][1];
			//	valType SubFactor05 = m[2][0] * m[3][1] - m[3][0] * m[2][1];
			//	valType SubFactor12 = m[1][0] * m[3][1] - m[3][0] * m[1][1];
			//	valType SubFactor18 = m[1][0] * m[2][1] - m[2][0] * m[1][1];
			
			__m128 Swp0a = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(1, 1, 1, 1));
			__m128 Swp0b = _mm_shuffle_ps(in[3], in[2], _MM_SHUFFLE(0, 0, 0, 0));
			
			__m128 Swp00 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(0, 0, 0, 0));
			__m128 Swp01 = _mm_shuffle_ps(Swp0a, Swp0a, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp02 = _mm_shuffle_ps(Swp0b, Swp0b, _MM_SHUFFLE(2, 0, 0, 0));
			__m128 Swp03 = _mm_shuffle_ps(in[2], in[1], _MM_SHUFFLE(1, 1, 1, 1));
			
			__m128 Mul00 = _mm_mul_ps(Swp00, Swp01);
			__m128 Mul01 = _mm_mul_ps(Swp02, Swp03);
			Fac5 = _mm_sub_ps(Mul00, Mul01);
		}
		
		__m128 SignA = _mm_set_ps( 1.0f,-1.0f, 1.0f,-1.0f);
		__m128 SignB = _mm_set_ps(-1.0f, 1.0f,-1.0f, 1.0f);
		
		// m[1][0]
		// m[0][0]
		// m[0][0]
		// m[0][0]
		__m128 Temp0 = _mm_shuffle_ps(in[1], in[0], _MM_SHUFFLE(0, 0, 0, 0));
		__m128 Vec0 = _mm_shuffle_ps(Temp0, Temp0, _MM_SHUFFLE(2, 2, 2, 0));
		
		// m[1][1]
		// m[0][1]
		// m[0][1]
		// m[0][1]
		__m128 Temp1 = _mm_shuffle_ps(in[1], in[0], _MM_SHUFFLE(1, 1, 1, 1));
		__m128 Vec1 = _mm_shuffle_ps(Temp1, Temp1, _MM_SHUFFLE(2, 2, 2, 0));
		
		// m[1][2]
		// m[0][2]
		// m[0][2]
		// m[0][2]
		__m128 Temp2 = _mm_shuffle_ps(in[1], in[0], _MM_SHUFFLE(2, 2, 2, 2));
		__m128 Vec2 = _mm_shuffle_ps(Temp2, Temp2, _MM_SHUFFLE(2, 2, 2, 0));
		
		// m[1][3]
		// m[0][3]
		// m[0][3]
		// m[0][3]
		__m128 Temp3 = _mm_shuffle_ps(in[1], in[0], _MM_SHUFFLE(3, 3, 3, 3));
		__m128 Vec3 = _mm_shuffle_ps(Temp3, Temp3, _MM_SHUFFLE(2, 2, 2, 0));
		
		// col0
		// + (Vec1[0] * Fac0[0] - Vec2[0] * Fac1[0] + Vec3[0] * Fac2[0]),
		// - (Vec1[1] * Fac0[1] - Vec2[1] * Fac1[1] + Vec3[1] * Fac2[1]),
		// + (Vec1[2] * Fac0[2] - Vec2[2] * Fac1[2] + Vec3[2] * Fac2[2]),
		// - (Vec1[3] * Fac0[3] - Vec2[3] * Fac1[3] + Vec3[3] * Fac2[3]),
		__m128 Mul00 = _mm_mul_ps(Vec1, Fac0);
		__m128 Mul01 = _mm_mul_ps(Vec2, Fac1);
		__m128 Mul02 = _mm_mul_ps(Vec3, Fac2);
		__m128 Sub00 = _mm_sub_ps(Mul00, Mul01);
		__m128 Add00 = _mm_add_ps(Sub00, Mul02);
		__m128 Inv0 = _mm_mul_ps(SignB, Add00);
		
		// col1
		// - (Vec0[0] * Fac0[0] - Vec2[0] * Fac3[0] + Vec3[0] * Fac4[0]),
		// + (Vec0[0] * Fac0[1] - Vec2[1] * Fac3[1] + Vec3[1] * Fac4[1]),
		// - (Vec0[0] * Fac0[2] - Vec2[2] * Fac3[2] + Vec3[2] * Fac4[2]),
		// + (Vec0[0] * Fac0[3] - Vec2[3] * Fac3[3] + Vec3[3] * Fac4[3]),
		__m128 Mul03 = _mm_mul_ps(Vec0, Fac0);
		__m128 Mul04 = _mm_mul_ps(Vec2, Fac3);
		__m128 Mul05 = _mm_mul_ps(Vec3, Fac4);
		__m128 Sub01 = _mm_sub_ps(Mul03, Mul04);
		__m128 Add01 = _mm_add_ps(Sub01, Mul05);
		__m128 Inv1 = _mm_mul_ps(SignA, Add01);
		
		// col2
		// + (Vec0[0] * Fac1[0] - Vec1[0] * Fac3[0] + Vec3[0] * Fac5[0]),
		// - (Vec0[0] * Fac1[1] - Vec1[1] * Fac3[1] + Vec3[1] * Fac5[1]),
		// + (Vec0[0] * Fac1[2] - Vec1[2] * Fac3[2] + Vec3[2] * Fac5[2]),
		// - (Vec0[0] * Fac1[3] - Vec1[3] * Fac3[3] + Vec3[3] * Fac5[3]),
		__m128 Mul06 = _mm_mul_ps(Vec0, Fac1);
		__m128 Mul07 = _mm_mul_ps(Vec1, Fac3);
		__m128 Mul08 = _mm_mul_ps(Vec3, Fac5);
		__m128 Sub02 = _mm_sub_ps(Mul06, Mul07);
		__m128 Add02 = _mm_add_ps(Sub02, Mul08);
		__m128 Inv2 = _mm_mul_ps(SignB, Add02);
		
		// col3
		// - (Vec1[0] * Fac2[0] - Vec1[0] * Fac4[0] + Vec2[0] * Fac5[0]),
		// + (Vec1[0] * Fac2[1] - Vec1[1] * Fac4[1] + Vec2[1] * Fac5[1]),
		// - (Vec1[0] * Fac2[2] - Vec1[2] * Fac4[2] + Vec2[2] * Fac5[2]),
		// + (Vec1[0] * Fac2[3] - Vec1[3] * Fac4[3] + Vec2[3] * Fac5[3]));
		__m128 Mul09 = _mm_mul_ps(Vec0, Fac2);
		__m128 Mul10 = _mm_mul_ps(Vec1, Fac4);
		__m128 Mul11 = _mm_mul_ps(Vec2, Fac5);
		__m128 Sub03 = _mm_sub_ps(Mul09, Mul10);
		__m128 Add03 = _mm_add_ps(Sub03, Mul11);
		__m128 Inv3 = _mm_mul_ps(SignA, Add03);
		
		__m128 Row0 = _mm_shuffle_ps(Inv0, Inv1, _MM_SHUFFLE(0, 0, 0, 0));
		__m128 Row1 = _mm_shuffle_ps(Inv2, Inv3, _MM_SHUFFLE(0, 0, 0, 0));
		__m128 Row2 = _mm_shuffle_ps(Row0, Row1, _MM_SHUFFLE(2, 0, 2, 0));
		
		//	valType Determinant = m[0][0] * Inverse[0][0]
		//						+ m[0][1] * Inverse[1][0]
		//						+ m[0][2] * Inverse[2][0]
		//						+ m[0][3] * Inverse[3][0];
		__m128 Det0 = SIMD_vec4_dot(in[0], Row2);
		__m128 Rcp0 = _mm_rcp_ps(Det0);
		//__m128 Rcp0 = _mm_div_ps(one, Det0);
		//	Inverse /= Determinant;
		out[0] = _mm_mul_ps(Inv0, Rcp0);
		out[1] = _mm_mul_ps(Inv1, Rcp0);
		out[2] = _mm_mul_ps(Inv2, Rcp0);
		out[3] = _mm_mul_ps(Inv3, Rcp0);
	}
	/*
	 FORCE_INLINE void SIMD_mat4_rotate(__m128 const in[4], float Angle, float const v[3], __m128 out[4])
	 {
	 float a = SML::radians(Angle);
	 float c = cos(a);
	 float s = sin(a);
	 
	 SML::vec4 AxisA(v[0], v[1], v[2], float(0));
	 __m128 AxisB = _mm_set_ps(AxisA.w, AxisA.z, AxisA.y, AxisA.x);
	 __m128 AxisC = detail::sse_nrm_ps(AxisB);
	 
	 __m128 Cos0 = _mm_set_ss(c);
	 __m128 CosA = _mm_shuffle_ps(Cos0, Cos0, _MM_SHUFFLE(0, 0, 0, 0));
	 __m128 Sin0 = _mm_set_ss(s);
	 __m128 SinA = _mm_shuffle_ps(Sin0, Sin0, _MM_SHUFFLE(0, 0, 0, 0));
	 
	 // tvec3<T, P> temp = (valType(1) - c) * axis;
	 __m128 Temp0 = _mm_sub_ps(one, CosA);
	 __m128 Temp1 = _mm_mul_ps(Temp0, AxisC);
	 
	 //Rotate[0][0] = c + temp[0] * axis[0];
	 //Rotate[0][1] = 0 + temp[0] * axis[1] + s * axis[2];
	 //Rotate[0][2] = 0 + temp[0] * axis[2] - s * axis[1];
	 __m128 Axis0 = _mm_shuffle_ps(AxisC, AxisC, _MM_SHUFFLE(0, 0, 0, 0));
	 __m128 TmpA0 = _mm_mul_ps(Axis0, AxisC);
	 __m128 CosA0 = _mm_shuffle_ps(Cos0, Cos0, _MM_SHUFFLE(1, 1, 1, 0));
	 __m128 TmpA1 = _mm_add_ps(CosA0, TmpA0);
	 __m128 SinA0 = SinA;//_mm_set_ps(0.0f, s, -s, 0.0f);
	 __m128 TmpA2 = _mm_shuffle_ps(AxisC, AxisC, _MM_SHUFFLE(3, 1, 2, 3));
	 __m128 TmpA3 = _mm_mul_ps(SinA0, TmpA2);
	 __m128 TmpA4 = _mm_add_ps(TmpA1, TmpA3);
	 
	 //Rotate[1][0] = 0 + temp[1] * axis[0] - s * axis[2];
	 //Rotate[1][1] = c + temp[1] * axis[1];
	 //Rotate[1][2] = 0 + temp[1] * axis[2] + s * axis[0];
	 __m128 Axis1 = _mm_shuffle_ps(AxisC, AxisC, _MM_SHUFFLE(1, 1, 1, 1));
	 __m128 TmpB0 = _mm_mul_ps(Axis1, AxisC);
	 __m128 CosA1 = _mm_shuffle_ps(Cos0, Cos0, _MM_SHUFFLE(1, 1, 0, 1));
	 __m128 TmpB1 = _mm_add_ps(CosA1, TmpB0);
	 __m128 SinB0 = SinA;//_mm_set_ps(-s, 0.0f, s, 0.0f);
	 __m128 TmpB2 = _mm_shuffle_ps(AxisC, AxisC, _MM_SHUFFLE(3, 0, 3, 2));
	 __m128 TmpB3 = _mm_mul_ps(SinA0, TmpB2);
	 __m128 TmpB4 = _mm_add_ps(TmpB1, TmpB3);
	 
	 //Rotate[2][0] = 0 + temp[2] * axis[0] + s * axis[1];
	 //Rotate[2][1] = 0 + temp[2] * axis[1] - s * axis[0];
	 //Rotate[2][2] = c + temp[2] * axis[2];
	 __m128 Axis2 = _mm_shuffle_ps(AxisC, AxisC, _MM_SHUFFLE(2, 2, 2, 2));
	 __m128 TmpC0 = _mm_mul_ps(Axis2, AxisC);
	 __m128 CosA2 = _mm_shuffle_ps(Cos0, Cos0, _MM_SHUFFLE(1, 0, 1, 1));
	 __m128 TmpC1 = _mm_add_ps(CosA2, TmpC0);
	 __m128 SinC0 = SinA;//_mm_set_ps(s, -s, 0.0f, 0.0f);
	 __m128 TmpC2 = _mm_shuffle_ps(AxisC, AxisC, _MM_SHUFFLE(3, 3, 0, 1));
	 __m128 TmpC3 = _mm_mul_ps(SinA0, TmpC2);
	 __m128 TmpC4 = _mm_add_ps(TmpC1, TmpC3);
	 
	 __m128 Result[4];
	 Result[0] = TmpA4;
	 Result[1] = TmpB4;
	 Result[2] = TmpC4;
	 Result[3] = _mm_set_ps(1, 0, 0, 0);
	 
	 //tmat4x4<valType> Result(uninitialize);
	 //Result[0] = m[0] * Rotate[0][0] + m[1] * Rotate[0][1] + m[2] * Rotate[0][2];
	 //Result[1] = m[0] * Rotate[1][0] + m[1] * Rotate[1][1] + m[2] * Rotate[1][2];
	 //Result[2] = m[0] * Rotate[2][0] + m[1] * Rotate[2][1] + m[2] * Rotate[2][2];
	 //Result[3] = m[3];
	 //return Result;
	 sse_mul_ps(in, Result, out);
	 }
	 */
	SIMD_STATIC_INLINE void SIMD_mat4_outerProduct(__m128 const & c, __m128 const & r, __m128 out[4])
	{
		out[0] = _mm_mul_ps(c, _mm_shuffle_ps(r, r, _MM_SHUFFLE(0, 0, 0, 0)));
		out[1] = _mm_mul_ps(c, _mm_shuffle_ps(r, r, _MM_SHUFFLE(1, 1, 1, 1)));
		out[2] = _mm_mul_ps(c, _mm_shuffle_ps(r, r, _MM_SHUFFLE(2, 2, 2, 2)));
		out[3] = _mm_mul_ps(c, _mm_shuffle_ps(r, r, _MM_SHUFFLE(3, 3, 3, 3)));
	}
	
#endif//CPU_ARCH & CPU_ARCH_SSE2_BIT
	
	
	
	
	
	
	
	
	template<typename T>
	mat_t<4, 4, T> translate(mat_t<4, 4, T> const& m, pack_t<3, T> const& v)
	{
		mat_t<4, 4, T> Result(m);
		Result[3] = m[0] * v[0] + m[1] * v[1] + m[2] * v[2] + m[3];
		return Result;
	}
	
	template<typename T>
	mat_t<4, 4, T> rotate(mat_t<4, 4, T> const& m, T angle, pack_t<3, T> const& v)
	{
		T const a = angle;
		T const c = cos(a);
		T const s = sin(a);
		
		pack_t<3, T> axis(normalize(v));
		pack_t<3, T> temp((T(1) - c) * axis);
		
		mat_t<4, 4, T> Rotate;
		Rotate[0][0] = c + temp[0] * axis[0];
		Rotate[0][1] = temp[0] * axis[1] + s * axis[2];
		Rotate[0][2] = temp[0] * axis[2] - s * axis[1];
		
		Rotate[1][0] = temp[1] * axis[0] - s * axis[2];
		Rotate[1][1] = c + temp[1] * axis[1];
		Rotate[1][2] = temp[1] * axis[2] + s * axis[0];
		
		Rotate[2][0] = temp[2] * axis[0] + s * axis[1];
		Rotate[2][1] = temp[2] * axis[1] - s * axis[0];
		Rotate[2][2] = c + temp[2] * axis[2];
		
		mat_t<4, 4, T> Result;
		Result[0] = m[0] * Rotate[0][0] + m[1] * Rotate[0][1] + m[2] * Rotate[0][2];
		Result[1] = m[0] * Rotate[1][0] + m[1] * Rotate[1][1] + m[2] * Rotate[1][2];
		Result[2] = m[0] * Rotate[2][0] + m[1] * Rotate[2][1] + m[2] * Rotate[2][2];
		Result[3] = m[3];
		return Result;
	}
	
	template<typename T>
	mat_t<4, 4, T> rotate_slow(mat_t<4, 4, T> const& m, T angle, pack_t<3, T> const& v)
	{
		T const a = angle;
		T const c = cos(a);
		T const s = sin(a);
		mat_t<4, 4, T> Result;
		
		pack_t<3, T> axis = normalize(v);
		
		Result[0][0] = c + (static_cast<T>(1) - c)      * axis.x     * axis.x;
		Result[0][1] = (static_cast<T>(1) - c) * axis.x * axis.y + s * axis.z;
		Result[0][2] = (static_cast<T>(1) - c) * axis.x * axis.z - s * axis.y;
		Result[0][3] = static_cast<T>(0);
		
		Result[1][0] = (static_cast<T>(1) - c) * axis.y * axis.x - s * axis.z;
		Result[1][1] = c + (static_cast<T>(1) - c) * axis.y * axis.y;
		Result[1][2] = (static_cast<T>(1) - c) * axis.y * axis.z + s * axis.x;
		Result[1][3] = static_cast<T>(0);
		
		Result[2][0] = (static_cast<T>(1) - c) * axis.z * axis.x + s * axis.y;
		Result[2][1] = (static_cast<T>(1) - c) * axis.z * axis.y - s * axis.x;
		Result[2][2] = c + (static_cast<T>(1) - c) * axis.z * axis.z;
		Result[2][3] = static_cast<T>(0);
		
		Result[3] = pack_t<4, T>(0, 0, 0, 1);
		return m * Result;
	}
	
	template<typename T>
	mat_t<4, 4, T> scale(mat_t<4, 4, T> const& m, pack_t<3, T> const& v)
	{
		mat_t<4, 4, T> Result;
		Result[0] = m[0] * v[0];
		Result[1] = m[1] * v[1];
		Result[2] = m[2] * v[2];
		Result[3] = m[3];
		return Result;
	}
	
	template<typename T>
	mat_t<4, 4, T> scale_slow(mat_t<4, 4, T> const& m, pack_t<3, T> const& v)
	{
		mat_t<4, 4, T> Result(T(1));
		Result[0][0] = v.x;
		Result[1][1] = v.y;
		Result[2][2] = v.z;
		return m * Result;
	}

	template<typename T, typename U>
	pack_t<3, T> projectZO(pack_t<3, T> const& obj, mat_t<4, 4, T> const& model, mat_t<4, 4, T> const& proj, pack_t<4, U> const& viewport)
	{
		pack_t<4, T> tmp = pack_t<4, T>(obj, static_cast<T>(1));
		tmp = model * tmp;
		tmp = proj * tmp;
		
		tmp /= tmp.w;
		tmp.x = tmp.x * static_cast<T>(0.5) + static_cast<T>(0.5);
		tmp.y = tmp.y * static_cast<T>(0.5) + static_cast<T>(0.5);
		
		tmp[0] = tmp[0] * T(viewport[2]) + T(viewport[0]);
		tmp[1] = tmp[1] * T(viewport[3]) + T(viewport[1]);
		
		return pack_t<3, T>(tmp);
	}
	
	template<typename T, typename U>
	pack_t<3, T> projectNO(pack_t<3, T> const& obj, mat_t<4, 4, T> const& model, mat_t<4, 4, T> const& proj, pack_t<4, U> const& viewport)
	{
		pack_t<4, T> tmp = pack_t<4, T>(obj, static_cast<T>(1));
		tmp = model * tmp;
		tmp = proj * tmp;
		
		tmp /= tmp.w;
		tmp = tmp * static_cast<T>(0.5) + static_cast<T>(0.5);
		tmp[0] = tmp[0] * T(viewport[2]) + T(viewport[0]);
		tmp[1] = tmp[1] * T(viewport[3]) + T(viewport[1]);
		
		return pack_t<3, T>(tmp);
	}
	
	template<typename T, typename U>
	pack_t<3, T> project(pack_t<3, T> const& obj, mat_t<4, 4, T> const& model, mat_t<4, 4, T> const& proj, pack_t<4, U> const& viewport)
	{
#		if SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_ZERO_TO_ONE
		return projectZO(obj, model, proj, viewport);
#		else
		return projectNO(obj, model, proj, viewport);
#		endif
	}
	
	template<typename T, typename U>
	pack_t<3, T> unProjectZO(pack_t<3, T> const& win, mat_t<4, 4, T> const& model, mat_t<4, 4, T> const& proj, pack_t<4, U> const& viewport)
	{
		mat_t<4, 4, T> Inverse = inverse(proj * model);
		
		pack_t<4, T> tmp = pack_t<4, T>(win, T(1));
		tmp.x = (tmp.x - T(viewport[0])) / T(viewport[2]);
		tmp.y = (tmp.y - T(viewport[1])) / T(viewport[3]);
		tmp.x = tmp.x * static_cast<T>(2) - static_cast<T>(1);
		tmp.y = tmp.y * static_cast<T>(2) - static_cast<T>(1);
		
		pack_t<4, T> obj = Inverse * tmp;
		obj /= obj.w;
		
		return pack_t<3, T>(obj);
	}
	
	template<typename T, typename U>
	pack_t<3, T> unProjectNO(pack_t<3, T> const& win, mat_t<4, 4, T> const& model, mat_t<4, 4, T> const& proj, pack_t<4, U> const& viewport)
	{
		mat_t<4, 4, T> Inverse = inverse(proj * model);
		
		pack_t<4, T> tmp = pack_t<4, T>(win, T(1));
		tmp.x = (tmp.x - T(viewport[0])) / T(viewport[2]);
		tmp.y = (tmp.y - T(viewport[1])) / T(viewport[3]);
		tmp = tmp * static_cast<T>(2) - static_cast<T>(1);
		
		pack_t<4, T> obj = Inverse * tmp;
		obj /= obj.w;
		
		return pack_t<3, T>(obj);
	}
	
	template<typename T, typename U>
	pack_t<3, T> unProject(pack_t<3, T> const& win, mat_t<4, 4, T> const& model, mat_t<4, 4, T> const& proj, pack_t<4, U> const& viewport)
	{
#		if SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_ZERO_TO_ONE
		return unProjectZO(win, model, proj, viewport);
#		else
		return unProjectNO(win, model, proj, viewport);
#		endif
	}
	
	template<typename T, bool N, typename U>
	mat_t<4, 4, T> pickMatrix(pack_t<2, T> const& center, pack_t<2, T> const& delta, pack_t<4, U> const& viewport)
	{
		assert(delta.x > static_cast<T>(0) && delta.y > static_cast<T>(0));
		mat_t<4, 4, T> Result(static_cast<T>(1));
		
		if(!(delta.x > static_cast<T>(0) && delta.y > static_cast<T>(0)))
			return Result; // Error
		
		pack_t<3, T> Temp(
						  (static_cast<T>(viewport[2]) - static_cast<T>(2) * (center.x - static_cast<T>(viewport[0]))) / delta.x,
						  (static_cast<T>(viewport[3]) - static_cast<T>(2) * (center.y - static_cast<T>(viewport[1]))) / delta.y,
						  static_cast<T>(0));
		
		// Translate and scale the picked region to the entire window
		Result = translate(Result, Temp);
		return scale(Result, pack_t<3, T>(static_cast<T>(viewport[2]) / delta.x, static_cast<T>(viewport[3]) / delta.y, static_cast<T>(1)));
	}
	

}



namespace sml
{
	
	// -- Constructors --

//
//	// -- Matrix conversions --
//
//	template<typename T>
//	mat_t<4, 4, T>::mat_t(mat_t<2, 2, T> const& m)
//	{
//		this->value[0] = col_type(m[0], 0, 0);
//		this->value[1] = col_type(m[1], 0, 0);
//		this->value[2] = col_type(0, 0, 1, 0);
//		this->value[3] = col_type(0, 0, 0, 1);
//	}
//
//	template<typename T>
//	mat_t<4, 4, T>::mat_t(mat_t<3, 3, T> const& m)
//	{
//		this->value[0] = col_type(m[0], 0);
//		this->value[1] = col_type(m[1], 0);
//		this->value[2] = col_type(m[2], 0);
//		this->value[3] = col_type(0, 0, 0, 1);
//	}
//
//	template<typename T>
//	mat_t<4, 4, T>::mat_t(mat_t<2, 3, T> const& m)
//	{
//		this->value[0] = col_type(m[0], 0);
//		this->value[1] = col_type(m[1], 0);
//		this->value[2] = col_type(0, 0, 1, 0);
//		this->value[3] = col_type(0, 0, 0, 1);
//	}
//
//	template<typename T>
//	mat_t<4, 4, T>::mat_t(mat_t<3, 2, T> const& m)
//	{
//		this->value[0] = col_type(m[0], 0, 0);
//		this->value[1] = col_type(m[1], 0, 0);
//		this->value[2] = col_type(m[2], 1, 0);
//		this->value[3] = col_type(0, 0, 0, 1);
//	}
//
//	template<typename T>
//	mat_t<4, 4, T>::mat_t(mat_t<2, 4, T> const& m)
//	{
//		this->value[0] = m[0];
//		this->value[1] = m[1];
//		this->value[2] = col_type(0, 0, 1, 0);
//		this->value[3] = col_type(0, 0, 0, 1);
//	}
//
//	template<typename T>
//	mat_t<4, 4, T>::mat_t(mat_t<4, 2, T> const& m)
//	{
//		this->value[0] = col_type(m[0], 0, 0);
//		this->value[1] = col_type(m[1], 0, 0);
//		this->value[2] = col_type(0, 0, 1, 0);
//		this->value[3] = col_type(0, 0, 0, 1);
//	}
//
//	template<typename T>
//	mat_t<4, 4, T>::mat_t(mat_t<3, 4, T> const& m)
//	{
//		this->value[0] = m[0];
//		this->value[1] = m[1];
//		this->value[2] = m[2];
//		this->value[3] = col_type(0, 0, 0, 1);
//	}
//
//	template<typename T>
//	mat_t<4, 4, T>::mat_t(mat_t<4, 3, T> const& m)
//	{
//		this->value[0] = col_type(m[0], 0);
//		this->value[1] = col_type(m[1], 0);
//		this->value[2] = col_type(m[2], 0);
//		this->value[3] = col_type(m[3], 1);
//	}
	
	// -- Accesses --
	
//	template<typename T>
//	typename mat_t<4, 4, T>::col_type & mat_t<4, 4, T>::operator[](typename mat_t<4, 4, T>::length_type i)
//	{
//		assert(i < this->length());
//		return this->value[i];
//	}
//
//	template<typename T>
//	typename mat_t<4, 4, T>::col_type const& mat_t<4, 4, T>::operator[](typename mat_t<4, 4, T>::length_type i) const
//	{
//		assert(i < this->length());
//		return this->value[i];
//	}
	
	// -- Unary arithmetic operators --
	
//#	if !SIMD_HAS_DEFAULTED_FUNCTIONS
//#	endif//!SIMD_HAS_DEFAULTED_FUNCTIONS
//	
//	template<typename T>
//	template<typename U>
//	SIMD_FUNC_QUALIFIER SIMD_CONSTEXPR_CXX14 mat_t<4, 4, T>& mat_t<4, 4, T>::operator=(mat_t<4, 4, U> const& m)
//	{
//		//memcpy could be faster
//		//memcpy(&this->value, &m.value, 16 * sizeof(valType));
//		this->value[0] = m[0];
//		this->value[1] = m[1];
//		this->value[2] = m[2];
//		this->value[3] = m[3];
//		return *this;
//	}
//	
//	template<typename T>
//	template<typename U>
//	SIMD_FUNC_QUALIFIER mat_t<4, 4, T>& mat_t<4, 4, T>::operator+=(U s)
//	{
//		this->value[0] += s;
//		this->value[1] += s;
//		this->value[2] += s;
//		this->value[3] += s;
//		return *this;
//	}
//	
//	template<typename T>
//	template<typename U>
//	SIMD_FUNC_QUALIFIER mat_t<4, 4, T>& mat_t<4, 4, T>::operator+=(mat_t<4, 4, U> const& m)
//	{
//		this->value[0] += m[0];
//		this->value[1] += m[1];
//		this->value[2] += m[2];
//		this->value[3] += m[3];
//		return *this;
//	}
//	
//	template<typename T>
//	template<typename U>
//	DLL_APICALL	mat_t<4, 4, T>	DLL_CALLING_CONVENTION & mat_t<4, 4, T>::operator-=(U s)
//	{
//		this->value[0] -= s;
//		this->value[1] -= s;
//		this->value[2] -= s;
//		this->value[3] -= s;
//		return *this;
//	}
//	
//	template<typename T>
//	template<typename U>
//	DLL_APICALL	mat_t<4, 4, T>	DLL_CALLING_CONVENTION & mat_t<4, 4, T>::operator-=(mat_t<4, 4, U> const& m)
//	{
//		this->value[0] -= m[0];
//		this->value[1] -= m[1];
//		this->value[2] -= m[2];
//		this->value[3] -= m[3];
//		return *this;
//	}
//	
//	template<typename T>
//	template<typename U>
//	DLL_APICALL	mat_t<4, 4, T>	DLL_CALLING_CONVENTION & mat_t<4, 4, T>::operator*=(U s)
//	{
//		this->value[0] *= s;
//		this->value[1] *= s;
//		this->value[2] *= s;
//		this->value[3] *= s;
//		return *this;
//	}
//	
//	template<typename T>
//	template<typename U>
//	DLL_APICALL	mat_t<4, 4, T>	DLL_CALLING_CONVENTION & mat_t<4, 4, T>::operator*=(mat_t<4, 4, U> const& m)
//	{
//		return (*this = *this * m);
//	}
//	
//	template<typename T>
//	template<typename U>
//	DLL_APICALL	mat_t<4, 4, T>	DLL_CALLING_CONVENTION & mat_t<4, 4, T>::operator/=(U s)
//	{
//		this->value[0] /= s;
//		this->value[1] /= s;
//		this->value[2] /= s;
//		this->value[3] /= s;
//		return *this;
//	}
//	
//	template<typename T>
//	template<typename U>
//	DLL_APICALL	mat_t<4, 4, T>	DLL_CALLING_CONVENTION & mat_t<4, 4, T>::operator/=(mat_t<4, 4, U> const& m)
//	{
//		return *this *= inverse(m);
//	}
//	
//	// -- Increment and decrement operators --
//	
//	template<typename T>
//	DLL_APICALL	mat_t<4, 4, T>	DLL_CALLING_CONVENTION & mat_t<4, 4, T>::operator++()
//	{
//		++this->value[0];
//		++this->value[1];
//		++this->value[2];
//		++this->value[3];
//		return *this;
//	}
//	
//	template<typename T>
//	DLL_APICALL	mat_t<4, 4, T>	DLL_CALLING_CONVENTION & mat_t<4, 4, T>::operator--()
//	{
//		--this->value[0];
//		--this->value[1];
//		--this->value[2];
//		--this->value[3];
//		return *this;
//	}
//	
//	template<typename T>
//	DLL_APICALL	mat_t<4, 4, T>	DLL_CALLING_CONVENTION mat_t<4, 4, T>::operator++(int)
//	{
//		mat_t<4, 4, T> Result(*this);
//		++*this;
//		return Result;
//	}
//	
//	template<typename T>
//	DLL_APICALL	mat_t<4, 4, T>	DLL_CALLING_CONVENTION mat_t<4, 4, T>::operator--(int)
//	{
//		mat_t<4, 4, T> Result(*this);
//		--*this;
//		return Result;
//	}
//	
//	// -- Unary constant operators --
//	
//	template<typename T>
//	DLL_APICALL	mat_t<4, 4, T>	DLL_CALLING_CONVENTION operator+(mat_t<4, 4, T> const& m)
//	{
//		return m;
//	}
//	
//	template<typename T>
//	DLL_APICALL	mat_t<4, 4, T>	DLL_CALLING_CONVENTION operator-(mat_t<4, 4, T> const& m)
//	{
//		return mat_t<4, 4, T>(
//							   -m[0],
//							   -m[1],
//							   -m[2],
//							   -m[3]);
//	}
//	
//	// -- Binary arithmetic operators --
//	
//	template<typename T>
//	DLL_APICALL	mat_t<4, 4, T>	DLL_CALLING_CONVENTION operator+(mat_t<4, 4, T> const& m, T const& s)
//	{
//		return mat_t<4, 4, T>(
//							   m[0] + s,
//							   m[1] + s,
//							   m[2] + s,
//							   m[3] + s);
//	}
//	
//	template<typename T>
//	DLL_APICALL	mat_t<4, 4, T>	DLL_CALLING_CONVENTION operator+(T const& s, mat_t<4, 4, T> const& m)
//	{
//		return mat_t<4, 4, T>(
//							   m[0] + s,
//							   m[1] + s,
//							   m[2] + s,
//							   m[3] + s);
//	}
//	
//	template<typename T>
//	DLL_APICALL	mat_t<4, 4, T>	DLL_CALLING_CONVENTION operator+(mat_t<4, 4, T> const& m1, mat_t<4, 4, T> const& m2)
//	{
//		return mat_t<4, 4, T>(
//							   m1[0] + m2[0],
//							   m1[1] + m2[1],
//							   m1[2] + m2[2],
//							   m1[3] + m2[3]);
//	}
//	
//	template<typename T>
//	DLL_APICALL	mat_t<4, 4, T>	DLL_CALLING_CONVENTION operator-(mat_t<4, 4, T> const& m, T const& s)
//	{
//		return mat_t<4, 4, T>(
//							   m[0] - s,
//							   m[1] - s,
//							   m[2] - s,
//							   m[3] - s);
//	}
//	
//	template<typename T>
//	DLL_APICALL	mat_t<4, 4, T>	DLL_CALLING_CONVENTION operator-(T const& s, mat_t<4, 4, T> const& m)
//	{
//		return mat_t<4, 4, T>(
//							   s - m[0],
//							   s - m[1],
//							   s - m[2],
//							   s - m[3]);
//	}
//	
//	template<typename T>
//	DLL_APICALL	mat_t<4, 4, T>	DLL_CALLING_CONVENTION operator-(mat_t<4, 4, T> const& m1, mat_t<4, 4, T> const& m2)
//	{
//		return mat_t<4, 4, T>(
//							   m1[0] - m2[0],
//							   m1[1] - m2[1],
//							   m1[2] - m2[2],
//							   m1[3] - m2[3]);
//	}
//	
//	template<typename T>
//	DLL_APICALL	mat_t<4, 4, T>	DLL_CALLING_CONVENTION operator*(mat_t<4, 4, T> const& m, T const  & s)
//	{
//		return mat_t<4, 4, T>(
//							   m[0] * s,
//							   m[1] * s,
//							   m[2] * s,
//							   m[3] * s);
//	}
//	
//	template<typename T>
//	DLL_APICALL	mat_t<4, 4, T>	DLL_CALLING_CONVENTION operator*(T const& s, mat_t<4, 4, T> const& m)
//	{
//		return mat_t<4, 4, T>(
//							   m[0] * s,
//							   m[1] * s,
//							   m[2] * s,
//							   m[3] * s);
//	}
//	
//	template<typename T>
//	SIMD_FUNC_QUALIFIER typename mat_t<4, 4, T>::col_type operator*
//	(
//	 mat_t<4, 4, T> const& m,
//	 typename mat_t<4, 4, T>::row_type const& v
//	 )
//	{
//		/*
//		 __m128 v0 = _mm_shuffle_ps(v.data, v.data, _MM_SHUFFLE(0, 0, 0, 0));
//		 __m128 v1 = _mm_shuffle_ps(v.data, v.data, _MM_SHUFFLE(1, 1, 1, 1));
//		 __m128 v2 = _mm_shuffle_ps(v.data, v.data, _MM_SHUFFLE(2, 2, 2, 2));
//		 __m128 v3 = _mm_shuffle_ps(v.data, v.data, _MM_SHUFFLE(3, 3, 3, 3));
//		 
//		 __m128 m0 = _mm_mul_ps(m[0].data, v0);
//		 __m128 m1 = _mm_mul_ps(m[1].data, v1);
//		 __m128 a0 = _mm_add_ps(m0, m1);
//		 
//		 __m128 m2 = _mm_mul_ps(m[2].data, v2);
//		 __m128 m3 = _mm_mul_ps(m[3].data, v3);
//		 __m128 a1 = _mm_add_ps(m2, m3);
//		 
//		 __m128 a2 = _mm_add_ps(a0, a1);
//		 
//		 return typename mat_t<4, 4, T>::col_type(a2);
//		 */
//		
//		typename mat_t<4, 4, T>::col_type const Mov0(v[0]);
//		typename mat_t<4, 4, T>::col_type const Mov1(v[1]);
//		typename mat_t<4, 4, T>::col_type const Mul0 = m[0] * Mov0;
//		typename mat_t<4, 4, T>::col_type const Mul1 = m[1] * Mov1;
//		typename mat_t<4, 4, T>::col_type const Add0 = Mul0 + Mul1;
//		typename mat_t<4, 4, T>::col_type const Mov2(v[2]);
//		typename mat_t<4, 4, T>::col_type const Mov3(v[3]);
//		typename mat_t<4, 4, T>::col_type const Mul2 = m[2] * Mov2;
//		typename mat_t<4, 4, T>::col_type const Mul3 = m[3] * Mov3;
//		typename mat_t<4, 4, T>::col_type const Add1 = Mul2 + Mul3;
//		typename mat_t<4, 4, T>::col_type const Add2 = Add0 + Add1;
//		return Add2;
//		
//		/*
//		 return typename mat_t<4, 4, T>::col_type(
//		 m[0][0] * v[0] + m[1][0] * v[1] + m[2][0] * v[2] + m[3][0] * v[3],
//		 m[0][1] * v[0] + m[1][1] * v[1] + m[2][1] * v[2] + m[3][1] * v[3],
//		 m[0][2] * v[0] + m[1][2] * v[1] + m[2][2] * v[2] + m[3][2] * v[3],
//		 m[0][3] * v[0] + m[1][3] * v[1] + m[2][3] * v[2] + m[3][3] * v[3]);
//		 */
//	}
//	
//	template<typename T>
//	SIMD_FUNC_QUALIFIER typename mat_t<4, 4, T>::row_type operator*
//	(
//	 typename mat_t<4, 4, T>::col_type const& v,
//	 mat_t<4, 4, T> const& m
//	 )
//	{
//		return typename mat_t<4, 4, T>::row_type(
//												  m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2] + m[0][3] * v[3],
//												  m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2] + m[1][3] * v[3],
//												  m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2] + m[2][3] * v[3],
//												  m[3][0] * v[0] + m[3][1] * v[1] + m[3][2] * v[2] + m[3][3] * v[3]);
//	}
//	
//	template<typename T>
//	SIMD_FUNC_QUALIFIER mat_t<2, 4, T> operator*(mat_t<4, 4, T> const& m1, mat_t<2, 4, T> const& m2)
//	{
//		return mat_t<2, 4, T>(
//							   m1[0][0] * m2[0][0] + m1[1][0] * m2[0][1] + m1[2][0] * m2[0][2] + m1[3][0] * m2[0][3],
//							   m1[0][1] * m2[0][0] + m1[1][1] * m2[0][1] + m1[2][1] * m2[0][2] + m1[3][1] * m2[0][3],
//							   m1[0][2] * m2[0][0] + m1[1][2] * m2[0][1] + m1[2][2] * m2[0][2] + m1[3][2] * m2[0][3],
//							   m1[0][3] * m2[0][0] + m1[1][3] * m2[0][1] + m1[2][3] * m2[0][2] + m1[3][3] * m2[0][3],
//							   m1[0][0] * m2[1][0] + m1[1][0] * m2[1][1] + m1[2][0] * m2[1][2] + m1[3][0] * m2[1][3],
//							   m1[0][1] * m2[1][0] + m1[1][1] * m2[1][1] + m1[2][1] * m2[1][2] + m1[3][1] * m2[1][3],
//							   m1[0][2] * m2[1][0] + m1[1][2] * m2[1][1] + m1[2][2] * m2[1][2] + m1[3][2] * m2[1][3],
//							   m1[0][3] * m2[1][0] + m1[1][3] * m2[1][1] + m1[2][3] * m2[1][2] + m1[3][3] * m2[1][3]);
//	}
//	
//	template<typename T>
//	SIMD_FUNC_QUALIFIER mat_t<3, 4, T> operator*(mat_t<4, 4, T> const& m1, mat_t<3, 4, T> const& m2)
//	{
//		return mat_t<3, 4, T>(
//							   m1[0][0] * m2[0][0] + m1[1][0] * m2[0][1] + m1[2][0] * m2[0][2] + m1[3][0] * m2[0][3],
//							   m1[0][1] * m2[0][0] + m1[1][1] * m2[0][1] + m1[2][1] * m2[0][2] + m1[3][1] * m2[0][3],
//							   m1[0][2] * m2[0][0] + m1[1][2] * m2[0][1] + m1[2][2] * m2[0][2] + m1[3][2] * m2[0][3],
//							   m1[0][3] * m2[0][0] + m1[1][3] * m2[0][1] + m1[2][3] * m2[0][2] + m1[3][3] * m2[0][3],
//							   m1[0][0] * m2[1][0] + m1[1][0] * m2[1][1] + m1[2][0] * m2[1][2] + m1[3][0] * m2[1][3],
//							   m1[0][1] * m2[1][0] + m1[1][1] * m2[1][1] + m1[2][1] * m2[1][2] + m1[3][1] * m2[1][3],
//							   m1[0][2] * m2[1][0] + m1[1][2] * m2[1][1] + m1[2][2] * m2[1][2] + m1[3][2] * m2[1][3],
//							   m1[0][3] * m2[1][0] + m1[1][3] * m2[1][1] + m1[2][3] * m2[1][2] + m1[3][3] * m2[1][3],
//							   m1[0][0] * m2[2][0] + m1[1][0] * m2[2][1] + m1[2][0] * m2[2][2] + m1[3][0] * m2[2][3],
//							   m1[0][1] * m2[2][0] + m1[1][1] * m2[2][1] + m1[2][1] * m2[2][2] + m1[3][1] * m2[2][3],
//							   m1[0][2] * m2[2][0] + m1[1][2] * m2[2][1] + m1[2][2] * m2[2][2] + m1[3][2] * m2[2][3],
//							   m1[0][3] * m2[2][0] + m1[1][3] * m2[2][1] + m1[2][3] * m2[2][2] + m1[3][3] * m2[2][3]);
//	}
//	
//	template<typename T>
//	DLL_APICALL	mat_t<4, 4, T>	DLL_CALLING_CONVENTION operator*(mat_t<4, 4, T> const& m1, mat_t<4, 4, T> const& m2)
//	{
//		typename mat_t<4, 4, T>::col_type const SrcA0 = m1[0];
//		typename mat_t<4, 4, T>::col_type const SrcA1 = m1[1];
//		typename mat_t<4, 4, T>::col_type const SrcA2 = m1[2];
//		typename mat_t<4, 4, T>::col_type const SrcA3 = m1[3];
//		
//		typename mat_t<4, 4, T>::col_type const SrcB0 = m2[0];
//		typename mat_t<4, 4, T>::col_type const SrcB1 = m2[1];
//		typename mat_t<4, 4, T>::col_type const SrcB2 = m2[2];
//		typename mat_t<4, 4, T>::col_type const SrcB3 = m2[3];
//		
//		mat_t<4, 4, T> Result;
//		Result[0] = SrcA0 * SrcB0[0] + SrcA1 * SrcB0[1] + SrcA2 * SrcB0[2] + SrcA3 * SrcB0[3];
//		Result[1] = SrcA0 * SrcB1[0] + SrcA1 * SrcB1[1] + SrcA2 * SrcB1[2] + SrcA3 * SrcB1[3];
//		Result[2] = SrcA0 * SrcB2[0] + SrcA1 * SrcB2[1] + SrcA2 * SrcB2[2] + SrcA3 * SrcB2[3];
//		Result[3] = SrcA0 * SrcB3[0] + SrcA1 * SrcB3[1] + SrcA2 * SrcB3[2] + SrcA3 * SrcB3[3];
//		return Result;
//	}
//
//	template<typename T>
//	DLL_APICALL	mat_t<4, 4, T>	DLL_CALLING_CONVENTION operator/(mat_t<4, 4, T> const& m, T const& s)
//	{
//		return mat_t<4, 4, T>(
//							   m[0] / s,
//							   m[1] / s,
//							   m[2] / s,
//							   m[3] / s);
//	}
//	
//	template<typename T>
//	DLL_APICALL	mat_t<4, 4, T>	DLL_CALLING_CONVENTION operator/(T const& s,	mat_t<4, 4, T> const& m)
//	{
//		return mat_t<4, 4, T>(
//							   s / m[0],
//							   s / m[1],
//							   s / m[2],
//							   s / m[3]);
//	}
//	
//	template<typename T>
//	SIMD_FUNC_QUALIFIER typename mat_t<4, 4, T>::col_type operator/(mat_t<4, 4, T> const& m, typename mat_t<4, 4, T>::row_type const& v)
//	{
//		return inverse(m) * v;
//	}
//	
//	template<typename T>
//	SIMD_FUNC_QUALIFIER typename mat_t<4, 4, T>::row_type operator/(typename mat_t<4, 4, T>::col_type const& v, mat_t<4, 4, T> const& m)
//	{
//		return v * inverse(m);
//	}
//	
//	template<typename T>
//	DLL_APICALL	mat_t<4, 4, T>	DLL_CALLING_CONVENTION operator/(mat_t<4, 4, T> const& m1, mat_t<4, 4, T> const& m2)
//	{
//		mat_t<4, 4, T> m1_copy(m1);
//		return m1_copy /= m2;
//	}
//	
//	// -- Boolean operators --
//	
//	template<typename T>
//	SIMD_FUNC_QUALIFIER bool operator==(mat_t<4, 4, T> const& m1, mat_t<4, 4, T> const& m2)
//	{
//		return (m1[0] == m2[0]) && (m1[1] == m2[1]) && (m1[2] == m2[2]) && (m1[3] == m2[3]);
//	}
//	
//	template<typename T>
//	SIMD_FUNC_QUALIFIER bool operator!=(mat_t<4, 4, T> const& m1, mat_t<4, 4, T> const& m2)
//	{
//		return (m1[0] != m2[0]) || (m1[1] != m2[1]) || (m1[2] != m2[2]) || (m1[3] != m2[3]);
//	}
}//namespace glm

















////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// operators
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace simd
{
#define MATRIX_PERFORM_PARTIAL_OPERATOR(_O, _F) \
	template<int C, int R, typename T> \
	SIMD_STATIC_INLINE	mat_t<C, R, T>	DLL_CALLING_CONVENTION perform_operator_##_F(const mat_t<C, R, T>& a, const mat_t<C, R, T>& b) \
	{ \
		mat_t<C, R, T> ret; \
		for (int i = 0; i < C; i++) \
			ret[i] = a[i] _O b[i]; \
		return ret; \
	}

#define MATRIX_PERFORM_PARTIAL_OPERATOR_ESCALAR(_O, _F) \
	template<int C, int R, typename T> \
	SIMD_STATIC_INLINE	mat_t<C, R, T>	DLL_CALLING_CONVENTION perform_operator_##_F(const mat_t<C, R, T>& a, T& b) \
	{ \
		mat_t<C, R, T> ret; \
		for (int i = 0; i < C; i++) \
			ret[i] = a[i] _O b[i]; \
		return ret; \
	} \
	template<int C, int R, typename T> \
	SIMD_STATIC_INLINE	mat_t<C, R, T>	DLL_CALLING_CONVENTION perform_operator_##_F(const T& a, const mat_t<C, R, T>& b) \
	{ \
		mat_t<C, R, T> ret; \
		for (int i = 0; i < C; i++) \
			ret[i] = a[i] _O b[i]; \
		return ret; \
	}
	
	template<int C, int R, typename T>
	SIMD_STATIC_INLINE	mat_t<C, R, T>	DLL_CALLING_CONVENTION perform_operator_cross(const mat_t<C, R, T>& m1, const mat_t<C, R, T>& m2)
	{
		mat_t<C, R, T> Result;
		assert(0);
		return Result;
	}

	template<typename T>
	SIMD_STATIC_INLINE	mat_t<4, 4, T>	DLL_CALLING_CONVENTION perform_operator_cross(const mat_t<4, 4, T>& m1, const mat_t<4, 4, T>& m2)
	{
		typename mat_t<4, 4, T>::col_type const SrcA0 = m1[0];
		typename mat_t<4, 4, T>::col_type const SrcA1 = m1[1];
		typename mat_t<4, 4, T>::col_type const SrcA2 = m1[2];
		typename mat_t<4, 4, T>::col_type const SrcA3 = m1[3];
		
		typename mat_t<4, 4, T>::col_type const SrcB0 = m2[0];
		typename mat_t<4, 4, T>::col_type const SrcB1 = m2[1];
		typename mat_t<4, 4, T>::col_type const SrcB2 = m2[2];
		typename mat_t<4, 4, T>::col_type const SrcB3 = m2[3];
		
		mat_t<4, 4, T> Result;
		Result[0] = SrcA0 * SrcB0[0] + SrcA1 * SrcB0[1] + SrcA2 * SrcB0[2] + SrcA3 * SrcB0[3];
		Result[1] = SrcA0 * SrcB1[0] + SrcA1 * SrcB1[1] + SrcA2 * SrcB1[2] + SrcA3 * SrcB1[3];
		Result[2] = SrcA0 * SrcB2[0] + SrcA1 * SrcB2[1] + SrcA2 * SrcB2[2] + SrcA3 * SrcB2[3];
		Result[3] = SrcA0 * SrcB3[0] + SrcA1 * SrcB3[1] + SrcA2 * SrcB3[2] + SrcA3 * SrcB3[3];
		return Result;
	}

	
#if CPU_ARCH & CPU_ARCH_SSE2_BIT
	
	SIMD_STATIC_INLINE	mat_t<4, 4, float>	DLL_CALLING_CONVENTION perform_operator_add(const mat_t<4, 4, float>& m1, const mat_t<4, 4, float>& m2)
	{
		mat_t<4, 4, float> ret;
		SIMD_mat4_add((const SIMD_vec4*)m1.internal_vector, (const SIMD_vec4*)m2.internal_vector, (SIMD_vec4*)ret.internal_vector);
		return ret;
	}
	
	
	SIMD_STATIC_INLINE	mat_t<4, 4, float>	DLL_CALLING_CONVENTION perform_operator_sub(const mat_t<4, 4, float>& m1, const mat_t<4, 4, float>& m2)
	{
		mat_t<4, 4, float> ret;
		SIMD_mat4_sub((const SIMD_vec4*)m1.internal_vector, (const SIMD_vec4*)m2.internal_vector, (SIMD_vec4*)ret.internal_vector);
		return ret;
	}
	
	
	SIMD_STATIC_INLINE	mat_t<4, 4, float>	DLL_CALLING_CONVENTION perform_operator_cross(const mat_t<4, 4, float>& m1, const mat_t<4, 4, float>& m2)
	{
		mat_t<4, 4, float> ret;
		SIMD_mat4_mul((const SIMD_vec4*)m1.internal_vector, (const SIMD_vec4*)m2.internal_vector, (SIMD_vec4*)ret.internal_vector);
		return ret;
	}
#endif

	MATRIX_PERFORM_PARTIAL_OPERATOR(+, add)
	MATRIX_PERFORM_PARTIAL_OPERATOR(-, sub)

	MATRIX_PERFORM_PARTIAL_OPERATOR_ESCALAR(+, add)
	MATRIX_PERFORM_PARTIAL_OPERATOR_ESCALAR(-, sub)
	MATRIX_PERFORM_PARTIAL_OPERATOR_ESCALAR(*, mul)
	MATRIX_PERFORM_PARTIAL_OPERATOR_ESCALAR(/, div)
	
	
	MATRIX_M_M_M_O_F(+, add)
	MATRIX_M_M_M_O_F(-, sub)
	MATRIX_M_M_M_O_F(*, cross)

//	MATRIX_M_M_S_O_F(+, add)
//	MATRIX_M_M_S_O_F(-, sub)
//	MATRIX_M_M_S_O_F(*, mul)
//	MATRIX_M_M_S_O_F(/, div)
}














////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// transpose
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace simd
{
	template<int C, int R, typename T>
	SIMD_STATIC_INLINE	mat_t<C, R, T>	DLL_CALLING_CONVENTION perform_transpose(const mat_t<C, R, T>& m)
	{
		mat_t<R, C, T> Result;
		for (int c = 0; c < C; c++)
		{
			for (int r = 0; r < R; r++)
			{
				Result[r][c] = m[c][r];
			}
		}
		return Result;
	}
	
#if CPU_ARCH & CPU_ARCH_SSE2_BIT
	
	SIMD_STATIC_INLINE	mat_t<4, 4, float>	DLL_CALLING_CONVENTION perform_transpose(const mat_t<4, 4, float>& m)
	{
		mat_t<4, 4, float> Result;
		SIMD_mat4_transpose((const SIMD_vec4*)&m.internal_vector[0], (SIMD_vec4*)&Result.internal_vector[0]);
		return Result;
	}
#endif
	
	MATRIX_M_M_F(transpose)
}








////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// mulComponents
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace simd
{
	template<int C, int R, typename T1, typename T2>
	SIMD_STATIC_INLINE	mat_t<C, R, T1>	DLL_CALLING_CONVENTION perform_mulComponents(const mat_t<C, R, T1>& a, const mat_t<C, R, T2>& b)
	{
		mat_t<R, C, T1> Result;
		for (int c = 0; c < C; c++)
			Result[c] = a[c] * b[c];
		return Result;
	}
	
#if CPU_ARCH & CPU_ARCH_SSE2_BIT
	SIMD_STATIC_INLINE	mat_t<4, 4, float>	DLL_CALLING_CONVENTION perform_mulComponents(const mat_t<4, 4, float>& a, const mat_t<4, 4, float>& b)
	{
		mat_t<4, 4, float> Result;
		SIMD_mat4_matrixCompMult((const SIMD_vec4*)&a.internal_vector[0], (const SIMD_vec4*)&b.internal_vector[0], (SIMD_vec4*)&Result.internal_vector[0]);
		return Result;
	}
#endif
	
	//TYPE_M_M_F(mulComponents)
}

















////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// inverse
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace simd
{
	template<int C, int R, typename T>
	SIMD_STATIC_INLINE	mat_t<C, R, T>	DLL_CALLING_CONVENTION perform_inverse(const mat_t<C, R, T>& m)
	{
		mat_t<R, C, T> Result;
		assert(0);
		return Result;
	}
	
	template<int C, int R, typename T>
	SIMD_STATIC_INLINE	mat_t<C, R, T>	DLL_CALLING_CONVENTION perform_inverse_fast(const mat_t<C, R, T>& m)
	{
		mat_t<R, C, T> Result;
		assert(0);
		return Result;
	}
	
	template<int C, int R, typename T>
	SIMD_STATIC_INLINE	mat_t<C, R, T>	DLL_CALLING_CONVENTION perform_inverse_homogeneous(const mat_t<C, R, T>& m)
	{
		mat_t<R, C, T> Result;
		assert(0);
		return Result;
	}

	template<typename T>
	SIMD_STATIC_INLINE	mat_t<4, 4, T>	DLL_CALLING_CONVENTION perform_inverse_homogeneous(const mat_t<4, 4, T>& m)
	{
		mat_t<4, 4, T> Result;
		
		for (int i = 0; i < 3; i++)
		{
			pre::move(m[i][0], Result[0][i]);
			pre::move(m[i][1], Result[1][i]);
			pre::move(m[i][2], Result[2][i]);
			Result[3][i] = -dot3(m[i], m[3]);
		}
		pre::zero(Result[0][3]);
		pre::zero(Result[1][3]);
		pre::zero(Result[2][3]);
		Result[3][3] = T(1);
		return Result;
	}

#if CPU_ARCH & CPU_ARCH_SSE2_BIT
	
	SIMD_STATIC_INLINE	mat_t<4, 4, float>	DLL_CALLING_CONVENTION perform_inverse(const mat_t<4, 4, float>& m)
	{
		mat_t<4, 4, float> Result;
		SIMD_mat4_inverse((const SIMD_vec4*)&m.internal_vector[0], (SIMD_vec4*)&Result.internal_vector[0]);
		return Result;
	}
	
	SIMD_STATIC_INLINE	mat_t<4, 4, float>	DLL_CALLING_CONVENTION perform_inverse_fast(const mat_t<4, 4, float>& m)
	{
		mat_t<4, 4, float> Result;
		SIMD_mat4_inverse_lowp((const SIMD_vec4*)&m.internal_vector[0], (SIMD_vec4*)&Result.internal_vector[0]);
		return Result;
	}
#endif
	
	MATRIX_M_M_F(inverse)
	MATRIX_M_M_F(inverse_fast)
	MATRIX_M_M_F(inverse_homogeneous)
}























namespace simd
{
//
//#	if !SIMD_HAS_DEFAULTED_FUNCTIONS || defined(PRE_FORCE_CTOR_INIT)
//	template<typename T>
//	mat_t<4, 4, T>::mat_t()
//	{
//#			ifdef PRE_FORCE_CTOR_INIT
//		this->value[0] = col_type(1, 0, 0, 0);
//		this->value[1] = col_type(0, 1, 0, 0);
//		this->value[2] = col_type(0, 0, 1, 0);
//		this->value[3] = col_type(0, 0, 0, 1);
//#			endif
//	}
//#	endif
	
//#	if !SIMD_HAS_DEFAULTED_FUNCTIONS
//	template<typename T>
//	mat_t<4, 4, T>::mat_t(mat_t<4, 4, T> const& m)
//	{
//		this->value[0] = m[0];
//		this->value[1] = m[1];
//		this->value[2] = m[2];
//		this->value[3] = m[3];
//	}
//#	endif//!SIMD_HAS_DEFAULTED_FUNCTIONS
//
//	template<typename T>
//	template<qualifier P>
//	mat_t<4, 4, T>::mat_t(mat_t<4, 4, T, P> const& m)
//	{
//		this->value[0] = m[0];
//		this->value[1] = m[1];
//		this->value[2] = m[2];
//		this->value[3] = m[3];
//	}

	template<int C, int R, typename T>
	inline	void mat_t_ctor(mat_t<C, R, T>& m, const T& s)
	{
	}
	
	template<typename T>
	inline	void mat_t_ctor(mat_t<4, 4, T>& m, const T& s)
	{
		m.internal_vector[0] = typename mat_t<4, 4, T>::col_type(s, 0, 0, 0);
		m.internal_vector[1] = typename mat_t<4, 4, T>::col_type(0, s, 0, 0);
		m.internal_vector[2] = typename mat_t<4, 4, T>::col_type(0, 0, s, 0);
		m.internal_vector[3] = typename mat_t<4, 4, T>::col_type(0, 0, 0, s);
	}
	
	template<int C, int R, typename T>
	mat_t<C, R, T>::mat_t(const T& s)
	{
		mat_t_ctor(*this, s);
	}
	

//
//	template<typename T>
//	mat_t<4, 4, T>::mat_t
//	(
//	 T const& x0, T const& y0, T const& z0, T const& w0,
//	 T const& x1, T const& y1, T const& z1, T const& w1,
//	 T const& x2, T const& y2, T const& z2, T const& w2,
//	 T const& x3, T const& y3, T const& z3, T const& w3
//	 )
//	{
//		this->value[0] = col_type(x0, y0, z0, w0);
//		this->value[1] = col_type(x1, y1, z1, w1);
//		this->value[2] = col_type(x2, y2, z2, w2);
//		this->value[3] = col_type(x3, y3, z3, w3);
//	}
//
//	template<typename T>
//	mat_t<4, 4, T>::mat_t
//	(
//	 col_type const& v0,
//	 col_type const& v1,
//	 col_type const& v2,
//	 col_type const& v3
//	 )
//	{
//		this->value[0] = v0;
//		this->value[1] = v1;
//		this->value[2] = v2;
//		this->value[3] = v3;
//	}
//
//	template<typename T>
//	template<typename U, Nualifier P>
//	mat_t<4, 4, T>::mat_t(mat_t<4, 4, U, P> const& m)
//	{
//		this->value[0] = col_type(m[0]);
//		this->value[1] = col_type(m[1]);
//		this->value[2] = col_type(m[2]);
//		this->value[3] = col_type(m[3]);
//	}
//
//	// -- Conversions --
//
//	template<typename T>
//	template<
//	typename X1, typename Y1, typename Z1, typename W1,
//	typename X2, typename Y2, typename Z2, typename W2,
//	typename X3, typename Y3, typename Z3, typename W3,
//	typename X4, typename Y4, typename Z4, typename W4>
//	mat_t<4, 4, T>::mat_t
//	(
//	 X1 const& x1, Y1 const& y1, Z1 const& z1, W1 const& w1,
//	 X2 const& x2, Y2 const& y2, Z2 const& z2, W2 const& w2,
//	 X3 const& x3, Y3 const& y3, Z3 const& z3, W3 const& w3,
//	 X4 const& x4, Y4 const& y4, Z4 const& z4, W4 const& w4
//	 )
//	{
////		SIMD_STATIC_ASSERT(std::numeric_limits<X1>::is_iec559 || std::numeric_limits<X1>::is_integer || SIMD_UNRESTRICTED_GENTYPE, "*mat4x4 constructor only takes float and integer types, 1st parameter type invalid.");
////		SIMD_STATIC_ASSERT(std::numeric_limits<Y1>::is_iec559 || std::numeric_limits<Y1>::is_integer || SIMD_UNRESTRICTED_GENTYPE, "*mat4x4 constructor only takes float and integer types, 2nd parameter type invalid.");
////		SIMD_STATIC_ASSERT(std::numeric_limits<Z1>::is_iec559 || std::numeric_limits<Z1>::is_integer || SIMD_UNRESTRICTED_GENTYPE, "*mat4x4 constructor only takes float and integer types, 3rd parameter type invalid.");
////		SIMD_STATIC_ASSERT(std::numeric_limits<W1>::is_iec559 || std::numeric_limits<W1>::is_integer || SIMD_UNRESTRICTED_GENTYPE, "*mat4x4 constructor only takes float and integer types, 4th parameter type invalid.");
////
////		SIMD_STATIC_ASSERT(std::numeric_limits<X2>::is_iec559 || std::numeric_limits<X2>::is_integer || SIMD_UNRESTRICTED_GENTYPE, "*mat4x4 constructor only takes float and integer types, 5th parameter type invalid.");
////		SIMD_STATIC_ASSERT(std::numeric_limits<Y2>::is_iec559 || std::numeric_limits<Y2>::is_integer || SIMD_UNRESTRICTED_GENTYPE, "*mat4x4 constructor only takes float and integer types, 6th parameter type invalid.");
////		SIMD_STATIC_ASSERT(std::numeric_limits<Z2>::is_iec559 || std::numeric_limits<Z2>::is_integer || SIMD_UNRESTRICTED_GENTYPE, "*mat4x4 constructor only takes float and integer types, 7th parameter type invalid.");
////		SIMD_STATIC_ASSERT(std::numeric_limits<W2>::is_iec559 || std::numeric_limits<W2>::is_integer || SIMD_UNRESTRICTED_GENTYPE, "*mat4x4 constructor only takes float and integer types, 8th parameter type invalid.");
////
////		SIMD_STATIC_ASSERT(std::numeric_limits<X3>::is_iec559 || std::numeric_limits<X3>::is_integer || SIMD_UNRESTRICTED_GENTYPE, "*mat4x4 constructor only takes float and integer types, 9th parameter type invalid.");
////		SIMD_STATIC_ASSERT(std::numeric_limits<Y3>::is_iec559 || std::numeric_limits<Y3>::is_integer || SIMD_UNRESTRICTED_GENTYPE, "*mat4x4 constructor only takes float and integer types, 10th parameter type invalid.");
////		SIMD_STATIC_ASSERT(std::numeric_limits<Z3>::is_iec559 || std::numeric_limits<Z3>::is_integer || SIMD_UNRESTRICTED_GENTYPE, "*mat4x4 constructor only takes float and integer types, 11th parameter type invalid.");
////		SIMD_STATIC_ASSERT(std::numeric_limits<W3>::is_iec559 || std::numeric_limits<W3>::is_integer || SIMD_UNRESTRICTED_GENTYPE, "*mat4x4 constructor only takes float and integer types, 12th parameter type invalid.");
////
////		SIMD_STATIC_ASSERT(std::numeric_limits<X4>::is_iec559 || std::numeric_limits<X4>::is_integer || SIMD_UNRESTRICTED_GENTYPE, "*mat4x4 constructor only takes float and integer types, 13th parameter type invalid.");
////		SIMD_STATIC_ASSERT(std::numeric_limits<Y4>::is_iec559 || std::numeric_limits<Y4>::is_integer || SIMD_UNRESTRICTED_GENTYPE, "*mat4x4 constructor only takes float and integer types, 14th parameter type invalid.");
////		SIMD_STATIC_ASSERT(std::numeric_limits<Z4>::is_iec559 || std::numeric_limits<Z4>::is_integer || SIMD_UNRESTRICTED_GENTYPE, "*mat4x4 constructor only takes float and integer types, 15th parameter type invalid.");
////		SIMD_STATIC_ASSERT(std::numeric_limits<W4>::is_iec559 || std::numeric_limits<W4>::is_integer || SIMD_UNRESTRICTED_GENTYPE, "*mat4x4 constructor only takes float and integer types, 16th parameter type invalid.");
//
//		this->value[0] = col_type(static_cast<T>(x1), value_type(y1), value_type(z1), value_type(w1));
//		this->value[1] = col_type(static_cast<T>(x2), value_type(y2), value_type(z2), value_type(w2));
//		this->value[2] = col_type(static_cast<T>(x3), value_type(y3), value_type(z3), value_type(w3));
//		this->value[3] = col_type(static_cast<T>(x4), value_type(y4), value_type(z4), value_type(w4));
//	}
//
//	template<typename T>
//	template<typename V1, typename V2, typename V3, typename V4>
//	mat_t<4, 4, T>::mat_t
//	(
//	 pack_t<4, V1> const& v1,
//	 pack_t<4, V2> const& v2,
//	 pack_t<4, V3> const& v3,
//	 pack_t<4, V4> const& v4
//	 )
//	{
////		SIMD_STATIC_ASSERT(std::numeric_limits<V1>::is_iec559 || std::numeric_limits<V1>::is_integer || SIMD_UNRESTRICTED_GENTYPE, "*mat4x4 constructor only takes float and integer types, 1st parameter type invalid.");
////		SIMD_STATIC_ASSERT(std::numeric_limits<V2>::is_iec559 || std::numeric_limits<V2>::is_integer || SIMD_UNRESTRICTED_GENTYPE, "*mat4x4 constructor only takes float and integer types, 2nd parameter type invalid.");
////		SIMD_STATIC_ASSERT(std::numeric_limits<V3>::is_iec559 || std::numeric_limits<V3>::is_integer || SIMD_UNRESTRICTED_GENTYPE, "*mat4x4 constructor only takes float and integer types, 3rd parameter type invalid.");
////		SIMD_STATIC_ASSERT(std::numeric_limits<V4>::is_iec559 || std::numeric_limits<V4>::is_integer || SIMD_UNRESTRICTED_GENTYPE, "*mat4x4 constructor only takes float and integer types, 4th parameter type invalid.");
//
//		this->value[0] = col_type(v1);
//		this->value[1] = col_type(v2);
//		this->value[2] = col_type(v3);
//		this->value[3] = col_type(v4);
//	}
}

namespace simd
{
//	template<int C, int R, typename T>
//	mat_t<C, R, T>& mat_t<C, R, T>::operator = (const mat_t<C, R, T>& m)
//	{
//		for (int i = 0; i < C; i++)
//			column[i] = m.column[i];
//		return *this;
//	}

	template<int C, int R, typename T>
	auto	mat_t<C, R, T>::row(int Index) const -> mat_t::row_type
	{
		assert(0 <= Index && Index < R);
		
		row_type Result;
		
		for(int i = 0; i < C; i++)
			Result[i] = internal_vector[i][Index];
			return Result;
	}

}


namespace simd
{
#define MATRIX_DEFINITIONS_C_R_T(_C, _R, _T) \
	template	mat_t<_C, _R, _T>::mat_t(const _T&);\
	template	mat_t<_C, _R, _T>::row_type	mat_t<_C, _R, _T>::row(int) const;\

//	template	mat_t<_C, _R, _T, _N>&	mat_t<_C, _R, _T, _N>::operator = (const mat_t<_C, _R, _T, _N>&);

#define MATRIX_DEFINITIONS_C_R(_C, _R) \
	MATRIX_DEFINITIONS_C_R_T(_C, _R, int8_t) \
	MATRIX_DEFINITIONS_C_R_T(_C, _R, uint8_t) \
	MATRIX_DEFINITIONS_C_R_T(_C, _R, int16_t) \
	MATRIX_DEFINITIONS_C_R_T(_C, _R, uint16_t) \
	MATRIX_DEFINITIONS_C_R_T(_C, _R, int32_t) \
	MATRIX_DEFINITIONS_C_R_T(_C, _R, uint32_t) \
	MATRIX_DEFINITIONS_C_R_T(_C, _R, int64_t) \
	MATRIX_DEFINITIONS_C_R_T(_C, _R, uint64_t) \
	MATRIX_DEFINITIONS_C_R_T(_C, _R, float) \
	MATRIX_DEFINITIONS_C_R_T(_C, _R, double) \

	
	MATRIX_DEFINITIONS_C_R(4, 4)
}



























namespace simd
{
	template<int C, int R, typename T>
	FORCE_INLINE	auto	perform_operator_mul(const mat_t<C, R, T>& m, const pack_t<C, T>& v) -> pack_t<R, T>
	{
		assert(0);
	}

	template<typename T>
	FORCE_INLINE	auto	perform_operator_mul(const mat_t<4, 4, T>& m, const pack_t<4, T>& v) -> pack_t<4, T>
	{
		pack_t<4, T> const Mov0(v[0]);
		pack_t<4, T> const Mov1(v[1]);
		pack_t<4, T> const Mul0 = m[0] * Mov0;
		pack_t<4, T> const Mul1 = m[1] * Mov1;
		pack_t<4, T> const Add0 = Mul0 + Mul1;
		pack_t<4, T> const Mov2(v[2]);
		pack_t<4, T> const Mov3(v[3]);
		pack_t<4, T> const Mul2 = m[2] * Mov2;
		pack_t<4, T> const Mul3 = m[3] * Mov3;
		pack_t<4, T> const Add1 = Mul2 + Mul3;
		pack_t<4, T> const Add2 = Add0 + Add1;
		return Add2;
	}
	
#if CPU_ARCH & CPU_ARCH_SSE2_BIT
	FORCE_INLINE	auto	perform_operator_mul(const mat_t<4, 4, float>& m, const pack_t<4, float>& v) -> pack_t<4, float>
	{
		return SIMD_mat4_mul_vec4((const SIMD_vec4*)&m, v.m);
	}
#endif

	template<int C, int R, typename T>
	DLL_FUNCTION(auto)	operator *(const mat_t<C, R, T>& m, const pack_t<C, T>& v) -> pack_t<R, T>
	{
		return perform_operator_mul(m, v);
	}
	
#define MAT_VEC_MUL_C_R_T(_C, _R, _T) \
	template	DLL_FUNCTION(auto)	operator * <_C, _R, _T>(const mat_t<_C, _R, _T>& m, const pack_t<_C, _T>& v) -> pack_t<_R, _T>;
	
	MAT_VEC_MUL_C_R_T(4, 4, float)

}

