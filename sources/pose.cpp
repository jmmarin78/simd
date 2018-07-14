// This software is MIT licensed (see LICENSE)

#define BUILD_DLL

#include <math.h>

#include <pre/number.h>
#include <simd/shapes.h>

#include "defines.h"


namespace simd
{
	
	template<typename T1, typename T2>
	mat_t<4, 4, T1> lookAtRH(pack_t<3, T1> const& eye, pack_t<3, T1> const& center, pack_t<3, T2> const& up)
	{
		auto	f(normalize(center - eye));
		auto	s(normalize(cross(f, up)));
		auto	u(cross(s, f));
		
		mat_t<4, 4, T1> Result(1);
		Result[0][0] = s.x;
		Result[1][0] = s.y;
		Result[2][0] = s.z;
		Result[0][1] = u.x;
		Result[1][1] = u.y;
		Result[2][1] = u.z;
		Result[0][2] =-f.x;
		Result[1][2] =-f.y;
		Result[2][2] =-f.z;
		Result[3][0] =-dot(s, eye);
		Result[3][1] =-dot(u, eye);
		Result[3][2] = dot(f, eye);
		return Result;
	}
	
	template<typename T1, typename T2>
	mat_t<4, 4, T1> lookAtLH(pack_t<3, T1> const& eye, pack_t<3, T1> const& center, pack_t<3, T2> const& up)
	{
		auto	f(normalize(center - eye));
		auto	s(normalize(cross(up, f)));
		auto	u(cross(f, s));
		
		mat_t<4, 4, T1> Result(1);
		Result[0][0] = s.x;
		Result[1][0] = s.y;
		Result[2][0] = s.z;
		Result[0][1] = u.x;
		Result[1][1] = u.y;
		Result[2][1] = u.z;
		Result[0][2] = f.x;
		Result[1][2] = f.y;
		Result[2][2] = f.z;
		Result[3][0] = -dot(s, eye);
		Result[3][1] = -dot(u, eye);
		Result[3][2] = -dot(f, eye);
		return Result;
	}
	
	template<typename T1, typename T2>
	mat_t<4, 4, T1> lookAt(pack_t<3, T1> const& eye, pack_t<3, T1> const& center, pack_t<3, T2> const& up)
	{
#		if SIMD_COORDINATE_SYSTEM == SIMD_LEFT_HANDED
		return lookAtLH(eye, center, up);
#		else
		return lookAtRH(eye, center, up);
#		endif
	}
	
	template<typename T1, typename T2>
	mat_t<4, 4, T1> lookAtDirect(pack_t<3, T1> const& eye, pack_t<3, T1> const& center, pack_t<3, T2> const& up)
	{
#		if SIMD_COORDINATE_SYSTEM == SIMD_LEFT_HANDED
		auto	f(normalize(center - eye));
		auto	s(normalize(cross(f, up)));
		auto	u(cross(s, f));
#		else
		auto	f(normalize(eye - center));
		auto	s(normalize(cross(up, f)));
		auto	u(cross(f, s));
#		endif

		mat_t<4, 4, T1> Result;
		Result[0][0] = s.x;
		Result[0][1] = s.y;
		Result[0][2] = s.z;
		pre::zero(Result[0][3]);
		Result[1][0] = u.x;
		Result[1][1] = u.y;
		Result[1][2] = u.z;
		pre::zero(Result[1][3]);
		Result[2][0] = f.x;
		Result[2][1] = f.y;
		Result[2][2] = f.z;
		pre::zero(Result[2][3]);
		Result[3][0] = eye.x;
		Result[3][1] = eye.y;
		Result[3][2] = eye.z;
		Result[3][3] = 1;
		return Result;
	}

}





namespace simd
{
//	mat_t<4, 4, float> lookAt<float>(pack_t<3, float> const&, pack_t<3, float> const&, pack_t<3, float> const&);

#define MAT4_DEFINITIONS_T(_T) \
template	mat_t<4, 4, _T> lookAt<_T, _T>(pack_t<3, _T> const& eye, pack_t<3, _T> const& center, pack_t<3, _T> const& up);\
template	mat_t<4, 4, _T> lookAtDirect<_T, _T>(pack_t<3, _T> const& eye, pack_t<3, _T> const& center, pack_t<3, _T> const& up);\

#define MAT4_DEFINITIONS_T2(_T1, _T2) \
template	mat_t<4, 4, _T1> lookAt<_T1, _T2>(pack_t<3, _T1> const& eye, pack_t<3, _T1> const& center, pack_t<3, _T2> const& up);\
template	mat_t<4, 4, _T1> lookAtDirect<_T1, _T2>(pack_t<3, _T1> const& eye, pack_t<3, _T1> const& center, pack_t<3, _T2> const& up);\


	MAT4_DEFINITIONS_T(float)
	MAT4_DEFINITIONS_T(double)
	
	MAT4_DEFINITIONS_T2(double, float)

}


