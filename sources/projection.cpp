// This software is MIT licensed (see LICENSE)

#define BUILD_DLL

#include <math.h>

#include <simd/shapes.h>

#include "defines.h"


namespace simd
{
	
	template<typename T>
	mat_t<4, 4, T> orthoMatrix(T left, T right, T bottom, T top)
	{
		mat_t<4, 4, T> Result(static_cast<T>(1));
		Result[0][0] = static_cast<T>(2) / (right - left);
		Result[1][1] = static_cast<T>(2) / (top - bottom);
		Result[3][0] = - (right + left) / (right - left);
		Result[3][1] = - (top + bottom) / (top - bottom);
		return Result;
	}
	
	template<typename T>
	mat_t<4, 4, T> orthoLH_ZO(T left, T right, T bottom, T top, T zNear, T zFar)
	{
		mat_t<4, 4, T> Result(1);
		Result[0][0] = static_cast<T>(2) / (right - left);
		Result[1][1] = static_cast<T>(2) / (top - bottom);
		Result[2][2] = static_cast<T>(1) / (zFar - zNear);
		Result[3][0] = - (right + left) / (right - left);
		Result[3][1] = - (top + bottom) / (top - bottom);
		Result[3][2] = - zNear / (zFar - zNear);
		return Result;
	}
	
	template<typename T>
	mat_t<4, 4, T> orthoLH_NO(T left, T right, T bottom, T top, T zNear, T zFar)
	{
		mat_t<4, 4, T> Result(1);
		Result[0][0] = static_cast<T>(2) / (right - left);
		Result[1][1] = static_cast<T>(2) / (top - bottom);
		Result[2][2] = static_cast<T>(2) / (zFar - zNear);
		Result[3][0] = - (right + left) / (right - left);
		Result[3][1] = - (top + bottom) / (top - bottom);
		Result[3][2] = - (zFar + zNear) / (zFar - zNear);
		return Result;
	}
	
	template<typename T>
	mat_t<4, 4, T> orthoRH_ZO(T left, T right, T bottom, T top, T zNear, T zFar)
	{
		mat_t<4, 4, T> Result(1);
		Result[0][0] = static_cast<T>(2) / (right - left);
		Result[1][1] = static_cast<T>(2) / (top - bottom);
		Result[2][2] = - static_cast<T>(1) / (zFar - zNear);
		Result[3][0] = - (right + left) / (right - left);
		Result[3][1] = - (top + bottom) / (top - bottom);
		Result[3][2] = - zNear / (zFar - zNear);
		return Result;
	}
	
	template<typename T>
	mat_t<4, 4, T> orthoRH_NO(T left, T right, T bottom, T top, T zNear, T zFar)
	{
		mat_t<4, 4, T> Result(1);
		Result[0][0] = static_cast<T>(2) / (right - left);
		Result[1][1] = static_cast<T>(2) / (top - bottom);
		Result[2][2] = - static_cast<T>(2) / (zFar - zNear);
		Result[3][0] = - (right + left) / (right - left);
		Result[3][1] = - (top + bottom) / (top - bottom);
		Result[3][2] = - (zFar + zNear) / (zFar - zNear);
		return Result;
	}
	
	template<typename T>
	mat_t<4, 4, T> orthoZO(T left, T right, T bottom, T top, T zNear, T zFar)
	{
#		if SIMD_COORDINATE_SYSTEM == SIMD_LEFT_HANDED
		return orthoLH_ZO(left, right, bottom, top, zNear, zFar);
#		else
		return orthoRH_ZO(left, right, bottom, top, zNear, zFar);
#		endif
	}
	
	template<typename T>
	mat_t<4, 4, T> orthoNO(T left, T right, T bottom, T top, T zNear, T zFar)
	{
#		if SIMD_COORDINATE_SYSTEM == SIMD_LEFT_HANDED
		return orthoLH_NO(left, right, bottom, top, zNear, zFar);
#		else
		return orthoRH_NO(left, right, bottom, top, zNear, zFar);
#		endif
	}
	
	template<typename T>
	mat_t<4, 4, T> orthoLH(T left, T right, T bottom, T top, T zNear, T zFar)
	{
#		if SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_ZERO_TO_ONE
		return orthoLH_ZO(left, right, bottom, top, zNear, zFar);
#		else
		return orthoLH_NO(left, right, bottom, top, zNear, zFar);
#		endif
	}
	
	template<typename T>
	mat_t<4, 4, T> orthoRH(T left, T right, T bottom, T top, T zNear, T zFar)
	{
#		if SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_ZERO_TO_ONE
		return orthoRH_ZO(left, right, bottom, top, zNear, zFar);
#		else
		return orthoRH_NO(left, right, bottom, top, zNear, zFar);
#		endif
	}
	
	template<typename T>
	mat_t<4, 4, T> orthoMatrix(T left, T right, T bottom, T top, T zNear, T zFar)
	{
#		if SIMD_COORDINATE_SYSTEM == SIMD_LEFT_HANDED && SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_ZERO_TO_ONE
		return orthoLH_ZO(left, right, bottom, top, zNear, zFar);
#		elif SIMD_COORDINATE_SYSTEM == SIMD_LEFT_HANDED && SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_NEGATIVE_ONE_TO_ONE
		return orthoLH_NO(left, right, bottom, top, zNear, zFar);
#		elif SIMD_COORDINATE_SYSTEM == SIMD_RIGHT_HANDED && SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_ZERO_TO_ONE
		return orthoRH_ZO(left, right, bottom, top, zNear, zFar);
#		elif SIMD_COORDINATE_SYSTEM == SIMD_RIGHT_HANDED && SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_NEGATIVE_ONE_TO_ONE
		return orthoRH_NO(left, right, bottom, top, zNear, zFar);
#		endif
	}
	
	template<typename T>
	mat_t<4, 4, T> frustumLH_ZO(T left, T right, T bottom, T top, T nearVal, T farVal)
	{
		mat_t<4, 4, T> Result(0);
		Result[0][0] = (static_cast<T>(2) * nearVal) / (right - left);
		Result[1][1] = (static_cast<T>(2) * nearVal) / (top - bottom);
		Result[2][0] = (right + left) / (right - left);
		Result[2][1] = (top + bottom) / (top - bottom);
		Result[2][2] = farVal / (farVal - nearVal);
		Result[2][3] = static_cast<T>(1);
		Result[3][2] = -(farVal * nearVal) / (farVal - nearVal);
		return Result;
	}
	
	template<typename T>
	mat_t<4, 4, T> frustumLH_NO(T left, T right, T bottom, T top, T nearVal, T farVal)
	{
		mat_t<4, 4, T> Result(0);
		Result[0][0] = (static_cast<T>(2) * nearVal) / (right - left);
		Result[1][1] = (static_cast<T>(2) * nearVal) / (top - bottom);
		Result[2][0] = (right + left) / (right - left);
		Result[2][1] = (top + bottom) / (top - bottom);
		Result[2][2] = (farVal + nearVal) / (farVal - nearVal);
		Result[2][3] = static_cast<T>(1);
		Result[3][2] = - (static_cast<T>(2) * farVal * nearVal) / (farVal - nearVal);
		return Result;
	}
	
	template<typename T>
	mat_t<4, 4, T> frustumRH_ZO(T left, T right, T bottom, T top, T nearVal, T farVal)
	{
		mat_t<4, 4, T> Result(0);
		Result[0][0] = (static_cast<T>(2) * nearVal) / (right - left);
		Result[1][1] = (static_cast<T>(2) * nearVal) / (top - bottom);
		Result[2][0] = (right + left) / (right - left);
		Result[2][1] = (top + bottom) / (top - bottom);
		Result[2][2] = farVal / (nearVal - farVal);
		Result[2][3] = static_cast<T>(-1);
		Result[3][2] = -(farVal * nearVal) / (farVal - nearVal);
		return Result;
	}
	
	template<typename T>
	mat_t<4, 4, T> frustumRH_NO(T left, T right, T bottom, T top, T nearVal, T farVal)
	{
		mat_t<4, 4, T> Result(0);
		Result[0][0] = (static_cast<T>(2) * nearVal) / (right - left);
		Result[1][1] = (static_cast<T>(2) * nearVal) / (top - bottom);
		Result[2][0] = (right + left) / (right - left);
		Result[2][1] = (top + bottom) / (top - bottom);
		Result[2][2] = - (farVal + nearVal) / (farVal - nearVal);
		Result[2][3] = static_cast<T>(-1);
		Result[3][2] = - (static_cast<T>(2) * farVal * nearVal) / (farVal - nearVal);
		return Result;
	}
	
	template<typename T>
	mat_t<4, 4, T> frustumZO(T left, T right, T bottom, T top, T nearVal, T farVal)
	{
#		if SIMD_COORDINATE_SYSTEM == SIMD_LEFT_HANDED
		return frustumLH_ZO(left, right, bottom, top, nearVal, farVal);
#		else
		return frustumRH_ZO(left, right, bottom, top, nearVal, farVal);
#		endif
	}
	
	template<typename T>
	mat_t<4, 4, T> frustumNO(T left, T right, T bottom, T top, T nearVal, T farVal)
	{
#		if SIMD_COORDINATE_SYSTEM == SIMD_LEFT_HANDED
		return frustumLH_NO(left, right, bottom, top, nearVal, farVal);
#		else
		return frustumRH_NO(left, right, bottom, top, nearVal, farVal);
#		endif
	}
	
	template<typename T>
	mat_t<4, 4, T> frustumLH(T left, T right, T bottom, T top, T nearVal, T farVal)
	{
#		if SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_ZERO_TO_ONE
		return frustumLH_ZO(left, right, bottom, top, nearVal, farVal);
#		else
		return frustumLH_NO(left, right, bottom, top, nearVal, farVal);
#		endif
	}
	
	template<typename T>
	mat_t<4, 4, T> frustumRH(T left, T right, T bottom, T top, T nearVal, T farVal)
	{
#		if SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_ZERO_TO_ONE
		return frustumRH_ZO(left, right, bottom, top, nearVal, farVal);
#		else
		return frustumRH_NO(left, right, bottom, top, nearVal, farVal);
#		endif
	}
	
	template<typename T>
	mat_t<4, 4, T> frustumMatrix(T left, T right, T bottom, T top, T nearVal, T farVal)
	{
#		if SIMD_COORDINATE_SYSTEM == SIMD_LEFT_HANDED && SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_ZERO_TO_ONE
		return frustumLH_ZO(left, right, bottom, top, nearVal, farVal);
#		elif SIMD_COORDINATE_SYSTEM == SIMD_LEFT_HANDED && SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_NEGATIVE_ONE_TO_ONE
		return frustumLH_NO(left, right, bottom, top, nearVal, farVal);
#		elif SIMD_COORDINATE_SYSTEM == SIMD_RIGHT_HANDED && SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_ZERO_TO_ONE
		return frustumRH_ZO(left, right, bottom, top, nearVal, farVal);
#		elif SIMD_COORDINATE_SYSTEM == SIMD_RIGHT_HANDED && SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_NEGATIVE_ONE_TO_ONE
		return frustumRH_NO(left, right, bottom, top, nearVal, farVal);
#		endif
	}
	
	template<typename T>
	mat_t<4, 4, T> perspectiveRH_ZO(T fovy, T aspect, T zNear, T zFar)
	{
		assert(abs(aspect - std::numeric_limits<T>::epsilon()) > static_cast<T>(0));
		
		T const tanHalfFovy = tan(fovy / static_cast<T>(2));
		
		mat_t<4, 4, T> Result(static_cast<T>(0));
		Result[0][0] = static_cast<T>(1) / (aspect * tanHalfFovy);
		Result[1][1] = static_cast<T>(1) / (tanHalfFovy);
		Result[2][2] = zFar / (zNear - zFar);
		Result[2][3] = - static_cast<T>(1);
		Result[3][2] = -(zFar * zNear) / (zFar - zNear);
		return Result;
	}
	
	template<typename T>
	mat_t<4, 4, T> perspectiveRH_NO(T fovy, T aspect, T zNear, T zFar)
	{
		assert(abs(aspect - std::numeric_limits<T>::epsilon()) > static_cast<T>(0));
		
		T const tanHalfFovy = tan(fovy / static_cast<T>(2));
		
		mat_t<4, 4, T> Result(static_cast<T>(0));
		Result[0][0] = static_cast<T>(1) / (aspect * tanHalfFovy);
		Result[1][1] = static_cast<T>(1) / (tanHalfFovy);
		Result[2][2] = - (zFar + zNear) / (zFar - zNear);
		Result[2][3] = - static_cast<T>(1);
		Result[3][2] = - (static_cast<T>(2) * zFar * zNear) / (zFar - zNear);
		return Result;
	}
	
	template<typename T>
	mat_t<4, 4, T> perspectiveLH_ZO(T fovy, T aspect, T zNear, T zFar)
	{
		assert(abs(aspect - std::numeric_limits<T>::epsilon()) > static_cast<T>(0));
		
		T const tanHalfFovy = tan(fovy / static_cast<T>(2));
		
		mat_t<4, 4, T> Result(static_cast<T>(0));
		Result[0][0] = static_cast<T>(1) / (aspect * tanHalfFovy);
		Result[1][1] = static_cast<T>(1) / (tanHalfFovy);
		Result[2][2] = zFar / (zFar - zNear);
		Result[2][3] = static_cast<T>(1);
		Result[3][2] = -(zFar * zNear) / (zFar - zNear);
		return Result;
	}
	
	template<typename T>
	mat_t<4, 4, T> perspectiveLH_NO(T fovy, T aspect, T zNear, T zFar)
	{
		assert(abs(aspect - std::numeric_limits<T>::epsilon()) > static_cast<T>(0));
		
		T const tanHalfFovy = tan(fovy / static_cast<T>(2));
		
		mat_t<4, 4, T> Result(static_cast<T>(0));
		Result[0][0] = static_cast<T>(1) / (aspect * tanHalfFovy);
		Result[1][1] = static_cast<T>(1) / (tanHalfFovy);
		Result[2][2] = (zFar + zNear) / (zFar - zNear);
		Result[2][3] = static_cast<T>(1);
		Result[3][2] = - (static_cast<T>(2) * zFar * zNear) / (zFar - zNear);
		return Result;
	}
	
	template<typename T>
	mat_t<4, 4, T> perspectiveZO(T fovy, T aspect, T zNear, T zFar)
	{
#		if SIMD_COORDINATE_SYSTEM == SIMD_LEFT_HANDED
		return perspectiveLH_ZO(fovy, aspect, zNear, zFar);
#		else
		return perspectiveRH_ZO(fovy, aspect, zNear, zFar);
#		endif
	}
	
	template<typename T>
	mat_t<4, 4, T> perspectiveNO(T fovy, T aspect, T zNear, T zFar)
	{
#		if SIMD_COORDINATE_SYSTEM == SIMD_LEFT_HANDED
		return perspectiveLH_NO(fovy, aspect, zNear, zFar);
#		else
		return perspectiveRH_NO(fovy, aspect, zNear, zFar);
#		endif
	}
	
	template<typename T>
	mat_t<4, 4, T> perspectiveLH(T fovy, T aspect, T zNear, T zFar)
	{
#		if SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_ZERO_TO_ONE
		return perspectiveLH_ZO(fovy, aspect, zNear, zFar);
#		else
		return perspectiveLH_NO(fovy, aspect, zNear, zFar);
#		endif
	}
	
	template<typename T>
	mat_t<4, 4, T> perspectiveRH(T fovy, T aspect, T zNear, T zFar)
	{
#		if SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_ZERO_TO_ONE
		return perspectiveRH_ZO(fovy, aspect, zNear, zFar);
#		else
		return perspectiveRH_NO(fovy, aspect, zNear, zFar);
#		endif
	}
	
	template<typename T>
	mat_t<4, 4, T> perspectiveMatrix(T fovy, T aspect, T zNear, T zFar)
	{
#		if SIMD_COORDINATE_SYSTEM == SIMD_LEFT_HANDED && SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_ZERO_TO_ONE
		return perspectiveLH_ZO<T>(fovy, aspect, zNear, zFar);
#		elif SIMD_COORDINATE_SYSTEM == SIMD_LEFT_HANDED && SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_NEGATIVE_ONE_TO_ONE
		return perspectiveLH_NO<T>(fovy, aspect, zNear, zFar);
#		elif SIMD_COORDINATE_SYSTEM == SIMD_RIGHT_HANDED && SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_ZERO_TO_ONE
		return perspectiveRH_ZO<T>(fovy, aspect, zNear, zFar);
#		elif SIMD_COORDINATE_SYSTEM == SIMD_RIGHT_HANDED && SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_NEGATIVE_ONE_TO_ONE
		return perspectiveRH_NO<T>(fovy, aspect, zNear, zFar);
#		endif
	}
	
	template<typename T>
	mat_t<4, 4, T> perspectiveFovRH_ZO(T fov, T width, T height, T zNear, T zFar)
	{
		assert(width > static_cast<T>(0));
		assert(height > static_cast<T>(0));
		assert(fov > static_cast<T>(0));
		
		T const rad = fov;
		T const h = simd::cos(static_cast<T>(0.5) * rad) / simd::sin(static_cast<T>(0.5) * rad);
		T const w = h * height / width; ///todo max(width , Height) / min(width , Height)?
		
		mat_t<4, 4, T> Result(static_cast<T>(0));
		Result[0][0] = w;
		Result[1][1] = h;
		Result[2][2] = zFar / (zNear - zFar);
		Result[2][3] = - static_cast<T>(1);
		Result[3][2] = -(zFar * zNear) / (zFar - zNear);
		return Result;
	}
	
	template<typename T>
	mat_t<4, 4, T> perspectiveFovRH_NO(T fov, T width, T height, T zNear, T zFar)
	{
		assert(width > static_cast<T>(0));
		assert(height > static_cast<T>(0));
		assert(fov > static_cast<T>(0));
		
		T const rad = fov;
		T const h = simd::cos(static_cast<T>(0.5) * rad) / simd::sin(static_cast<T>(0.5) * rad);
		T const w = h * height / width; ///todo max(width , Height) / min(width , Height)?
		
		mat_t<4, 4, T> Result(static_cast<T>(0));
		Result[0][0] = w;
		Result[1][1] = h;
		Result[2][2] = - (zFar + zNear) / (zFar - zNear);
		Result[2][3] = - static_cast<T>(1);
		Result[3][2] = - (static_cast<T>(2) * zFar * zNear) / (zFar - zNear);
		return Result;
	}
	
	template<typename T>
	mat_t<4, 4, T> perspectiveFovLH_ZO(T fov, T width, T height, T zNear, T zFar)
	{
		assert(width > static_cast<T>(0));
		assert(height > static_cast<T>(0));
		assert(fov > static_cast<T>(0));
		
		T const rad = fov;
		T const h = simd::cos(static_cast<T>(0.5) * rad) / simd::sin(static_cast<T>(0.5) * rad);
		T const w = h * height / width; ///todo max(width , Height) / min(width , Height)?
		
		mat_t<4, 4, T> Result(static_cast<T>(0));
		Result[0][0] = w;
		Result[1][1] = h;
		Result[2][2] = zFar / (zFar - zNear);
		Result[2][3] = static_cast<T>(1);
		Result[3][2] = -(zFar * zNear) / (zFar - zNear);
		return Result;
	}
	
	template<typename T>
	mat_t<4, 4, T> perspectiveFovLH_NO(T fov, T width, T height, T zNear, T zFar)
	{
		assert(width > static_cast<T>(0));
		assert(height > static_cast<T>(0));
		assert(fov > static_cast<T>(0));
		
		T const rad = fov;
		T const h = simd::cos(static_cast<T>(0.5) * rad) / simd::sin(static_cast<T>(0.5) * rad);
		T const w = h * height / width; ///todo max(width , Height) / min(width , Height)?
		
		mat_t<4, 4, T> Result(static_cast<T>(0));
		Result[0][0] = w;
		Result[1][1] = h;
		Result[2][2] = (zFar + zNear) / (zFar - zNear);
		Result[2][3] = static_cast<T>(1);
		Result[3][2] = - (static_cast<T>(2) * zFar * zNear) / (zFar - zNear);
		return Result;
	}
	
	template<typename T>
	mat_t<4, 4, T> perspectiveFovZO(T fov, T width, T height, T zNear, T zFar)
	{
#		if SIMD_COORDINATE_SYSTEM == SIMD_LEFT_HANDED
		return perspectiveFovLH_ZO(fov, width, height, zNear, zFar);
#		else
		return perspectiveFovRH_ZO(fov, width, height, zNear, zFar);
#		endif
	}
	
	template<typename T>
	mat_t<4, 4, T> perspectiveFovNO(T fov, T width, T height, T zNear, T zFar)
	{
#		if SIMD_COORDINATE_SYSTEM == SIMD_LEFT_HANDED
		return perspectiveFovLH_NO(fov, width, height, zNear, zFar);
#		else
		return perspectiveFovRH_NO(fov, width, height, zNear, zFar);
#		endif
	}
	
	template<typename T>
	mat_t<4, 4, T> perspectiveFovLH(T fov, T width, T height, T zNear, T zFar)
	{
#		if SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_ZERO_TO_ONE
		return perspectiveFovLH_ZO(fov, width, height, zNear, zFar);
#		else
		return perspectiveFovLH_NO(fov, width, height, zNear, zFar);
#		endif
	}
	
	template<typename T>
	mat_t<4, 4, T> perspectiveFovRH(T fov, T width, T height, T zNear, T zFar)
	{
#		if SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_ZERO_TO_ONE
		return perspectiveFovRH_ZO(fov, width, height, zNear, zFar);
#		else
		return perspectiveFovRH_NO(fov, width, height, zNear, zFar);
#		endif
	}
	
	template<typename T>
	mat_t<4, 4, T> perspectiveFovMatrix(T fov, T width, T height, T zNear, T zFar)
	{
#		if SIMD_COORDINATE_SYSTEM == SIMD_LEFT_HANDED && SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_ZERO_TO_ONE
		return perspectiveFovLH_ZO(fov, width, height, zNear, zFar);
#		elif SIMD_COORDINATE_SYSTEM == SIMD_LEFT_HANDED && SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_NEGATIVE_ONE_TO_ONE
		return perspectiveFovLH_NO(fov, width, height, zNear, zFar);
#		elif SIMD_COORDINATE_SYSTEM == SIMD_RIGHT_HANDED && SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_ZERO_TO_ONE
		return perspectiveFovRH_ZO(fov, width, height, zNear, zFar);
#		elif SIMD_COORDINATE_SYSTEM == SIMD_RIGHT_HANDED && SIMD_DEPTH_CLIP_SPACE == SIMD_DEPTH_NEGATIVE_ONE_TO_ONE
		return perspectiveFovRH_NO(fov, width, height, zNear, zFar);
#		endif
	}
	
	template<typename T>
	mat_t<4, 4, T> infinitePerspectiveRH(T fovy, T aspect, T zNear)
	{
		T const range = tan(fovy / static_cast<T>(2)) * zNear;
		T const left = -range * aspect;
		T const right = range * aspect;
		T const bottom = -range;
		T const top = range;
		
		mat_t<4, 4, T> Result(static_cast<T>(0));
		Result[0][0] = (static_cast<T>(2) * zNear) / (right - left);
		Result[1][1] = (static_cast<T>(2) * zNear) / (top - bottom);
		Result[2][2] = - static_cast<T>(1);
		Result[2][3] = - static_cast<T>(1);
		Result[3][2] = - static_cast<T>(2) * zNear;
		return Result;
	}
	
	template<typename T>
	mat_t<4, 4, T> infinitePerspectiveLH(T fovy, T aspect, T zNear)
	{
		T const range = tan(fovy / static_cast<T>(2)) * zNear;
		T const left = -range * aspect;
		T const right = range * aspect;
		T const bottom = -range;
		T const top = range;
		
		mat_t<4, 4, T> Result(T(0));
		Result[0][0] = (static_cast<T>(2) * zNear) / (right - left);
		Result[1][1] = (static_cast<T>(2) * zNear) / (top - bottom);
		Result[2][2] = static_cast<T>(1);
		Result[2][3] = static_cast<T>(1);
		Result[3][2] = - static_cast<T>(2) * zNear;
		return Result;
	}
	
	template<typename T>
	mat_t<4, 4, T> infinitePerspective(T fovy, T aspect, T zNear)
	{
#		if SIMD_COORDINATE_SYSTEM == SIMD_LEFT_HANDED
		return infinitePerspectiveLH(fovy, aspect, zNear);
#		else
		return infinitePerspectiveRH(fovy, aspect, zNear);
#		endif
	}
	
	template<typename T>
	mat_t<4, 4, T> tweakedInfinitePerspectiveMatrix(T fovy, T aspect, T zNear, T ep)
	{
		T const range = tan(fovy / static_cast<T>(2)) * zNear;
		T const left = -range * aspect;
		T const right = range * aspect;
		T const bottom = -range;
		T const top = range;
		
		mat_t<4, 4, T> Result(static_cast<T>(0));
		Result[0][0] = (static_cast<T>(2) * zNear) / (right - left);
		Result[1][1] = (static_cast<T>(2) * zNear) / (top - bottom);
		Result[2][2] = ep - static_cast<T>(1);
		Result[2][3] = static_cast<T>(-1);
		Result[3][2] = (ep - static_cast<T>(2)) * zNear;
		return Result;
	}
	
	template<typename T>
	mat_t<4, 4, T> tweakedInfinitePerspectiveMatrix(T fovy, T aspect, T zNear)
	{
		return tweakedInfinitePerspective(fovy, aspect, zNear, std::numeric_limits<T>::epsilon());
	}
	
	
}


namespace simd
{
	template<typename T>
	frustum_t<T>	frustum_t<T>::perspective(T fovy, T aspect, T zNear, T zFar)
	{
		frustum_t<T>	ret;
		ret.m_mode = PERSPECTIVE_FOVY_ASPECT;
		ret.m_fovy = fovy;
		ret.m_aspect = aspect;
		ret.m_zNear = zNear;
		ret.m_zFar = zFar;
		return ret;
	}

	template<typename T>
	mat_t<4, 4, T> frustum_t<T>::to_matrix() const
	{
		if (m_mode == PERSPECTIVE_FOVY_ASPECT)
			return perspectiveMatrix<T>(m_fovy, m_aspect, m_zNear, m_zFar);
		return mat_t<4, 4, T>();
	}
	
}





namespace simd
{
#define PROJECTION_DEFINITIONS_T(_T) \
template	mat_t<4, 4, _T> perspectiveMatrix(_T fovy, _T aspect, _T zNear, _T zFar); \
template	frustum_t<_T>	frustum_t<_T>::perspective(_T fovy, _T aspect, _T zNear, _T zFar); \
template	mat_t<4, 4, _T>	to_matrix(const frustum_t<_T>& AFrustum); \


	PROJECTION_DEFINITIONS_T(float)
	PROJECTION_DEFINITIONS_T(double)

}


