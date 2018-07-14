// This software is MIT licensed (see LICENSE)

#pragma once

#include <simd/core.h>

namespace simd
{
	
	// Denavit-Hartenberg
	template<typename T>
	struct denavit_hartenberg_t
	{
		T	o, z, d, a;
	};
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Euler
	// Brief:
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	template<typename T, char A1, char A2, char A3>
	struct euler_t
	{
		T	angle[3];
		FORCE_INLINE	euler_t()							{}
		FORCE_INLINE	euler_t(T i)						{ angle[0] = angle[1] = angle[2] = i; }
		FORCE_INLINE	euler_t(T a0,T a1,T a2)  			{ angle[0] = a0; angle[1] = a1; angle[2] = a2; }
		FORCE_INLINE	auto	operator[](int i) -> T&		{return angle[i];}
	};
	
	template<typename T>
	class axis_t
	{
		pack_t<4, T>	mData;
	public:
		FORCE_INLINE	axis_t()										{}
		FORCE_INLINE	auto	axis() -> pack_t<3, T>& 			{return *(pack_t<3, T>*)&mData;}
		FORCE_INLINE	auto	axis() const -> const pack_t<3, T>&	{return *(const pack_t<3, T>*)&mData;}
		FORCE_INLINE	auto	angle() -> T&						{return mData.w;}
		FORCE_INLINE	auto	angle() const -> const T&			{return mData.w;}
	};
	
	template<typename T, char _Axis>
	struct aligned_axis_t
	{
		T		mAngle;
		
		FORCE_INLINE			aligned_axis_t()									{}
		FORCE_INLINE			aligned_axis_t(const T& iAngle):mAngle(iAngle)		{assert(_Axis=='x' ||_Axis=='y' ||_Axis=='z');}
		static	char		axis_name()											{return _Axis;}
		DLL_FUNCTION(auto)		axis() const -> pack_t<3, T>;
		FORCE_INLINE	auto	angle() -> T&										{return mAngle;}
		FORCE_INLINE	auto	angle() const -> const T&							{return mAngle;}
		
	};
	
	template<typename T>
	struct quaternion_t
	{
		typedef quaternion_t<T> type;
		typedef T value_type;

		union
		{
			struct
			{
				T x, y, z, w;
			};
			typename reg<T, 4>::type_t	m;
		};
		
		/// Return the count of components of a quaternion
		FORCE_INLINE static constexpr int component_count(){return 4;}
		FORCE_INLINE T & operator[](int i);
		FORCE_INLINE T const& operator[](int i) const;
		// -- Implicit basic constructors --
		FORCE_INLINE constexpr quaternion_t();
		FORCE_INLINE constexpr quaternion_t(quaternion_t<T> const& q);
		// -- Explicit basic constructors --
		FORCE_INLINE constexpr quaternion_t(T s, pack_t<3, T> const& v);
		FORCE_INLINE constexpr quaternion_t(T w, T x, T y, T z);
		// -- Conversion constructors --
		template<typename U>
		FORCE_INLINE constexpr explicit quaternion_t(quaternion_t<U> const& q);
		/// Create a quaternion from two normalized axis
		///
		/// @param u A first normalized axis
		/// @param v A second normalized axis
		/// @see gtc_quaternion
		/// @see http://lolengine.net/blog/2013/09/18/beautiful-maths-quaternion-from-vectors
		FORCE_INLINE quaternion_t(pack_t<3, T> const& u, pack_t<3, T> const& v);
		/// Build a quaternion from euler angles (pitch, yaw, roll), in radians.
		FORCE_INLINE explicit quaternion_t(pack_t<3, T> const& eulerAngles);
		FORCE_INLINE explicit quaternion_t(mat_t<3, 3, T> const& q);
		FORCE_INLINE explicit quaternion_t(mat_t<4, 4, T> const& q);
		// -- Unary arithmetic operators --
		FORCE_INLINE quaternion_t<T> & operator=(quaternion_t<T> const& q);
		template<typename U>
		FORCE_INLINE quaternion_t<T> & operator=(quaternion_t<U> const& q);
		template<typename U>
		FORCE_INLINE quaternion_t<T> & operator+=(quaternion_t<U> const& q);
		template<typename U>
		FORCE_INLINE quaternion_t<T> & operator-=(quaternion_t<U> const& q);
		template<typename U>
		FORCE_INLINE quaternion_t<T> & operator*=(quaternion_t<U> const& q);
		template<typename U>
		FORCE_INLINE quaternion_t<T> & operator*=(U s);
		template<typename U>
		FORCE_INLINE quaternion_t<T> & operator/=(U s);
	};
	
	template<typename TP, typename TV>
	struct lookat_t
	{
		pack_t<3, TP>			from, to;
		pack_t<3, TV>			up;

		FORCE_INLINE	lookat_t()																										{}
		lookat_t(const pack_t<3, TP>& FromPosition, const pack_t<3, TP>& ToTarget, const pack_t<3, TV>& UpVector, bool iCorrectUp);
		FORCE_INLINE	auto	vector() const -> pack_t<3, TP> 																		{return to - from;}
		FORCE_INLINE	auto	right() const -> pack_t<3, TV>																			{return cross(vector(), up);}
	};
	
	template<typename _T_Position,typename _T_Spin,typename _T_Scale>
	struct position_t
	{
		_T_Position		translation;
		_T_Spin			spin;
		_T_Scale		scaling;
		FORCE_INLINE	position_t()																					{}
		FORCE_INLINE	position_t(const _T_Position& iPosition, const _T_Spin& iSpin, const _T_Scale& iScaling)		{translation=iPosition;spin=iSpin;scaling=iScaling;}
	};
	

	template<typename T1, typename T2>	mat_t<4, 4, T1> lookAt(pack_t<3, T1> const& eye, pack_t<3, T1> const& center, pack_t<3, T2> const& up);
	template<typename T1, typename T2>	mat_t<4, 4, T1> lookAtDirect(pack_t<3, T1> const& eye, pack_t<3, T1> const& center, pack_t<3, T2> const& up);
}

namespace simd
{
#define DEFINE_EXTERN_EULER_T(_T) \
	extern template struct euler_t<_T, 'x', 'y', 'z'>; \

	DEFINE_EXTERN_EULER_T(float)
	DEFINE_EXTERN_EULER_T(double)

	typedef	euler_t<float, 'x', 'y', 'z'>	euler;
	typedef	euler_t<float, 'x', 'y', 'z'>	f32euler;
	typedef	euler_t<double, 'x', 'y', 'z'>	f64euler;

	
	
//	typedef aligned_axis_t<float,'x'>							axis_x_f32;
//	typedef aligned_axis_t<float,'y'>							axis_y_f32;
//	typedef aligned_axis_t<float,'z'>							axis_z_f32;
//	typedef aligned_axis_t<double,'x'>							axis_x_f64;
//	typedef aligned_axis_t<double,'y'>							axis_y_f64;
//	typedef aligned_axis_t<double,'z'>							axis_z_f64;
//	typedef euler_t<float, 'x', 'y', 'z'>						euler_f32;
//	typedef euler_t<double, 'x', 'y', 'z'>					euler_f64;
//	typedef	lookat_t<float>									lookat_f32;
//	typedef	lookat_t<double>									lookat_f64;
//	typedef	position_t<simd::TPack<float, 3, false>,euler_f32,simd::TPack<float, 3, false>>	position3d_f32;
//	typedef	position_t<simd::TPack<double, 3, false>,euler_f32,simd::TPack<double, 3, false>>	position3d_f64;
//
//	typedef	denavit_hartenberg_t<float>						denavit_hartenberg_f32;
//	typedef	denavit_hartenberg_t<double>						denavit_hartenberg_f64;
//	typedef	denavit_hartenberg_t<double>						denavit_hartenberg_t;
}

#include "inline/pose.h"
