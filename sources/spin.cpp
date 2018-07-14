// This software is MIT licensed (see LICENSE)

#define BUILD_DLL

#include <simd/pose.h>

//namespace glm
//{
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> quat_identity()
//	{
//		return tquat<T, Q>(static_cast<T>(1), static_cast<T>(0), static_cast<T>(0), static_cast<T>(0));
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER vec<3, T, Q> cross(vec<3, T, Q> const& v, tquat<T, Q> const& q)
//	{
//		return inverse(q) * v;
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER vec<3, T, Q> cross(tquat<T, Q> const& q, vec<3, T, Q> const& v)
//	{
//		return q * v;
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> squad
//	(
//	 tquat<T, Q> const& q1,
//	 tquat<T, Q> const& q2,
//	 tquat<T, Q> const& s1,
//	 tquat<T, Q> const& s2,
//	 T const& h)
//	{
//		return mix(mix(q1, q2, h), mix(s1, s2, h), static_cast<T>(2) * (static_cast<T>(1) - h) * h);
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> intermediate
//	(
//	 tquat<T, Q> const& prev,
//	 tquat<T, Q> const& curr,
//	 tquat<T, Q> const& next
//	 )
//	{
//		tquat<T, Q> invQuat = inverse(curr);
//		return exp((log(next * invQuat) + log(prev * invQuat)) / static_cast<T>(-4)) * curr;
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> exp(tquat<T, Q> const& q)
//	{
//		vec<3, T, Q> u(q.x, q.y, q.z);
//		T const Angle = glm::length(u);
//		if (Angle < epsilon<T>())
//			return tquat<T, Q>();
//		
//		vec<3, T, Q> const v(u / Angle);
//		return tquat<T, Q>(cos(Angle), sin(Angle) * v);
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> log(tquat<T, Q> const& q)
//	{
//		vec<3, T, Q> u(q.x, q.y, q.z);
//		T Vec3Len = length(u);
//		
//		if (Vec3Len < epsilon<T>())
//		{
//			if(q.w > static_cast<T>(0))
//				return tquat<T, Q>(log(q.w), static_cast<T>(0), static_cast<T>(0), static_cast<T>(0));
//			else if(q.w < static_cast<T>(0))
//				return tquat<T, Q>(log(-q.w), pi<T>(), static_cast<T>(0), static_cast<T>(0));
//			else
//				return tquat<T, Q>(std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity());
//		}
//		else
//		{
//			T t = atan(Vec3Len, T(q.w)) / Vec3Len;
//			T QuatLen2 = Vec3Len * Vec3Len + q.w * q.w;
//			return tquat<T, Q>(static_cast<T>(0.5) * log(QuatLen2), t * q.x, t * q.y, t * q.z);
//		}
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> pow(tquat<T, Q> const& x, T const& y)
//	{
//		//Raising to the power of 0 should yield 1
//		//Needed to prevent a division by 0 error later on
//		if(y > -epsilon<T>() && y < epsilon<T>())
//			return tquat<T, Q>(1,0,0,0);
//		
//		//To deal with non-unit quaternions
//		T magnitude = sqrt(x.x * x.x + x.y * x.y + x.z * x.z + x.w *x.w);
//		
//		//Equivalent to raising a real number to a power
//		//Needed to prevent a division by 0 error later on
//		if(abs(x.w / magnitude) > static_cast<T>(1) - epsilon<T>() && abs(x.w / magnitude) < static_cast<T>(1) + epsilon<T>())
//			return tquat<T, Q>(pow(x.w, y),0,0,0);
//		
//		T Angle = acos(x.w / magnitude);
//		T NewAngle = Angle * y;
//		T Div = sin(NewAngle) / sin(Angle);
//		T Mag = pow(magnitude, y - static_cast<T>(1));
//		
//		return tquat<T, Q>(cos(NewAngle) * magnitude * Mag, x.x * Div * Mag, x.y * Div * Mag, x.z * Div * Mag);
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER vec<3, T, Q> rotate(tquat<T, Q> const& q, vec<3, T, Q> const& v)
//	{
//		return q * v;
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER vec<4, T, Q> rotate(tquat<T, Q> const& q, vec<4, T, Q> const& v)
//	{
//		return q * v;
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER T extractRealComponent(tquat<T, Q> const& q)
//	{
//		T w = static_cast<T>(1) - q.x * q.x - q.y * q.y - q.z * q.z;
//		if(w < T(0))
//			return T(0);
//		else
//			return -sqrt(w);
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER T length2(tquat<T, Q> const& q)
//	{
//		return q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w;
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> shortMix(tquat<T, Q> const& x, tquat<T, Q> const& y, T const& a)
//	{
//		if(a <= static_cast<T>(0)) return x;
//		if(a >= static_cast<T>(1)) return y;
//		
//		T fCos = dot(x, y);
//		tquat<T, Q> y2(y); //BUG!!! tquat<T> y2;
//		if(fCos < static_cast<T>(0))
//		{
//			y2 = -y;
//			fCos = -fCos;
//		}
//		
//		//if(fCos > 1.0f) // problem
//		T k0, k1;
//		if(fCos > (static_cast<T>(1) - epsilon<T>()))
//		{
//			k0 = static_cast<T>(1) - a;
//			k1 = static_cast<T>(0) + a; //BUG!!! 1.0f + a;
//		}
//		else
//		{
//			T fSin = sqrt(T(1) - fCos * fCos);
//			T fAngle = atan(fSin, fCos);
//			T fOneOverSin = static_cast<T>(1) / fSin;
//			k0 = sin((static_cast<T>(1) - a) * fAngle) * fOneOverSin;
//			k1 = sin((static_cast<T>(0) + a) * fAngle) * fOneOverSin;
//		}
//		
//		return tquat<T, Q>(
//						   k0 * x.w + k1 * y2.w,
//						   k0 * x.x + k1 * y2.x,
//						   k0 * x.y + k1 * y2.y,
//						   k0 * x.z + k1 * y2.z);
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> fastMix(tquat<T, Q> const& x, tquat<T, Q> const& y, T const& a)
//	{
//		return glm::normalize(x * (static_cast<T>(1) - a) + (y * a));
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> rotation(vec<3, T, Q> const& orig, vec<3, T, Q> const& dest)
//	{
//		T cosTheta = dot(orig, dest);
//		vec<3, T, Q> rotationAxis;
//		
//		if(cosTheta >= static_cast<T>(1) - epsilon<T>()) {
//			// orig and dest point in the same direction
//			return quat_identity<T,Q>();
//		}
//		
//		if(cosTheta < static_cast<T>(-1) + epsilon<T>())
//		{
//			// special case when vectors in opposite directions :
//			// there is no "ideal" rotation axis
//			// So guess one; any will do as long as it's perpendicular to start
//			// This implementation favors a rotation around the Up axis (Y),
//			// since it's often what you want to do.
//			rotationAxis = cross(vec<3, T, Q>(0, 0, 1), orig);
//			if(length2(rotationAxis) < epsilon<T>()) // bad luck, they were parallel, try again!
//				rotationAxis = cross(vec<3, T, Q>(1, 0, 0), orig);
//			
//			rotationAxis = normalize(rotationAxis);
//			return angleAxis(pi<T>(), rotationAxis);
//		}
//		
//		// Implementation from Stan Melax's Game Programming Gems 1 article
//		rotationAxis = cross(orig, dest);
//		
//		T s = sqrt((T(1) + cosTheta) * static_cast<T>(2));
//		T invs = static_cast<T>(1) / s;
//		
//		return tquat<T, Q>(
//						   s * static_cast<T>(0.5f),
//						   rotationAxis.x * invs,
//						   rotationAxis.y * invs,
//						   rotationAxis.z * invs);
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> quatLookAt(vec<3, T, Q> const& direction, vec<3, T, Q> const& up)
//	{
//#		if GLM_COORDINATE_SYSTEM == GLM_LEFT_HANDED
//		return quatLookAtLH(direction, up);
//#		else
//		return quatLookAtRH(direction, up);
//# 		endif
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> quatLookAtRH(vec<3, T, Q> const& direction, vec<3, T, Q> const& up)
//	{
//		mat<3, 3, T, Q> Result;
//		
//		Result[2] = -direction;
//		Result[0] = normalize(cross(up, Result[2]));
//		Result[1] = cross(Result[2], Result[0]);
//		
//		return quat_cast(Result);
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> quatLookAtLH(vec<3, T, Q> const& direction, vec<3, T, Q> const& up)
//	{
//		mat<3, 3, T, Q> Result;
//		
//		Result[2] = direction;
//		Result[0] = normalize(cross(up, Result[2]));
//		Result[1] = cross(Result[2], Result[0]);
//		
//		return quat_cast(Result);
//	}
//	
//}//namespace glm
//
//
//
//
//
//
//namespace glm{
//	namespace detail
//	{
//		template<typename T, qualifier Q, bool Aligned>
//		struct compute_dot<tquat<T, Q>, T, Aligned>
//		{
//			static GLM_FUNC_QUALIFIER T call(tquat<T, Q> const& a, tquat<T, Q> const& b)
//			{
//				vec<4, T, Q> tmp(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w);
//				return (tmp.x + tmp.y) + (tmp.z + tmp.w);
//			}
//		};
//		
//		template<typename T, qualifier Q, bool Aligned>
//		struct compute_quat_add
//		{
//			static tquat<T, Q> call(tquat<T, Q> const& q, tquat<T, Q> const& p)
//			{
//				return tquat<T, Q>(q.w + p.w, q.x + p.x, q.y + p.y, q.z + p.z);
//			}
//		};
//		
//		template<typename T, qualifier Q, bool Aligned>
//		struct compute_quat_sub
//		{
//			static tquat<T, Q> call(tquat<T, Q> const& q, tquat<T, Q> const& p)
//			{
//				return tquat<T, Q>(q.w - p.w, q.x - p.x, q.y - p.y, q.z - p.z);
//			}
//		};
//		
//		template<typename T, qualifier Q, bool Aligned>
//		struct compute_quat_mul_scalar
//		{
//			static tquat<T, Q> call(tquat<T, Q> const& q, T s)
//			{
//				return tquat<T, Q>(q.w * s, q.x * s, q.y * s, q.z * s);
//			}
//		};
//		
//		template<typename T, qualifier Q, bool Aligned>
//		struct compute_quat_div_scalar
//		{
//			static tquat<T, Q> call(tquat<T, Q> const& q, T s)
//			{
//				return tquat<T, Q>(q.w / s, q.x / s, q.y / s, q.z / s);
//			}
//		};
//		
//		template<typename T, qualifier Q, bool Aligned>
//		struct compute_quat_mul_vec4
//		{
//			static vec<4, T, Q> call(tquat<T, Q> const& q, vec<4, T, Q> const& v)
//			{
//				return vec<4, T, Q>(q * vec<3, T, Q>(v), v.w);
//			}
//		};
//	}//namespace detail
//	
//	// -- Component accesses --
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER T & tquat<T, Q>::operator[](typename tquat<T, Q>::length_type i)
//	{
//		assert(i >= 0 && i < this->length());
//		return (&x)[i];
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER T const& tquat<T, Q>::operator[](typename tquat<T, Q>::length_type i) const
//	{
//		assert(i >= 0 && i < this->length());
//		return (&x)[i];
//	}
//	
//	// -- Implicit basic constructors --
//	
//#	if !GLM_HAS_DEFAULTED_FUNCTIONS || defined(GLM_FORCE_CTOR_INIT)
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER GLM_CONSTEXPR tquat<T, Q>::tquat()
//#			ifdef GLM_FORCE_CTOR_INIT
//	: x(0), y(0), z(0), w(1)
//#			endif
//	{}
//#	endif
//	
//#	if !GLM_HAS_DEFAULTED_FUNCTIONS
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER GLM_CONSTEXPR tquat<T, Q>::tquat(tquat<T, Q> const& q)
//	: x(q.x), y(q.y), z(q.z), w(q.w)
//	{}
//#	endif//!GLM_HAS_DEFAULTED_FUNCTIONS
//	
//	template<typename T, qualifier Q>
//	template<qualifier P>
//	GLM_FUNC_QUALIFIER GLM_CONSTEXPR tquat<T, Q>::tquat(tquat<T, P> const& q)
//	: x(q.x), y(q.y), z(q.z), w(q.w)
//	{}
//	
//	// -- Explicit basic constructors --
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER GLM_CONSTEXPR tquat<T, Q>::tquat(T s, vec<3, T, Q> const& v)
//	: x(v.x), y(v.y), z(v.z), w(s)
//	{}
//	
//	template <typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER GLM_CONSTEXPR tquat<T, Q>::tquat(T _w, T _x, T _y, T _z)
//	: x(_x), y(_y), z(_z), w(_w)
//	{}
//	
//	// -- Conversion constructors --
//	
//	template<typename T, qualifier Q>
//	template<typename U, qualifier P>
//	GLM_FUNC_QUALIFIER GLM_CONSTEXPR tquat<T, Q>::tquat(tquat<U, P> const& q)
//	: x(static_cast<T>(q.x))
//	, y(static_cast<T>(q.y))
//	, z(static_cast<T>(q.z))
//	, w(static_cast<T>(q.w))
//	{}
//	
//	//template<typename valType>
//	//GLM_FUNC_QUALIFIER tquat<valType>::tquat
//	//(
//	//	valType const& pitch,
//	//	valType const& yaw,
//	//	valType const& roll
//	//)
//	//{
//	//	vec<3, valType> eulerAngle(pitch * valType(0.5), yaw * valType(0.5), roll * valType(0.5));
//	//	vec<3, valType> c = glm::cos(eulerAngle * valType(0.5));
//	//	vec<3, valType> s = glm::sin(eulerAngle * valType(0.5));
//	//
//	//	this->w = c.x * c.y * c.z + s.x * s.y * s.z;
//	//	this->x = s.x * c.y * c.z - c.x * s.y * s.z;
//	//	this->y = c.x * s.y * c.z + s.x * c.y * s.z;
//	//	this->z = c.x * c.y * s.z - s.x * s.y * c.z;
//	//}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q>::tquat(vec<3, T, Q> const& u, vec<3, T, Q> const& v)
//	{
//		T norm_u_norm_v = sqrt(dot(u, u) * dot(v, v));
//		T real_part = norm_u_norm_v + dot(u, v);
//		vec<3, T, Q> t;
//		
//		if(real_part < static_cast<T>(1.e-6f) * norm_u_norm_v)
//		{
//			// If u and v are exactly opposite, rotate 180 degrees
//			// around an arbitrary orthogonal axis. Axis normalisation
//			// can happen later, when we normalise the quaternion.
//			real_part = static_cast<T>(0);
//			t = abs(u.x) > abs(u.z) ? vec<3, T, Q>(-u.y, u.x, static_cast<T>(0)) : vec<3, T, Q>(static_cast<T>(0), -u.z, u.y);
//		}
//		else
//		{
//			// Otherwise, build quaternion the standard way.
//			t = cross(u, v);
//		}
//		
//		*this = normalize(tquat<T, Q>(real_part, t.x, t.y, t.z));
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q>::tquat(vec<3, T, Q> const& eulerAngle)
//	{
//		vec<3, T, Q> c = glm::cos(eulerAngle * T(0.5));
//		vec<3, T, Q> s = glm::sin(eulerAngle * T(0.5));
//		
//		this->w = c.x * c.y * c.z + s.x * s.y * s.z;
//		this->x = s.x * c.y * c.z - c.x * s.y * s.z;
//		this->y = c.x * s.y * c.z + s.x * c.y * s.z;
//		this->z = c.x * c.y * s.z - s.x * s.y * c.z;
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q>::tquat(mat<3, 3, T, Q> const& m)
//	{
//		*this = quat_cast(m);
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q>::tquat(mat<4, 4, T, Q> const& m)
//	{
//		*this = quat_cast(m);
//	}
//	
//#	if GLM_HAS_EXPLICIT_CONVERSION_OPERATORS
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q>::operator mat<3, 3, T, Q>()
//	{
//		return mat3_cast(*this);
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q>::operator mat<4, 4, T, Q>()
//	{
//		return mat4_cast(*this);
//	}
//#	endif//GLM_HAS_EXPLICIT_CONVERSION_OPERATORS
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> conjugate(tquat<T, Q> const& q)
//	{
//		return tquat<T, Q>(q.w, -q.x, -q.y, -q.z);
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> inverse(tquat<T, Q> const& q)
//	{
//		return conjugate(q) / dot(q, q);
//	}
//	
//	// -- Unary arithmetic operators --
//	
//#	if !GLM_HAS_DEFAULTED_FUNCTIONS
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> & tquat<T, Q>::operator=(tquat<T, Q> const& q)
//	{
//		this->w = q.w;
//		this->x = q.x;
//		this->y = q.y;
//		this->z = q.z;
//		return *this;
//	}
//#	endif//!GLM_HAS_DEFAULTED_FUNCTIONS
//	
//	template<typename T, qualifier Q>
//	template<typename U>
//	GLM_FUNC_QUALIFIER tquat<T, Q> & tquat<T, Q>::operator=(tquat<U, Q> const& q)
//	{
//		this->w = static_cast<T>(q.w);
//		this->x = static_cast<T>(q.x);
//		this->y = static_cast<T>(q.y);
//		this->z = static_cast<T>(q.z);
//		return *this;
//	}
//	
//	template<typename T, qualifier Q>
//	template<typename U>
//	GLM_FUNC_QUALIFIER tquat<T, Q> & tquat<T, Q>::operator+=(tquat<U, Q> const& q)
//	{
//		return (*this = detail::compute_quat_add<T, Q, detail::is_aligned<Q>::value>::call(*this, tquat<T, Q>(q)));
//	}
//	
//	template<typename T, qualifier Q>
//	template<typename U>
//	GLM_FUNC_QUALIFIER tquat<T, Q> & tquat<T, Q>::operator-=(tquat<U, Q> const& q)
//	{
//		return (*this = detail::compute_quat_sub<T, Q, detail::is_aligned<Q>::value>::call(*this, tquat<T, Q>(q)));
//	}
//	
//	template<typename T, qualifier Q>
//	template<typename U>
//	GLM_FUNC_QUALIFIER tquat<T, Q> & tquat<T, Q>::operator*=(tquat<U, Q> const& r)
//	{
//		tquat<T, Q> const p(*this);
//		tquat<T, Q> const q(r);
//		
//		this->w = p.w * q.w - p.x * q.x - p.y * q.y - p.z * q.z;
//		this->x = p.w * q.x + p.x * q.w + p.y * q.z - p.z * q.y;
//		this->y = p.w * q.y + p.y * q.w + p.z * q.x - p.x * q.z;
//		this->z = p.w * q.z + p.z * q.w + p.x * q.y - p.y * q.x;
//		return *this;
//	}
//	
//	template<typename T, qualifier Q>
//	template<typename U>
//	GLM_FUNC_QUALIFIER tquat<T, Q> & tquat<T, Q>::operator*=(U s)
//	{
//		return (*this = detail::compute_quat_mul_scalar<T, Q, detail::is_aligned<Q>::value>::call(*this, static_cast<U>(s)));
//	}
//	
//	template<typename T, qualifier Q>
//	template<typename U>
//	GLM_FUNC_QUALIFIER tquat<T, Q> & tquat<T, Q>::operator/=(U s)
//	{
//		return (*this = detail::compute_quat_div_scalar<T, Q, detail::is_aligned<Q>::value>::call(*this, static_cast<U>(s)));
//	}
//	
//	// -- Unary bit operators --
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> operator+(tquat<T, Q> const& q)
//	{
//		return q;
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> operator-(tquat<T, Q> const& q)
//	{
//		return tquat<T, Q>(-q.w, -q.x, -q.y, -q.z);
//	}
//	
//	// -- Binary operators --
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> operator+(tquat<T, Q> const& q, tquat<T, Q> const& p)
//	{
//		return tquat<T, Q>(q) += p;
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> operator-(tquat<T, Q> const& q, tquat<T, Q> const& p)
//	{
//		return tquat<T, Q>(q) -= p;
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> operator*(tquat<T, Q> const& q, tquat<T, Q> const& p)
//	{
//		return tquat<T, Q>(q) *= p;
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER vec<3, T, Q> operator*(tquat<T, Q> const& q, vec<3, T, Q> const& v)
//	{
//		vec<3, T, Q> const QuatVector(q.x, q.y, q.z);
//		vec<3, T, Q> const uv(glm::cross(QuatVector, v));
//		vec<3, T, Q> const uuv(glm::cross(QuatVector, uv));
//		
//		return v + ((uv * q.w) + uuv) * static_cast<T>(2);
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER vec<3, T, Q> operator*(vec<3, T, Q> const& v, tquat<T, Q> const& q)
//	{
//		return glm::inverse(q) * v;
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER vec<4, T, Q> operator*(tquat<T, Q> const& q, vec<4, T, Q> const& v)
//	{
//		return detail::compute_quat_mul_vec4<T, Q, detail::is_aligned<Q>::value>::call(q, v);
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER vec<4, T, Q> operator*(vec<4, T, Q> const& v, tquat<T, Q> const& q)
//	{
//		return glm::inverse(q) * v;
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> operator*(tquat<T, Q> const& q, T const& s)
//	{
//		return tquat<T, Q>(
//						   q.w * s, q.x * s, q.y * s, q.z * s);
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> operator*(T const& s, tquat<T, Q> const& q)
//	{
//		return q * s;
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> operator/(tquat<T, Q> const& q, T const& s)
//	{
//		return tquat<T, Q>(
//						   q.w / s, q.x / s, q.y / s, q.z / s);
//	}
//	
//	// -- Boolean operators --
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER bool operator==(tquat<T, Q> const& q1, tquat<T, Q> const& q2)
//	{
//		return all(epsilonEqual(q1, q2, epsilon<T>()));
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER bool operator!=(tquat<T, Q> const& q1, tquat<T, Q> const& q2)
//	{
//		return any(epsilonNotEqual(q1, q2, epsilon<T>()));
//	}
//	
//	// -- Operations --
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER T dot(tquat<T, Q> const& x, tquat<T, Q> const& y)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'dot' accepts only floating-point inputs");
//		return detail::compute_dot<tquat<T, Q>, T, detail::is_aligned<Q>::value>::call(x, y);
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER T length(tquat<T, Q> const& q)
//	{
//		return glm::sqrt(dot(q, q));
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> normalize(tquat<T, Q> const& q)
//	{
//		T len = length(q);
//		if(len <= T(0)) // Problem
//			return tquat<T, Q>(static_cast<T>(1), static_cast<T>(0), static_cast<T>(0), static_cast<T>(0));
//		T oneOverLen = T(1) / len;
//		return tquat<T, Q>(q.w * oneOverLen, q.x * oneOverLen, q.y * oneOverLen, q.z * oneOverLen);
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> cross(tquat<T, Q> const& q1, tquat<T, Q> const& q2)
//	{
//		return tquat<T, Q>(
//						   q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z,
//						   q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y,
//						   q1.w * q2.y + q1.y * q2.w + q1.z * q2.x - q1.x * q2.z,
//						   q1.w * q2.z + q1.z * q2.w + q1.x * q2.y - q1.y * q2.x);
//	}
//	/*
//	 // (x * sin(1 - a) * angle / sin(angle)) + (y * sin(a) * angle / sin(angle))
//	 template<typename T, qualifier Q>
//	 GLM_FUNC_QUALIFIER tquat<T, Q> mix(tquat<T, Q> const& x, tquat<T, Q> const& y, T const& a)
//	 {
//	 if(a <= T(0)) return x;
//	 if(a >= T(1)) return y;
//	 
//	 float fCos = dot(x, y);
//	 tquat<T, Q> y2(y); //BUG!!! tquat<T, Q> y2;
//	 if(fCos < T(0))
//	 {
//	 y2 = -y;
//	 fCos = -fCos;
//	 }
//	 
//	 //if(fCos > 1.0f) // problem
//	 float k0, k1;
//	 if(fCos > T(0.9999))
//	 {
//	 k0 = T(1) - a;
//	 k1 = T(0) + a; //BUG!!! 1.0f + a;
//	 }
//	 else
//	 {
//	 T fSin = sqrt(T(1) - fCos * fCos);
//	 T fAngle = atan(fSin, fCos);
//	 T fOneOverSin = static_cast<T>(1) / fSin;
//	 k0 = sin((T(1) - a) * fAngle) * fOneOverSin;
//	 k1 = sin((T(0) + a) * fAngle) * fOneOverSin;
//	 }
//	 
//	 return tquat<T, Q>(
//	 k0 * x.w + k1 * y2.w,
//	 k0 * x.x + k1 * y2.x,
//	 k0 * x.y + k1 * y2.y,
//	 k0 * x.z + k1 * y2.z);
//	 }
//	 
//	 template<typename T, qualifier Q>
//	 GLM_FUNC_QUALIFIER tquat<T, Q> mix2
//	 (
//	 tquat<T, Q> const& x,
//	 tquat<T, Q> const& y,
//	 T const& a
//	 )
//	 {
//	 bool flip = false;
//	 if(a <= static_cast<T>(0)) return x;
//	 if(a >= static_cast<T>(1)) return y;
//	 
//	 T cos_t = dot(x, y);
//	 if(cos_t < T(0))
//	 {
//	 cos_t = -cos_t;
//	 flip = true;
//	 }
//	 
//	 T alpha(0), beta(0);
//	 
//	 if(T(1) - cos_t < 1e-7)
//	 beta = static_cast<T>(1) - alpha;
//	 else
//	 {
//	 T theta = acos(cos_t);
//	 T sin_t = sin(theta);
//	 beta = sin(theta * (T(1) - alpha)) / sin_t;
//	 alpha = sin(alpha * theta) / sin_t;
//	 }
//	 
//	 if(flip)
//	 alpha = -alpha;
//	 
//	 return normalize(beta * x + alpha * y);
//	 }
//	 */
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> mix(tquat<T, Q> const& x, tquat<T, Q> const& y, T a)
//	{
//		T cosTheta = dot(x, y);
//		
//		// Perform a linear interpolation when cosTheta is close to 1 to avoid side effect of sin(angle) becoming a zero denominator
//		if(cosTheta > T(1) - epsilon<T>())
//		{
//			// Linear interpolation
//			return tquat<T, Q>(
//							   mix(x.w, y.w, a),
//							   mix(x.x, y.x, a),
//							   mix(x.y, y.y, a),
//							   mix(x.z, y.z, a));
//		}
//		else
//		{
//			// Essential Mathematics, page 467
//			T angle = acos(cosTheta);
//			return (sin((T(1) - a) * angle) * x + sin(a * angle) * y) / sin(angle);
//		}
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> lerp(tquat<T, Q> const& x, tquat<T, Q> const& y, T a)
//	{
//		// Lerp is only defined in [0, 1]
//		assert(a >= static_cast<T>(0));
//		assert(a <= static_cast<T>(1));
//		
//		return x * (T(1) - a) + (y * a);
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> slerp(tquat<T, Q> const& x,	tquat<T, Q> const& y, T a)
//	{
//		tquat<T, Q> z = y;
//		
//		T cosTheta = dot(x, y);
//		
//		// If cosTheta < 0, the interpolation will take the long way around the sphere.
//		// To fix this, one quat must be negated.
//		if (cosTheta < T(0))
//		{
//			z        = -y;
//			cosTheta = -cosTheta;
//		}
//		
//		// Perform a linear interpolation when cosTheta is close to 1 to avoid side effect of sin(angle) becoming a zero denominator
//		if(cosTheta > T(1) - epsilon<T>())
//		{
//			// Linear interpolation
//			return tquat<T, Q>(
//							   mix(x.w, z.w, a),
//							   mix(x.x, z.x, a),
//							   mix(x.y, z.y, a),
//							   mix(x.z, z.z, a));
//		}
//		else
//		{
//			// Essential Mathematics, page 467
//			T angle = acos(cosTheta);
//			return (sin((T(1) - a) * angle) * x + sin(a * angle) * z) / sin(angle);
//		}
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> rotate(tquat<T, Q> const& q, T const& angle, vec<3, T, Q> const& v)
//	{
//		vec<3, T, Q> Tmp = v;
//		
//		// Axis of rotation must be normalised
//		T len = glm::length(Tmp);
//		if(abs(len - T(1)) > T(0.001))
//		{
//			T oneOverLen = static_cast<T>(1) / len;
//			Tmp.x *= oneOverLen;
//			Tmp.y *= oneOverLen;
//			Tmp.z *= oneOverLen;
//		}
//		
//		T const AngleRad(angle);
//		T const Sin = sin(AngleRad * T(0.5));
//		
//		return q * tquat<T, Q>(cos(AngleRad * T(0.5)), Tmp.x * Sin, Tmp.y * Sin, Tmp.z * Sin);
//		//return gtc::quaternion::cross(q, tquat<T, Q>(cos(AngleRad * T(0.5)), Tmp.x * fSin, Tmp.y * fSin, Tmp.z * fSin));
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER vec<3, T, Q> eulerAngles(tquat<T, Q> const& x)
//	{
//		return vec<3, T, Q>(pitch(x), yaw(x), roll(x));
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER T roll(tquat<T, Q> const& q)
//	{
//		return static_cast<T>(atan(static_cast<T>(2) * (q.x * q.y + q.w * q.z), q.w * q.w + q.x * q.x - q.y * q.y - q.z * q.z));
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER T pitch(tquat<T, Q> const& q)
//	{
//		//return T(atan(T(2) * (q.y * q.z + q.w * q.x), q.w * q.w - q.x * q.x - q.y * q.y + q.z * q.z));
//		const T y = static_cast<T>(2) * (q.y * q.z + q.w * q.x);
//		const T x = q.w * q.w - q.x * q.x - q.y * q.y + q.z * q.z;
//		
//		if(detail::compute_equal<T>::call(y, static_cast<T>(0)) && detail::compute_equal<T>::call(x, static_cast<T>(0))) //avoid atan2(0,0) - handle singularity - Matiis
//			return static_cast<T>(static_cast<T>(2) * atan(q.x,q.w));
//		
//		return static_cast<T>(atan(y,x));
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER T yaw(tquat<T, Q> const& q)
//	{
//		return asin(clamp(static_cast<T>(-2) * (q.x * q.z - q.w * q.y), static_cast<T>(-1), static_cast<T>(1)));
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER mat<3, 3, T, Q> mat3_cast(tquat<T, Q> const& q)
//	{
//		mat<3, 3, T, Q> Result(T(1));
//		T qxx(q.x * q.x);
//		T qyy(q.y * q.y);
//		T qzz(q.z * q.z);
//		T qxz(q.x * q.z);
//		T qxy(q.x * q.y);
//		T qyz(q.y * q.z);
//		T qwx(q.w * q.x);
//		T qwy(q.w * q.y);
//		T qwz(q.w * q.z);
//		
//		Result[0][0] = T(1) - T(2) * (qyy +  qzz);
//		Result[0][1] = T(2) * (qxy + qwz);
//		Result[0][2] = T(2) * (qxz - qwy);
//		
//		Result[1][0] = T(2) * (qxy - qwz);
//		Result[1][1] = T(1) - T(2) * (qxx +  qzz);
//		Result[1][2] = T(2) * (qyz + qwx);
//		
//		Result[2][0] = T(2) * (qxz + qwy);
//		Result[2][1] = T(2) * (qyz - qwx);
//		Result[2][2] = T(1) - T(2) * (qxx +  qyy);
//		return Result;
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER mat<4, 4, T, Q> mat4_cast(tquat<T, Q> const& q)
//	{
//		return mat<4, 4, T, Q>(mat3_cast(q));
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> quat_cast(mat<3, 3, T, Q> const& m)
//	{
//		T fourXSquaredMinus1 = m[0][0] - m[1][1] - m[2][2];
//		T fourYSquaredMinus1 = m[1][1] - m[0][0] - m[2][2];
//		T fourZSquaredMinus1 = m[2][2] - m[0][0] - m[1][1];
//		T fourWSquaredMinus1 = m[0][0] + m[1][1] + m[2][2];
//		
//		int biggestIndex = 0;
//		T fourBiggestSquaredMinus1 = fourWSquaredMinus1;
//		if(fourXSquaredMinus1 > fourBiggestSquaredMinus1)
//		{
//			fourBiggestSquaredMinus1 = fourXSquaredMinus1;
//			biggestIndex = 1;
//		}
//		if(fourYSquaredMinus1 > fourBiggestSquaredMinus1)
//		{
//			fourBiggestSquaredMinus1 = fourYSquaredMinus1;
//			biggestIndex = 2;
//		}
//		if(fourZSquaredMinus1 > fourBiggestSquaredMinus1)
//		{
//			fourBiggestSquaredMinus1 = fourZSquaredMinus1;
//			biggestIndex = 3;
//		}
//		
//		T biggestVal = sqrt(fourBiggestSquaredMinus1 + static_cast<T>(1)) * static_cast<T>(0.5);
//		T mult = static_cast<T>(0.25) / biggestVal;
//		
//		switch(biggestIndex)
//		{
//			case 0:
//				return tquat<T, Q>(biggestVal, (m[1][2] - m[2][1]) * mult, (m[2][0] - m[0][2]) * mult, (m[0][1] - m[1][0]) * mult);
//			case 1:
//				return tquat<T, Q>((m[1][2] - m[2][1]) * mult, biggestVal, (m[0][1] + m[1][0]) * mult, (m[2][0] + m[0][2]) * mult);
//			case 2:
//				return tquat<T, Q>((m[2][0] - m[0][2]) * mult, (m[0][1] + m[1][0]) * mult, biggestVal, (m[1][2] + m[2][1]) * mult);
//			case 3:
//				return tquat<T, Q>((m[0][1] - m[1][0]) * mult, (m[2][0] + m[0][2]) * mult, (m[1][2] + m[2][1]) * mult, biggestVal);
//			default: // Silence a -Wswitch-default warning in GCC. Should never actually get here. Assert is just for sanity.
//				assert(false);
//				return tquat<T, Q>(1, 0, 0, 0);
//		}
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> quat_cast(mat<4, 4, T, Q> const& m4)
//	{
//		return quat_cast(mat<3, 3, T, Q>(m4));
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER T angle(tquat<T, Q> const& x)
//	{
//		return acos(x.w) * static_cast<T>(2);
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER vec<3, T, Q> axis(tquat<T, Q> const& x)
//	{
//		T tmp1 = static_cast<T>(1) - x.w * x.w;
//		if(tmp1 <= static_cast<T>(0))
//			return vec<3, T, Q>(0, 0, 1);
//		T tmp2 = static_cast<T>(1) / sqrt(tmp1);
//		return vec<3, T, Q>(x.x * tmp2, x.y * tmp2, x.z * tmp2);
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER tquat<T, Q> angleAxis(T const& angle, vec<3, T, Q> const& v)
//	{
//		tquat<T, Q> Result;
//		
//		T const a(angle);
//		T const s = glm::sin(a * static_cast<T>(0.5));
//		
//		Result.w = glm::cos(a * static_cast<T>(0.5));
//		Result.x = v.x * s;
//		Result.y = v.y * s;
//		Result.z = v.z * s;
//		return Result;
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER vec<4, bool, Q> lessThan(tquat<T, Q> const& x, tquat<T, Q> const& y)
//	{
//		vec<4, bool, Q> Result;
//		for(length_t i = 0; i < x.length(); ++i)
//			Result[i] = x[i] < y[i];
//		return Result;
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER vec<4, bool, Q> lessThanEqual(tquat<T, Q> const& x, tquat<T, Q> const& y)
//	{
//		vec<4, bool, Q> Result;
//		for(length_t i = 0; i < x.length(); ++i)
//			Result[i] = x[i] <= y[i];
//		return Result;
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER vec<4, bool, Q> greaterThan(tquat<T, Q> const& x, tquat<T, Q> const& y)
//	{
//		vec<4, bool, Q> Result;
//		for(length_t i = 0; i < x.length(); ++i)
//			Result[i] = x[i] > y[i];
//		return Result;
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER vec<4, bool, Q> greaterThanEqual(tquat<T, Q> const& x, tquat<T, Q> const& y)
//	{
//		vec<4, bool, Q> Result;
//		for(length_t i = 0; i < x.length(); ++i)
//			Result[i] = x[i] >= y[i];
//		return Result;
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER vec<4, bool, Q> equal(tquat<T, Q> const& x, tquat<T, Q> const& y)
//	{
//		vec<4, bool, Q> Result;
//		for(length_t i = 0; i < x.length(); ++i)
//			Result[i] = detail::compute_equal<T>::call(x[i], y[i]);
//		return Result;
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER vec<4, bool, Q> notEqual(tquat<T, Q> const& x, tquat<T, Q> const& y)
//	{
//		vec<4, bool, Q> Result;
//		for(length_t i = 0; i < x.length(); ++i)
//			Result[i] = !detail::compute_equal<T>::call(x[i], y[i]);
//		return Result;
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER vec<4, bool, Q> isnan(tquat<T, Q> const& q)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'isnan' only accept floating-point inputs");
//		
//		return vec<4, bool, Q>(isnan(q.x), isnan(q.y), isnan(q.z), isnan(q.w));
//	}
//	
//	template<typename T, qualifier Q>
//	GLM_FUNC_QUALIFIER vec<4, bool, Q> isinf(tquat<T, Q> const& q)
//	{
//		GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'isinf' only accept floating-point inputs");
//		
//		return vec<4, bool, Q>(isinf(q.x), isinf(q.y), isinf(q.z), isinf(q.w));
//	}
//}//namespace glm
//
//
//
//
//
//#if GLM_ARCH & GLM_ARCH_SSE2_BIT
//
//namespace glm{
//	namespace detail
//	{
//		/*
//		 template<qualifier Q>
//		 struct compute_quat_mul<float, Q, true>
//		 {
//		 static tquat<float, Q> call(tquat<float, Q> const& q1, tquat<float, Q> const& q2)
//		 {
//		 // SSE2 STATS: 11 shuffle, 8 mul, 8 add
//		 // SSE4 STATS: 3 shuffle, 4 mul, 4 dpps
//		 
//		 __m128 const mul0 = _mm_mul_ps(q1.Data, _mm_shuffle_ps(q2.Data, q2.Data, _MM_SHUFFLE(0, 1, 2, 3)));
//		 __m128 const mul1 = _mm_mul_ps(q1.Data, _mm_shuffle_ps(q2.Data, q2.Data, _MM_SHUFFLE(1, 0, 3, 2)));
//		 __m128 const mul2 = _mm_mul_ps(q1.Data, _mm_shuffle_ps(q2.Data, q2.Data, _MM_SHUFFLE(2, 3, 0, 1)));
//		 __m128 const mul3 = _mm_mul_ps(q1.Data, q2.Data);
//		 
//		 #			if GLM_ARCH & GLM_ARCH_SSE41_BIT
//		 __m128 const add0 = _mm_dp_ps(mul0, _mm_set_ps(1.0f, -1.0f,  1.0f,  1.0f), 0xff);
//		 __m128 const add1 = _mm_dp_ps(mul1, _mm_set_ps(1.0f,  1.0f,  1.0f, -1.0f), 0xff);
//		 __m128 const add2 = _mm_dp_ps(mul2, _mm_set_ps(1.0f,  1.0f, -1.0f,  1.0f), 0xff);
//		 __m128 const add3 = _mm_dp_ps(mul3, _mm_set_ps(1.0f, -1.0f, -1.0f, -1.0f), 0xff);
//		 #			else
//		 __m128 const mul4 = _mm_mul_ps(mul0, _mm_set_ps(1.0f, -1.0f,  1.0f,  1.0f));
//		 __m128 const add0 = _mm_add_ps(mul0, _mm_movehl_ps(mul4, mul4));
//		 __m128 const add4 = _mm_add_ss(add0, _mm_shuffle_ps(add0, add0, 1));
//		 
//		 __m128 const mul5 = _mm_mul_ps(mul1, _mm_set_ps(1.0f,  1.0f,  1.0f, -1.0f));
//		 __m128 const add1 = _mm_add_ps(mul1, _mm_movehl_ps(mul5, mul5));
//		 __m128 const add5 = _mm_add_ss(add1, _mm_shuffle_ps(add1, add1, 1));
//		 
//		 __m128 const mul6 = _mm_mul_ps(mul2, _mm_set_ps(1.0f,  1.0f, -1.0f,  1.0f));
//		 __m128 const add2 = _mm_add_ps(mul6, _mm_movehl_ps(mul6, mul6));
//		 __m128 const add6 = _mm_add_ss(add2, _mm_shuffle_ps(add2, add2, 1));
//		 
//		 __m128 const mul7 = _mm_mul_ps(mul3, _mm_set_ps(1.0f, -1.0f, -1.0f, -1.0f));
//		 __m128 const add3 = _mm_add_ps(mul3, _mm_movehl_ps(mul7, mul7));
//		 __m128 const add7 = _mm_add_ss(add3, _mm_shuffle_ps(add3, add3, 1));
//		 #endif
//		 
//		 // This SIMD code is a politically correct way of doing this, but in every test I've tried it has been slower than
//		 // the final code below. I'll keep this here for reference - maybe somebody else can do something better...
//		 //
//		 //__m128 xxyy = _mm_shuffle_ps(add4, add5, _MM_SHUFFLE(0, 0, 0, 0));
//		 //__m128 zzww = _mm_shuffle_ps(add6, add7, _MM_SHUFFLE(0, 0, 0, 0));
//		 //
//		 //return _mm_shuffle_ps(xxyy, zzww, _MM_SHUFFLE(2, 0, 2, 0));
//		 
//		 tquat<float, Q> Result;
//		 _mm_store_ss(&Result.x, add4);
//		 _mm_store_ss(&Result.y, add5);
//		 _mm_store_ss(&Result.z, add6);
//		 _mm_store_ss(&Result.w, add7);
//		 return Result;
//		 }
//		 };
//		 */
//		
//		template<qualifier Q>
//		struct compute_dot<tquat<float, Q>, float, true>
//		{
//			static GLM_FUNC_QUALIFIER float call(tquat<float, Q> const& x, tquat<float, Q> const& y)
//			{
//				return _mm_cvtss_f32(glm_vec1_dot(x.data, y.data));
//			}
//		};
//		
//		template<qualifier Q>
//		struct compute_quat_add<float, Q, true>
//		{
//			static tquat<float, Q> call(tquat<float, Q> const& q, tquat<float, Q> const& p)
//			{
//				tquat<float, Q> Result;
//				Result.data = _mm_add_ps(q.data, p.data);
//				return Result;
//			}
//		};
//		
//#	if GLM_ARCH & GLM_ARCH_AVX_BIT
//		template<qualifier Q>
//		struct compute_quat_add<double, Q, true>
//		{
//			static tquat<double, Q> call(tquat<double, Q> const& a, tquat<double, Q> const& b)
//			{
//				tquat<double, Q> Result;
//				Result.data = _mm256_add_pd(a.data, b.data);
//				return Result;
//			}
//		};
//#	endif
//		
//		template<qualifier Q>
//		struct compute_quat_sub<float, Q, true>
//		{
//			static tquat<float, Q> call(tquat<float, Q> const& q, tquat<float, Q> const& p)
//			{
//				vec<4, float, Q> Result;
//				Result.data = _mm_sub_ps(q.data, p.data);
//				return Result;
//			}
//		};
//		
//#	if GLM_ARCH & GLM_ARCH_AVX_BIT
//		template<qualifier Q>
//		struct compute_quat_sub<double, Q, true>
//		{
//			static tquat<double, Q> call(tquat<double, Q> const& a, tquat<double, Q> const& b)
//			{
//				tquat<double, Q> Result;
//				Result.data = _mm256_sub_pd(a.data, b.data);
//				return Result;
//			}
//		};
//#	endif
//		
//		template<qualifier Q>
//		struct compute_quat_mul_scalar<float, Q, true>
//		{
//			static tquat<float, Q> call(tquat<float, Q> const& q, float s)
//			{
//				vec<4, float, Q> Result;
//				Result.data = _mm_mul_ps(q.data, _mm_set_ps1(s));
//				return Result;
//			}
//		};
//		
//#	if GLM_ARCH & GLM_ARCH_AVX_BIT
//		template<qualifier Q>
//		struct compute_quat_mul_scalar<double, Q, true>
//		{
//			static tquat<double, Q> call(tquat<double, Q> const& q, double s)
//			{
//				tquat<double, Q> Result;
//				Result.data = _mm256_mul_pd(q.data, _mm_set_ps1(s));
//				return Result;
//			}
//		};
//#	endif
//		
//		template<qualifier Q>
//		struct compute_quat_div_scalar<float, Q, true>
//		{
//			static tquat<float, Q> call(tquat<float, Q> const& q, float s)
//			{
//				vec<4, float, Q> Result;
//				Result.data = _mm_div_ps(q.data, _mm_set_ps1(s));
//				return Result;
//			}
//		};
//		
//#	if GLM_ARCH & GLM_ARCH_AVX_BIT
//		template<qualifier Q>
//		struct compute_quat_div_scalar<double, Q, true>
//		{
//			static tquat<double, Q> call(tquat<double, Q> const& q, double s)
//			{
//				tquat<double, Q> Result;
//				Result.data = _mm256_div_pd(q.data, _mm_set_ps1(s));
//				return Result;
//			}
//		};
//#	endif
//		
//		template<qualifier Q>
//		struct compute_quat_mul_vec4<float, Q, true>
//		{
//			static vec<4, float, Q> call(tquat<float, Q> const& q, vec<4, float, Q> const& v)
//			{
//				__m128 const q_wwww = _mm_shuffle_ps(q.data, q.data, _MM_SHUFFLE(3, 3, 3, 3));
//				__m128 const q_swp0 = _mm_shuffle_ps(q.data, q.data, _MM_SHUFFLE(3, 0, 2, 1));
//				__m128 const q_swp1 = _mm_shuffle_ps(q.data, q.data, _MM_SHUFFLE(3, 1, 0, 2));
//				__m128 const v_swp0 = _mm_shuffle_ps(v.data, v.data, _MM_SHUFFLE(3, 0, 2, 1));
//				__m128 const v_swp1 = _mm_shuffle_ps(v.data, v.data, _MM_SHUFFLE(3, 1, 0, 2));
//				
//				__m128 uv      = _mm_sub_ps(_mm_mul_ps(q_swp0, v_swp1), _mm_mul_ps(q_swp1, v_swp0));
//				__m128 uv_swp0 = _mm_shuffle_ps(uv, uv, _MM_SHUFFLE(3, 0, 2, 1));
//				__m128 uv_swp1 = _mm_shuffle_ps(uv, uv, _MM_SHUFFLE(3, 1, 0, 2));
//				__m128 uuv     = _mm_sub_ps(_mm_mul_ps(q_swp0, uv_swp1), _mm_mul_ps(q_swp1, uv_swp0));
//				
//				__m128 const two = _mm_set1_ps(2.0f);
//				uv  = _mm_mul_ps(uv, _mm_mul_ps(q_wwww, two));
//				uuv = _mm_mul_ps(uuv, two);
//				
//				vec<4, float, Q> Result;
//				Result.data = _mm_add_ps(v.Data, _mm_add_ps(uv, uuv));
//				return Result;
//			}
//		};
//	}//namespace detail
//}//namespace glm
//
//#endif//GLM_ARCH & GLM_ARCH_SSE2_BIT



