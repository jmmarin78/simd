// This software is MIT licensed (see LICENSE)

namespace simd
{
	// -- Unary bit operators --
	template<typename T>	FORCE_INLINE quaternion_t<T> operator + (quaternion_t<T> const& q);
	template<typename T>	FORCE_INLINE quaternion_t<T> operator - (quaternion_t<T> const& q);
	// -- Binary operators --
	template<typename T>	FORCE_INLINE quaternion_t<T> operator + (quaternion_t<T> const& q, quaternion_t<T> const& p);
	template<typename T>	FORCE_INLINE quaternion_t<T> operator - (quaternion_t<T> const& q, quaternion_t<T> const& p);
	template<typename T>	FORCE_INLINE quaternion_t<T> operator * (quaternion_t<T> const& q, quaternion_t<T> const& p);
	template<typename T>	FORCE_INLINE pack_t<3, T> operator * (quaternion_t<T> const& q, pack_t<3, T> const& v);
	template<typename T>	FORCE_INLINE pack_t<3, T> operator * (pack_t<3, T> const& v, quaternion_t<T> const& q);
	template<typename T>	FORCE_INLINE pack_t<4, T> operator * (quaternion_t<T> const& q, pack_t<4, T> const& v);
	template<typename T>	FORCE_INLINE pack_t<4, T> operator * (pack_t<4, T> const& v, quaternion_t<T> const& q);
	template<typename T>	FORCE_INLINE quaternion_t<T> operator * (quaternion_t<T> const& q, T const& s);
	template<typename T>	FORCE_INLINE quaternion_t<T> operator * (T const& s, quaternion_t<T> const& q);
	template<typename T>	FORCE_INLINE quaternion_t<T> operator / (quaternion_t<T> const& q, T const& s);
	// -- Boolean operators --
	template<typename T>	FORCE_INLINE bool operator==(quaternion_t<T> const& q1, quaternion_t<T> const& q2);
	template<typename T>	FORCE_INLINE bool operator!=(quaternion_t<T> const& q1, quaternion_t<T> const& q2);
	/// Returns the length of the quaternion.
	///
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE T length(quaternion_t<T> const& q);
	/// Returns the normalized quaternion.
	///
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE quaternion_t<T> normalize(quaternion_t<T> const& q);
	/// Returns dot product of q1 and q2, i.e., q1[0] * q2[0] + q1[1] * q2[1] + ...
	///
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE T dot(quaternion_t<T> const& x, quaternion_t<T> const& y);
	
	/// Spherical linear interpolation of two quaternions.
	/// The interpolation is oriented and the rotation is performed at constant speed.
	/// For short path spherical linear interpolation, use the slerp function.
	///
	/// @param x A quaternion
	/// @param y A quaternion
	/// @param a Interpolation factor. The interpolation is defined beyond the range [0, 1].
	/// @tparam T Floating-point scalar types.
	///
	/// @see - slerp(quaternion_t<T> const& x, quaternion_t<T> const& y, T const& a)
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE quaternion_t<T> mix(quaternion_t<T> const& x, quaternion_t<T> const& y, T a);
	
	/// Linear interpolation of two quaternions.
	/// The interpolation is oriented.
	///
	/// @param x A quaternion
	/// @param y A quaternion
	/// @param a Interpolation factor. The interpolation is defined in the range [0, 1].
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE quaternion_t<T> lerp(quaternion_t<T> const& x, quaternion_t<T> const& y, T a);
	
	/// Spherical linear interpolation of two quaternions.
	/// The interpolation always take the short path and the rotation is performed at constant speed.
	///
	/// @param x A quaternion
	/// @param y A quaternion
	/// @param a Interpolation factor. The interpolation is defined beyond the range [0, 1].
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE quaternion_t<T> slerp(quaternion_t<T> const& x, quaternion_t<T> const& y, T a);
	
	/// Returns the q conjugate.
	///
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE quaternion_t<T> conjugate(quaternion_t<T> const& q);
	
	/// Returns the q inverse.
	///
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE quaternion_t<T> inverse(quaternion_t<T> const& q);
	
	/// Rotates a quaternion from a vector of 3 components axis and an angle.
	///
	/// @param q Source orientation
	/// @param angle Angle expressed in radians.
	/// @param axis Axis of the rotation
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE quaternion_t<T> rotate(quaternion_t<T> const& q, T const& angle, pack_t<3, T> const& axis);
	
	/// Returns euler angles, pitch as x, yaw as y, roll as z.
	/// The result is expressed in radians.
	///
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE pack_t<3, T> eulerAngles(quaternion_t<T> const& x);
	
	/// Returns roll value of euler angles expressed in radians.
	///
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE T roll(quaternion_t<T> const& x);
	
	/// Returns pitch value of euler angles expressed in radians.
	///
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE T pitch(quaternion_t<T> const& x);
	
	/// Returns yaw value of euler angles expressed in radians.
	///
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE T yaw(quaternion_t<T> const& x);
	
	/// Converts a quaternion to a 3 * 3 matrix.
	///
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE mat_t<3, 3, T> mat3_cast(quaternion_t<T> const& x);
	
	/// Converts a quaternion to a 4 * 4 matrix.
	///
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE mat_t<4, 4, T> mat4_cast(quaternion_t<T> const& x);
	
	/// Converts a pure rotation 3 * 3 matrix to a quaternion.
	///
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE quaternion_t<T> quat_cast(mat_t<3, 3, T> const& x);
	
	/// Converts a pure rotation 4 * 4 matrix to a quaternion.
	///
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE quaternion_t<T> quat_cast(mat_t<4, 4, T> const& x);
	
	/// Returns the quaternion rotation angle.
	///
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE T angle(quaternion_t<T> const& x);
	
	/// Returns the q rotation axis.
	///
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE pack_t<3, T> axis(quaternion_t<T> const& x);
	
	/// Build a quaternion from an angle and a normalized axis.
	///
	/// @param angle Angle expressed in radians.
	/// @param axis Axis of the quaternion, must be normalized.
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE quaternion_t<T> angleAxis(T const& angle, pack_t<3, T> const& axis);
	
	/// Returns the component-wise comparison result of x < y.
	///
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE pack_t<4, bool> lessThan(quaternion_t<T> const& x, quaternion_t<T> const& y);
	
	/// Returns the component-wise comparison of result x <= y.
	///
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE pack_t<4, bool> lessThanEqual(quaternion_t<T> const& x, quaternion_t<T> const& y);
	
	/// Returns the component-wise comparison of result x > y.
	///
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE pack_t<4, bool> greaterThan(quaternion_t<T> const& x, quaternion_t<T> const& y);
	
	/// Returns the component-wise comparison of result x >= y.
	///
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE pack_t<4, bool> greaterThanEqual(quaternion_t<T> const& x, quaternion_t<T> const& y);
	
	/// Returns the component-wise comparison of result x == y.
	///
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE pack_t<4, bool> equal(quaternion_t<T> const& x, quaternion_t<T> const& y);
	
	/// Returns the component-wise comparison of result x != y.
	///
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE pack_t<4, bool> notEqual(quaternion_t<T> const& x, quaternion_t<T> const& y);
	
	/// Returns true if x holds a NaN (not a number)
	/// representation in the underlying implementation's set of
	/// floating point representations. Returns false otherwise,
	/// including for implementations with no NaN
	/// representations.
	///
	/// /!\ When using compiler fast math, this function may fail.
	///
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE pack_t<4, bool> isnan(quaternion_t<T> const& x);
	
	/// Returns true if x holds a positive infinity or negative
	/// infinity representation in the underlying implementation's
	/// set of floating point representations. Returns false
	/// otherwise, including for implementations with no infinity
	/// representations.
	///
	/// @tparam T Floating-point scalar types.
	///
	/// @see gtc_quaternion
	template<typename T>
	FORCE_INLINE pack_t<4, bool> isinf(quaternion_t<T> const& x);
	
	//namespace glm
	//{
	//	/// @addtogroup gtx_quaternion
	//	/// @{
	//
	//	/// Create an identity quaternion.
	//	///
	//	/// @see gtx_quaternion
	//	template<typename T>
	//	GLM_FUNC_DECL quaternion_t<T> quat_identity();
	//
	//	/// Compute a cross product between a quaternion and a vector.
	//	///
	//	/// @see gtx_quaternion
	//	template<typename T>
	//	GLM_FUNC_DECL pack_t<3, T> cross(
	//									 quaternion_t<T> const& q,
	//									 pack_t<3, T> const& v);
	//
	//	//! Compute a cross product between a vector and a quaternion.
	//	///
	//	/// @see gtx_quaternion
	//	template<typename T>
	//	GLM_FUNC_DECL pack_t<3, T> cross(
	//									 pack_t<3, T> const& v,
	//									 quaternion_t<T> const& q);
	//
	//	//! Compute a point on a path according squad equation.
	//	//! q1 and q2 are control points; s1 and s2 are intermediate control points.
	//	///
	//	/// @see gtx_quaternion
	//	template<typename T>
	//	GLM_FUNC_DECL quaternion_t<T> squad(
	//									quaternion_t<T> const& q1,
	//									quaternion_t<T> const& q2,
	//									quaternion_t<T> const& s1,
	//									quaternion_t<T> const& s2,
	//									T const& h);
	//
	//	//! Returns an intermediate control point for squad interpolation.
	//	///
	//	/// @see gtx_quaternion
	//	template<typename T>
	//	GLM_FUNC_DECL quaternion_t<T> intermediate(
	//										   quaternion_t<T> const& prev,
	//										   quaternion_t<T> const& curr,
	//										   quaternion_t<T> const& next);
	//
	//	//! Returns a exp of a quaternion.
	//	///
	//	/// @see gtx_quaternion
	//	template<typename T>
	//	GLM_FUNC_DECL quaternion_t<T> exp(
	//								  quaternion_t<T> const& q);
	//
	//	//! Returns a log of a quaternion.
	//	///
	//	/// @see gtx_quaternion
	//	template<typename T>
	//	GLM_FUNC_DECL quaternion_t<T> log(
	//								  quaternion_t<T> const& q);
	//
	//	/// Returns x raised to the y power.
	//	///
	//	/// @see gtx_quaternion
	//	template<typename T>
	//	GLM_FUNC_DECL quaternion_t<T> pow(
	//								  quaternion_t<T> const& x,
	//								  T const& y);
	//
	//	//! Returns quarternion square root.
	//	///
	//	/// @see gtx_quaternion
	//	//template<typename T>
	//	//quaternion_t<T> sqrt(
	//	//	quaternion_t<T> const& q);
	//
	//	//! Rotates a 3 components vector by a quaternion.
	//	///
	//	/// @see gtx_quaternion
	//	template<typename T>
	//	GLM_FUNC_DECL pack_t<3, T> rotate(
	//									  quaternion_t<T> const& q,
	//									  pack_t<3, T> const& v);
	//
	//	/// Rotates a 4 components vector by a quaternion.
	//	///
	//	/// @see gtx_quaternion
	//	template<typename T>
	//	GLM_FUNC_DECL pack_t<4, T> rotate(
	//									  quaternion_t<T> const& q,
	//									  pack_t<4, T> const& v);
	//
	//	/// Extract the real component of a quaternion.
	//	///
	//	/// @see gtx_quaternion
	//	template<typename T>
	//	GLM_FUNC_DECL T extractRealComponent(
	//										 quaternion_t<T> const& q);
	//
	//	/// Converts a quaternion to a 3 * 3 matrix.
	//	///
	//	/// @see gtx_quaternion
	//	template<typename T>
	//	GLM_FUNC_DECL mat_t<3, 3, T> toMat3(
	//										 quaternion_t<T> const& x){return mat3_cast(x);}
	//
	//	/// Converts a quaternion to a 4 * 4 matrix.
	//	///
	//	/// @see gtx_quaternion
	//	template<typename T>
	//	GLM_FUNC_DECL mat_t<4, 4, T> toMat4(
	//										 quaternion_t<T> const& x){return mat4_cast(x);}
	//
	//	/// Converts a 3 * 3 matrix to a quaternion.
	//	///
	//	/// @see gtx_quaternion
	//	template<typename T>
	//	GLM_FUNC_DECL quaternion_t<T> toQuat(
	//									 mat_t<3, 3, T> const& x){return quat_cast(x);}
	//
	//	/// Converts a 4 * 4 matrix to a quaternion.
	//	///
	//	/// @see gtx_quaternion
	//	template<typename T>
	//	GLM_FUNC_DECL quaternion_t<T> toQuat(
	//									 mat_t<4, 4, T> const& x){return quat_cast(x);}
	//
	//	/// Quaternion interpolation using the rotation short path.
	//	///
	//	/// @see gtx_quaternion
	//	template<typename T>
	//	GLM_FUNC_DECL quaternion_t<T> shortMix(
	//									   quaternion_t<T> const& x,
	//									   quaternion_t<T> const& y,
	//									   T const& a);
	//
	//	/// Quaternion normalized linear interpolation.
	//	///
	//	/// @see gtx_quaternion
	//	template<typename T>
	//	GLM_FUNC_DECL quaternion_t<T> fastMix(
	//									  quaternion_t<T> const& x,
	//									  quaternion_t<T> const& y,
	//									  T const& a);
	//
	//	/// Compute the rotation between two vectors.
	//	/// param orig vector, needs to be normalized
	//	/// param dest vector, needs to be normalized
	//	///
	//	/// @see gtx_quaternion
	//	template<typename T>
	//	GLM_FUNC_DECL quaternion_t<T> rotation(
	//									   pack_t<3, T> const& orig,
	//									   pack_t<3, T> const& dest);
	//
	//	/// Build a look at quaternion based on the default handedness.
	//	///
	//	/// @param direction Desired forward direction. Needs to be normalized.
	//	/// @param up Up vector, how the camera is oriented. Typically (0, 1, 0).
	//	template<typename T>
	//	GLM_FUNC_DECL quaternion_t<T> quatLookAt(
	//										 pack_t<3, T> const& direction,
	//										 pack_t<3, T> const& up);
	//
	//	/// Build a right-handed look at quaternion.
	//	///
	//	/// @param direction Desired forward direction onto which the -z-axis gets mapped. Needs to be normalized.
	//	/// @param up Up vector, how the camera is oriented. Typically (0, 1, 0).
	//	template<typename T>
	//	GLM_FUNC_DECL quaternion_t<T> quatLookAtRH(
	//										   pack_t<3, T> const& direction,
	//										   pack_t<3, T> const& up);
	//
	//	/// Build a left-handed look at quaternion.
	//	///
	//	/// @param direction Desired forward direction onto which the +z-axis gets mapped. Needs to be normalized.
	//	/// @param up Up vector, how the camera is oriented. Typically (0, 1, 0).
	//	template<typename T>
	//	GLM_FUNC_DECL quaternion_t<T> quatLookAtLH(
	//										   pack_t<3, T> const& direction,
	//										   pack_t<3, T> const& up);
	//
	//	/// Returns the squared length of x.
	//	///
	//	/// @see gtx_quaternion
	//	template<typename T>
	//	GLM_FUNC_DECL T length2(quaternion_t<T> const& q);
	//
	//	/// @}
	//}//namespace glm
	
}
