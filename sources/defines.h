// This software is MIT licensed (see LICENSE)

#pragma once


#define IMPLEMENT_VECTOR_OP_EQ(_dim) \
template<typename T> \
DLL_FUNCTION(void)	pack_t<_dim, T>::operator += (const pack_t<_dim, T>& v) \
{\
	for (int i = 0; i < _dim; i++) ((T*)this)[i] += v[i];\
}\
template<typename T> \
DLL_FUNCTION(void)	pack_t<_dim, T>::operator += (const T& v) \
{\
	for (int i = 0; i < _dim; i++) ((T*)this)[i] += v;\
}\
template<typename T> \
DLL_FUNCTION(void)	pack_t<_dim, T>::operator -= (const pack_t<_dim, T>& v) \
{\
	for (int i = 0; i < _dim; i++) ((T*)this)[i] -= v[i];\
}\
template<typename T> \
DLL_FUNCTION(void)	pack_t<_dim, T>::operator -= (const T& v) \
{\
	for (int i = 0; i < _dim; i++) ((T*)this)[i] -= v;\
}\
template<typename T> \
DLL_FUNCTION(void)	pack_t<_dim, T>::operator *= (const pack_t<_dim, T>& v) \
{\
	for (int i = 0; i < _dim; i++) ((T*)this)[i] *= v[i];\
}\
template<typename T> \
DLL_FUNCTION(void)	pack_t<_dim, T>::operator *= (const T& v) \
{\
	for (int i = 0; i < _dim; i++) ((T*)this)[i] *= v;\
}\
template<typename T> \
DLL_FUNCTION(void)	pack_t<_dim, T>::operator /= (const pack_t<_dim, T>& v) \
{\
	for (int i = 0; i < _dim; i++) ((T*)this)[i] /= v[i];\
}\
template<typename T> \
DLL_FUNCTION(void)	pack_t<_dim, T>::operator /= (const T& v)\
{\
	for (int i = 0; i < _dim; i++) ((T*)this)[i] /= v;\
}\



#define VECTOR_L_V_V_F_D_T(_F, _D, _T) \
	template<> DLL_APICALL typename length_t<_T>::type DLL_CALLING_CONVENTION _F<pack_t<_D, _T>>(pack_t<_D, _T> const& a, pack_t<_D, _T> const& b) \
	{ \
		return perform_##_F(a, b); \
	}
#define VECTOR_L_V_V_F_T(_F, _T) \
	VECTOR_L_V_V_F_D_T(_F, 2, _T) \
	VECTOR_L_V_V_F_D_T(_F, 3, _T) \
	VECTOR_L_V_V_F_D_T(_F, 4, _T)
#define VECTOR_L_V_V_F(_F) \
	VECTOR_L_V_V_F_T(_F, int8_t) \
	VECTOR_L_V_V_F_T(_F, uint8_t) \
	VECTOR_L_V_V_F_T(_F, int16_t) \
	VECTOR_L_V_V_F_T(_F, uint16_t) \
	VECTOR_L_V_V_F_T(_F, int32_t) \
	VECTOR_L_V_V_F_T(_F, uint32_t) \
	VECTOR_L_V_V_F_T(_F, int64_t) \
	VECTOR_L_V_V_F_T(_F, uint64_t) \
	VECTOR_L_V_V_F_T(_F, float) \
	VECTOR_L_V_V_F_T(_F, double)

#define VECTOR_L_V_F_D_T(_F, _D, _T) \
	template<> DLL_APICALL typename length_t<_T>::type DLL_CALLING_CONVENTION _F<pack_t<_D, _T>>(pack_t<_D, _T> const& a) \
	{ \
		return perform_##_F(a); \
	}
#define VECTOR_L_V_F_T(_F, _T) \
	VECTOR_L_V_F_D_T(_F, 2, _T) \
	VECTOR_L_V_F_D_T(_F, 3, _T) \
	VECTOR_L_V_F_D_T(_F, 4, _T)
#define VECTOR_L_V_F(_F) \
	VECTOR_L_V_F_T(_F, int8_t) \
	VECTOR_L_V_F_T(_F, uint8_t) \
	VECTOR_L_V_F_T(_F, int16_t) \
	VECTOR_L_V_F_T(_F, uint16_t) \
	VECTOR_L_V_F_T(_F, int32_t) \
	VECTOR_L_V_F_T(_F, uint32_t) \
	VECTOR_L_V_F_T(_F, int64_t) \
	VECTOR_L_V_F_T(_F, uint64_t) \
	VECTOR_L_V_F_T(_F, float) \
	VECTOR_L_V_F_T(_F, double)


#define VECTOR_IL_V_F_D_T(_F, _D, _T) \
	template<> DLL_APICALL typename lengthi_t<_T>::type DLL_CALLING_CONVENTION _F<pack_t<_D, _T>>(pack_t<_D, _T> const& a) \
	{ \
		return perform_##_F(a); \
	}
#define VECTOR_IL_V_F_T(_F, _T) \
	VECTOR_IL_V_F_D_T(_F, 2, _T) \
	VECTOR_IL_V_F_D_T(_F, 3, _T) \
	VECTOR_IL_V_F_D_T(_F, 4, _T)
#define VECTOR_IL_V_F(_F) \
	VECTOR_IL_V_F_T(_F, int8_t) \
	VECTOR_IL_V_F_T(_F, uint8_t) \
	VECTOR_IL_V_F_T(_F, int16_t) \
	VECTOR_IL_V_F_T(_F, uint16_t) \
	VECTOR_IL_V_F_T(_F, int32_t) \
	VECTOR_IL_V_F_T(_F, uint32_t) \
	VECTOR_IL_V_F_T(_F, int64_t) \
	VECTOR_IL_V_F_T(_F, uint64_t) \
	VECTOR_IL_V_F_T(_F, float) \
	VECTOR_IL_V_F_T(_F, double)



#define VECTOR_V_V_F_D_T(_F, _D, _T) \
	template<> DLL_APICALL pack_t<_D, _T> DLL_CALLING_CONVENTION _F<pack_t<_D, _T>>(pack_t<_D, _T> const& a) \
	{ \
		return perform_##_F(a); \
	}
#define VECTOR_V_V_F_T(_F, _T) \
	VECTOR_V_V_F_D_T(_F, 2, _T) \
	VECTOR_V_V_F_D_T(_F, 3, _T) \
	VECTOR_V_V_F_D_T(_F, 4, _T)
#define VECTOR_V_V_F(_F) \
	VECTOR_V_V_F_T(_F, int8_t) \
	VECTOR_V_V_F_T(_F, uint8_t) \
	VECTOR_V_V_F_T(_F, int16_t) \
	VECTOR_V_V_F_T(_F, uint16_t) \
	VECTOR_V_V_F_T(_F, int32_t) \
	VECTOR_V_V_F_T(_F, uint32_t) \
	VECTOR_V_V_F_T(_F, int64_t) \
	VECTOR_V_V_F_T(_F, uint64_t) \
	VECTOR_V_V_F_T(_F, float) \
	VECTOR_V_V_F_T(_F, double)

#define VECTOR_V_V_V_F_D_T(_F, _D, _T) \
	template<> DLL_APICALL pack_t<_D, _T> DLL_CALLING_CONVENTION _F<pack_t<_D, _T>>(pack_t<_D, _T> const& a, pack_t<_D, _T> const& b) \
	{ \
		return perform_##_F(a, b); \
	}
#define VECTOR_V_V_V_F_T(_F, _T) \
	VECTOR_V_V_V_F_D_T(_F, 2, _T) \
	VECTOR_V_V_V_F_D_T(_F, 3, _T) \
	VECTOR_V_V_V_F_D_T(_F, 4, _T)
#define VECTOR_V_V_V_F(_F) \
	VECTOR_V_V_V_F_T(_F, int8_t) \
	VECTOR_V_V_V_F_T(_F, uint8_t) \
	VECTOR_V_V_V_F_T(_F, int16_t) \
	VECTOR_V_V_V_F_T(_F, uint16_t) \
	VECTOR_V_V_V_F_T(_F, int32_t) \
	VECTOR_V_V_V_F_T(_F, uint32_t) \
	VECTOR_V_V_V_F_T(_F, int64_t) \
	VECTOR_V_V_V_F_T(_F, uint64_t) \
	VECTOR_V_V_V_F_T(_F, float) \
	VECTOR_V_V_V_F_T(_F, double)


//simd::pack_t<3, double, (simd::qualifier)2>
//simd::operator-<3, double, (simd::qualifier)2>
//(simd::pack_t<3, double, (simd::qualifier)2> const&, simd::pack_t<3, double, (simd::qualifier)2> const&)


#define VECTOR_V_V_V_O_F_D_T(_O, _F, _D, _T) \
	template<> DLL_APICALL pack_t<_D, _T> DLL_CALLING_CONVENTION operator _O <_D, _T>(pack_t<_D, _T> const& a, pack_t<_D, _T> const& b) \
	{ \
		return perform_operator_##_F(a, b); \
	} \
	template<> DLL_APICALL pack_t<_D, _T> DLL_CALLING_CONVENTION operator _O <_D, _T>(_T const& a, pack_t<_D, _T> const& b) \
	{ \
		return perform_operator_##_F(a, b); \
	} \
	template<> DLL_APICALL pack_t<_D, _T> DLL_CALLING_CONVENTION operator _O <_D, _T>(pack_t<_D, _T> const& a, _T const& b) \
	{ \
		return perform_operator_##_F(a, b); \
	}
#define VECTOR_V_V_V_O_F_T(_O, _F, _T) \
	VECTOR_V_V_V_O_F_D_T(_O, _F, 2, _T) \
	VECTOR_V_V_V_O_F_D_T(_O, _F, 3, _T) \
	VECTOR_V_V_V_O_F_D_T(_O, _F, 4, _T)
#define VECTOR_V_V_V_O_F(_O, _F) \
	VECTOR_V_V_V_O_F_T(_O, _F, int8_t) \
	VECTOR_V_V_V_O_F_T(_O, _F, uint8_t) \
	VECTOR_V_V_V_O_F_T(_O, _F, int16_t) \
	VECTOR_V_V_V_O_F_T(_O, _F, uint16_t) \
	VECTOR_V_V_V_O_F_T(_O, _F, int32_t) \
	VECTOR_V_V_V_O_F_T(_O, _F, uint32_t) \
	VECTOR_V_V_V_O_F_T(_O, _F, int64_t) \
	VECTOR_V_V_V_O_F_T(_O, _F, uint64_t) \
	VECTOR_V_V_V_O_F_T(_O, _F, float) \
	VECTOR_V_V_V_O_F_T(_O, _F, double)


#define MATRIX_M_M_M_O_F_C_R_T(_O, _F, _C, _R, _T) \
template<>	DLL_APICALL	mat_t<_C, _R, _T>	DLL_CALLING_CONVENTION	operator _O <_C, _R, _T>(mat_t<_C, _R, _T> const& a, mat_t<_C, _R, _T> const& b) \
{ \
	return perform_operator_##_F(a, b); \
}
#define MATRIX_M_M_M_O_F_T(_O, _F, _T) \
	MATRIX_M_M_M_O_F_C_R_T(_O, _F, 2, 2, _T) \
	MATRIX_M_M_M_O_F_C_R_T(_O, _F, 3, 3, _T) \
	MATRIX_M_M_M_O_F_C_R_T(_O, _F, 4, 4, _T)
#define MATRIX_M_M_M_O_F(_O, _F) \
	MATRIX_M_M_M_O_F_T(_O, _F, int8_t) \
	MATRIX_M_M_M_O_F_T(_O, _F, uint8_t) \
	MATRIX_M_M_M_O_F_T(_O, _F, int16_t) \
	MATRIX_M_M_M_O_F_T(_O, _F, uint16_t) \
	MATRIX_M_M_M_O_F_T(_O, _F, int32_t) \
	MATRIX_M_M_M_O_F_T(_O, _F, uint32_t) \
	MATRIX_M_M_M_O_F_T(_O, _F, int64_t) \
	MATRIX_M_M_M_O_F_T(_O, _F, uint64_t) \
	MATRIX_M_M_M_O_F_T(_O, _F, float) \
	MATRIX_M_M_M_O_F_T(_O, _F, double)

#define MATRIX_M_M_S_O_F_C_R_T(_O, _F, _C, _R, _T) \
template<>	DLL_APICALL	mat_t<_C, _R, _T>	DLL_CALLING_CONVENTION	operator _O <_C, _R, _T>(mat_t<_C, _R, _T> const& a, _T const& b) \
{ \
	return perform_operator_##_F(a, b); \
} \
template<>	DLL_APICALL	mat_t<_C, _R, _T>	DLL_CALLING_CONVENTION	operator _O <_C, _R, _T>(_T const& a, mat_t<_C, _R, _T> const& b) \
{ \
	return perform_operator_##_F(a, b); \
}
#define MATRIX_M_M_S_O_F_T(_O, _F, _T) \
	MATRIX_M_M_S_O_F_C_R_T(_O, _F, 2, 2, _T) \
	MATRIX_M_M_S_O_F_C_R_T(_O, _F, 3, 3, _T) \
	MATRIX_M_M_S_O_F_C_R_T(_O, _F, 4, 4, _T)
#define MATRIX_M_M_S_O_F(_O, _F) \
	MATRIX_M_M_S_O_F_T(_O, _F, int8_t) \
	MATRIX_M_M_S_O_F_T(_O, _F, uint8_t) \
	MATRIX_M_M_S_O_F_T(_O, _F, int16_t) \
	MATRIX_M_M_S_O_F_T(_O, _F, uint16_t) \
	MATRIX_M_M_S_O_F_T(_O, _F, int32_t) \
	MATRIX_M_M_S_O_F_T(_O, _F, uint32_t) \
	MATRIX_M_M_S_O_F_T(_O, _F, int64_t) \
	MATRIX_M_M_S_O_F_T(_O, _F, uint64_t) \
	MATRIX_M_M_S_O_F_T(_O, _F, float) \
	MATRIX_M_M_S_O_F_T(_O, _F, double)



#define MATRIX_M_M_F_C_R_T(_F, _C, _R, _T) \
	template<> DLL_APICALL mat_t<_C, _R, _T> DLL_CALLING_CONVENTION _F <_C, _R, _T> (mat_t<_C, _R, _T> const& a) \
	{ \
		return perform_##_F(a); \
	}
#define MATRIX_M_M_F_T(_F, _T) \
	MATRIX_M_M_F_C_R_T(_F, 2, 2, _T) \
	MATRIX_M_M_F_C_R_T(_F, 3, 3, _T) \
	MATRIX_M_M_F_C_R_T(_F, 4, 4, _T)
#define MATRIX_M_M_F(_F) \
	MATRIX_M_M_F_T(_F, int8_t) \
	MATRIX_M_M_F_T(_F, uint8_t) \
	MATRIX_M_M_F_T(_F, int16_t) \
	MATRIX_M_M_F_T(_F, uint16_t) \
	MATRIX_M_M_F_T(_F, int32_t) \
	MATRIX_M_M_F_T(_F, uint32_t) \
	MATRIX_M_M_F_T(_F, int64_t) \
	MATRIX_M_M_F_T(_F, uint64_t) \
	MATRIX_M_M_F_T(_F, float) \
	MATRIX_M_M_F_T(_F, double)


#define TYPE_M_M_F_C_R_T(_F, _C, _R, _T) \
	template<> DLL_APICALL mat_t<_C, _R, _T> DLL_CALLING_CONVENTION _F <mat_t<_C, _R, _T>> (mat_t<_C, _R, _T> const& a) \
	{ \
		return perform_##_F(a); \
	}
#define TYPE_M_M_F_T(_F, _T) \
	TYPE_M_M_F_C_R_T(_F, 2, 2, _T) \
	TYPE_M_M_F_C_R_T(_F, 3, 3, _T) \
	TYPE_M_M_F_C_R_T(_F, 4, 4, _T)
#define TYPE_M_M_F(_F) \
	TYPE_M_M_F_T(_F, int8_t) \
	TYPE_M_M_F_T(_F, uint8_t) \
	TYPE_M_M_F_T(_F, int16_t) \
	TYPE_M_M_F_T(_F, uint16_t) \
	TYPE_M_M_F_T(_F, int32_t) \
	TYPE_M_M_F_T(_F, uint32_t) \
	TYPE_M_M_F_T(_F, int64_t) \
	TYPE_M_M_F_T(_F, uint64_t) \
	TYPE_M_M_F_T(_F, float) \
	TYPE_M_M_F_T(_F, double)
