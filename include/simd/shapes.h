// This software is MIT licensed (see LICENSE)

#pragma once

#include <simd/core.h>

namespace Maths
{
#define SHAPES_FUNC_DEF
#define SHAPES_FUNC_INL	FORCE_INLINE
#define SHAPES_FUNC_SIN	static FORCE_INLINE
#define SHAPES_FUNC_STT	static
#define SHAPES_FUNC_TM	template<typename O>
#define SHAPES_FUNC_TMI	template<typename O> FORCE_INLINE

	template<typename TD, typename TS>
	void	shapes_operator_assign(TD* o, const TS&);
}

// Boxes
namespace simd
{
	template <int D, typename T>
	struct aabox_t
	{
		simd::pack_t<D, T>	min, max;

		FORCE_INLINE	aabox_t(){}
		template<int D2, typename T2>	inline	aabox_t(const aabox_t<D2, T2>& i);
		aabox_t(const pack_t<D, T>& AMinValue, const pack_t<D, T>& AMaxValue);
		aabox_t(const T& AMinValueX, const T& AMinValueY, const T& AMaxValueX, const T& AMaxValueY);
		aabox_t(const pack_t<D, T>& _Center, const T& _Radius);
		aabox_t(const T& AMinValue, const T& AMaxValue);
		aabox_t(const pack_t<D, T>* Points, int PointsCount);
		aabox_t(const aabox_t<D, T> Boxes[], int BoxCount);
		
		SHAPES_FUNC_SIN	auto	dimension() -> int						{return D;}
		SHAPES_FUNC_SIN	auto	invalid_aabox() -> aabox_t<D, T>		{aabox_t<D, T> ret; ret.invalidate(); return ret;}

		DLL_FUNCTION(void)	invalidate();
		DLL_FUNCTION(bool)	is_valid() const;
		DLL_FUNCTION(auto)	center() const -> pack_t<D, T>;
		DLL_FUNCTION(auto)	size() const -> pack_t<D, T>;
		DLL_FUNCTION(auto)	width() const -> T;
		DLL_FUNCTION(auto)	height() const -> T;
		DLL_FUNCTION(auto)	depth() const -> T;
		DLL_FUNCTION(auto)	inner_radius() const -> T;
		DLL_FUNCTION(auto)	outer_radius() const -> T;
		DLL_FUNCTION(auto)	vertex(bool AMaxValueX, bool AMaxValueY, bool AMaxValueZ) const -> pack_t<D, T>;
		DLL_FUNCTION(auto)	vertex(int i) const -> pack_t<D, T>;
		DLL_FUNCTION(void)	vertices(pack_t<D, T>* o) const;
		
		DLL_FUNCTION(void)	set(const aabox_t<D, T>& _Box, int iSplit);
		DLL_FUNCTION(void)	set(const pack_t<D, T>* _Points, int _PointsCount);
		DLL_FUNCTION(void)	set(const pack_t<D, T>& _Value);
		DLL_FUNCTION(void)	expand(const pack_t<D, T>*, int iCount);
		DLL_FUNCTION(void)	expand(const pack_t<D, T>&);
		DLL_FUNCTION(void)	expand(const aabox_t<D, T>&);
		DLL_FUNCTION(bool)	contains(const pack_t<D, T>& _Value) const;
		DLL_FUNCTION(bool)	intersect(T* hitDist, const pack_t<D, T>* origPt, const pack_t<D, T>* dir);
		DLL_FUNCTION(int)	dominant_axis() const;
	};
	
	template <int D, typename T, bool A, int VC>
	struct TBox
	{
		pack_t<D, T>	vertex[VC];

		SHAPES_FUNC_INL	TBox(){};
	};
}

// Plane, triangle, sphere
namespace simd
{
	template <int D, typename T>
	struct TPlane
	{
	private:
		pack_t<D + 1, T>	mData;
	public:
		SHAPES_FUNC_INL	auto	axis() -> pack_t<3, T>&						{return *(pack_t<3, T>*)&mData;}
		SHAPES_FUNC_INL	auto	axis() const -> const pack_t<3, T>&			{return *(const pack_t<3, T>*)&mData;}
		SHAPES_FUNC_INL	auto	origin_distance() -> T&								{return mData.w;}
		SHAPES_FUNC_INL	auto	origin_distance() const -> const T&					{return mData.w;}
	};
	
	template <int D, typename T>
	struct TSphere
	{
	private:
		pack_t<D + 1, T>	mData;
	public:
		SHAPES_FUNC_INL	auto	center() -> pack_t<3, T>&					{return *(pack_t<3, T>*)&mData;}
		SHAPES_FUNC_INL	auto	center() const -> const	pack_t<3, T>&		{return *(const pack_t<3, T>*)&mData;}
		SHAPES_FUNC_INL	auto	radius() -> T&										{return mData.w;}
		SHAPES_FUNC_INL	auto	radius() const -> const T&							{return mData.w;}
	};
	
	template <int D, typename T>
	struct TTriangle
	{
		pack_t<D, T>	vertex[3];
	};
}

// Lines
namespace simd
{
	template<int D, typename T, bool _HasVector>
	struct TLinearShapeBase
	{
	};
	
	template<int D, typename T>
	struct TLinearShapeBase<D, T, true>
	{
		typedef pack_t<D, T>	vector_t;
		pack_t<D, T>			Vertex;
		pack_t<D, T>			Vector;

		SHAPES_FUNC_INL	void	set(const TLinearShapeBase& i)							{Vector=i.Vector;Vertex=i.Vertex;}
		SHAPES_FUNC_INL	auto	vertex(int i) const -> pack_t<D, T>				{if(i==0)return Vertex;else return Vertex+Vector;}
		SHAPES_FUNC_INL	auto	vector() const -> const	pack_t<D, T>&				{return Vector;}
		SHAPES_FUNC_INL	auto	origin() const -> const	pack_t<D, T>&				{return Vertex;}
		SHAPES_FUNC_INL	auto	extrem() const -> pack_t<D, T>					{return Vertex+Vector;}
		SHAPES_FUNC_INL	auto	component0() -> pack_t<D, T>&						{return Vertex;}
		SHAPES_FUNC_INL	auto	component1() -> pack_t<D, T>&						{return Vector;}
		SHAPES_FUNC_INL	auto	component0() const -> const	pack_t<D, T>&			{return Vertex;}
		SHAPES_FUNC_INL	auto	component1() const -> const	pack_t<D, T>&			{return Vector;}
	};
	
	template<int D, typename T>
	struct TLinearShapeBase<D, T, false>
	{
		typedef	pack_t<D, T>	vector_t;
		pack_t<D, T>			Vertex[2];
		
		SHAPES_FUNC_INL	void	set(const TLinearShapeBase& i)					{Vertex[0]=i.Vertex[0];Vertex[1]=i.Vertex[1];}
		SHAPES_FUNC_INL	auto	vertex(int i) const->const	pack_t<D, T>&		{assert(0<=i && i<=1);return Vertex[i];}
		SHAPES_FUNC_INL	auto	vector() const->pack_t<D, T>					{return Vertex[1]-Vertex[0];}
		SHAPES_FUNC_INL	auto	origin() const->const pack_t<D, T>&				{return Vertex[0];}
		SHAPES_FUNC_INL	auto	extrem() const->const pack_t<D, T>&				{return Vertex[1];}
		SHAPES_FUNC_INL	auto	component0()->pack_t<D, T>&						{return Vertex[0];}
		SHAPES_FUNC_INL	auto	component1()->pack_t<D, T>&						{return Vertex[1];}
		SHAPES_FUNC_INL	auto	component0() const->const pack_t<D, T>&			{return Vertex[0];}
		SHAPES_FUNC_INL	auto	component1() const->const pack_t<D, T>&			{return Vertex[1];}
	};
	
	template<int D, typename T, bool _ForwardsInfinite, bool _BackwardsInfinite>
	struct TLinearShape : public TLinearShapeBase<D, T, _ForwardsInfinite || _BackwardsInfinite>
	{
	};
	
	
	
	template<int D, typename T>
	struct TSegment : public TLinearShape<D, T, false, false>
	{
		FORCE_INLINE	TSegment(){};
	};
	
	template<int D, typename T>
	struct TRay : public TLinearShape<D, T, true, false>
	{
		FORCE_INLINE	TRay(){};
	};
	
	template<int D, typename T>
	struct TLine : public TLinearShape<D, T, true, true>
	{
		FORCE_INLINE	TLine(){};
	};
}

// Volumes
namespace simd
{
	template<typename T>
	class frustum_t
	{
		enum mode_t {PERSPECTIVE_FOVY_ASPECT};
		mode_t	m_mode;
		T	m_fovy, m_aspect, m_zNear, m_zFar;
		friend 	mat_t<4, 4, T> to_matrix(const frustum_t<T>& AFrustum);
	public:
		DLL_FUNCTION_STATIC(auto)	perspective(T fovy, T aspect, T zNear, T zFar) -> frustum_t;
		DLL_FUNCTION(auto)	 		to_matrix() const -> mat_t<4, 4, T>;
	};
	
	template<typename T>	mat_t<4, 4, T> orthoMatrix(T left, T right, T bottom, T top);
	template<typename T>	mat_t<4, 4, T> orthoMatrix(T left, T right, T bottom, T top, T zNear, T zFar);
	template<typename T>	mat_t<4, 4, T> frustumMatrix(T left, T right, T bottom, T top, T nearVal, T farVal);
	template<typename T>	mat_t<4, 4, T> perspectiveFovMatrix(T fov, T width, T height, T zNear, T zFar);
	template<typename T>	mat_t<4, 4, T> perspectiveMatrix(T fovy, T aspect, T zNear, T zFar);
	template<typename T>	mat_t<4, 4, T> infinitePerspectiveMatrix(T fovy, T aspect, T zNear);
	// Infinite projection matrix: http://www.terathon.com/gdc07_lengyel.pdf
	template<typename T>	mat_t<4, 4, T> tweakedInfinitePerspectiveMatrix(T fovy, T aspect, T zNear, T ep);
	
	template<typename T>	FORCE_INLINE	mat_t<4, 4, T> to_matrix(const frustum_t<T>& F) {return F.to_matrix();}
	
	
	
	template <typename T>
	struct TFrustum
	{
	private:
		enum	EFrustumSemanthic{fsUnknown,fsBoard,fsFOV,fsFocus,};
		EFrustumSemanthic	Semanthic;
	public:
		const T	near_plane, far_plane, aspect, offset;	// Always valid
		const T	top, bottom, right, left; // Board data
		const T	focus_distance, focus_radius; // Focus data
		const T	vfov, hfov; // FOV data

		SHAPES_FUNC_DEF	TFrustum();
		template<typename T2, int D, bool A>
		SHAPES_FUNC_DEF	TFrustum(const aabox_t<D, T2>&);

		DLL_FUNCTION(auto)	set_offset(const T&) -> TFrustum&;
		DLL_FUNCTION(void)	set_board(const T& iWidth, const T& iHeight, const T& iNearPlaneDistance, const T& iFarPlaneDistance);
		DLL_FUNCTION(void)	set_board(const T& AMinValueX, const T& AMaxValueX, const T& AMinValueY, const T& AMaxValueY, const T& iNearPlaneDistance, const T& iFarPlaneDistance);
		DLL_FUNCTION(void)	set_focus(const T& iDistance, const T& iRadius, const T& iAspectRatio, const T& iNearPlaneDistance, const T& iFarPlaneDistance);
		DLL_FUNCTION(void)	set_vfov(const T& iVerticalFov, const T& iAspectRatio, const T& iNearPlaneDistance, const T& iFarPlaneDistance);
		DLL_FUNCTION(void)	set_fovs(const T& iVerticalFov, const T& iHorizontalFov, const T& iNearPlaneDistance, const T& iFarPlaneDistance);
		DLL_FUNCTION(auto)	subfrustum(int _NumSections, int _CurrentSection,T *_pNearPlane,T *_pFarPlane, const T& _NearPlane, const T& _FarPlane) const -> TFrustum;

		SHAPES_FUNC_STT auto	make_board(const T& iWidth, const T& iHeight, const T& iNearPlaneDistance, const T& iFarPlaneDistance) -> TFrustum;
		SHAPES_FUNC_STT auto	make_board(const T& AMinValueX, const T& AMaxValueX, const T& AMinValueY, const T& AMaxValueY, const T& iNearPlaneDistance, const T& iFarPlaneDistance) -> TFrustum;
		SHAPES_FUNC_STT auto	make_focus(const T& iDistance, const T& iRadius, const T& iAspectRatio, const T& iNearPlaneDistance, const T& iFarPlaneDistance) -> TFrustum;
		SHAPES_FUNC_STT auto	make_vfov(const T& iVerticalFov, const T& iAspectRatio, const T& iNearPlaneDistance, const T& iFarPlaneDistance) -> TFrustum;
		SHAPES_FUNC_STT auto	make_fovs(const T& iVerticalFov, const T& iHorizontalFov, const T& iNearPlaneDistance, const T& iFarPlaneDistance) -> TFrustum;

		SHAPES_FUNC_INL	bool	is_focus() const	{return Semanthic==fsFocus;}
		SHAPES_FUNC_INL	bool	is_board() const	{return Semanthic==fsBoard;}
		SHAPES_FUNC_INL	bool	is_fov() const	{return Semanthic==fsFOV;}
	};
	
	template<int D, typename T>
	struct TSixPlanes
	{
		TPlane<D, T>	Plane[6];
		
		inline	TSixPlanes(){};
		template<int D2, typename T2>
		TSixPlanes& operator=(const aabox_t<D2, T2>&);
		template<typename T2>
		TSixPlanes& operator=(const TFrustum<T2>&);
	};

	template <typename M>
	struct TVolume
	{
		inline	TVolume(){};
		M	Matrix;
		
		TVolume& operator=(const M&);
		template<int D2, typename T2>
		TVolume& operator=(const aabox_t<D2, T2>&);
		template<typename T2>
		TVolume& operator=(const TFrustum<T2>&);
		
		typename M::Vertex		vertex(bool _Up, bool _Right, bool _Near) const;
		typename M::Vertex		local_vertex(bool AMaxValueX, bool AMaxValueY, bool AMaxValueZ) const;
	};


	template<
	int D1, typename T1, bool _ForwardsInfinite1, bool _BackwardsInfinite1,
	int D2, typename T2, bool _ForwardsInfinite2, bool _BackwardsInfinite2
	>
	bool	intersect(const TLinearShape<D1, T1, _ForwardsInfinite1, _BackwardsInfinite1>&, const TLinearShape<D2, T2, _ForwardsInfinite2, _BackwardsInfinite2>&);

	
	#define	DECLARE_SHAPE(_Shape, _Type, _NewType) \
	extern template struct _Shape<_Type>; \
	typedef _Shape<_Type> _NewType;

	#define	DECLARE_SHAPE2(_Shape, _Type, _Dimension, _NewType) \
	extern template struct _Shape<_Dimension, _Type>; \
	typedef _Shape<_Dimension, _Type> _NewType;

	#define	DECLARE_SHAPE2F(_Shape, _NewType) \
	DECLARE_SHAPE2(_Shape, float, 2, _NewType##2d_f32) \
	DECLARE_SHAPE2(_Shape, double, 2, _NewType##2d_f64) \
	DECLARE_SHAPE2(_Shape, float, 3, _NewType##3d_f32) \
	DECLARE_SHAPE2(_Shape, double, 3, _NewType##3d_f64)

	DECLARE_SHAPE(TFrustum, float, frustum3d_f32)
	DECLARE_SHAPE(TFrustum, double, frustum3d_f64)

	DECLARE_SHAPE2F(TTriangle, triangle)
	DECLARE_SHAPE2F(TSegment, segment)
	DECLARE_SHAPE2F(aabox_t, aabox)
	DECLARE_SHAPE2F(TSphere, sphere)

	typedef aabox3d_f32	aabox3d;
	
	DECLARE_SHAPE2(aabox_t, int64_t, 2, aabox2d_i64)
	DECLARE_SHAPE2(aabox_t, int64_t, 3, aabox3d_i64)
	DECLARE_SHAPE2(aabox_t, int32_t, 3, aabox3d_i32)
	DECLARE_SHAPE2(aabox_t, int32_t, 2, aabox2d_i32)
	DECLARE_SHAPE2(aabox_t, int16_t, 3, aabox3d_i16)
	DECLARE_SHAPE2(aabox_t, int16_t, 2, aabox2d_i16)
}




