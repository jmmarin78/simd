// This software is MIT licensed (see LICENSE)

#define BUILD_DLL

#include <math.h>

#include <pre/number.h>
#include <simd/shapes.h>

#include "defines.h"


namespace simd
{
	template<int D, typename T>
	template<int D2, typename T2>
	aabox_t<D, T>::aabox_t(const aabox_t<D2, T2>& i)
	{
		for(int a = 0; a < D && a < D2; a++)
		{
			min[a] = i.min[a];
			max[a] = i.max[a];
		}
		
	}

	//	template <int D, typename T>
	//	V	GetRelativeCoordinates(const aabox_t<D, T>& _Box,const V& _Coordinates)
	//	{
	//		return (_Coordinates-_Box.min) / _Box.GetSize();
	//	}
	//
	//	template <int D, typename T>
	//	aabox_t<D, T>	ScaleFromCenter(const aabox_t<D, T>& _Box,const V& _Factor)
	//	{
	//		V	mid=_Box.GetCenter();
	//		V	vec=vector(mid,_Box.max) * _Factor;
	//		return aabox_t<D, T>(mid-vec,mid+vec);
	//	}
	
	template <int D, typename T>
	aabox_t<D, T>::aabox_t(const T& iMinX,const T& iMinY,const T& iMaxX,const T& iMaxY)
	{
		min.x = iMinX;
		min.y = iMinY;
		max.x = iMaxX;
		max.y = iMaxY;
	}
	
	
	template <int D, typename T>
	aabox_t<D, T>::aabox_t(const pack_t<D, T>& _Min,const pack_t<D, T>& _Max)
	{
		min = _Min;
		max = _Max;
	}
	
	template <int D, typename T>
	aabox_t<D, T>::aabox_t(const pack_t<D, T>& _Center,const T& _Radius)
	{
		min = _Center-pack_t<D, T>(_Radius);
		max = _Center+pack_t<D, T>(_Radius);
	}
	
	template <int D, typename T>
	aabox_t<D, T>::aabox_t(const T& _Min,const T& _Max)
	{
		min = _Min;
		max = _Max;
	}
	
	template <int D, typename T>
	aabox_t<D, T>::aabox_t(const pack_t<D, T>* _Points, int _PointsCount)
	{
		set(_Points, _PointsCount);
	}
	
	template <int D, typename T>
	aabox_t<D, T>::aabox_t(const aabox_t<D, T> _Boxes[], int _BoxCount)
	{
		const aabox_t<D, T>*	boxes=_Boxes;
		if(_BoxCount > 0)
		{
			min=boxes->min;
			max=boxes->max;
			_BoxCount--;
			
			while(_BoxCount--)
			{
				boxes++;
				min.x = pre::min(min.x, boxes->min.x);
				min.y = pre::min(min.y, boxes->min.y);
				
				max.x = pre::max(max.x, boxes->max.x);
				max.y = pre::max(max.y, boxes->max.y);
				
				if(D == 3)
				{
					min[2] = pre::min(min[2], boxes->min[2]);
					max[2] = pre::max(max[2], boxes->max[2]);
				}
			}
		}
		else
		{
			for(int a = 0; a < D; a++)
			{
				min[a] = std::numeric_limits<T>::max();
				min[a] = -std::numeric_limits<T>::max();
			}
		}
	}
	
	template <int D, typename T>
	pack_t<D, T>
	aabox_t<D, T>::vertex(bool _MaxX,bool _MaxY,bool _MaxZ)const
	{
		pack_t<D, T> ret;
		ret[0] = _MaxX ? max.x:min.x;
		ret[1] = _MaxY ? max.y:min.y;
		if(D == 3)
			ret[2] = _MaxZ ? max[2] : min[2];
		return ret;
	}
	
	template <int D, typename T>
	pack_t<D, T>
	aabox_t<D, T>::vertex(int i)const
	{
		pack_t<D, T> ret;
		for(int a = 0; a < D; a++)
			ret[a] = ((i >> a) & 1) ? max[a] : min[a];
		return ret;
	}
	
	template <int D, typename T>
	void
	aabox_t<D, T>::vertices(pack_t<D, T>* o)const
	{
		if(D == 2)
		{
			for(int a = 0; a < 4; a++)
				o[a] = vertex(a);
		}
		if(D == 3)
		{
			for(int a = 0; a < 8; a++)
				o[a] = vertex(a);
		}
	}
	
	template <int D, typename T>
	void
	aabox_t<D, T>::set(const pack_t<D, T>* _Points, int _PointsCount)
	{
		if(_PointsCount)
		{
			min = _Points[0];
			max = _Points[0];
			for(int a = 1; a < _PointsCount; a++)
				expand(_Points[a]);
		}
		else
		{
			min = pack_t<D, T>(T(0));
			max = pack_t<D, T>(T(0));
		}
	}
	
	
	template <int D, typename T>
	void
	aabox_t<D, T>::invalidate()
	{
		min = std::numeric_limits<T>::max();
		max = -std::numeric_limits<T>::max();
		// 	min = gMax_f64;
		// 	max = -gMax_f64;
	}
	
	template <int D, typename T>
	bool
	aabox_t<D, T>::is_valid() const
	{
		if (min.x == std::numeric_limits<T>::max())
			return false;
		return true;
	}
	
	//	template <int D, typename T>
	//	template <typename M>
	//	void
	//	aabox_t<D, T>::operator=(const TVolume<M>& _volume)
	//	{
	//
	//		set(	_volume.GetVertex(true,true,true));
	//		expand(	_volume.GetVertex(true,true,false));
	//		expand(	_volume.GetVertex(true,false,true));
	//		expand(	_volume.GetVertex(true,false,false));
	//		expand(	_volume.GetVertex(false,true,true));
	//		expand(	_volume.GetVertex(false,true,false));
	//		expand(	_volume.GetVertex(false,false,true));
	//		expand(	_volume.GetVertex(false,false,false));
	//	}
	//
	//	template <int D, typename T>
	//	template<int _VertexCount>
	//	void
	//	aabox_t<D, T>::operator=(const TBox<V,_VertexCount>& _Box)
	//	{
	//		Set(	_Box.Vertex[0]);
	//		for(int a=1;a<8;a++)
	//			Expand(	_Box.Vertex[a]);
	//	}
	
	
	//	template <int D, typename T>
	//	bool
	//	aabox_t<D, T>::intersect(T* hitDist, const pack_t<D, T>* origPt, const pack_t<D, T>* dir)
	//	{
	//		if(D==3)
	//		{
	//			TPlane<D, T> sides[6] = {	TPlane<D, T>( 1, 0, 0,-min[0]), TPlane<D, T>(-1, 0, 0, max[0]),
	//				TPlane<D, T>( 0, 1, 0,-min[1]), TPlane<D, T>( 0,-1, 0, max[1]),
	//				TPlane<D, T>( 0, 0, 1,-min[2]), TPlane<D, T>( 0, 0,-1, max[2])
	//			};
	//
	//			*hitDist = T(0);  // safe initial value
	//			pack_t<D, T> hitPt = *origPt;
	//
	//			bool inside = false;
	//
	//			for ( int i=0; (i<6) && !inside; i++ )
	//			{
	//				T cosTheta = dot( sides[i].Normal, *dir );
	//				T dist = dot ( sides[i].Normal, *origPt )+sides[i].DistanceToOrigin;
	//
	//				//  if we're nearly intersecting, just punt and call it an intersection
	//				if ( IsAlmostZero(dist) ) return true;
	//				//  skip nearly (&actually) parallel rays
	//				if ( IsAlmostZero(cosTheta) ) continue;
	//				//  only interested in intersections along the ray, not before it.
	//				*hitDist = -dist / cosTheta;
	//				if ( *hitDist < 0.f ) continue;
	//
	//				hitPt = (*hitDist)*(*dir) + (*origPt);
	//
	//				inside = true;
	//
	//				for ( int j=0; (j<6) && inside; j++ )
	//				{
	//					if ( j==i )
	//						continue;
	//					T d = dot( sides[j].Normal, hitPt )+sides[j].DistanceToOrigin;
	//
	//					inside = ((d + 0.00015) >= 0.f);
	//				}
	//			}
	//
	//			return inside;
	//		}
	//		else
	//		{
	//			return false;
	//		}
	//	}
	
	template<int D, typename T1, typename T2>
	bool	overlaps(const aabox_t<D, T1>& ABox, const pack_t<D, T2>& APoint)
	{
		for (int i = 0; i < D; i++)
		{
			if (APoint[i] < ABox.min[i])
				return false;
			if (APoint[i] > ABox.max[i])
				return false;
		}
		return true;
	}
	
	template <typename T>
	static bool RaySlabIntersect(T slabmin, T slabmax, T raystart, T rayend, T& tbenter, T& tbexit)
	{
		T raydir = rayend - raystart;
		
		// ray parallel to the slab
		if (fabs(raydir) < 1.0E-9f)
		{
			// ray parallel to the slab, but ray not inside the slab planes
			if(raystart < slabmin || raystart > slabmax)
			{
				return false;
			}
			// ray parallel to the slab, but ray inside the slab planes
			else
			{
				return true;
			}
		}
		
		// slab's enter and exit parameters
		T tsenter = (slabmin - raystart) / raydir;
		T tsexit = (slabmax - raystart) / raydir;
		
		// order the enter / exit values.
		if(tsenter > tsexit)
		{
			auto v = tsenter;
			tsenter = tsexit;
			tsexit = v;
			//swapf(tsenter, tsexit);
		}
		
		// make sure the slab interval and the current box intersection interval overlap
		if (tbenter > tsexit || tsenter > tbexit)
		{
			// nope. Ray missed the box.
			return false;
		}
		// yep, the slab and current intersection interval overlap
		else
		{
			// update the intersection interval
			tbenter = max(tbenter, tsenter);
			tbexit = min(tbexit, tsexit);
			return true;
		}
	}
	
	template <typename T>
	static bool SegmentIntersectRectangle(T a_rectangleMinX,
								   T a_rectangleMinY,
								   T a_rectangleMaxX,
								   T a_rectangleMaxY,
								   T a_p1x,
								   T a_p1y,
								   T a_p2x,
								   T a_p2y)
	{
		// Find min and max X for the segment
		
		T minX = a_p1x;
		T maxX = a_p2x;
		
		if(a_p1x > a_p2x)
		{
			minX = a_p2x;
			maxX = a_p1x;
		}
		
		// Find the intersection of the segment's and rectangle's x-projections
		
		if(maxX > a_rectangleMaxX)
		{
			maxX = a_rectangleMaxX;
		}
		
		if(minX < a_rectangleMinX)
		{
			minX = a_rectangleMinX;
		}
		
		if(minX > maxX) // If their projections do not intersect return false
		{
			return false;
		}
		
		// Find corresponding min and max Y for min and max X we found before
		
		T minY = a_p1y;
		T maxY = a_p2y;
		
		T dx = a_p2x - a_p1x;
		
		if(fabs(dx) > 0.0000001)
		{
			T a = (a_p2y - a_p1y) / dx;
			T b = a_p1y - a * a_p1x;
			minY = a * minX + b;
			maxY = a * maxX + b;
		}
		
		if(minY > maxY)
		{
			T tmp = maxY;
			maxY = minY;
			minY = tmp;
		}
		
		// Find the intersection of the segment's and rectangle's y-projections
		
		if(maxY > a_rectangleMaxY)
		{
			maxY = a_rectangleMaxY;
		}
		
		if(minY < a_rectangleMinY)
		{
			minY = a_rectangleMinY;
		}
		
		if(minY > maxY) // If Y-projections do not intersect return false
		{
			return false;
		}
		
		return true;
	}
	
	template<int D, typename T1, bool _ForwardsInfinite, bool _BackwardsInfinite, typename T2>
	bool	intersect(const TLinearShape<D, T1, _ForwardsInfinite, _BackwardsInfinite>& line, const aabox_t<D, T2>& rect)
	{
//		if (D == 2)
//		{
//			return SegmentIntersectRectangle(rect.min.x, rect.min.y, rect.max.x, rect.max.y,
//											 line.Vertex[0].x, line.Vertex[0].y, line.Vertex[1].x, line.Vertex[1].y);
//		}
//		return false;
		
		T1	tenter = 0, texit = 1;

		// test X slab

		for (int i = 0; i < D; i++)
		{
			if (!RaySlabIntersect(rect.min[i], rect.max[i], line.Vertex[0][i], line.Vertex[1][i], tenter, texit))
			{
				return false;
			}
		}
		return true;

		
//		if (D == 2)
//		{
//			// initialise to the segment's boundaries.
//
//			// test Y slab
//
//			if (!RaySlabIntersect(rect.min.y, rect.max.y, Seg.Start.y, Seg.End.y, tenter, texit))
//			{
//				return false;
//			}
//
//			// test Z slab
//			if (!RaySlabIntersect(rect.min.z, rect.max.z, Seg.Start.z, Seg.End.z, tenter, texit))
//			{
//				return false;
//			}
//
//			// all intersections in the green. Return the first time of intersection, tenter.
//			tinter = tenter;
//			return  true;
//
//			/*
//			if (overlaps(line.Vertex[0], rect) || overlaps(line.Vertex[1], rect)) {
//				return true;
//			}
//
//			auto dir = line.Vertex[0] - line.Vertex[1];
//			auto len2 = dot(dir, dir);
//			auto norm = dir * rsq(len2);
//			for (int i = 0; i < D; i++)
//				norm[i] = (norm[i] != 0) ? (1 / norm[i]) : 0;
//
//			auto min = (rect.min - line.Vertex[0]) * norm;
//			auto max = (rect.max - line.Vertex[0]) * norm;
//
//			auto tmin = fmaxf(fminf(min.x, max.x), fminf(min.y, max.y));
//			auto tmax = fminf(fmaxf(min.x, max.x), fmaxf(min.y, max.y));
//			if (tmax < 0 || tmin > tmax) {
//				return false;
//			}
//			float t = (tmin < 0.0f) ? tmax : tmin;
//			return t > 0.0f && t * t < len2;
//			 */
//		}
//		return false;
		
		
//		high_precission_t<T1, T2>::type tmin, tmax, tymin, tymax, tzmin, tzmax;
//
//		tmin = (parameters[r.sign[0]].x() - r.origin.x()) * r.inv_direction.x();
//		tmax = (parameters[1-r.sign[0]].x() - r.origin.x()) * r.inv_direction.x();
//		tymin = (parameters[r.sign[1]].y() - r.origin.y()) * r.inv_direction.y();
//		tymax = (parameters[1-r.sign[1]].y() - r.origin.y()) * r.inv_direction.y();
//		if ( (tmin > tymax) || (tymin > tmax) )
//			return false;
//		if (tymin > tmin)
//			tmin = tymin;
//		if (tymax < tmax)
//			tmax = tymax;
//		tzmin = (parameters[r.sign[2]].z() - r.origin.z()) * r.inv_direction.z();
//		tzmax = (parameters[1-r.sign[2]].z() - r.origin.z()) * r.inv_direction.z();
//		if ( (tmin > tzmax) || (tzmin > tmax) )
//			return false;
//		if (tzmin > tmin)
//			tmin = tzmin;
//		if (tzmax < tmax)
//			tmax = tzmax;
//		return ( (tmin < t1) && (tmax > t0) );

	}
	
	template	bool			intersect<2, float, false, false, float>(const TLinearShape<2, float, false, false>& line, const aabox_t<2, float>& rect);

	
//	template<
//	int D1, typename T1,
//	int D2, typename T2, bool _ForwardsInfinite, bool _BackwardsInfinite
//	>
//	inline	bool	intersect(const aabox_t<D2, T2>& b, const TLinearShape<D1, T1, _ForwardsInfinite, _BackwardsInfinite>& l) {return intersect(l, b);}
//
//	bool Box::intersect(const Ray &r, float t0, float t1) const {
//	}
//
//
//
	
	/*
	 bool
	 aabox_t<D, T>::Contains(const V &_Vertex,EConceptLocation _Mask,EConceptLocation *_pOut) const
	 {
	 if(_pOut==NULL)
	 {
	 if(_Mask & clOutside)
	 {
	 if (	min.x > _Vertex.x ||
	 min.y > _Vertex.y ||
	 min.z > _Vertex.z ||
	 max.x < _Vertex.x ||
	 max.y < _Vertex.y ||
	 max.z < _Vertex.z )
	 return true;
	 }
	 else
	 {
	 if (	min.x > _Vertex.x ||
	 min.y > _Vertex.y ||
	 min.z > _Vertex.z ||
	 max.x < _Vertex.x ||
	 max.y < _Vertex.y ||
	 max.z < _Vertex.z )
	 return false;
	 if(_Mask & clTouch)
	 {
	 if (	min.x == _Vertex.x ||
	 min.y == _Vertex.y ||
	 min.z == _Vertex.z ||
	 max.x == _Vertex.x ||
	 max.y == _Vertex.y ||
	 max.z == _Vertex.z )
	 return true;
	 }
	 if
	 }
	 }
	 return false;
	 }
	 
	 int
	 aabox_t<D, T>::Contains(const aabox_t<D, T>& _box) const
	 {
	 if ( min.x > _box.max.x ) return 0;
	 if ( min.y > _box.max.y ) return 0;
	 if ( min.z > _box.max.z ) return 0;
	 if ( max.x < _box.min.x ) return 0;
	 if ( max.y < _box.min.y ) return 0;
	 if ( max.z < _box.min.z ) return 0;
	 
	 return 3;
	 }
	 */
	template <int D, typename T>
	pack_t<D, T>
	aabox_t<D, T>::size() const
	{
		return max - min;
	}
	
	template <int D, typename T>
	pack_t<D, T>
	aabox_t<D, T>::center() const
	{
		return (max+min) * pack_t<D, T>(T(0.5));
	}
	
	//	template <int D, typename T>
	//	void
	//	aabox_t<D, T>::GetCenter(pack_t<D, T>* _pData) const
	//	{
	//		*_pData = (max+min) * V::make(T(0.5));
	//	}
	//
	//	template<typename T>
	//	static	inline	void	InternalGetSize(const T* iMin,const T* iMax,T* o, int iCount)
	//	{
	//		rep(a,iCount)
	//		o[a]=iMax[a]-iMin[a];
	//	}
	//
	//#define MACRO_InternalGetSizeInteger(_Type)
	//static	inline	void	InternalGetSize(const _Type* iMin,const _Type* iMax,_Type* o, int iCount)
	//{
	//rep(a,iCount)
	//o[a]=iMax[a]-iMin[a]+1;
	//}
	//
	//	MACRO_InternalGetSizeInteger(int8_t)
	//	MACRO_InternalGetSizeInteger(uint8_t)
	//	MACRO_InternalGetSizeInteger(int16_t)
	//	MACRO_InternalGetSizeInteger(uint16_t)
	//	MACRO_InternalGetSizeInteger(int32_t)
	//	MACRO_InternalGetSizeInteger(uint32_t)
	//	MACRO_InternalGetSizeInteger(int64_t)
	//	MACRO_InternalGetSizeInteger(uint64_t)
	
	
	template <int D, typename T>
	T
	aabox_t<D, T>::width() const
	{
		return max[0] - min[0];
	}
	
	template <int D, typename T>
	T
	aabox_t<D, T>::height() const
	{
		if (D < 2)
			return T(0);
		return max[1] - min[1];
	}
	
	
	template <int D, typename T>
	T
	aabox_t<D, T>::depth() const
	{
		if (D < 3)
			return T(0);
		
		return max[2] - min[2];
	}
	
	
	template <int D, typename T>
	T
	aabox_t<D, T>::outer_radius() const
	{
		return distance(min, max) * T(0.5);
	}
	
	template <int D, typename T>
	T sqdistance(const pack_t<D, T>& p, const aabox_t<D, T>& b)
	{
		T sqDist = T(0);
		
		for (int i = 0; i < D; i++)
		{
			// For each axis count any excess distance outside box extents
			T v = p[i];
			if (v < b.min[i])
				sqDist += (b.min[i] - v) * (b.min[i] - v);
			if (v > b.max[i])
				sqDist += (v - b.max[i]) * (v - b.max[i]);
		}
		return sqDist;
	}
	
	template <int D, typename T>
	bool
	aabox_t<D, T>::contains(const pack_t<D, T>& _Value)const
	{
		for (int i = 0; i < D; i++)
		{
			if (_Value[i] < min[i] || _Value[i] > max[i])
				return false;
		}
		return true;
	}
	
	
	template <int D, typename T>
	T distance(const pack_t<D, T>& p, const aabox_t<D, T>& b)
	{
		return xsqrt(sqdistance(p, b));
	}
	
	template <int D, typename T>
	int
	aabox_t<D, T>::dominant_axis() const
	{
		pack_t<D, T> dim = size();
		
		if ( dim.y > dim.x )
		{
			if (D==3 && dim.z > dim.y)
				return 2;
			return 1;
		}
		
		if (D==3 && dim.z > dim.x )
			return 2;
		return 0;
	}
	
	template <int D, typename T>
	void
	aabox_t<D, T>::set(const aabox_t<D, T>& _Box,int iSplit)
	{
		pack_t<D, T> hs = size()/T(2);
		
		if(D == 2)
		{
			if(iSplit == 0)
			{
				min = _Box.min;
				max = min + hs;
			}
			else if(iSplit == 1)
			{
				min = _Box.min;
				min.x += hs.x;
				max = min + hs;
				max.x += hs.x;
			}
			else if(iSplit==2)
			{
				min = _Box.min;
				min.y += hs.y;
				max = min + hs;
				max.y += hs.y;
			}
			else
			{
				min = _Box.min + hs;
				max=max;
			}
		}
		else
		{
			assert(0);
		}
	}
	
	template <int D, typename T>
	void
	aabox_t<D, T>::set(const pack_t<D, T>& _Value)
	{
		min = _Value;
		max = _Value;
	}
	
	template <int D, typename T>
	void
	aabox_t<D, T>::expand(const pack_t<D, T>* ipVertices,int iCount)
	{
		for(int i = 0; i < iCount; i++)
		{
			for(int a = 0; a < D; a++)
			{
				min[a] = pre::min(min[a], ipVertices[i][a]);
				max[a] = pre::max(max[a], ipVertices[i][a]);
			}
		}
	}
	
	template <int D, typename T>
	void
	aabox_t<D, T>::expand(const pack_t<D, T>& _Value)
	{
		for(int a = 0; a < D; a++)
		{
			min[a] = pre::min(min[a], _Value[a]);
			max[a] = pre::max(max[a], _Value[a]);
		}
	}
	
	template <int D, typename T>
	void
	aabox_t<D, T>::expand(const aabox_t<D, T>& _Box)
	{
		for(int a = 0; a < D; a++)
		{
			min[a] = pre::min(min[a], _Box.min[a]);
			max[a] = pre::max(max[a], _Box.max[a]);
		}
	}
	
	template <int D, typename T>
	void	Normalize(pack_t<D, T> _Source[],pack_t<D, T> _Detination[], int _NumVertices,const aabox_t<D, T>& _NormalizationBox)
	{
		pack_t<D, T>	size = (_NormalizationBox.max-_NormalizationBox.min)*T(0.5);
		pack_t<D, T>	isize = T(1)/size;
		
		for(int a = 0; a < _NumVertices; a++)
		{
			_Detination[a] = (_Source[a]-_NormalizationBox.min) * isize;
		}
	}
	
	//	template <int D, typename T>
	//	ERelativeLocation	GetRelativeLocation(const aabox_t<D, T>& _A,const pack_t<D, T>& _B,ERelativeLocation _Mask)
	//	{
	//		rep(a,D)
	//			if(_B[a]<_A.min[a])		return rlOutside;
	//
	//		rep(a,V::Dimension)
	//		if(_B[a]>_A.max[a])		return rlOutside;
	//
	//		rep(a,V::Dimension)
	//		if(_B[a]==_A.min[a])	return ERelativeLocation(rlContact & rlBInsideToA);
	//
	//		rep(a,V::Dimension)
	//		if(_B[a]==_A.max[a])	return ERelativeLocation(rlContact & rlBInsideToA);
	//
	//		return rlBInsideToA;
	//	}
	
	
	
	//  Transform an axis-aligned bounding box by the specified matrix, and compute a n ew axis-aligned bounding box
	//	template<typename T,typename V,int _Dimension,typename M>
	//	void XFormBoundingBox( aabox_t<D, T>& _Result, const aabox_t<D, T>& _Source, const M& _Matrix )
	//	{
	//		V  pts[8];
	//		for ( int i=0; i<8; i++ )
	//			pts[i] = _Source.GetVertex(i);
	//
	//		_Result.min = V(T(3.3e33), T(3.3e33), T(3.3e33));
	//		_Result.max = V(-T(3.3e33), -T(3.3e33), -T(3.3e33));
	//
	//		for (int i=0; i<8; i++)
	//		{
	//			V tmp;
	//			transform(_Matrix,tmp,pts[i]);
	//			_Result.Expand(tmp);
	//		}
	//	}
}



//namespace simd
//{
//
//
//	template<typename T>
//	TFrustum<T>	TFrustum<T>::make_board(const T& iWidth,const T& iHeight,const T& iNearPlaneDistance,const T& iFarPlaneDistance)
//	{
//		TFrustum<T> ret;
//		ret.set_board(iWidth,iHeight,iNearPlaneDistance,iFarPlaneDistance);
//		return ret;
//	}
//	template<typename T>
//	TFrustum<T>	TFrustum<T>::make_board(const T& iMinX,const T& iMaxX,const T& iMinY,const T& iMaxY,const T& iNearPlaneDistance,const T& iFarPlaneDistance)
//	{
//		TFrustum<T> ret;
//		ret.set_board(iMinX,iMaxX,iMinY,iMaxY,iNearPlaneDistance,iFarPlaneDistance);
//		return ret;
//	}
//	template<typename T>
//	TFrustum<T>	TFrustum<T>::make_focus(const T& iDistance,const T& iRadius,const T& iAspectRatio,const T& iNearPlaneDistance,const T& iFarPlaneDistance)
//	{
//		TFrustum<T> ret;
//		ret.set_focus(iDistance,iRadius,iAspectRatio,iNearPlaneDistance,iFarPlaneDistance);
//		return ret;
//	}
//	template<typename T>
//	TFrustum<T>	TFrustum<T>::make_vfov(const T& iVerticalFov,const T& iAspectRatio,const T& iNearPlaneDistance,const T& iFarPlaneDistance)
//	{
//		TFrustum<T> ret;
//		ret.set_vfov(iVerticalFov,iAspectRatio,iNearPlaneDistance,iFarPlaneDistance);
//		return ret;
//	}
//	template<typename T>
//	TFrustum<T>	TFrustum<T>::make_fovs(const T& iVerticalFov,const T& iHorizontalFov,const T& iNearPlaneDistance,const T& iFarPlaneDistance)
//	{
//		TFrustum<T> ret;
//		ret.set_fovs(iVerticalFov,iHorizontalFov,iNearPlaneDistance,iFarPlaneDistance);
//		return ret;
//	}
//
//	template<typename T>
//	TFrustum<T>::TFrustum()
//	: near_plane(T(0)),far_plane(T(0)),aspect(T(0)),offset(T(0))
//	, top(T(0)),bottom(T(0)),right(T(0)),left(T(0))
//	, focus_distance(T(0)), focus_radius(T(0))
//	, vfov(T(0)), hfov(T(0))
//
//	{
//		Semanthic=fsUnknown;
//	}
//
//	template<typename T>
//	TFrustum<T>&
//	TFrustum<T>::set_offset(const T& i)
//	{
//		num::set(i,pre::remove_const(offset));
//		return *this;
//	}
//	template<typename T>
//	void
//	TFrustum<T>::set_board(const T& iWidth,const T& iHeight,const T& iNearPlaneDistance,const T& iFarPlaneDistance)
//	{
//		pre::remove_const(right)=iWidth*0.5f;
//		pre::remove_const(left)=iWidth*(-0.5f);
//		pre::remove_const(top)=iHeight*0.5f;
//		pre::remove_const(bottom)=iHeight*(-0.5f);
//		num::set(iNearPlaneDistance,pre::remove_const(near_plane));
//		num::set(iFarPlaneDistance,pre::remove_const(far_plane));
//		pre::remove_const(aspect)=iHeight/iWidth;
//		Semanthic=fsBoard;
//	}
//	template<typename T>
//	void
//	TFrustum<T>::set_board(const T& iMinX,const T& iMaxX,const T& iMinY,const T& iMaxY,const T& iNearPlaneDistance,const T& iFarPlaneDistance)
//	{
//		num::set(iMaxX,pre::remove_const(right));
//		num::set(iMinX,pre::remove_const(left));
//		num::set(iMaxY,pre::remove_const(top));
//		num::set(iMinY,pre::remove_const(bottom));
//		num::set(iNearPlaneDistance,pre::remove_const(near_plane));
//		num::set(iFarPlaneDistance,pre::remove_const(far_plane));
//		pre::remove_const(aspect)=(iMaxX-iMinX)/(iMaxY-iMinY);
//		Semanthic=fsBoard;
//	}
//	template<typename T>
//	void
//	TFrustum<T>::set_focus(const T& iDistance,const T& iRadius,const T& iAspectRatio,const T& iNearPlaneDistance,const T& iFarPlaneDistance)
//	{
//		num::set(iDistance,pre::remove_const(focus_distance));
//		num::set(iRadius,pre::remove_const(focus_radius));
//		num::set(iAspectRatio,pre::remove_const(aspect));
//		num::set(iNearPlaneDistance,pre::remove_const(near_plane));
//		num::set(iFarPlaneDistance,pre::remove_const(far_plane));
//		Semanthic=fsFocus;
//	}
//	template<typename T>
//	void
//	TFrustum<T>::set_vfov(const T& iVerticalFov,const T& iAspectRatio,const T& iNearPlaneDistance,const T& iFarPlaneDistance)
//	{
//		num::set(iVerticalFov,pre::remove_const(vfov));
//		num::zero(pre::remove_const(hfov));
//		num::set(iAspectRatio,pre::remove_const(aspect));
//		num::set(iNearPlaneDistance,pre::remove_const(near_plane));
//		num::set(iFarPlaneDistance,pre::remove_const(far_plane));
//		Semanthic=fsFOV;
//	}
//	template<typename T>
//	void
//	TFrustum<T>::set_fovs(const T& iVerticalFov,const T& iHorizontalFov,const T& iNearPlaneDistance,const T& iFarPlaneDistance)
//	{
//		num::set(iVerticalFov,pre::remove_const(vfov));
//		num::set(iHorizontalFov,pre::remove_const(hfov));
//		pre::remove_const(aspect)=iVerticalFov/iHorizontalFov;
//		num::set(iNearPlaneDistance,pre::remove_const(near_plane));
//		num::set(iFarPlaneDistance,pre::remove_const(far_plane));
//		Semanthic=fsFOV;
//	}
//}

namespace simd
{

	
}

namespace simd
{
#define INSTANCIATE_FUNCTIONS_FRUSTUM(T) \
template					TFrustum<T>::TFrustum(); \
template	TFrustum<T>&	TFrustum<T>::set_offset(const T& i); \
template	void			TFrustum<T>::set_board(const T& iWidth,const T& iHeight,const T& iNearPlaneDistance,const T& iFarPlaneDistance); \
template	void			TFrustum<T>::set_board(const T& iMinX,const T& iMaxX,const T& iMinY,const T& iMaxY,const T& iNearPlaneDistance,const T& iFarPlaneDistance); \
template	void			TFrustum<T>::set_focus(const T& iDistance,const T& iRadius,const T& iAspectRatio,const T& iNearPlaneDistance,const T& iFarPlaneDistance); \
template	void			TFrustum<T>::set_vfov(const T& iVerticalFov,const T& iAspectRatio,const T& iNearPlaneDistance,const T& iFarPlaneDistance); \
template	void			TFrustum<T>::set_fovs(const T& iVerticalFov,const T& iHorizontalFov,const T& iNearPlaneDistance,const T& iFarPlaneDistance); \
template	TFrustum<T>		TFrustum<T>::make_board(const T& iWidth,const T& iHeight,const T& iNearPlaneDistance,const T& iFarPlaneDistance); \
template	TFrustum<T>		TFrustum<T>::make_board(const T& iMinX,const T& iMaxX,const T& iMinY,const T& iMaxY,const T& iNearPlaneDistance,const T& iFarPlaneDistance); \
template	TFrustum<T>		TFrustum<T>::make_focus(const T& iDistance,const T& iRadius,const T& iAspectRatio,const T& iNearPlaneDistance,const T& iFarPlaneDistance); \
template	TFrustum<T>		TFrustum<T>::make_vfov(const T& iVerticalFov,const T& iAspectRatio,const T& iNearPlaneDistance,const T& iFarPlaneDistance); \
template	TFrustum<T>		TFrustum<T>::make_fovs(const T& iVerticalFov,const T& iHorizontalFov,const T& iNearPlaneDistance,const T& iFarPlaneDistance);
	
	
	
	
//	INSTANCIATE_FUNCTIONS_FRUSTUM(float)
//	INSTANCIATE_FUNCTIONS_FRUSTUM(double)
	
	
	
#define INSTANCIATE_FUNCTIONS(T,D,A) \
template					aabox_t<D, T>::aabox_t(const pack_t<D, T>& _Min,const pack_t<D, T>& _Max); \
template					aabox_t<D, T>::aabox_t(const T& iMinX,const T& iMinY,const T& iMaxX,const T& iMaxY); \
template					aabox_t<D, T>::aabox_t(const pack_t<D, T>& _Center,const T& _Radius); \
template					aabox_t<D, T>::aabox_t(const T& _Min,const T& _Max); \
template					aabox_t<D, T>::aabox_t(const pack_t<D, T>* _Points, int _PointsCount); \
template					aabox_t<D, T>::aabox_t(const aabox_t<D, T> _Boxes[], int _BoxCount); \
template	void			aabox_t<D, T>::invalidate(); \
template	bool			aabox_t<D, T>::is_valid() const; \
template	pack_t<D, T>	aabox_t<D, T>::center()const; \
template	pack_t<D, T>	aabox_t<D, T>::size()const; \
template	T				aabox_t<D, T>::width()const; \
template	T				aabox_t<D, T>::height()const; \
template	T				aabox_t<D, T>::depth()const; \
template	T				aabox_t<D, T>::outer_radius()const; \
template	pack_t<D, T>	aabox_t<D, T>::vertex(bool _MaxX,bool _MaxY,bool _MaxZ)const; \
template	pack_t<D, T>	aabox_t<D, T>::vertex(int i)const; \
template	void			aabox_t<D, T>::vertices(pack_t<D, T>* o)const; \
template	void			aabox_t<D, T>::set(const aabox_t<D, T>& _Box,int iSplit); \
template	void			aabox_t<D, T>::set(const pack_t<D, T>* _Points, int _PointsCount); \
template	void			aabox_t<D, T>::expand(const pack_t<D, T>*,int); \
template	void			aabox_t<D, T>::expand(const pack_t<D, T>&); \
template	void			aabox_t<D, T>::expand(const aabox_t<D, T>& _Box); \
template	bool			aabox_t<D, T>::contains(const pack_t<D, T>& _Value)const; \

	
	
	
	INSTANCIATE_FUNCTIONS(float,2,false)
	INSTANCIATE_FUNCTIONS(double,2,false)
	INSTANCIATE_FUNCTIONS(float,3,false)
	INSTANCIATE_FUNCTIONS(double,3,false)
	//	INSTANCIATE_FUNCTIONS(float,2,true)
	//	INSTANCIATE_FUNCTIONS(double,2,true)
	//	INSTANCIATE_FUNCTIONS(float,3,true)
	//	INSTANCIATE_FUNCTIONS(double,3,true)
	
}

namespace simd
{
	template<typename TD, typename TS>
	void	shapes_operator_assign(TD* o,const TS& i)
	{
		shapes_operator_assign_procedure(o, i);
	}
}
