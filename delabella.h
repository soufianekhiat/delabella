/*
* Simple DelaBella: sdb
* Aim to be Closer to C-like less C++
Froked from:
DELABELLA - Delaunay triangulation library
Copyright (C) 2018-2022 GUMIX - Marcin Sokalski
*/

#ifndef DELABELLA_H
#define DELABELLA_H

#ifdef __cplusplus
namespace sdb {
#endif

#ifndef SCALAR_TYPE
#define SCALAR_TYPE double
#endif
#ifndef INTEGER_TYPE
#define INTEGER_TYPE int
#endif
typedef SCALAR_TYPE		Scalar;
typedef INTEGER_TYPE	Integer;

struct SDB
{
	struct Vertex;
	struct Simplex;
	struct Iterator;

	struct Vertex
	{
		Vertex* next; // next in internal / boundary set of vertices
		Simplex* sew; // one of triangles sharing this vertex
		Scalar x, y; // coordinates (input copy)
		Integer i; // index of original point

		inline const Simplex* StartIterator(Iterator* it/*not_null*/) const;
	};

	struct Simplex
	{
		Vertex* v[3];  // 3 vertices spanning this triangle
		Simplex* f[3]; // 3 adjacent faces, f[i] is at the edge opposite to vertex v[i]
		Simplex* next; // next triangle (of delaunay set or hull set)
		Integer index;       // list index
		
		unsigned char flags; 

		bool IsDelaunay() const
		{
			return !(flags & 0b10000000);
		}

		bool IsInterior(int at) const
		{
			return flags & 0b01000000;
		}

		bool IsEdgeFixed(int at) const
		{
			return flags & (1 << at);
		}

		inline const Simplex* StartIterator(Iterator* it/*not_null*/, int around/*0,1,2*/) const;
	};

	struct Iterator
	{
		const Simplex* current;
		int around;

		const Simplex* Next()
		{
			int pivot = around + 1;
			if (pivot == 3)
				pivot = 0;

			Simplex* next = current->f[pivot];
			Vertex* v = current->v[around];

			if (next->v[0] == v)
				around = 0;
			else
			if (next->v[1] == v)
				around = 1;
			else
				around = 2;

			current = next;
			return current;
		}

		const Simplex* Prev()
		{
			int pivot = around - 1;
			if (pivot == -1)
				pivot = 2;

			Simplex* prev = current->f[pivot];
			Vertex* v = current->v[around];

			if (prev->v[0] == v)
				around = 0;
			else
			if (prev->v[1] == v)
				around = 1;
			else
				around = 2;

			current = prev;
			return current;
		}
	};

	static SDB* Create();
	virtual void Destroy() = 0;

	virtual void SetErrLog(int(*proc)(void* stream, const char* fmt, ...), void* stream) = 0;

	// return 0: no output 
	// negative: all points are colinear, output hull vertices form colinear segment list, no triangles on output
	// positive: output hull vertices form counter-clockwise ordered segment contour, delaunay and hull triangles are available
	// if 'y' pointer is null, y coords are treated to be located immediately after every x
	// if advance_bytes is less than 2*sizeof coordinate type, it is treated as 2*sizeof coordinate type  
	virtual Integer Triangulate(Integer points, const Scalar* x, const Scalar* y = 0, size_t advance_bytes = 0, Integer stop = -1) = 0;

	// num of points passed to last call to Triangulate()
	virtual Integer GetNumInputPoints() const = 0;

	// num of indices returned from last call to Triangulate()
	virtual Integer GetNumOutputIndices() const = 0;

	// num of hull faces (non delaunay triangles)
	virtual Integer GetNumOutputHullFaces() const = 0;

	// num of boundary vertices
	virtual Integer GetNumBoundaryVerts() const = 0;

	// num of internal vertices
	virtual Integer GetNumInternalVerts() const = 0;

	// when called right after Triangulate() / Constrain() / FloodFill() it returns number of triangles,
	// but if called after Polygonize() it returns number of polygons
	virtual Integer GetNumPolygons() const = 0;

	virtual const Simplex* GetFirstDelaunaySimplex() const = 0; // valid only if Triangulate() > 0
	virtual const Simplex* GetFirstHullSimplex() const = 0; // valid only if Triangulate() > 0
	virtual const Vertex*  GetFirstBoundaryVertex() const = 0; // if Triangulate() < 0 it is list, otherwise closed contour! 
	virtual const Vertex*  GetFirstInternalVertex() const = 0; // valid only if Triangulate() > 0

	// given input point index, returns corresponding vertex pointer 
	virtual const Vertex*  GetVertexByIndex(Integer i) const = 0;

	// insert constraint edges into triangulation, valid only if Triangulate() > 0
	virtual Integer ConstrainEdges(Integer edges, const Integer* pa, const Integer* pb, size_t advance_bytes) = 0;

	// assigns interior / exterior flags to all faces, valid only if Triangulate() > 0
	// returns number of 'land' faces (they start at GetFirstDelaunaySimplex)
	// optionally <exterior> pointer is set to the first 'sea' face
	// if invert is set, outer-most faces will become 'land' (instead of 'sea')
	// depth controls how many times (0=INF) wave front can pass through constraints rings
	virtual Integer FloodFill(bool invert, const Simplex** exterior = 0, int depth = 0) = 0;

	// groups adjacent faces, not separated by constraint edges, built on concyclic vertices into polygons
	// first 3 vertices of a polygon are all 3 vertices of first face Simplex::v[0], v[1], v[2]
	// every next face in polygon defines additional 1 polygon vertex at its Simplex::v[0]
	// usefull as preprocessing step before genereating voronoi diagrams 
	// and as unification step before comparing 2 or more triangulations
	// valid only if Triangulate() > 0
	virtual Integer Polygonize(const Simplex* poly[/*GetNumOutputIndices()/3*/] = 0) = 0; 

	// GenVoronoiDiagramVerts(), valid only if Triangulate() > 0
	// it makes sense to call it prior to constraining only
	// generates VD vertices (for use with VD edges or VD polys indices)
	// <x>,<y> can be null if only number of vertices is needed
	// assuming:
	//   N = GetNumBoundaryVerts()
	//   M = GetNumInternalVerts()
	//   P = GetNumPolygons()
	//   V = P + N
	// <x>,<y> array must be (at least) V elements long
	// first P <x>,<y> elements will contain internal points (divisor W=1)
	// next N <x>,<y> elements will contain edge normals (divisor W=0)
	// function returns number vertices filled (V) on success, otherwise 0
	virtual Integer GenVoronoiDiagramVerts(Scalar* x, Scalar* y, size_t advance_bytes = 0) const = 0;


	// GenVoronoiDiagramEdges(), valid only if Triangulate() > 0
	// it makes sense to call it prior to constraining only
	// generates unidirected VD edges (without ones in opposite direction)
	// assuming:
	//   N = GetNumBoundaryVerts()
	//   M = GetNumInternalVerts()
	//   P = GetNumPolygons()
	//   Integer = 2 * (N + M + P - 1)
	// <indices> must be (at least) Integer elements long or must be null
	// every pair of consecutive values in <indices> represent VD edge
	// there is no guaranteed correspondence between edges order and other data
	// function returns number of indices filled (I) on success, otherwise 0
	virtual Integer GenVoronoiDiagramEdges(Integer* indices, size_t advance_bytes = 0) const = 0;

	// GenVoronoiDiagramPolys() valid only if Triangulate() > 0
	// it makes sense to call it prior to constraining only
	// generates VD polygons
	// assuming:
	//   N = GetNumBoundaryVerts()
	//   M = GetNumInternalVerts()
	//   P = GetNumPolygons()
	//   Integer = 3 * (N + M) + 2 * (P - 1) + N
	// <indices> must be (at least) Integer elements long or must be null
	// first M polys in <indices> represent closed VD cells, thay are in order
	// and corresponding to vertices from GetFirstInternalVertex() -> Vertex::next list
	// next N polys in <indices> represent open VD cells, thay are in order
	// and corresponding to vertices from GetFirstBoundaryVertex() -> Vertex::next list
	// every poly written to <indices> is terminated with ~0 value
	// if both <indices> and <closed_indices> are not null, 
	// number of closed VD cells indices is written to <closed_indices>
	// function returns number of indices filled (I) on success, otherwise 0
	virtual Integer GenVoronoiDiagramPolys(Integer* indices, size_t advance_bytes=0, Integer* closed_indices=0) const = 0;

	virtual void CheckTopology() const = 0;
};

inline const typename SDB::Simplex* SDB::Simplex::StartIterator( SDB::Iterator* it/*not_null*/, int around/*0,1,2*/) const
{
	it->current = this;
	it->around = around;
	return this;
}

inline const typename SDB::Simplex* SDB::Vertex::StartIterator( SDB::Iterator* it/*not_null*/) const
{
	it->current = sew;
	if (sew->v[0] == this)
		it->around = 0;
	else
	if (sew->v[1] == this)
		it->around = 1;
	else
		it->around = 2;
	return sew;
}

#ifdef __cplusplus
}
#endif

#endif
