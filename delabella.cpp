/*
* Simple DelaBella: sdb
* Aim to be Closer to C-like less C++
Froked from:
DELABELLA - Delaunay triangulation library
Copyright (C) 2018-2022 GUMIX - Marcin Sokalski
*/

//#define DELABELLA_AUTOTEST

// in case of troubles, allows to see if any assert pops up.
// define it globally (like with -DDELABELLA_AUTOTEST)

// this is to work around bug in the DekkersProduct()
// note Splitter constant is also wrong but not used outside DekkersProduct()
#undef FP_FAST_FMAF
#define FP_FAST_FMAF
#undef FP_FAST_FMA
#define FP_FAST_FMA
#undef FP_FAST_FMAL
#define FP_FAST_FMAL

#include <random>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <stdint.h>
#include <string.h>

#include "delabella.h"
#include "predicates.h"

#ifdef __cplusplus
namespace sdb{
#endif

// benching hack, fixme!
uint64_t sorting_bench = 0;

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif

static uint64_t uSec()
{
#ifdef _WIN32
	LARGE_INTEGER c;
	static LARGE_INTEGER f;
	static BOOL bf = QueryPerformanceFrequency(&f);
	QueryPerformanceCounter(&c);
	uint64_t n = c.QuadPart;
	uint64_t d = f.QuadPart;
	uint64_t m = 1000000;
	// calc microseconds = n*m/d carefully!
	// naive mul/div would work only for upto 5h on 1GHz freq
	// we exploit fact that m*d fits in uint64 (upto 18THz freq)
	// so n%d*m fits as well,
	return n / d * m + n % d * m / d;
#else
	timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return (uint64_t)ts.tv_sec * 1000000 + ts.tv_nsec / 1000;
#endif
}

struct SDBImpl
{
	SDBImpl() : vert_map(0),
				vert_alloc(0),
				face_alloc(0),
				max_verts(0),
				max_faces(0),
				first_dela_face(0),
				first_hull_face(0),
				first_boundary_vert(0),
				first_internal_vert(0),
				inp_verts(0),
				out_verts(0),
				out_boundary_verts(0),
				polygons(0),
				out_hull_faces(0),
				unique_points(0),
				errlog_proc(0),
				errlog_file(0)
	{}

	struct Face;

	struct Iter : SDBIterator
	{};

	struct Vert : SDBVertex
	{
		static const Scalar resulterrbound;

		static bool overlap(const Vert *v1, const Vert *v2)
		{
			return v1->x == v2->x && v1->y == v2->y;
		}

		/*
		bool operator<(const Vert &v) const
		{
			// assuming no dups
			// return this->x < v.x;

			// seams to work, but results in super slow triangulation on gamma+
			// return this->x < v.x || this->x == v.x && this->y < v.y;

			// reversing?
			// return this->x > v.x || this->x == v.x && this->y > v.y;

			// reversing paraboloid, somewhat faster without predicate
			{
				Scalar ax = this->x+Scalar(1.1);
				Scalar ay = this->y+Scalar(0.9);
				Scalar bx = v.x+Scalar(1.1);
				Scalar by = v.y+Scalar(0.9);
				Scalar a = ax*ax+ay*ay;
				Scalar b = bx*bx+by*by;
				if (a == b)
				{
					if (this->x > v.x || this->x == v.x && this->y > v.y)
						return true;
					return false;
				}
				return a > b;
			}

			// reversing paraboloid, speeds up: gam,sym,hex, bias (1.1,0.9) added mostly for sym to de-sym it
			{
				Scalar dif = predicates::adaptive::sqrlendif2d(this->x + Scalar(1.1), this->y + Scalar(0.9), v.x + Scalar(1.1), v.y + Scalar(0.9));
				if (dif > 0)
					return true;
				if (dif < 0)
					return false;
				if (this->x > v.x || this->x == v.x && this->y > v.y)
					return true;
				return false;
			}


			{	// somewhat faster, poor compiler inlining?
				const Scalar a = this->x * this->x + this->y * this->y;
				const Scalar b = v.x * v.x + v.y * v.y;
				const Scalar c = a - b;
				if (std::abs(c) > (a + b) * resulterrbound)
					return c < 0;
			}

			Scalar dif = predicates::adaptive::sqrlendif2d(this->x, this->y, v.x, v.y);

			if (dif < 0)
				return true;
			if (dif > 0)
				return false;
			if (this->x < v.x || this->x == v.x && this->y < v.y)
				return true;
			return false;
		}
		*/
	};


	struct Face : SDBSimplex
	{
		static const Scalar iccerrboundA;

		void RotateEdgeFlagsCCW() // <<
		{
			uint8_t f = this->flags;
			this->flags = ((f << 1) & 0b00110110) | ((f >> 2) & 0b00001001) | (f & 0b11000000);
		}

		void RotateEdgeFlagsCW() // >>
		{
			uint8_t f = this->flags;
			this->flags = ((f >> 1) & 0b00011011) | ((f << 2) & 0b00100100) | (f & 0b11000000);
		}

		void ToggleEdgeFixed(int at)
		{
			this->flags = (this->flags ^ (0b00001000 << at)) | (0b00000001 << at);
		}

		uint8_t GetEdgeBits(int at) const
		{
			return (this->flags >> at) & 0b00001001;
		}

		void SetEdgeBits(int at, uint8_t bits)
		{
			this->flags = (bits << at) | (this->flags & ~(0b00001001 << at));
		}

		static Face *Alloc(Face **from)
		{
			Face *f = *from;
			*from = (Face *)f->next;
			f->next = 0;
			return f;
		}

		void Free(Face **to)
		{
			this->next = *to;
			*to = this;
		}

		Face *Next(const Vert *p) const
		{
			if (this->v[0] == p)
				return (Face *)this->f[1];
			if (this->v[1] == p)
				return (Face *)this->f[2];
			if (this->v[2] == p)
				return (Face *)this->f[0];
			return 0;
		}

		bool signN() const
		{
			return 0 > predicates::adaptive::orient2d(
						   this->v[0]->x, this->v[0]->y,
						   this->v[1]->x, this->v[1]->y,
						   this->v[2]->x, this->v[2]->y);
		}

		bool sign0() const
		{
			return 0 == predicates::adaptive::orient2d(
							this->v[0]->x, this->v[0]->y,
							this->v[1]->x, this->v[1]->y,
							this->v[2]->x, this->v[2]->y);
		}

		bool dot0(const Vert &p) const
		{
			return predicates::adaptive::incircle(
					   p.x, p.y,
					   this->v[0]->x, this->v[0]->y,
					   this->v[1]->x, this->v[1]->y,
					   this->v[2]->x, this->v[2]->y) == 0;
		}

		bool dotP(const Vert &p) const
		{
			return predicates::adaptive::incircle(
					   p.x, p.y,
					   this->v[0]->x, this->v[0]->y,
					   this->v[1]->x, this->v[1]->y,
					   this->v[2]->x, this->v[2]->y) > 0;
		}

		bool dotN(const Vert &p) const
		{
			return predicates::adaptive::incircle(
					   p.x, p.y,
					   this->v[0]->x, this->v[0]->y,
					   this->v[1]->x, this->v[1]->y,
					   this->v[2]->x, this->v[2]->y) < 0;
		}

		bool dotNP(const Vert &p) /*const*/
		{
			{ // somewhat faster, poor compiler inlining?

				const Scalar dx = this->v[2]->x;
				const Scalar dy = this->v[2]->y;

				const Scalar adx = p.x - dx;
				const Scalar ady = p.y - dy;
				const Scalar bdx = this->v[0]->x - dx;
				const Scalar bdy = this->v[0]->y - dy;
				const Scalar cdx = this->v[1]->x - dx;
				const Scalar cdy = this->v[1]->y - dy;

				const Scalar adxcdy = adx * cdy;
				const Scalar adxbdy = adx * bdy;
				const Scalar bdxcdy = bdx * cdy;
				const Scalar bdxady = bdx * ady;
				const Scalar cdxbdy = cdx * bdy;
				const Scalar cdxady = cdx * ady;

				const Scalar alift = adx * adx + ady * ady;
				const Scalar blift = bdx * bdx + bdy * bdy;
				const Scalar clift = cdx * cdx + cdy * cdy;

				const Scalar dif_bdxcdy_cdxbdy = bdxcdy - cdxbdy;
				const Scalar sum_abs_bdxcdy_cdxbdy = std::abs(bdxcdy) + std::abs(cdxbdy);

				const Scalar det_a = alift * dif_bdxcdy_cdxbdy;
				const Scalar det_b = blift * (cdxady - adxcdy);
				const Scalar det_c = clift * (adxbdy - bdxady);

				const Scalar det = det_a + det_b + det_c;

				const Scalar permanent = sum_abs_bdxcdy_cdxbdy * alift + (std::abs(cdxady) + std::abs(adxcdy)) * blift + (std::abs(adxbdy) + std::abs(bdxady)) * clift;

				Scalar errbound = iccerrboundA * permanent;
				if (std::abs(det) >= std::abs(errbound))
					return det <= 0;
			}

			return predicates::adaptive::incircle(
					   p.x, p.y,
					   this->v[0]->x, this->v[0]->y,
					   this->v[1]->x, this->v[1]->y,
					   this->v[2]->x, this->v[2]->y) <= 0;
		}
	};

	Vert *vert_alloc;
	Face *face_alloc;
	Integer *vert_map;
	Integer max_verts;
	Integer max_faces;

	Face *first_dela_face;
	Face *first_hull_face;
	Vert *first_boundary_vert;
	Vert *first_internal_vert;

	Integer inp_verts;
	Integer out_verts;
	Integer polygons;
	Integer out_hull_faces;
	Integer out_boundary_verts;
	Integer unique_points;

	Scalar trans[2]; // experimental

	int (*errlog_proc)(void *file, const char *fmt, ...);
	void *errlog_file;

	Integer Prepare(Integer *start, Face **hull, Integer *out_hull_faces, Face **cache, uint64_t *sort_stamp, Integer stop)
	{
		//uint64_t time0 = uSec();

		//if (errlog_proc)
		//	errlog_proc(errlog_file, "[...] sorting vertices");

		Integer points = inp_verts;

		if (stop >= 0 && points > stop)
			points = stop;

		#if 0
		struct CMP
		{
			bool operator () (const Vert& l, const Vert& r) const
			{
				// reversing paraboloid, somewhat faster without predicate
				Scalar ax = l.x + tx;
				Scalar ay = l.y + ty;
				Scalar bx = r.x + tx;
				Scalar by = r.y + ty;
				Scalar a = ax * ax + ay * ay;
				Scalar b = bx * bx + by * by;
				if (a == b)
				{
					if (l.x > r.x || l.x == r.x && l.y > r.y)
						return true;
					return false;
				}
				return a > b;
			}

			const Scalar tx, ty;
		} cmp = { trans[0], trans[1] };

		std::sort(vert_alloc, vert_alloc + points, cmp);
		#endif

		// rmove dups
		{
			vert_map[vert_alloc[0].i] = 0;

			Integer w = 0, r = 1; // skip initial no-dups block
			while (r < points && !Vert::overlap(vert_alloc + r, vert_alloc + w))
			{
				vert_map[vert_alloc[r].i] = r;
				w++;
				r++;
			}

			Integer d = w; // dup map
			w++;

			while (r < points)
			{
				vert_map[vert_alloc[r].i] = d; // add first dup in run
				r++;

				// skip dups
				while (r < points && Vert::overlap(vert_alloc + r, vert_alloc + r - 1))
				{
					vert_map[vert_alloc[r].i] = d; // add next dup in run
					r++;
				}

				// copy next no-dups block (in percent chunks?)
				while (r < points && !Vert::overlap(vert_alloc + r, vert_alloc + r - 1))
				{
					vert_map[vert_alloc[r].i] = w;
					vert_alloc[w++] = vert_alloc[r++];
				}

				d = w - 1;
			}

			// MAKE CRACK IN THE CRYSTAL
			// - move last point into front
			#if 1
			{
				if (w > 3)
				{
					const Integer tail = w - 1;
					Vert tmp = vert_alloc[tail];
					memmove(vert_alloc + 1, vert_alloc, sizeof(Vert)*tail);
					vert_alloc[0] = tmp;
					for (int i = 0; i < points; i++)
					{
						if (vert_map[i] == tail)
							vert_map[i] = 0;
						else
							vert_map[i]++;
					}
				}
			}
			#endif

			#if 0
			{
				// more cracks
				const int cracks = 1000;
				if (w > cracks)
				{
					Vert crack[cracks];
					for (Integer c = 0; c < cracks; c++)
					{
						Integer d = w - 1 - w * c / cracks;
						Integer b = w - 1 - w * (c + 1) / cracks;

						crack[c] = vert_alloc[d];
						memmove(vert_alloc + b + c + 2, vert_alloc + b + 1, sizeof(Vert)*(d - b - 1));
					}
					memcpy(vert_alloc, crack, sizeof(Vert)*cracks);

					const Integer ofs = cracks * (w - 1) / w;
					for (int i = 0; i < points; i++)
					{
						Integer m = vert_map[i];
						Integer c  = ofs - cracks * m / w;
						Integer cc = ofs - cracks * (m + 1) / w;

						if (cc < c)
							vert_map[i] = c;
						else
							vert_map[i] += c + 1;
					}
					
				}
			}
			#endif


			uint64_t time1 = uSec();

			sorting_bench = time1 - *sort_stamp;

			if (errlog_proc)
				errlog_proc(errlog_file, "\r[100] sorting vertices (%lld ms)\n", (time1 - *sort_stamp) / 1000);
			
			*sort_stamp = time1;

			if (points - w)
			{
				if (errlog_proc)
					errlog_proc(errlog_file, "[WRN] detected %d duplicates in xy array!\n", points - w);
				points = w;
			}
		}

		if (points < 3)
		{
			if (points == 2)
			{
				if (errlog_proc)
					errlog_proc(errlog_file, "[WRN] all input points are colinear, returning single segment!\n");
				first_boundary_vert = vert_alloc + 0;
				first_internal_vert = 0;
				vert_alloc[0].next = vert_alloc + 1;
				vert_alloc[1].next = 0;
			}
			else
			{
				if (errlog_proc)
					errlog_proc(errlog_file, "[WRN] all input points are identical, returning signle point!\n");
				first_boundary_vert = vert_alloc + 0;
				first_internal_vert = 0;
				vert_alloc[0].next = 0;
			}

			out_boundary_verts = points;
			return -points;
		}

		Integer i;
		Face f; // tmp
		f.v[0] = vert_alloc + 0;

		Scalar lo_x = vert_alloc[0].x, hi_x = lo_x;
		Scalar lo_y = vert_alloc[0].y, hi_y = lo_y;
		Integer lower_left = 0;
		Integer upper_right = 0;
		for (i = 1; i < 3; i++)
		{
			Vert *v = vert_alloc + i;
			f.v[i] = v;

			if (v->x < lo_x)
			{
				lo_x = v->x;
				lower_left = i;
			}
			else if (v->x == lo_x && v->y < vert_alloc[lower_left].y)
				lower_left = i;

			if (v->x > hi_x)
			{
				hi_x = v->x;
				upper_right = i;
			}
			else if (v->x == hi_x && v->y > vert_alloc[upper_right].y)
				upper_right = i;

			lo_y = v->y < lo_y ? v->y : lo_y;
			hi_y = v->y > hi_y ? v->y : hi_y;
		}

		int pro = 0;
		// skip until points are coplanar
		while (i < points && f.dot0(vert_alloc[i]))
		{
			if (i >= pro)
			{
				uint64_t p = (int)((uint64_t)100 * i / points);
				pro = (int)((p + 1) * points / 100);
				if (pro >= points)
					pro = (int)points - 1;
				if (i == points - 1)
				{
					p = 100;
				}
				if (errlog_proc)
					errlog_proc(errlog_file, "\r[%2d%s] convex hull triangulation ", p, p >= 100 ? "" : "%");
			}

			Vert *v = vert_alloc + i;

			if (v->x < lo_x)
			{
				lo_x = v->x;
				lower_left = i;
			}
			else if (v->x == lo_x && v->y < vert_alloc[lower_left].y)
				lower_left = i;

			if (v->x > hi_x)
			{
				hi_x = v->x;
				upper_right = i;
			}
			else if (v->x == hi_x && v->y > vert_alloc[upper_right].y)
				upper_right = i;

			lo_y = v->y < lo_y ? v->y : lo_y;
			hi_y = v->y > hi_y ? v->y : hi_y;

			i++;
		}

		Integer *vert_sub = 0;

		bool colinear = f.sign0(); // hybrid
		if (colinear)
		{
			vert_sub = (Integer *)malloc(sizeof(Integer) * ((size_t)i + 1));
			if (!vert_sub)
			{
				if (errlog_proc)
					errlog_proc(errlog_file, "[ERR] Not enough memory, shop for some more RAM. See you!\n");
				return 0;
			}

			for (Integer s = 0; s <= i; s++)
				vert_sub[s] = s;

			// choose x or y axis to sort verts (no need to be exact)
			if (hi_x - lo_x > hi_y - lo_y)
			{
				struct
				{
					bool operator()(const Integer &a, const Integer &b)
					{
						return less(vert_alloc[a], vert_alloc[b]);
					}

					bool less(const Vert &a, const Vert &b) const
					{
						return a.x < b.x;
					}

					Vert *vert_alloc;
				} c;

				c.vert_alloc = vert_alloc;
				std::sort(vert_sub, vert_sub + i, c);
			}
			else
			{
				struct
				{
					bool operator()(const Integer &a, const Integer &b)
					{
						return less(vert_alloc[a], vert_alloc[b]);
					}

					bool less(const Vert &a, const Vert &b) const
					{
						return a.y < b.y;
					}

					Vert *vert_alloc;
				} c;

				c.vert_alloc = vert_alloc;
				std::sort(vert_sub, vert_sub + i, c);
			}
		}
		else
		{
			// split to lower (below diag) and uppert (above diag) parts,
			// sort parts separately in opposite directions
			// mark part with Vert::sew temporarily

			for (Integer j = 0; j < i; j++)
			{
				if (j == lower_left)
				{
					// lower
					vert_alloc[j].sew = &f;
				}
				else if (j == upper_right)
				{
					// upper
					vert_alloc[j].sew = 0;
				}
				else
				{
					Vert *ll = vert_alloc + lower_left;
					Vert *ur = vert_alloc + upper_right;

					Scalar dot = predicates::adaptive::orient2d(ll->x, ll->y, ur->x, ur->y, vert_alloc[j].x, vert_alloc[j].y);

					if (dot < 0)
					{
						// lower
						vert_alloc[j].sew = &f;
					}
					else
					{
						if (dot == 0)
						{
							int dbg = 1;
						}
#ifdef DELABELLA_AUTOTEST
						assert(dot > 0);
#endif
						// upper
						vert_alloc[j].sew = 0;
					}
				}
			}

			struct
			{
				// default to CW order (unlikely, if wrong, we will reverse using one=2 and two=1)

				bool operator()(const Integer &a, const Integer &b) const
				{
					return less(vert_alloc[a], vert_alloc[b]);
				}

				bool less(const Vert &a, const Vert &b) const
				{
					// if a is lower and b is upper, return true
					if (a.sew && !b.sew)
						return false;

					// if a is upper and b is lower, return false
					if (!a.sew && b.sew)
						return true;

					// actually we can compare coords directly
					if (a.sew)
					{
						// lower
						if (a.x > b.x)
							return true;
						if (a.x == b.x)
							return a.y > b.y;
						return false;
					}
					else
					{
						// upper
						if (a.x < b.x)
							return true;
						if (a.x == b.x)
							return a.y < b.y;
						return false;
					}

// otherwise
#ifdef DELABELLA_AUTOTEST
					assert(0);
#endif
					return false;
				}

				Vert *vert_alloc;
			} c;

			vert_sub = (Integer *)malloc(sizeof(Integer) * ((size_t)i + 1));
			if (!vert_sub)
			{
				if (errlog_proc)
					errlog_proc(errlog_file, "[ERR] Not enough memory, shop for some more RAM. See you!\n");
				return 0;
			}

			for (Integer s = 0; s <= i; s++)
				vert_sub[s] = s;

			c.vert_alloc = vert_alloc;
			std::sort(vert_sub, vert_sub + i, c);
		}

		*out_hull_faces = 0;

		// alloc faces only if we're going to create them
		if (i < points || !colinear)
		{
			Integer hull_faces = 2 * points - 4;
			*out_hull_faces = hull_faces;

			if (max_faces < hull_faces)
			{
				if (max_faces)
				{
					free(face_alloc);
					// delete[] face_alloc;
				}
				max_faces = 0;

				/*
				try
				{
					face_alloc = new Face[(size_t)hull_faces];
				}
				catch (...)
				{
					face_alloc = 0;
				}
				*/

				face_alloc = (Face *)malloc(sizeof(Face) * hull_faces);

				if (face_alloc)
					max_faces = hull_faces;
				else
				{
					if (errlog_proc)
						errlog_proc(errlog_file, "[ERR] Not enough memory, shop for some more RAM. See you!\n");
					return 0;
				}
			}

			face_alloc[0].flags = 0;
			for (Integer i = 1; i < hull_faces; i++)
			{
				face_alloc[i - 1].next = face_alloc + i;
				face_alloc[i].flags = 0;
			}
			face_alloc[hull_faces - 1].next = 0;

			*cache = face_alloc;
		}

		if (i == points)
		{
			if (errlog_proc)
			{
				// TODO:
				// this should be moved to end of Triangulate()
				// after printing final "[100] ... \n"
				if (colinear)
					errlog_proc(errlog_file, "[WRN] all input points are colinear\n");
				else
					errlog_proc(errlog_file, "[WRN] all input points are cocircular\n");
			}

			if (colinear)
			{
				// link verts into open list
				first_boundary_vert = vert_alloc + vert_sub[0];
				first_internal_vert = 0;
				out_boundary_verts = points;

				for (Integer j = 1; j < points; j++)
					vert_alloc[j - 1].next = vert_alloc + vert_sub[j];
				vert_alloc[points - 1].next = 0;

				free(vert_sub);
				return -points;
			}

			// we're almost done!
			// time to build final flat hull
			Face *next_p = Face::Alloc(cache);
			Face *next_q = Face::Alloc(cache);
			Face *prev_p = 0;
			Face *prev_q = 0;

			for (Integer j = 2; j < i; j++)
			{
				Face *p = next_p;
				p->v[0] = vert_alloc + vert_sub[0];
				p->v[1] = vert_alloc + vert_sub[j - 1];
				p->v[2] = vert_alloc + vert_sub[j];

				// mirrored
				Face *q = next_q;
				q->v[0] = vert_alloc + vert_sub[0];
				q->v[1] = vert_alloc + vert_sub[j];
				q->v[2] = vert_alloc + vert_sub[j - 1];

				if (j < i - 1)
				{
					next_p = Face::Alloc(cache);
					next_q = Face::Alloc(cache);
				}
				else
				{
					next_p = 0;
					next_q = 0;
				}

				p->f[0] = q;
				p->f[1] = next_p ? next_p : q;
				p->f[2] = prev_p ? prev_p : q;

				q->f[0] = p;
				q->f[1] = prev_q ? prev_q : p;
				q->f[2] = next_q ? next_q : p;

				prev_p = p;
				prev_q = q;
			}

			*hull = prev_q;
		}
		else
		{
			// time to build cone hull with i'th vertex at the tip and 0..i-1 verts in the base
			// build cone's base in direction it is invisible to cone's tip!

			f.v[0] = vert_alloc + vert_sub[0];
			f.v[1] = vert_alloc + vert_sub[1];
			f.v[2] = vert_alloc + vert_sub[2];

			int one = 1, two = 2;
			if (!f.dotNP(vert_alloc[vert_sub[i]]))
			{
				// if i-th vert can see the contour we will flip every face
				one = 2;
				two = 1;
			}

			Face *next_p = Face::Alloc(cache);
			Face *prev_p = 0;

			Face *first_q = 0;
			Face *next_q = Face::Alloc(cache);
			Face *prev_q = 0;

			for (Integer j = 2; j < i; j++)
			{
				Face *p = next_p;
				p->v[0] = vert_alloc + vert_sub[0];
				p->v[one] = vert_alloc + vert_sub[j - 1];
				p->v[two] = vert_alloc + vert_sub[j];

				Face *q;

				if (j == 2)
				{
					// first base triangle also build extra tip face
					q = Face::Alloc(cache);
					q->v[0] = vert_alloc + vert_sub[i];
					q->v[one] = vert_alloc + vert_sub[1];
					q->v[two] = vert_alloc + vert_sub[0];

					q->f[0] = p;
					q->f[one] = 0; // LAST_Q;
					q->f[two] = next_q;

					first_q = q;
					prev_q = q;
				}

				q = next_q;
				q->v[0] = vert_alloc + vert_sub[i];
				q->v[one] = vert_alloc + vert_sub[j];
				q->v[two] = vert_alloc + vert_sub[j - 1];

				next_q = Face::Alloc(cache);

				q->f[0] = p;
				q->f[one] = prev_q;
				q->f[two] = next_q;
				prev_q = q;

				p->f[0] = q;

				if (j < i - 1)
				{
					next_p = Face::Alloc(cache);
				}
				else
				{
					// last base triangle also build extra tip face
					q = next_q;
					q->v[0] = vert_alloc + vert_sub[i];
					q->v[one] = vert_alloc + vert_sub[0];
					q->v[two] = vert_alloc + vert_sub[i - 1];

					q->f[0] = p;
					q->f[one] = prev_q;
					q->f[two] = first_q;

					first_q->f[one] = q;

					next_p = 0;
					prev_q = q;
				}

				p->f[one] = next_p ? next_p : q;
				p->f[two] = prev_p ? prev_p : first_q;

				prev_p = p;
			}

			*hull = prev_q;

			i++;
		}

		free(vert_sub);

		*start = i;
		return points;
	}

	Integer Triangulate(Integer *other_faces, uint64_t* sort_stamp, Integer stop)
	{
		Integer i = 0;
		Face *hull = 0;
		Integer hull_faces = 0;
		Face *cache = 0;

		Integer points = Prepare(&i, &hull, &hull_faces, &cache, sort_stamp, stop);
		unique_points = points < 0 ? -points : points;
		if (points <= 0)
		{
			return points;
		}

		/////////////////////////////////////////////////////////////////////////
		// ACTUAL ALGORITHM

		/*
		// perf meters
		int seeds = 0;
		int grows = 0;
		*/

		int pro = 0;
		for (; i < points; i++)
		{
			if (i >= pro)
			{
				uint64_t p = (int)((uint64_t)100 * i / points);
				pro = (int)((p + 1) * points / 100);
				if (pro >= points)
					pro = (int)points - 1;
				if (i == points - 1)
					p = 100;
				if (errlog_proc)
					errlog_proc(errlog_file, "\r[%2d%s] convex hull triangulation ", p, p >= 100 ? "" : "%");
			}

			Vert *q = vert_alloc + i;
			Vert *p = vert_alloc + i - 1;
			Face *f = hull;

			// 1. FIND FIRST VISIBLE FACE
			//    simply iterate around last vertex using last added triange adjecency info
			while (/*++seeds &&*/ f->dotNP(*q))
			{
				f = f->Next(p);
				if (f == hull)
				{
					// printf(".");
					//  if no visible face can be located at last vertex,
					//  let's run through all faces (approximately last to first),
					//  yes this is emergency fallback and should not ever happen.
					f = face_alloc + (intptr_t)2 * i - 4 - 1;
					while (/*++seeds &&*/ f->dotNP(*q))
					{
#ifdef DELABELLA_AUTOTEST

						if (f == face_alloc)
						{
							//dups test
							struct
							{
								bool operator () (const Vert& a, const Vert& b) const
								{
									return a.x == b.x ? a.y < b.y : a.x < b.x;
								}
							} dup;

							std::sort(vert_alloc, vert_alloc + points, dup);
							for (Integer d = 1; d < points; d++)
							{
								assert(vert_alloc[d - 1].x != vert_alloc[d].x ||
									vert_alloc[d - 1].y != vert_alloc[d].y);
							}
						}

						assert(f != face_alloc); // no face is visible? you must be kidding!
#endif
						f--;
					}
				}
			}

			// 2. DELETE VISIBLE FACES & ADD NEW ONES
			//    (we also build silhouette (vertex loop) between visible & invisible faces)

			Integer del = 0;
			Integer add = 0;

			// push first visible face onto stack (of visible faces)
			Face *stack = f;
			f->next = f; // old trick to use list pointers as 'on-stack' markers
			while (stack)
			{
				// pop, take care of last item ptr (it's not null!)
				f = stack;
				stack = (Face *)f->next;
				if (stack == f)
					stack = 0;
				f->next = 0;

				// copy parts of old face that we still need after removal
				Vert *fv[3] = {(Vert *)f->v[0], (Vert *)f->v[1], (Vert *)f->v[2]};
				Face *ff[3] = {(Face *)f->f[0], (Face *)f->f[1], (Face *)f->f[2]};

				// delete visible face
				f->Free(&cache);
				del++;

				// check all 3 neighbors
				for (int e = 0; e < 3; e++)
				{
					Face *n = ff[e];
					if (n && !n->next) // ensure neighbor is not processed yet & isn't on stack
					{
						// if neighbor is not visible we have slihouette edge
						if (/*++grows &&*/ n->dotNP(*q))
						{
							// build face
							add++;

							// ab: given face adjacency [index][],
							// it provides [][2] vertex indices on shared edge (CCW order)
							const static int ab[3][2] = {{1, 2}, {2, 0}, {0, 1}};

							Vert *a = fv[ab[e][0]];
							Vert *b = fv[ab[e][1]];

							Face *s = Face::Alloc(&cache);
							s->v[0] = a;
							s->v[1] = b;
							s->v[2] = q;

							s->f[2] = n;

							// change neighbour's adjacency from old visible face to cone side
							if (n->f[0] == f)
								n->f[0] = s;
							else if (n->f[1] == f)
								n->f[1] = s;
							else if (n->f[2] == f)
								n->f[2] = s;
#ifdef DELABELLA_AUTOTEST
							else
								assert(0);
#endif

							// build silhouette needed for sewing sides in the second pass
							a->sew = s;
							a->next = b;
						}
						else
						{
							// disjoin visible faces
							// so they won't be processed more than once

							if (n->f[0] == f)
								n->f[0] = 0;
							else if (n->f[1] == f)
								n->f[1] = 0;
							else if (n->f[2] == f)
								n->f[2] = 0;
#ifdef DELABELLA_AUTOTEST
							else
								assert(0);
#endif

							// push neighbor face, it's visible and requires processing
							n->next = stack ? stack : n;
							stack = n;
						}
					}
				}
			}

#ifdef DELABELLA_AUTOTEST
			// if add<del+2 hungry hull has consumed some point
			// that means we can't do delaunay for some under precission reasons
			// althought convex hull would be fine with it
			assert(add == del + 2);
#endif

			// 3. SEW SIDES OF CONE BUILT ON SLIHOUTTE SEGMENTS

			hull = face_alloc + (intptr_t)2 * i - 4 + 1; // last added face

			// last face must contain part of the silhouette
			// (edge between its v[0] and v[1])
			Vert *entry = (Vert *)hull->v[0];

			Vert *pr = entry;
			do
			{
				// sew pr<->nx
				Vert *nx = (Vert *)pr->next;
				pr->sew->f[0] = nx->sew;
				nx->sew->f[1] = pr->sew;
				pr = nx;
			} while (pr != entry);
		}

		// printf("seeds: %d, grows: %d\n", seeds, grows);

#ifdef DELABELLA_AUTOTEST
		assert(2 * i - 4 == hull_faces);
#endif

		for (Integer j = 0; j < points; j++)
		{
			vert_alloc[j].next = 0;
			vert_alloc[j].sew = 0;
		}

		Integer others = 0;

		i = 0;
		Face **prev_dela = &first_dela_face;
		Face **prev_hull = &first_hull_face;
		for (Integer j = 0; j < hull_faces; j++)
		{
			Face *f = face_alloc + j;

			// back-link all verts to some_face
			// yea, ~6x times, sorry
			((Vert *)f->v[0])->sew = f;
			((Vert *)f->v[1])->sew = f;
			((Vert *)f->v[2])->sew = f;

			if (f->signN())
			{
				f->index = i;		   // store index in dela list
				f->flags = 0b01000000; // interior
				*prev_dela = f;
				prev_dela = (Face **)&f->next;
				i++;
			}
			else
			{
				f->index = others;	   // store index in hull list (~ to mark it is not dela)
				f->flags = 0b10000000; // hull
				*prev_hull = f;
				prev_hull = (Face **)&f->next;
				others++;
			}
		}

		if (other_faces)
			*other_faces = others;

		*prev_dela = 0;
		*prev_hull = 0;

		// let's trace boudary contour, at least one vertex of first_hull_face
		// must be shared with dela face, find that dela face
		Iter it;
		Vert *v = (Vert *)first_hull_face->v[0];
		Face *t = (Face *)v->StartIterator(&it);
		Face *e = t; // end

		first_boundary_vert = (Vert *)v;
		out_boundary_verts = 1;

		while (1)
		{
			if (t->IsDelaunay())
			{
				int pr = it.around - 1;
				if (pr < 0)
					pr = 2;
				int nx = it.around + 1;
				if (nx > 2)
					nx = 0;

				Face *fpr = (Face *)t->f[pr];
				if (!fpr->IsDelaunay())
				{
					// let's move from: v to t->v[nx]
					v->next = t->v[nx];
					v = (Vert *)v->next;
					if (v == first_boundary_vert)
						break; // lap finished
					out_boundary_verts++;
					t = (Face *)t->StartIterator(&it, nx);
					e = t;
					continue;
				}
			}
			t = (Face *)it.Next();

#ifdef DELABELLA_AUTOTEST
			assert(t != e);
#endif
		}

		// link all other verts into internal list
		first_internal_vert = 0;
		Vert **prev_inter = &first_internal_vert;
		for (Integer j = 0; j < points; j++)
		{
			if (!vert_alloc[j].next)
			{
				Vert *next = vert_alloc + j;
				*prev_inter = next;
				prev_inter = (Vert **)&next->next;
			}
		}

		if (errlog_proc)
			errlog_proc(errlog_file, "\r[100] convex hull triangulation (%lld ms)\n", (uSec() - *sort_stamp) / 1000);

		return 3 * i;
	}

	bool ReallocVerts(Integer points)
	{
		inp_verts = points;
		out_verts = 0;
		polygons = 0;

		first_dela_face = 0;
		first_hull_face = 0;
		first_boundary_vert = 0;

		if (max_verts < points)
		{
			if (max_verts)
			{
				free(vert_map);
				vert_map = 0;

				free(vert_alloc);
				// delete [] vert_alloc;
				vert_alloc = 0;
				max_verts = 0;
			}

			/*
			try
			{
				vert_alloc = new Vert[(size_t)points];
			}
			catch (...)
			{
				vert_alloc = 0;
			}
			*/

			vert_alloc = (Vert *)malloc(sizeof(Vert) * points);

			if (vert_alloc)
				vert_map = (Integer *)malloc(sizeof(Integer) * (size_t)points);

			if (vert_alloc && vert_map)
				max_verts = points;
			else
			{
				if (errlog_proc)
					errlog_proc(errlog_file, "[ERR] Not enough memory, shop for some more RAM. See you!\n");
				return false;
			}
		}

		return true;
	}

	// private
	Face *FindConstraintOffenders(Vert *va, Vert *vb, Face ***ptail, Vert **restart)
	{
		static const int rotate[3][3] = {{0, 1, 2}, {1, 2, 0}, {2, 0, 1}};
		static const int other_vert[3][2] = {{1, 2}, {2, 0}, {0, 1}};

		// returns list of faces!
		// with Face::index replaced with vertex indice opposite to the offending edge
		// (we are about to change triangulation after all so indexes will change as well)
		// and it's safe for detrimination of dela/hull (sign is preserved)

		Face *list = 0;
		Face **tail = &list;

		Iter it;

		Face *first = (Face *)va->StartIterator(&it);
		Face *face = first;

		// find first face around va, containing offending edge
		Vert *v0;
		Vert *v1;
		Face *N = 0;
		int a, b, c;

		while (1)
		{
			if (!face->IsDelaunay())
			{
				face = (Face *)it.Next();
#ifdef DELABELLA_AUTOTEST
				assert(face != first);
#endif
				continue;
			}

			a = it.around;
			b = other_vert[a][0];
			v0 = (Vert *)(face->v[b]);
			if (v0 == vb)
			{
				*tail = 0;
				*ptail = list ? tail : 0;
				return list; // ab is already there
			}

			c = other_vert[a][1];
			v1 = (Vert *)(face->v[c]);
			if (v1 == vb)
			{
				*tail = 0;
				*ptail = list ? tail : 0;
				return list; // ab is already there
			}

			Scalar a0b = predicates::adaptive::orient2d(va->x, va->y, v0->x, v0->y, vb->x, vb->y);
			Scalar a1b = predicates::adaptive::orient2d(va->x, va->y, v1->x, v1->y, vb->x, vb->y);

			if (a0b <= 0 && a1b >= 0)
			{
				// note:
				// check co-linearity only if v0,v1 are pointing
				// to the right direction (from va to vb)

				if (a0b == 0)
				{
					*restart = v0;
					*tail = 0;
					*ptail = list ? tail : 0;
					return list;
				}

				if (a1b == 0)
				{
					*restart = v1;
					*tail = 0;
					*ptail = list ? tail : 0;
					return list;
				}

				// offending edge!
				N = (Face *)face;
				break;
			}

			face = (Face *)it.Next();

#ifdef DELABELLA_AUTOTEST
			assert(face != first);
#endif
		}

		while (1)
		{
			if (a)
			{
				// rotate N->v[] and N->f 'a' times 'backward' such offending edge appears opposite to v[0]
				const int *r = rotate[a];

				Vert *v[3] = {(Vert *)N->v[0], (Vert *)N->v[1], (Vert *)N->v[2]};
				N->v[0] = v[r[0]];
				N->v[1] = v[r[1]];
				N->v[2] = v[r[2]];

				Face *f[3] = {(Face *)N->f[0], (Face *)N->f[1], (Face *)N->f[2]};
				N->f[0] = f[r[0]];
				N->f[1] = f[r[1]];
				N->f[2] = f[r[2]];

				// if (classify)
				{
					// rotate face_bits!
					if (a == 1)
						N->RotateEdgeFlagsCW();
					else // a==2
						N->RotateEdgeFlagsCCW();
				}
			}

			// add edge
			*tail = N;
			tail = (Face **)&N->next;

			// what is our next face?
			Face *F = (Face *)(N->f[0]);
			int d, e, f;

			if (F->f[0] == N)
			{
				d = 0;
				e = 1;
				f = 2;
			}
			else if (F->f[1] == N)
			{
				d = 1;
				e = 2;
				f = 0;
			}
			else
			{
				d = 2;
				e = 0;
				f = 1;
			}

			Vert *vr = (Vert *)(F->v[d]);

			if (vr == vb)
			{
				*restart = 0;
				*tail = 0;
				*ptail = list ? tail : 0;
				return list;
			}

			// is vr above or below ab ?
			Scalar abr = predicates::adaptive::orient2d(va->x, va->y, vb->x, vb->y, vr->x, vr->y);

			if (abr == 0)
			{
				*restart = vr;
				*tail = 0;
				*ptail = list ? tail : 0;
				return list;
			}

			if (abr > 0)
			{
				// above: de edge (a' = f vert)
				a = f;
				v0 = (Vert *)F->v[d];
				v1 = (Vert *)F->v[e];
				N = F;
			}
			else
			{
				// below: fd edge (a' = e vert)
				a = e;
				v0 = (Vert *)F->v[f];
				v1 = (Vert *)F->v[d];
				N = F;
			}
		}

#ifdef DELABELLA_AUTOTEST
		assert(0);
#endif
		*restart = 0;
		*ptail = 0;
		return 0;
	}

	// private:
	bool FindEdgeFaces(const Vert *v1, const Vert *v2, Face *twin[2], int opposed[2])
	{
		static const int prev[3] = {2, 0, 1};
		static const int next[3] = {1, 2, 0};

		Iter it;
		Face *f = (Face *)v1->StartIterator(&it);
		Face *e = f;
		do
		{
			if (f->v[prev[it.around]] == v2)
			{
				// at least one of them must be dela, it will appear as twin[0]
				if (f->IsDelaunay())
				{
					twin[0] = f;
					opposed[0] = next[it.around];
					twin[1] = (Face *)it.Next();
					opposed[1] = prev[it.around];

#ifdef DELABELLA_AUTOTEST
					assert(twin[0]->f[opposed[0]] == twin[1]);
					assert(twin[1]->f[opposed[1]] == twin[0]);
#endif

					return true;
				}

				Face *tmp_f = f;
				int tmp_op = next[it.around];
				f = (Face *)it.Next();

				if (f->IsDelaunay())
				{
					twin[0] = f;
					opposed[0] = prev[it.around];
					twin[1] = tmp_f;
					opposed[1] = tmp_op;

#ifdef DELABELLA_AUTOTEST
					assert(twin[0]->f[opposed[0]] == twin[1]);
					assert(twin[1]->f[opposed[1]] == twin[0]);
#endif

					return true;
				}
			}
			else
				f = (Face *)it.Next();
		} while (f != e);

		return false;
	}

	Integer ConstrainEdges(Integer edges, const Integer *pa, const Integer *pb, size_t advance_bytes)
	{
		if (advance_bytes == 0)
			advance_bytes = 2 * sizeof(Integer);

		int unfixed = 0;

		uint64_t time0 = uSec();

		int pro = 0;
		for (Integer con = 0; con < edges; con++)
		{
			if (con >= pro)
			{
				uint64_t p = (int)((uint64_t)100 * con / edges);
				pro = (int)((p + 1) * edges / 100);
				if (pro >= edges)
					pro = (int)edges - 1;
				if (con == edges - 1)
					p = 100;
				if (errlog_proc)
					errlog_proc(errlog_file, "\r[%2d%s] constraining ", p, p >= 100 ? "" : "%");
			}

			Integer a = *(const Integer *)((const char *)pa + (intptr_t)con * advance_bytes);
			Integer b = *(const Integer *)((const char *)pb + (intptr_t)con * advance_bytes);

			if (!first_dela_face || a == b)
				continue;

			// current
			Vert *va = (Vert *)GetVertexByIndex(a);
			if (!va)
				continue;

			// destination
			Vert *vc = (Vert *)GetVertexByIndex(b);
			if (!vc)
				continue;

			do
			{
				// 1. Make list of offenders
				Face **tail = 0;
				Vert *restart = 0;
				Face *list = FindConstraintOffenders(va, vc, &tail, &restart);

				if (!list)
				{
					if (restart)
					{
						// if (classify)
						{
							Face *nf[2];
							int op[2];

#ifdef DELABELLA_AUTOTEST
							assert(FindEdgeFaces(va, restart, nf, op));
#else
							FindEdgeFaces(va, restart, nf, op);
#endif

							nf[0]->ToggleEdgeFixed(op[0]);
							nf[1]->ToggleEdgeFixed(op[1]);
						}

						va = restart;
						continue;
					}

					// if (classify)
					{
						Face *nf[2];
						int op[2];

#ifdef DELABELLA_AUTOTEST
						assert(FindEdgeFaces(va, vc, nf, op));
#else
						FindEdgeFaces(va, vc, nf, op);
#endif

						nf[0]->ToggleEdgeFixed(op[0]);
						nf[1]->ToggleEdgeFixed(op[1]);
					}
				}

				Face *flipped = 0;

				Vert *vb = restart ? restart : vc;

				// will we relink'em back to dela and replace indexes?

				// 2. Repeatedly until there are no offenders on the list
				// - remove first edge from the list of offenders
				// - if it forms concave quad diagonal, push it back to the list of offenders
				// - otherwise do flip then:
				//   - if flipped diagonal still intersects ab, push it back to the list of offenders
				//   - otherwise push it to the list of flipped edges
				while (list)
				{
					const int a = 0, b = 1, c = 2;

					Face *N = list;
					list = (Face *)N->next;
					N->next = 0;

					Vert *v0 = (Vert *)(N->v[b]);
					Vert *v1 = (Vert *)(N->v[c]);

					Face *F = (Face *)(N->f[a]);
					int d, e, f;

					if (F->f[0] == N)
					{
						d = 0;
						e = 1;
						f = 2;
					}
					else if (F->f[1] == N)
					{
						d = 1;
						e = 2;
						f = 0;
					}
					else
					{
						d = 2;
						e = 0;
						f = 1;
					}

					Vert *v = (Vert *)N->v[0]; // may be not same as global va (if we have had a skip)
					Vert *vr = (Vert *)(F->v[d]);

					// is v,v0,vr,v1 a convex quad?
					Scalar v0r = predicates::adaptive::orient2d(v->x, v->y, v0->x, v0->y, vr->x, vr->y);
					Scalar v1r = predicates::adaptive::orient2d(v->x, v->y, v1->x, v1->y, vr->x, vr->y);

					if (v0r >= 0 || v1r <= 0)
					{
						// CONCAVE CUNT!
						*tail = N;
						tail = (Face **)&N->next;
						continue;
					}

#ifdef DELABELLA_AUTOTEST
					assert(v0r < 0 && v1r > 0);
#endif

					// if (classify)
					{
						if (N->IsEdgeFixed(a))
						{
							unfixed++;
#ifdef DELABELLA_AUTOTEST
							assert(F->IsEdgeFixed(d));
						}
						else
						{
							assert(!F->IsEdgeFixed(d));
#endif
						}
					}

					// it's convex, xa already checked
					// if (l_a0r < r_a0r && l_a1r > r_a1r)
					{
						if (f == 0)
						{
							/*           *                             *
								   va*  / \                      va*  / \
									  \/   \                        \/   \
									  /\ O  \                       /\ O  \
								   v /  \    \v0                 v /  \    \v0
									*----\----*---------*	      *----\----*---------*
								   / \a   \ b/f\      q/	     / \'-,c\   a\      q/
								  /   \  N \/   \  Q  /   -->   /   \f '-, N  \  Q  /
								 /  P  \   /\ F  \   /		   /  P  \  F '-, b\   /
								/p      \c/e \   d\ /		  /p      \e   \d'-,\ /
							   *---------*----+----*		 *---------*----\----*
										v1     \     vr		           v1    \    vr
												*vb                           *vb
							*/

							Face *O = (Face *)N->f[c];
							Face *P = (Face *)(N->f[b]);
							Face *Q = (Face *)(F->f[e]);

							// bad if O==P
							// int p = P->f[0] == N ? 0 : P->f[1] == N ? 1 : 2;
							// we have to look for Nac edge
							int p = P->v[0] == v ? 2 : P->v[1] == v ? 0 : 1;

							// int q = Q->f[0] == F ? 0 : Q->f[1] == F ? 1 : 2;
							// look for Fdf edge
							int q = Q->v[0] == vr ? 2 : Q->v[1] == vr ? 0 : 1;

							// do flip
							N->v[a] = v0;
							N->v[b] = vr;
							N->v[c] = v;
							N->f[a] = F;
							N->f[b] = O;
							N->f[c] = Q;
							F->v[f] = v; // from v0
							F->f[d] = P; // from N
							F->f[e] = N; // from Q
							P->f[p] = F; // from N
							Q->f[q] = N; // from F

							v0->sew = N;
							v1->sew = F;

							// if (classify)
							{
								// TRANSFER bits:
								// N->b -> F->d
								// N->c -> N->b
								// F->e -> N->c

								F->SetEdgeBits(d, N->GetEdgeBits(b));
								N->SetEdgeBits(b, N->GetEdgeBits(c));
								N->SetEdgeBits(c, F->GetEdgeBits(e));

								if (v == va && vr == vb)
								{
									// set N->a & F->e
									N->SetEdgeBits(a, 0b00001001);
									F->SetEdgeBits(e, 0b00001001);
								}
								else
								{
									// clear N->a & F->e
									N->SetEdgeBits(a, 0);
									F->SetEdgeBits(e, 0);
								}
							}
						}
						else // e==0, (or d==0 but we don't care about it)
						{
							/*           *						       *
										/p\						      /p\
									   /   \					     /   \
									  /  P  \					    /  P  \
								   v /       \ v0                v /       \ v0
									*---------*          	      *---------*
								   / \a  N  b/f\         	     / \'-,e    f\
							 va*__/___\_____/___\__*vb --> va*__/___\b '-, F__\__*vb
								 /     \   /  F  \             /     \  N '-, d\
								/   O   \c/e     d\  		  /   O   \a    c'-,\
							   *---------*---------*		 *---------*---------*
									   v1 \       / vr		         v1 \       / vr
										   \  Q  /                       \  Q  /
											\   /						  \   /
											 \q/						   \q/
											  *							    *
							*/

							Face *O = (Face *)N->f[b];

							Face *P = (Face *)(N->f[c]);
							Face *Q = (Face *)(F->f[f]);
							
							// int p = P->f[0] == N ? 0 : P->f[1] == N ? 1 : 2;
							// look for Nba edge
							int p = P->v[0] == v0 ? 2 : P->v[1] == v0 ? 0 : 1;

							//int q = Q->f[0] == F ? 0 : Q->f[1] == F ? 1 : 2;
							// look for Fed edge
							int q = Q->v[0] == v1 ? 2 : Q->v[1] == v1 ? 0 : 1;

							// do flip
							N->v[a] = v1;
							N->v[b] = v;
							N->v[c] = vr;
							N->f[a] = F;
							N->f[b] = Q;
							N->f[c] = O;
							F->v[e] = v; // from v1
							F->f[d] = P; // from N
							F->f[f] = N; // from Q
							P->f[p] = F; // from N
							Q->f[q] = N; // from F

							v0->sew = F;
							v1->sew = N;

							// if (classify)
							{
								// TRANSFER bits:
								// N->c -> F->d
								// N->b -> N->c
								// F->f -> N->b
								// F->e -> F->e

								F->SetEdgeBits(d, N->GetEdgeBits(c));
								N->SetEdgeBits(c, N->GetEdgeBits(b));
								N->SetEdgeBits(b, F->GetEdgeBits(f));

								if (v == va && vr == vb)
								{
									// set N->a & F->f
									N->SetEdgeBits(a, 0b00001001);
									F->SetEdgeBits(f, 0b00001001);
								}
								else
								{
									N->SetEdgeBits(a, 0);
									F->SetEdgeBits(f, 0);
								}
							}
						}

						// if v and vr are on strongly opposite sides of the edge
						// push N's edge back to offenders otherwise push to the new edges

						if (va == v || vr == vb)
						{

							// resolved!
							N->next = flipped;
							flipped = N;
						}
						else
						{
							// check if v and vr are on the same side of a--b
							Scalar abv = predicates::adaptive::orient2d(va->x, va->y, vb->x, vb->y, v->x, v->y);
							Scalar abr = predicates::adaptive::orient2d(va->x, va->y, vb->x, vb->y, vr->x, vr->y);

							if (abv >= 0 && abr >= 0 || abv <= 0 && abr <= 0)
							{
								// resolved
								N->next = flipped;
								flipped = N;
							}
							else
							{
								// unresolved
								*tail = N;
								tail = (Face **)&N->next;
							}
						}
					}
				}

				// 3. Repeatedly until no flip occurs
				// for every edge from new edges list,
				// if 2 triangles sharing the edge violates delaunay criterion
				// do diagonal flip

				while (1)
				{
					bool no_flips = true;
					Face *N = flipped;
					while (N)
					{
						if (N->v[1] == va && N->v[2] == vb || N->v[1] == vb && N->v[2] == va)
						{
							N = (Face *)N->next;
							continue;
						}

						const int a = 0, b = 1, c = 2;

						Vert *v0 = (Vert *)(N->v[b]);
						Vert *v1 = (Vert *)(N->v[c]);

						Face *F = (Face *)(N->f[a]);
						int d, e, f;

						if (F->f[0] == N)
						{
							d = 0;
							e = 1;
							f = 2;
						}
						else if (F->f[1] == N)
						{
							d = 1;
							e = 2;
							f = 0;
						}
						else
						{
							d = 2;
							e = 0;
							f = 1;
						}

						Vert *v = (Vert *)N->v[0];
						Vert *vr = (Vert *)(F->v[d]);

						// fixed by calling F->cross()
						// bool np = N->dotP(*vr);
						// bool fp = F->dotP(*v);
						// assert(np && fp || !np && !fp);

						// can we check if it was flipped last time?
						// if ((N->index & 0x40000000) == 0)
						if (N->dotP(*vr) /* || F->dotP(*v)*/)
						{
							// if (classify)
							{
								// check if N->a or F->d is marked as fixed
								// it should not ever happen -- i think

								// ... but if it does happen
								// remember we've broken fixed edge
								// then display warning right before we return

								if (N->IsEdgeFixed(a))
								{
									printf("IT DOES HAPPEN!\n");
									unfixed++;
#ifdef DELABELLA_AUTOTEST
									assert(F->IsEdgeFixed(d));
								}
								else
								{
									assert(!F->IsEdgeFixed(d));
#endif
								}
							}

							no_flips = false;

							if (f == 0)
							{
								Face *O = (Face *)N->f[c];
								Face *P = (Face *)(N->f[b]);
								Face *Q = (Face *)(F->f[e]);

								// bad if O==P
								// int p = P->f[0] == N ? 0 : P->f[1] == N ? 1 : 2;
								// we have to look for Nac edge
								int p = P->v[0] == v ? 2 : P->v[1] == v ? 0 : 1;

								// int q = Q->f[0] == F ? 0 : Q->f[1] == F ? 1 : 2;
								// look for Fdf edge
								int q = Q->v[0] == vr ? 2 : Q->v[1] == vr ? 0 : 1;

								// do flip
								N->v[a] = v0;
								N->v[b] = vr;
								N->v[c] = v;
								N->f[a] = F;
								N->f[b] = O;
								N->f[c] = Q;
								F->v[f] = v; // from v0
								F->f[d] = P; // from N
								F->f[e] = N; // from Q
								P->f[p] = F; // from N
								Q->f[q] = N; // from F

								v0->sew = N;
								v1->sew = F;

								// if (classify)
								{
									// TRANSFER bits:
									// N->b -> F->d
									// N->c -> N->b
									// F->f -> F->f
									// F->e -> N->c

									F->SetEdgeBits(d, N->GetEdgeBits(b));
									N->SetEdgeBits(b, N->GetEdgeBits(c));
									N->SetEdgeBits(c, F->GetEdgeBits(e));

									// clear N->a & F->e
									N->SetEdgeBits(a, 0);
									F->SetEdgeBits(e, 0);
								}
							}
							else
							{
								Face *O = (Face *)N->f[b];
								Face *P = (Face *)(N->f[c]);
								Face *Q = (Face *)(F->f[f]);

								// int p = P->f[0] == N ? 0 : P->f[1] == N ? 1 : 2;
								// look for Nba edge
								int p = P->v[0] == v0 ? 2 : P->v[1] == v0 ? 0 : 1;

								//int q = Q->f[0] == F ? 0 : Q->f[1] == F ? 1 : 2;
								// look for Fed edge
								int q = Q->v[0] == v1 ? 2 : Q->v[1] == v1 ? 0 : 1;																		

								// do flip
								N->v[a] = v1;
								N->v[b] = v;
								N->v[c] = vr;
								N->f[a] = F;
								N->f[b] = Q;
								N->f[c] = O;
								F->v[e] = v; // from v1
								F->f[d] = P; // from N
								F->f[f] = N; // from Q
								P->f[p] = F; // from N
								Q->f[q] = N; // from F

								v0->sew = F;
								v1->sew = N;

								// if (classify)
								{
									// TRANSFER bits:
									// N->c -> F->d
									// N->b -> N->c
									// F->f -> N->b
									// F->e -> F->e

									F->SetEdgeBits(d, N->GetEdgeBits(c));
									N->SetEdgeBits(c, N->GetEdgeBits(b));
									N->SetEdgeBits(b, F->GetEdgeBits(f));

									// clear N->a & F->f
									N->SetEdgeBits(a, 0);
									F->SetEdgeBits(f, 0);
								}
							}

							// can we un-mark not flipped somehow?
							// N->index &= 0x3fffffff;
							// F->index &= 0x3fffffff;
						}
						else
						{
							// can we mark it as not flipped somehow?
							// N->index |= 0x40000000;
						}

						N = (Face *)N->next;
					}

					if (no_flips)
						break;
				}

				va = restart;
			} while (va);
		}

		// clean up the mess we've made with dela faces list !!!
		Integer hull_faces = 2 * unique_points - 4;
		Face **tail = &first_dela_face;
		Integer index = 0;
		for (Integer i = 0; i < hull_faces; i++)
		{
			if (face_alloc[i].IsDelaunay())
			{
				*tail = face_alloc + i;
				tail = (Face **)&face_alloc[i].next;
				face_alloc[i].index = index++;
			}
		}
		*tail = 0;

		polygons = index;

		if (errlog_proc)
		{
			errlog_proc(errlog_file, "\r[100] constraining (%lld ms)\n", (uSec() - time0) / 1000);
			if (unfixed)
				errlog_proc(errlog_file, "\r[WRN] %d crossing edges detected - unfixed\n", unfixed);
		}

		return polygons;
	}

	Integer FloodFill(bool invert, const typename SDBSimplex **exterior, int depth)
	{
		if (!first_dela_face)
			return 0;

		uint64_t time0 = uSec();

		if (errlog_proc)
			errlog_proc(errlog_file, "[...] flood filling ");

		if (depth <= 0)
			depth = -1;

		const Integer marker = -1;
		const Integer seeded = -2;
		uint8_t fill = invert ? 0b01000000 : 0;

		// 0. seed list with boundary faces
		//    excluding ones with all their boundary edges marked as fixed odd times
		//    - iterate boundary verts, iterate faces around, check boundary and fixed edges

		Face *seed = 0;
		Face *flip = 0;
		Face* flip_tail = 0;

		Vert *v = first_boundary_vert;
		Vert *e = v;
		do
		{
			Iter it;
			Face *hull = (Face *)v->StartIterator(&it);
			Face *dela = (Face *)it.Next();

			while (hull->IsDelaunay())
			{
				hull = dela;
				dela = (Face *)it.Next();
			}

			while (!dela->IsDelaunay())
			{
				hull = dela;
				dela = (Face *)it.Next();
			}

			// ok we're on the boundary dela(around,around+1)

			static const int next[3] = {1, 2, 0};
			static const int inside[3] = {2, 0, 1};

			int odds = 0, bounds = 1;
			if (dela->flags & (0b1000 << inside[it.around]))
				odds++;

			if (!dela->f[it.around]->IsDelaunay())
			{
				if (dela->flags & (0b1000 << it.around))
					odds++;
				bounds++;
			}

			if (!dela->f[next[it.around]]->IsDelaunay())
			{
				if (dela->flags & (0b1000 << next[it.around]))
					odds++;
				bounds++;
			}

			if (odds < bounds)
			{
				// make sure we dont seed same dela twice
				if (dela->index != seeded)
				{
					dela->index = seeded;
					dela->next = seed;
					seed = dela;
				}
			}
			else
			{
				// make sure we dont seed same dela twice
				if (dela->index != seeded)
				{
					dela->index = seeded;
					dela->next = flip;
					flip = dela;
					if (!flip_tail)
						flip_tail = flip;
				}
			}

			v = (Vert *)v->next;
		} while (v != e);

		if (!seed)
		{
			// invert if boundary is fully
			// covered with odd-fixed edges
			seed = flip;
			flip = 0;
			fill ^= 0b01000000;
			flip_tail = 0;

			depth--;
		}

		int acc = 0;
		int pro = 0;
		int tot = (int)out_verts / 3;

		while (seed && depth)
		{
			Face *stack = seed;
			seed = 0;
			// 1. convert seed list to stack
			//    - for every face in list
			//      - overwrite interior/exterior flag with current fill
			//      - mark as proceed index = -1
			Face *f = stack;
			while (f)
			{
				f->flags = (f->flags & 0b00111111) | fill;
				f->index = marker;
				f = (Face *)f->next;
				acc++;

				if (acc >= pro)
				{
					uint64_t p = (int)((uint64_t)100 * acc / tot);
					pro = (int)((p + 1) * tot / 100);
					if (pro >= tot)
						pro = (int)tot - 1;
					if (acc == tot - 1)
					{
						p = 100;
					}
					if (errlog_proc)
						errlog_proc(errlog_file, "\r[%2d%s] flood filling ", p, p >= 100 ? "" : "%");
				}
			}

			if (flip)
			{
				flip_tail->next = seed;
				seed = flip;
				flip = 0;
				flip_tail = 0;
			}			

			// 2. until stack is not empty
			//    - pop 1 face from stack
			//    - check its every neighbor
			//      - if not processed yet (index != -1) and
			//      - if dela face (flag & 0x80 == 0) and
			//      - if shared edge is not odd-fixed then:
			//        - remove from seed list if it is there! (maybe its not frequent, so iterate seed stack)
			//        - push it to stack
			//        - overwrite interior/exterior flag with current fill
			//        - mark as processed index = -1
			//      else
			//        - push face to seed list (it must be possible to remove it, maybe it is not frequent so we can stack it)
			while (stack)
			{
				f = stack;
				stack = (Face *)f->next;

				for (int i = 0; i < 3; i++)
				{
					Face *n = (Face *)f->f[i];
					if (n->index != marker && n->IsDelaunay())
					{
						if ((f->flags & (0b1000 << i)) == 0)
						{
							if (n->index == seeded)
							{
								// can be slow - check / rethink
								// try reverse pushing !!!
								Face **p = &seed;
								Face *s = seed;
								while (s != n)
								{
									p = (Face **)&s->next;
									s = (Face *)s->next;
								}
								*p = (Face *)n->next;
							}

							n->next = stack;
							stack = n;

							n->flags = (n->flags & 0b00111111) | fill;
							n->index = marker;

							acc++;
							if (acc >= pro)
							{
								uint64_t p = (int)((uint64_t)100 * acc / tot);
								pro = (int)((p + 1) * tot / 100);
								if (pro >= tot)
									pro = (int)tot - 1;
								if (acc == tot - 1)
								{
									p = 100;
								}
								if (errlog_proc)
									errlog_proc(errlog_file, "\r[%2d%s] flood filling ", p, p >= 100 ? "" : "%");
							}
						}
						else if (n->index != seeded)
						{
							n->index = seeded;
							n->next = seed;
							seed = n;
						}
					}
				}
			}

			// 3. if seed list is not empty
			//    - flip current fill mode
			//    - jump to 1.

			depth--;
			fill ^= 0b01000000;
		}

		// 4. rebuild face list and assign indexes

		// - run through face_alloc make 3 lists:
		//   - interior
		//   - exterior
		//   - hull

		polygons = out_verts / 3;

		Integer faces = polygons + out_hull_faces;

		first_hull_face = 0;
		first_dela_face = 0;

		Face *last_interior_face = 0;
		Face *first_exterior_face = 0;

		Integer interior = 0;

		// index
		Integer dela = 0;
		Integer hull = 0;

		for (int i = 0; i < faces; i++)
		{
			Face *f = face_alloc + i;

			if (f->IsDelaunay())
			{

				if (f->index != marker)
				{
					f->flags = (f->flags & ~0b01000000) | fill;
				}

				if (f->flags & 0b01000000)
				{
					if (!interior)
					{
						last_interior_face = f;
						f->next = first_exterior_face;
					}
					else
						f->next = first_dela_face;

					interior++;
					first_dela_face = f;
				}
				else
				{
					f->next = first_exterior_face;
					first_exterior_face = f;
					if (interior)
						last_interior_face->next = f;
				}

				f->index = dela++;
			}
			else
			{
				f->next = first_hull_face;
				first_hull_face = f;
				f->index = hull++;
			}
		}

		if (!first_dela_face)
			first_dela_face = first_exterior_face;

		if (exterior)
			*exterior = first_exterior_face;

		if (errlog_proc)
			errlog_proc(errlog_file, "\r[100] flood filling (%lld ms)\n", (uSec() - time0) / 1000);

		return interior;
	}

	Integer Polygonize(const typename SDBSimplex *poly[])
	{
		const Integer marker = (Integer)(-1);

		uint64_t time0 = uSec();
		Face **buf = 0;
		if (!poly)
		{
			buf = (Face **)malloc(sizeof(Face *) * (size_t)out_verts / 3);
			poly = (const typename SDBSimplex **)buf;
			if (!poly)
				return -1;
		}

		// clear poly indices;
		Face *f = first_dela_face;
		while (f)
		{
			f->index = marker;
			f = (Face *)f->next;
		}

		Integer num = 0;
		f = first_dela_face;
		int pro = 0;
		Integer faces = out_verts / 3;
		Integer i = 0;
		while (f)
		{
			if (i >= pro)
			{
				uint64_t p = (int)((uint64_t)100 * i / faces);
				pro = (int)((p + 1) * faces / 100);
				if (pro >= faces)
					pro = (int)faces - 1;
				if (i == faces - 1)
					p = 100;
				if (errlog_proc)
					errlog_proc(errlog_file, "\r[%2d%s] polygonizing ", p, p >= 100 ? "" : "%");
			}
			i++;

			bool new_poly = true;
			Face *next = (Face *)f->next;
			for (int i = 0; i < 3; i++)
			{
				// if (classify)
				{
					if (f->IsEdgeFixed(i))
						continue;
				}

				Face *a = (Face *)f->f[i];

				Integer index = a->index;

				if (index != marker && a->IsDelaunay())
				{
					int j = 0;
					for (; j < 3; j++)
					{
						Vert *v = (Vert *)a->v[j];
						if (v != f->v[0] && v != f->v[1] && v != f->v[2] && !f->dot0(*v))
							break;
					}

					if (j == 3)
					{
						Integer dest = f->index;
						if (dest != marker)
						{
							// merging polys !!!
							Face *m = (Face *)poly[index];
							while (m)
							{
								Face *n = (Face *)m->next;
								m->index = dest;
								m->next = (Face *)poly[dest];
								poly[dest] = m;
								m = n;
							}

							// fill the gap
							if (index < num - 1)
							{
								dest = index;
								index = num - 1;
								poly[dest] = 0;

								m = (Face *)poly[index];
								while (m)
								{
									Face *n = (Face *)m->next;
									m->index = dest;
									m->next = (Face *)poly[dest];
									poly[dest] = m;
									m = n;
								}
							}

							num--;
						}
						else
						{
							new_poly = false;
							// merge f with existing poly
							f->index = index;
							f->next = (Face *)poly[index];
							poly[index] = f;
						}
					}
				}
			}

			if (new_poly)
			{
				// start new poly
				f->next = 0;
				f->index = num;
				poly[num++] = f;
			}

			f = next;
		}

		polygons = num;

#if 1
		// ALTER POST PROC:
		// re-order triangles in every polygon and re-order indices in faces:
		// - first face defines first 3 verts in contour with v[0],v[1],v[2]
		// - every next face define just 1 additional vertex at its v[0]
		// thay can form fan or strip or arbitraty mix of both
		// without changing existing edges!

		for (Integer p = num - 1; p >= 0; p--)
		{
			Face *f = (Face *)(poly[p]);
			if (!f->next)
			{
				// single triangle
				// just link into list of polys
				if (p < num - 1)
					f->next = (Face *)poly[p + 1];
				continue;
			}

			Face *first = 0;
			// break list appart
			// so we can unmark all poly faces

			while (f)
			{
				Face *n = (Face *)f->next;

				if (!first)
				{
					// lookup good starting face (it will appear as last in the poly after all)
					// having exactly 2 poly boundary edges, exactly 1 inner edge

					int inner_edges =
						(f->f[0]->index == p) +
						(f->f[1]->index == p) +
						(f->f[2]->index == p);

					if (inner_edges == 1)
						first = f;
				}

				f->next = 0; // unmark
				f = n;
			}

#ifdef DELABELLA_AUTOTEST
			assert(first);
#endif

			Face *list = first; // here list is empty, first is used as sentinel

			f = first; // current face

			Face *last = 0; // will be one inserted right after first

			bool step_on = false; // is current vertex inserted

			Iter it;
			f->StartIterator(&it, f->f[0]->index == p ? 2 : f->f[1]->index == p ? 0
																				: 1);

			int dbg = 0;

			while (1)
			{
				if (!step_on && !f->next)
				{
					if (list == first)
						last = f;
					step_on = true;
					f->next = list;
					list = f;
					dbg++;

					// rotate such it.around becomes 0
					if (it.around != 0)
					{
						Face *fr = (Face *)f->f[0];
						Vert *vr = (Vert *)f->v[0];

						if (it.around == 1)
						{
							f->f[0] = f->f[1];
							f->f[1] = f->f[2];
							f->f[2] = fr;
							f->v[0] = f->v[1];
							f->v[1] = f->v[2];
							f->v[2] = vr;

							// if (classify)
							{
								f->RotateEdgeFlagsCW();
							}
						}
						else // it.around == 2
						{
							f->f[0] = f->f[2];
							f->f[2] = f->f[1];
							f->f[1] = fr;
							f->v[0] = f->v[2];
							f->v[2] = f->v[1];
							f->v[1] = vr;

							// if (classify)
							{
								f->RotateEdgeFlagsCCW();
							}
						}

						// adjust iterator after rot
						it.around = 0;
					}
				}

				static const int next_probe[] = {1, 2, 0};
				if (f->f[next_probe[it.around]]->index != p)
				{
					// step on other leg:
					// and restart iterator
					static const int other_leg[3] = {2, 0, 1};
					f->StartIterator(&it, other_leg[it.around]);
					step_on = false;
					continue;
				}

				f = (Face *)it.Next();
				if (f == first)
					break;
			}

			// rotate list so first will be wrapped back to head of the face list
			first->next = list;
			list = first;

			// link last face with first face in next poly
			last->next = (p < num - 1) ? (Face *)poly[p + 1] : 0;

			// store ordered list in poly
			poly[p] = list;
		}

#else
		// merge polys into single list
		// they are separated by indexes
		for (int i = 0; i < num - 1; i++)
		{
			f = (Face *)poly[i];
			while (f->next)
				f = (Face *)f->next;
			f->next = (Face *)poly[i + 1];
		}
#endif

		first_dela_face = (Face *)poly[0];

		if (buf)
			free(buf);

		if (errlog_proc)
			errlog_proc(errlog_file, "\r[100] polygonizing (%lld ms)\n", (uSec() - time0) / 1000);

		return num;
	}

	Integer Triangulate(Integer points, const Scalar *x, const Scalar *y, size_t advance_bytes, Integer stop)
	{
		uint64_t sort_stamp = uSec();

		// const size_t max_triangulate_indices = (size_t)points * 6 - 15;
		// const size_t max_voronoi_edge_indices = (size_t)points * 6 - 12;
		const size_t max_voronoi_poly_indices = (size_t)points * 7 - 9; // winner of shame!

		if (max_voronoi_poly_indices > (size_t)std::numeric_limits<Integer>::max())
		{
			if (errlog_proc)
				errlog_proc(errlog_file, "[ERR] index type too small for provided number of points!\n");
			return 0;
		}

		if (!x)
			return 0;

		if (!y)
			y = x + 1;

		if (advance_bytes < sizeof(Scalar) * 2)
			advance_bytes = sizeof(Scalar) * 2;

		if (!ReallocVerts(points))
			return 0;

		if (errlog_proc)
			errlog_proc(errlog_file, "[...] sorting vertices ");

		for (Integer i = 0; i < points; i++)
		{
			Vert* v = vert_alloc + i;
			v->i = i;
			v->x = *(const Scalar*)((const char*)x + i * advance_bytes);
			v->y = *(const Scalar*)((const char*)y + i * advance_bytes);
		}

		struct KD
		{
			const Scalar xx, xy;
			const Scalar yx, yy;

			Integer pro, acc, tot;
			int (* const errlog_proc)(void *file, const char *fmt, ...);
			void * const errlog_file;

			Scalar bbox[4];

			Scalar X(const Vert& v) const
			{
				return v.x * xx + v.y * xy;
			}

			Scalar Y(const Vert& v) const
			{
				return v.x * yx + v.y * yy;
			}

			void Progress(Integer n)
			{
				acc += n;
				if (acc >= pro)
				{
					uint64_t p = (int)((uint64_t)100 * acc / tot);
					pro = (int)((p + 1) * tot / 100);
					if (pro >= tot)
						pro = (int)tot - 1;
					if (acc == tot - 1)
					{
						p = 100;
					}
					if (errlog_proc)
						errlog_proc(errlog_file, "\r[%2d%s] sorting vertices ", p, p >= 100 ? "" : "%");
				}
			}

			bool Split(Vert* v, Integer n)
			{
				const int limit = 256;

				struct Stk
				{
					struct
					{
						Vert* v;
						Integer n;
					} sub[limit];
					Stk* pop;
					Stk* push;
				};

				int depth = 1;
				Stk base = { {{v,n}}, 0,0 };
				Stk* stack = &base;

				while (1)
				{
					depth--;
					if (depth < 0)
					{
						stack = stack->pop;
						if (!stack)
							break;
						depth = limit - 1;
					}
					v = stack->sub[depth].v;
					n = stack->sub[depth].n;

					if (n < 2)
					{
						Progress(n);
						//return;
						continue;
					}

					bbox[0] = bbox[1] = X(v[0]);
					bbox[2] = bbox[3] = Y(v[0]);
					for (Integer i = 1; i < n; i++)
					{
						Scalar x = X(v[i]);
						Scalar y = Y(v[i]);
						bbox[0] = std::min(bbox[0], x);
						bbox[1] = std::max(bbox[1], x);
						bbox[2] = std::min(bbox[2], y);
						bbox[3] = std::max(bbox[3], y);
					}

					if (bbox[1] - bbox[0] >= bbox[3] - bbox[2] && bbox[0] < bbox[1])
					{
						Scalar half = (bbox[1] + bbox[0]) / 2;
						if (half == bbox[0])
							half = bbox[1];

						Vert* u = v;
						Scalar temp = X(v[0]);
						for (int i = 0; i < n; i++)
						{
							Scalar vx = X(v[i]);
							if (vx >= half)
							{
								u = v + i;
								temp = vx;
							}
						}

						half = temp;

						struct
						{
							Scalar X(const Vert& v) const
							{
								return v.x * xx + v.y * xy;
							}
							bool operator () (const Vert& a, const Vert& b) const
							{
								const Scalar ax = X(a);
								const Scalar bx = X(b);
								if (ax == bx)
								{
									if (a.x == b.x)
										return a.y < b.y;
									return a.x < b.x;
								}
								return ax < bx;
							}
							bool operator () (const Vert& v) const
							{
								const Scalar vx = X(v);
								if (vx == t)
								{
									if (v.x == u.x)
										return v.y < u.y;
									return v.x < u.x;
								}
								return vx < t;

								// return X(v) < t;
							}
							const Vert& u;
							const Scalar t;
							const Scalar xx, xy;
						} p = { *u, half,xx,xy };

						if (bbox[2] == bbox[3] || n == 2)
						{
							if (n>2)
								std::sort(v, v + n, p);
							Progress(n);

							#ifdef DELABELLA_AUTOTEST
							/*
							for (int i = 0; i < n; i++)
							{
								bool dup_allowed = true;
								for (int j = i+1; j < n; j++)
								{
									// all dups can occur only right after i
									bool diff = v[i].x != v[j].x || v[i].y != v[j].y;
									assert(dup_allowed || diff);
									if (dup_allowed && diff)
										dup_allowed = false;
								}
							}
							*/
							#endif
							//return;
							continue;
						}

						u = std::partition(v, v + n, p);

						#ifdef DELABELLA_AUTOTEST
						assert(u != v && u != v + n);
						/*
						for (int i = 0; i < (int)(u - v); i++)
						{
							for (int j = (int)(u - v); j < n; j++)
							{
								// all dups of i must be in the same partition!
								bool diff = v[i].x != v[j].x || v[i].y != v[j].y;
								assert(diff);
							}
						}
						for (int i = (int)(u - v); i < n; i++)
						{
							for (int j = 0; j < (int)(u - v); j++)
							{
								// all dups of i must be in the same partition!
								bool diff = v[i].x != v[j].x || v[i].y != v[j].y;
								assert(diff);
							}
						}
						*/
						#endif

						if (depth == limit)
						{
							depth = 0;
							if (!stack->push)
							{
								// oops
								stack->push = (Stk*)malloc(sizeof(Stk));
								if (!stack->push)
									break;
								stack->push->push = 0;
								stack->push->pop = stack;
							}
							stack = stack->push;
						}
					
						//Split(v, (Integer)(u - v));
						stack->sub[depth].v = v;
						stack->sub[depth].n = (Integer)(u - v);
						depth++;

						if (depth == limit)
						{
							depth = 0;
							if (!stack->push)
							{
								// oops
								stack->push = (Stk*)malloc(sizeof(Stk));
								if (!stack->push)
									break;
								stack->push->push = 0;
								stack->push->pop = stack;
							}
							stack = stack->push;
						}

						//Split(u, n - (Integer)(u - v));
						stack->sub[depth].v = u;
						stack->sub[depth].n = n - (Integer)(u - v);
						depth++;

						continue;
					}
					else
					if (bbox[2] < bbox[3])
					{
						Scalar half = (bbox[3] + bbox[2]) / 2;
						if (half == bbox[2])
							half = bbox[3];

						Vert* u = v;
						Scalar temp = Y(v[0]);
						for (int i = 0; i < n; i++)
						{
							Scalar vy = Y(v[i]);
							if (vy >= half)
							{
								u = v + i;
								temp = vy;
							}
						}

						half = temp;

						struct
						{
							Scalar Y(const Vert& v) const
							{
								return v.x * yx + v.y * yy;
							}
							bool operator () (const Vert& a, const Vert& b) const
							{
								const Scalar ay = Y(a);
								const Scalar by = Y(b);
								if (ay == by)
								{
									if (a.y == b.y)
										return a.x < b.x;
									return a.y < b.y;
								}
								return ay < by;
							}
							bool operator () (const Vert& v) const
							{
								const Scalar vy = Y(v);
								if (vy == t)
								{
									if (v.y == u.y)
										return v.x < u.x;
									return v.y < u.y;
								}
								return vy < t;

								//return Y(v) < t;
							}
							const Vert& u;
							const Scalar t;
							const Scalar yx, yy;
						} p = { *u,half,yx,yy };

						if (bbox[0] == bbox[1] || n == 2)
						{
							if (n > 2)
								std::sort(v, v + n, p);
							Progress(n);

							#ifdef DELABELLA_AUTOTEST
							/*
							for (int i = 0; i < n; i++)
							{
								bool dup_allowed = true;
								for (int j = i+1; j < n; j++)
								{
									// all dups can occur only right after i
									bool diff = v[i].x != v[j].x || v[i].y != v[j].y;
									assert(dup_allowed || diff);
									if (dup_allowed && diff)
										dup_allowed = false;
								}
							}
							*/
							#endif

							//return;
							continue;
						}

						u = std::partition(v, v + n, p);

						#ifdef DELABELLA_AUTOTEST
						assert(u != v && u != v + n);
						/*
						for (int i = 0; i < (int)(u - v); i++)
						{
							for (int j = (int)(u - v); j < n; j++)
							{
								// all dups of i must be in the same partition!
								bool diff = v[i].x != v[j].x || v[i].y != v[j].y;
								if (!diff)
								{
									int a = 0;
								}
								assert(diff);
							}
						}
						for (int i = (int)(u - v); i < n; i++)
						{
							for (int j = 0; j < (int)(u - v); j++)
							{
								// all dups of i must be in the same partition!
								bool diff = v[i].x != v[j].x || v[i].y != v[j].y;
								assert(diff);
							}
						}
						*/
						#endif

						if (depth == limit)
						{
							depth = 0;
							if (!stack->push)
							{
								// oops
								stack->push = (Stk*)malloc(sizeof(Stk));
								if (!stack->push)
									break;
								stack->push->push = 0;
								stack->push->pop = stack;
							}
							stack = stack->push;
						}

						//Split(v, (Integer)(u - v));
						stack->sub[depth].v = v;
						stack->sub[depth].n = (Integer)(u - v);
						depth++;

						if (depth == limit)
						{
							depth = 0;
							if (!stack->push)
							{
								// oops
								stack->push = (Stk*)malloc(sizeof(Stk));
								if (!stack->push)
									break;
								stack->push->push = 0;
								stack->push->pop = stack;
							}
							stack = stack->push;
						}

						//Split(u, n - (Integer)(u - v));
						stack->sub[depth].v = u;
						stack->sub[depth].n = n - (Integer)(u - v);
						depth++;

						continue;
					}
					else
					{
						// all verts are dups
						// dont touch

						// ehm, actually they can appear as dups (inexact rot) 
						// but there may be a non dup hidden in the middle
						struct
						{
							bool operator () (const Vert& a, const Vert& b) const
							{
								if (a.x == b.x)
									return a.y < b.y;
								return a.x < b.x;
							}
						} p;

						if (n > 2)
							std::sort(v, v + n, p);
						Progress(n);

						#ifdef DELABELLA_AUTOTEST
						/*
						for (int i = 0; i < n; i++)
						{
							bool dup_allowed = true;
							for (int j = i+1; j < n; j++)
							{
								// all dups can occur only right after i
								bool diff = v[i].x != v[j].x || v[i].y != v[j].y;
								assert(dup_allowed || diff);
								if (dup_allowed && diff)
									dup_allowed = false;
							}
						}
						*/
						#endif
					}
				} // stack while

				bool ret = stack == 0;

				stack = base.push;
				while (stack)
				{
					Stk* s = stack->push;
					free(stack);
					stack = s;
				}

				return ret;
			}
		};

		KD kd = { (Scalar)2, (Scalar)1, (Scalar)-1, (Scalar)2, 0,0,points, errlog_proc, errlog_file };

		// KD kd = { (Scalar)1.1, (Scalar)0.1, (Scalar)-0.1, (Scalar)1.1 }; // robustenss test
		if (!kd.Split(vert_alloc, points))
		{
			if (errlog_proc)
				errlog_proc(errlog_file, "\n[ERR] Not enough memory, shop for some more RAM. See you!\n");
			return 0;
		}

		#if 0
		exit(0);

		Scalar xlo, xhi, ylo, yhi;
		xlo = xhi = *x;
		ylo = yhi = *y;

		for (Integer i = 0; i < points; i++)
		{
			Vert* v = vert_alloc + i;
			v->i = i;
			v->x = *(const Scalar *)((const char *)x + i * advance_bytes);
			v->y = *(const Scalar *)((const char *)y + i * advance_bytes);

			xlo = std::min(xlo,v->x);
			ylo = std::min(ylo,v->y);
			xhi = std::max(xhi,v->x);
			yhi = std::max(yhi,v->y);
		}

		// todo: scale up more robust!
		Scalar xr = (xhi - xlo) / 2;
		Scalar yr = (yhi - ylo) / 2;

		xlo -= xr;
		xhi += xr;
		ylo -= yr;
		yhi += yr;

		trans[0] = -xlo;
		trans[1] = -ylo;
		#endif


		out_hull_faces = 0;
		unique_points = 0;
		out_verts = Triangulate(&out_hull_faces, &sort_stamp, stop);
		polygons = out_verts / 3;
		return out_verts;
	}

	void Destroy()
	{
		if (vert_map)
			free(vert_map);

		if (face_alloc)
		{
			// delete [] face_alloc;
			free(face_alloc);
		}

		if (vert_alloc)
		{
			// delete [] vert_alloc;
			free(vert_alloc);
		}

		delete this;
	}

	// num of points passed to last call to Triangulate()
	Integer GetNumInputPoints() const
	{
		return inp_verts;
	}

	// num of verts returned from last call to Triangulate()
	Integer GetNumOutputIndices() const
	{
		return out_verts;
	}

	Integer GetNumOutputHullFaces() const
	{
		return out_hull_faces;
	}

	Integer GetNumBoundaryVerts() const
	{
		return out_verts < 0 ? -out_verts : out_boundary_verts;
	}

	Integer GetNumInternalVerts() const
	{
		return out_verts < 0 ? 0 : unique_points - out_boundary_verts;
	}

	// num of polygons
	Integer GetNumPolygons() const
	{
		return polygons;
	}

	const typename SDBSimplex *GetFirstDelaunaySimplex() const
	{
		return first_dela_face;
	}

	const typename SDBSimplex *GetFirstHullSimplex() const
	{
		return first_hull_face;
	}

	const typename SDBVertex *GetFirstBoundaryVertex() const
	{
		return first_boundary_vert;
	}

	const typename SDBVertex *GetFirstInternalVertex() const
	{
		return first_internal_vert;
	}

	const typename SDBVertex *GetVertexByIndex(Integer i) const
	{
		if (i < 0 || i >= inp_verts)
			return 0;
		return vert_alloc + vert_map[i];
	}

	void SetErrLog(int (*proc)(void *stream, const char *fmt, ...), void *stream)
	{
		errlog_proc = proc;
		errlog_file = stream;
	}

	Integer GenVoronoiDiagramVerts(Scalar *x, Scalar *y, size_t advance_bytes) const
	{
		if (!first_dela_face)
			return 0;

		const Integer polys = polygons;
		const Integer contour = out_boundary_verts;
		Integer ret = polys + contour;

		if (!x || !y)
			return ret;

		if (advance_bytes < sizeof(Scalar) * 2)
			advance_bytes = sizeof(Scalar) * 2;

		const Face *f = first_dela_face;
		while (f)
		{
			const Scalar v1x = f->v[1]->x - f->v[0]->x;
			const Scalar v1y = f->v[1]->y - f->v[0]->y;
			const Scalar v2x = f->v[2]->x - f->v[0]->x;
			const Scalar v2y = f->v[2]->y - f->v[0]->y;

			const Scalar v11 = v1x * v1x + v1y * v1y;
			const Scalar v22 = v2x * v2x + v2y * v2y;
			const Scalar v12 = (v1x * v2y - v1y * v2x) * 2;

			Scalar cx = f->v[0]->x + (v2y * v11 - v1y * v22) / v12;
			Scalar cy = f->v[0]->y + (v1x * v22 - v2x * v11) / v12;

			// yes, for polys < tris
			// we possibly calc it multiple times
			// and overwrite already calculated centers
			size_t offs = advance_bytes * (size_t)f->index;
			*(Scalar *)((char *)x + offs) = cx;
			*(Scalar *)((char *)y + offs) = cy;

			f = (Face *)f->next;
		}

		{
			size_t offs = advance_bytes * (size_t)polys;
			x = (Scalar *)((char *)x + offs);
			y = (Scalar *)((char *)y + offs);
		}

		Vert *prev = first_boundary_vert;
		Vert *vert = (Vert *)prev->next;
		for (Integer i = 0; i < contour; i++)
		{
			Scalar nx = prev->y - vert->y;
			Scalar ny = vert->x - prev->x;

			Scalar nn = 1 / sqrt(nx * nx + ny * ny);
			nx *= nn;
			ny *= nn;

			size_t offs = advance_bytes * (size_t)i;
			*(Scalar *)((char *)x + offs) = nx;
			*(Scalar *)((char *)y + offs) = ny;

			prev = vert;
			vert = (Vert *)vert->next;
		}

		return ret;
	}

	Integer GenVoronoiDiagramEdges(Integer *indices, size_t advance_bytes) const
	{
		if (!first_dela_face)
			return 0;

		const Integer polys = polygons;
		const Integer verts = unique_points;
		Integer ret = 2 * (verts + polys - 1);

		if (!indices)
			return ret;

		if (advance_bytes < sizeof(Integer))
			advance_bytes = sizeof(Integer);

		const Integer contour = out_boundary_verts;
		const Integer inter = verts - contour;

		Integer *idx = indices;

		Vert *vert = first_internal_vert;
		for (Integer i = 0; i < inter; i++)
		{
			Iter it;
			Face *t = (Face *)vert->StartIterator(&it);
			Face *e = t;

			Integer a = t->index; // begin
			do
			{
				Integer b = t->index;
				if (a < b)
				{
					*idx = a;
					idx = (Integer *)((char *)idx + advance_bytes);
					*idx = b;
					idx = (Integer *)((char *)idx + advance_bytes);
				}
				a = b;
				t = (Face *)it.Next();
			} while (t != e);

			// loop closing seg
			Integer b = t->index;
			if (a < b)
			{
				*idx = a;
				idx = (Integer *)((char *)idx + advance_bytes);
				*idx = b;
				idx = (Integer *)((char *)idx + advance_bytes);
			}
			a = b;

			vert = (Vert *)vert->next;
		}

		Vert *prev = first_boundary_vert;
		vert = (Vert *)prev->next;
		for (Integer i = 0; i < contour; i++)
		{
			Integer a = i + polys; // begin

			// iterate all dela faces around prev
			// add their voro-vert index == dela face index
			Iter it;
			Face *t = (Face *)prev->StartIterator(&it);

			// it starts at random face, so lookup the prev->vert edge
			while (1)
			{
				if (t->IsDelaunay())
				{
					if (t->v[0] == prev && t->v[1] == vert ||
						t->v[1] == prev && t->v[2] == vert ||
						t->v[2] == prev && t->v[0] == vert)
						break;
				}
				t = (Face *)it.Next();
			}

			// now iterate around, till we're inside the boundary
			while (t->IsDelaunay())
			{
				Integer b = t->index;

				if (a < b)
				{
					*idx = a;
					idx = (Integer *)((char *)idx + advance_bytes);
					*idx = b;
					idx = (Integer *)((char *)idx + advance_bytes);
				}
				a = b;

				t = (Face *)it.Next();
			}

			Integer b = (i == 0 ? contour - 1 : i - 1) + polys; // loop-wrapping!
			if (a < b)
			{
				*idx = a;
				idx = (Integer *)((char *)idx + advance_bytes);
				*idx = b;
				idx = (Integer *)((char *)idx + advance_bytes);
			}
			a = b;

			prev = vert;
			vert = (Vert *)vert->next;
		}

#ifdef DELABELLA_AUTOTEST
		assert(((char *)idx - (char *)indices) / advance_bytes == ret);
#endif

		return ret;
	}

	Integer GenVoronoiDiagramPolys(Integer *indices, size_t advance_bytes, Integer *closed_indices) const
	{
		if (!first_dela_face)
			return 0;

		const Integer polys = polygons;
		const Integer contour = out_boundary_verts;
		const Integer verts = unique_points;
		Integer ret = 3 * verts + 2 * (polys - 1) + contour;

		if (!indices)
			return ret;

		if (advance_bytes < sizeof(Integer))
			advance_bytes = sizeof(Integer);

		const Integer inter = verts - contour;
		const Integer poly_ending = ~0;

		Integer *idx = indices;

		Integer closed = 0;
		Vert *vert = first_internal_vert;
		for (Integer i = 0; i < inter; i++)
		{
			Iter it;
			Face *t = (Face *)vert->StartIterator(&it);
			Face *e = t;

			do
			{
				Integer a = t->index;
				*idx = a;
				idx = (Integer *)((char *)idx + advance_bytes);
				closed++;

				t = (Face *)it.Next();
			} while (t != e);

			// loop closing seg
			Integer a = poly_ending;
			*idx = a;
			idx = (Integer *)((char *)idx + advance_bytes);
			closed++;

			vert = (Vert *)vert->next;
		}

		if (closed_indices)
			*closed_indices = closed;

		Vert *prev = first_boundary_vert;
		vert = (Vert *)prev->next;
		for (Integer i = 0; i < contour; i++)
		{
			Integer a = i + polys; // begin
			*idx = a;
			idx = (Integer *)((char *)idx + advance_bytes);

			// iterate all dela faces around prev
			// add their voro-vert index == dela face index
			Iter it;
			Face *t = (Face *)prev->StartIterator(&it);

			// it starts at random face, so lookup the prev->vert edge
			while (1)
			{
				if (t->IsDelaunay())
				{
					if (t->v[0] == prev && t->v[1] == vert ||
						t->v[1] == prev && t->v[2] == vert ||
						t->v[2] == prev && t->v[0] == vert)
						break;
				}
				t = (Face *)it.Next();
			}

			// now iterate around, till we're inside the boundary
			while (t->IsDelaunay())
			{
				Integer a = t->index;
				*idx = a;
				idx = (Integer *)((char *)idx + advance_bytes);

				t = (Face *)it.Next();
			}

			a = (i == 0 ? contour - 1 : i - 1) + polys; // loop-wrapping!
			*idx = a;
			idx = (Integer *)((char *)idx + advance_bytes);

			a = poly_ending;
			*idx = a;
			idx = (Integer *)((char *)idx + advance_bytes);

			prev = vert;
			vert = (Vert *)vert->next;
		}

#ifdef DELABELLA_AUTOTEST
		assert(((char *)idx - (char *)indices) / advance_bytes == ret);
#endif

		return verts;
	}

	// private:
	void CheckVert(Vert *v) const
	{
		Integer all_faces = out_hull_faces + out_verts / 3;

		assert(v - vert_alloc >= 0);
		assert(v - vert_alloc < unique_points);

		Face *f = (Face *)v->sew;
		assert(f);
		assert(f - face_alloc >= 0);
		assert(f - face_alloc < all_faces);

		int refs = 0;
		for (int i = 0; i < 3; i++)
		{
			if (f->v[i] == v)
				refs++;
		}

		assert(refs == 1);
	}

	void CheckFace(Face *f) const
	{
		Integer all_faces = out_hull_faces + out_verts / 3;

		assert(f - face_alloc >= 0);
		assert(f - face_alloc < all_faces);

		for (int i = 0; i < 3; i++)
		{
			assert(f->v[i]);
			assert((Vert *)f->v[i] - vert_alloc >= 0);
			assert((Vert *)f->v[i] - vert_alloc < unique_points);
		}

		assert(f->v[0] != f->v[1] && f->v[1] != f->v[2] && f->v[2] != f->v[0]);

		for (int i = 0; i < 3; i++)
		{
			assert(f->f[i]);
			assert((Face *)f->f[i] - face_alloc >= 0);
			assert((Face *)f->f[i] - face_alloc < all_faces);
		}

		assert(f->f[0] != f && f->f[1] != f && f->f[2] != f);

		// check back refs & edge flags
		Vert *fv[][2] =
		{
			{(Vert *)f->v[1], (Vert *)f->v[2]},
			{(Vert *)f->v[2], (Vert *)f->v[0]},
			{(Vert *)f->v[0], (Vert *)f->v[1]}
		};

		for (int i = 0; i < 3; i++)
		{

			uint8_t nf = f->GetEdgeBits(i);

			Face *h = (Face *)f->f[i];
			Vert *a = fv[i][0];
			Vert *b = fv[i][1];
			if (h->v[0] == b && h->v[1] == a)
			{
				assert(h->f[2] == f);
				assert(nf == h->GetEdgeBits(2));
			}	
			else 
			if (h->v[1] == b && h->v[2] == a)
			{
				assert(h->f[0] == f);
				assert(nf == h->GetEdgeBits(0));
			}
			else
			{
				assert(h->v[2] == b && h->v[0] == a);
				assert(h->f[1] == f);
				assert(nf == h->GetEdgeBits(1));
			}
		}
	}

	void CheckTopology() const
	{
		assert(first_boundary_vert);
		if (unique_points > out_boundary_verts)
			assert(first_internal_vert);
		assert(first_internal_vert != first_boundary_vert);

		{
			int check = 0;
			Vert *v = first_boundary_vert;
			while (v)
			{
				CheckVert(v);
				check++;
				v = (Vert *)v->next;
				if (v == first_boundary_vert)
					break;
			}
			assert(check == out_boundary_verts);
		}

		{
			int check = 0;
			Vert *v = first_internal_vert;
			while (v)
			{
				CheckVert(v);
				check++;
				v = (Vert *)v->next;
			}
			assert(check == unique_points - out_boundary_verts);
		}

		assert(first_dela_face);
		assert(first_hull_face);
		assert(first_dela_face != first_hull_face);

		{
			int check = 0;
			Face *f = first_dela_face;
			while (f)
			{
				CheckFace(f);
				check++;
				f = (Face *)f->next;
			}
			assert(check == out_verts / 3);
		}

		{
			int check = 0;
			Face *f = first_hull_face;
			while (f)
			{
				CheckFace(f);
				check++;
				f = (Face *)f->next;
			}
			assert(check == out_hull_faces);
		}
	}
};


void SDB::Create()
{
	m_pImpl = new SDBImpl;
}

void SDB::Destroy()
{
	m_pImpl->Destroy();
}

void SDB::SetErrLog( int( *proc )( void* stream, const char* fmt, ... ), void* stream )
{
	m_pImpl->SetErrLog( proc, stream );
}

Integer SDB::Triangulate( Integer points, const Scalar* x, const Scalar* y, size_t advance_bytes, Integer stop )
{
	return m_pImpl->Triangulate( points, x, y, advance_bytes, stop );
}

Integer SDB::GetNumInputPoints() const
{
	return m_pImpl->GetNumInputPoints();
}

Integer SDB::GetNumOutputIndices() const
{
	return m_pImpl->GetNumOutputIndices();
}

Integer SDB::GetNumOutputHullFaces() const
{
	return m_pImpl->GetNumOutputHullFaces();
}

Integer SDB::GetNumBoundaryVerts() const
{
	return m_pImpl->GetNumBoundaryVerts();
}

Integer SDB::GetNumInternalVerts() const
{
	return m_pImpl->GetNumInternalVerts();
}

Integer SDB::GetNumPolygons() const
{
	return m_pImpl->GetNumPolygons();
}

const SDBSimplex* SDB::GetFirstDelaunaySimplex() const
{
	return m_pImpl->GetFirstDelaunaySimplex();
}

const SDBSimplex* SDB::GetFirstHullSimplex() const
{
	return m_pImpl->GetFirstHullSimplex();
}

const SDBVertex* SDB::GetFirstBoundaryVertex() const
{
	return m_pImpl->GetFirstBoundaryVertex();
}

const SDBVertex* SDB::GetFirstInternalVertex() const
{
	return m_pImpl->GetFirstInternalVertex();
}

const SDBVertex* SDB::GetVertexByIndex( Integer i ) const
{
	return m_pImpl->GetVertexByIndex( i );
}

Integer SDB::ConstrainEdges( Integer edges, const Integer* pa, const Integer* pb, size_t advance_bytes )
{
	return m_pImpl->ConstrainEdges( edges, pa, pb, advance_bytes );
}

Integer SDB::FloodFill( bool invert, const SDBSimplex** exterior, int depth )
{
	return m_pImpl->FloodFill( invert, exterior, depth );
}

Integer SDB::Polygonize( const SDBSimplex* poly[/*GetNumOutputIndices()/3*/ ] )
{
	return m_pImpl->Polygonize( poly );
}

Integer SDB::GenVoronoiDiagramVerts( Scalar* x, Scalar* y, size_t advance_bytes ) const
{
	return m_pImpl->GenVoronoiDiagramVerts( x, y, advance_bytes );
}

Integer SDB::GenVoronoiDiagramEdges( Integer* indices, size_t advance_bytes ) const
{
	return m_pImpl->GenVoronoiDiagramEdges( indices, advance_bytes );
}

Integer SDB::GenVoronoiDiagramPolys( Integer* indices, size_t advance_bytes, Integer* closed_indices ) const
{
	return m_pImpl->GenVoronoiDiagramPolys( indices, advance_bytes, closed_indices );
}

void SDB::CheckTopology() const
{
	m_pImpl->CheckTopology();
}

const Scalar SDBImpl::Face::iccerrboundA = ((Scalar(10) + Scalar(96) * std::exp2(-(Scalar)std::numeric_limits<Scalar>::digits)) * std::exp2(-(Scalar)std::numeric_limits<Scalar>::digits));

const Scalar SDBImpl::Vert::resulterrbound = (Scalar(3) + Scalar(8) * std::exp2(-(Scalar)std::numeric_limits<Scalar>::digits)) * std::exp2(-(Scalar)std::numeric_limits<Scalar>::digits);

#ifdef __cplusplus
}
#endif
