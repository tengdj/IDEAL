/*
 * Poly2Tri Copyright (c) 2009-2010, Poly2Tri Contributors
 * http://code.google.com/p/poly2tri/
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 * * Neither the name of Poly2Tri nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// Include guard
#ifndef SHAPES_H
#define SHAPES_H

#include <vector>
#include <cstddef>
#include <assert.h>
#include <cmath>
#include "../geometry/Point.h"
namespace p2t {

// Triangle-based data structures are know to have better performance than quad-edge structures
// See: J. Shewchuk, "Triangle: Engineering a 2D Quality Mesh Generator and Delaunay Triangulator"
//      "Triangulations in CGAL"
class Triangle {
private:

	/// Triangle points
	Point* points_[3];
	/// Neighbor list
	Triangle* neighbors_[3];

	/// Has this triangle been marked as an interior triangle?
	bool interior_;
public:

	/// Constructor
	Triangle(Point& a, Point& b, Point& c);

	/// Flags to determine if an edge is a Constrained edge
	bool constrained_edge[3];
	/// Flags to determine if an edge is a Delauney edge
	bool delaunay_edge[3];

	Point* GetPoint(const int& index);
	Point* PointCW(Point& point);
	Point* PointCCW(Point& point);
	Point* OppositePoint(Triangle& t, Point& p);

	Triangle* GetNeighbor(const int& index);
	void MarkNeighbor(Point* p1, Point* p2, Triangle* t);
	void MarkNeighbor(Triangle& t);

	void MarkConstrainedEdge(const int index);
	void MarkConstrainedEdge(Edge& edge);
	void MarkConstrainedEdge(Point* p, Point* q);

	int Index(const Point* p);
	int EdgeIndex(const Point* p1, const Point* p2);

	Triangle* NeighborCW(Point& point);
	Triangle* NeighborCCW(Point& point);
	bool GetConstrainedEdgeCCW(Point& p);
	bool GetConstrainedEdgeCW(Point& p);
	void SetConstrainedEdgeCCW(Point& p, bool ce);
	void SetConstrainedEdgeCW(Point& p, bool ce);
	bool GetDelunayEdgeCCW(Point& p);
	bool GetDelunayEdgeCW(Point& p);
	void SetDelunayEdgeCCW(Point& p, bool e);
	void SetDelunayEdgeCW(Point& p, bool e);

	bool Contains(Point* p);
	bool Contains(const Edge& e);
	bool Contains(Point* p, Point* q);
	void Legalize(Point& point);
	void Legalize(Point& opoint, Point& npoint);
	/**
	 * Clears all references to all other triangles and points
	 */
	void Clear();
	void ClearNeighbor(Triangle *triangle );
	void ClearNeighbors();
	void ClearDelunayEdges();

	inline bool IsInterior();
	inline void IsInterior(bool b);

	Triangle& NeighborAcross(Point& opoint);


	void DebugPrint();

	Point *point(int p){
		return points_[p];
	}
	void print(){
		printf("((%f %f, %f %f, %f %f, %f %f))",
				points_[0]->x,points_[0]->y,
				points_[1]->x,points_[1]->y,
				points_[2]->x,points_[2]->y,
				points_[0]->x,points_[0]->y);
	}

};

inline bool cmp(const Point* a, const Point* b)
{
  if (a->y < b->y) {
	return true;
  } else if (a->y == b->y) {
	// Make sure q is point with greater x value
	if (a->x < b->x) {
	  return true;
	}
  }
  return false;
}

/// Add two points_ component-wise.
inline Point operator +(const Point& a, const Point& b)
{
  return Point(a.x + b.x, a.y + b.y);
}

/// Subtract two points_ component-wise.
inline Point operator -(const Point& a, const Point& b)
{
  return Point(a.x - b.x, a.y - b.y);
}

/// Multiply point by scalar
inline Point operator *(double s, const Point& a)
{
  return Point(s * a.x, s * a.y);
}

inline bool operator ==(const Point& a, const Point& b)
{
  return a.x == b.x && a.y == b.y;
}

inline bool operator !=(const Point& a, const Point& b)
{
  return !(a.x == b.x) && !(a.y == b.y);
}

/// Peform the dot product on two vectors.
inline double Dot(const Point& a, const Point& b)
{
  return a.x * b.x + a.y * b.y;
}

/// Perform the cross product on two vectors. In 2D this produces a scalar.
inline double Cross(const Point& a, const Point& b)
{
  return a.x * b.y - a.y * b.x;
}

/// Perform the cross product on a point and a scalar. In 2D this produces
/// a point.
inline Point Cross(const Point& a, double s)
{
  return Point(s * a.y, -s * a.x);
}

/// Perform the cross product on a scalar and a point. In 2D this produces
/// a point.
inline Point Cross(const double s, const Point& a)
{
  return Point(-s * a.y, s * a.x);
}

inline Point* Triangle::GetPoint(const int& index)
{
  return points_[index];
}

inline Triangle* Triangle::GetNeighbor(const int& index)
{
  return neighbors_[index];
}

inline bool Triangle::Contains(Point* p)
{
  return p == points_[0] || p == points_[1] || p == points_[2];
}

inline bool Triangle::Contains(const Edge& e)
{
  return Contains(e.p) && Contains(e.q);
}

inline bool Triangle::Contains(Point* p, Point* q)
{
  return Contains(p) && Contains(q);
}

inline bool Triangle::IsInterior()
{
  return interior_;
}

inline void Triangle::IsInterior(bool b)
{
  interior_ = b;
}

}

#endif


