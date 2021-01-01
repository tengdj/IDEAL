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

namespace p2t {

struct Edge;

struct TrPoint {

  double x, y;

  /// Default constructor does nothing (for performance).
  TrPoint()
  {
    x = 0.0;
    y = 0.0;
  }

  /// The edges this point constitutes an upper ending point
  std::vector<Edge*> edge_list;

  /// Construct using coordinates.
  TrPoint(double x, double y) : x(x), y(y) {}

  /// Set this point to all zeros.
  void set_zero()
  {
    x = 0.0;
    y = 0.0;
  }

  /// Set this point to some specified coordinates.
  void set(double x_, double y_)
  {
    x = x_;
    y = y_;
  }

  /// Negate this point.
  TrPoint operator -() const
  {
    TrPoint v;
    v.set(-x, -y);
    return v;
  }

  /// Add a point to this point.
  void operator +=(const TrPoint& v)
  {
    x += v.x;
    y += v.y;
  }

  /// Subtract a point from this point.
  void operator -=(const TrPoint& v)
  {
    x -= v.x;
    y -= v.y;
  }

  /// Multiply this point by a scalar.
  void operator *=(double a)
  {
    x *= a;
    y *= a;
  }

  /// Get the length of this point (the norm).
  double Length() const
  {
    return sqrt(x * x + y * y);
  }

  /// Convert this point into a unit point. Returns the Length.
  double Normalize()
  {
    double len = Length();
    x /= len;
    y /= len;
    return len;
  }

};

// Represents a simple polygon's edge
struct Edge {

  TrPoint* p, *q;

  /// Constructor
  Edge(TrPoint& p1, TrPoint& p2) : p(&p1), q(&p2)
  {
    if (p1.y > p2.y) {
      q = &p1;
      p = &p2;
    } else if (p1.y == p2.y) {
      if (p1.x > p2.x) {
        q = &p1;
        p = &p2;
      } else if (p1.x == p2.x) {
        // Repeat points
        assert(false);
      }
    }

    q->edge_list.push_back(this);
  }
};

// Triangle-based data structures are know to have better performance than quad-edge structures
// See: J. Shewchuk, "Triangle: Engineering a 2D Quality Mesh Generator and Delaunay Triangulator"
//      "Triangulations in CGAL"
class Triangle {
public:

/// Constructor
Triangle(TrPoint& a, TrPoint& b, TrPoint& c);

/// Flags to determine if an edge is a Constrained edge
bool constrained_edge[3];
/// Flags to determine if an edge is a Delauney edge
bool delaunay_edge[3];

TrPoint* GetPoint(const int& index);
TrPoint* PointCW(TrPoint& point);
TrPoint* PointCCW(TrPoint& point);
TrPoint* OppositePoint(Triangle& t, TrPoint& p);

Triangle* GetNeighbor(const int& index);
void MarkNeighbor(TrPoint* p1, TrPoint* p2, Triangle* t);
void MarkNeighbor(Triangle& t);

void MarkConstrainedEdge(const int index);
void MarkConstrainedEdge(Edge& edge);
void MarkConstrainedEdge(TrPoint* p, TrPoint* q);

int Index(const TrPoint* p);
int EdgeIndex(const TrPoint* p1, const TrPoint* p2);

Triangle* NeighborCW(TrPoint& point);
Triangle* NeighborCCW(TrPoint& point);
bool GetConstrainedEdgeCCW(TrPoint& p);
bool GetConstrainedEdgeCW(TrPoint& p);
void SetConstrainedEdgeCCW(TrPoint& p, bool ce);
void SetConstrainedEdgeCW(TrPoint& p, bool ce);
bool GetDelunayEdgeCCW(TrPoint& p);
bool GetDelunayEdgeCW(TrPoint& p);
void SetDelunayEdgeCCW(TrPoint& p, bool e);
void SetDelunayEdgeCW(TrPoint& p, bool e);

bool Contains(TrPoint* p);
bool Contains(const Edge& e);
bool Contains(TrPoint* p, TrPoint* q);
void Legalize(TrPoint& point);
void Legalize(TrPoint& opoint, TrPoint& npoint);
/**
 * Clears all references to all other triangles and points
 */
void Clear();
void ClearNeighbor(Triangle *triangle );
void ClearNeighbors();
void ClearDelunayEdges();

inline bool IsInterior();
inline void IsInterior(bool b);

Triangle& NeighborAcross(TrPoint& opoint);

void DebugPrint();

TrPoint *point(int p){
	return points_[p];
}
private:

/// Triangle points
TrPoint* points_[3];
/// Neighbor list
Triangle* neighbors_[3];

/// Has this triangle been marked as an interior triangle?
bool interior_;
};

inline bool cmp(const TrPoint* a, const TrPoint* b)
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
inline TrPoint operator +(const TrPoint& a, const TrPoint& b)
{
  return TrPoint(a.x + b.x, a.y + b.y);
}

/// Subtract two points_ component-wise.
inline TrPoint operator -(const TrPoint& a, const TrPoint& b)
{
  return TrPoint(a.x - b.x, a.y - b.y);
}

/// Multiply point by scalar
inline TrPoint operator *(double s, const TrPoint& a)
{
  return TrPoint(s * a.x, s * a.y);
}

inline bool operator ==(const TrPoint& a, const TrPoint& b)
{
  return a.x == b.x && a.y == b.y;
}

inline bool operator !=(const TrPoint& a, const TrPoint& b)
{
  return !(a.x == b.x) && !(a.y == b.y);
}

/// Peform the dot product on two vectors.
inline double Dot(const TrPoint& a, const TrPoint& b)
{
  return a.x * b.x + a.y * b.y;
}

/// Perform the cross product on two vectors. In 2D this produces a scalar.
inline double Cross(const TrPoint& a, const TrPoint& b)
{
  return a.x * b.y - a.y * b.x;
}

/// Perform the cross product on a point and a scalar. In 2D this produces
/// a point.
inline TrPoint Cross(const TrPoint& a, double s)
{
  return TrPoint(s * a.y, -s * a.x);
}

/// Perform the cross product on a scalar and a point. In 2D this produces
/// a point.
inline TrPoint Cross(const double s, const TrPoint& a)
{
  return TrPoint(-s * a.y, s * a.x);
}

inline TrPoint* Triangle::GetPoint(const int& index)
{
  return points_[index];
}

inline Triangle* Triangle::GetNeighbor(const int& index)
{
  return neighbors_[index];
}

inline bool Triangle::Contains(TrPoint* p)
{
  return p == points_[0] || p == points_[1] || p == points_[2];
}

inline bool Triangle::Contains(const Edge& e)
{
  return Contains(e.p) && Contains(e.q);
}

inline bool Triangle::Contains(TrPoint* p, TrPoint* q)
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


