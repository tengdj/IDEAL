/*
 * geometry_computation.h
 *
 *  Created on: Apr 28, 2022
 *      Author: teng
 */

#ifndef SRC_GEOMETRY_GEOMETRY_COMPUTATION_H_
#define SRC_GEOMETRY_GEOMETRY_COMPUTATION_H_

#include <math.h>

/*
 *
 * some utility functions for geometry computations
 *
 * */

/*
 *
 * calculating distance between point and polygon
 *
 * */
inline double point_to_segment_distance(const double x, const double y, const double x1, double y1, const double x2, const double y2, bool geography = false) {

  double A = x - x1;
  double B = y - y1;
  double C = x2 - x1;
  double D = y2 - y1;

  double dot = A * C + B * D;
  double len_sq = C * C + D * D;
  double param = -1;
  if (len_sq != 0) //in case of 0 length line
      param = dot / len_sq;

  double xx, yy;

  if (param < 0) {
    xx = x1;
    yy = y1;
  } else if (param > 1) {
    xx = x2;
    yy = y2;
  } else {
    xx = x1 + param * C;
    yy = y1 + param * D;
  }

  double dx = x - xx;
  double dy = y - yy;
  if(geography){
      dy = dy/degree_per_kilometer_latitude;
      dx = dx/degree_per_kilometer_longitude(y);
  }
  return sqrt(dx * dx + dy * dy);
}

inline double point_to_segment_distance_batch(Point &p, Point *vs, size_t seq_len, bool geography = false){
    double mindist = DBL_MAX;
    for (int i = 0; i < seq_len-1; i++) {
        double dist = point_to_segment_distance(p.x, p.y, vs[i].x, vs[i].y, vs[i+1].x, vs[i+1].y, geography);
        if(dist<mindist){
            mindist = dist;
        }
    }
    return mindist;
}

/*
 * calculating distance between segments, for distance calculation between polygons
 * */

inline double segment_to_segment_distance_batch(Point *vs1, Point *vs2, size_t s1, size_t s2, bool geography = false){
    double mindist = DBL_MAX;
	for(int i=0;i<s1;i++){
		double dist = point_to_segment_distance_batch(vs1[i], vs2, s2, geography);
		if(dist<mindist){
			mindist = dist;
		}
	}
	for(int i=0;i<s2;i++){
		double dist = point_to_segment_distance_batch(vs2[i], vs1, s1, geography);
		if(dist<mindist){
			mindist = dist;
		}
	}
	return mindist;
}

/*
 *
 * checking whether two polygons intersect
 *
 * */

inline int sgn(const double& x) {
	return x >= 0 ? x ? 1 : 0 : -1;
}

inline bool inter1(double a, double b, double c, double d) {

	double tmp;
    if (a > b){
    	tmp = a;
    	a = b;
    	b = tmp;
    }
    if (c > d){
    	tmp = c;
    	c = d;
    	d = tmp;
    }
    return max(a, c) <= min(b, d);
}

inline bool segment_intersect(const Point& a, const Point& b, const Point& c, const Point& d) {
    if (c.cross(a, d) == 0 && c.cross(b, d) == 0)
        return inter1(a.x, b.x, c.x, d.x) && inter1(a.y, b.y, c.y, d.y);
    return sgn(a.cross(b, c)) != sgn(a.cross(b, d)) &&
           sgn(c.cross(d, a)) != sgn(c.cross(d, b));
}


// checking whether two segments intersect
inline bool segment_intersect_batch(Point *p1, Point *p2, int s1, int s2){
	for(int i=0;i<s1;i++){
		for(int j=0;j<s2;j++){
			if(segment_intersect(p1[i],p1[(i+1)%s1],p2[j],p2[(j+1)%s2])){
//				p1[i].print();
//				p1[(i+1)%s1].print();
//				p2[j].print();
//				p2[(j+1)%s2].print();
//				log("%d %d intersect",i,j);
				return true;
			}
		}
	}
	return false;
}



#endif /* SRC_GEOMETRY_GEOMETRY_COMPUTATION_H_ */
