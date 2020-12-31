/*
 * triangulate.h
 *
 *  Created on: Dec 29, 2020
 *      Author: teng
 */

#ifndef SRC_GEOMETRY_TRIANGULATE_H_
#define SRC_GEOMETRY_TRIANGULATE_H_


#include "MyPolygon.h"

/* Segment attributes */

class segment_t{
public:
  Point v0, v1;		/* two endpoints */
  int is_inserted;		/* inserted in trapezoidation yet ? */
  int root0, root1;		/* root nodes in Q */
  int next;			/* Next logical segment */
  int prev;			/* Previous segment */
} ;


/* Trapezoid attributes */

class trap_t{
public:
  int lseg, rseg;		/* two adjoining segments */
  Point hi, lo;		/* max/min y-values */
  int u0, u1;
  int d0, d1;
  int sink;			/* pointer to corresponding in Q */
  int usave, uside;		/* I forgot what this means */
  int state;
} ;


/* Node attributes for every node in the query structure */

class node_t{
public:
  int nodetype;			/* Y-node or S-node */
  int segnum;
  Point yval;
  int trnum;
  int parent;			/* doubly linked DAG */
  int left, right;		/* children */
} ;


class monchain_t{
public:
  int vnum;
  int next;			/* Circularly linked list  */
  int prev;			/* describing the monotone */
  int marked;			/* polygon */
};


class vertexchain_t {
public:
  Point pt;
  int vnext[4];			/* next vertices for the 4 chains */
  int vpos[4];			/* position of v in the 4 chains */
  int nextfree;
} ;


/* Node types */

#define T_X     1
#define T_Y     2
#define T_SINK  3


//#define SEGSIZE 200		/* max# of segments. Determines how */
//				/* many points can be specified as */
//				/* input. If your datasets have large */
//				/* number of points, increase this */
//				/* value accordingly. */
//
//#define QSIZE   8*SEGSIZE	/* maximum table sizes */
//#define TRSIZE  4*SEGSIZE	/* max# trapezoids */


#define TRUE  1
#define FALSE 0


#define FIRSTPT 1		/* checking whether pt. is inserted */
#define LASTPT  2


//#define INFINITY 1<<30
#define C_EPS 1.0e-7		/* tolerance value: Used for making */
				/* all decisions about collinearity or */
				/* left/right of segment. Decrease */
				/* this value if the input points are */
				/* spaced very close together */


#define S_LEFT 1		/* for merge-direction */
#define S_RIGHT 2


#define ST_VALID 1		/* for trapezium state */
#define ST_INVALID 2


#define SP_SIMPLE_LRUP 1	/* for splitting trapezoids */
#define SP_SIMPLE_LRDN 2
#define SP_2UP_2DN     3
#define SP_2UP_LEFT    4
#define SP_2UP_RIGHT   5
#define SP_2DN_LEFT    6
#define SP_2DN_RIGHT   7
#define SP_NOSPLIT    -1

#define TR_FROM_UP 1		/* for traverse-direction */
#define TR_FROM_DN 2

#define TRI_LHS 1
#define TRI_RHS 2


#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

#define CROSS(v0, v1, v2) (((v1).x - (v0).x)*((v2).y - (v0).y) - \
			   ((v1).y - (v0).y)*((v2).x - (v0).x))

#define DOT(v0, v1) ((v0).x * (v1).x + (v0).y * (v1).y)

#define FP_EQUAL(s, t) (fabs(s - t) <= C_EPS)

#define CROSS_SINE(v0, v1) ((v0).x * (v1).y - (v1).x * (v0).y)
#define LENGTH(v0) (sqrt((v0).x * (v0).x + (v0).y * (v0).y))

/* utility functions */
inline int _greater_than(Point *v0, Point *v1)
{
  if (v0->y > v1->y + C_EPS)
    return TRUE;
  else if (v0->y < v1->y - C_EPS)
    return FALSE;
  else
    return (v0->x > v1->x);
}


inline int _equal_to(Point *v0, Point *v1)
{
  return (FP_EQUAL(v0->y, v1->y) && FP_EQUAL(v0->x, v1->x));
}

inline int _greater_than_equal_to(Point *v0, Point *v1)
{
  if (v0->y > v1->y + C_EPS)
    return TRUE;
  else if (v0->y < v1->y - C_EPS)
    return FALSE;
  else
    return (v0->x >= v1->x);
}

inline int _less_than(Point *v0, Point *v1)
{
  if (v0->y < v1->y - C_EPS)
    return TRUE;
  else if (v0->y > v1->y + C_EPS)
    return FALSE;
  else
    return (v0->x < v1->x);
}


/* Get log*n for given n */
inline int math_logstar_n(int n)
{
  register int i;
  double v;

  for (i = 0, v = (double) n; v >= 1; i++)
    v = log2(v);

  return (i - 1);
}


inline int math_N(int n, int h)
{
  register int i;
  double v;

  for (i = 0, v = (int) n; i < h; i++)
    v = log2(v);

  return (int) ceil((double) 1.0*n/v);
}

/* Return the maximum of the two points into the yval structure */
inline int _max(Point *yval, Point *v0, Point *v1){
	if (v0->y > v1->y + C_EPS){
		*yval = *v0;
	}else if(FP_EQUAL(v0->y, v1->y)){
		if (v0->x > v1->x + C_EPS)
			*yval = *v0;
		else
			*yval = *v1;
	} else{
		*yval = *v1;
	}
	return 0;
}


/* Return the minimum of the two points into the yval structure */
inline int _min(Point *yval, Point *v0, Point *v1){
	if (v0->y < v1->y - C_EPS){
		*yval = *v0;
	}else if(FP_EQUAL(v0->y, v1->y)){
		if (v0->x < v1->x)
			*yval = *v0;
		else
			*yval = *v1;
	}else{
		*yval = *v1;
	}

	return 0;
}

inline double get_angle(Point *vp0, Point *vpnext, Point *vp1)
{
  Point v0, v1;

  v0.x = vpnext->x - vp0->x;
  v0.y = vpnext->y - vp0->y;

  v1.x = vp1->x - vp0->x;
  v1.y = vp1->y - vp0->y;

  if (CROSS_SINE(v0, v1) >= 0)	/* sine is positive */
    return DOT(v0, v1)/LENGTH(v0)/LENGTH(v1);
  else
    return (-1.0 * DOT(v0, v1)/LENGTH(v0)/LENGTH(v1) - 2);
}


class triangulator{
	int num_edges = 0;
	int num_trap = 0;
	int num_query = 0;

	node_t *qs = NULL;		/* Query structure */
	trap_t *tr = NULL;		/* Trapezoid structure */
	int q_idx = 0;
	int tr_idx = 0;
	segment_t *seg = NULL;		/* Segment table */

	int choose_idx = 0;
	int *permute = NULL;

	monchain_t *mchain = NULL; /* Table to hold all the monotone */
					  /* polygons . Each monotone polygon */
					  /* is a circularly linked list */

	vertexchain_t *vert = NULL; /* chain init. information. This */
					    /* is used to decide which */
					    /* monotone polygon to split if */
					    /* there are several other */
					    /* polygons touching at the same */
					    /* vertex  */

	int *mon = NULL;	/* contains position of any vertex in */
					/* the monotone chain for the polygon */
	int *visited = NULL;
	int chain_idx = 0, op_idx = 0, mon_idx = 0;

	/* Return a new node to be added into the query tree */
	int newnode(){
		return q_idx++;
	}

	/* Return a free trapezoid */
	int newtrap(){
		tr[tr_idx].lseg = -1;
		tr[tr_idx].rseg = -1;
		tr[tr_idx].state = ST_VALID;
		return tr_idx++;
	}

	/* return a new mon structure from the table */
	int newmon(){
	  return ++mon_idx;
	}


	/* return a new chain element from the table */
	int new_chain_element(){
	  return ++chain_idx;
	}

	int locate_endpoint(Point *, Point *, int);

	int generate_random_ordering();
	int choose_segment(void);

	int init_query_structure(int segnum);
	int is_left_of(int segnum, Point *v);
	int inserted(int segnum, int whichpt);
	int merge_trapezoids(int segnum, int tfirst, int tlast, int side);
	int add_segment(int segnum);
	int find_new_roots(int segnum);
	int inside_polygon(trap_t *t);
	int get_vertex_positions(int v0, int v1, int *ip, int *iq);
	int make_new_monotone_poly(int mcur, int v0, int v1);
	int traverse_polygon(int mcur, int trnum, int from, int dir);
	int triangulate_single_polygon(int posmax, int side, int *op);
	int initialise(int n);
public:
	triangulator(){
	}


	int construct_trapezoids();
	int monotonate_trapezoids();
	int triangulate_monotone_polygons(int nmonpoly, int *op);

	int is_point_inside_polygon(Point *v);
	int triangulate_polygon(VertexSequence *vc, int *triangles);

	void print_trapezoids();


};


#endif /* SRC_GEOMETRY_TRIANGULATE_H_ */
