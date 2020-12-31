/*
 * triangulate.cpp
 *
 *  Created on: Dec 29, 2020
 *      Author: teng
 */


# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <fstream>
# include <iomanip>
# include <iostream>
#include <math.h>

using namespace std;

double angle_degree ( double x1, double y1, double x2, double y2, double x3, double y3 );
bool between ( double xa, double ya, double xb, double yb, double xc, double yc );
char ch_cap ( char ch );
bool ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
bool collinear ( double xa, double ya, double xb, double yb, double xc,  double yc );
bool diagonal ( int im1, int ip1, int n, int prev[], int next[], double x[], double y[] );
bool diagonalie ( int im1, int ip1, int n, int next[], double x[], double y[] );
bool in_cone ( int im1, int ip1, int n, int prev[], int next[], double x[], double y[] );
bool intersect ( double xa, double ya, double xb, double yb, double xc,  double yc, double xd, double yd );
bool intersect_prop ( double xa, double ya, double xb, double yb, double xc, double yc, double xd, double yd );
bool l4_xor ( bool l1, bool l2 );
double polygon_area ( int n, double x[], double y[] );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double triangle_area ( double xa, double ya, double xb, double yb, double xc, double yc );

//****************************************************************************80

double angle_degree ( double x1, double y1, double x2, double y2, double x3,
  double y3 )

//****************************************************************************80
//
//  Purpose:
//
//    ANGLE_DEGREE returns the degree angle defined by three points.
//
//  Discussion:
//
//        P1
//        /
//       /
//      /
//     /
//    P2--------->P3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2016
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X1, Y1, X2, Y2, X3, Y3, the coordinates of the points
//    P1, P2, P3.
//
//    Output, double ANGLE_DEGREE, the angle swept out by the rays, measured
//    in degrees.  0 <= VALUE < 360.  If either ray has zero length,
//    then VALUE is set to 0.
//
{
  double r8_pi = 3.141592653589793;
  double value;
  double x;
  double y;

  x = ( x3 - x2 ) * ( x1 - x2 ) + ( y3 - y2 ) * ( y1 - y2 );

  y = ( x3 - x2 ) * ( y1 - y2 ) - ( y3 - y2 ) * ( x1 - x2 );

  if ( x == 0.0 && y == 0.0 )
  {
    value = 0.0;
  }
  else
  {
    value = atan2 ( y, x );

    if ( value < 0.0 )
    {
      value = value + 2.0 * r8_pi;
    }
    value = 180.0 * value / r8_pi;
  }

  return value;
}
//****************************************************************************80

bool between ( double xa, double ya, double xb, double yb, double xc,
  double yc )

//****************************************************************************80
//
//  Purpose:
//
//    BETWEEN is TRUE if vertex C is between vertices A and B.
//
//  Discussion:
//
//    The points must be (numerically) collinear.
//
//    Given that condition, we take the greater of XA - XB and YA - YB
//    as a "scale" and check where C's value lies.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 May 2014
//
//  Author:
//
//    Original C version by Joseph ORourke.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Joseph ORourke,
//    Computational Geometry in C,
//    Cambridge, 1998,
//    ISBN: 0521649765,
//    LC: QA448.D38.
//
//  Parameters:
//
//    Input, double XA, YA, XB, YB, XC, YC, the coordinates of
//    the vertices.
//
//    Output, bool BETWEEN, is TRUE if C is between A and B.
//
{
  bool value;
  double xmax;
  double xmin;
  double ymax;
  double ymin;

  if ( ! collinear ( xa, ya, xb, yb, xc, yc ) )
  {
    value = false;
  }
  else if ( fabs ( ya - yb ) < fabs ( xa - xb ) )
  {
    xmax = r8_max ( xa, xb );
    xmin = r8_min ( xa, xb );
    value = ( xmin <= xc && xc <= xmax );
  }
  else
  {
    ymax = r8_max ( ya, yb );
    ymin = r8_min ( ya, yb );
    value = ( ymin <= yc && yc <= ymax );
  }

  return value;
}
//****************************************************************************80

char ch_cap ( char ch )

//****************************************************************************80
//
//  Purpose:
//
//    CH_CAP capitalizes a single character.
//
//  Discussion:
//
//    This routine should be equivalent to the library "toupper" function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 July 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH, the character to capitalize.
//
//    Output, char CH_CAP, the capitalized character.
//
{
  if ( 97 <= ch && ch <= 122 )
  {
    ch = ch - 32;
  }

  return ch;
}
//****************************************************************************80

bool ch_eqi ( char ch1, char ch2 )

//****************************************************************************80
//
//  Purpose:
//
//    CH_EQI is true if two characters are equal, disregarding case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH1, CH2, the characters to compare.
//
//    Output, bool CH_EQI, is true if the two characters are equal,
//    disregarding case.
//
{
  if ( 97 <= ch1 && ch1 <= 122 )
  {
    ch1 = ch1 - 32;
  }
  if ( 97 <= ch2 && ch2 <= 122 )
  {
    ch2 = ch2 - 32;
  }

  return ( ch1 == ch2 );
}
//****************************************************************************80

int ch_to_digit ( char ch )

//****************************************************************************80
//
//  Purpose:
//
//    CH_TO_DIGIT returns the integer value of a base 10 digit.
//
//  Example:
//
//     CH  DIGIT
//    ---  -----
//    '0'    0
//    '1'    1
//    ...  ...
//    '9'    9
//    ' '    0
//    'X'   -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH, the decimal digit, '0' through '9' or blank are legal.
//
//    Output, int CH_TO_DIGIT, the corresponding integer value.  If the
//    character was 'illegal', then DIGIT is -1.
//
{
  int digit;

  if ( '0' <= ch && ch <= '9' )
  {
    digit = ch - '0';
  }
  else if ( ch == ' ' )
  {
    digit = 0;
  }
  else
  {
    digit = -1;
  }

  return digit;
}
//****************************************************************************80

bool collinear ( double xa, double ya, double xb, double yb, double xc,
  double yc )

//****************************************************************************80
//
//  Purpose:
//
//    COLLINEAR returns a measure of collinearity for three points.
//
//  Discussion:
//
//    In order to deal with collinear points whose coordinates are not
//    numerically exact, we compare the area of the largest square
//    that can be created by the line segment between two of the points
//    to (twice) the area of the triangle formed by the points.
//
//    If the points are collinear, their triangle has zero area.
//    If the points are close to collinear, then the area of this triangle
//    will be small relative to the square of the longest segment.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2016
//
//  Author:
//
//    Original C version by Joseph ORourke.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Joseph ORourke,
//    Computational Geometry in C,
//    Cambridge, 1998,
//    ISBN: 0521649765,
//    LC: QA448.D38.
//
//  Parameters:
//
//    Input, double XA, YA, XB, YB, XC, YC, the vertex coordinates.
//
//    Output, bool COLLINEAR, is TRUE if the points are judged
//    to be collinear.
//
{
  double area;
  const double r8_eps = 2.220446049250313E-016;
  double side_ab_sq;
  double side_bc_sq;
  double side_ca_sq;
  double side_max_sq;
  bool value;

  area = triangle_area ( xa, ya, xb, yb, xc, yc );

  side_ab_sq = pow ( xa - xb, 2 ) + pow ( ya - yb, 2 );
  side_bc_sq = pow ( xb - xc, 2 ) + pow ( yb - yc, 2 );
  side_ca_sq = pow ( xc - xa, 2 ) + pow ( yc - ya, 2 );

  side_max_sq = r8_max ( side_ab_sq, r8_max ( side_bc_sq, side_ca_sq ) );

  if ( side_max_sq <= r8_eps )
  {
    value = true;
  }
  else if ( 2.0 * fabs ( area ) <= r8_eps * side_max_sq )
  {
    value = true;
  }
  else
  {
    value = false;
  }

  return value;
}
//****************************************************************************80

bool diagonal ( int im1, int ip1, int n, int prev_node[], int next_node[],
  double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIAGONAL: VERTEX(IM1) to VERTEX(IP1) is a proper internal diagonal.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 May 2014
//
//  Author:
//
//    Original C version by Joseph ORourke.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Joseph ORourke,
//    Computational Geometry in C,
//    Cambridge, 1998,
//    ISBN: 0521649765,
//    LC: QA448.D38.
//
//  Parameters:
//
//    Input, int IM1, IP1, the indices of two vertices.
//
//    Input, int N, the number of vertices.
//
//    Input, int PREV_NODE[N], the previous neighbor of each vertex.
//
//    Input, int NEXT_NODE[N], the next neighbor of each vertex.
//
//    Input, double X[N], Y[N], the coordinates of each vertex.
//
//    Output, bool DIAGONAL, the value of the test.
//
{
  bool value;
  bool value1;
  bool value2;
  bool value3;

  value1 = in_cone ( im1, ip1, n, prev_node, next_node, x, y );
  value2 = in_cone ( ip1, im1, n, prev_node, next_node, x, y );
  value3 = diagonalie ( im1, ip1, n, next_node, x, y );

  value = ( value1 && value2 && value3 );

  return value;
}
//****************************************************************************80

bool diagonalie ( int im1, int ip1, int n, int next_node[], double x[],
  double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIAGONALIE is true if VERTEX(IM1):VERTEX(IP1) is a proper diagonal.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 May 2014
//
//  Author:
//
//    Original C version by Joseph ORourke.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Joseph ORourke,
//    Computational Geometry in C,
//    Cambridge, 1998,
//    ISBN: 0521649765,
//    LC: QA448.D38.
//
//  Parameters:
//
//    Input, int IM1, IP1, the indices of two vertices.
//
//    Input, int N, the number of vertices.
//
//    Input, int NEXT_NODE[N], the next neighbor of each vertex.
//
//    Input, double X[N], Y[N], the coordinates of each vertex.
//
//    Output, bool DIAGONALIE, the value of the test.
//
{
  int first;
  int j;
  int jp1;
  bool value;
  bool value2;

  first = im1;
  j = first;
  jp1 = next_node[first];

  value = true;
//
//  For each edge VERTEX(J):VERTEX(JP1) of the polygon:
//
  while ( 1 )
  {
//
//  Skip any edge that includes vertex IM1 or IP1.
//
    if ( j == im1 || j == ip1 || jp1 == im1 || jp1 == ip1 )
    {
    }
    else
    {
      value2 = intersect ( x[im1], y[im1], x[ip1], y[ip1], x[j], y[j],
        x[jp1], y[jp1] );

      if ( value2 )
      {
        value = false;
        break;
      }
    }
    j = jp1;
    jp1 = next_node[j];

    if ( j == first )
    {
      break;
    }
  }

  return value;
}
//****************************************************************************80









bool in_cone ( int im1, int ip1, int n, int prev_node[], int next_node[],
  double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    IN_CONE is TRUE if the diagonal VERTEX(IM1):VERTEX(IP1) is strictly internal.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 May 2014
//
//  Author:
//
//    Original C version by Joseph ORourke.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Joseph ORourke,
//    Computational Geometry in C,
//    Cambridge, 1998,
//    ISBN: 0521649765,
//    LC: QA448.D38.
//
//  Parameters:
//
//    Input, int IM1, IP1, the indices of two vertices.
//
//    Input, int N, the number of vertices.
//
//    Input, int PREV_NODE[N], the previous neighbor of each vertex.
//
//    Input, int NEXT_NODE[N], the next neighbor of each vertex.
//
//    Input, double X[N], Y[N], the coordinates of each vertex.
//
//    Output, bool IN_CONE, the value of the test.
//
{
  int i;
  int im2;
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;
  bool value;

  im2 = prev_node[im1];
  i = next_node[im1];

  t1 = triangle_area ( x[im1], y[im1], x[i], y[i], x[im2], y[im2] );

  if ( 0.0 <= t1 )
  {
    t2 = triangle_area ( x[im1], y[im1], x[ip1], y[ip1], x[im2], y[im2] );
    t3 = triangle_area ( x[ip1], y[ip1], x[im1], y[im1], x[i], y[i] );
    value = ( ( 0.0 < t2 ) && ( 0.0 < t3 ) );
  }
  else
  {
    t4 = triangle_area ( x[im1], y[im1], x[ip1], y[ip1], x[i], y[i] );
    t5 = triangle_area ( x[ip1], y[ip1], x[im1], y[im1], x[im2], y[im2] );
    value = ! ( ( 0.0 <= t4 ) && ( 0.0 <= t5 ) );
  }
  return value;
}
//****************************************************************************80

bool intersect ( double xa, double ya, double xb, double yb, double xc,
  double yc, double xd, double yd )

//****************************************************************************80
//
//  Purpose:
//
//    INTERSECT is true if lines VA:VB and VC:VD intersect.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 May 2014
//
//  Author:
//
//    Original C version by Joseph ORourke.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Joseph ORourke,
//    Computational Geometry in C,
//    Cambridge, 1998,
//    ISBN: 0521649765,
//    LC: QA448.D38.
//
//  Parameters:
//
//    Input, double XA, YA, XB, YB, XC, YC, XD, YD, the X and Y
//    coordinates of the four vertices.
//
//    Output, bool INTERSECT, the value of the test.
//
{
  bool value;

  if ( intersect_prop ( xa, ya, xb, yb, xc, yc, xd, yd ) )
  {
    value = true;
  }
  else if ( between ( xa, ya, xb, yb, xc, yc ) )
  {
    value = true;
  }
  else if ( between ( xa, ya, xb, yb, xd, yd ) )
  {
    value = true;
  }
  else if ( between ( xc, yc, xd, yd, xa, ya ) )
  {
    value = true;
  }
  else if ( between ( xc, yc, xd, yd, xb, yb ) )
  {
    value = true;
  }
  else
  {
    value = false;
  }
  return value;
}
//****************************************************************************80

bool intersect_prop ( double xa, double ya, double xb, double yb, double xc,
  double yc, double xd, double yd )

//****************************************************************************80
//
//  Purpose:
//
//    INTERSECT_PROP is TRUE if lines VA:VB and VC:VD have a proper intersection.
//
//  Discussion:
//
//    Thanks to Gene Dial for correcting the call to intersect_prop(),
//    08 September 2016.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 September 2016
//
//  Author:
//
//    Original C version by Joseph ORourke.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Joseph ORourke,
//    Computational Geometry in C,
//    Cambridge, 1998,
//    ISBN: 0521649765,
//    LC: QA448.D38.
//
//  Parameters:
//
//    Input, double XA, YA, XB, YB, XC, YC, XD, YD, the X and Y
//    coordinates of the four vertices.
//
//    Output, bool INTERSECT_PROP, the result of the test.
//
{
  double t1;
  double t2;
  double t3;
  double t4;
  bool value;
  bool value1;
  bool value2;
  bool value3;
  bool value4;

  if ( collinear ( xa, ya, xb, yb, xc, yc ) )
  {
    value = false;
  }
  else if ( collinear ( xa, ya, xb, yb, xd, yd ) )
  {
    value = false;
  }
  else if ( collinear ( xc, yc, xd, yd, xa, ya ) )
  {
    value = false;
  }
  else if ( collinear ( xc, yc, xd, yd, xb, yb ) )
  {
    value = false;
  }
  else
  {
    t1 = triangle_area ( xa, ya, xb, yb, xc, yc );
    t2 = triangle_area ( xa, ya, xb, yb, xd, yd );
    t3 = triangle_area ( xc, yc, xd, yd, xa, ya );
    t4 = triangle_area ( xc, yc, xd, yd, xb, yb );

    value1 = ( 0.0 < t1 );
    value2 = ( 0.0 < t2 );
    value3 = ( 0.0 < t3 );
    value4 = ( 0.0 < t4 );

    value = ( l4_xor ( value1, value2 ) ) && ( l4_xor ( value3, value4 ) );
  }
  return value;
}
//****************************************************************************80

bool l4_xor ( bool l1, bool l2 )

//****************************************************************************80
//
//  Purpose:
//
//    L4_XOR returns the exclusive OR of two L4's.
//
//  Discussion:
//
//    An L4 is a logical value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 May 2014
//
//  Author:
//
//   John Burkardt
//
//  Parameters:
//
//    Input, bool L1, L2, two values whose exclusive OR is needed.
//
//    Output, bool L4_XOR, the exclusive OR of L1 and L2.
//
{
  bool value;
  bool value1;
  bool value2;

  value1 = (     l1   && ( ! l2 ) );
  value2 = ( ( ! l1 ) &&     l2   );

  value = ( value1 || value2 );

  return value;
}
//****************************************************************************80

double polygon_area ( int n, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYGON_AREA returns the area of a polygon.
//
//  Discussion:
//
//    The vertices should be listed in counter-clockwise order so that
//    the area will be positive.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2016
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, int N, the number of vertices.
//
//    Input, double X[N], Y[N], the vertex coordinates.
//
//    Output, double POLYGON_AREA, the area of the polygon.
//
{
  double area;
  int i;
  int im1;

  area = 0.0;
  im1 = n - 1;
  for ( i = 0; i < n; i++ )
  {
    area = area + x[im1] * y[i] - x[i] * y[im1];
    im1 = i;
  }
  area = 0.5 * area;

  return area;
}
//****************************************************************************80

int *polygon_triangulate ( int n, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYGON_TRIANGULATE determines a triangulation of a polygon.
//
//  Discussion:
//
//    There are N-3 triangles in the triangulation.
//
//    For the first N-2 triangles, the first edge listed is always an
//    internal diagonal.
//
//    Thanks to Gene Dial for pointing out a mistake in the area calculation,
//    10 September 2016.
//
//    Gene Dial requested an angle tolerance of about 1 millionth radian or
//    5.7E-05 degrees, 26 June 2018.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2018
//
//  Author:
//
//    Original C version by Joseph ORourke.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Joseph ORourke,
//    Computational Geometry in C,
//    Cambridge, 1998,
//    ISBN: 0521649765,
//    LC: QA448.D38.
//
//  Parameters:
//
//    Input, int N, the number of vertices.
//
//    Input, double X[N], Y[N], the coordinates of each vertex.
//
//    Output, int TRIANGLES[3*(N-2)], the triangles of the
//    triangulation.
//
{

  double angle;
  const double angle_tol = 5.7E-05;
  double area;
  bool *ear;
  int *next_node;
  int *prev_node;
  int triangle_num;
  int *triangles;

  /*
   *
   * evaluating polygons
   *
   * */
//
//  We must have at least 3 vertices.
//
  if ( n < 3 )
  {
    cerr << "\n";
    cerr << "POLYGON_TRIANGULATE - Fatal error!\n";
    cerr << "  N < 3.\n";
    return NULL;
  }
//
//  Consecutive vertices cannot be equal.
//
  for (int node = 0, node_m1=n-1; node < n; node_m1 = node++ )
  {
    if ( x[node_m1] == x[node] && y[node_m1] == y[node] )
    {
      cerr << "\n";
      cerr << "POLYGON_TRIANGULATE - Fatal error!\n";
      cerr << "  Two consecutive nodes are identical.\n";
      return NULL;
    }
  }
//
//  No node can be the vertex of an angle less than 1 degree
//  in absolute value.
//

  for (int node2 = 0, node1 = n-1; node2 < n; node1=node2++ )
  {
    int node3 = ( ( node2 + 1 ) % n );

    angle = angle_degree (
      x[node1], y[node1],
      x[node2], y[node2],
      x[node3], y[node3] );

    if ( fabs ( angle ) <= angle_tol )
    {
      cerr << "\n";
      cerr << "POLYGON_TRIANGULATE - Fatal error!\n";
      cerr << "  Polygon has an angle smaller than " << angle_tol << "\n";
      cerr << "  occurring at node " << node2 << "\n";
      return NULL;
    }
  }
//
//  Area must be positive.
//
  area = polygon_area ( n, x, y );

  if ( area <= 0.0 )
  {
    cerr << area<<"\n";
    cerr << "POLYGON_TRIANGULATE - Fatal error!\n";
    cerr << "  Polygon has zero or negative area.\n";
    return NULL;
  }

  /*
   *
   * start triangulating
   *
   * */

  triangles = new int[3*(n-2)];
//
//  PREV_NODE and NEXT_NODE point to the previous and next nodes.
//
  prev_node = new int[n];
  next_node = new int[n];

  prev_node[0] = n - 1;
  next_node[0] = 1;

  for (int i = 1; i < n - 1; i++ )
  {
    prev_node[i] = i - 1;
    next_node[i] = i + 1;
  }

  prev_node[n - 1] = n - 2;
  next_node[n - 1] = 0;
//
//  EAR indicates whether the node and its immediate neighbors form an ear
//  that can be sliced off immediately.
//
  ear = new bool[n];
  for (int i = 0; i < n; i++ )
  {
    ear[i] = diagonal ( prev_node[i], next_node[i], n, prev_node, next_node, x, y );
  }

  triangle_num = 0;

  int i0,i1,i2 = 0,i3,i4;

  while ( triangle_num < n - 3 )
  {
//	 cout<< triangle_num<<" "<<i2<<" "<<n-3<<endl;
//
//  If I2 is an ear, gather information necessary to carry out
//  the slicing operation and subsequent "healing".
//
    if ( ear[i2] )
    {
      i3 = next_node[i2];
      i4 = next_node[i3];
      i1 = prev_node[i2];
      i0 = prev_node[i1];
//
//  Make vertex I2 disappear.
//
      next_node[i1] = i3;
      prev_node[i3] = i1;
//
//  Update the earity of I1 and I3, because I2 disappeared.
//
      ear[i1] = diagonal ( i0, i3, n, prev_node, next_node, x, y );
      ear[i3] = diagonal ( i1, i4, n, prev_node, next_node, x, y );
//
//  Add the diagonal [I3, I1, I2] to the list.
//
      triangles[0+triangle_num*3] = i3;
      triangles[1+triangle_num*3] = i1;
      triangles[2+triangle_num*3] = i2;
      triangle_num = triangle_num + 1;
    }
//
//  Try the next vertex.
//
    i2 = next_node[i2];
  }
//
//  The last triangle is formed from the three remaining vertices.
//
  i3 = next_node[i2];
  i1 = prev_node[i2];

  triangles[0+triangle_num*3] = i3;
  triangles[1+triangle_num*3] = i1;
  triangles[2+triangle_num*3] = i2;
  triangle_num = triangle_num + 1;
//
//  Free memory.
//
  delete [] ear;
  delete [] next_node;
  delete [] prev_node;

  return triangles;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = y;
  }
  else
  {
    value = x;
  }
  return value;
}
//****************************************************************************80



double triangle_area ( double xa, double ya, double xb, double yb, double xc, double yc )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_AREA computes the signed area of a triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double XA, YA, XB, YB, XC, YC, the coordinates of
//    the vertices of the triangle, given in counterclockwise order.
//
//    Output, double TRIANGLE_AREA, the signed area of the triangle.
//
{
  double value;

  value = 0.5 * (
      ( xb - xa ) * ( yc - ya )
    - ( xc - xa ) * ( yb - ya ) );

  return value;
}

