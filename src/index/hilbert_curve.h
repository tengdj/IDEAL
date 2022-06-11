/*
 * hilbert_curve.hpp
 *
 *  Created on: Jun 8, 2022
 *      Author: teng
 */

#ifndef SRC_INDEX_HILBERT_CURVE_H_
#define SRC_INDEX_HILBERT_CURVE_H_

# include <cstdlib>
# include <ctime>
# include <iomanip>
# include <iostream>

using namespace std;

void rot ( size_t n, size_t &x, size_t &y, size_t rx, size_t ry )

//****************************************************************************80
//
//  Purpose:
//
//    ROT rotates and flips a quadrant appropriately
//
//  Modified:
//
//    24 December 2015
//
//  Parameters:
//
//    Input, size_t N, the length of a side of the square.  N must be a power of 2.
//
//    Input/output, size_t &X, &Y, the input and output coordinates of a point.
//
//    Input, size_t RX, RY, ???
//
{
  size_t t;

  if ( ry == 0 )
  {
//
//  Reflect.
//
    if ( rx == 1 )
    {
      x = n - 1 - x;
      y = n - 1 - y;
    }
//
//  Flip.
//
    t = x;
    x = y;
    y = t;
  }
  return;
}

size_t i4_power ( size_t i, size_t j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_POWER returns the value of I^J.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, size_t I, J, the base and the power.  J should be nonnegative.
//
//    Output, size_t I4_POWER, the value of I^J.
//
{
  size_t k;
  size_t value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J negative.\n";
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J = 0.\n";
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}

//****************************************************************************80

void d2xy ( size_t m, size_t d, size_t &x, size_t &y )

//****************************************************************************80
//
//  Purpose:
//
//    D2XY converts a 1D Hilbert coordinate to a 2D Cartesian coordinate.
//
//  Modified:
//
//    24 December 2015
//
//  Parameters:
//
//    Input, size_t M, the index of the Hilbert curve.
//    The number of cells is N=2^M.
//    0 < M.
//
//    Input, size_t D, the Hilbert coordinate of the cell.
//    0 <= D < N * N.
//
//    Output, size_t &X, &Y, the Cartesian coordinates of the cell.
//    0 <= X, Y < N.
//
{
  size_t n;
  size_t rx;
  size_t ry;
  size_t s;
  size_t t = d;

  n = i4_power ( 2, m );

  x = 0;
  y = 0;
  for ( s = 1; s < n; s = s * 2 )
  {
    rx = ( 1 & ( t / 2 ) );
    ry = ( 1 & ( t ^ rx ) );
    rot ( s, x, y, rx, ry );
    x = x + s * rx;
    y = y + s * ry;
    t = t / 4;
  }
  return;
}

size_t xy2d ( size_t m, size_t x, size_t y )

//****************************************************************************80
//
//  Purpose:
//
//    XY2D converts a 2D Cartesian coordinate to a 1D Hilbert coordinate.
//
//  Discussion:
//
//    It is assumed that a square has been divided into an NxN array of cells,
//    where N is a power of 2.
//
//    Cell (0,0) is in the lower left corner, and (N-1,N-1) in the upper
//    right corner.
//
//  Modified:
//
//    24 December 2015
//
//  Parameters:
//
//    Input, size_t M, the index of the Hilbert curve.
//    The number of cells is N=2^M.
//    0 < M.
//
//    Input, size_t X, Y, the Cartesian coordinates of a cell.
//    0 <= X, Y < N.
//
//    Output, size_t XY2D, the Hilbert coordinate of the cell.
//    0 <= D < N * N.
//
{
  size_t d = 0;
  size_t n;
  size_t rx;
  size_t ry;
  size_t s;

  n = i4_power ( 2, m );

  for ( s = n / 2; s > 0; s = s / 2 )
  {
    rx = ( x & s ) > 0;
    ry = ( y & s ) > 0;
    d = d + s * s * ( ( 3 * rx ) ^ ry );
    rot ( s, x, y, rx, ry );
  }
  return d;
}


#endif /* SRC_INDEX_HILBERT_CURVE_H_ */
