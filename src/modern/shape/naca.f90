!> naca — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

subroutine naca4_cambered ( m, p, t, c, n, xc, xu, yu, xl, yl )

!*****************************************************************************80
!
!! NACA4_CAMBERED: (xu,yu), (xl,yl) for a NACA cambered 4-digit airfoil.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eastman Jacobs, Kenneth Ward, Robert Pinkerton,
!    "The characteristics of 78 related airfoil sections from tests in
!    the variable-density wind tunnel",
!    NACA Report 460, 1933.
!
!  Parameters:
!
!    Input, double precision M, the maximum camber.
!    0.0 < M.
!
!    Input, double precision P, the location of maximum camber.
!    0.0 < P < 1.0
!
!    Input, double precision T, the maximum relative thickness.
!    0.0 < T <= 1.0
!
!    Input, double precision C, the chord length.
!    0.0 < C.
!
!    Input, integer N, the number of sample points.
!
!    Input, double precision XC(N), points along the chord length.  
!    0.0 <= XC(*) <= C.
!
!    Output, double precision XU(N), YU(N), XL(N), YL(N), for each value of 
!    XC, measured along the camber line, the corresponding values (XU,YU) 
!    on the upper airfoil surface and (XL,YL) on the lower airfoil surface.
!
  implicit none

  integer n

  double precision c
  double precision divisor
  double precision dycdx
  integer i
  double precision m
  double precision p
  double precision t
  double precision theta
  double precision xc(n)
  double precision xl(n)
  double precision xu(n)
  double precision yc
  double precision yl(n)
  double precision yt
  double precision yu(n)

  do i = 1, n

    if ( 0.0D+00 <= xc(i) / c .and. xc(i) / c <= p ) then
      divisor = p ** 2
    else if ( p <= xc(i) / c .and. xc(i) / c <= 1.0D+00 ) then
      divisor = ( 1.0D+00 - p ) ** 2
    else
      divisor = 1.0D+00
    end if

    dycdx = 2.0D+00 * m * ( p - xc(i) / c ) / divisor

    theta = atan ( dycdx )
   
    yt = 5.0D+00 * t * c * ( &
       0.2969D+00 * sqrt ( xc(i) / c ) &
       + (((( &
         - 0.1015D+00 ) * ( xc(i) / c ) &
         + 0.2843D+00 ) * ( xc(i) / c ) &
         - 0.3516D+00 ) * ( xc(i) / c ) &
         - 0.1260D+00 ) * ( xc(i) / c ) )

    if ( 0.0D+00 <= xc(i) / c .and. xc(i) / c <= p ) then
      yc = m * xc(i) * ( 2.0D+00 * p - xc(i) / c ) / p ** 2
    else if ( p <= xc(i) / c .and. xc(i) / c <= 1.0D+00 ) then
      yc = m * ( xc(i) - c ) * ( 2.0D+00 * p - xc(i) / c - 1.0D+00 ) &
        / ( 1.0D+00 - p ) ** 2
    else
      yc = 0.0D+00
    end if

    xu(i) = xc(i) - yt * sin ( theta )
    yu(i) = yc + yt * cos ( theta )
    xl(i) = xc(i) + yt * sin ( theta )
    yl(i) = yc - yt * cos ( theta )

  end do
end

subroutine naca4_symmetric ( t, c, n, x, y )

!*****************************************************************************80
!
!! NACA4_SYMMETRIC evaluates y(x) for a NACA symmetric 4-digit airfoil.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eastman Jacobs, Kenneth Ward, Robert Pinkerton,
!    "The characteristics of 78 related airfoil sections from tests in
!    the variable-density wind tunnel",
!    NACA Report 460, 1933.
!
!  Parameters:
!
!    Input, double precision T, the maximum relative thickness.
!
!    Input, double precision C, the chord length.
!
!    Input, integer N, the number of sample points.
!
!    Input, double precision X(N), points along the chord length.  
!    0.0 <= X(*) <= C.
!
!    Output, double precision Y(N), for each value of X, the corresponding
!    value of Y so that (X,Y) is on the upper wing surface, and (X,-Y) is on the
!    lower wing surface.
!
  implicit none

  integer n

  double precision c
  double precision t
  double precision x(n)
  double precision y(n)

  y(1:n) = 5.0D+00 * t * c * ( &
    0.2969D+00 * sqrt ( x(1:n) / c ) &
    + (((( &
      - 0.1015D+00 ) * ( x(1:n) / c ) &
      + 0.2843D+00 ) * ( x(1:n) / c ) &
      - 0.3516D+00 ) * ( x(1:n) / c ) &
      - 0.1260D+00 ) * ( x(1:n) / c ) )
end

subroutine r8vec_linspace ( n, a, b, x )

!*****************************************************************************80
!
!! R8VEC_LINSPACE creates a vector of linearly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
!
!    In other words, the interval is divided into N-1 even subintervals,
!    and the endpoints of intervals are used as the points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, double precision A, B, the first and last entries.
!
!    Output, double precision X(N), a vector of linearly spaced data.
!
  implicit none

  integer n

  double precision a
  double precision b
  integer i
  double precision x(n)

  if ( n == 1 ) then

    x(1) = ( a + b ) / 2.0D+00

  else

    do i = 1, n
      x(i) = ( real ( n - i) * a   &
             + real (     i - 1) * b ) &
             / real ( n     - 1)
    end do

  end if
end

function r8vec_max ( n, a )

!*****************************************************************************80
!
!! R8VEC_MAX returns the maximum value in an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, double precision A(N), the array.
!
!    Output, double precision R8VEC_MAX, the value of the largest entry.
!
  implicit none

  integer n

  double precision a(n)
  double precision r8vec_max
  double precision value

  value = maxval ( a(1:n) )

  r8vec_max = value
end

function r8vec_min ( n, a )

!*****************************************************************************80
!
!! R8VEC_MIN returns the minimum value of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, double precision A(N), the array.
!
!    Output, double precision R8VEC_MIN, the value of the smallest entry.
!
  implicit none

  integer n

  double precision a(n)
  double precision r8vec_min
  double precision value

  value = minval ( a(1:n) )

  r8vec_min = value
end
