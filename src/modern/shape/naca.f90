!> naca — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module naca_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  implicit none
  private

  public :: naca4_cambered, naca4_symmetric, r8vec_linspace, r8vec_max, r8vec_min

contains

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
  !    Input, real(real64) M, the maximum camber.
  !    0.0 < M.
  !
  !    Input, real(real64) P, the location of maximum camber.
  !    0.0 < P < 1.0
  !
  !    Input, real(real64) T, the maximum relative thickness.
  !    0.0 < T <= 1.0
  !
  !    Input, real(real64) C, the chord length.
  !    0.0 < C.
  !
  !    Input, integer(int32) N, the number of sample points.
  !
  !    Input, real(real64) XC(N), points along the chord length.  
  !    0.0 <= XC(*) <= C.
  !
  !    Output, real(real64) XU(N), YU(N), XL(N), YL(N), for each value of 
  !    XC, measured along the camber line, the corresponding values (XU,YU) 
  !    on the upper airfoil surface and (XL,YL) on the lower airfoil surface.
  !

    integer(int32) n

    real(real64) c
    real(real64) divisor
    real(real64) dycdx
    integer(int32) i
    real(real64) m
    real(real64) p
    real(real64) t
    real(real64) theta
    real(real64) xc(n)
    real(real64) xl(n)
    real(real64) xu(n)
    real(real64) yc
    real(real64) yl(n)
    real(real64) yt
    real(real64) yu(n)

    do i = 1, n

      if ( 0.0e+00_real64 <= xc(i) / c .and. xc(i) / c <= p ) then
        divisor = p ** 2
      else if ( p <= xc(i) / c .and. xc(i) / c <= 1.0e+00_real64 ) then
        divisor = ( 1.0e+00_real64 - p ) ** 2
      else
        divisor = 1.0e+00_real64
      end if

      dycdx = 2.0e+00_real64 * m * ( p - xc(i) / c ) / divisor

      theta = atan ( dycdx )

      yt = 5.0e+00_real64 * t * c * ( &
         0.2969e+00_real64 * sqrt ( xc(i) / c ) &
         + (((( &
           - 0.1015e+00_real64 ) * ( xc(i) / c ) &
           + 0.2843e+00_real64 ) * ( xc(i) / c ) &
           - 0.3516e+00_real64 ) * ( xc(i) / c ) &
           - 0.1260e+00_real64 ) * ( xc(i) / c ) )

      if ( 0.0e+00_real64 <= xc(i) / c .and. xc(i) / c <= p ) then
        yc = m * xc(i) * ( 2.0e+00_real64 * p - xc(i) / c ) / p ** 2
      else if ( p <= xc(i) / c .and. xc(i) / c <= 1.0e+00_real64 ) then
        yc = m * ( xc(i) - c ) * ( 2.0e+00_real64 * p - xc(i) / c - 1.0e+00_real64 ) &
          / ( 1.0e+00_real64 - p ) ** 2
      else
        yc = 0.0e+00_real64
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
  !    Input, real(real64) T, the maximum relative thickness.
  !
  !    Input, real(real64) C, the chord length.
  !
  !    Input, integer(int32) N, the number of sample points.
  !
  !    Input, real(real64) X(N), points along the chord length.  
  !    0.0 <= X(*) <= C.
  !
  !    Output, real(real64) Y(N), for each value of X, the corresponding
  !    value of Y so that (X,Y) is on the upper wing surface, and (X,-Y) is on the
  !    lower wing surface.
  !

    integer(int32) n

    real(real64) c
    real(real64) t
    real(real64) x(n)
    real(real64) y(n)

    y(1:n) = 5.0e+00_real64 * t * c * ( &
      0.2969e+00_real64 * sqrt ( x(1:n) / c ) &
      + (((( &
        - 0.1015e+00_real64 ) * ( x(1:n) / c ) &
        + 0.2843e+00_real64 ) * ( x(1:n) / c ) &
        - 0.3516e+00_real64 ) * ( x(1:n) / c ) &
        - 0.1260e+00_real64 ) * ( x(1:n) / c ) )
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
  !    Input, integer(int32) N, the number of entries in the vector.
  !
  !    Input, real(real64) A, B, the first and last entries.
  !
  !    Output, real(real64) X(N), a vector of linearly spaced data.
  !

    integer(int32) n

    real(real64) a
    real(real64) b
    integer(int32) i
    real(real64) x(n)

    if ( n == 1 ) then

      x(1) = ( a + b ) / 2.0e+00_real64

    else

      do i = 1, n
        x(i) = ( real ( n - i, real64) * a   &
               + real (     i - 1, real64) * b ) &
               / real ( n     - 1, real64)
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
  !    Input, integer(int32) N, the number of entries in the array.
  !
  !    Input, real(real64) A(N), the array.
  !
  !    Output, real(real64) R8VEC_MAX, the value of the largest entry.
  !

    integer(int32) n

    real(real64) a(n)
    real(real64) r8vec_max
    real(real64) value

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
  !    Input, integer(int32) N, the number of entries in the array.
  !
  !    Input, real(real64) A(N), the array.
  !
  !    Output, real(real64) R8VEC_MIN, the value of the smallest entry.
  !

    integer(int32) n

    real(real64) a(n)
    real(real64) r8vec_min
    real(real64) value

    value = minval ( a(1:n) )

    r8vec_min = value
  end

end module naca_mod
