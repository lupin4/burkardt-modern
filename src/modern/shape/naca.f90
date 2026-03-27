!> naca -- Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module naca_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: naca4_cambered, naca4_symmetric, r8vec_linspace, r8vec_max, r8vec_min

contains

  pure subroutine naca4_cambered ( m, p, t, c, n, xc, xu, yu, xl, yl ) &
        bind(C, name="naca4_cambered")

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
  !    Input, real(dp) M, the maximum camber.
  !    0.0 < M.
  !
  !    Input, real(dp) P, the location of maximum camber.
  !    0.0 < P < 1.0
  !
  !    Input, real(dp) T, the maximum relative thickness.
  !    0.0 < T <= 1.0
  !
  !    Input, real(dp) C, the chord length.
  !    0.0 < C.
  !
  !    Input, integer(ip) N, the number of sample points.
  !
  !    Input, real(dp) XC(N), points along the chord length.
  !    0.0 <= XC(*) <= C.
  !
  !    Output, real(dp) XU(N), YU(N), XL(N), YL(N), for each value of
  !    XC, measured along the camber line, the corresponding values (XU,YU)
  !    on the upper airfoil surface and (XL,YL) on the lower airfoil surface.
  !

    integer(ip), intent(in), value :: n
    real(dp), intent(in), value :: c
    real(dp), intent(in), value :: m
    real(dp), intent(in), value :: p
    real(dp), intent(in), value :: t
    real(dp), intent(in) :: xc(n)
    real(dp), intent(out) :: xl(n)
    real(dp), intent(out) :: xu(n)
    real(dp), intent(out) :: yl(n)
    real(dp), intent(out) :: yu(n)
    real(dp) :: divisor
    real(dp) :: dycdx
    integer(ip) :: i
    real(dp) :: theta
    real(dp) :: yc
    real(dp) :: yt

    do i = 1, n

      if ( 0.0_dp <= xc(i) / c .and. xc(i) / c <= p ) then
        divisor = p ** 2
      else if ( p <= xc(i) / c .and. xc(i) / c <= 1.0_dp ) then
        divisor = ( 1.0_dp - p ) ** 2
      else
        divisor = 1.0_dp
      end if

      dycdx = 2.0_dp * m * ( p - xc(i) / c ) / divisor

      theta = atan ( dycdx )

      yt = 5.0_dp * t * c * ( &
         0.2969_dp * sqrt ( xc(i) / c ) &
         + (((( &
           - 0.1015_dp ) * ( xc(i) / c ) &
           + 0.2843_dp ) * ( xc(i) / c ) &
           - 0.3516_dp ) * ( xc(i) / c ) &
           - 0.1260_dp ) * ( xc(i) / c ) )

      if ( 0.0_dp <= xc(i) / c .and. xc(i) / c <= p ) then
        yc = m * xc(i) * ( 2.0_dp * p - xc(i) / c ) / p ** 2
      else if ( p <= xc(i) / c .and. xc(i) / c <= 1.0_dp ) then
        yc = m * ( xc(i) - c ) * ( 2.0_dp * p - xc(i) / c - 1.0_dp ) &
          / ( 1.0_dp - p ) ** 2
      else
        yc = 0.0_dp
      end if

      xu(i) = xc(i) - yt * sin ( theta )
      yu(i) = yc + yt * cos ( theta )
      xl(i) = xc(i) + yt * sin ( theta )
      yl(i) = yc - yt * cos ( theta )

    end do
  end subroutine naca4_cambered

  pure subroutine naca4_symmetric ( t, c, n, x, y ) &
        bind(C, name="naca4_symmetric")

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
  !    Input, real(dp) T, the maximum relative thickness.
  !
  !    Input, real(dp) C, the chord length.
  !
  !    Input, integer(ip) N, the number of sample points.
  !
  !    Input, real(dp) X(N), points along the chord length.
  !    0.0 <= X(*) <= C.
  !
  !    Output, real(dp) Y(N), for each value of X, the corresponding
  !    value of Y so that (X,Y) is on the upper wing surface, and (X,-Y) is on the
  !    lower wing surface.
  !

    integer(ip), intent(in), value :: n
    real(dp), intent(in), value :: c
    real(dp), intent(in), value :: t
    real(dp), intent(in) :: x(n)
    real(dp), intent(out) :: y(n)

    y(1:n) = 5.0_dp * t * c * ( &
      0.2969_dp * sqrt ( x(1:n) / c ) &
      + (((( &
        - 0.1015_dp ) * ( x(1:n) / c ) &
        + 0.2843_dp ) * ( x(1:n) / c ) &
        - 0.3516_dp ) * ( x(1:n) / c ) &
        - 0.1260_dp ) * ( x(1:n) / c ) )
  end subroutine naca4_symmetric

  pure subroutine r8vec_linspace ( n, a, b, x ) &
        bind(C, name="r8vec_linspace")

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
  !    Input, integer(ip) N, the number of entries in the vector.
  !
  !    Input, real(dp) A, B, the first and last entries.
  !
  !    Output, real(dp) X(N), a vector of linearly spaced data.
  !

    integer(ip), intent(in), value :: n
    real(dp), intent(in), value :: a
    real(dp), intent(in), value :: b
    real(dp), intent(out) :: x(n)
    integer(ip) :: i

    if ( n == 1 ) then

      x(1) = ( a + b ) / 2.0_dp

    else

      do i = 1, n
        x(i) = ( real ( n - i, dp) * a   &
               + real (     i - 1, dp) * b ) &
               / real ( n     - 1, dp)
      end do

    end if
  end subroutine r8vec_linspace

  pure function r8vec_max ( n, a ) &
        bind(C, name="r8vec_max")

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
  !    Input, integer(ip) N, the number of entries in the array.
  !
  !    Input, real(dp) A(N), the array.
  !
  !    Output, real(dp) R8VEC_MAX, the value of the largest entry.
  !

    integer(ip), intent(in), value :: n
    real(dp), intent(in) :: a(n)
    real(dp) :: r8vec_max
    real(dp) :: value

    value = maxval ( a(1:n) )

    r8vec_max = value
  end function r8vec_max

  pure function r8vec_min ( n, a ) &
        bind(C, name="r8vec_min")

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
  !    Input, integer(ip) N, the number of entries in the array.
  !
  !    Input, real(dp) A(N), the array.
  !
  !    Output, real(dp) R8VEC_MIN, the value of the smallest entry.
  !

    integer(ip), intent(in), value :: n
    real(dp), intent(in) :: a(n)
    real(dp) :: r8vec_min
    real(dp) :: value

    value = minval ( a(1:n) )

    r8vec_min = value
  end function r8vec_min

end module naca_mod
