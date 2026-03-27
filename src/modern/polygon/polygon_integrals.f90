!> polygon_integrals — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module polygon_integrals_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: moment, moment_central, moment_normalized, r8_choose, r8_mop

contains

  subroutine moment ( n, x, y, p, q, nu_pq ) &
        bind(C, name="moment")

  !*****************************************************************************80
  !
  !! MOMENT computes an unnormalized moment of a polygon.
  !
  !  Discussion:
  !
  !    Nu(P,Q) = Integral ( x, y in polygon ) x^p y^q dx dy
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    03 October 2012
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Carsten Steger,
  !    On the calculation of arbitrary moments of polygons,
  !    Technical Report FGBV-96-05,
  !    Forschungsgruppe Bildverstehen, Informatik IX,
  !    Technische Universitaet Muenchen, October 1996.
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of vertices of the polygon.
  !
  !    Input, real(dp) X(N), Y(N), the vertex coordinates.
  !
  !    Input, integer(ip) P, Q, the indices of the moment.
  !
  !    Output, real(dp) NU_PQ, the unnormalized moment Nu(P,Q).
  !

    integer(ip), intent(in), value :: n
    real(dp), intent(in) :: x(n)
    real(dp), intent(in) :: y(n)
    integer(ip), intent(in), value :: p
    integer(ip), intent(in), value :: q
    real(dp), intent(out) :: nu_pq

    integer(ip) :: i
    integer(ip) :: k
    integer(ip) :: l
    real(dp) :: r8_choose
    real(dp) :: s_pq
    real(dp) :: xi
    real(dp) :: xj
    real(dp) :: yi
    real(dp) :: yj

    nu_pq = 0.0_dp

    xj = x(n)
    yj = y(n)

    do i = 1, n

      xi = x(i)
      yi = y(i)

      s_pq = 0.0_dp
      do k = 0, p
        do l = 0, q
          s_pq = s_pq &
            + r8_choose ( k + l, l ) * r8_choose ( p + q - k - l, q - l ) &
            * xi ** k * xj ** ( p - k ) &
            * yi ** l * yj ** ( q - l )
        end do
      end do

      nu_pq = nu_pq + ( xj * yi - xi * yj ) * s_pq

      xj = xi
      yj = yi

    end do

    nu_pq = nu_pq / real ( p + q + 2, dp) &
      / real ( p + q + 1, dp) &
      / r8_choose ( p + q, p )
  end subroutine moment

  subroutine moment_central ( n, x, y, p, q, mu_pq ) &
        bind(C, name="moment_central")

  !*****************************************************************************80
  !
  !! MOMENT_CENTRAL computes central moments of a polygon.
  !
  !  Discussion:
  !
  !    The central moment Mu(P,Q) is defined by
  !
  !      Mu(P,Q) = Integral ( polygon ) (x-Alpha(1,0))^p (y-Alpha(0,1))^q dx dy
  !              / Area ( polygon )
  !
  !    where
  !
  !      Alpha(1,0) = Integral ( polygon ) x dx dy / Area ( polygon )
  !      Alpha(0,1) = Integral ( polygon ) y dx dy / Area ( polygon )
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    03 October 2012
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Carsten Steger,
  !    On the calculation of arbitrary moments of polygons,
  !    Technical Report FGBV-96-05,
  !    Forschungsgruppe Bildverstehen, Informatik IX,
  !    Technische Universitaet Muenchen, October 1996.
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of vertices of the polygon.
  !
  !    Input, real(dp) X(N), Y(N), the vertex coordinates.
  !
  !    Input, integer(ip) P, Q, the indices of the moment.
  !
  !    Output, real(dp) MU_PQ, the unnormalized moment Mu(P,Q).
  !

    integer(ip), intent(in), value :: n
    real(dp), intent(in) :: x(n)
    real(dp), intent(in) :: y(n)
    integer(ip), intent(in), value :: p
    integer(ip), intent(in), value :: q
    real(dp), intent(out) :: mu_pq

    real(dp) :: alpha_01
    real(dp) :: alpha_10
    real(dp) :: alpha_ij
    integer(ip) :: i
    integer(ip) :: j
    real(dp) :: r8_choose
    real(dp) :: r8_mop

    call moment_normalized ( n, x, y, 1, 0, alpha_10 )
    call moment_normalized ( n, x, y, 0, 1, alpha_01 )

    mu_pq = 0.0_dp

    do i = 0, p
      do j = 0, q

        call moment_normalized ( n, x, y, i, j, alpha_ij )

        mu_pq = mu_pq + r8_mop ( p + q - i - j ) &
          * r8_choose ( p, i ) * r8_choose ( q, j ) &
          * alpha_10 ** ( p - i ) * alpha_01 ** ( q - j ) * alpha_ij

      end do
    end do
  end subroutine moment_central

  subroutine moment_normalized ( n, x, y, p, q, alpha_pq ) &
        bind(C, name="moment_normalized")

  !*****************************************************************************80
  !
  !! MOMENT_NORMALIZED computes a normalized moment of a polygon.
  !
  !  Discussion:
  !
  !    Alpha(P,Q) = Integral ( x, y in polygon ) x^p y^q dx dy / Area ( polygon )
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    03 October 2012
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Carsten Steger,
  !    On the calculation of arbitrary moments of polygons,
  !    Technical Report FGBV-96-05,
  !    Forschungsgruppe Bildverstehen, Informatik IX,
  !    Technische Universitaet Muenchen, October 1996.
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of vertices of the polygon.
  !
  !    Input, real(dp) X(N), Y(N), the vertex coordinates.
  !
  !    Input, integer(ip) P, Q, the indices of the moment.
  !
  !    Output, real(dp) ALPHA_PQ, the normalized moment Alpha(P,Q).
  !

    integer(ip), intent(in), value :: n
    real(dp), intent(in) :: x(n)
    real(dp), intent(in) :: y(n)
    integer(ip), intent(in), value :: p
    integer(ip), intent(in), value :: q
    real(dp), intent(out) :: alpha_pq

    real(dp) :: nu_00
    real(dp) :: nu_pq

    call moment ( n, x, y, p, q, nu_pq )
    call moment ( n, x, y, 0, 0, nu_00 )

    alpha_pq = nu_pq / nu_00
  end subroutine moment_normalized

  pure function r8_choose ( n, k ) &
        bind(C, name="r8_choose")

  !*****************************************************************************80
  !
  !! R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
  !
  !  Discussion:
  !
  !    The value is calculated in such a way as to avoid overflow and
  !    roundoff.  The calculation is done in R8 arithmetic.
  !
  !    The formula used is:
  !
  !      C(N,K) = N! / ( K! * (N-K)! )
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    24 March 2008
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    ML Wolfson, HV Wright,
  !    Algorithm 160:
  !    Combinatorial of M Things Taken N at a Time,
  !    Communications of the ACM,
  !    Volume 6, Number 4, April 1963, page 161.
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, K, are the values of N and K.
  !
  !    Output, real(dp) R8_CHOOSE, the number of combinations of N
  !    things taken K at a time.
  !

    integer(ip), intent(in), value :: n
    integer(ip), intent(in), value :: k
    real(dp) :: r8_choose

    integer(ip) :: i
    integer(ip) :: mn
    integer(ip) :: mx
    real(dp) :: value

    mn = min ( k, n - k )

    if ( mn < 0 ) then

      value = 0.0_dp

    else if ( mn == 0 ) then

      value = 1.0_dp

    else

      mx = max ( k, n - k )
      value = real ( mx + 1, dp)

      do i = 2, mn
        value = ( value * real ( mx + i, dp) ) / real ( i, dp)
      end do

    end if

    r8_choose = value
  end function r8_choose

  pure function r8_mop ( i ) &
        bind(C, name="r8_mop")

  !*****************************************************************************80
  !
  !! R8_MOP returns the I-th power of -1 as an R8.
  !
  !  Discussion:
  !
  !    An R8 is a real(dp) value.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 November 2007
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) I, the power of -1.
  !
  !    Output, real(dp) R8_MOP, the I-th power of -1.
  !

    integer(ip), intent(in), value :: i
    real(dp) :: r8_mop

    if ( mod ( i, 2 ) == 0 ) then
      r8_mop = + 1.0_dp
    else
      r8_mop = - 1.0_dp
    end if
  end function r8_mop

end module polygon_integrals_mod
