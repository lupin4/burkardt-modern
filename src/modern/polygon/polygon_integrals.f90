!> polygon_integrals — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

subroutine moment ( n, x, y, p, q, nu_pq )

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
!    Input, integer N, the number of vertices of the polygon.
!
!    Input, double precision X(N), Y(N), the vertex coordinates.
!
!    Input, integer P, Q, the indices of the moment.
!
!    Output, double precision NU_PQ, the unnormalized moment Nu(P,Q).
!
  implicit none

  integer n

  integer i
  integer k
  integer l
  double precision nu_pq
  integer p
  integer q
  double precision r8_choose
  double precision s_pq
  double precision x(n)
  double precision xi
  double precision xj
  double precision y(n)
  double precision yi
  double precision yj

  nu_pq = 0.0D+00

  xj = x(n)
  yj = y(n)

  do i = 1, n

    xi = x(i)
    yi = y(i)

    s_pq = 0.0D+00
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

  nu_pq = nu_pq / real ( p + q + 2) &
    / real ( p + q + 1) &
    / r8_choose ( p + q, p )
end

subroutine moment_central ( n, x, y, p, q, mu_pq )

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
!    Input, integer N, the number of vertices of the polygon.
!
!    Input, double precision X(N), Y(N), the vertex coordinates.
!
!    Input, integer P, Q, the indices of the moment.
!
!    Output, double precision MU_PQ, the unnormalized moment Mu(P,Q).
!
  implicit none

  integer n

  double precision alpha_01
  double precision alpha_10
  double precision alpha_ij
  integer i
  integer j
  double precision mu_pq
  integer p
  integer q
  double precision r8_choose
  double precision r8_mop
  double precision x(n)
  double precision y(n)

  call moment_normalized ( n, x, y, 1, 0, alpha_10 )
  call moment_normalized ( n, x, y, 0, 1, alpha_01 )

  mu_pq = 0.0D+00

  do i = 0, p
    do j = 0, q

      call moment_normalized ( n, x, y, i, j, alpha_ij )

      mu_pq = mu_pq + r8_mop ( p + q - i - j ) &
        * r8_choose ( p, i ) * r8_choose ( q, j ) &
        * alpha_10 ** ( p - i ) * alpha_01 ** ( q - j ) * alpha_ij

    end do
  end do
end

subroutine moment_normalized ( n, x, y, p, q, alpha_pq )

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
!    Input, integer N, the number of vertices of the polygon.
!
!    Input, double precision X(N), Y(N), the vertex coordinates.
!
!    Input, integer P, Q, the indices of the moment.
!
!    Output, double precision ALPHA_PQ, the normalized moment Alpha(P,Q).
!
  implicit none

  integer n

  double precision alpha_pq
  double precision nu_00
  double precision nu_pq
  integer p
  integer q
  double precision x(n)
  double precision y(n)

  call moment ( n, x, y, p, q, nu_pq )
  call moment ( n, x, y, 0, 0, nu_00 )

  alpha_pq = nu_pq / nu_00
end

function r8_choose ( n, k )

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
!    Input, integer N, K, are the values of N and K.
!
!    Output, double precision R8_CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer i
  integer k
  integer mn
  integer mx
  integer n
  double precision r8_choose
  double precision value

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    value = 0.0D+00

  else if ( mn == 0 ) then

    value = 1.0D+00

  else

    mx = max ( k, n - k )
    value = real ( mx + 1)

    do i = 2, mn
      value = ( value * real ( mx + i) ) / real ( i)
    end do

  end if

  r8_choose = value
end

function r8_mop ( i )

!*****************************************************************************80
!
!! R8_MOP returns the I-th power of -1 as an R8.
!
!  Discussion:
!
!    An R8 is a double precision value.
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
!    Input, integer I, the power of -1.
!
!    Output, double precision R8_MOP, the I-th power of -1.
!
  implicit none

  integer i
  double precision r8_mop

  if ( mod ( i, 2 ) == 0 ) then
    r8_mop = + 1.0D+00
  else
    r8_mop = - 1.0D+00
  end if
end
