!> triangle_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

subroutine triangle_grid ( n, t, tg )  

!*****************************************************************************80
!
!! TRIANGLE_GRID computes points on a triangular grid.
!
!  Discussion:
!
!    The grid is defined by specifying the coordinates of an enclosing
!    triangle T, and the number of subintervals each side of the triangle
!    should be divided into.
!
!    Choosing N = 10, for instance, breaks each side into 10 subintervals,
!    and produces a grid of ((10+1)*(10+2))/2 = 66 points.
!
!              X
!             9 X
!            8 9 X
!           7 8 9 X
!          6 7 8 9 X
!         5 6 7 8 9 X
!        4 5 6 7 8 9 X
!       3 4 5 6 7 8 9 X
!      2 3 4 5 6 7 8 9 X
!     1 2 3 4 5 6 7 8 9 X
!    0 1 2 3 4 5 6 7 8 9 X
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of subintervals.
!
!    Input, double precision T(2,3), the coordinates of the points
!    defining the triangle.
!
!    Output, double precision TG(2,((N+1)*(N+2))/2), the coordinates
!    of the points in the triangle.
!
  implicit none

  integer n

  integer i
  double precision ir
  integer j
  double precision jr
  integer k
  double precision kr
  double precision nr
  integer p
  double precision t(2,3)
  double precision tg(2,((n+1)*(n+2))/2)

  p = 0
  nr = real ( n)

  do i = 0, n
    ir = real ( i)
    do j = 0, n - i
      jr = real ( j)
      k = n - i - j
      kr = real ( k)
      p = p + 1
      tg(1:2,p) = ( ir * t(1:2,1) + jr * t(1:2,2) + kr * t(1:2,3) ) / nr
    end do
  end do
end
