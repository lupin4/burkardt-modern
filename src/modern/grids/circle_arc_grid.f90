!> circle_arc_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

subroutine circle_arc_grid ( r, c, a, n, xy )

!*****************************************************************************80
!
!! CIRCLE_ARC_GRID computes grid points along a circular arc.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision R, the radius of the circle.
!
!    Input, double precision C(2), the coordinates of the center.
!
!    Input, double precision A(2), the angle of the first and last
!    points, in DEGREES.
!
!    Input, integer N, the number of points.
!
!    Output, double precision XY(2,N), the grid points.
!
  implicit none

  integer n

  double precision a(2)
  double precision aj
  double precision c(2)
  integer j
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision r
  double precision xy(2,n)

  do j = 1, n

    aj = ( real ( n - j) * a(1)   &
         + real (     j - 1) * a(2) ) &
         / real ( n     - 1)

    xy(1,j) = c(1) + r * cos ( aj * pi / 180.0D+00 )
    xy(2,j) = c(2) + r * sin ( aj * pi / 180.0D+00 )

  end do
end
