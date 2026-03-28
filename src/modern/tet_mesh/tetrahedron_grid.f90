!> tetrahedron_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

subroutine tetrahedron_grid ( n, t, ng, tg )

!*****************************************************************************80
!
!! TETRAHEDRON_GRID computes points on a tetrahedral grid.
!
!  Discussion:
!
!    The grid is defined by specifying the coordinates of an enclosing
!    tetrahedron T, and the number of subintervals each edge of the 
!    tetrahedron should be divided into.
!
!    Choosing N = 10, for instance, breaks each side into 10 subintervals,
!    and produces a grid of ((10+1)*(10+2)*(10+3))/6 = 286 points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of subintervals.
!
!    Input, double precision T(3,4), the vertices of the tetrahedron.
!
!    Input, integer NG, the number of grid points.
!
!    Output, double precision TG(3,NG), the tetrahedron grid points.
!
  implicit none

  integer ng

  integer i
  integer j
  integer k
  integer l
  integer n
  integer p
  double precision t(3,4)
  double precision tg(3,ng)

  p = 0

  do i = 0, n
    do j = 0, n - i
      do k = 0, n - i - j
        l = n - i - j - k
        p = p + 1
        tg(1:3,p) = ( real ( i) * t(1:3,1) &
                    + real ( j) * t(1:3,2) &
                    + real ( k) * t(1:3,3) &
                    + real ( l) * t(1:3,4) ) / real ( n)
      end do
    end do
  end do
end

subroutine tetrahedron_grid_count ( n, ng )

!*****************************************************************************80
!
!! TETRAHEDRON_GRID_COUNT counts the grid points inside a tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of subintervals.
!
!    Output, integer NG, the number of grid points.
!
  implicit none

  integer n
  integer ng

  ng = ( ( n + 1 ) * ( n + 2 ) * ( n + 3 ) ) / 6
end
