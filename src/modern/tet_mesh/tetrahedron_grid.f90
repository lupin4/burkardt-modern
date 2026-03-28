!> tetrahedron_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module tetrahedron_grid_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  implicit none
  private

  public :: tetrahedron_grid, tetrahedron_grid_count

contains

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
  !    Input, integer(int32) N, the number of subintervals.
  !
  !    Input, real(real64) T(3,4), the vertices of the tetrahedron.
  !
  !    Input, integer(int32) NG, the number of grid points.
  !
  !    Output, real(real64) TG(3,NG), the tetrahedron grid points.
  !

    integer(int32) ng

    integer(int32) i
    integer(int32) j
    integer(int32) k
    integer(int32) l
    integer(int32) n
    integer(int32) p
    real(real64) t(3,4)
    real(real64) tg(3,ng)

    p = 0

    do i = 0, n
      do j = 0, n - i
        do k = 0, n - i - j
          l = n - i - j - k
          p = p + 1
          tg(1:3,p) = ( real ( i, real64) * t(1:3,1) &
                      + real ( j, real64) * t(1:3,2) &
                      + real ( k, real64) * t(1:3,3) &
                      + real ( l, real64) * t(1:3,4) ) / real ( n, real64)
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
  !    Input, integer(int32) N, the number of subintervals.
  !
  !    Output, integer(int32) NG, the number of grid points.
  !

    integer(int32) n
    integer(int32) ng

    ng = ( ( n + 1 ) * ( n + 2 ) * ( n + 3 ) ) / 6
  end

end module tetrahedron_grid_mod
