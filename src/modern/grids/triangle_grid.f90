!> triangle_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module triangle_grid_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  implicit none
  private

  public :: triangle_grid

contains

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
  !    Input, integer(int32) N, the number of subintervals.
  !
  !    Input, real(real64) T(2,3), the coordinates of the points
  !    defining the triangle.
  !
  !    Output, real(real64) TG(2,((N+1)*(N+2))/2), the coordinates
  !    of the points in the triangle.
  !

    integer(int32) n

    integer(int32) i
    real(real64) ir
    integer(int32) j
    real(real64) jr
    integer(int32) k
    real(real64) kr
    real(real64) nr
    integer(int32) p
    real(real64) t(2,3)
    real(real64) tg(2,((n+1)*(n+2))/2)

    p = 0
    nr = real ( n, real64)

    do i = 0, n
      ir = real ( i, real64)
      do j = 0, n - i
        jr = real ( j, real64)
        k = n - i - j
        kr = real ( k, real64)
        p = p + 1
        tg(1:2,p) = ( ir * t(1:2,1) + jr * t(1:2,2) + kr * t(1:2,3) ) / nr
      end do
    end do
  end

end module triangle_grid_mod
