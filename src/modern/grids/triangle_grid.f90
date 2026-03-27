!> triangle_grid -- Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module triangle_grid_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: triangle_grid

contains

  pure subroutine triangle_grid ( n, t, tg ) &
        bind(C, name="triangle_grid")

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
  !    Input, integer(ip) N, the number of subintervals.
  !
  !    Input, real(dp) T(2,3), the coordinates of the points
  !    defining the triangle.
  !
  !    Output, real(dp) TG(2,((N+1)*(N+2))/2), the coordinates
  !    of the points in the triangle.
  !

    integer(ip), intent(in), value :: n
    integer(ip) :: i
    real(dp) :: ir
    integer(ip) :: j
    real(dp) :: jr
    integer(ip) :: k
    real(dp) :: kr
    real(dp) :: nr
    integer(ip) :: p
    real(dp), intent(in) :: t(2,3)
    real(dp), intent(out) :: tg(2,((n+1)*(n+2))/2)

    p = 0
    nr = real ( n, dp)

    do i = 0, n
      ir = real ( i, dp)
      do j = 0, n - i
        jr = real ( j, dp)
        k = n - i - j
        kr = real ( k, dp)
        p = p + 1
        tg(1:2,p) = ( ir * t(1:2,1) + jr * t(1:2,2) + kr * t(1:2,3) ) / nr
      end do
    end do
  end subroutine triangle_grid

end module triangle_grid_mod
