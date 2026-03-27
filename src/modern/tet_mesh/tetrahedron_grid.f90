!> tetrahedron_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module tetrahedron_grid_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: tetrahedron_grid, tetrahedron_grid_count

contains

  pure subroutine tetrahedron_grid ( n, t, ng, tg ) &
        bind(C, name="tetrahedron_grid")

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
  !    Input, integer(ip) N, the number of subintervals.
  !
  !    Input, real(dp) T(3,4), the vertices of the tetrahedron.
  !
  !    Input, integer(ip) NG, the number of grid points.
  !
  !    Output, real(dp) TG(3,NG), the tetrahedron grid points.
  !

    integer(ip), intent(in), value :: ng

    integer(ip) :: i
    integer(ip) :: j
    integer(ip) :: k
    integer(ip) :: l
    integer(ip), intent(in), value :: n
    integer(ip) :: p
    real(dp), intent(in) :: t(3,4)
    real(dp), intent(out) :: tg(3,ng)

    p = 0

    do i = 0, n
      do j = 0, n - i
        do k = 0, n - i - j
          l = n - i - j - k
          p = p + 1
          tg(1:3,p) = ( real ( i, dp) * t(1:3,1) &
                      + real ( j, dp) * t(1:3,2) &
                      + real ( k, dp) * t(1:3,3) &
                      + real ( l, dp) * t(1:3,4) ) / real ( n, dp)
        end do
      end do
    end do
  end subroutine tetrahedron_grid

  pure subroutine tetrahedron_grid_count ( n, ng ) &
        bind(C, name="tetrahedron_grid_count")

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
  !    Input, integer(ip) N, the number of subintervals.
  !
  !    Output, integer(ip) NG, the number of grid points.
  !

    integer(ip), intent(in), value :: n
    integer(ip), intent(out) :: ng

    ng = ( ( n + 1 ) * ( n + 2 ) * ( n + 3 ) ) / 6
  end subroutine tetrahedron_grid_count

end module tetrahedron_grid_mod
