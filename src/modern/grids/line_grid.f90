!> line_grid -- Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module line_grid_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: line_grid

contains

  pure subroutine line_grid ( n, a, b, c, x ) &
        bind(C, name="line_grid")

  !*****************************************************************************80
  !
  !! LINEE_GRID: grid points over the interior of a line segment in 1D.
  !
  !  Discussion:
  !
  !    In 1D, a grid is to be created, using N points.
  !
  !    Over the interval [A,B], we have 5 choices for grid centering:
  !      1: 0,   1/3, 2/3, 1
  !      2: 1/5, 2/5, 3/5, 4/5
  !      3: 0,   1/4, 2/4, 3/4
  !      4: 1/4, 2/4, 3/4, 1
  !      5: 1/8, 3/8, 5/8, 7/8
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    31 August 2014
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of points.
  !
  !    Input, real(dp) A, B, the endpoints.
  !
  !    Input, integer(ip) C, the grid centering.
  !    1 <= C <= 5.
  !
  !    Output, real(dp) X(N), the points.
  !

    integer(ip), intent(in), value :: n
    real(dp), intent(in), value :: a
    real(dp), intent(in), value :: b
    integer(ip), intent(in), value :: c
    integer(ip) :: j
    real(dp), intent(out) :: x(n)

    do j = 1, n

      if ( c == 1 ) then

        if ( n == 1 ) then
          x(j) = 0.5_dp * ( a + b )
        else
          x(j) = (   real ( n - j, dp) * a   &
                   + real (     j - 1, dp) * b ) &
                   / real ( n     - 1, dp)
        end if
      else if ( c == 2 ) then
        x(j) = (   real ( n - j + 1, dp) * a   &
                 + real (     j, dp) * b ) &
                 / real ( n     + 1, dp)
      else if ( c == 3 ) then
        x(j) = (   real ( n - j + 1, dp) * a   &
                 + real (     j - 1, dp) * b ) &
                 / real ( n, dp)
      else if ( c == 4 ) then
        x(j) = (   real ( n - j, dp) * a   &
                 + real (     j, dp) * b ) &
                 / real ( n, dp)
      else if ( c == 5 ) then
        x(j) = (   real ( 2 * n - 2 * j + 1, dp) * a   &
                 + real (         2 * j - 1, dp) * b ) &
                 / real ( 2 * n, dp)
      end if

    end do
  end subroutine line_grid

end module line_grid_mod
