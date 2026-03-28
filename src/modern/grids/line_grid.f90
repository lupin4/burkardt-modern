!> line_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module line_grid_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  implicit none
  private

  public :: line_grid

contains

  subroutine line_grid ( n, a, b, c, x )

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
  !    Input, integer(int32) N, the number of points.
  !
  !    Input, real(real64) A, B, the endpoints.
  !
  !    Input, integer(int32) C, the grid centering.
  !    1 <= C <= 5.
  !
  !    Output, real(real64) X(N), the points.
  !

    integer(int32) n

    real(real64) a
    real(real64) b
    integer(int32) c
    integer(int32) j
    real(real64) x(n)

    do j = 1, n

      if ( c == 1 ) then

        if ( n == 1 ) then
          x(j) = 0.5e+00_real64 * ( a + b )
        else
          x(j) = (   real ( n - j, real64) * a   &
                   + real (     j - 1, real64) * b ) & 
                   / real ( n     - 1, real64)
        end if
      else if ( c == 2 ) then
        x(j) = (   real ( n - j + 1, real64) * a   &
                 + real (     j, real64) * b ) & 
                 / real ( n     + 1, real64)
      else if ( c == 3 ) then
        x(j) = (   real ( n - j + 1, real64) * a   &
                 + real (     j - 1, real64) * b ) & 
                 / real ( n, real64)
      else if ( c == 4 ) then
        x(j) = (   real ( n - j, real64) * a   &
                 + real (     j, real64) * b ) & 
                 / real ( n, real64)
      else if ( c == 5 ) then
        x(j) = (   real ( 2 * n - 2 * j + 1, real64) * a   &
                 + real (         2 * j - 1, real64) * b ) & 
                 / real ( 2 * n, real64)
      end if

    end do
  end

end module line_grid_mod
