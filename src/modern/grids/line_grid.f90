!> line_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

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
!    Input, integer N, the number of points.
!
!    Input, double precision A, B, the endpoints.
!
!    Input, integer C, the grid centering.
!    1 <= C <= 5.
!
!    Output, double precision X(N), the points.
!
  implicit none

  integer n

  double precision a
  double precision b
  integer c
  integer j
  double precision x(n)

  do j = 1, n

    if ( c == 1 ) then

      if ( n == 1 ) then
        x(j) = 0.5D+00 * ( a + b )
      else
        x(j) = (   real ( n - j) * a   &
                 + real (     j - 1) * b ) & 
                 / real ( n     - 1)
      end if
    else if ( c == 2 ) then
      x(j) = (   real ( n - j + 1) * a   &
               + real (     j) * b ) & 
               / real ( n     + 1)
    else if ( c == 3 ) then
      x(j) = (   real ( n - j + 1) * a   &
               + real (     j - 1) * b ) & 
               / real ( n)
    else if ( c == 4 ) then
      x(j) = (   real ( n - j) * a   &
               + real (     j) * b ) & 
               / real ( n)
    else if ( c == 5 ) then
      x(j) = (   real ( 2 * n - 2 * j + 1) * a   &
               + real (         2 * j - 1) * b ) & 
               / real ( 2 * n)
    end if

  end do
end
