!> cube_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

subroutine cube_grid ( n, ns, a, b, c, x )

!*****************************************************************************80
!
!! CUBE_GRID: grid points over the interior of a cube in 3D.
!
!  Discussion:
!
!    In 3D, a logically rectangular grid is to be created.
!    In the I-th dimension, the grid will use S(I) points.
!    The total number of grid points is 
!      N = product ( 1 <= I <= 3 ) S(I)
!
!    Over the interval [A(i),B(i)], we have 5 choices for grid centering:
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
!    N = product ( 1 <= I <= 3 ) NS(I).
!
!    Input, integer NS(3), the number of points along 
!    each dimension.
!
!    Input, double precision A(3), B(3), the endpoints for each dimension.
!
!    Input, integer C(3), the grid centering for each dimension.
!    1 <= C(*) <= 5.
!
!    Output, double precision X(3,N) = X(3,S(1)*S(2)*S(3)), the points.
!
  implicit none

  integer , parameter :: m = 3
  integer n

  double precision a(m)
  double precision b(m)
  integer c(m)
  integer i
  integer j
  integer ns(m)
  integer s
  double precision x(m,n)
  double precision , allocatable :: xs(:)
!
!  Create the 1D grids in each dimension.
!
  do i = 1, m

    s = ns(i)

    allocate ( xs(1:s) )

    do j = 1, s

      if ( c(i) == 1 ) then

        if ( s == 1 ) then
          xs(j) = 0.5D+00 * ( a(i) + b(i) )
        else
          xs(j) = (   real ( s - j) * a(i)   &
                    + real (     j - 1) * b(i) ) & 
                    / real ( s     - 1)
        end if
      else if ( c(i) == 2 ) then
        xs(j) = (   real ( s - j + 1) * a(i)   &
                  + real (     j) * b(i) ) & 
                  / real ( s     + 1)
      else if ( c(i) == 3 ) then
        xs(j) = (   real ( s - j + 1) * a(i)   &
                  + real (     j - 1) * b(i) ) & 
                  / real ( s)
      else if ( c(i) == 4 ) then
        xs(j) = (   real ( s - j) * a(i)   &
                  + real (     j) * b(i) ) & 
                  / real ( s)
      else if ( c(i) == 5 ) then
        xs(j) = (   real ( 2 * s - 2 * j + 1) * a(i)   &
                  + real (         2 * j - 1) * b(i) ) & 
                  / real ( 2 * s)
      end if

    end do

    call r8vec_direct_product ( i, s, xs, m, n, x )

    deallocate ( xs )

  end do
end
