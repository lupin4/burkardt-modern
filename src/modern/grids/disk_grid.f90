!> disk_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

subroutine disk_grid ( n, r, c, ng, cg )

!*****************************************************************************80
!
!! DISK_GRID computes grid points inside a disk.
!
!  Discussion:
!
!    The grid is defined by specifying the radius and center of the circle,
!    and the number of subintervals N into which the horizontal radius
!    should be divided.  Thus, a value of N = 2 will result in 5 points
!    along that horizontal line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of subintervals.
!
!    Input, double precision R, the radius of the circle.
!
!    Input, double precision C(2), the coordinates of the center of the circle.
!
!    Input, integer NG, the number of grid points, as determined by
!    DISK_GRID_COUNT.
!
!    Output, double precision CG(2,NG), the grid points inside the circle.
!
  implicit none

  integer ng

  double precision c(2)
  double precision cg(2,ng)
  integer i
  integer j
  integer n
  integer p
  double precision r
  double precision x
  double precision y

  p = 0

  do j = 0, n

    i = 0
    x = c(1)
    y = c(2) + r * real ( 2 * j) / real ( 2 * n + 1)
    p = p + 1
    cg(1,p) = x
    cg(2,p) = y

    if ( 0 < j ) then
      p = p + 1
      cg(1,p) = x
      cg(2,p) = 2.0D+00 * c(2) - y
    end if

    do

      i = i + 1
      x = c(1) + r * real ( 2 * i) / real ( 2 * n + 1)

      if ( r * r < ( x - c(1) )**2 + ( y - c(2) )**2 ) then
        exit
      end if

      p = p + 1
      cg(1,p) = x
      cg(2,p) = y
      p = p + 1
      cg(1,p) = 2.0D+00 * c(1) - x
      cg(2,p) = y

      if ( 0 < j ) then
        p = p + 1
        cg(1,p) = x
        cg(2,p) = 2.0D+00 * c(2) - y
        p = p + 1
        cg(1,p) = 2.0D+00 * c(1) - x
        cg(2,p) = 2.0D+00 * c(2) - y
      end if

    end do

  end do
end

subroutine disk_grid_count ( n, r, c, ng )

!*****************************************************************************80
!
!! DISK_GRID_COUNT counts the grid points inside a disk.
!
!  Discussion:
!
!    The grid is defined by specifying the radius and center of the circle,
!    and the number of subintervals N into which the horizontal radius
!    should be divided.  Thus, a value of N = 2 will result in 5 points
!    along that horizontal line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of subintervals.
!
!    Input, double precision R, the radius of the circle.
!
!    Input, double precision C(2), the coordinates of the center of the circle.
!
!    Output, integer NG, the number of grid points inside 
!    the circle.
!
  implicit none

  double precision c(2)
  integer i
  integer j
  integer n
  integer ng
  double precision r
  double precision x
  double precision y

  ng = 0

  do j = 0, n

    i = 0
    x = c(1)
    y = c(2) + r * real ( 2 * j) / real ( 2 * n + 1)
    ng = ng + 1

    if ( 0 < j ) then
      ng = ng + 1
    end if

    do

      i = i + 1
      x = c(1) + r * real ( 2 * i) / real ( 2 * n + 1)

      if ( r * r < ( x - c(1) )**2 + ( y - c(2) )**2 ) then
        exit
      end if

      ng = ng + 1
      ng = ng + 1
      if ( 0 < j ) then
        ng = ng + 1
        ng = ng + 1
      end if

    end do

  end do
end

subroutine disk_grid_fibonacci ( n, r, c, g )

!*****************************************************************************80
!
!! DISK_GRID_FIBONACCI computes Fibonacci grid points inside a disk.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Swinbank, James Purser,
!    Fibonacci grids: A novel approach to global modelling,
!    Quarterly Journal of the Royal Meteorological Society,
!    Volume 132, Number 619, July 2006 Part B, pages 1769-1793.
!
!  Parameters:
!
!    Input, integer N, the number of points desired.
!
!    Input, double precision R, the radius of the circle.
!
!    Input, double precision C(2), the coordinates of the center of the circle.
!
!    Output, double precision G(2,N), the grid points.
!
  implicit none

  integer n

  double precision c(2)
  double precision g(2,n)
  double precision gr
  double precision gt
  integer i
  double precision phi
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision r
  double precision r0

  r0 = r / sqrt ( real ( n) - 0.5D+00 )
  phi = ( 1.0D+00 + sqrt ( 5.0D+00 ) ) / 2.0D+00

  do i = 1, n
    gr = r0 * sqrt ( real ( i) - 0.5D+00 )
    gt = 2.0D+00 * pi * real ( i) / phi
    g(1,i) = c(1) + gr * cos ( gt )
    g(2,i) = c(2) + gr * sin ( gt )
  end do
end

subroutine r82vec_print_part ( n, a, max_print, title )

!*****************************************************************************80
!
!! R82VEC_PRINT_PART prints "part" of an R82VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vector.
!
!    Input, double precision A(2,N), the vector to be printed.
!
!    Input, integer MAX_PRINT, the maximum number of lines
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer n

  double precision a(2,n)
  integer i
  integer max_print
  character ( len = * ) title

  if ( max_print <= 0 ) then
  end if

  if ( n <= 0 ) then
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(1:2,i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(1:2,i)
    end do
    write ( *, '(a)' ) '  ........  ..............  ..............'
    i = n
    write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(1:2,i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(1:2,i)
    end do
    i = max_print
    write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,a)' ) i, ':', a(1:2,i), &
      '...more entries...'

  end if
end
