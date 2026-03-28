!> circle_arc_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

subroutine circle_arc_grid ( r, c, a, n, xy )

!*****************************************************************************80
!
!! CIRCLE_ARC_GRID computes grid points along a circular arc.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision R, the radius of the circle.
!
!    Input, double precision C(2), the coordinates of the center.
!
!    Input, double precision A(2), the angle of the first and last
!    points, in DEGREES.
!
!    Input, integer N, the number of points.
!
!    Output, double precision XY(2,N), the grid points.
!
  implicit none

  integer n

  double precision a(2)
  double precision aj
  double precision c(2)
  integer j
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision r
  double precision xy(2,n)

  do j = 1, n

    aj = ( real ( n - j) * a(1)   &
         + real (     j - 1) * a(2) ) &
         / real ( n     - 1)

    xy(1,j) = c(1) + r * cos ( aj * pi / 180.0D+00 )
    xy(2,j) = c(2) + r * sin ( aj * pi / 180.0D+00 )

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
