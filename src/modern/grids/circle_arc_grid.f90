!> circle_arc_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module circle_arc_grid_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  implicit none
  private

  public :: circle_arc_grid, r82vec_print_part

contains

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
  !    Input, real(real64) R, the radius of the circle.
  !
  !    Input, real(real64) C(2), the coordinates of the center.
  !
  !    Input, real(real64) A(2), the angle of the first and last
  !    points, in DEGREES.
  !
  !    Input, integer(int32) N, the number of points.
  !
  !    Output, real(real64) XY(2,N), the grid points.
  !

    integer(int32) n

    real(real64) a(2)
    real(real64) aj
    real(real64) c(2)
    integer(int32) j
    real(real64), parameter :: pi = 3.141592653589793e+00_real64
    real(real64) r
    real(real64) xy(2,n)

    do j = 1, n

      aj = ( real ( n - j, real64) * a(1)   &
           + real (     j - 1, real64) * a(2) ) &
           / real ( n     - 1, real64)

      xy(1,j) = c(1) + r * cos ( aj * pi / 180.0e+00_real64 )
      xy(2,j) = c(2) + r * sin ( aj * pi / 180.0e+00_real64 )

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
  !    Input, integer(int32) N, the number of entries of the vector.
  !
  !    Input, real(real64) A(2,N), the vector to be printed.
  !
  !    Input, integer(int32) MAX_PRINT, the maximum number of lines
  !    to print.
  !
  !    Input, character ( len = * ) TITLE, a title.
  !

    integer(int32) n

    real(real64) a(2,n)
    integer(int32) i
    integer(int32) max_print
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

end module circle_arc_grid_mod
