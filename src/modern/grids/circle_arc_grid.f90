!> circle_arc_grid -- Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module circle_arc_grid_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: circle_arc_grid, r82vec_print_part

contains

  pure subroutine circle_arc_grid ( r, c, a, n, xy ) &
        bind(C, name="circle_arc_grid")

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
  !    Input, real(dp) R, the radius of the circle.
  !
  !    Input, real(dp) C(2), the coordinates of the center.
  !
  !    Input, real(dp) A(2), the angle of the first and last
  !    points, in DEGREES.
  !
  !    Input, integer(ip) N, the number of points.
  !
  !    Output, real(dp) XY(2,N), the grid points.
  !

    integer(ip), intent(in), value :: n
    real(dp), intent(in) :: a(2)
    real(dp) :: aj
    real(dp), intent(in) :: c(2)
    integer(ip) :: j
    real(dp), parameter :: pi = 3.141592653589793e+00_dp
    real(dp), intent(in), value :: r
    real(dp), intent(out) :: xy(2,n)

    do j = 1, n

      aj = ( real ( n - j, dp) * a(1)   &
           + real (     j - 1, dp) * a(2) ) &
           / real ( n     - 1, dp)

      xy(1,j) = c(1) + r * cos ( aj * pi / 180.0_dp )
      xy(2,j) = c(2) + r * sin ( aj * pi / 180.0_dp )

    end do
  end subroutine circle_arc_grid

  subroutine r82vec_print_part ( n, a, max_print, title ) &
        bind(C, name="r82vec_print_part")

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
  !    Input, integer(ip) N, the number of entries of the vector.
  !
  !    Input, real(dp) A(2,N), the vector to be printed.
  !
  !    Input, integer(ip) MAX_PRINT, the maximum number of lines
  !    to print.
  !
  !    Input, character ( len = * ) TITLE, a title.
  !

    integer(ip), intent(in), value :: n
    real(dp), intent(in) :: a(2,n)
    integer(ip) :: i
    integer(ip), intent(in), value :: max_print
    character ( len = * ), intent(in) :: title

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
  end subroutine r82vec_print_part

end module circle_arc_grid_mod
