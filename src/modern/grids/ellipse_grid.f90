!> ellipse_grid -- Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module ellipse_grid_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: ellipse_grid, ellipse_grid_count, i4_ceiling, r82vec_print_part

contains

  subroutine ellipse_grid ( n, r, c, ng, xy ) &
        bind(C, name="ellipse_grid")

  !*****************************************************************************80
  !
  !! ELLIPSE_GRID generates grid points inside an ellipse.
  !
  !  Discussion:
  !
  !    The ellipse is specified as
  !
  !      ( ( X - C1 ) / R1 )^2 + ( ( Y - C2 ) / R2 )^2 = 1
  !
  !    The user supplies a number N.  There will be N+1 grid points along
  !    the shorter axis.
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
  !    Input, real(dp) R(2), the half axis lengths.
  !
  !    Input, real(dp) C(2), the center of the ellipse.
  !
  !    Input, integer(int64) NG, the number of grid points
  !    inside the ellipse.
  !
  !    Output, real(dp) XY(2,NG), the grid points.
  !

    integer(ip), intent(in), value :: ng
    real(dp), intent(in) :: c(2)
    real(dp) :: h
    integer(ip) :: i
    integer(ip) :: i4_ceiling
    integer(ip) :: j
    integer(ip), intent(in), value :: n
    integer(ip) :: ni
    integer(ip) :: nj
    integer(ip) :: p
    real(dp), intent(in) :: r(2)
    real(dp) :: x
    real(dp), intent(out) :: xy(2,ng)
    real(dp) :: y

    if ( r(1) < r(2) ) then
      h = 2.0_dp * r(1) / real ( 2 * n + 1, dp)
      ni = n
      nj = i4_ceiling ( r(2) / r(1) ) * n
    else
      h = 2.0_dp * r(2) / real ( 2 * n + 1, dp)
      nj = n
      ni = i4_ceiling ( r(1) / r(2) ) * n
    end if

    p = 0

    do j = 0, nj

      i = 0
      x = c(1)
      y = c(2) + real ( j, dp) * h
      p = p + 1
      xy(1,p) = x
      xy(2,p) = y

      if ( j > 0 ) then
        p = p + 1
        xy(1,p) = x
        xy(2,p) = 2.0_dp * c(2) - y
      end if

      do

        i = i + 1
        x = c(1) + real ( i, dp) * h

        if ( 1.0_dp < ( ( x - c(1) ) / r(1) ) ** 2 &
                     + ( ( y - c(2) ) / r(2) ) ** 2 ) then
          exit
        end if

        p = p + 1
        xy(1,p) = x
        xy(2,p) = y
        p = p + 1
        xy(1,p) = 2.0_dp * c(1) - x
        xy(2,p) = y

        if ( j > 0 ) then
          p = p + 1
          xy(1,p) = x
          xy(2,p) = 2.0_dp * c(2) - y
          p = p + 1
          xy(1,p) = 2.0_dp * c(1) - x
          xy(2,p) = 2.0_dp * c(2) - y
        end if

      end do

    end do
  end subroutine ellipse_grid

  subroutine ellipse_grid_count ( n, r, c, ng ) &
        bind(C, name="ellipse_grid_count")

  !*****************************************************************************80
  !
  !! ELLIPSE_GRID_COUNT counts the grid points inside an ellipse.
  !
  !  Discussion:
  !
  !    The ellipse is specified as
  !
  !      ( ( X - C1 ) / R1 )^2 + ( ( Y - C2 ) / R2 )^2 = 1
  !
  !    The user supplies a number N.  There will be N+1 grid points along
  !    the shorter axis.
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
  !    Input, real(dp) R(2), the half axis lengths.
  !
  !    Input, real(dp) C(2), the center of the ellipse.
  !
  !    Output, integer(ip) NG, the number of grid points inside
  !    the ellipse.
  !

    real(dp), intent(in) :: c(2)
    real(dp) :: h
    integer(ip) :: i
    integer(ip) :: i4_ceiling
    integer(ip) :: j
    integer(ip), intent(in), value :: n
    integer(ip), intent(out) :: ng
    integer(ip) :: ni
    integer(ip) :: nj
    integer(ip) :: p
    real(dp), intent(in) :: r(2)
    real(dp) :: x
    real(dp) :: y

    if ( r(1) < r(2) ) then
      h = 2.0_dp * r(1) / real ( 2 * n + 1, dp)
      ni = n
      nj = i4_ceiling ( r(2) / r(1) ) * n
    else
      h = 2.0_dp * r(2) / real ( 2 * n + 1, dp)
      nj = n
      ni = i4_ceiling ( r(1) / r(2) ) * n
    end if

    p = 0

    do j = 0, nj

      i = 0
      x = c(1)
      y = c(2) + real ( j, dp) * h
      p = p + 1

      if ( j > 0 ) then
        p = p + 1
      end if

      do

        i = i + 1
        x = c(1) + real ( i, dp) * h

        if ( 1.0_dp < ( ( x - c(1) ) / r(1) ) ** 2 &
                     + ( ( y - c(2) ) / r(2) ) ** 2 ) then
          exit
        end if

        p = p + 1
        p = p + 1

        if ( j > 0 ) then
          p = p + 1
          p = p + 1
        end if

      end do

    end do

    ng = p
  end subroutine ellipse_grid_count

  pure function i4_ceiling ( r ) &
        bind(C, name="i4_ceiling")

  !*****************************************************************************80
  !
  !! I4_CEILING rounds an R8 "up" (towards +oo) to the next I4.
  !
  !  Example:
  !
  !    R     Value
  !
  !   -1.1  -1
  !   -1.0  -1
  !   -0.9   0
  !    0.0   0
  !    5.0   5
  !    5.1   6
  !    5.9   6
  !    6.0   6
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
  !    Input, real(dp) R, the value to be rounded up.
  !
  !    Output, integer(ip) I4_CEILING, the rounded value.
  !

    integer(ip) :: i4_ceiling
    real(dp), intent(in), value :: r
    integer(ip) :: value

    value = int ( r )
    if ( real ( value, dp) < r ) then
      value = value + 1
    end if

    i4_ceiling = value
  end function i4_ceiling

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

end module ellipse_grid_mod
