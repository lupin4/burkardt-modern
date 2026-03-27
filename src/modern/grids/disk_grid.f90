!> disk_grid -- Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module disk_grid_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: disk_grid, disk_grid_count, disk_grid_fibonacci, r82vec_print_part

contains

  subroutine disk_grid ( n, r, c, ng, cg ) &
        bind(C, name="disk_grid")

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
  !    Input, integer(ip) N, the number of subintervals.
  !
  !    Input, real(dp) R, the radius of the circle.
  !
  !    Input, real(dp) C(2), the coordinates of the center of the circle.
  !
  !    Input, integer(ip) NG, the number of grid points, as determined by
  !    DISK_GRID_COUNT.
  !
  !    Output, real(dp) CG(2,NG), the grid points inside the circle.
  !

    integer(ip), intent(in), value :: ng
    real(dp), intent(in) :: c(2)
    real(dp), intent(out) :: cg(2,ng)
    integer(ip) :: i
    integer(ip) :: j
    integer(ip), intent(in), value :: n
    integer(ip) :: p
    real(dp), intent(in), value :: r
    real(dp) :: x
    real(dp) :: y

    p = 0

    do j = 0, n

      i = 0
      x = c(1)
      y = c(2) + r * real ( 2 * j, dp) / real ( 2 * n + 1, dp)
      p = p + 1
      cg(1,p) = x
      cg(2,p) = y

      if ( j > 0 ) then
        p = p + 1
        cg(1,p) = x
        cg(2,p) = 2.0_dp * c(2) - y
      end if

      do

        i = i + 1
        x = c(1) + r * real ( 2 * i, dp) / real ( 2 * n + 1, dp)

        if ( r * r < ( x - c(1) )**2 + ( y - c(2) )**2 ) then
          exit
        end if

        p = p + 1
        cg(1,p) = x
        cg(2,p) = y
        p = p + 1
        cg(1,p) = 2.0_dp * c(1) - x
        cg(2,p) = y

        if ( j > 0 ) then
          p = p + 1
          cg(1,p) = x
          cg(2,p) = 2.0_dp * c(2) - y
          p = p + 1
          cg(1,p) = 2.0_dp * c(1) - x
          cg(2,p) = 2.0_dp * c(2) - y
        end if

      end do

    end do
  end subroutine disk_grid

  subroutine disk_grid_count ( n, r, c, ng ) &
        bind(C, name="disk_grid_count")

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
  !    Input, integer(ip) N, the number of subintervals.
  !
  !    Input, real(dp) R, the radius of the circle.
  !
  !    Input, real(dp) C(2), the coordinates of the center of the circle.
  !
  !    Output, integer(ip) NG, the number of grid points inside
  !    the circle.
  !

    real(dp), intent(in) :: c(2)
    integer(ip) :: i
    integer(ip) :: j
    integer(ip), intent(in), value :: n
    integer(ip), intent(out) :: ng
    real(dp), intent(in), value :: r
    real(dp) :: x
    real(dp) :: y

    ng = 0

    do j = 0, n

      i = 0
      x = c(1)
      y = c(2) + r * real ( 2 * j, dp) / real ( 2 * n + 1, dp)
      ng = ng + 1

      if ( j > 0 ) then
        ng = ng + 1
      end if

      do

        i = i + 1
        x = c(1) + r * real ( 2 * i, dp) / real ( 2 * n + 1, dp)

        if ( r * r < ( x - c(1) )**2 + ( y - c(2) )**2 ) then
          exit
        end if

        ng = ng + 1
        ng = ng + 1
        if ( j > 0 ) then
          ng = ng + 1
          ng = ng + 1
        end if

      end do

    end do
  end subroutine disk_grid_count

  pure subroutine disk_grid_fibonacci ( n, r, c, g ) &
        bind(C, name="disk_grid_fibonacci")

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
  !    Input, integer(ip) N, the number of points desired.
  !
  !    Input, real(dp) R, the radius of the circle.
  !
  !    Input, real(dp) C(2), the coordinates of the center of the circle.
  !
  !    Output, real(dp) G(2,N), the grid points.
  !

    integer(ip), intent(in), value :: n
    real(dp), intent(in) :: c(2)
    real(dp), intent(out) :: g(2,n)
    real(dp) :: gr
    real(dp) :: gt
    integer(ip) :: i
    real(dp) :: phi
    real(dp), parameter :: pi = 3.141592653589793e+00_dp
    real(dp), intent(in), value :: r
    real(dp) :: r0

    r0 = r / sqrt ( real ( n, dp) - 0.5_dp )
    phi = ( 1.0_dp + sqrt ( 5.0_dp ) ) / 2.0_dp

    do i = 1, n
      gr = r0 * sqrt ( real ( i, dp) - 0.5_dp )
      gt = 2.0_dp * pi * real ( i, dp) / phi
      g(1,i) = c(1) + gr * cos ( gt )
      g(2,i) = c(2) + gr * sin ( gt )
    end do
  end subroutine disk_grid_fibonacci

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

end module disk_grid_mod
