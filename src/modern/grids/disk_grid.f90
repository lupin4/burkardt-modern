!> disk_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module disk_grid_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  implicit none
  private

  public :: disk_grid, disk_grid_count, disk_grid_fibonacci, r82vec_print_part

contains

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
  !    Input, integer(int32) N, the number of subintervals.
  !
  !    Input, real(real64) R, the radius of the circle.
  !
  !    Input, real(real64) C(2), the coordinates of the center of the circle.
  !
  !    Input, integer(int32) NG, the number of grid points, as determined by
  !    DISK_GRID_COUNT.
  !
  !    Output, real(real64) CG(2,NG), the grid points inside the circle.
  !

    integer(int32) ng

    real(real64) c(2)
    real(real64) cg(2,ng)
    integer(int32) i
    integer(int32) j
    integer(int32) n
    integer(int32) p
    real(real64) r
    real(real64) x
    real(real64) y

    p = 0

    do j = 0, n

      i = 0
      x = c(1)
      y = c(2) + r * real ( 2 * j, real64) / real ( 2 * n + 1, real64)
      p = p + 1
      cg(1,p) = x
      cg(2,p) = y

      if ( 0 < j ) then
        p = p + 1
        cg(1,p) = x
        cg(2,p) = 2.0e+00_real64 * c(2) - y
      end if

      do

        i = i + 1
        x = c(1) + r * real ( 2 * i, real64) / real ( 2 * n + 1, real64)

        if ( r * r < ( x - c(1) )**2 + ( y - c(2) )**2 ) then
          exit
        end if

        p = p + 1
        cg(1,p) = x
        cg(2,p) = y
        p = p + 1
        cg(1,p) = 2.0e+00_real64 * c(1) - x
        cg(2,p) = y

        if ( 0 < j ) then
          p = p + 1
          cg(1,p) = x
          cg(2,p) = 2.0e+00_real64 * c(2) - y
          p = p + 1
          cg(1,p) = 2.0e+00_real64 * c(1) - x
          cg(2,p) = 2.0e+00_real64 * c(2) - y
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
  !    Input, integer(int32) N, the number of subintervals.
  !
  !    Input, real(real64) R, the radius of the circle.
  !
  !    Input, real(real64) C(2), the coordinates of the center of the circle.
  !
  !    Output, integer(int32) NG, the number of grid points inside 
  !    the circle.
  !

    real(real64) c(2)
    integer(int32) i
    integer(int32) j
    integer(int32) n
    integer(int32) ng
    real(real64) r
    real(real64) x
    real(real64) y

    ng = 0

    do j = 0, n

      i = 0
      x = c(1)
      y = c(2) + r * real ( 2 * j, real64) / real ( 2 * n + 1, real64)
      ng = ng + 1

      if ( 0 < j ) then
        ng = ng + 1
      end if

      do

        i = i + 1
        x = c(1) + r * real ( 2 * i, real64) / real ( 2 * n + 1, real64)

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
  !    Input, integer(int32) N, the number of points desired.
  !
  !    Input, real(real64) R, the radius of the circle.
  !
  !    Input, real(real64) C(2), the coordinates of the center of the circle.
  !
  !    Output, real(real64) G(2,N), the grid points.
  !

    integer(int32) n

    real(real64) c(2)
    real(real64) g(2,n)
    real(real64) gr
    real(real64) gt
    integer(int32) i
    real(real64) phi
    real(real64), parameter :: pi = 3.141592653589793e+00_real64
    real(real64) r
    real(real64) r0

    r0 = r / sqrt ( real ( n, real64) - 0.5e+00_real64 )
    phi = ( 1.0e+00_real64 + sqrt ( 5.0e+00_real64 ) ) / 2.0e+00_real64

    do i = 1, n
      gr = r0 * sqrt ( real ( i, real64) - 0.5e+00_real64 )
      gt = 2.0e+00_real64 * pi * real ( i, real64) / phi
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

end module disk_grid_mod
