!> ellipse_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module ellipse_grid_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  implicit none
  private

  public :: ellipse_grid, ellipse_grid_count, i4_ceiling, r82vec_print_part

contains

  subroutine ellipse_grid ( n, r, c, ng, xy )

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
  !    Input, integer(int32) N, the number of subintervals.
  !
  !    Input, real(real64) R(2), the half axis lengths.
  !
  !    Input, real(real64) C(2), the center of the ellipse.
  !
  !    Input, integer(int64) NG, the number of grid points 
  !    inside the ellipse.
  !
  !    Output, real(real64) XY(2,NG), the grid points.
  !

    integer(int32) ng

    real(real64) c(2)
    real(real64) h
    integer(int32) i
    integer(int32) i4_ceiling
    integer(int32) j
    integer(int32) n
    integer(int32) ni
    integer(int32) nj
    integer(int32) p
    real(real64) r(2)
    real(real64) x
    real(real64) xy(2,ng)
    real(real64) y

    if ( r(1) < r(2) ) then
      h = 2.0e+00_real64 * r(1) / real ( 2 * n + 1, real64)
      ni = n
      nj = i4_ceiling ( r(2) / r(1) ) * n
    else
      h = 2.0e+00_real64 * r(2) / real ( 2 * n + 1, real64)
      nj = n
      ni = i4_ceiling ( r(1) / r(2) ) * n
    end if

    p = 0

    do j = 0, nj

      i = 0
      x = c(1)
      y = c(2) + real ( j, real64) * h
      p = p + 1
      xy(1,p) = x
      xy(2,p) = y

      if ( 0 < j ) then
        p = p + 1
        xy(1,p) = x
        xy(2,p) = 2.0e+00_real64 * c(2) - y
      end if

      do

        i = i + 1
        x = c(1) + real ( i, real64) * h

        if ( 1.0e+00_real64 < ( ( x - c(1) ) / r(1) ) ** 2 &
                     + ( ( y - c(2) ) / r(2) ) ** 2 ) then
          exit
        end if

        p = p + 1
        xy(1,p) = x
        xy(2,p) = y
        p = p + 1
        xy(1,p) = 2.0e+00_real64 * c(1) - x
        xy(2,p) = y

        if ( 0 < j ) then
          p = p + 1
          xy(1,p) = x
          xy(2,p) = 2.0e+00_real64 * c(2) - y
          p = p + 1
          xy(1,p) = 2.0e+00_real64 * c(1) - x
          xy(2,p) = 2.0e+00_real64 * c(2) - y
        end if

      end do

    end do
  end

  subroutine ellipse_grid_count ( n, r, c, ng )

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
  !    Input, integer(int32) N, the number of subintervals.
  !
  !    Input, real(real64) R(2), the half axis lengths.
  !
  !    Input, real(real64) C(2), the center of the ellipse.
  !
  !    Output, integer(int32) NG, the number of grid points inside 
  !    the ellipse.
  !

    real(real64) c(2)
    real(real64) h
    integer(int32) i
    integer(int32) i4_ceiling
    integer(int32) j
    integer(int32) n
    integer(int32) ng
    integer(int32) ni
    integer(int32) nj
    integer(int32) p
    real(real64) r(2)
    real(real64) x
    real(real64) y

    if ( r(1) < r(2) ) then
      h = 2.0e+00_real64 * r(1) / real ( 2 * n + 1, real64)
      ni = n
      nj = i4_ceiling ( r(2) / r(1) ) * n
    else
      h = 2.0e+00_real64 * r(2) / real ( 2 * n + 1, real64)
      nj = n
      ni = i4_ceiling ( r(1) / r(2) ) * n
    end if

    p = 0

    do j = 0, nj

      i = 0
      x = c(1)
      y = c(2) + real ( j, real64) * h
      p = p + 1

      if ( 0 < j ) then
        p = p + 1
      end if

      do

        i = i + 1
        x = c(1) + real ( i, real64) * h

        if ( 1.0e+00_real64 < ( ( x - c(1) ) / r(1) ) ** 2 &
                     + ( ( y - c(2) ) / r(2) ) ** 2 ) then
          exit
        end if

        p = p + 1
        p = p + 1

        if ( 0 < j ) then
          p = p + 1
          p = p + 1
        end if

      end do

    end do

    ng = p
  end

  function i4_ceiling ( r )

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
  !    Input, real(real64) R, the value to be rounded up.
  !
  !    Output, integer(int32) I4_CEILING, the rounded value.
  !

    integer(int32) i4_ceiling
    real(real64) r
    integer(int32) value

    value = int ( r )
    if ( real ( value, real64) < r ) then
      value = value + 1
    end if

    i4_ceiling = value
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
    character ( len = * )  title

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

end module ellipse_grid_mod
