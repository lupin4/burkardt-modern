!> ellipsoid_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module ellipsoid_grid_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  implicit none
  private

  public :: ellipsoid_grid, ellipsoid_grid_count, i4_ceiling

contains

  subroutine ellipsoid_grid ( n, r, c, ng, xyz )

  !*****************************************************************************80
  !
  !! ELLIPSOID_GRID generates the grid points inside an ellipsoid.
  !
  !  Discussion:
  !
  !    The ellipsoid is specified as
  !
  !      ( ( X - C1 ) / R1 )^2 
  !    + ( ( Y - C2 ) / R2 )^2 
  !    + ( ( Z - C3 ) / R3 )^2 = 1
  !
  !    The user supplies a number N.  There will be N+1 grid points along
  !    the shortest axis.
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
  !    Input, real(real64) R(3), the half axis lengths.
  !
  !    Input, real(real64) C(3), the center of the ellipsoid.
  !
  !    Input, integer(int32) NG, the number of grid points.
  !
  !    Output, real(real64) XYZ(3,NG), the grid point coordinates.
  !

    integer(int32) ng

    real(real64) c(3)
    real(real64) h
    integer(int32) i
    integer(int32) i4_ceiling
    integer(int32) j
    integer(int32) k
    integer(int32) m
    integer(int32) n
    integer(int32) ng2
    integer(int32) ni
    integer(int32) nj
    integer(int32) nk
    integer(int32) np
    real(real64) p(3,8)
    real(real64) r(3)
    real(real64) x
    real(real64) xyz(3,ng)
    real(real64) y
    real(real64) z

    if ( r(1) == minval ( r(1:3) ) ) then
      h = 2.0e+00_real64 * r(1) / real ( 2 * n + 1, real64)
      ni = n
      nj = i4_ceiling ( r(2) / r(1) ) * n
      nk = i4_ceiling ( r(3) / r(1) ) * n
    else if ( r(2) == minval ( r(1:3) ) ) then
      h = 2.0e+00_real64 * r(2) / real ( 2 * n + 1, real64)
      nj = n
      ni = i4_ceiling ( r(1) / r(2) ) * n
      nk = i4_ceiling ( r(3) / r(2) ) * n
    else
      h = 2.0e+00_real64 * r(3) / real ( 2 * n + 1, real64)
      nk = n
      ni = i4_ceiling ( r(1) / r(3) ) * n
      nj = i4_ceiling ( r(2) / r(3) ) * n
    end if

    ng2 = 0

    do k = 0, nk
      z = c(3) + real ( k, real64) * h
      do j = 0, nj
        y = c(2) + real ( j, real64) * h
        do i = 0, ni
          x = c(1) + real ( i, real64) * h
  !
  !  If we have left the ellipsoid, the I loop is completed.
  !
          if ( 1.0e+00_real64 < ( ( x - c(1) ) / r(1) ) ** 2 &
                       + ( ( y - c(2) ) / r(2) ) ** 2 &
                       + ( ( z - c(3) ) / r(3) ) ** 2 ) then
            exit
          end if
  !
  !  At least one point is generated, but more possible by symmetry.
  !
          np = 1
          p(1,np) = x
          p(2,np) = y
          p(3,np) = z

          if ( 0 < i ) then
            do m = 1, np
              p(1,m+np) = 2.0e+00_real64 * c(1) - p(1,m)
              p(2,m+np) = p(2,m)
              p(3,m+np) = p(3,m)
            end do
            np = 2 * np
          end if

          if ( 0 < j ) then
            do m = 1, np
              p(1,m+np) = p(1,m)
              p(2,m+np) = 2.0e+00_real64 * c(2) - p(2,m)
              p(3,m+np) = p(3,m)
            end do
            np = 2 * np
          end if

          if ( 0 < k ) then
            do m = 1, np
              p(1,m+np) = p(1,m)
              p(2,m+np) = p(2,m)
              p(3,m+np) = 2.0e+00_real64 * c(3) - p(3,m)
            end do
            np = 2 * np
          end if

          xyz(1:3,ng2+1:ng2+np) = p(1:3,1:np)
          ng2 = ng2 + np

        end do
      end do
    end do
  end

  subroutine ellipsoid_grid_count ( n, r, c, ng )

  !*****************************************************************************80
  !
  !! ELLIPSOID_GRID_COUNT counts the grid points inside an ellipsoid.
  !
  !  Discussion:
  !
  !    The ellipsoid is specified as
  !
  !      ( ( X - C1 ) / R1 )^2 
  !    + ( ( Y - C2 ) / R2 )^2 
  !    + ( ( Z - C3 ) / R3 )^2 = 1
  !
  !    The user supplies a number N.  There will be N+1 grid points along
  !    the shortest axis.
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
  !    Input, real(real64) R(3), the half axis lengths.
  !
  !    Input, real(real64) C(3), the center of the ellipsoid.
  !
  !    Output, integer(int32) NG, the number of grid points.
  !

    real(real64) c(3)
    real(real64) h
    integer(int32) i
    integer(int32) i4_ceiling
    integer(int32) j
    integer(int32) k
    integer(int32) n
    integer(int32) ng
    integer(int32) ni
    integer(int32) nj
    integer(int32) nk
    integer(int32) np
    real(real64) r(3)
    real(real64) x
    real(real64) y
    real(real64) z

    if ( r(1) == minval ( r(1:3) ) ) then
      h = 2.0e+00_real64 * r(1) / real ( 2 * n + 1, real64)
      ni = n
      nj = i4_ceiling ( r(2) / r(1) ) * n
      nk = i4_ceiling ( r(3) / r(1) ) * n
    else if ( r(2) == minval ( r(1:3) ) ) then
      h = 2.0e+00_real64 * r(2) / real ( 2 * n + 1, real64)
      nj = n
      ni = i4_ceiling ( r(1) / r(2) ) * n
      nk = i4_ceiling ( r(3) / r(2) ) * n
    else
      h = 2.0e+00_real64 * r(3) / real ( 2 * n + 1, real64)
      nk = n
      ni = i4_ceiling ( r(1) / r(3) ) * n
      nj = i4_ceiling ( r(2) / r(3) ) * n
    end if

    ng = 0

    do k = 0, nk
      z = c(3) + real ( k, real64) * h
      do j = 0, nj
        y = c(2) + real ( j, real64) * h
        do i = 0, ni
          x = c(1) + real ( i, real64) * h
  !
  !  If we have left the ellipsoid, the I loop is completed.
  !
          if ( 1.0e+00_real64 < ( ( x - c(1) ) / r(1) ) ** 2 &
                       + ( ( y - c(2) ) / r(2) ) ** 2 &
                       + ( ( z - c(3) ) / r(3) ) ** 2 ) then
            exit
          end if
  !
  !  At least one point is generated, but more possible by symmetry.
  !
          np = 1

          if ( 0 < i ) then
            np = 2 * np
          end if

          if ( 0 < j ) then
            np = 2 * np
          end if

          if ( 0 < k ) then
            np = 2 * np
          end if

          ng = ng + np

        end do
      end do
    end do
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

end module ellipsoid_grid_mod
