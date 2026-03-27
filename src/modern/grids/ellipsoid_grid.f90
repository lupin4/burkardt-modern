!> ellipsoid_grid -- Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module ellipsoid_grid_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: ellipsoid_grid, ellipsoid_grid_count, i4_ceiling

contains

  subroutine ellipsoid_grid ( n, r, c, ng, xyz ) &
        bind(C, name="ellipsoid_grid")

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
  !    Input, integer(ip) N, the number of subintervals.
  !
  !    Input, real(dp) R(3), the half axis lengths.
  !
  !    Input, real(dp) C(3), the center of the ellipsoid.
  !
  !    Input, integer(ip) NG, the number of grid points.
  !
  !    Output, real(dp) XYZ(3,NG), the grid point coordinates.
  !

    integer(ip), intent(in), value :: ng
    real(dp), intent(in) :: c(3)
    real(dp) :: h
    integer(ip) :: i
    integer(ip) :: i4_ceiling
    integer(ip) :: j
    integer(ip) :: k
    integer(ip) :: m
    integer(ip), intent(in), value :: n
    integer(ip) :: ng2
    integer(ip) :: ni
    integer(ip) :: nj
    integer(ip) :: nk
    integer(ip) :: np
    real(dp) :: p(3,8)
    real(dp), intent(in) :: r(3)
    real(dp) :: x
    real(dp), intent(out) :: xyz(3,ng)
    real(dp) :: y
    real(dp) :: z

    if ( r(1) == minval ( r(1:3) ) ) then
      h = 2.0_dp * r(1) / real ( 2 * n + 1, dp)
      ni = n
      nj = i4_ceiling ( r(2) / r(1) ) * n
      nk = i4_ceiling ( r(3) / r(1) ) * n
    else if ( r(2) == minval ( r(1:3) ) ) then
      h = 2.0_dp * r(2) / real ( 2 * n + 1, dp)
      nj = n
      ni = i4_ceiling ( r(1) / r(2) ) * n
      nk = i4_ceiling ( r(3) / r(2) ) * n
    else
      h = 2.0_dp * r(3) / real ( 2 * n + 1, dp)
      nk = n
      ni = i4_ceiling ( r(1) / r(3) ) * n
      nj = i4_ceiling ( r(2) / r(3) ) * n
    end if

    ng2 = 0

    do k = 0, nk
      z = c(3) + real ( k, dp) * h
      do j = 0, nj
        y = c(2) + real ( j, dp) * h
        do i = 0, ni
          x = c(1) + real ( i, dp) * h
  !
  !  If we have left the ellipsoid, the I loop is completed.
  !
          if ( 1.0_dp < ( ( x - c(1) ) / r(1) ) ** 2 &
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

          if ( i > 0 ) then
            do m = 1, np
              p(1,m+np) = 2.0_dp * c(1) - p(1,m)
              p(2,m+np) = p(2,m)
              p(3,m+np) = p(3,m)
            end do
            np = 2 * np
          end if

          if ( j > 0 ) then
            do m = 1, np
              p(1,m+np) = p(1,m)
              p(2,m+np) = 2.0_dp * c(2) - p(2,m)
              p(3,m+np) = p(3,m)
            end do
            np = 2 * np
          end if

          if ( k > 0 ) then
            do m = 1, np
              p(1,m+np) = p(1,m)
              p(2,m+np) = p(2,m)
              p(3,m+np) = 2.0_dp * c(3) - p(3,m)
            end do
            np = 2 * np
          end if

          xyz(1:3,ng2+1:ng2+np) = p(1:3,1:np)
          ng2 = ng2 + np

        end do
      end do
    end do
  end subroutine ellipsoid_grid

  subroutine ellipsoid_grid_count ( n, r, c, ng ) &
        bind(C, name="ellipsoid_grid_count")

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
  !    Input, integer(ip) N, the number of subintervals.
  !
  !    Input, real(dp) R(3), the half axis lengths.
  !
  !    Input, real(dp) C(3), the center of the ellipsoid.
  !
  !    Output, integer(ip) NG, the number of grid points.
  !

    real(dp), intent(in) :: c(3)
    real(dp) :: h
    integer(ip) :: i
    integer(ip) :: i4_ceiling
    integer(ip) :: j
    integer(ip) :: k
    integer(ip), intent(in), value :: n
    integer(ip), intent(out) :: ng
    integer(ip) :: ni
    integer(ip) :: nj
    integer(ip) :: nk
    integer(ip) :: np
    real(dp), intent(in) :: r(3)
    real(dp) :: x
    real(dp) :: y
    real(dp) :: z

    if ( r(1) == minval ( r(1:3) ) ) then
      h = 2.0_dp * r(1) / real ( 2 * n + 1, dp)
      ni = n
      nj = i4_ceiling ( r(2) / r(1) ) * n
      nk = i4_ceiling ( r(3) / r(1) ) * n
    else if ( r(2) == minval ( r(1:3) ) ) then
      h = 2.0_dp * r(2) / real ( 2 * n + 1, dp)
      nj = n
      ni = i4_ceiling ( r(1) / r(2) ) * n
      nk = i4_ceiling ( r(3) / r(2) ) * n
    else
      h = 2.0_dp * r(3) / real ( 2 * n + 1, dp)
      nk = n
      ni = i4_ceiling ( r(1) / r(3) ) * n
      nj = i4_ceiling ( r(2) / r(3) ) * n
    end if

    ng = 0

    do k = 0, nk
      z = c(3) + real ( k, dp) * h
      do j = 0, nj
        y = c(2) + real ( j, dp) * h
        do i = 0, ni
          x = c(1) + real ( i, dp) * h
  !
  !  If we have left the ellipsoid, the I loop is completed.
  !
          if ( 1.0_dp < ( ( x - c(1) ) / r(1) ) ** 2 &
                       + ( ( y - c(2) ) / r(2) ) ** 2 &
                       + ( ( z - c(3) ) / r(3) ) ** 2 ) then
            exit
          end if
  !
  !  At least one point is generated, but more possible by symmetry.
  !
          np = 1

          if ( i > 0 ) then
            np = 2 * np
          end if

          if ( j > 0 ) then
            np = 2 * np
          end if

          if ( k > 0 ) then
            np = 2 * np
          end if

          ng = ng + np

        end do
      end do
    end do
  end subroutine ellipsoid_grid_count

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

end module ellipsoid_grid_mod
