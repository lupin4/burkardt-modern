!> ellipsoid_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

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
!    Input, integer N, the number of subintervals.
!
!    Input, double precision R(3), the half axis lengths.
!
!    Input, double precision C(3), the center of the ellipsoid.
!
!    Input, integer NG, the number of grid points.
!
!    Output, double precision XYZ(3,NG), the grid point coordinates.
!
  implicit none

  integer ng

  double precision c(3)
  double precision h
  integer i
  integer i4_ceiling
  integer j
  integer k
  integer m
  integer n
  integer ng2
  integer ni
  integer nj
  integer nk
  integer np
  double precision p(3,8)
  double precision r(3)
  double precision x
  double precision xyz(3,ng)
  double precision y
  double precision z

  if ( r(1) == minval ( r(1:3) ) ) then
    h = 2.0D+00 * r(1) / real ( 2 * n + 1)
    ni = n
    nj = i4_ceiling ( r(2) / r(1) ) * n
    nk = i4_ceiling ( r(3) / r(1) ) * n
  else if ( r(2) == minval ( r(1:3) ) ) then
    h = 2.0D+00 * r(2) / real ( 2 * n + 1)
    nj = n
    ni = i4_ceiling ( r(1) / r(2) ) * n
    nk = i4_ceiling ( r(3) / r(2) ) * n
  else
    h = 2.0D+00 * r(3) / real ( 2 * n + 1)
    nk = n
    ni = i4_ceiling ( r(1) / r(3) ) * n
    nj = i4_ceiling ( r(2) / r(3) ) * n
  end if

  ng2 = 0

  do k = 0, nk
    z = c(3) + real ( k) * h
    do j = 0, nj
      y = c(2) + real ( j) * h
      do i = 0, ni
        x = c(1) + real ( i) * h
!
!  If we have left the ellipsoid, the I loop is completed.
!
        if ( 1.0D+00 < ( ( x - c(1) ) / r(1) ) ** 2 &
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
            p(1,m+np) = 2.0D+00 * c(1) - p(1,m)
            p(2,m+np) = p(2,m)
            p(3,m+np) = p(3,m)
          end do
          np = 2 * np
        end if

        if ( 0 < j ) then
          do m = 1, np
            p(1,m+np) = p(1,m)
            p(2,m+np) = 2.0D+00 * c(2) - p(2,m)
            p(3,m+np) = p(3,m)
          end do
          np = 2 * np
        end if

        if ( 0 < k ) then
          do m = 1, np
            p(1,m+np) = p(1,m)
            p(2,m+np) = p(2,m)
            p(3,m+np) = 2.0D+00 * c(3) - p(3,m)
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
!    Input, integer N, the number of subintervals.
!
!    Input, double precision R(3), the half axis lengths.
!
!    Input, double precision C(3), the center of the ellipsoid.
!
!    Output, integer NG, the number of grid points.
!
  implicit none

  double precision c(3)
  double precision h
  integer i
  integer i4_ceiling
  integer j
  integer k
  integer n
  integer ng
  integer ni
  integer nj
  integer nk
  integer np
  double precision r(3)
  double precision x
  double precision y
  double precision z

  if ( r(1) == minval ( r(1:3) ) ) then
    h = 2.0D+00 * r(1) / real ( 2 * n + 1)
    ni = n
    nj = i4_ceiling ( r(2) / r(1) ) * n
    nk = i4_ceiling ( r(3) / r(1) ) * n
  else if ( r(2) == minval ( r(1:3) ) ) then
    h = 2.0D+00 * r(2) / real ( 2 * n + 1)
    nj = n
    ni = i4_ceiling ( r(1) / r(2) ) * n
    nk = i4_ceiling ( r(3) / r(2) ) * n
  else
    h = 2.0D+00 * r(3) / real ( 2 * n + 1)
    nk = n
    ni = i4_ceiling ( r(1) / r(3) ) * n
    nj = i4_ceiling ( r(2) / r(3) ) * n
  end if

  ng = 0

  do k = 0, nk
    z = c(3) + real ( k) * h
    do j = 0, nj
      y = c(2) + real ( j) * h
      do i = 0, ni
        x = c(1) + real ( i) * h
!
!  If we have left the ellipsoid, the I loop is completed.
!
        if ( 1.0D+00 < ( ( x - c(1) ) / r(1) ) ** 2 &
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
