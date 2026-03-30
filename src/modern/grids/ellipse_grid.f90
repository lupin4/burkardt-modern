!> ellipse_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

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
!    Input, integer N, the number of subintervals.
!
!    Input, double precision R(2), the half axis lengths.
!
!    Input, double precision C(2), the center of the ellipse.
!
!    Input, integer(8) NG, the number of grid points 
!    inside the ellipse.
!
!    Output, double precision XY(2,NG), the grid points.
!
  implicit none

  integer ng

  double precision c(2)
  double precision h
  integer i
  integer i4_ceiling
  integer j
  integer n
  integer ni
  integer nj
  integer p
  double precision r(2)
  double precision x
  double precision xy(2,ng)
  double precision y

  if ( r(1) < r(2) ) then
    h = 2.0D+00 * r(1) / real ( 2 * n + 1)
    ni = n
    nj = i4_ceiling ( r(2) / r(1) ) * n
  else
    h = 2.0D+00 * r(2) / real ( 2 * n + 1)
    nj = n
    ni = i4_ceiling ( r(1) / r(2) ) * n
  end if

  p = 0

  do j = 0, nj

    i = 0
    x = c(1)
    y = c(2) + real ( j) * h
    p = p + 1
    xy(1,p) = x
    xy(2,p) = y

    if ( 0 < j ) then
      p = p + 1
      xy(1,p) = x
      xy(2,p) = 2.0D+00 * c(2) - y
    end if

    do

      i = i + 1
      x = c(1) + real ( i) * h

      if ( 1.0D+00 < ( ( x - c(1) ) / r(1) ) ** 2 &
                   + ( ( y - c(2) ) / r(2) ) ** 2 ) then
        exit
      end if

      p = p + 1
      xy(1,p) = x
      xy(2,p) = y
      p = p + 1
      xy(1,p) = 2.0D+00 * c(1) - x
      xy(2,p) = y

      if ( 0 < j ) then
        p = p + 1
        xy(1,p) = x
        xy(2,p) = 2.0D+00 * c(2) - y
        p = p + 1
        xy(1,p) = 2.0D+00 * c(1) - x
        xy(2,p) = 2.0D+00 * c(2) - y
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
!    Input, integer N, the number of subintervals.
!
!    Input, double precision R(2), the half axis lengths.
!
!    Input, double precision C(2), the center of the ellipse.
!
!    Output, integer NG, the number of grid points inside 
!    the ellipse.
!
  implicit none

  double precision c(2)
  double precision h
  integer i
  integer i4_ceiling
  integer j
  integer n
  integer ng
  integer ni
  integer nj
  integer p
  double precision r(2)
  double precision x
  double precision y

  if ( r(1) < r(2) ) then
    h = 2.0D+00 * r(1) / real ( 2 * n + 1)
    ni = n
    nj = i4_ceiling ( r(2) / r(1) ) * n
  else
    h = 2.0D+00 * r(2) / real ( 2 * n + 1)
    nj = n
    ni = i4_ceiling ( r(1) / r(2) ) * n
  end if

  p = 0

  do j = 0, nj

    i = 0
    x = c(1)
    y = c(2) + real ( j) * h
    p = p + 1

    if ( 0 < j ) then
      p = p + 1
    end if

    do

      i = i + 1
      x = c(1) + real ( i) * h

      if ( 1.0D+00 < ( ( x - c(1) ) / r(1) ) ** 2 &
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
