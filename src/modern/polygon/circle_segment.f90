!> circle_segment — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

subroutine circle_segment_angle_from_chord ( r, c, p1, p2, theta )

!*****************************************************************************80
!
!! CIRCLE_SEGMENT_ANGLE_FROM_CHORD computes the angle of a circle segment.
!
!  Discussion:
!
!    Begin with a circle of radius R.  Choose two points P1 and P2 on the
!    circle, and draw the chord P1:P2.  This chord divides the circle
!    into two pieces, each of which is called a circle segment.
!    Consider one of the pieces.  The "angle" of this segment is the angle 
!    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
!    the chord P1:P2 which is closest to C.  The "height" of the segment
!    is the distance from Q to the perimeter of the circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision R, the radius of the circle.
!    0 < R.
!
!    Input, double precision C(2), the center of the circle.
!
!    Input, double precision P1(2), P2(2), the ends of the chord.
!
!    Output, double precision THETA, the angle of the circle segment.
!    0 <= THETA < 2 * PI.
!
  implicit none

  double precision c(2)
  double precision p1(2)
  double precision p2(2)
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision r
  double precision r8_atan
  double precision theta
  double precision v1(2)
  double precision v2(2)
!
!  Compute the radial vectors V1 and V2.
!
  v1(1:2) = p1(1:2) - c(1:2)
  v2(1:2) = p2(1:2) - c(1:2)
!
!  The arc cosine will only give us an answer between 0 and PI.
!
  theta = r8_atan ( v2(2), v2(1) ) - r8_atan ( v1(2), v1(1) )
!
!  Force 0 <= THETA < 2 * PI.
!
  do while ( theta < 0.0D+00 )
    theta = theta + 2.0D+00 * pi
  end do

  do while ( 2.0D+00 * pi <= theta )
    theta = theta - 2.0D+00 * pi
  end do
end

subroutine circle_segment_angle_from_chord_angles ( omega1, omega2, theta )

!*****************************************************************************80
!
!! CIRCLE_SEGMENT_ANGLE_FROM_CHORD_ANGLES computes angle of a circle segment.
!
!  Discussion:
!
!    Begin with a circle of radius R.  Choose two points P1 and P2 on the
!    circle, and draw the chord P1:P2.  This chord divides the circle
!    into two pieces, each of which is called a circle segment.
!    Consider one of the pieces.  The "angle" of this segment is the angle 
!    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
!    the chord P1:P2 which is closest to C.  The "height" of the segment
!    is the distance from Q to the perimeter of the circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision OMEGA1, OMEGA2, the angles of the points P1 
!    and P2.  OMEGA1 <= OMEGA2.
!
!    Output, double precision THETA, the angle of the circle segment.
!    Essentially, THETA = OMEGA2 - OMEGA1.
!
  implicit none

  double precision , parameter :: pi = 3.141592653589793D+00
  double precision omega1
  double precision omega2
  double precision theta

  do while ( omega2 < omega1 )
    omega2 = omega2 + 2.0D+00 * pi
  end do

  theta = omega2 - omega1
end

subroutine circle_segment_angle_from_height ( r, h, theta )

!*****************************************************************************80
!
!! CIRCLE_SEGMENT_ANGLE_FROM_HEIGHT computes the angle of a circle segment.
!
!  Discussion:
!
!    Begin with a circle of radius R.  Choose two points P1 and P2 on the
!    circle, and draw the chord P1:P2.  This chord divides the circle
!    into two pieces, each of which is called a circle segment.
!    Consider one of the pieces.  The "angle" of this segment is the angle 
!    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
!    the chord P1:P2 which is closest to C.  The "height" of the segment
!    is the distance from Q to the perimeter of the circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision R, the radius of the circle.
!    0 < R.
!
!    Input, double precision H, the "height" of the circle segment.
!    0 <= H <= 2 * R.
!
!    Output, double precision THETA, the angle of the circle segment.
!
  implicit none

  double precision h
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision r
  double precision r8_asin
! double precision r8_acos
  double precision theta

  if ( h <= 0.0D+00 ) then

    theta = 0.0D+00

  else if ( h <= r ) then

!   theta = 2.0D+00 * r8_acos ( ( r - h ) / r )
    theta = 2.0D+00 * r8_asin ( sqrt ( r * r - ( r - h ) * ( r - h ) ) / r )

  else if ( h <= 2.0D+00 * r ) then

!   theta = 2.0D+00 * r8_acos ( ( r - h ) / r )
    theta = 2.0D+00 * r8_asin ( sqrt ( r * r - ( r - h ) * ( r - h ) ) / r )
    theta = 2.0D+00 * pi - theta

  else

    theta = 2.0D+00 * pi

  end if
end

subroutine circle_segment_area_from_angle ( r, theta, area )

!*****************************************************************************80
!
!! CIRCLE_SEGMENT_AREA_FROM_ANGLE computes the area of a circle segment.
!
!  Discussion:
!
!    Begin with a circle of radius R.  Choose two points P1 and P2 on the
!    circle, and draw the chord P1:P2.  This chord divides the circle
!    into two pieces, each of which is called a circle segment.
!    Consider one of the pieces.  The "angle" of this segment is the angle 
!    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
!    the chord P1:P2 which is closest to C.  The "height" of the segment
!    is the distance from Q to the perimeter of the circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision R, the radius of the circle.
!    0 < R.
!
!    Input, double precision THETA, the angle of the circle segment.
!
!    Output, double precision AREA, the area of the circle segment.
!
  implicit none

  double precision area
  double precision r
  double precision theta

  area = r * r * ( theta - sin ( theta ) ) / 2.0D+00
end

subroutine circle_segment_area_from_chord ( r, c, p1, p2, area )

!*****************************************************************************80
!
!! CIRCLE_SEGMENT_AREA_FROM_CHORD computes the area of a circle segment.
!
!  Discussion:
!
!    Begin with a circle of radius R.  Choose two points P1 and P2 on the
!    circle, and draw the chord P1:P2.  This chord divides the circle
!    into two pieces, each of which is called a circle segment.
!    Consider one of the pieces.  The "angle" of this segment is the angle 
!    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
!    the chord P1:P2 which is closest to C.  The "height" of the segment
!    is the distance from Q to the perimeter of the circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision R, the radius of the circle.
!    0 < R.
!
!    Input, double precision C(2), the center of the circle.
!
!    Input, double precision P1(2), P2(2), the ends of the chord.
!
!    Output, double precision AREA, the area of the circle segment.
!
  implicit none

  double precision area
  double precision c(2)
  double precision p1(2)
  double precision p2(2)
  double precision r
  double precision theta

  call circle_segment_angle_from_chord ( r, c, p1, p2, theta )

  area = r * r * ( theta - sin ( theta ) ) / 2.0D+00
end

subroutine circle_segment_area_from_height ( r, h, area )

!*****************************************************************************80
!
!! CIRCLE_SEGMENT_AREA_FROM_HEIGHT computes the area of a circle segment.
!
!  Discussion:
!
!    Begin with a circle of radius R.  Choose two points P1 and P2 on the
!    circle, and draw the chord P1:P2.  This chord divides the circle
!    into two pieces, each of which is called a circle segment.
!    Consider one of the pieces.  The "angle" of this segment is the angle 
!    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
!    the chord P1:P2 which is closest to C.  The "height" of the segment
!    is the distance from Q to the perimeter of the circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision R, the radius of the circle.
!    0 < R.
!
!    Input, double precision H, the height of the circle segment.
!    0 <= H <= 2 * R.
!
!    Output, double precision AREA, the area of the circle segment.
!
  implicit none

  double precision area
  double precision h
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision r
  double precision r8_asin
  double precision theta

  if ( h <= 0.0D+00 ) then

    area = 0.0D+00

  else if ( h <= r ) then

    theta = 2.0D+00 * r8_asin ( sqrt ( r * r - ( r - h ) * ( r - h ) ) / r )
    area = r * r * ( theta - sin ( theta ) ) / 2.0D+00

  else if ( h <= 2.0D+00 * r ) then

    theta = 2.0D+00 * r8_asin ( sqrt ( r * r - ( r - h ) * ( r - h ) ) / r )
    theta = 2.0D+00 * pi - theta
    area = r * r * ( theta - sin ( theta ) ) / 2.0D+00

  else

    area = pi * r * r

  end if
end

subroutine circle_segment_area_from_sample ( r, c, p1, p2, n, seed, area )

!*****************************************************************************80
!
!! CIRCLE_SEGMENT_AREA_FROM_SAMPLE computes the area of a circle segment.
!
!  Discussion:
!
!    Begin with a circle of radius R.  Choose two points P1 and P2 on the
!    circle, and draw the chord P1:P2.  This chord divides the circle
!    into two pieces, each of which is called a circle segment.
!    Consider one of the pieces.  The "angle" of this segment is the angle 
!    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
!    the chord P1:P2 which is closest to C.  The "height" of the segment
!    is the distance from Q to the perimeter of the circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision R, the radius of the circle.
!    0 < R.
!
!    Input, double precision C(2,1), the center of the circle.
!
!    Input, double precision P1(2,1), P2(2,1), the ends of the chord.
!
!    Input, integer N, the number of sample points.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, double precision AREA, the area of the circle segment.
!
  implicit none

  integer n

  double precision angle(n)
  double precision area
  double precision c(2)
  integer i
  integer m
  double precision omega1
  double precision omega2
  double precision p(2)
  double precision p1(2)
  double precision p2(2)
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision r
  double precision r2(n)
  double precision r8_atan
  double precision rmh
  integer seed
  double precision vdotp(n)
  double precision x(n)
  double precision y(n)
!
!  Determine the angles of the chord endpoints.
!
  omega1 = r8_atan ( p1(2) - c(2), p1(1) - c(1) )
  do while ( omega1 < 0.0D+00 )
    omega1 = omega1 + 2.0D+00 * pi
  end do

  omega2 = r8_atan ( p2(2) - c(2), p2(1) - c(1) );
  do while ( omega2 < omega1 )
    omega2 = omega2 + 2.0D+00 * pi
  end do
!
!  Get N random points in the circle.
!  To simplify angle measurement, take OMEGA1 as your smallest angle.
!  That way, the check OMEGA1 <= ANGLE <= OMEGA2 will be legitimate.
!
  call r8vec_uniform_01 ( n, seed, angle )
  angle(1:n) = omega1 + 2.0D+00 * pi * angle(1:n)
  
  call r8vec_uniform_01 ( n, seed, r2 )
  r2(1:n) = sqrt ( r2(1:n) )

  x(1:n) = c(1) + r2(1:n) * cos ( angle(1:n) )
  y(1:n) = c(2) + r2(1:n) * sin ( angle(1:n) )
!
!  Determine the vector that touches the circle segment base.
!
  p(1:2) = 0.5D+00 * ( p1(1:2) + p2(1:2) ) - c(1:2)

  rmh = sqrt ( p(1)**2 + p(2)**2 )
  p(1:2) = p(1:2) / rmh

  if ( pi < omega2 - omega1 ) then
    p(1:2) = - p(1:2)
    rmh =  - rmh
  end if
!
!  Compute the projection of each point onto P.
!
  vdotp(1:n) = ( x(1:n) - c(1) ) * p(1) + ( y(1:n) - c(2) ) * p(2)
!
!  Points in the segment lie in the sector, and project at least RMH onto P.
!
  m = 0
  do i = 1, n
    if ( omega1 < angle(i) .and. &
                  angle(i) < omega2 .and. &
         rmh < vdotp(i) ) then
      m = m + 1
    end if
  end do
!
!  The area of the segment is its relative share of the circle area.
!
  area = pi * r**2 * real ( m) / real ( n)
end

subroutine circle_segment_cdf ( r, h, h2, cdf )

!*****************************************************************************80
!
!! CIRCLE_SEGMENT_CDF computes a CDF related to a circle segment.
!
!  Discussion:
!
!    Begin with a circle of radius R.  Choose two points P1 and P2 on the
!    circle, and draw the chord P1:P2.  This chord divides the circle
!    into two pieces, each of which is called a circle segment.
!    Consider one of the pieces.  The "angle" of this segment is the angle 
!    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
!    the chord P1:P2 which is closest to C.  The "height" of the segment
!    is the distance from Q to the perimeter of the circle.
!
!    Now, suppose we want to assign a cumulative density function or CDF
!    based on a variable H2 which measures the height of the circle segment
!    formed by an arbitrary point in the circle segment.  CDF(H2) will
!    measure the probability that a point drawn uniformly at random from
!    the circle segment defines a (smaller) circle segment of height H2.
!
!    If we can define this CDF, then we will be able to sample uniformly
!    from the circle segment, since our "Y" value can be determined from H2,
!    and our X value is chosen uniformly over the horizontal chord 
!    through (0,Y).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision R, the radius of the circle.
!    0 < R.
!
!    Input, double precision H, the "height" of the circle segment.
!    0 <= H <= 2 * R.
!
!    Input, double precision H2, the "height" of the new circle segment 
!    defined by a given point in the circle segment.  0 <= H2 <= H.
!
!    Output, double precision CDF, the cumulative density function for H2, 
!    the probability that a point chosen at random in the circle segment 
!    would define a smaller circle segment of height H2 or less.
!
  implicit none

  double precision a
  double precision a2
  double precision cdf
  double precision h
  double precision h2
  double precision r

  if ( h2 <= 0.0D+00 ) then
    cdf = 0.0D+00
  else if ( h <= h2 ) then
    cdf = 1.0D+00
  else
    call circle_segment_area_from_height ( r, h,  a  )
    call circle_segment_area_from_height ( r, h2, a2 )
    cdf = a2 / a
  end if
end

subroutine circle_segment_centroid_from_chord ( r, c, p1, p2, d )

!*****************************************************************************80
!
!! CIRCLE_SEGMENT_CENTROID_FROM_CHORD computes the centroid of a circle segment.
!
!  Discussion:
!
!    Begin with a circle of radius R.  Choose two points P1 and P2 on the
!    circle, and draw the chord P1:P2.  This chord divides the circle
!    into two pieces, each of which is called a circle segment.
!    Consider one of the pieces.  The "angle" of this segment is the angle 
!    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
!    the chord P1:P2 which is closest to C.  The "height" of the segment
!    is the distance from Q to the perimeter of the circle.
!
!    For this function, we assume that the center of the circle is at (0,0),
!    that the chord is horizontal, and that the circle segment is at the top.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision R, the radius of the circle.
!    0 < R.
!
!    Input, double precision C(2), the center of the circle.
!
!    Input, double precision P1(2), P2(2), the coordinates of the endpoints 
!    of the chord.
!
!    Output, double precision D(2), the coordinates of the centroid.
!
  implicit none

  double precision c(2)
  double precision d(2)
  double precision p1(2)
  double precision p2(2)
  double precision r
  double precision s
  double precision theta
  double precision thetah
  double precision v1(2)
!
!  Get the angle subtended by P1:P2.
!
  call circle_segment_angle_from_chord ( r, c, p1, p2, theta )
!
!  Construct V1, the vector from C to P1.
!
  v1(1:2) = p1(1:2) - c(1:2)
!
!  Rotate V1 through THETA / 2.
!
  thetah = theta / 2.0D+00

  d(1) = cos ( thetah ) * v1(1) - sin ( thetah ) * v1(2)
  d(2) = sin ( thetah ) * v1(1) + cos ( thetah ) * v1(2)
!
!  Scale this vector so it represents the distance to the centroid
!  relative to R.
!
  s = 4.0D+00 * ( sin ( theta / 2.0 ) ) ** 3 &
    / 3.0D+00 / ( theta - sin ( theta ) )

  d(1) = s * d(1)
  d(2) = s * d(2)
!
!  Add the center.
!
  d(1) = d(1) + c(1)
  d(2) = d(2) + c(2)
end

subroutine circle_segment_centroid_from_height ( r, h, d )

!*****************************************************************************80
!
!! CIRCLE_SEGMENT_CENTROID_FROM_HEIGHT computes centroid of a circle segment.
!
!  Discussion:
!
!    Begin with a circle of radius R.  Choose two points P1 and P2 on the
!    circle, and draw the chord P1:P2.  This chord divides the circle
!    into two pieces, each of which is called a circle segment.
!    Consider one of the pieces.  The "angle" of this segment is the angle 
!    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
!    the chord P1:P2 which is closest to C.  The "height" of the segment
!    is the distance from Q to the perimeter of the circle.
!
!    For this function, we assume that the center of the circle is at (0,0).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision R, the radius of the circle.
!    0 < R.
!
!    Input, double precision H, the "height" of the circle segment.
!    0 <= H <= 2 * R.
!
!    Output, double precision D(2), the coordinates of the centroid.
!
  implicit none

  double precision d(2)
  double precision h
  double precision r
  double precision theta
  double precision x
  double precision y

  call circle_segment_angle_from_height ( r, h, theta )

  d(1) = 0.0D+00
  d(2) = 4.0D+00 * r * ( sin ( theta / 2.0D+00 ) ) ** 3 / 3.0D+00 &
    / ( theta - sin ( theta ) )
end

subroutine circle_segment_centroid_from_sample ( r, c, p1, p2, n, seed, d )

!*****************************************************************************80
!
!! CIRCLE_SEGMENT_CENTROID_FROM_SAMPLE estimates a circle segment centroid.
!
!  Discussion:
!
!    Begin with a circle of radius R.  Choose two points P1 and P2 on the
!    circle, and draw the chord P1:P2.  This chord divides the circle
!    into two pieces, each of which is called a circle segment.
!    Consider one of the pieces.  The "angle" of this segment is the angle 
!    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
!    the chord P1:P2 which is closest to C.  The "height" of the segment
!    is the distance from Q to the perimeter of the circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision R, the radius of the circle.
!    0 < R.
!
!    Input, double precision C(2), the center of the circle.
!
!    Input, double precision P1(2), P2(2), the ends of the chord.
!
!    Input, integer N, the number of sample points.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, double precision D(2), the estimated centroid of the 
!    circle segment.
!
  implicit none

  integer n

  double precision c(2)
  double precision d(2)
  double precision p1(2)
  double precision p2(2)
  double precision r
  integer seed
  double precision x(n)
  double precision y(n)

  call circle_segment_sample_from_chord ( r, c, p1, p2, n, seed, x, y )

  d(1) = sum ( x(1:n) ) / real ( n)
  d(2) = sum ( y(1:n) ) / real ( n)
end

subroutine circle_segment_contains_point ( r, c, omega1, omega2, xy, value )

!*****************************************************************************80
!
!! CIRCLE_SEGMENT_CONTAINS_POINT reports whether a point is in a circle segment.
!
!  Discussion:
!
!    Begin with a circle of radius R.  Choose two points P1 and P2 on the
!    circle, and draw the chord P1:P2.  This chord divides the circle
!    into two pieces, each of which is called a circle segment.
!    Consider one of the pieces.  The "angle" of this segment is the angle 
!    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
!    the chord P1:P2 which is closest to C.  The "height" of the segment
!    is the distance from Q to the perimeter of the circle.
!
!    In this function, we allow the circle to have an arbitrary center C,
!    arbitrary radius R, and we describe the points P1 and P2 by specifying
!    their angles OMEGA1 and OMEGA2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision R, the radius of the circle.
!    0 < R.
!
!    Input, double precision C(2), the center of the circle.
!
!    Input, double precision OMEGA1, OMEGA2, the angles of the two points on 
!    the circumference of the circle that define the circle segment.
!    OMEGA1 < OMEGA2 <= OMEGA1 + 2 * PI
!
!    Input, double precision XY(2), a point.
!
!    Output, integer VALUE, is TRUE if the point is inside 
!    the circle segment.
!
  implicit none

  double precision c(2)
  double precision h
  double precision omega1
  double precision omega2
  double precision omegah
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision r
  double precision theta
  double precision v(2)
  double precision v_omega
  double precision v_project
  double precision v_r
  integer value
  double precision xy(2)

  if ( r <= 0.0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CIRCLE_SEGMENT_CONTAINS_POINT - Fatal error!'
    write ( *, '(a)' ) '  R <= 0.0.'
    stop
  end if

  do while ( omega2 < omega1 )
    omega2 = omega2 + 2.0D+00 * pi
  end do
!
!  Compute the vector V = XY - C:
!
  v(1:2) = xy(1:2) - c(1:2)
!
!  a: Point must be inside the circle, so ||V|| <= R.
!
  v_r = sqrt ( v(1) ** 2 + v(2) ** 2 )

  if ( r < v_r ) then
    value = 0
  end if
!
!  b: Angle made by the vector V must be between OMEGA1 and OMEGA2.
!
  v_omega = atan2 ( v(2), v(1) )

  do while ( omega1 <= v_omega + 2.0D+00 * pi )
    v_omega = v_omega - 2.0D+00 * pi
  end do

  do while ( v_omega + 2.0D+00 * pi <= omega1 )
    v_omega = v_omega + 2.0D+00 * pi
  end do

  if ( omega2 < v_omega ) then
    value = 0
  end if
!
!  c: Projection of V onto unit centerline must be at least R-H.
!
  omegah = 0.5D+00 * ( omega1 + omega2 )
  v_project = v(1) * cos ( omegah ) + v(2) * sin ( omegah )

  theta = omega2 - omega1
  call circle_segment_height_from_angle ( r, theta, h )

  if ( v_project < r - h ) then
    value = 0
  end if

  value = 1
end

subroutine circle_segment_height_from_angle ( r, angle, h )

!*****************************************************************************80
!
!! CIRCLE_SEGMENT_HEIGHT_FROM_ANGLE: height of a circle segment from angle.
!
!  Discussion:
!
!    Begin with a circle of radius R.  Choose two points P1 and P2 on the
!    circle, and draw the chord P1:P2.  This chord divides the circle
!    into two pieces, each of which is called a circle segment.
!    Consider one of the pieces.  The "angle" of this segment is the angle 
!    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
!    the chord P1:P2 which is closest to C.  The "height" of the segment
!    is the distance from Q to the perimeter of the circle.
!
!    This function is given the radius R and angle of the segment, and
!    determines the corresponding height.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision R, the radius of the circle.
!    0 < R.
!
!    Input, double precision ANGLE, the angle of the circle segment.
!    0 <= ANGLE <= 2.0 * PI.
!
!    Output, double precision H, the height of the circle segment.
!
  implicit none

  double precision angle
  double precision h
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision r

  if ( angle < 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'CIRCLE_SEGMENT_HEIGHT_FROM_ANGLE - Fatal error!'
    write ( *, '(a)' ) '  ANGLE < 0.0.'
    stop
  end if

  if ( angle == 0.0D+00 ) then
    h = 0.0D+00
  end if

  if ( angle == 2.0D+00 * pi ) then
    h = 2.0D+00 * r
  end if

  if ( 2.0D+00 * pi < angle ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'CIRCLE_SEGMENT_HEIGHT_FROM_ANGLE - Fatal error!'
    write ( *, '(a)' ) '  2.0 * pi < ANGLE.'
    stop
  end if

  if ( r <= 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'CIRCLE_SEGMENT_HEIGHT_FROM_ANGLE - Fatal error!'
    write ( *, '(a)' ) '  R <= 0.0.'
    stop
  end if

  if ( angle <= pi ) then
    h = r * ( 1.0D+00 - cos (                  angle   / 2.0D+00 ) )
  else
    h = r * ( 1.0D+00 + cos ( ( 2.0D+00 * pi - angle ) / 2.0D+00 ) )
  end if
end

subroutine circle_segment_height_from_area ( r, area, h )

!*****************************************************************************80
!
!! CIRCLE_SEGMENT_HEIGHT_FROM_AREA: height of a circle segment from area.
!
!  Discussion:
!
!    Begin with a circle of radius R.  Choose two points P1 and P2 on the
!    circle, and draw the chord P1:P2.  This chord divides the circle
!    into two pieces, each of which is called a circle segment.
!    Consider one of the pieces.  The "angle" of this segment is the angle 
!    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
!    the chord P1:P2 which is closest to C.  The "height" of the segment
!    is the distance from Q to the perimeter of the circle.
!
!    This function is given the radius R and area of the segment, and
!    determines the corresponding height.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision R, the radius of the circle.
!    0 < R.
!
!    Input, double precision AREA, the area of the circle segment.
!    0 <= AREA <= 2.0 * PI * R^2.
!
!    Output, double precision H, the height of the circle segment.
!
  implicit none

  double precision a
  double precision a1
  double precision a2
  double precision area
  double precision area_circle
  double precision eps
  double precision h
  double precision h1
  double precision h2
  integer it
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision r
  double precision r8_epsilon

  if ( area < 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'CIRCLE_SEGMENT_HEIGHT_FROM_AREA - Fatal error!'
    write ( *, '(a)' ) '  AREA < 0.0.'
    stop
  end if

  area_circle = 2.0D+00 * pi * r ** 2

  if ( area == 0.0D+00 ) then
    h = 0.0D+00
  end if

  if ( area == area_circle ) then
    h = 2.0D+00 * r
  end if

  if ( area_circle < area ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'CIRCLE_SEGMENT_HEIGHT_FROM_AREA - Fatal error!'
    write ( *, '(a)' ) '  2.0 * pi * r^2 < AREA.'
    stop
  end if

  if ( r <= 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'CIRCLE_SEGMENT_HEIGHT_FROM_AREA - Fatal error!'
    write ( *, '(a)' ) '  R <= 0.0.'
    stop
  end if

  h1 = 0.0D+00
  call circle_segment_area_from_height ( r, h1, a1 )
  h2 = 2.0D+00 * r
  call circle_segment_area_from_height ( r, h2, a2 )

  it = 0
  eps = r8_epsilon ( )

  do while ( it < 30 )

    h = 0.5D+00 * ( h1 + h2 )
    call circle_segment_area_from_height ( r, h, a )
    it = it + 1

    if ( abs ( a - area ) < sqrt ( eps ) * area ) then
      exit
    end if

    if ( a < area ) then
      h1 = h
      a1 = a
    else
      h2 = h
      a2 = a
    end if

  end do
end

subroutine circle_segment_height_from_chord ( r, c, p1, p2, h )

!*****************************************************************************80
!
!! CIRCLE_SEGMENT_HEIGHT_FROM_CHORD: height of a circle segment from chord.
!
!  Discussion:
!
!    Begin with a circle of radius R.  Choose two points P1 and P2 on the
!    circle, and draw the chord P1:P2.  This chord divides the circle
!    into two pieces, each of which is called a circle segment.
!    Consider one of the pieces.  The "angle" of this segment is the angle 
!    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
!    the chord P1:P2 which is closest to C.  The "height" of the segment
!    is the distance from Q to the perimeter of the circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision R, the radius of the circle.
!    0 < R.
!
!    Input, double precision C(2), the coordinates of the circle center.
!
!    Input, double precision P1(2), P2(2), the coordinates of the 
!    chord endpoints.
!
!    Output, double precision H, the height of the circle segment.
!
  implicit none

  double precision c(2)
  double precision h
  double precision p1(2)
  double precision p2(2)
  double precision r
  double precision theta

  call circle_segment_angle_from_chord ( r, c, p1, p2, theta )

  call circle_segment_height_from_angle ( r, theta, h )
end

subroutine circle_segment_rotation_from_chord ( r, c, p1, p2, alpha )

!*****************************************************************************80
!
!! CIRCLE_SEGMENT_ROTATION_FROM_CHORD computes the rotation of a circle segment.
!
!  Discussion:
!
!    Begin with a circle of radius R.  Choose two points P1 and P2 on the
!    circle, and draw the chord P1:P2.  This chord divides the circle
!    into two pieces, each of which is called a circle segment.
!    Consider one of the pieces.  The "angle" of this segment is the angle 
!    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
!    the chord P1:P2 which is closest to C.  The "height" of the segment
!    is the distance from Q to the perimeter of the circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 July 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision R, the radius of the circle.
!    0 < R.
!
!    Input, double precision C(2), the center of the circle.
!
!    Input, double precision P1(2), P2(2), the ends of the chord.
!    Warning! If P1 = P2, we can't tell whether the segment is the whole
!    circle or none of it!
!
!    Output, double precision ALPHA, the rotation of the circle segment.
!    0 <= ALPHA < 2 * PI.
!
  implicit none

  double precision alpha
  double precision c(2)
  double precision p1(2)
  double precision p2(2)
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision r
  double precision r8_atan
  double precision rho1
  double precision rho2
  double precision theta
  double precision v1(2)
  double precision v2(2)
!
!  Compute the radial vectors V1 and V2.
!
  v1(1:2) = p1(1:2) - c(1:2)
  v2(1:2) = p2(1:2) - c(1:2)
!
!  Use R8_ATAN to guarantee that 0 <= RHO1, RHO2 <= 2 * PI.
!
  rho1 = r8_atan ( v1(2), v1(1) )
  rho2 = r8_atan ( v2(2), v2(1) )
!
!  Force RHO2 to be bigger than RHO1.
!
  do while ( rho2 <= rho1 )
    rho2 = rho2 + 2.0D+00 * pi
  end do
!
!  Compute THETA.
!
  theta = rho2 - rho1
!
!  ALPHA is RHO1, plus half of the angular distance between P1 and P2.
!
  alpha = rho1 + 0.5D+00 * theta

  do while ( 2.0D+00 * pi <= alpha )
    alpha = alpha - 2.0D+00 * pi
  end do
end

subroutine circle_segment_sample_from_chord ( r, c, p1, p2, n, seed, x, y )

!*****************************************************************************80
!
!! CIRCLE_SEGMENT_SAMPLE_FROM_CHORD samples points from a circle segment.
!
!  Discussion:
!
!    Begin with a circle of radius R.  Choose two points P1 and P2 on the
!    circle, and draw the chord P1:P2.  This chord divides the circle
!    into two pieces, each of which is called a circle segment.
!    Consider one of the pieces.  The "angle" of this segment is the angle 
!    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
!    the chord P1:P2 which is closest to C.  The "height" of the segment
!    is the distance from Q to the perimeter of the circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision R, the radius of the circle.
!    0 < R.
!
!    Input, double precision C(2), the center of the circle.
!
!    Input, double precision P1(2), P2(2), the endpoints of the chord.
!
!    Input, integer N, the number of sample points.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, double precision X(N), Y(N), the sample points.
!
  implicit none

  integer n

  double precision c(2)
  double precision c2(2)
  double precision eta(n)
  double precision h
  double precision p1(2)
  double precision p2(2)
  double precision r
  integer seed
  double precision vc(2)
  double precision vr(2)
  double precision x(n)
  double precision xi(n)
  double precision y(n)
!
!  Determine unit vectors VR and VC.
!  VR points to the center of the chord from the radius.
!  VC points along the chord, from P1 to P2.
!
  vr(1:2) = 0.5D+00 * ( p1(1:2) + p2(1:2) ) - c(1:2)
  vr(1:2) = vr(1:2) / sqrt ( vr(1)**2 + vr(2)**2 )
  vc(1:2) = p2(1:2) - p1(1:2)
  vc(1:2) = vc(1:2) / sqrt ( vc(1)**2 + vc(2)**2 )
!
!  Get the height of the circle segment.
!
  c2 = (/ 0.0D+00, 0.0D+00 /)
  call circle_segment_height_from_chord ( r, c2, p1, p2, h )
!
!  Sample (xi,eta) in the reference coordinates, where the chord
!  is horizontal.
!
  call circle_segment_sample_from_height ( r, h, n, seed, xi, eta )
!
!  XI is the left/right coordinate along VC.
!  ETA is the distance along VR.
!
  x(1:n) = c(1) + eta(1:n) * vr(1) + xi(1:n) * vc(1)
  y(1:n) = c(2) + eta(1:n) * vr(2) + xi(1:n) * vc(2)
end

subroutine circle_segment_sample_from_height ( r, h, n, seed, x, y )

!*****************************************************************************80
!
!! CIRCLE_SEGMENT_SAMPLE_FROM_HEIGHT samples points from a circle segment.
!
!  Discussion:
!
!    Begin with a circle of radius R.  Choose two points P1 and P2 on the
!    circle, and draw the chord P1:P2.  This chord divides the circle
!    into two pieces, each of which is called a circle segment.
!    Consider one of the pieces.  The "angle" of this segment is the angle 
!    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
!    the chord P1:P2 which is closest to C.  The "height" of the segment
!    is the distance from Q to the perimeter of the circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision R, the radius of the circle.
!    0 < R.
!
!    Input, double precision H, the height of the circle segment.
!    0 <= H <= 2 * R.
!
!    Input, integer N, the number of sample points.
!
!    Input/output, integer SEED, a seed for the random
!    number generator.
!
!    Output, double precision X(N), Y(N), the sample points.
!
  implicit none

  integer n

  double precision area
  double precision area2(n)
  double precision h
  double precision h2(n)
  integer i
  double precision r
  integer seed
  double precision u(n)
  double precision wh(n)
  double precision x(n)
  double precision y(n)

  call circle_segment_area_from_height ( r, h, area )
!
!  Pick CDF's randomly.
!
  call r8vec_uniform_01 ( n, seed, u )
!
!  Choose points randomly by choosing ordered areas randomly.
!
  area2(1:n) = u(1:n) * area
!
!  Each area corresponds to a height H2.  Find it.
!
  do i = 1, n
    call circle_segment_height_from_area ( r, area2(i), h2(i) )
  end do
!
!  Determine the half-width WH of the segment for each H2.
!
  wh(1:n) = sqrt ( h2(1:n) * ( 2.0 * r - h2(1:n) ) )
!
!  Choose an X position randomly in [-WH,+WH].
!
  call r8vec_uniform_01 ( n, seed, u )
  
  x(1:n) = ( 2.0D+00 * u(1:n) - 1.0 ) * wh(1:n)
!
!  Our circle center is at (0,0).  Our height of H2 is subtracted
!  from the height R at the peak of the circle.  Determine the Y
!  coordinate using this fact.
!
  y(1:n) = r - h2(1:n)
end

subroutine circle_segment_width_from_height ( r, h, w )

!*****************************************************************************80
!
!! CIRCLE_SEGMENT_WIDTH_FROM_HEIGHT computes the width of a circle segment.
!
!  Discussion:
!
!    Begin with a circle of radius R.  Choose two points P1 and P2 on the
!    circle, and draw the chord P1:P2.  This chord divides the circle
!    into two pieces, each of which is called a circle segment.
!    Consider one of the pieces.  The "angle" of this segment is the angle 
!    P1:C:P2, where C is the center of the circle.  Let Q be the point on 
!    the chord P1:P2 which is closest to C.  The "height" of the segment
!    is the distance from Q to the perimeter of the circle.  The "width"
!    of the circle segment is the length of P1:P2.
!
!    This function is given the radius R and height H of the segment, and
!    determines the corresponding width W.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision R, the radius of the circle.
!    0 < R.
!
!    Input, double precision H, the height of the circle segment.
!    0 <= H <= 2 * R.
!
!    Output, double precision W, the width of the circle segment.
!
  implicit none

  double precision h
  double precision r
  double precision w

  w = 2.0D+00 * sqrt ( h * ( 2.0D+00 * r - h ) )
end

subroutine filename_inc ( filename )

!*****************************************************************************80
!
!! FILENAME_INC increments a partially numeric filename.
!
!  Discussion:
!
!    It is assumed that the digits in the name, whether scattered or
!    connected, represent a number that is to be increased by 1 on
!    each call.  If this number is all 9's on input, the output number
!    is all 0's.  Non-numeric letters of the name are unaffected.
!
!    If the name is empty, then the routine stops.
!
!    If the name contains no digits, the empty string is returned.
!
!  Example:
!
!      Input            Output
!      -----            ------
!      'a7to11.txt'     'a7to12.txt'
!      'a7to99.txt'     'a8to00.txt'
!      'a9to99.txt'     'a0to00.txt'
!      'cat.txt'        ' '
!      ' '              STOP!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILENAME.
!    On input, a character string to be incremented.
!    On output, the incremented string.
!
  implicit none

  character c
  integer change
  integer digit
  character ( len = * ) filename
  integer i
  integer lens

  lens = len_trim ( filename )

  if ( lens <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILENAME_INC - Fatal error!'
    write ( *, '(a)' ) '  The input string is empty.'
    stop
  end if

  change = 0

  do i = lens, 1, -1

    c = filename(i:i)

    if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

      change = change + 1

      digit = ichar ( c ) - 48
      digit = digit + 1

      if ( digit == 10 ) then
        digit = 0
      end if

      c = char ( digit + 48 )

      filename(i:i) = c

      if ( c /= '0' ) then
      end if

    end if

  end do
!
!  No digits were found.  Return blank.
!
  if ( change == 0 ) then
    filename = ' '
  end if
end

subroutine gauss ( n, alpha, beta, x, w )

!*****************************************************************************80
!
!! GAUSS computes a Gauss quadrature rule.
!
!  Discussion:
!
!    Given a weight function W encoded by the first N recurrence coefficients 
!    ALPHA and BETA for the associated orthogonal polynomials, the call 
!      call gauss ( n, alpha, beta, x, w ) 
!    generates the nodes and weights of the N-point Gauss quadrature rule 
!    for the weight function W.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 July 2013
!
!  Author:
!
!    Original MATLAB version by Walter Gautschi.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Walter Gautschi,
!    Orthogonal Polynomials: Computation and Approximation,
!    Oxford, 2004,
!    ISBN: 0-19-850672-4,
!    LC: QA404.5 G3555.
!
!  Parameters:
!
!    Input, integer N, the order of the desired quadrature rule.
!
!    Input, double precision ALPHA(N), BETA(N), the alpha and beta recurrence 
!    coefficients for the othogonal polynomials associated with the
!    weight function.
!
!    Output, double precision X(N), W(N), the nodes and  weights of the desired 
!    quadrature rule.  The nodes are listed in increasing order.
!
  implicit none

  integer n

  double precision a(n,n)
  double precision alpha(n)
  double precision beta(n)
  integer i
  integer it_max
  integer it_num
  integer rot_num
  double precision t
  double precision v(n,n)
  double precision w(n)
  double precision x(n)
!
!  Define the tridiagonal Jacobi matrix.
!
  a(1:n,1:n) = 0.0D+00

  do i = 1, n
    a(i,i) = alpha(i)
  end do

  do i = 2, n
    t = sqrt ( beta(i) )
    a(i,i-1) = t
    a(i-1,i) = t
  end do
!
!  Get the eigenvectors and eigenvalues.
!
  it_max = 100

  call jacobi_eigenvalue ( n, a, it_max, v, x, it_num, rot_num )

  w(1:n) = beta(1) * v(1,1:n)**2
end

subroutine jacobi_eigenvalue ( n, a, it_max, v, d, it_num, rot_num )

!*****************************************************************************80
!
!! JACOBI_EIGENVALUE carries out the Jacobi eigenvalue iteration.
!
!  Discussion:
!
!    This function computes the eigenvalues and eigenvectors of a
!    real symmetric matrix, using Rutishauser's modfications of the classical
!    Jacobi rotation method with threshold pivoting. 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 July 2013
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, double precision A(N,N), the matrix, which must be square, real,
!    and symmetric.
!
!    Input, integer IT_MAX, the maximum number of iterations.
!
!    Output, double precision V(N,N), the matrix of eigenvectors.
!
!    Output, double precision D(N), the eigenvalues, in descending order.
!
!    Output, integer IT_NUM, the total number of iterations.
!
!    Output, integer ROT_NUM, the total number of rotations.
!
  implicit none

  integer n

  double precision a(n,n)
  double precision bw(n)
  double precision c
  double precision d(n)
  double precision g
  double precision gapq
  double precision h
  integer i
  integer it_max
  integer it_num
  integer j
  integer k
  integer l
  integer m
  integer p
  integer q
  integer rot_num
  double precision s
  double precision t
  double precision tau
  double precision term
  double precision termp
  double precision termq
  double precision theta
  double precision thresh
  double precision v(n,n)
  double precision w(n)
  double precision zw(n)

  do j = 1, n
    do i = 1, n
      if ( i == j ) then
        v(i,j) = 1.0D+00
      else
        v(i,j) = 0.0D+00
      end if
    end do
  end do

  do i = 1, n
    d(i) = a(i,i)
  end do

  bw(1:n) = d(1:n)
  zw(1:n) = 0.0D+00
  it_num = 0
  rot_num = 0

  do while ( it_num < it_max )

    it_num = it_num + 1
!
!  The convergence threshold is based on the size of the elements in
!  the strict upper triangle of the matrix.
!
    thresh = 0.0D+00
    do j = 1, n
      do i = 1, j - 1
        thresh = thresh + a(i,j) ** 2
      end do
    end do

    thresh = sqrt ( thresh ) / real ( 4 * n)

    if ( thresh == 0.0D+00 ) then
      exit 
    end if

    do p = 1, n
      do q = p + 1, n

        gapq = 10.0D+00 * abs ( a(p,q) )
        termp = gapq + abs ( d(p) )
        termq = gapq + abs ( d(q) )
!
!  Annihilate tiny offdiagonal elements.
!
        if ( 4 < it_num .and. &
             termp == abs ( d(p) ) .and. &
             termq == abs ( d(q) ) ) then

          a(p,q) = 0.0D+00
!
!  Otherwise, apply a rotation.
!
        else if ( thresh <= abs ( a(p,q) ) ) then

          h = d(q) - d(p)
          term = abs ( h ) + gapq

          if ( term == abs ( h ) ) then
            t = a(p,q) / h
          else
            theta = 0.5D+00 * h / a(p,q)
            t = 1.0D+00 / ( abs ( theta ) + sqrt ( 1.0D+00 + theta * theta ) )
            if ( theta < 0.0D+00 ) then 
              t = - t
            end if
          end if

          c = 1.0D+00 / sqrt ( 1.0D+00 + t * t )
          s = t * c
          tau = s / ( 1.0D+00 + c )
          h = t * a(p,q)
!
!  Accumulate corrections to diagonal elements.
!
          zw(p) = zw(p) - h                  
          zw(q) = zw(q) + h
          d(p) = d(p) - h
          d(q) = d(q) + h

          a(p,q) = 0.0D+00
!
!  Rotate, using information from the upper triangle of A only.
!
          do j = 1, p - 1
            g = a(j,p)
            h = a(j,q)
            a(j,p) = g - s * ( h + g * tau )
            a(j,q) = h + s * ( g - h * tau )
          end do

          do j = p + 1, q - 1
            g = a(p,j)
            h = a(j,q)
            a(p,j) = g - s * ( h + g * tau )
            a(j,q) = h + s * ( g - h * tau )
          end do

          do j = q + 1, n
            g = a(p,j)
            h = a(q,j)
            a(p,j) = g - s * ( h + g * tau )
            a(q,j) = h + s * ( g - h * tau )
          end do
!
!  Accumulate information in the eigenvector matrix.
!
          do j = 1, n
            g = v(j,p)
            h = v(j,q)
            v(j,p) = g - s * ( h + g * tau )
            v(j,q) = h + s * ( g - h * tau )
          end do

          rot_num = rot_num + 1

        end if

      end do
    end do

    bw(1:n) = bw(1:n) + zw(1:n)
    d(1:n) = bw(1:n)
    zw(1:n) = 0.0D+00

  end do
!
!  Restore upper triangle of input matrix.
!
  do j = 1, n
    do i = 1, j - 1
      a(i,j) = a(j,i)
    end do
  end do
!
!  Ascending sort the eigenvalues and eigenvectors.
!
  do k = 1, n - 1

    m = k

    do l = k + 1, n
      if ( d(l) < d(m) ) then
        m = l
      end if
    end do

    if ( m /= k ) then

      t    = d(m)
      d(m) = d(k)
      d(k) = t

      w(1:n)   = v(1:n,m)
      v(1:n,m) = v(1:n,k)
      v(1:n,k) = w(1:n)

    end if

  end do
end

subroutine r_jacobi ( n, a, b, alpha, beta )

!*****************************************************************************80
!
!! R_JACOBI computes recurrence coefficients for monic Jacobi polynomials.
!
!  Discussion:
!
!    This function generates the first N recurrence coefficients for monic 
!    Jacobi polynomials with parameters A and B. 
!
!    These polynomials are orthogonal on [-1,1] relative to the weight
!
!      w(x) = (1.0-x)^A * (1.0+x)^B. 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2013
!
!  Author:
!
!    Original MATLAB version by Dirk Laurie, Walter Gautschi.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Walter Gautschi,
!    Orthogonal Polynomials: Computation and Approximation,
!    Oxford, 2004,
!    ISBN: 0-19-850672-4,
!    LC: QA404.5 G3555.
!
!  Parameters:
!
!    Input, integer N, the number of coefficients desired.
!
!    Input, double precision A, B, the parameters for the Jacobi polynomial.
!    -1.0 < A, -1.0 < B.
!
!    Output, double precision ALPHA(N), BETA(N), the first N recurrence
!    coefficients.
!
  implicit none

  integer n

  double precision a
  double precision alpha(n)
  double precision b
  double precision beta(n)
  integer i
  double precision i_r8
  double precision mu
  double precision nab
  double precision nu
  double precision r8_gamma

  if ( a <= -1.0D+00 ) then 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R_JACOBI - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of A.'
    stop
  end if

  if ( b <= -1.0D+00 ) then 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R_JACOBI - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of B.'
    stop
  end if

  nu = ( b - a ) / ( a + b + 2.0D+00 )

  mu = 2.0D+00 ** ( a + b + 1.0D+00 ) &
    * r8_gamma ( a + 1.0D+00 ) &
    * r8_gamma ( b + 1.0D+00 ) &
    / r8_gamma ( a + b + 2.0D+00 )

  alpha(1) = nu
  beta(1) = mu 

  if ( n == 1 ) then
  end if

  do i = 2, n
    i_r8 = real ( i)
    alpha(i) = ( b - a ) * ( b + a ) & 
      / ( 2.0D+00 * ( i_r8 - 1.0D+00 ) + a + b ) &
      / ( 2.0D+00 * i_r8 + a + b )
  end do

  beta(2) = 4.0D+00 * ( a + 1.0D+00 ) * ( b + 1.0D+00 ) &
    / ( a + b + 2.0D+00 ) ** 2 &
    / ( a + b + 3.0D+00 )

  do i = 3, n
    i_r8 = real ( i)
    nab = 2.0D+00 * ( i_r8 - 1.0D+00 ) + a +  b
    beta(i) = 4.0D+00 * ( i_r8 - 1.0D+00 + a ) * ( i_r8 - 1.0D+00 + b ) &
      * ( i_r8 - 1.0D+00 ) * ( i_r8 - 1.0D+00 + a + b ) &
      / nab ** 2 &
      / ( nab + 1.0D+00 ) &
      / ( nab - 1.0D+00 )
  end do
end

subroutine tridisolve ( n, a, b, c, d, x )

!*****************************************************************************80
!
!! TRIDISOLVE solves a tridiagonal system of linear equations.
!
!  Discussion:
!
!    We can describe an NxN tridiagonal matrix by vectors A, B, and C, where
!    A and C are of length N-1.  In that case, a linear system can be
!    represented as
!                        b(1) * x(1) + c(1) * x(2)   = d(1),
!      a(j-1) * x(j-1) + b(j) * x(j) + c(j) * x(j+1) = d(j), j = 2:n-1,
!      a(n-1) * x(n-1) + b(n) * x(n)                 = d(n)
!
!    This function produces the solution vector X.
!
!    This function is derived from Cleve Moler's Matlab suite.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2013
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer N, the order of the linear system.
!
!    Input, double precision A(N-1), B(N), C(N-1), the matrix entries.
!
!    Input, double precision D(N), the right hand side.
!
!    Output, double precision X(N), the solution.
!
  implicit none

  integer n

  double precision a(n-1)
  double precision b(n)
  double precision bi(n)
  double precision c(n-1)
  double precision d(n)
  integer j
  double precision mu
  double precision x(n)

  x(1:n) = d(1:n)

  bi(1:n) = 1.0D+00 / b(1:n)

  do j = 1, n - 1
    mu = a(j) * bi(j)
    b(j+1) = b(j+1) - mu * c(j)
    x(j+1) = x(j+1) - mu * x(j)
  end do

  x(n) = x(n) * bi(n)
  do j = n - 1, 1, -1
    x(j) = ( x(j) - c(j) * x(j+1) ) * bi(j)
  end do
end
