!> dutch — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

function angle_deg_2d ( p1, p2, p3 )

!*****************************************************************************80
!
!! ANGLE_DEG_2D returns the angle swept out between two rays in 2D.
!
!  Discussion:
!
!    Except for the zero angle case, it should be true that
!
!      ANGLE_DEG_2D ( P1, P2, P3 ) + ANGLE_DEG_2D ( P3, P2, P1 ) = 360.0
!
!        P1
!        /
!       /
!      /
!     /
!    P2--------->P3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision P1(2), P2(2), P3(2), define the rays
!    P1 - P2 and P3 - P2 which define the angle.
!
!    Output, double precision ANGLE_DEG_2D, the angle swept out by the
!    rays, measured in degrees.  0 <= ANGLE_DEG_2D < 360.  If either ray
!    has zero length, then ANGLE_DEG_2D is set to 0.
!
  implicit none

  integer , parameter :: dim_num = 2

  double precision angle_deg_2d
  double precision angle_rad_2d
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision radians_to_degrees
  double precision p(dim_num)
  double precision p1(dim_num)
  double precision p2(dim_num)
  double precision p3(dim_num)

  p(1) = ( p3(1) - p2(1) ) * ( p1(1) - p2(1) ) &
       + ( p3(2) - p2(2) ) * ( p1(2) - p2(2) )

  p(2) = ( p3(1) - p2(1) ) * ( p1(2) - p2(2) ) &
       - ( p3(2) - p2(2) ) * ( p1(1) - p2(1) )

  if ( p(1) == 0.0D+00 .and. p(2) == 0.0D+00 ) then
    angle_deg_2d = 0.0D+00
  end if

  angle_rad_2d = atan2 ( p(2), p(1) )

  if ( angle_rad_2d < 0.0D+00 ) then
    angle_rad_2d = angle_rad_2d + 2.0D+00 * pi
  end if

  angle_deg_2d = radians_to_degrees ( angle_rad_2d )
end

function angle_rad_2d ( p1, p2, p3 )

!*****************************************************************************80
!
!! ANGLE_RAD_2D returns the angle in radians swept out between two rays in 2D.
!
!  Discussion:
!
!    Except for the zero angle case, it should be true that
!
!      ANGLE_RAD_2D ( P1, P2, P3 ) + ANGLE_RAD_2D ( P3, P2, P1 ) = 2 * PI
!
!        P1
!        /
!       /
!      /
!     /
!    P2--------->P3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision P1(2), P2(2), P3(2), define the rays
!    P1 - P2 and P3 - P2 which define the angle.
!
!    Output, double precision ANGLE_RAD_2D, the angle swept out by the rays,
!    in radians.  0 <= ANGLE_RAD_2D < 2 * PI.  If either ray has zero
!    length, then ANGLE_RAD_2D is set to 0.
!
  implicit none

  integer , parameter :: dim_num = 2

  double precision angle_rad_2d
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision p(dim_num)
  double precision p1(dim_num)
  double precision p2(dim_num)
  double precision p3(dim_num)

  p(1) = ( p3(1) - p2(1) ) * ( p1(1) - p2(1) ) &
       + ( p3(2) - p2(2) ) * ( p1(2) - p2(2) )

  p(2) = ( p3(1) - p2(1) ) * ( p1(2) - p2(2) ) &
       - ( p3(2) - p2(2) ) * ( p1(1) - p2(1) )

  if ( p(1) == 0.0D+00 .and. p(2) == 0.0D+00 ) then
    angle_rad_2d = 0.0D+00
  end if

  angle_rad_2d = atan2 ( p(2), p(1) )

  if ( angle_rad_2d < 0.0D+00 ) then
    angle_rad_2d = angle_rad_2d + 2.0D+00 * pi
  end if
end

subroutine circle_dia2imp_2d ( x1, y1, x2, y2, r, cx, cy )

!*****************************************************************************80
!
!! CIRCLE_DIA2IMP_2D converts a diameter to an implicit circle in 2D.
!
!  Discussion:
!
!    The diameter form of a circle is:
!
!      (X1,Y1) and (X2,Y2) are endpoints of a diameter.
!
!    The implicit form of a circle in 2D is:
!
!      (X-CX)**2 + (Y-CY)**2 = R**2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision X1, Y1, X2, Y2, are the X and Y coordinates
!    of two points which form a diameter of the circle.
!
!    Output, double precision R, the computed radius of the circle.
!
!    Output, double precision CX, CY, the computed center of the circle.
!
  implicit none

  double precision r
  double precision x1
  double precision x2
  double precision cx
  double precision y1
  double precision y2
  double precision cy

  r = 0.5D+00 * sqrt ( ( x1 - x2 )**2 + ( y1 - y2 )**2 )

  cx = 0.5D+00 * ( x1 + x2 )
  cy = 0.5D+00 * ( y1 + y2 )
end

subroutine circle_exp2imp_2d ( x1, y1, x2, y2, x3, y3, r, cx, cy )

!*****************************************************************************80
!
!! CIRCLE_EXP2IMP_2D converts a circle from explicit to implicit form in 2D.
!
!  Formula:
!
!    The explicit form of a circle in 2D is:
!
!      The circle passing through (X1,Y1), (X2,Y2), (X3,Y3).
!
!    The implicit form of a circle in 2D is:
!
!      (X-CX)**2 + (Y-CY)**2 = R**2
!
!  Discussion:
!
!    Any three points define a circle, as long as they don't lie on a straight
!    line.  (If the points do lie on a straight line, we could stretch the
!    definition of a circle to allow an infinite radius and a center at
!    some infinite point.)
!
!    Instead of the formulas used here, you can use the linear system
!    approach in the routine TRIANGLE_OUTCIRCLE_2D.
!
!    The diameter of the circle can be found by solving a 2 by 2 linear system.
!    This is because the vectors P2 - P1 and P3 - P1 are secants of the circle,
!    and each forms a right triangle with the diameter.  Hence, the dot product
!    of P2 - P1 with the diameter is equal to the square of the length
!    of P2 - P1, and similarly for P3 - P1.  These two equations determine the
!    diameter vector originating at P1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph O'Rourke,
!    Computational Geometry,
!    Cambridge University Press,
!    Second Edition, 1998, page 187.
!
!  Parameters:
!
!    Input, double precision X1, Y1, X2, Y2, X3, Y3, are the X and Y
!    coordinates of three points that lie on the circle.  These points should be
!    distinct, and not collinear.
!
!    Output, double precision R, the radius of the circle.  Normally, R will
!    be positive.  R is returned as -1 in the unlikely event that the points are
!    numerically collinear.
!
!    Output, double precision CX, CY, the center of the circle.
!
  implicit none

  double precision a
  double precision b
  double precision c
  double precision d
  double precision e
  double precision f
  double precision g
  double precision r
  double precision x1
  double precision x2
  double precision x3
  double precision cx
  double precision y1
  double precision y2
  double precision y3
  double precision cy

  a = x2 - x1
  b = y2 - y1
  c = x3 - x1
  d = y3 - y1

  e = a * ( x1 + x2 ) + b * ( y1 + y2 )
  f = c * ( x1 + x3 ) + d * ( y1 + y3 )
!
!  Our formula is:
!
!    G = a * ( d - b ) - b * ( c - a )
!
!  but we get slightly better results using the original data.
!
  g = a * ( y3 - y2 ) - b * ( x3 - x2 )
!
!  We check for collinearity.  A more useful check would compare the
!  absolute value of G to a small quantity.
!
  if ( g == 0.0D+00 ) then
    cx = 0.0D+00
    cy = 0.0D+00
    r = -1.0D+00
  end if
!
!  The center is halfway along the diameter vector from (X1,Y1).
!
  cx = 0.5D+00 * ( d * e - b * f ) / g
  cy = 0.5D+00 * ( a * f - c * e ) / g
!
!  Knowing the center, the radius is now easy to compute.
!
  r = sqrt ( ( x1 - cx )**2 + ( y1 - cy )**2 )
end

subroutine circle_imp_contains_point_2d ( r, cx, cy, x, y, inside )

!*****************************************************************************80
!
!! CIRCLE_IMP_CONTAINS_POINT_2D: does an implicit circle contains a point in 2D?
!
!  Formula:
!
!    An implicit circle in 2D satisfies the equation:
!
!      ( X - CX )^2 + ( Y - CY )^2 = R^2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision R, the radius of the circle.
!
!    Input, double precision CX, CY, the coordinates of the center
!    of the circle.
!
!    Input, double precision X, Y, the point to be checked.
!
!    Output, logical INSIDE, is TRUE if the point is inside or on the circle,
!    FALSE otherwise.
!
  implicit none

  logical inside
  double precision r
  double precision x
  double precision cx
  double precision y
  double precision cy

  if ( ( x - cx )**2 + ( y - cy )**2 <= r * r ) then
    inside = .true.
  else
    inside = .false.
  end if
end

function cross0_2d ( p0, p1, p2 )

!*****************************************************************************80
!
!! CROSS0_2D finds the cross product of (P1-P0) and (P2-P0) in 2D.
!
!  Discussion:
!
!    Strictly speaking, the vectors lie in the (X,Y) plane, and
!    the cross product here is a vector in the Z direction.
!
!    The vectors are specified with respect to a basis point P0.
!    We are computing the normal to the triangle (P0,P1,P2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision P0(2), P1(2), P2(2), the coordinates of
!    the three points.
!
!    Output, double precision CROSS0_2D, the Z component of the cross product
!    (P1-P0) x (P2-P0).
!
  implicit none

  double precision cross0_2d
  double precision p0(2)
  double precision p1(2)
  double precision p2(2)

  cross0_2d = ( p1(1) - p0(1) ) * ( p2(2) - p0(2) ) &
            - ( p1(2) - p0(2) ) * ( p2(1) - p0(1) )
end

subroutine ij_next ( i, j, n )

!*****************************************************************************80
!
!! IJ_NEXT returns the next matrix index.
!
!  Discussion:
!
!    For N = 3, the sequence of indices returned is:
!
!      (1,1), (1,2), (1,3), (2,1), (2,2), (2,3), (3,1), (3,2), (3,3), (0,0).
!
!    Note that once the value (N,N) is returned, the next value returned
!    will be (0,0).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer I, J.  On input, the current pair of
!    indices.  On output, the next pair of indices.  If either index is illegal
!    on input, the output value of (I,J) will be (1,1).
!
!    Input, integer N, the maximum value for I and J.
!
  implicit none

  integer i
  integer j
  integer n

  if ( n < 1 ) then
    i = 0
    j = 0
  end if

  if ( i < 1 .or. n < i .or. j < 1 .or. n < j ) then
    i = 1
    j = 1
  end if

  if ( j < n ) then
    j = j + 1
  else if ( i < n ) then
    i = i + 1
    j = 1
  else
    i = 0
    j = 0
  end if
end

subroutine ij_next_gt ( i, j, n )

!*****************************************************************************80
!
!! IJ_NEXT_GT returns the next matrix index, with the constraint that I < J.
!
!  Discussion:
!
!    For N = 5, the sequence of indices returned is:
!
!      (1,2), (1,3), (1,4), (1,5), (2,3), (2,4), (2,5), (3,4), (3,5), (4,5).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer I, J.  On input, the current pair of
!    indices.  On output, the next pair of indices.  If either index is illegal
!    on input, the output value of (I,J) will be (1,2).
!
!    Input, integer N, the maximum value for I and J.
!    A value of N less than 2 is nonsense.
!
  implicit none

  integer i
  integer j
  integer n

  if ( n < 2 ) then
    i = 0
    j = 0
  end if

  if ( i < 1 .or. n < i .or. j < 1 .or. n < j .or. j <= i ) then
    i = 1
    j = 2
  end if

  if ( j < n ) then
    j = j + 1
  else if ( i < n - 1 ) then
    i = i + 1
    j = i + 1
  else
    i = 0
    j = 0
  end if
end

subroutine line_exp2imp_2d ( x1, y1, x2, y2, a, b, c )

!*****************************************************************************80
!
!! LINE_EXP2IMP_2D converts an explicit line to implicit form in 2D.
!
!  Discussion:
!
!    The explicit form of a line in 2D is:
!
!      (X1,Y1), (X2,Y2).
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision X1, Y1, X2, Y2.  (X1,Y1) and (X2,Y2) are
!    two points on the line. (X1,Y1) must be different
!    from (X2,Y2).
!
!    Output, double precision A, B, C, three coefficients which describe
!    the line that passes through (X1,Y1) and (X2,Y2).
!
  implicit none

  double precision a
  double precision b
  double precision c
  double precision x1
  double precision x2
  double precision y1
  double precision y2
!
!  Take care of degenerate cases.
!
  if ( x1 == x2 .and. y1 == y2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE_EXP2IMP_2D - Fatal error!'
    write ( *, '(a)' ) '  (X1,Y1) = (X2,Y2)'
    write ( *, '(a,g14.6,2x,g14.6)' ) '  (X1,Y1) = ', x1, y1
    write ( *, '(a,g14.6,2x,g14.6)' ) '  (X2,Y2) = ', x2, y2
    stop
  end if

  a = y2 - y1
  b = x1 - x2
  c = x2 * y1 - x1 * y2
end

subroutine line_exp_point_dist_2d ( x1, y1, x2, y2, x, y, dist )

!*****************************************************************************80
!
!! LINE_EXP_POINT_DIST_2D: distance ( explicit line, point ) in 2D.
!
!  Discussion:
!
!    The explicit form of a line in 2D is:
!
!      (X1,Y1), (X2,Y2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision X1, Y1, X2, Y2.  (X1,Y1) and (X2,Y2) are two
!    points on the line.
!
!    Input, double precision X, Y, the point whose distance from the line is
!    to be measured.
!
!    Output, double precision DIST, the distance from the point to the line.
!
  implicit none

  double precision bot
  double precision dist
  double precision dot
  double precision t
  double precision x
  double precision xn
  double precision x1
  double precision x2
  double precision y
  double precision yn
  double precision y1
  double precision y2

  bot = ( x1 - x2 )**2 + ( y1 - y2 )**2

  if ( bot == 0.0D+00 ) then

    xn = x1
    yn = y1
!
!  (P-P1) dot (P2-P1) = Norm(P-P1) * Norm(P2-P1) * Cos(Theta).
!
!  (P-P1) dot (P2-P1) / Norm(P-P1)**2 = normalized coordinate T
!  of the projection of (P-P1) onto (P2-P1).
!
  else

    dot = ( x - x1 ) * ( x2 - x1 ) + ( y - y1 ) * ( y2 - y1 )

    t = dot / bot

    xn = x1 + t * ( x2 - x1 )
    yn = y1 + t * ( y2 - y1 )

  end if

  dist = sqrt ( ( xn - x )**2 + ( yn - y )**2 )
end

subroutine line_exp_point_dist_signed_2d ( x1, y1, x2, y2, x, y, dist_signed )

!*****************************************************************************80
!
!! LINE_EXP_POINT_DIST_SIGNED_2D: signed distance ( explicit line, point ).
!
!  Discussion:
!
!    The signed distance has two interesting properties:
!
!    *  The absolute value of the signed distance is the
!        usual (Euclidean) distance.
!
!    *  Points with signed distance 0 lie on the line,
!       points with a negative signed distance lie on one side
!         of the line,
!       points with a positive signed distance lie on the
!         other side of the line.
!
!    Assuming that C is nonnegative, then if a point is a positive
!    distance away from the line, it is on the same side of the
!    line as the point (0,0), and if it is a negative distance
!    from the line, it is on the opposite side from (0,0).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision X1, Y1, X2, Y2, define the two points
!    (X1,Y1) and (X2,Y2) that determine the line.
!
!    Input, double precision X, Y, the point (X,Y) whose signed distance
!    is desired.
!
!    Output, double precision DIST_SIGNED, the signed distance from the
!    point to the line.
!
  implicit none

  double precision a
  double precision b
  double precision c
  double precision dist_signed
  double precision x
  double precision x1
  double precision x2
  double precision y
  double precision y1
  double precision y2
!
!  Convert the line to implicit form.
!
  call line_exp2imp_2d ( x1, y1, x2, y2, a, b, c )
!
!  Compute the signed distance from the point to the line.
!
  dist_signed = ( a * x + b * y + c ) / sqrt ( a * a + b * b )
end

subroutine line_seg_contains_point_2d ( x1, y1, x2, y2, x3, y3, u, v )

!*****************************************************************************80
!
!! LINE_SEG_CONTAINS_POINT_2D reports if a line segment contains a point in 2D.
!
!  Discussion:
!
!    In exact arithmetic, point P3 = (X3,Y3) is on the line segment between
!    P1=(X1,Y1) and P2=(X2,Y2) if and only if 0 <= U <= 1 and V = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision X1, Y1, X2, Y2, the endpoints P1 and P2 of
!    a line segment.
!
!    Input, double precision X3, Y3, a point P3 to be tested.
!
!    Output, double precision U, the coordinate of (X3,Y3) along the axis from
!    with origin at P1 and unit at P2.
!
!    Output, double precision V, the magnitude of the off-axis portion of the
!    vector P3-P1, measured in units of (P2-P1).
!
  implicit none

  double precision u
  double precision unit
  double precision v
  double precision x1
  double precision x2
  double precision x3
  double precision y1
  double precision y2
  double precision y3

  unit = sqrt ( ( x2 - x1 )**2 + ( y2 - y1 )**2 )

  if ( unit == 0.0D+00 ) then

    if ( x3 == x1 .and. y3 == y1 ) then
      u = 0.5D+00
      v = 0.0D+00
    else
      u = 0.5D+00
      v = huge ( 1.0D+00 )
    end if

  else

    u = ( ( x3 - x1 ) * ( x2 - x1 ) + ( y3 - y1 ) * ( y2 - y1 ) ) / unit**2

    v = sqrt ( ( ( u - 1.0D+00 ) * x1 - u * x2 + x3 )**2 &
             + ( ( u - 1.0D+00 ) * y1 - u * y2 + y3 )**2 ) / unit

  end if
end

subroutine line_seg_vec_int_2d ( n, x1, y1, x2, y2, i, j, flag, xint, yint )

!*****************************************************************************80
!
!! LINE_SEG_VEC_INT_2D computes intersections of a set of line segments.
!
!  Discussion:
!
!    This is an implementation of the relatively inefficient algorithm
!    for computing all intersections of a set of line segments.
!
!    This is an "incremental" code, which returns as soon as it finds
!    a single intersection.  To find the next intersection, simply call
!    again.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of line segments.
!
!    Input, double precision X1(N), Y1(N), X2(N), Y2(N), the coordinates of the
!    endpoints of the line segments.
!
!    Input/output, integer I, J, used to keep track of the
!    computation.  On first call with a given set of line segments,
!    set I = J = 0.  On return with FLAG = 1, I and J record the indices of the
!    line segments whose intersection has just been found.  To find the
!    next intersection, simply call again, but do not alter I and J.
!
!    Output, integer FLAG:
!    0, no more intersections, the computation is done.
!    1, an intersection was detected between segments I and J.
!
!    Output, double precision XINT, YINT, the location of an intersection
!    of line segments I and J, if FLAG is 1.
!
  implicit none

  integer n

  integer flag
  integer i
  integer j
  double precision x1(n)
  double precision x2(n)
  double precision xint
  double precision y1(n)
  double precision y2(n)
  double precision yint

  do

    call ij_next_gt ( i, j, n )

    if ( i == 0 ) then
      flag = 0
      exit
    end if

    call lines_seg_int_2d ( x1(i), y1(i), x2(i), y2(i), x1(j), y1(j), &
      x2(j), y2(j), flag, xint, yint )

    if ( flag == 1 ) then
    end if

  end do
end

subroutine lines_exp_int_2d ( ival, x, y, x1, y1, x2, y2, x3, y3, x4, y4 )

!*****************************************************************************80
!
!! LINES_EXP_INT_2D determines where two explicit lines intersect in 2D.
!
!  Discussion:
!
!    The explicit form of a line in 2D is:
!
!      (X1,Y1), (X2,Y2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IVAL, reports on the intersection.
!
!     0, no intersection, the lines may be parallel or degenerate.
!     1, one intersection point, returned in X, Y.
!     2, infinitely many intersections, the lines are identical.
!
!    Output, double precision X, Y, if IVAl = 1, then X, Y contains
!    the intersection point.  Otherwise, X = 0, Y = 0.
!
!    Input, double precision X1, Y1, X2, Y2, define the first line.
!
!    Input, double precision X3, Y3, X4, Y4, define the second line.
!
  implicit none

  double precision a1
  double precision a2
  double precision b1
  double precision b2
  double precision c1
  double precision c2
  integer ival
  logical point_1
  logical point_2
  double precision x
  double precision x1
  double precision x2
  double precision x3
  double precision x4
  double precision y
  double precision y1
  double precision y2
  double precision y3
  double precision y4

  ival = 0
  x = 0.0D+00
  y = 0.0D+00
!
!  Check whether either line is a point.
!
  if ( x1 == x2 .and. y1 == y2 ) then
    point_1 = .true.
  else
    point_1 = .false.
  end if

  if ( x3 == x4 .and. y3 == y4 ) then
    point_2 = .true.
  else
    point_2 = .false.
  end if
!
!  Convert the lines to ABC format.
!
  if ( .not. point_1 ) then
    call line_exp2imp_2d ( x1, y1, x2, y2, a1, b1, c1 )
  end if

  if ( .not. point_2 ) then
    call line_exp2imp_2d ( x3, y3, x4, y4, a2, b2, c2 )
  end if
!
!  Search for intersection of the lines.
!
  if ( point_1 .and. point_2 ) then
    if ( x1 == x3 .and. y1 == y3 ) then
      ival = 1
      x = x1
      y = y1
    end if
  else if ( point_1 ) then
    if ( a2 * x1 + b2 * y1 == c2 ) then
      ival = 1
      x = x1
      y = y1
    end if
  else if ( point_2 ) then
    if ( a1 * x3 + b1 * y3 == c1 ) then
      ival = 1
      x = x3
      y = y3
    end if
  else
    call lines_imp_int_2d ( a1, b1, c1, a2, b2, c2, ival, x, y )
  end if
end

subroutine lines_imp_int_2d ( a1, b1, c1, a2, b2, c2, ival, x, y )

!*****************************************************************************80
!
!! LINES_IMP_INT_2D determines where two implicit lines intersect in 2D.
!
!  Discussion:
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision A1, B1, C1, define the first line.
!    At least one of A1 and B1 must be nonzero.
!
!    Input, double precision A2, B2, C2, define the second line.
!    At least one of A2 and B2 must be nonzero.
!
!    Output, integer IVAL, reports on the intersection.
!
!    -1, both A1 and B1 were zero.
!    -2, both A2 and B2 were zero.
!     0, no intersection, the lines are parallel.
!     1, one intersection point, returned in X, Y.
!     2, infinitely many intersections, the lines are identical.
!
!    Output, double precision X, Y, if IVAL = 1, then X, Y contains
!    the intersection point.  Otherwise, X = 0, Y = 0.
!
  implicit none

  double precision a(2,2)
  double precision a1
  double precision a2
  double precision b(2,2)
  double precision b1
  double precision b2
  double precision c1
  double precision c2
  double precision det
  integer ival
  double precision x
  double precision y

  x = 0.0D+00
  y = 0.0D+00
!
!  Refuse to handle degenerate lines.
!
  if ( a1 == 0.0D+00 .and. b1 == 0.0D+00 ) then
    ival = -1
    return
  else if ( a2 == 0.0D+00 .and. b2 == 0.0D+00 ) then
    ival = -2
  end if
!
!  Set up a linear system, and compute its inverse.
!
  a(1,1) = a1
  a(1,2) = b1
  a(2,1) = a2
  a(2,2) = b2

  call r8mat2_inverse ( a, b, det )
!
!  If the inverse exists, then the lines intersect.
!  Multiply the inverse times -C to get the intersection point.
!
  if ( det /= 0.0D+00 ) then

    ival = 1
    x = - b(1,1) * c1 - b(1,2) * c2
    y = - b(2,1) * c1 - b(2,2) * c2
!
!  If the inverse does not exist, then the lines are parallel
!  or coincident.  Check for parallelism by seeing if the
!  C entries are in the same ratio as the A or B entries.
!
  else

    ival = 0

    if ( a1 == 0.0D+00 ) then
      if ( b2 * c1 == c2 * b1 ) then
        ival = 2
      end if
    else
      if ( a2 * c1 == c2 * a1 ) then
        ival = 2
      end if
    end if

  end if
end

subroutine lines_seg_dist_2d ( x1, y1, x2, y2, x3, y3, x4, y4, flag, x5, y5 )

!*****************************************************************************80
!
!! LINES_SEG_DIST_2D computes the distance of two line segments in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision X1, Y1, X2, Y2, the endpoints of the first
!    segment.
!
!    Input, double precision X3, Y3, X4, Y4, the endpoints of the second
!    segment.
!
!    Output, integer FLAG, records the results.
!    0, the line segments do not intersect.
!    1, the line segments intersect.
!
!    Output, double precision X5, Y5.
!    If FLAG = 0, X5 = Y5 = 0.
!    If FLAG = 1, then (X5,Y5) is a point of intersection.
!
  implicit none

  integer flag
  integer ival
  double precision u
  double precision v
  double precision x1
  double precision x2
  double precision x3
  double precision x4
  double precision x5
  double precision y1
  double precision y2
  double precision y3
  double precision y4
  double precision y5

  x5 = 0.0D+00
  y5 = 0.0D+00
!
!  Find the intersection of the two lines.
!
  call lines_exp_int_2d ( ival, x5, y5, x1, y1, x2, y2, x3, y3, x4, y4 )

  if ( ival == 0 ) then
    flag = 0
  end if
!
!  Is the point on the first segment?
!
  call line_seg_contains_point_2d ( x1, y1, x2, y2, x5, y5, u, v )

  if ( u < 0 .or. 1.0D+00 < u .or. 0.001D+00 < v ) then
    flag = 0
  end if
!
!  Is the point on the second segment?
!
  call line_seg_contains_point_2d ( x3, y3, x4, y4, x5, y5, u, v )

  if ( u < 0 .or. 1.0D+00 < u .or. 0.001D+00 < v ) then
    flag = 0
  end if

  flag = 1
end

subroutine lines_seg_int_1d ( x1, x2, x3, x4, flag, x5, x6 )

!*****************************************************************************80
!
!! LINES_SEG_INT_1D computes the intersection of two line segments in 1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision X1, X2, the endpoints of the first segment.
!
!    Input, double precision X3, X4, the endpoints of the second segment.
!
!    Output, integer FLAG, records the results.
!    0, the line segments do not intersect.
!    1, the line segments intersect.
!
!    Output, double precision X5, X6, the endpoints of the intersection
!    segment.  If FLAG = 0, X5 = X6 = 0.
!
  implicit none

  integer flag
  double precision x1
  double precision x2
  double precision x3
  double precision x4
  double precision x5
  double precision x6
  double precision y1
  double precision y2
  double precision y3
  double precision y4

  y1 = min ( x1, x2 )
  y2 = max ( x1, x2 )
  y3 = min ( x3, x4 )
  y4 = max ( x3, x4 )

  flag = 0
  x5 = 0.0D+00
  x6 = 0.0D+00

  if ( y4 < y1 ) then
    return
  else if ( y2 < y3 ) then
  end if

  flag = 1
  x5 = max ( y1, y3 )
  x6 = min ( y2, y4 )
end

subroutine lines_seg_int_2d ( x1, y1, x2, y2, x3, y3, x4, y4, flag, x5, y5 )

!*****************************************************************************80
!
!! LINES_SEG_INT_2D computes the intersection of two line segments in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision X1, Y1, X2, Y2, the endpoints of the first
!    segment.
!
!    Input, double precision X3, Y3, X4, Y4, the endpoints of the second
!    segment.
!
!    Output, integer FLAG, records the results.
!    0, the line segments do not intersect.
!    1, the line segments intersect.
!
!    Output, double precision X5, Y5.
!    If FLAG = 0, X5 = Y5 = 0.
!    If FLAG = 1, then (X5,Y5) is a point of intersection.
!
  implicit none

  integer flag
  integer ival
  double precision u
  double precision v
  double precision x1
  double precision x2
  double precision x3
  double precision x4
  double precision x5
  double precision y1
  double precision y2
  double precision y3
  double precision y4
  double precision y5

  x5 = 0.0D+00
  y5 = 0.0D+00
!
!  Find the intersection of the two lines.
!
  call lines_exp_int_2d ( ival, x5, y5, x1, y1, x2, y2, x3, y3, x4, y4 )

  if ( ival == 0 ) then
    flag = 0
  end if
!
!  Is the point on the first segment?
!
  call line_seg_contains_point_2d ( x1, y1, x2, y2, x5, y5, u, v )

  if ( u < 0 .or. 1.0D+00 < u .or. 0.001D+00 < v ) then
    flag = 0
  end if
!
!  Is the point on the second segment?
!
  call line_seg_contains_point_2d ( x3, y3, x4, y4, x5, y5, u, v )

  if ( u < 0 .or. 1.0D+00 < u .or. 0.001D+00 < v ) then
    flag = 0
  end if

  flag = 1
end

subroutine points_convex_hull_cubic_2d ( node_num, node_xy, hull_num, hull )

!*****************************************************************************80
!
!! POINTS_CONVEX_HULL_CUBIC_2D computes the convex hull of 2D points.
!
!  Discussion:
!
!    The algorithm used requires time that is cubic in the number of input
!    points.  Algorithms are available which are much more efficient!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marc de Berg, Marc van Kreveld, Mark Overmars, Otfried Schwarzkopf,
!    Computational Geometry,
!    Second Edition,
!    Springer, 2000, pages 2-8.
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of points.
!
!    Input, double precision NODE_XY(2,NODE_NUM), the coordinates of
!    the points.
!
!    Output, integer HULL_NUM, the number of vertices in the
!    convex hull.
!
!    Output, integer HULL(NODE_NUM), the HULL_NUM vertices that
!    form the convex hull, in counter clockwise order.
!
  implicit none

  integer node_num

  integer b(node_num)
  double precision cross
  double precision cross0_2d
  integer e(node_num)
  integer hull(node_num)
  integer hull_num
  integer i
  integer j
  integer k
  integer match
  double precision node_xy(2,node_num)
  logical valid

  hull_num = 0

  do i = 1, node_num

    do j = 1,node_num
      if ( i /= j ) then
        valid = .true.
        do k = 1, node_num
          if ( k /= i .and. k /= j ) then
!
!  Compute the cross product: P(J)-P(I) x P(K)-P(I)
!  to determine if P(K) is strictly to the left of the line from P(I) to P(J).
!
!  I'll have to think if CROSS should be positive or negative...
!
            cross = cross0_2d ( node_xy(1:2,i), node_xy(1:2,j), node_xy(1:2,k) )

            if ( 0.0D+00 < cross ) then
              valid = .false.
              exit
            end if

          end if

        end do
!
!  Add the line from P(I) to P(J) to the list.
!
        if ( valid ) then
          hull_num = hull_num + 1
          b(hull_num) = i
          e(hull_num) = j
        end if

      end if

    end do

  end do
!
!  From the unordered edge list of vertex pairs,
!  construct an ordered list of vertices.
!
  hull(1) = b(1)
  match = e(1)

  do k = 2, hull_num

    hull(k) = 0

    do j = k, hull_num
      if ( b(j) == match ) then
        hull(k) = b(j)
        match = e(j)
        call i4_swap ( b(j), b(k) )
        call i4_swap ( e(j), e(k) )
        exit
      end if
    end do

    if ( hull(k) == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POINTS_2D_CONVEX_HULL_CUBIC - Fatal error!'
      write ( *, '(a,i8)' ) '  Algorithm failed for K = ', k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  B, E arrays:'
      write ( *, '(a)' ) ' '
      do i = 1, hull_num
        write ( *, '(2x,i8,2x,i8)' ) b(i), e(i)
      end do
      stop
    end if

  end do

  if ( match /= hull(1) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POINTS_2D_CONVEX_HULL_CUBIC - Fatal error!'
    write ( *, '(a)' ) '  Algorithm failed to link last node to first.'
    stop
  end if
!
!  Reverse the order.
!
  hull(1:hull_num) = hull(hull_num:1:-1)
end

subroutine points_convex_hull_nlogh_2d ( node_num, node_xy, hull_num, hull )

!*****************************************************************************80
!
!! POINTS_CONVEX_HULL_NLOGH_2D computes the convex hull of 2D points.
!
!  Discussion:
!
!    The work involved is N*log(H), where N is the number of points, and H is
!    the number of points that are on the hull.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, double precision NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Output, integer HULL_NUM, the number of nodes that lie
!    on the convex hull.
!
!    Output, integer HULL(NODE_NUM).  Entries 1 through HULL_NUM
!    contain the indices of the nodes that form the convex hull, in order.
!
  implicit none

  integer node_num

  double precision angle
  double precision angle_max
  double precision angle_rad_2d
  double precision di
  double precision dr
  integer first
  integer hull(node_num)
  integer hull_num
  integer i
  double precision node_xy(2,node_num)
  double precision p_xy(2)
  integer q
  double precision q_xy(2)
  integer r
  double precision r_xy(2)

  if ( node_num < 1 ) then
    hull_num = 0
  end if
!
!  If NODE_NUM = 1, the hull is the point.
!
  if ( node_num == 1 ) then
    hull_num = 1
    hull(1) = 1
  end if
!
!  If NODE_NUM = 2, then the convex hull is either the two distinct points,
!  or possibly a single (repeated) point.
!
  if ( node_num == 2 ) then

    if ( node_xy(1,1) /= node_xy(1,2) .or. node_xy(2,1) /= node_xy(2,2) ) then
      hull_num = 2
      hull(1) = 1
      hull(2) = 2
    else
      hull_num = 1
      hull(1) = 1
    end if
  end if
!
!  Find the leftmost point and call it "Q".
!  In case of ties, take the bottom-most.
!
  q = 1
  do i = 2, node_num
    if ( node_xy(1,i) < node_xy(1,q) .or. &
       ( node_xy(1,i) == node_xy(1,q) .and. &
         node_xy(2,i) < node_xy(2,q) ) ) then
      q = i
    end if
  end do

  q_xy(1:2) = node_xy(1:2,q)
!
!  Remember the starting point, so we know when to stop!
!
  first = q
  hull_num = 1
  hull(1) = q
!
!  For the first point, make a dummy previous point, 1 unit south,
!  and call it "P".
!
  p_xy(1) = q_xy(1)
  p_xy(2) = q_xy(2) - 1.0D+00
!
!  Now, having old point P, and current point Q, find the new point R
!  so the angle PQR is maximal.
!
!  Watch out for the possibility that the two nodes are identical.
!
  do

    r = 0
    angle_max = 0.0D+00

    do i = 1, node_num

      if ( i /= q .and. &
           ( node_xy(1,i) /= q_xy(1) .or. node_xy(2,i) /= q_xy(2) ) ) then

        angle = angle_rad_2d ( p_xy, q_xy, node_xy(1:2,i) )

        if ( r == 0 .or. angle_max < angle ) then

          r = i
          r_xy(1:2) = node_xy(1:2,r)
          angle_max = angle
!
!  In case of ties, choose the nearer point.
!
        else if ( r /= 0 .and. angle == angle_max ) then

          di = ( node_xy(1,i) - q_xy(1) )**2 + ( node_xy(2,i) - q_xy(2) )**2
          dr = ( r_xy(1)      - q_xy(1) )**2 + ( r_xy(2)      - q_xy(2) )**2

          if ( di < dr ) then
            r = i
            r_xy(1:2) = node_xy(1:2,r)
            angle_max = angle
          end if

        end if

      end if

    end do
!
!  We are done when we have returned to the first point on the convex hull.
!
    if ( r == first ) then
      exit
    end if

    hull_num = hull_num + 1

    if ( node_num < hull_num ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POINTS_CONVEX_HULL_NLOGH_2D - Fatal error!'
      write ( *, '(a)' ) '  The algorithm has failed.'
      stop
    end if
!
!  Add point R to convex hull.
!
    hull(hull_num) = r
!
!  Set P := Q, Q := R, and prepare to search for next point R.
!
    q = r

    p_xy(1:2) = q_xy(1:2)
    q_xy(1:2) = r_xy(1:2)

  end do
!
!  Reverse the order of the points.
!
  hull(1:hull_num) = hull(hull_num:1:-1)
end

subroutine points_convex_hull_nlogn_2d ( node_num, node_xy, hull_num, hull )

!*****************************************************************************80
!
!! POINTS_CONVEX_HULL_NLOGN_2D computes the convex hull of 2D points.
!
!  Discussion:
!
!    The algorithm used requires time that is of order N * Log ( N ) in
!    N. the number of input points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marc de Berg, Marc van Kreveld, Mark Overmars, Otfried Schwarzkopf,
!    Computational Geometry,
!    Second Edition,
!    Springer, 2000, pages 2-8.
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of points.
!
!    Input, double precision NODE_XY(2,NODE_NUM), the coordinates of
!    the points.
!
!    Output, integer HULL_NUM, the number of vertices in the
!    convex hull.
!
!    Output, integer HULL(NODE_NUM), the HULL_NUM vertices that
!    form the convex hull, in counter clockwise order.
!
  implicit none

  integer node_num

  double precision angle_rad_2d
  integer hull(node_num)
  integer hull_num
  integer i
  integer l(node_num+100)
  integer n1
  integer n2
  integer n3
  double precision node_xy(2,node_num)
  double precision p1(2)
  double precision p2(2)
  double precision p3(2)
  double precision , parameter :: pi = 3.14159265D+00
!
!  Sort the (X,Y) points lexicographically.
!
  call r82vec_sort_quick_a ( node_num, node_xy )
!
!  Start with L1 = ( 1, 2 )
!
!  List L1 is stored in indices N1 through N2 of L.
!
  n1 = 1
  n2 = 2

  l(1) = 1
  p2(1:2) = node_xy(1:2,1)
  l(2) = 2
  p3(1:2) = node_xy(1:2,2)
!
!  Append I to L1.
!    Then, while L1 contains at least 3 points
!     if the last 3 points do not make a strict right turn, delete middle one.
!
  do i = 3, node_num

    n2 = n2 + 1
    l(n2) = i
    p1(1:2) = p2(1:2)
    p2(1:2) = p3(1:2)
    p3(1:2) = node_xy(1:2,i)

    do while ( 3 <= n2 + 1 - n1 )

      if ( angle_rad_2d ( p1, p2, p3 ) < pi ) then
        exit
      end if

      l(n2-1) = l(n2)
      n2 = n2 - 1
      p2(1:2) = p1(1:2)

      if ( 1 <= n2 - 2 ) then
        p1(1:2) = node_xy(1:2,l(n2-2))
      end if

    end do

  end do
!
!  Start with L2 = ( N, N-1 )
!
!  List L2 is stored in indices N2 through N3 of L.
!
  n3 = n2
  n3 = n3 + 1
  p2(1:2) = node_xy(1:2,node_num)
  l(n3) = node_num - 1
  p3(1:2) = node_xy(1:2,node_num-1)
!
!  Append I to L2.
!    Then, while L2 contains at least 3 points,
!    if the last 3 points do not make a strict right turn, delete middle one.
!
  do i = node_num-2, 1, -1

    n3 = n3 + 1
    l(n3) = i
    p1(1:2) = p2(1:2)
    p2(1:2) = p3(1:2)
    p3(1:2) = node_xy(1:2,i)

    do while ( 3 <= n3 + 1 - n2 )

      if ( angle_rad_2d ( p1, p2, p3 ) < pi ) then
        exit
      end if

      l(n3-1) = l(n3)
      n3 = n3 - 1
      p2(1:2) = p1(1:2)

      if ( n2 <= n3 - 2 ) then
        p1(1:2) = node_xy(1:2,l(n3-2))
      end if

    end do

  end do
!
!  Last entry of L is 1, so knock it off.
!
  hull_num = n3 - 1
!
!  Reverse the order.
!
  hull(1:hull_num) = l(hull_num:1:-1)
end

subroutine points_minidisc1_2d ( n, px, py, qx, qy, r, cx, cy )

!*****************************************************************************80
!
!! POINTS_MINIDISC1_2D finds the smallest circle through Q containing points P.
!
!  Discussion:
!
!    It is assumed that there IS at least one disk which has the
!    points Q on its boundary, and which contains all the
!    points P.  For arbitrary points, this need not be the case.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marc de Berg, Marc van Kreveld, Mark Overmars, Otfried Schwarzkopf,
!    Computational Geometry,
!    Second Edition,
!    Springer, 2000, pages 86-90.
!
!  Parameters:
!
!    Input, integer N, the number of points in the set P.
!
!    Input, double precision PX(N), PY(N), the X and Y coordinates of a set of
!    points in the plane.
!
!    Input, double precision QX, QY, a point in the plane.
!
!    Output, double precision R, CX, CY, the radius and center of the smallest
!    disk that encloses P and has Q on its boundary.
!
  implicit none

  integer n

  double precision cx
  double precision cy
  logical inside
  integer j
  double precision px(n)
  double precision py(n)
  double precision qx
  double precision qy
  double precision r
!
!  Determine the smallest disk with Q and P(1) on its boundary,
!  which is simply the circle whose diameter is the segment
!  from Q to P(1).
!
  call circle_dia2imp_2d ( qx, qy, px(1), py(1), r, cx, cy )
!
!  Now consider a point in the set P.  If it is not already
!  contained in the current circle, then expand the circle.
!
  do j = 2, n

    call circle_imp_contains_point_2d ( r, cx, cy, px(j), py(j), inside )

    if ( .not. inside ) then
      call points_minidisc2_2d ( j-1, px, py, px(j), py(j), qx, qy, r, cx, cy )
    end if

  end do
end

subroutine points_minidisc2_2d ( n, px, py, q1x, q1y, q2x, q2y, r, cx, cy )

!*****************************************************************************80
!
!! POINTS_MINIDISC2_2D: smallest circle through Q1 and Q2 containing points P.
!
!  Discussion:
!
!    It is assumed that there IS at least one disk which has the
!    points Q1 and Q2 on its boundary, and which contains all the
!    points P.  For arbitrary points, this need not be the case.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marc de Berg, Marc van Kreveld, Mark Overmars, Otfried Schwarzkopf,
!    Computational Geometry,
!    Second Edition,
!    Springer, 2000, pages 86-90.
!
!  Parameters:
!
!    Input, integer N, the number of points in the set P.
!
!    Input, double precision PX(N), PY(N), the X and Y coordinates of a set of
!    points in the plane.
!
!    Input, double precision Q1X, Q1Y, Q2X, Q2Y, two points in the plane.
!
!    Output, double precision R, CX, CY, the radius and center of the smallest
!    disk that encloses P and has Q1 and Q2 on its boundary.
!
  implicit none

  integer n

  double precision cx
  double precision cy
  logical inside
  integer k
  double precision px(n)
  double precision py(n)
  double precision q1x
  double precision q1y
  double precision q2x
  double precision q2y
  double precision r
!
!  Determine the smallest disk with Q1, Q2 on its boundary,
!  which is simply the circle whose diameter is the segment
!  from Q1 to Q2.
!
  call circle_dia2imp_2d ( q1x, q1y, q2x, q2y, r, cx, cy )
!
!  Now consider a point in the set P.  If it is not already
!  contained in the current circle, then expand the circle.
!
  do k = 1, n

    call circle_imp_contains_point_2d ( r, cx, cy, px(k), py(k), inside )

    if ( .not. inside ) then
      call circle_exp2imp_2d ( q1x, q1y, q2x, q2y, px(k), py(k), r, cx, cy )
    end if

  end do
end

subroutine points_minidisc_2d ( n, px, py, r, cx, cy )

!*****************************************************************************80
!
!! POINTS_MINIDISC_2D finds the smallest circle containing points P.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marc de Berg, Marc van Kreveld, Mark Overmars, Otfried Schwarzkopf,
!    Computational Geometry,
!    Second Edition,
!    Springer, 2000, pages 86-90.
!
!  Parameters:
!
!    Input, integer N, the number of points in the set P.
!
!    Input, double precision PX(N), PY(N), the X and Y coordinates of a set of
!    points in the plane.
!
!    Output, double precision R, CX, CY, the radius and center of the smallest
!    disk that encloses P.
!
  implicit none

  integer n

  double precision cx
  double precision cy
  integer i
  logical inside
  double precision px(n)
  double precision py(n)
  double precision r
!
!  N = 1
!
  if ( n == 1 ) then
    r = 0.0D+00
    cx = px(1)
    cy = py(1)
  end if
!
!  N = 2
!
  if ( n == 2 ) then
    r = 0.5D+00 * sqrt ( ( px(1) - px(2) )**2 + ( py(1) - py(2) )**2 )
    cx = 0.5D+00 * ( px(1) + px(2) )
    cy = 0.5D+00 * ( py(1) + py(2) )
  end if
!
!  Determine the smallest disk with P(1), P(2) on its boundary,
!
  call circle_dia2imp_2d ( px(1), py(1), px(2), py(2), r, cx, cy )
!
!  Now consider a point in the set P.  If it is not already
!  contained in the current circle, then expand the circle.
!
  do i = 3, n

    call circle_imp_contains_point_2d ( r, cx, cy, px(i), py(i), inside )

    if ( .not. inside ) then
      call points_minidisc1_2d ( i-1, px, py, px(i), py(i), r, cx, cy )
    end if

  end do
end

subroutine poly_triangulate_2d ( n, x, y, triang )

!*****************************************************************************80
!
!! POLY_TRIANGULATE_2D returns a triangulation of a polygon.
!
!  Discussion:
!
!    Given a polygon with N (distinct) nodes, a triangulation consists of
!    N-2 (distinct) triangles, with the property that each triangle is
!    constructed from vertices of the polygon, and that every point in the
!    polygon is contained in exactly one triangle, unless it lies on the
!    boundary of two or more triangles, and that no point exterior to the
!    polygon is contained in any triangle.  In other words, a triangulation
!    is a dissection of the area of the polygon into triangles.
!
!    This algorithm should be correct, but not efficient.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marc de Berg, Marc van Kreveld, Mark Overmars, Otfried Schwarzkopf,
!    Computational Geometry,
!    Second Edition,
!    Springer, 2000, pages 46-49.
!
!  Parameters:
!
!    Input, integer N, the number of nodes in the polygon.
!
!    Input, double precision X(N), Y(N), the coordinates of the nodes,
!    listed in counter-clockwise order.
!
!    Output, integer TRIANG(3,N-2), the triangulation of
!    the polygon.
!
  implicit none

  integer n
  integer , parameter :: stack2_max = 100

  integer best
  integer degree
  integer degree2
  double precision dist
  double precision dist_max
  logical inside
  integer j
  integer number
  integer poly(n)
  integer stack1(n-2)
  integer stack1_num
  integer stack2(stack2_max)
  integer stack2_num
  integer t
  integer triang(3,n-2)
  integer triang_num
  integer u
  integer v
  integer w
  double precision x(n)
  double precision x1
  double precision x2
  double precision x3
  double precision x4
  double precision y(n)
  double precision y1
  double precision y2
  double precision y3
  double precision y4

  if ( n <= 2 ) then
  end if

  degree = n
  call i4vec_indicator ( n, poly )

  call poly_reorder_nodes ( n, x, y, degree, poly )

  stack1_num = 0
  stack2_num = 0
  triang_num = 0

  do

    if ( degree == 3 ) then

      triang_num = triang_num + 1
      triang(1:3,triang_num) = poly(1:3)

      if ( stack1_num <= 0 ) then
        exit
      end if

      call i4vec_pop ( degree, poly, n-2, stack1_num, stack1, &
        stack2_max, stack2_num, stack2 )

      call poly_reorder_nodes ( n, x, y, degree, poly )

      cycle

    end if

    u = poly(degree)
    v = poly(1)
    w = poly(2)

    x1 = x(u)
    y1 = y(u)
    x2 = x(v)
    y2 = y(v)
    x3 = x(w)
    y3 = y(w)

    dist_max = 0.0D+00
    best = 0
    number = 0

    do j = 3, degree - 1

      t = j
      x4 = x(poly(t))
      y4 = y(poly(t))

      call triangle_contains_point_2d ( x1, y1, x2, y2, x3, y3, x4, y4, inside )

      if ( inside ) then

        number = number + 1
        call line_exp_point_dist_2d ( x1, y1, x3, y3, x4, y4, dist )

        if ( dist_max <= dist ) then
          dist_max = dist
          best = t
        end if

      end if

    end do
!
!  If there were no polygonal nodes inside the triangle (U,V,W),
!  then triangle (U,V,W) can be used as part of the triangulation,
!  and the remaining region is simply the polygon (W,...,U).
!
    if ( number == 0 ) then

      triang_num = triang_num + 1
      triang(1,triang_num) = poly(1)
      triang(2,triang_num) = poly(2)
      triang(3,triang_num) = poly(degree)

      poly(1:n-1) = poly(2:n)
      degree = degree - 1

      call poly_reorder_nodes ( n, x, y, degree, poly )

    else
!
!  If there were polygonal nodes inside the triangle (U,V,W),
!  then POLY(T) is the node that was furthest from the line (U,W).
!
!  Split the polygon (V,W,part1,T,part2,U)
!  into (V,W,part1,T) and (V,T,part2,U)
!
      degree2 = t

      call i4vec_push ( degree2, poly(1:t), n-2, stack1_num, stack1, &
        stack2_max, stack2_num, stack2 )

      degree2 = degree + 2 - t
      poly(2:degree2) = poly(t:degree)
      degree = degree2

      call poly_reorder_nodes ( n, x, y, degree, poly )

    end if

  end do
end

subroutine poly_reorder_nodes ( nxy, x, y, npoly, poly )

!*****************************************************************************80
!
!! POLY_REORDER_NODES reorders nodes of a polygon so node 1 is leftest lowest.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NXY, the number of nodes.
!
!    Input, double precision X(NXY), Y(NXY), the coordinates of the nodes.
!
!    Input, integer NPOLY, the number of nodes of the polygon.
!
!    Input/output, POLY(NPOLY), the indices of the nodes.
!
  implicit none

  integer npoly
  integer nxy

  integer i
  integer imin
  integer p
  integer pmin
  integer poly(npoly)
  integer poly2(npoly)
  double precision x(nxy)
  double precision y(nxy)

  imin = 1
  pmin = poly(imin)

  do i = 2, npoly
    p = poly(i)
    if ( &
      ( x(p) < x(pmin) ) .or. &
      ( x(p) == x(pmin) .and. y(p) < y(pmin) ) ) then
      imin = i
      pmin = p
    end if
  end do

  if ( imin /= 1 ) then

    poly2(1:npoly+1-imin) = poly(imin:npoly)
    poly2(npoly+2-imin:npoly) = poly(1:imin-1)

    poly(1:npoly) = poly2(1:npoly)

  end if
end

subroutine polycon_minkowski_sum_linear ( nu, ux, uy, nv, vx, vy, nw, wx, wy )

!*****************************************************************************80
!
!! POLYCON_MINKOWSKI_SUM_LINEAR: the Minkowski sum of two convex polygons.
!
!  Discussion:
!
!    For two geometric shapes U and V, the Minkowski sum W is the
!    set of all points
!      w = u + v
!    formed by adding points from the two shapes.
!
!    The Minkowski sum of two convex polygons is also a convex polygon.
!
!    The algorithm used here is only valid for convex polygons.
!
!    The vertices must be listed in counterclockwise order, with
!    the first vertex having the smallest Y coordinate.  If two
!    or more vertices have the same Y coordinate, then the one with
!    smallest X coordinate should be first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marc de Berg, Marc van Kreveld, Mark Overmars, Otfried Schwarzkopf,
!    Computational Geometry,
!    Second Edition,
!    Springer, 2000, pages 275-281.
!
!  Parameters:
!
!    Input, integer NU, the number of vertices of the 
!    first polygon.
!
!    Input, double precision UX(NU), UY(NU), the coordinates of the vertices
!    of the first polygon.
!
!    Input, integer NV, the number of vertices of the 
!    second polygon.
!
!    Input, double precision VX(NV), VY(NY), the coordinates of the vertices
!    of the second polygon.
!
!    Output, integer NW, the number of vertices of the sum polygon.
!    NW will be at most NV+NU.
!
!    Output, double precision WX(*), WY(*), the coordinates of the vertices
!    of the sum polygon.
!
  implicit none

  integer nu
  integer nv

  double precision angle_rad_2d
  integer i
  integer i4_wrap
  integer iw
  integer ip1w
  integer j
  integer jp1w
  integer jw
  integer nw
  double precision p1(2)
  double precision p2(2)
  double precision p3(2)
  double precision q1(2)
  double precision q2(2)
  double precision q3(2)
  double precision ux(nu)
  double precision uy(nu)
  double precision u_angle
  logical u_done
  double precision vx(nv)
  double precision vy(nv)
  double precision v_angle
  logical v_done
  double precision wx(nu+nv+2)
  double precision wy(nu+nv+2)

  i = 1
  iw = 1
  u_done = .false.
  j = 1
  jw = 1
  v_done = .false.
  nw = 0

  do

    nw = nw + 1

    wx(nw) = ux(iw) + vx(jw)
    wy(nw) = uy(iw) + vy(jw)

    ip1w = i4_wrap ( i+1, 1, nu )

    p1(1) = ux(iw) + 1.0D+00
    p1(2) = uy(iw)
    p2(1) = ux(iw)
    p2(2) = uy(iw)
    p3(1) = ux(ip1w)
    p3(2) = uy(ip1w)

    u_angle = angle_rad_2d ( p1, p2, p3 )
    u_angle = u_angle + real ( i - nu) / real ( nu) &
      * 2.0D+00 * 3.141592653589793D+00

    jp1w = i4_wrap ( j+1, 1, nv )

    q1(1) = vx(jw) + 1.0D+00
    q1(2) = vy(jw)
    q2(1) = vx(jw)
    q2(2) = vy(jw)
    q3(1) = vx(jp1w)
    q3(2) = vy(jp1w)

    v_angle = angle_rad_2d ( q1, q2, q3 )
    v_angle = v_angle + real ( j - nv) / real ( nv) &
      * 2.0D+00 * 3.141592653589793D+00

    if ( u_angle < v_angle ) then
      i = i + 1
    else if ( v_angle < u_angle ) then
      j = j + 1
    else
      i = i + 1
      j = j + 1
    end if

    if ( i == nu + 1 ) then
      u_done = .true.
    end if

    if ( j == nv + 1 ) then
      v_done = .true.
    end if

    if ( u_done .and. v_done ) then
      exit
    end if

    iw = i4_wrap ( i, 1, nu )
    jw = i4_wrap ( j, 1, nv )

  end do
end

subroutine polycon_minkowski_sum_n2logn2 ( nu, ux, uy, nv, vx, vy, nw, w_xy )

!*****************************************************************************80
!
!! POLYCON_MINKOWSKI_SUM_N2LOGN2 Minkowski sums two convex polygons.
!
!  Discussion:
!
!    For two geometric shapes U and V, the Minkowski sum W is the
!    set of all points
!      w = u + v
!    formed by adding points from the two shapes.
!
!    The Minkowski sum of two convex polygons is also a convex polygon.
!
!    The algorithm used here is only valid for convex polygons.
!
!    The vertices must be listed in counterclockwise order, with
!    the first vertex having the smallest Y coordinate.  If two
!    or more vertices have the same Y coordinate, then the one with
!    smallest X coordinate should be first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marc de Berg, Marc van Kreveld, Mark Overmars, Otfried Schwarzkopf,
!    Computational Geometry,
!    Second Edition,
!    Springer, 2000, pages 275-281.
!
!  Parameters:
!
!    Input, integer NU, the number of vertices of the
!    first polygon.
!
!    Input, double precision UX(NU), UY(NU), the coordinates of the vertices
!    of the first polygon.
!
!    Input, integer NV, the number of vertices of the
!    second polygon.
!
!    Input, double precision VX(NV), VY(NY), the coordinates of the vertices
!    of the second polygon.
!
!    Output, integer NW, the number of vertices of the sum polygon.
!    NW will be at most NU+NV.
!
!    Output, double precision W_XY(2,*), the coordinates of the vertices
!    of the sum polygon.
!
  implicit none

  integer nu
  integer nv

  integer i
  integer j
  integer nuv
  integer nw
  double precision uv_xy(2,nu*nv)
  double precision ux(nu)
  double precision uy(nu)
  double precision vx(nv)
  double precision vy(nv)
  integer w(nu+nv)
  double precision w_xy(2,nu+nv)
!
!  Generate points from all pairs.
!
  nuv = 0
  do i = 1, nu
    do j = 1, nv
      nuv = nuv + 1
      uv_xy(1,nuv) = ux(i) + vx(j)
      uv_xy(2,nuv) = uy(i) + vy(j)
    end do
  end do
!
!  Compute the convex hull.
!  The output value of NW should be no more than NU+NV.
!
  call points_convex_hull_nlogn_2d ( nuv, uv_xy, nw, w )
!
!  Collect the points.
!
  w_xy(1:2,1:nw) = uv_xy(1:2,w(1:nw))
end

function r4_uniform_01 ( seed )

!*****************************************************************************80
!
!! R4_UNIFORM_01 returns a unit pseudorandom R4.
!
!  Discussion:
!
!    An R4 is a real value.
!
!    This routine implements the recursion
!
!      seed = ( 16807 * seed )  mod ( 2^31 - 1 )
!      r4_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R4_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre LEcuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real R4_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer , parameter :: i4_huge = 2147483647
  integer k
  integer seed
  real r4_uniform_01

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r4_uniform_01 = real ( seed) * 4.656612875E-10
end

subroutine r8mat2_inverse ( a, b, det )

!*****************************************************************************80
!
!! R8MAT2_INVERSE inverts a 2 by 2 real matrix using Cramer's rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision A(2,2), the matrix to be inverted.
!
!    Output, double precision B(2,2), the inverse of the matrix A.
!
!    Output, double precision DET, the determinant of the matrix A.
!
!    If DET is zero, then A is singular, and does not have an
!    inverse.  In that case, B is simply set to zero, and a
!    message is printed.
!
!    If DET is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
  implicit none

  double precision a(2,2)
  double precision b(2,2)
  double precision det
!
!  Compute the determinant.
!
  det = a(1,1) * a(2,2) - a(1,2) * a(2,1)
!
!  If the determinant is zero, bail out.
!
  if ( det == 0.0D+00 ) then

    b(1:2,1:2) = 0.0D+00
  end if
!
!  Compute the entries of the inverse matrix using an explicit formula.
!
  b(1,1) = + a(2,2) / det
  b(1,2) = - a(1,2) / det
  b(2,1) = - a(2,1) / det
  b(2,2) = + a(1,1) / det
end

subroutine r8vec2_compare ( n, a1, a2, i, j, isgn )

!*****************************************************************************80
!
!! R8VEC2_COMPARE compares elements of an R8VEC2.
!
!  Discussion:
!
!    The lexicographic ordering is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data items.
!
!    Input, double precision A1(N), A2(N), contain the two components
!    of each item.
!
!    Input, integer I, J, the items to be compared.
!
!    Output, integer ISGN, the results of the comparison:
!    -1, item I < item J,
!     0, item I = item J,
!    +1, item I > item J.
!
  implicit none

  integer n

  double precision a1(n)
  double precision a2(n)
  integer i
  integer isgn
  integer j

  isgn = 0

       if ( a1(i) < a1(j) ) then

    isgn = -1

  else if ( a1(i) == a1(j) ) then

         if ( a2(i) < a2(j) ) then
      isgn = -1
    else if ( a2(i) < a2(j) ) then
      isgn = 0
    else if ( a2(i) > a2(j) ) then
      isgn = +1
    end if

  else if ( a1(i) > a1(j) ) then

    isgn = +1

  end if
end

subroutine r8vec2_print ( n, a1, a2, title )

!*****************************************************************************80
!
!! R8VEC2_PRINT prints a pair of real vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, double precision A1(N), A2(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer n

  double precision a1(n)
  double precision a2(n)
  integer i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i6,2g14.6)' ) i, a1(i), a2(i)
  end do
end

subroutine r8vec2_sort_a ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC2_SORT_A ascending sorts a vector of pairs of integers.
!
!  Discussion:
!
!    Each item to be sorted is a pair (I,J), with the I
!    and J values stored in separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items of data.
!
!    Input/output, double precision A1(N), A2(N), the data to be sorted.
!
  implicit none

  integer n

  double precision a1(n)
  double precision a2(n)
  integer i
  integer indx
  integer isgn
  integer j
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call r8_swap ( a1(i), a1(j) )
      call r8_swap ( a2(i), a2(j) )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call r8vec2_compare ( n, a1, a2, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do
end

subroutine rect_int_2d ( x1, y1, x2, y2, x3, y3, x4, y4, flag, x5, y5, x6, y6 )

!*****************************************************************************80
!
!! RECT_INT_2D computes the intersection of two rectangles in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision X1, Y1, X2, Y2, the corners of the first segment.
!
!    Input, double precision X3, Y3, X4, Y4, the corners of the second segment.
!
!    Output, integer FLAG, records the results.
!    0, the rectangles do not intersect.
!    1, the intersection is a point.
!    2, the intersection is a line.
!    3, the intersection is a rectangle.
!
!    Output, double precision X5, Y5, X6, Y6, the corners of the intersection.
!    If FLAG = 0, X5 = Y5 = X6 = Y6 = 0.
!
  implicit none

  integer flag
  double precision x1
  double precision x2
  double precision x3
  double precision x4
  double precision x5
  double precision x6
  double precision y1
  double precision y2
  double precision y3
  double precision y4
  double precision y5
  double precision y6

  x5 = 0.0D+00
  y5 = 0.0D+00
  x6 = 0.0D+00
  y6 = 0.0D+00

  call lines_seg_int_1d ( x1, x2, x3, x4, flag, x5, x6 )

  if ( flag == 0 ) then
  end if

  call lines_seg_int_1d ( y1, y2, y3, y4, flag, y5, y6 )
end

subroutine triangle_contains_point_2d ( x1, y1, x2, y2, x3, y3, x, y, inside )

!*****************************************************************************80
!
!! TRIANGLE_CONTAINS_POINT_2D finds if a point is inside a triangle in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision X1, Y1, X2, Y2, X3, Y3, the coordinates of
!    the corners of the triangle.
!
!    Input, double precision X, Y, the point to be checked.
!
!    Output, logical INSIDE, is .TRUE. if (X,Y) is inside
!    the triangle or on its boundary, and .FALSE. otherwise.
!
  implicit none

  integer , parameter :: N = 2
  integer , parameter :: NRHS = 1

  double precision a(N,N+NRHS)
  double precision c1
  double precision c2
  integer info
  logical inside
  double precision x
  double precision x1
  double precision x2
  double precision x3
  double precision y
  double precision y1
  double precision y2
  double precision y3
!
!  Set up the linear system
!
!    ( X2-X1  X3-X1 ) C1  = X-X1
!    ( Y2-Y1  Y3-Y1 ) C2    Y-Y1
!
!  which is satisfied by the barycentric coordinates of (X,Y).
!
  a(1,1) = x2 - x1
  a(1,2) = x3 - x1
  a(1,3) = x - x1

  a(2,1) = y2 - y1
  a(2,2) = y3 - y1
  a(2,3) = y - y1
!
!  Solve the linear system.
!
  call r8mat_solve ( a, N, NRHS, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE_CONTAINS_POINT - Fatal error!'
    write ( *, '(a)' ) '  The linear system is singular.'
    write ( *, '(a)' ) '  The input data does not form a proper triangle.'
    stop
  end if

  c1 = a(1,3)
  c2 = a(2,3)
!
!  If the point is in the triangle, its barycentric coordinates
!  must both be nonnegative, and sum to no more than 1.
!
  if ( c1 < 0.0D+00 .or. c2 < 0.0D+00 ) then
    inside = .false.
  else if ( 1.0D+00 < c1 + c2 ) then
    inside = .false.
  else
    inside = .true.
  end if
end

subroutine triangulate_tricolor ( node_num, triang, color )

!*****************************************************************************80
!
!! TRIANGULATE_TRICOLOR three-colors the nodes of a triangulated polygon.
!
!  Discussion:
!
!    While the data in TRIANG must represent a proper triangulation
!    of the polygon, no check for this is made.
!
!    The three-coloring is done in such a way that every triangle
!    of the triangulation has one node of each color.
!
!    This solves the art gallery problem, which asks for the minimum
!    number of guards necessary in a (polygonal) art gallery, so that
!    every point inside the gallery is seen by at least one guard.
!    The answer is that (N/3) (integer division, with truncation) guards
!    are sufficient.  Triangulate the polygon, tricolor the triangulation,
!    choose the color that shows up least often (in case N is not divisible
!    by 3) and place a guard at those points.  Since every triangle has
!    a node of that color, and the triangles make up the gallery, you're done.
!    This assumes that the gallery is a simple polygon.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marc de Berg, Marc van Kreveld, Mark Overmars, Otfried Schwarzkopf,
!    Computational Geometry,
!    Second Edition,
!    Springer, 2000, pages 47-48.
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes in the polygon.
!
!    Input, integer TRIANG(3,NODE_NUM-2), the triangulation of
!    the polygon.
!
!    Output, integer COLOR(NODE_NUM), an assignment of the
!    "colors" 1, 2 and 3 to the triangles in such a way that every triangle
!    in the triangulation has one node of each color.
!
  implicit none

  integer node_num

  integer color(node_num)
  integer color1
  integer color2
  integer color3
  integer node1
  integer node2
  integer node3
  integer stack(6*(node_num-3))
  integer stack_max
  integer stack_num
  integer t1
  integer t2
  integer triang(3,node_num-2)

  stack_max = 6 * ( node_num - 3 )
  t1 = 1

  node1 = triang(1,t1)
  node2 = triang(2,t1)
  node3 = triang(3,t1)

  color1 = 1
  color2 = 2
  color3 = 3

  color(node1) = color1
  color(node2) = color2
  color(node3) = color3

  stack_num = 0

  call triangulate_color_push ( t1, node1, color1, node2, color2, color3, &
    stack_max, stack_num, stack )

  call triangulate_color_push ( t1, node2, color2, node3, color3, color1, &
    stack_max, stack_num, stack )

  call triangulate_color_push ( t1, node3, color3, node1, color1, color2, &
    stack_max, stack_num, stack )

  do while ( 0 < stack_num )

    call triangulate_color_pop ( t1, node1, color1, node2, color2, color3, &
      stack_max, stack_num, stack )

    call triangulate_common_edge ( node_num-2, triang, node1, node2, t1, &
      t2, node3 )

    if ( t2 /= 0 ) then

      color(node3) = color3

      call triangulate_color_push ( t2, node1, color1, node3, color3, color2, &
        stack_max, stack_num, stack )

      call triangulate_color_push ( t2, node3, color3, node2, color2, color1, &
        stack_max, stack_num, stack )

    end if

  end do
end

subroutine triangulate_color_push ( t, node1, color1, node2, color2, color3, &
  stack_max, stack_num, stack )

!*****************************************************************************80
!
!! TRIANGULATE_COLOR_PUSH pushes a side of a colored triangle onto the stack.
!
!  Discussion:
!
!    This is a utility routine used by TRIANGULATE_TRICOLOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer T, the triangle that the side belongs to.
!
!    Input, integer NODE1, COLOR1, the starting node of the
!    edge, and its color.
!
!    Input, integer NODE2, COLOR2, the end node of the edge,
!    and its color.
!
!    Input, integer COLOR3, the remaining color.
!
!    Input, integer STACK_MAX, the maximum size of the stack.
!
!    Input/output, integer STACK_NUM, the current size of 
!    the stack.
!
!    Input/output, integer STACK(STACK_MAX), the stack.
!
  implicit none

  integer stack_max

  integer color1
  integer color2
  integer color3
  integer node1
  integer node2
  integer stack(stack_max)
  integer stack_num
  integer t

  stack_num = stack_num + 1
  stack(stack_num) = t
  stack_num = stack_num + 1
  stack(stack_num) = node1
  stack_num = stack_num + 1
  stack(stack_num) = color1
  stack_num = stack_num + 1
  stack(stack_num) = node2
  stack_num = stack_num + 1
  stack(stack_num) = color2
  stack_num = stack_num + 1
  stack(stack_num) = color3
end

subroutine triangulate_color_pop ( t, node1, color1, node2, color2, color3, &
  stack_max, stack_num, stack )

!*****************************************************************************80
!
!! TRIANGULATE_COLOR_POP pops a side of a colored triangle from the stack.
!
!  Discussion:
!
!    This is a utility routine used by TRIANGULATE_TRICOLOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer T, the triangle that the side belongs to.
!
!    Output, integer NODE1, COLOR1, the starting node of the
!    edge, and its color.
!
!    Output, integer NODE2, COLOR2, the end node of the edge,
!    and its color.
!
!    Output, integer COLOR3, the remaining color.
!
!    Input, integer STACK_MAX, the maximum size of the stack.
!
!    Input/output, integer STACK_NUM, the current size of 
!    the stack.
!
!    Input/output, integer STACK(STACK_MAX), the stack.
!
  implicit none

  integer stack_max

  integer color1
  integer color2
  integer color3
  integer node1
  integer node2
  integer stack(stack_max)
  integer stack_num
  integer t

  color3 = stack(stack_num)
  stack_num = stack_num - 1
  color2 = stack(stack_num)
  stack_num = stack_num - 1
  node2 = stack(stack_num)
  stack_num = stack_num - 1
  color1 = stack(stack_num)
  stack_num = stack_num - 1
  node1 = stack(stack_num)
  stack_num = stack_num - 1
  t = stack(stack_num)
  stack_num = stack_num - 1
end

subroutine triangulate_common_edge ( triang_num, triang, node1, node2, t1, &
  t2, node3 )

!*****************************************************************************80
!
!! TRIANGULATE_COMMON_EDGE seeks the other triangle that shares an edge.
!
!  Discussion:
!
!    This is a utility routine used by TRIANGULATE_TRICOLOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TRIANG_NUM, the number of triangles.
!
!    Input, integer TRIANG(3,TRIANG_NUM), the triangulation
!    of the polygon.
!
!    Input, integer NODE1, NODE2, the starting and ending nodes
!    of an edge.
!
!    Input, integer T1, the triangle to which the edge
!    (NODE1,NODE2) belongs.
!
!    Output, integer T2, is 0 if there is no triangle containing
!    the matching edge (NODE2,NODE1), or it is the index of a triangle
!    containing the edge (NODE2,NODE1).
!
!    Output, integer NODE3, the other node in triangle T2.
!
  implicit none

  integer triang_num

  integer i4_wrap
  integer j1
  integer j2
  integer j3
  integer node1
  integer node2
  integer node3
  integer t
  integer t1
  integer t2
  integer triang(3,triang_num)

  t2 = 0
  node3 = 0

  do t = 1, triang_num
    do j1 = 1, 3
      j2 = i4_wrap ( j1+1, 1, 3 )
      if ( triang(j2,t) == node1 .and. triang(j1,t) == node2 ) then
        t2 = t
        j3 = i4_wrap ( j1+2, 1, 3 )
        node3 = triang(j3,t)
      end if
    end do
  end do
end

subroutine triangulation_boundary_count ( point_num, tri_num, bound_num )

!*****************************************************************************80
!
!! TRIANGULATION_BOUNDARY_COUNT returns the number of boundary edges.
!
!  Discussion:
!
!    We assume we are given information about a legal, maximal triangulation
!    of a set of points in the plane.  In a maximal triangulation, no more
!    edges can be added without crossing an existing edge.
!
!    Given the number of points and triangles, we are going to apply
!    Euler's formula to determine the number of edges that lie on the
!    convex hull of the set of points.
!
!    The number of faces, including the infinite face, is TRI_NUM + 1.
!
!    Let BOUND_NUM denote the number of edges on the boundary of the convex
!    hull.  Each of the TRI_NUM triangles uses three edges.  Every edge
!    occurs in two different faces, so the number of edges must be
!    ( 3 * TRI_NUM + BOUND_NUM ) / 2.
!
!    The number of points is POINT_NUM.
!
!    Euler's formula asserts that, for a simple connected figure in the
!    plane with no edge crossings, POINT_NUM points, EDGE_NUM edges and
!    FACE_NUM faces:
!
!      POINT_NUM - EDGE_NUM + FACE_NUM = 2
!
!    In our context, this becomes
!
!      POINT_NUM - ( 3 * TRI_NUM + BOUND_NUM ) / 2 + TRI_NUM + 1 = 2
!
!    or
!
!      BOUND_NUM = 2 * POINT_NUM - TRI_NUM - 2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marc de Berg, Marc van Kreveld, Mark Overmars, Otfried Schwarzkopf,
!    Computational Geometry, Section 9.1,
!    Springer, 2000.
!
!  Parameters:
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, integer TRI_NUM, the number of triangles.
!
!    Output, integer BOUND_NUM, the number of edges that lie
!    on the convex hull of the triangulation.
!
  implicit none

  integer bound_num
  integer point_num
  integer tri_num

  bound_num = 2 * point_num - tri_num - 2
end
