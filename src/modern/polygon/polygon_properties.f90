!> polygon_properties — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module polygon_properties_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  implicit none
  private

  public :: angle_half, angle_rad, between, collinear, diagonal, diagonalie
  public :: i4_modp, i4_wrap, in_cone, intersect, intersect_prop, l4_xor
  public :: polygon_angles, polygon_area, polygon_area_2, polygon_centroid, polygon_centroid_2, polygon_contains_point
  public :: polygon_contains_point_2, polygon_contains_point_3, polygon_diameter, polygon_expand, polygon_inrad_data, polygon_integral_1
  public :: polygon_integral_x, polygon_integral_xx, polygon_integral_xy, polygon_integral_y, polygon_integral_yy, polygon_is_convex
  public :: polygon_lattice_area, polygon_outrad_data, polygon_perimeter, polygon_perimeter_quad, polygon_point_dist, polygon_point_near
  public :: polygon_sample, polygon_side_data, polygon_triangulate, r8_degrees, r8_uniform_01, r8mat_solve
  public :: r8vec_uniform_01, segment_point_dist, segment_point_near, triangle_area, triangle_barycentric, triangle_contains_point_1

contains

  subroutine angle_half ( p1, p2, p3, p4 )

  !*****************************************************************************80
  !
  !! ANGLE_HALF finds half an angle.
  !
  !  Discussion:
  !
  !    The original angle is defined by the sequence of points P1, P2 and P3.
  !
  !    The point P4 is calculated so that:
  !
  !      (P1,P2,P4) = (P1,P2,P3) / 2
  !
  !        P1
  !        /
  !       /   P4
  !      /  .  
  !     / .
  !    P2--------->P3
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    01 March 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real(real64) P1(2), P2(2), P3(2), points defining the angle. 
  !
  !    Input, real(real64) P4(2), a point defining the half angle.
  !    The vector P4 - P2 will have unit norm.
  !

    real(real64) p1(2)
    real(real64) p2(2)
    real(real64) p3(2)
    real(real64) p4(2)

    p4(1:2) = 0.5e+00_real64 * ( &
        ( p1(1:2) - p2(1:2) ) / sqrt ( sum ( ( p1(1:2) - p2(1:2) )**2 ) ) &
      + ( p3(1:2) - p2(1:2) ) / sqrt ( sum ( ( p3(1:2) - p2(1:2) )**2 ) ) )

     p4(1:2) = p2(1:2) + p4(1:2) / sqrt ( sum ( p4(1:2)**2 ) )
  end

  function angle_rad ( p1, p2, p3 )

  !*****************************************************************************80
  !
  !! ANGLE_RAD returns the angle in radians swept out between two rays.
  !
  !  Discussion:
  !
  !    Except for the zero angle case, it should be true that
  !
  !      ANGLE_RAD ( P1, P2, P3 ) + ANGLE_RAD ( P3, P2, P1 ) = 2 * PI
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
  !    Input, real(real64) P1(2), P2(2), P3(2), define the rays
  !    P1 - P2 and P3 - P2 which define the angle.
  !
  !    Output, real(real64) ANGLE_RAD, the angle swept out by the rays,
  !    in radians.  0 <= ANGLE_RAD < 2 * PI.  If either ray has zero
  !    length, then ANGLE_RAD is set to 0.
  !

    real(real64) angle_rad
    real(real64) p(2)
    real(real64) p1(2)
    real(real64) p2(2)
    real(real64) p3(2)
    real(real64), parameter :: r8_pi = 3.141592653589793e+00_real64

    p(1) = ( p3(1) - p2(1) ) * ( p1(1) - p2(1) ) &
         + ( p3(2) - p2(2) ) * ( p1(2) - p2(2) )

    p(2) = ( p3(1) - p2(1) ) * ( p1(2) - p2(2) ) &
         - ( p3(2) - p2(2) ) * ( p1(1) - p2(1) )

    if ( all ( p(1:2) == 0.0e+00_real64)  ) then
      angle_rad = 0.0e+00_real64
    end if

    angle_rad = atan2 ( p(2), p(1) )

    if ( angle_rad < 0.0e+00_real64 ) then
      angle_rad = angle_rad + 2.0e+00_real64 * r8_pi
    end if
  end

  function between ( xa, ya, xb, yb, xc, yc )

  !*****************************************************************************80
  !
  !! BETWEEN is TRUE if vertex C is between vertices A and B.
  !
  !  Discussion:
  !
  !    The points must be (numerically) collinear.
  !
  !    Given that condition, we take the greater of XA - XB and YA - YB
  !    as a "scale" and check where C's value lies.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 May 2014
  !
  !  Author:
  !
  !    Original C version by Joseph ORourke.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Joseph ORourke,
  !    Computational Geometry in C,
  !    Cambridge, 1998,
  !    ISBN: 0521649765,
  !    LC: QA448.e38_real64.
  !
  !  Parameters:
  !
  !    Input, real(real64) XA, YA, XB, YB, XC, YC, the coordinates of 
  !    the vertices.
  !
  !    Output, logical BETWEEN, is TRUE if C is between A and B.
  !

    logical between
    logical collinear
    logical value
    real(real64) xa
    real(real64) xb
    real(real64) xc
    real(real64) xmax
    real(real64) xmin
    real(real64) ya
    real(real64) yb
    real(real64) yc
    real(real64) ymax
    real(real64) ymin

    if ( .not. collinear ( xa, ya, xb, yb, xc, yc ) ) then
      value = .false.
    else if ( abs ( ya - yb ) < abs ( xa - xb ) ) then
      xmax = max ( xa, xb )
      xmin = min ( xa, xb )
      value = ( xmin <= xc .and. xc <= xmax )
    else
      ymax = max ( ya, yb )
      ymin = min ( ya, yb )
      value = ( ymin <= yc .and. yc <= ymax )
    end if

    between = value
  end

  function collinear ( xa, ya, xb, yb, xc, yc )

  !*****************************************************************************80
  !
  !! COLLINEAR returns a measure of collinearity for three points.
  !
  !  Discussion:
  !
  !    In order to deal with collinear points whose coordinates are not
  !    numerically exact, we compare the area of the largest square
  !    that can be created by the line segment between two of the points
  !    to (twice) the area of the triangle formed by the points.
  !
  !    If the points are collinear, their triangle has zero area.
  !    If the points are close to collinear, then the area of this triangle
  !    will be small relative to the square of the longest segment.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 May 2014
  !
  !  Author:
  !
  !    Original C version by Joseph ORourke.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Joseph ORourke,
  !    Computational Geometry in C,
  !    Cambridge, 1998,
  !    ISBN: 0521649765,
  !    LC: QA448.e38_real64.
  !
  !  Parameters:
  !
  !    Input, real(real64) XA, YA, XB, YB, XC, YC, the coordinates of 
  !    the vertices.
  !
  !    Output, logical COLLINEAR, is TRUE if the points are judged 
  !    to be collinear.
  !

    real(real64) area
    logical collinear
    real(real64), parameter :: r8_eps = 2.220446049250313e-016_real64
    real(real64) side_ab_sq
    real(real64) side_bc_sq
    real(real64) side_ca_sq
    real(real64) side_max_sq
    logical value
    real(real64) xa
    real(real64) xb
    real(real64) xc
    real(real64) ya
    real(real64) yb
    real(real64) yc

    area = 0.5e+00_real64 * ( &
        ( xb - xa ) * ( yc - ya ) &
      - ( xc - xa ) * ( yb - ya ) )

    side_ab_sq = ( xa - xb ) ** 2 + ( ya - yb ) ** 2
    side_bc_sq = ( xb - xc ) ** 2 + ( yb - yc ) ** 2
    side_ca_sq = ( xc - xa ) ** 2 + ( yc - ya ) ** 2

    side_max_sq = max ( side_ab_sq, max ( side_bc_sq, side_ca_sq ) )

    if ( side_max_sq <= r8_eps ) then
      value = .true.
    else if ( 2.0e+00_real64 * abs ( area ) <= r8_eps * side_max_sq ) then
      value = .true.
    else
      value = .false.
    end if

    collinear = value
  end

  function diagonal ( im1, ip1, n, prev, next, x, y )

  !*****************************************************************************80
  !
  !! DIAGONAL: VERTEX(IM1) to VERTEX(IP1) is a proper internal diagonal.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 May 2014
  !
  !  Author:
  !
  !    Original C version by Joseph ORourke.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Joseph ORourke,
  !    Computational Geometry in C,
  !    Cambridge, 1998,
  !    ISBN: 0521649765,
  !    LC: QA448.e38_real64.
  !
  !  Parameters:
  !
  !    Input, integer(int32) IM1, IP1, the indices of two vertices.
  !
  !    Input, integer(int32) N, the number of vertices.
  !
  !    Input, integer(int32) PREV(N), the previous neighbor of each vertex.
  !
  !    Input, integer(int32) NEXT(N), the next neighbor of each vertex.
  !
  !    Input, real(real64) X(N), Y(N), the coordinates of each vertex.
  !
  !    Output, logical DIAGONAL, the value of the test.
  !

    integer(int32) n

    logical diagonal
    logical diagonalie
    integer(int32) im1
    logical in_cone
    integer(int32) ip1
    integer(int32) next(n)
    integer(int32) prev(n)
    logical value1
    logical value2
    logical value3
    real(real64) x(n)
    real(real64) y(n)

    value1 = in_cone ( im1, ip1, n, prev, next, x, y )
    value2 = in_cone ( ip1, im1, n, prev, next, x, y )
    value3 = diagonalie ( im1, ip1, n, next, x, y )

    diagonal = ( value1 .and. value2 .and. value3 )
  end

  function diagonalie ( im1, ip1, n, next, x, y )

  !*****************************************************************************80
  !
  !! DIAGONALIE is true if VERTEX(IM1):VERTEX(IP1) is a proper diagonal.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 May 2014
  !
  !  Author:
  !
  !    Original C version by Joseph ORourke.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Joseph ORourke,
  !    Computational Geometry in C,
  !    Cambridge, 1998,
  !    ISBN: 0521649765,
  !    LC: QA448.e38_real64.
  !
  !  Parameters:
  !
  !    Input, integer(int32) IM1, IP1, the indices of two vertices.
  !
  !    Input, integer(int32) N, the number of vertices.
  !
  !    Input, integer(int32) NEXT(N), the next neighbor of each vertex.
  !
  !    Input, real(real64) X(N), Y(N), the coordinates of each vertex.
  !
  !    Output, logical DIAGONALIE, the value of the test.
  !

    integer(int32) n

    logical diagonalie
    integer(int32) first
    integer(int32) im1
    logical intersect
    integer(int32) ip1
    integer(int32) j
    integer(int32) jp1
    integer(int32) next(n)
    logical value
    logical value2
    real(real64) x(n)
    real(real64) y(n)

    first = im1
    j = first
    jp1 = next(first)

    value = .true.
  !
  !  For each edge VERTEX(J):VERTEX(JP1) of the polygon:
  !
    do
  !
  !  Skip any edge that includes vertex IM1 or IP1.
  !
      if ( j == im1 .or. j == ip1 .or. jp1 == im1 .or. jp1 == ip1 ) then

      else

        value2 = intersect ( x(im1), y(im1), x(ip1), y(ip1), x(j), y(j), &
          x(jp1), y(jp1) )

        if ( value2 ) then
          value = .false.
          exit
        end if

      end if

      j = jp1
      jp1 = next(j)

      if ( j == first ) then
        exit
      end if

    end do

    diagonalie = value
  end

  function i4_modp ( i, j )

  !*****************************************************************************80
  !
  !! I4_MODP returns the nonnegative remainder of I4 division.
  !
  !  Discussion:
  !
  !    If
  !      NREM = I4_MODP ( I, J )
  !      NMULT = ( I - NREM ) / J
  !    then
  !      I = J * NMULT + NREM
  !    where NREM is always nonnegative.
  !
  !    The MOD function computes a result with the same sign as the
  !    quantity being divided.  Thus, suppose you had an angle A,
  !    and you wanted to ensure that it was between 0 and 360.
  !    Then mod(A,360) would do, if A was positive, but if A
  !    was negative, your result would be between -360 and 0.
  !
  !    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
  !
  !    An I4 is an integer(int32) value.
  !
  !  Example:
  !
  !        I     J     MOD I4_MODP    Factorization
  !
  !      107    50       7       7    107 =  2 *  50 + 7
  !      107   -50       7       7    107 = -2 * -50 + 7
  !     -107    50      -7      43   -107 = -3 *  50 + 43
  !     -107   -50      -7      43   -107 =  3 * -50 + 43
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    02 March 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(int32) I, the number to be divided.
  !
  !    Input, integer(int32) J, the number that divides I.
  !
  !    Output, integer(int32) I4_MODP, the nonnegative remainder when I is
  !    divided by J.
  !

    integer(int32) i
    integer(int32) i4_modp
    integer(int32) j
    integer(int32) value

    if ( j == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'I4_MODP - Fatal error!'
      write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
      stop 1
    end if

    value = mod ( i, j )

    if ( value < 0 ) then
      value = value + abs ( j )
    end if

    i4_modp = value
  end

  function i4_wrap ( ival, ilo, ihi )

  !*****************************************************************************80
  !
  !! I4_WRAP forces an I4 to lie between given limits by wrapping.
  !
  !  Discussion:
  !
  !    An I4 is an integer(int32) value.
  !
  !    There appears to be a bug in the GFORTRAN compiler which can lead to
  !    erroneous results when the first argument of I4_WRAP is an expression.
  !    In particular:
  !
  !    do i = 1, 3
  !      if ( test ) then
  !        i4 = i4_wrap ( i + 1, 1, 3 )
  !      end if
  !    end do
  !
  !    was, when I = 3, returning I4 = 3.  So I had to replace this with
  !
  !    do i = 1, 3
  !      if ( test ) then
  !        i4 = i + 1
  !        i4 = i4_wrap ( i4, 1, 3 )
  !      end if
  !    end do
  !
  !  Example:
  !
  !    ILO = 4, IHI = 8
  !
  !    I  Value
  !
  !    -2     8
  !    -1     4
  !     0     5
  !     1     6
  !     2     7
  !     3     8
  !     4     4
  !     5     5
  !     6     6
  !     7     7
  !     8     8
  !     9     4
  !    10     5
  !    11     6
  !    12     7
  !    13     8
  !    14     4
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 September 2009
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(int32) IVAL, a value.
  !
  !    Input, integer(int32) ILO, IHI, the desired bounds.
  !
  !    Output, integer(int32) I4_WRAP, a "wrapped" version of the value.
  !

    integer(int32) i4_modp
    integer(int32) i4_wrap
    integer(int32) ihi
    integer(int32) ilo
    integer(int32) ival
    integer(int32) jhi
    integer(int32) jlo
    integer(int32) value
    integer(int32) wide

    jlo = min ( ilo, ihi )
    jhi = max ( ilo, ihi )

    wide = jhi - jlo + 1

    if ( wide == 1 ) then
      value = jlo
    else
      value = jlo + i4_modp ( ival - jlo, wide )
    end if

    i4_wrap = value
  end

  function in_cone ( im1, ip1, n, prev, next, x, y )

  !*****************************************************************************80
  !
  !! IN_CONE is TRUE if the diagonal VERTEX(IM1):VERTEX(IP1) is strictly internal.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 May 2014
  !
  !  Author:
  !
  !    Original C version by Joseph ORourke.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Joseph ORourke,
  !    Computational Geometry in C,
  !    Cambridge, 1998,
  !    ISBN: 0521649765,
  !    LC: QA448.e38_real64.
  !
  !  Parameters:
  !
  !    Input, integer(int32) IM1, IP1, the indices of two vertices.
  !
  !    Input, integer(int32) N, the number of vertices.
  !
  !    Input, integer(int32) PREV(N), the previous neighbor of each vertex.
  !
  !    Input, integer(int32) NEXT(N), the next neighbor of each vertex.
  !
  !    Input, real(real64) X(N), Y(N), the coordinates of each vertex.
  !
  !    Output, logical IN_CONE, the value of the test.
  !

    integer(int32) n

    integer(int32) i
    integer(int32) im1
    integer(int32) im2
    logical in_cone
    integer(int32) ip1
    integer(int32) next(n)
    integer(int32) prev(n)
    real(real64) t1
    real(real64) t2
    real(real64) t3
    real(real64) t4
    real(real64) t5
    real(real64) triangle_area
    logical value
    real(real64) x(n)
    real(real64) y(n)

    im2 = prev(im1)
    i = next(im1)

    t1 = triangle_area ( x(im1), y(im1), x(i), y(i), x(im2), y(im2) )

    if ( 0.0e+00_real64 <= t1 ) then

      t2 = triangle_area ( x(im1), y(im1), x(ip1), y(ip1), x(im2), y(im2) )
      t3 = triangle_area ( x(ip1), y(ip1), x(im1), y(im1), x(i), y(i) )
      value = ( ( 0.0e+00_real64 < t2 ) .and. ( 0.0e+00_real64 < t3 ) )

    else

      t4 = triangle_area ( x(im1), y(im1), x(ip1), y(ip1), x(i), y(i) )
      t5 = triangle_area ( x(ip1), y(ip1), x(im1), y(im1), x(im2), y(im2) )
      value = .not. ( ( 0.0e+00_real64 <= t4 ) .and. ( 0.0e+00_real64 <= t5 ) )

    end if

    in_cone = value
  end

  function intersect ( xa, ya, xb, yb, xc, yc, xd, yd )

  !*****************************************************************************80
  !
  !! INTERSECT is true if lines VA:VB and VC:VD intersect.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 May 2014
  !
  !  Author:
  !
  !    Original C version by Joseph ORourke.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Joseph ORourke,
  !    Computational Geometry in C,
  !    Cambridge, 1998,
  !    ISBN: 0521649765,
  !    LC: QA448.e38_real64.
  !
  !  Parameters:
  !
  !    Input, real(real64) XA, YA, XB, YB, XC, YC, XD, YD, the X and Y 
  !    coordinates of the four vertices.
  !
  !    Output, logical VALUE, the value of the test.
  !

    logical between
    logical intersect
    logical intersect_prop
    logical value
    real(real64) xa
    real(real64) xb
    real(real64) xc
    real(real64) xd
    real(real64) ya
    real(real64) yb
    real(real64) yc
    real(real64) yd

    if ( intersect_prop ( xa, ya, xb, yb, xc, yc, xc, yd ) ) then
      value = .true.
    else if ( between ( xa, ya, xb, yb, xc, yc ) ) then
      value = .true.
    else if ( between ( xa, ya, xb, yb, xd, yd ) ) then
      value = .true.
    else if ( between ( xc, yc, xd, yd, xa, ya ) ) then
      value = .true.
    else if ( between ( xc, yc, xd, yd, xb, yb ) ) then
      value = .true.
    else
      value = .false.
    end if

    intersect = value
  end

  function intersect_prop ( xa, ya, xb, yb, xc, yc, xd, yd )

  !*****************************************************************************80
  !
  !! INTERSECT_PROP is TRUE if lines VA:VB and VC:VD have a proper intersection.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 May 2014
  !
  !  Author:
  !
  !    Original C version by Joseph ORourke.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Joseph ORourke,
  !    Computational Geometry in C,
  !    Cambridge, 1998,
  !    ISBN: 0521649765,
  !    LC: QA448.e38_real64.
  !
  !  Parameters:
  !
  !    Input, real(real64) XA, YA, XB, YB, XC, YC, XD, YD, the X and Y 
  !    coordinates of the four vertices.
  !
  !    Output, logical INTERSECT_PROP, the result of the test.
  !

    logical collinear
    logical intersect_prop
    logical l4_xor
    real(real64) t1
    real(real64) t2
    real(real64) t3
    real(real64) t4
    real(real64) triangle_area
    logical value
    logical value1
    logical value2
    logical value3
    logical value4
    real(real64) xa
    real(real64) xb
    real(real64) xc
    real(real64) xd
    real(real64) ya
    real(real64) yb
    real(real64) yc
    real(real64) yd

    if ( collinear ( xa, ya, xb, yb, xc, yc ) ) then
      value = .false.
    else if ( collinear ( xa, ya, xb, yb, xd, yd ) ) then
      value = .false.
    else if ( collinear ( xc, yc, xd, yd, xa, ya ) ) then
      value = .false.
    else if ( collinear ( xc, yc, xd, yd, xb, yb ) ) then
      value = .false.
    else

      t1 = triangle_area ( xa, ya, xb, yb, xc, yc )
      t2 = triangle_area ( xa, ya, xb, yb, xd, yd )
      t3 = triangle_area ( xc, yc, xd, yd, xa, ya )
      t4 = triangle_area ( xc, yc, xd, yd, xb, yb )

      value1 = ( 0.0e+00_real64 < t1 )
      value2 = ( 0.0e+00_real64 < t2 )
      value3 = ( 0.0e+00_real64 < t3 )
      value4 = ( 0.0e+00_real64 < t4 )

      value = ( l4_xor ( value1, value2 ) ) .and. ( l4_xor ( value3, value4 ) )

    end if

    intersect_prop = value
  end

  function l4_xor ( l1, l2 )

  !*****************************************************************************80
  !
  !! L4_XOR returns the exclusive OR of two L4's.
  !
  !  Discussion:
  !
  !    An L4 is a logical value.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 May 2014
  !
  !  Author:
  !
  !   John Burkardt
  !
  !  Parameters:
  !
  !    Input, logical L1, L2, two values whose exclusive OR 
  !    is needed.
  !
  !    Output, logical L4_XOR, the exclusive OR of L1 and L2.
  !

    logical l1
    logical l2
    logical l4_xor
    logical value1
    logical value2

    value1 = (         l1   .and. ( .not. l2 ) )
    value2 = ( ( .not. l1 ) .and.         l2   )

    l4_xor = ( value1 .or. value2 )
  end

  subroutine polygon_angles ( n, v, angle )

  !*****************************************************************************80
  !
  !! POLYGON_ANGLES computes the interior angles of a polygon.
  !
  !  Discussion:
  !
  !    The vertices should be listed in counter clockwise order.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    14 March 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of vertices of the polygon.
  !
  !    Input, real(real64) V(2,N), the vertices.
  !
  !    Output, real(real64) ANGLE(N), the angles of the polygon,
  !    in radians.
  !

    integer(int32) n

    real(real64) angle(n)
    real(real64) angle_rad
    integer(int32) i
    integer(int32) i4_wrap
    integer(int32) im1
    integer(int32) ip1
    real(real64) v(2,n)

    if ( n < 3 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLYGON_ANGLES - Fatal error!'
      write ( *, '(a)' ) '  The number of vertices must be at least 3.'
      write ( *, '(a,i8)' ) '  The input value of N = ', n
      stop 1
    end if

    do i = 1, n

      im1 = i4_wrap ( i - 1, 1, n )
      ip1 = i4_wrap ( i + 1, 1, n )

      angle(i) = angle_rad ( v(1:2,im1), v(1:2,i), v(1:2,ip1) )

    end do
  end

  subroutine polygon_area ( n, v, area )

  !*****************************************************************************80
  !
  !! POLYGON_AREA computes the area of a polygon.
  !
  !  Discussion:
  !
  !    AREA = 1/2 * abs ( sum ( 1 <= I <= N ) X(I) * ( Y(I+1) - Y(I-1) ) )
  !    where Y(0) should be replaced by Y(N), and Y(N+1) by Y(1).
  !
  !    If the vertices are given in counter clockwise order, the area
  !    will be positive.  If the vertices are given in clockwise order,
  !    the area will be negative.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    17 October 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of vertices of the polygon.
  !
  !    Input, real(real64) V(2,N), the vertices.
  !
  !    Output, real(real64) AREA, the absolute area of the polygon.
  !

    integer(int32) n

    real(real64) area
    integer(int32) i
    integer(int32) i4_wrap
    integer(int32) im1
    integer(int32) ip1
    real(real64) v(2,n)

    area = 0.0e+00_real64

    do i = 1, n

      im1 = i4_wrap ( i-1, 1, n )
      ip1 = i4_wrap ( i+1, 1, n )

      area = area + v(1,i) * ( v(2,ip1) - v(2,im1) )

    end do

    area = 0.5e+00_real64 * area
  end

  subroutine polygon_area_2 ( n, v, area )

  !*****************************************************************************80
  !
  !! POLYGON_AREA_2 computes the area of a polygon.
  !
  !  Discussion:
  !
  !    The area is the sum of the areas of the triangles formed by
  !    node N with consecutive pairs of nodes.
  !
  !    If the vertices are given in counter clockwise order, the area
  !    will be positive.  If the vertices are given in clockwise order,
  !    the area will be negative.
  !
  !    Thanks to Martin Pineault for noticing that an earlier version
  !    of this routine would not correctly compute the area of a nonconvex
  !    polygon.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    17 October 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Adrian Bowyer, John Woodwark,
  !    A Programmer's Geometry,
  !    Butterworths, 1983,
  !    ISBN: 0408012420.
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of vertices of the polygon.
  !
  !    Input, real(real64) V(2,N), the vertices.
  !
  !    Output, real(real64) AREA, the absolute area of the polygon.
  !

    integer(int32) n

    real(real64) area
    real(real64) area_triangle
    integer(int32) i
    real(real64) triangle_area
    real(real64) v(2,n)

    area = 0.0e+00_real64

    do i = 1, n - 2

      area_triangle = triangle_area ( &
        v(1,i),   v(2,i), &
        v(1,i+1), v(2,i+1), &
        v(1,n),   v(2,n) )

      area = area + area_triangle

    end do
  end

  subroutine polygon_centroid ( n, v, centroid )

  !*****************************************************************************80
  !
  !! POLYGON_CENTROID computes the centroid of a polygon.
  !
  !  Discussion:
  !
  !    Denoting the centroid coordinates by CENTROID, then
  !
  !      CENTROID(1) = Integral ( Polygon interior ) x dx dy / Area ( Polygon )
  !      CENTROID(2) = Integral ( Polygon interior ) y dx dy / Area ( Polygon ).
  !
  !    Green's theorem states that for continuously differentiable functions
  !    M(x,y) and N(x,y),
  !
  !      Integral ( Polygon boundary ) ( M dx + N dy ) =
  !      Integral ( Polygon interior ) ( dN/dx - dM/dy ) dx dy.
  !
  !    Using M(x,y) = 0 and N(x,y) = x*x/2, we get:
  !
  !      CENTROID(1) = 0.5 * Integral ( Polygon boundary ) x*x dy 
  !                  / Area ( Polygon ),
  !
  !    which becomes
  !
  !      CENTROID(1) = 1/6 sum ( 1 <= I <= N )
  !        ( X(I+1) + X(I) ) * ( X(I) * Y(I+1) - X(I+1) * Y(I))
  !        / Area ( Polygon )
  !
  !    where, when I = N, the index "I+1" is replaced by 1.
  !
  !    A similar calculation gives us a formula for CENTROID(2).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    12 July 2003
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Gerard Bashein, Paul Detmer,
  !    Centroid of a Polygon,
  !    in Graphics Gems IV, 
  !    edited by Paul Heckbert,
  !    AP Professional, 1994,
  !    T385.G6974.
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of sides of the polygon.
  !
  !    Input, real(real64) V(2,N), the coordinates of the vertices.
  !
  !    Output, real(real64) CENTROID(2), the coordinates of the centroid.
  !

    integer(int32) n

    real(real64) area
    real(real64) centroid(2)
    integer(int32) i
    integer(int32) ip1
    real(real64) temp
    real(real64) v(2,n)

    area = 0.0e+00_real64
    centroid(1:2) = 0.0e+00_real64

    do i = 1, n

      if ( i < n ) then
        ip1 = i + 1
      else
        ip1 = 1
      end if

      temp = ( v(1,i) * v(2,ip1) - v(1,ip1) * v(2,i) )

      area = area + temp

      centroid(1:2) = centroid(1:2) + ( v(1:2,ip1) + v(1:2,i) ) * temp

    end do

    area = area / 2.0e+00_real64

    if ( area == 0.0e+00_real64 ) then
      centroid(1:2) = v(1:2,1)
    else
      centroid(1:2) = centroid(1:2) / ( 6.0e+00_real64 * area )
    end if
  end

  subroutine polygon_centroid_2 ( n, v, centroid )

  !*****************************************************************************80
  !
  !! POLYGON_CENTROID_2 computes the centroid of a polygon.
  !
  !  Discussion:
  !
  !    The centroid is the area-weighted sum of the centroids of
  !    disjoint triangles that make up the polygon.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    12 July 2003
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Adrian Bowyer, John Woodwark,
  !    A Programmer's Geometry,
  !    Butterworths, 1983,
  !    ISBN: 0408012420.
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of vertices of the polygon.
  !
  !    Input, real(real64) V(2,N), the coordinates of the vertices.
  !
  !    Output, real(real64) CENTROID(2), the coordinates of the centroid.
  !

    integer(int32) n

    real(real64) area_polygon
    real(real64) area_triangle
    real(real64) centroid(2)
    integer(int32) i
    real(real64) triangle_area
    real(real64) v(2,n)

    area_polygon = 0.0e+00_real64
    centroid(1:2) = 0.0e+00_real64

    do i = 1, n - 2

      area_triangle = triangle_area ( &
        v(1,i),   v(2,i), &
        v(1,i+1), v(2,i+1), &
        v(1,n),   v(2,n) )

      area_polygon = area_polygon + area_triangle

      centroid(1:2) = centroid(1:2) + area_triangle &
        * ( v(1:2,i) + v(1:2,i+1) + v(1:2,n) ) / 3.0e+00_real64

    end do

    if ( area_polygon == 0.0e+00_real64 ) then
      centroid(1:2) = v(1:2,1)
    else
      centroid(1:2) = centroid(1:2) / area_polygon
    end if
  end

  subroutine polygon_contains_point ( n, v, p, inside )

  !*****************************************************************************80
  !
  !! POLYGON_CONTAINS_POINT finds if a point is inside a polygon.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    06 November 2016
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of nodes or vertices in 
  !    the polygon.  N must be at least 3.
  !
  !    Input, real(real64) V(2,N), the vertices of the polygon.
  !
  !    Input, real(real64) P(2), the coordinates of the point to be tested.
  !
  !    Output, logical INSIDE, is TRUE if the point is inside 
  !    the polygon.
  !

    integer(int32) n

    integer(int32) i
    logical inside
    integer(int32) ip1
    real(real64) p(2)
    real(real64) px1
    real(real64) px2
    real(real64) py1
    real(real64) py2
    real(real64) v(2,n)
    real(real64) xints

    inside = .false.

    px1 = v(1,1)
    py1 = v(2,1)
    xints = p(1) - 1.0e+00_real64

    do i = 1, n

      px2 = v(1,mod(i,n)+1)
      py2 = v(2,mod(i,n)+1)

      if ( min ( py1, py2 ) < p(2) ) then
        if ( p(2) <= max ( py1, py2 ) ) then
          if ( p(1) <= max ( px1, px2 ) ) then
            if ( py1 /= py2 ) then
              xints = ( p(2) - py1 ) * ( px2 - px1 ) / ( py2 - py1 ) + px1
            end if
            if ( px1 == px2 .or. p(1) <= xints ) then
              inside = .not. inside
            end if
          end if
        end if
      end if

      px1 = px2
      py1 = py2

    end do
  end

  subroutine polygon_contains_point_2 ( n, v, p, inside )

  !*****************************************************************************80
  !
  !! POLYGON_CONTAINS_POINT_2: is a point inside a convex polygon.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    02 February 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of nodes or vertices in the 
  !    polygon.  N must be at least 3.
  !
  !    Input, real(real64) V(2,N), the vertices of the polygon.
  !
  !    Input, real(real64) P(2), the coordinates of the point to be tested.
  !
  !    Output, logical INSIDE, is TRUE if the point is inside
  !    the polygon or on its boundary.
  !

    integer(int32) n

    integer(int32) i
    logical inside
    real(real64) p(2)
    real(real64) t(2,3)
    real(real64) v(2,n)

    inside = .false.
  !
  !  A point is inside a convex polygon if and only if it is inside
  !  one of the triangles formed by X(1),Y(1) and any two consecutive
  !  points on the polygon's circumference.
  !
    t(1:2,1) = v(1:2,1)

    do i = 2, n - 1

      t(1:2,2) = v(1:2,i)
      t(1:2,3) = v(1:2,i+1)

      call triangle_contains_point_1 ( t, p, inside )

      if ( inside ) then
      end if

    end do
  end

  subroutine polygon_contains_point_3 ( n, v, p, inside )

  !*****************************************************************************80
  !
  !! POLYGON_CONTAINS_POINT_3 finds if a point is inside a simple polygon.
  !
  !  Discussion:
  !
  !    A simple polygon is one whose boundary never crosses itself.
  !    The polygon does not need to be convex.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    09 May 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Moshe Shimrat,
  !    ACM Algorithm 112,
  !    Position of Point Relative to Polygon,
  !    Communications of the ACM,
  !    Volume 5, Number 8, page 434, August 1962.
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of nodes or vertices in 
  !    the polygon.  N must be at least 3.
  !
  !    Input, real(real64) V(2,N), the vertices of the polygon.
  !
  !    Input, real(real64) P(2), the coordinates of the point to be tested.
  !
  !    Output, logical INSIDE, is TRUE if the point is inside 
  !    the polygon.
  !

    integer(int32) n

    integer(int32) i
    logical inside
    integer(int32) ip1
    real(real64) p(2)
    real(real64) v(2,n)

    inside = .false.

    do i = 1, n

      if ( i < n ) then
        ip1 = i + 1
      else
        ip1 = 1
      end if

      if ( ( v(2,i)   <  p(2) .and. p(2) <= v(2,ip1)   ) .or. &
           ( p(2) <= v(2,i)   .and. v(2,ip1)   < p(2) ) ) then
        if ( ( p(1) - v(1,i) ) - ( p(2) - v(2,i) ) &
           * ( v(1,ip1) - v(1,i) ) / ( v(2,ip1) - v(2,i) ) < 0.0e+00_real64 ) then
          inside = .not. inside
        end if
      end if

    end do
  end

  subroutine polygon_diameter ( n, v, diameter )

  !*****************************************************************************80
  !
  !! POLYGON_DIAMETER computes the diameter of a polygon.
  !
  !  Discussion:
  !
  !    The diameter of a polygon is the maximum distance between any
  !    two points on the polygon.  It is guaranteed that this maximum
  !    distance occurs between two vertices of the polygon.  It is
  !    sufficient to check the distance between all pairs of vertices.
  !    This is an N^2 algorithm.  There is an algorithm by Shamos which
  !    can compute this quantity in order N time instead.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    03 January 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of vertices of the polygon.
  !
  !    Input, real(real64) V(2,N), the vertices.
  !
  !    Output, real(real64) DIAMETER, the diameter of the polygon.
  !

    integer(int32) n

    real(real64) diameter
    integer(int32) i
    integer(int32) j
    real(real64) v(2,n)

    diameter = 0.0e+00_real64

    do i = 1, n

      do j = i + 1, n
        diameter = max ( diameter, &
          sqrt ( ( v(1,i) - v(1,j) ) ** 2 + ( v(2,i) - v(2,j) ) ** 2 ) )
      end do

    end do
  end

  subroutine polygon_expand ( n, v, h, w )

  !*****************************************************************************80
  !
  !! POLYGON_EXPAND expands a polygon.
  !
  !  Discussion:
  !
  !    This routine simple moves each vertex of the polygon outwards
  !    in such a way that the sides of the polygon advance by H.  
  !
  !    This approach should always work if the polygon is convex, or 
  !    star-shaped.  But for general polygons, it is possible
  !    that this procedure, for large enough H, will create a polygon
  !    whose sides intersect.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    03 March 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of sides of the polygon.
  !
  !    Input, real(real64) V(2,N), the coordinates of the vertices.
  !
  !    Input, real(real64) H, the expansion amount.
  !
  !    Output, real(real64) W(2,N), the "expanded" coordinates.
  !

    integer(int32) n

    real(real64) angle
    real(real64) angle_rad
    real(real64) h
    real(real64) h2
    integer(int32) i
    integer(int32) i4_wrap
    integer(int32) im1
    integer(int32) ip1
    real(real64) p4(2)
    real(real64) v(2,n)
    real(real64) w(2,n)
  !
  !  Consider each angle, formed by the nodes P(I-1), P(I), P(I+1).
  !
    do i = 1, n

      im1 = i4_wrap ( i-1, 1, n )
      ip1 = i4_wrap ( i+1, 1, n )
  !
  !        P1
  !        /
  !       /   P4
  !      /  .  
  !     / .
  !    P2--------->P3
  !
      call angle_half ( v(1:2,im1), v(1:2,i), v(1:2,ip1), p4 )
  !
  !  Compute the value of the half angle.
  !
      angle = angle_rad ( v(1:2,im1), v(1:2,i), p4(1:2) )
  !
  !  The stepsize along the ray must be adjusted so that the sides
  !  move out by H.
  !
      h2 = h / sin ( angle )

      w(1:2,i) = v(1:2,i) - h2 * ( p4(1:2) - v(1:2,i) )

    end do
  end

  subroutine polygon_inrad_data ( n, radin, area, radout, side )

  !*****************************************************************************80
  !
  !! POLYGON_INRAD_DATA determines polygonal data from its inner radius.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    09 September 2003
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of sides of the polygon.
  !    N must be at least 3.
  !
  !    Input, real(real64) RADIN, the inner radius of the polygon, that is,
  !    the radius of the largest circle that can be inscribed within
  !    the polygon.
  !
  !    Output, real(real64) AREA, the area of the regular polygon.
  !
  !    Output, real(real64) RADOUT, the outer radius of the polygon, that is,
  !    the radius of the smallest circle that can be described about
  !    the polygon.
  !
  !    Output, real(real64) SIDE, the length of one side of the polygon.
  !

    real(real64) angle
    real(real64) area
    integer(int32) n
    real(real64), parameter :: r8_pi = 3.141592653589793e+00_real64
    real(real64) radin
    real(real64) radout
    real(real64) side

    if ( n < 3 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLYGON_INRAD_DATA - Fatal error!'
      write ( *, '(a)' ) '  Input value of N must be at least 3'
      write ( *, '(a,i8)' ) '  but your input value was N = ', n
      stop 1
    end if

    angle = r8_pi / real ( n, real64)
    area = real ( n, real64) * radin * radin * tan ( angle )
    side = 2.0e+00_real64 * radin * tan ( angle )
    radout = 0.5e+00_real64 * side / sin ( angle )
  end

  subroutine polygon_integral_1 ( n, v, result )

  !*****************************************************************************80
  !
  !! POLYGON_INTEGRAL_1 integrates the function 1 over a polygon.
  !
  !  Discussion:
  !
  !    The polygon is bounded by the points (X(1:N), Y(1:N)).
  !
  !    INTEGRAL = 0.5 * sum ( 1 <= I <= N )
  !      ( X(I) + X(I-1) ) * ( Y(I) - Y(I-1) )
  !
  !    where X(0) and Y(0) should be replaced by X(N) and Y(N).
  !
  !    The integral of 1 over a polygon is the area of the polygon.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    02 February 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    SF Bockman,
  !    Generalizing the Formula for Areas of Polygons to Moments,
  !    American Mathematical Society Monthly,
  !    1989, pages 131-132.
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of vertices of the polygon.
  !    N should be at least 3 for a nonzero result.
  !
  !    Input, real(real64) V(2,N), the coordinates of the vertices
  !    of the polygon.  These vertices should be given in counter clockwise order.
  !
  !    Output, real(real64) RESULT, the value of the integral.
  !

    integer(int32) n

    integer(int32) i
    integer(int32) im1
    real(real64) result
    real(real64) v(2,n)

    result = 0.0e+00_real64

    if ( n < 3 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLYGON_INTEGRAL_1 - Fatal error!'
      write ( *, '(a)' ) '  The number of vertices must be at least 3.'
      write ( *, '(a,i8)' ) '  The input value of N = ', n
      stop 1
    end if

    do i = 1, n

      if ( i == 1 ) then
        im1 = n
      else
        im1 = i - 1
      end if

      result = result + 0.5e+00_real64 * ( v(1,i) + v(1,im1) ) * ( v(2,i) - v(2,im1) )

    end do
  end

  subroutine polygon_integral_x ( n, v, result )

  !*****************************************************************************80
  !
  !! POLYGON_INTEGRAL_X integrates the function X over a polygon.
  !
  !  Discussion:
  !
  !    The polygon is bounded by the points (X(1:N), Y(1:N)).
  !
  !    INTEGRAL = (1/6) * sum ( 1 <= I <= N )
  !      ( X(I)*X(I) + X(I) * X(I-1) + X(I-1)*X(I-1) ) * ( Y(I) - Y(I-1) )
  !
  !    where X(0) and Y(0) should be replaced by X(N) and Y(N).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    10 July 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    SF Bockman,
  !    Generalizing the Formula for Areas of Polygons to Moments,
  !    American Mathematical Society Monthly,
  !    1989, pages 131-132.
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of vertices of the polygon.
  !    N should be at least 3 for a nonzero result.
  !
  !    Input, real(real64) V(2,N), the coordinates of the vertices
  !    of the polygon.  These vertices should be given in counter clockwise order.
  !
  !    Output, real(real64) RESULT, the value of the integral.
  !

    integer(int32) n

    integer(int32) i
    integer(int32) im1
    real(real64) result
    real(real64) v(2,n)

    result = 0.0e+00_real64

    if ( n < 3 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLYGON_INTEGRAL_X - Fatal error!'
      write ( *, '(a)' ) '  The number of vertices must be at least 3.'
      write ( *, '(a,i8)' ) '  The input value of N = ', n
      stop 1
    end if

    do i = 1, n

      if ( i == 1 ) then
        im1 = n
      else
        im1 = i - 1
      end if

      result = result + ( v(1,i)**2 + v(1,i) * v(1,im1) + v(1,im1)**2 ) &
        * ( v(2,i) - v(2,im1) )

    end do

    result = result / 6.0e+00_real64
  end

  subroutine polygon_integral_xx ( n, v, result )

  !*****************************************************************************80
  !
  !! POLYGON_INTEGRAL_XX integrates the function X*X over a polygon.
  !
  !  Discussion:
  !
  !    The polygon is bounded by the points (X(1:N), Y(1:N)).
  !
  !    INTEGRAL = (1/12) * sum ( 1 <= I <= N )
  !      ( X(I)^3 + X(I)^2 * X(I-1) + X(I) * X(I-1)^2 + X(I-1)^3 )
  !      * ( Y(I) - Y(I-1) )
  !
  !    where X(0) and Y(0) should be replaced by X(N) and Y(N).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    10 July 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    SF Bockman,
  !    Generalizing the Formula for Areas of Polygons to Moments,
  !    American Mathematical Society Monthly,
  !    1989, pages 131-132.
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of vertices of the polygon.
  !    N should be at least 3 for a nonzero result.
  !
  !    Input, real(real64) V(2,N), the coordinates of the vertices
  !    of the polygon.  These vertices should be given in
  !    counter clockwise order.
  !
  !    Output, real(real64) RESULT, the value of the integral.
  !

    integer(int32) n

    integer(int32) i
    integer(int32) im1
    real(real64) result
    real(real64) v(2,n)

    result = 0.0e+00_real64

    if ( n < 3 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLYGON_INTEGRAL_XX - Fatal error!'
      write ( *, '(a)' ) '  The number of vertices must be at least 3.'
      write ( *, '(a,i8)' ) '  The input value of N = ', n
      stop 1
    end if

    do i = 1, n

      if ( i == 1 ) then
        im1 = n
      else
        im1 = i - 1
      end if

      result = result + ( v(1,i) ** 3 + v(1,i) ** 2 * v(1,im1) &
        + v(1,i) * v(1,im1) ** 2 + v(1,im1) ** 3 ) * ( v(2,i) - v(2,im1) )

    end do

    result = result / 12.0e+00_real64
  end

  subroutine polygon_integral_xy ( n, v, result )

  !*****************************************************************************80
  !
  !! POLYGON_INTEGRAL_XY integrates the function X*Y over a polygon.
  !
  !  Discussion:
  !
  !    The polygon is bounded by the points (X(1:N), Y(1:N)).
  !
  !    INTEGRAL = (1/24) * sum ( 1 <= I <= N )
  !      ( Y(I)   * ( 3 * X(I)^2 + 2 * X(I) * X(I-1) +     X(I-1)^2 )
  !      + Y(I-1) * (     X(I)^2 + 2 * X(I) * X(I-1) + 3 * X(I-1)^2 ) )
  !      * ( Y(I) - Y(I-1) )
  !
  !    where X(0) and Y(0) should be replaced by X(N) and Y(N).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    10 July 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    SF Bockman,
  !    Generalizing the Formula for Areas of Polygons to Moments,
  !    American Mathematical Society Monthly,
  !    1989, pages 131-132.
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of vertices of the polygon.
  !    N should be at least 3 for a nonzero result.
  !
  !    Input, real(real64) V(2,N), the coordinates of the vertices
  !    of the polygon.  These vertices should be given in
  !    counter clockwise order.
  !
  !    Output, real(real64) RESULT, the value of the integral.
  !

    integer(int32) n

    integer(int32) i
    integer(int32) im1
    real(real64) result
    real(real64) v(2,n)

    result = 0.0e+00_real64

    if ( n < 3 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLYGON_INTEGRAL_XY - Fatal error!'
      write ( *, '(a)' ) '  The number of vertices must be at least 3.'
      write ( *, '(a,i8)' ) '  The input value of N = ', n
      stop 1
    end if

    do i = 1, n

      if ( i == 1 ) then
        im1 = n
      else
        im1 = i - 1
      end if

      result = result + ( &
        v(2,i) * ( 3.0e+00_real64 * v(1,i)**2 + 2.0e+00_real64 * v(1,i) * v(1,im1) &
        + v(1,im1)**2 ) + v(2,im1) * ( v(1,i)**2 + 2.0e+00_real64 * v(1,i) * v(1,im1) &
        + 3.0e+00_real64 * v(1,im1)**2 ) ) * ( v(2,i) - v(2,im1) )

    end do

    result = result / 24.0e+00_real64
  end

  subroutine polygon_integral_y ( n, v, result )

  !*****************************************************************************80
  !
  !! POLYGON_INTEGRAL_Y integrates the function Y over a polygon.
  !
  !  Discussion:
  !
  !    The polygon is bounded by the points (X(1:N), Y(1:N)).
  !
  !    INTEGRAL = (1/6) * sum ( 1 <= I <= N )
  !      - ( Y(I)^2 + Y(I) * Y(I-1) + Y(I-1)^2 ) * ( X(I) - X(I-1) )
  !
  !    where X(0) and Y(0) should be replaced by X(N) and Y(N).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    10 July 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    SF Bockman,
  !    Generalizing the Formula for Areas of Polygons to Moments,
  !    American Mathematical Society Monthly,
  !    1989, pages 131-132.
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of vertices of the polygon.
  !    N should be at least 3 for a nonzero result.
  !
  !    Input, real(real64) V(2,N), the coordinates of the vertices
  !    of the polygon.  These vertices should be given in
  !    counter clockwise order.
  !
  !    Output, real(real64) RESULT, the value of the integral.
  !

    integer(int32) n

    integer(int32) i
    integer(int32) im1
    real(real64) result
    real(real64) v(2,n)

    result = 0.0e+00_real64

    if ( n < 3 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLYGON_INTEGRAL_Y - Fatal error!'
      write ( *, '(a)' ) '  The number of vertices must be at least 3.'
      write ( *, '(a,i8)' ) '  The input value of N = ', n
      stop 1
    end if

    do i = 1, n

      if ( i == 1 ) then
        im1 = n
      else
        im1 = i - 1
      end if

      result = result - ( v(2,i)**2 + v(2,i) * v(2,im1) + v(2,im1)**2 ) &
        * ( v(1,i) - v(1,im1) )

    end do

    result = result / 6.0e+00_real64
  end

  subroutine polygon_integral_yy ( n, v, result )

  !*****************************************************************************80
  !
  !! POLYGON_INTEGRAL_YY integrates the function Y*Y over a polygon.
  !
  !  Discussion:
  !
  !    The polygon is bounded by the points (X(1:N), Y(1:N)).
  !
  !    INTEGRAL = (1/12) * sum ( 1 <= I <= N )
  !      - ( Y(I)^3 + Y(I)^2 * Y(I-1) + Y(I) * Y(I-1)^2 + Y(I-1)^3 )
  !      * ( X(I) - X(I-1) )
  !
  !    where X(0) and Y(0) should be replaced by X(N) and Y(N).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    10 July 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    SF Bockman,
  !    Generalizing the Formula for Areas of Polygons to Moments,
  !    American Mathematical Society Monthly,
  !    1989, pages 131-132.
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of vertices of the polygon.
  !    N should be at least 3 for a nonzero result.
  !
  !    Input, real(real64) V(2,N), the coordinates of the vertices
  !    of the polygon.  These vertices should be given in
  !    counter clockwise order.
  !
  !    Output, real(real64) RESULT, the value of the integral.
  !

    integer(int32) n

    integer(int32) i
    integer(int32) im1
    real(real64) result
    real(real64) v(2,n)

    result = 0.0e+00_real64

    if ( n < 3 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLYGON_INTEGRAL_YY - Fatal error!'
      write ( *, '(a)' ) '  The number of polygonal vertices must be '
      write ( *, '(a,i8)' ) '  at least 3, but the input polygon has N = ', n
      stop 1
    end if

    do i = 1, n

      if ( i == 1 ) then
        im1 = n
      else
        im1 = i - 1
      end if

      result = result - ( v(2,i)**3 + v(2,i)**2 * v(2,im1) &
        + v(2,i) * v(2,im1)**2 + v(2,im1)**3 ) * ( v(1,i) - v(1,im1) )

    end do

    result = result / 12.0e+00_real64
  end

  function polygon_is_convex ( n, v )

  !*****************************************************************************80
  !
  !! POLYGON_IS_CONVEX determines whether a polygon is convex.
  !
  !  Discussion:
  !
  !    If the polygon has less than 3 distinct vertices, it is
  !    classified as convex degenerate.
  !
  !    If the polygon "goes around" more than once, it is classified
  !    as NOT convex.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    02 May 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Peter Schorn, Frederick Fisher,
  !    Testing the Convexity of a Polygon,
  !    in Graphics Gems IV, 
  !    edited by Paul Heckbert,
  !    AP Professional, 1994,
  !    T385.G6974.
  !
  !  Parameters
  !
  !    Input, integer(int32) N, the number of vertices.
  !
  !    Input/output, real(real64) V(2,N), the coordinates of the vertices 
  !    of the polygon.  On output, duplicate consecutive points have been 
  !    deleted, and the vertices have been reordered so that the 
  !    lexicographically least point comes first.
  !
  !    Output, integer(int32) POLYGON_IS_CONVEX:
  !    -1, the polygon is not convex;
  !     0, the polygon has less than 3 vertices; it is "degenerately" convex;
  !     1, the polygon is convex and counter clockwise;
  !     2, the polygon is convex and clockwise.
  !

    real(real64), parameter :: r8_pi = 3.141592653589793e+00_real64
    real(real64), parameter :: RAD_TO_DEG = 180.0e+00_real64 / r8_pi

    integer(int32) n

    real(real64) angle
    integer(int32), parameter :: CONVEX_CCW = 1
    integer(int32), parameter :: CONVEX_CW = 2
    real(real64) cross
    integer(int32), parameter :: DEGENERATE_CONVEX = 0
    real(real64) dot
    real(real64) exterior_total
    integer(int32) i
    integer(int32) ip1
    integer(int32) ip2
    integer(int32), parameter :: NOT_CONVEX = -1
    integer(int32) polygon_is_convex
    real(real64) sense
    real(real64), parameter :: tol = 1.0e+00_real64
    real(real64) v(2,n)

    exterior_total = 0.0e+00_real64
  !
  !  If there are not at least 3 distinct vertices, we are done.
  !
    if ( n < 3 ) then
      polygon_is_convex = DEGENERATE_CONVEX
    end if

    sense = 0.0e+00_real64
  !
  !  Consider each polygonal vertex I.
  !
    do i = 1, n

      ip1 = i + 1
      if ( n < ip1 ) then
        ip1 = ip1 - n
      end if

      ip2 = i + 2
      if ( n < ip2 ) then
        ip2 = ip2 - n
      end if

      dot =   ( v(1,ip2) - v(1,ip1) ) * ( v(1,i) - v(1,ip1) ) &
            + ( v(2,ip2) - v(2,ip1) ) * ( v(2,i) - v(2,ip1) )

      cross =   ( v(1,ip2) - v(1,ip1) ) * ( v(2,i) - v(2,ip1) ) &
              - ( v(1,i)   - v(1,ip1) ) * ( v(2,ip2) - v(2,ip1) )

      angle = atan2 ( cross, dot )
  !
  !  See if the turn defined by this vertex is our first indication of
  !  the "sense" of the polygon, or if it disagrees with the previously
  !  defined sense.
  !
      if ( sense == 0.0e+00_real64 ) then

        if ( angle < 0.0e+00_real64 ) then
          sense = -1.0e+00_real64
        else if ( 0.0e+00_real64 < angle ) then
          sense = +1.0e+00_real64
        end if

      else if ( sense == 1.0e+00_real64 ) then

        if ( angle < 0.0e+00_real64 ) then
          polygon_is_convex = NOT_CONVEX
        end if

      else if ( sense == -1.0e+00_real64 ) then

        if ( 0.0e+00_real64 < angle ) then
          polygon_is_convex = NOT_CONVEX
        end if

      end if
  !
  !  If the exterior total is greater than 360, then the polygon is
  !  going around again.
  !
      angle = atan2 ( -cross, -dot )

      exterior_total = exterior_total + angle

      if ( 360.0e+00_real64 + tol < abs ( exterior_total ) * RAD_TO_DEG ) then
        polygon_is_convex = NOT_CONVEX
      end if

    end do

    if ( sense == +1.0e+00_real64 ) then
      polygon_is_convex = CONVEX_CCW
    else if ( sense == -1.0e+00_real64 ) then
      polygon_is_convex = CONVEX_CW
    end if
  end

  subroutine polygon_lattice_area ( i, b, area )

  !*****************************************************************************80
  !
  !! POLYGON_LATTICE_AREA computes the area of a lattice polygon.
  !
  !  Discussion:
  !
  !    We define a lattice to be the 2D plane, in which the points
  !    whose (X,Y) coordinates are both integers are given a special
  !    status as "lattice points".
  !
  !    A lattice polygon is a polygon whose vertices are lattice points.
  !
  !    The area of a lattice polygon can be computed by Pick's Theorem:
  !
  !      Area = I + B / 2 - 1
  !
  !    where
  !
  !      I = the number of lattice points contained strictly inside the polygon;
  !
  !      B = the number of lattice points that lie exactly on the boundary.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    05 June 2002
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Branko Gruenbaum, Geoffrey Shephard,
  !    Pick's Theorem,
  !    The American Mathematical Monthly,
  !    Volume 100, Number 2, February 1993, pages 150-161.
  !
  !  Parameters:
  !
  !    Input, integer(int32) I, the number of interior lattice points.
  !
  !    Input, integer(int32) B, the number of boundary lattice points.
  !
  !    Output, real(real64) AREA, the area of the lattice polygon.
  !

    real(real64) area
    integer(int32) b
    integer(int32) i

    area = real ( i, real64) + real ( b, real64) / 2.0e+00_real64 - 1.0e+00_real64
  end

  subroutine polygon_outrad_data ( n, radout, area, radin, side )

  !*****************************************************************************80
  !
  !! POLYGON_OUTRAD_DATA determines polygonal data from its outer radius.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    09 September 2003
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of sides of the polygon.
  !    N must be at least 3.
  !
  !    Input, real(real64) RADOUT, the outer radius of the polygon, that is,
  !    the radius of the smallest circle that can be described
  !    around the polygon.
  !
  !    Output, real(real64) AREA, the area of the regular polygon.
  !
  !    Output, real(real64) RADIN, the inner radius of the polygon, that is,
  !    the radius of the largest circle that can be inscribed
  !    within the polygon.
  !
  !    Output, real(real64) SIDE, the length of one side of the polygon.
  !

    real(real64) angle
    real(real64) area
    integer(int32) n
    real(real64), parameter :: r8_pi = 3.141592653589793e+00_real64
    real(real64) radin
    real(real64) radout
    real(real64) side

    if ( n < 3 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLYGON_OUTRAD_DATA - Fatal error!'
      write ( *, '(a)' ) '  Input value of N must be at least 3'
      write ( *, '(a,i8)' ) '  but your input value was N = ', n
      stop 1
    end if

    angle = r8_pi / real ( n, real64)
    area = 0.5e+00_real64 * real ( n, real64) * radout * radout &
      * sin ( 2.0e+00_real64 * angle )
    side = 2.0e+00_real64 * radout * sin ( angle )
    radin = 0.5e+00_real64 * side / tan ( angle )
  end

  subroutine polygon_perimeter ( n, v, perimeter )

  !*****************************************************************************80
  !
  !! POLYGON_PERIMETER computes the perimeter of a polygon.
  !
  !  Discussion:
  !
  !    The perimeter is simply the sum of the side lengths.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    15 October 2015
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of vertices of the polygon.
  !
  !    Input, real(real64) V(2,N), the vertices.
  !
  !    Output, real(real64) PERIMETER, the perimeter.
  !

    integer(int32) n

    integer(int32) i
    integer(int32) im1
    real(real64) l
    real(real64) perimeter
    real(real64) v(2,n)

    perimeter = 0.0e+00_real64

    im1 = n

    do i = 1, n
      l = sqrt ( ( v(1,im1) - v(1,i) ) ** 2 + ( v(2,im1) - v(2,i) ) ** 2 )
      perimeter = perimeter + l
      im1 = i
    end do
  end

  subroutine polygon_perimeter_quad ( n, v, hmax, f, value )

  !*****************************************************************************80
  !
  !! POLYGON_PERIMETER_QUAD estimates an integral over the perimeter of a polygon.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    19 October 2015
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of vertices of the polygon.
  !
  !    Input, real(real64) V(2,N), the vertices.
  !
  !    Input, real(real64) HMAX, the maximum length of a quadrature interval.
  !
  !    Input, real(real64), external F ( X, Y ), a function whose integral 
  !    over the perimeter is desired.
  !
  !    Output, real(real64) VALUE, the estimated integral.
  !

    integer(int32) n

    real(real64) dxy
    real(real64), external :: f
    real(real64) fxy
    real(real64) hmax
    integer(int32) i
    integer(int32) i4_wrap
    integer(int32) ip1
    integer(int32) j
    real(real64) l
    integer(int32) m
    real(real64) v(2,n)
    real(real64) value
    real(real64) x
    real(real64) y

    value = 0.0e+00_real64

    do i = 1, n

      ip1 = i4_wrap ( i + 1, 1, n )
      l = sqrt ( ( v(1,ip1) - v(1,i) ) ** 2 + ( v(2,ip1) - v(2,i) ) ** 2 )
      m = ceiling ( l / hmax )
      dxy = l / real ( m, real64)

      do j = 1, 2 * m - 1, 2
        x = ( real ( 2 * m - j, real64) * v(1,i) &
            + real (         j, real64) * v(1,ip1) ) &
            / real ( 2 * m, real64)
        y = ( real ( 2 * m - j, real64) * v(2,i) &
            + real (         j, real64) * v(2,ip1) ) &
            / real ( 2 * m, real64)
        fxy = f ( x, y )
        value = value + fxy * dxy
      end do

    end do
  end

  subroutine polygon_point_dist ( n, v, p, dist )

  !*****************************************************************************80
  !
  !! POLYGON_POINT_DIST: distance ( polygon, point ).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    28 February 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of vertices.
  !
  !    Input, real(real64) V(2,N), the triangle vertices.
  !
  !    Input, real(real64) P(2), the point to be checked.
  !
  !    Output, real(real64) DIST, the distance from the point to the
  !    polygon.
  !

    integer(int32) n

    real(real64) dist
    real(real64) dist2
    integer(int32) i4_wrap
    integer(int32) j
    integer(int32) jp1
    real(real64) p(2)
    real(real64) v(2,n)
  !
  !  Find the distance to each of the line segments.
  !
    dist = huge ( dist )

    do j = 1, n

      jp1 = i4_wrap ( j+1, 1, n )

      call segment_point_dist ( v(1:2,j), v(1:2,jp1), p, dist2 )

      if ( dist2 < dist ) then
        dist = dist2
      end if

    end do
  end

  subroutine polygon_point_near ( n, v, p, pn, dist )

  !*****************************************************************************80
  !
  !! POLYGON_POINT_NEAR computes the nearest point on a polygon.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    28 February 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real(real64) V(2,N), the polygon vertices.
  !
  !    Input, real(real64) P(2), the point whose nearest polygon point
  !    is to be determined.
  !
  !    Output, real(real64) PN(2), the nearest point to P.
  !
  !    Output, real(real64) DIST, the distance from the point to the
  !    polygon.
  !

    integer(int32) n

    real(real64) dist
    real(real64) dist2
    integer(int32) i4_wrap
    integer(int32) j
    integer(int32) jp1
    real(real64) p(2)
    real(real64) pn(2)
    real(real64) pn2(2)
    real(real64) tval
    real(real64) v(2,n)
  !
  !  Find the distance to each of the line segments that make up the edges
  !  of the polygon.
  !
    dist = huge ( dist )
    pn(1:2) = 0.0e+00_real64

    do j = 1, n

      jp1 = i4_wrap ( j+1, 1, n )

      call segment_point_near ( v(1:2,j), v(1:2,jp1), p, &
        pn2, dist2, tval )

      if ( dist2 < dist ) then
        dist = dist2
        pn(1:2) = pn2(1:2)
      end if

    end do
  end

  subroutine polygon_sample ( nv, v, n, seed, s )

  !*****************************************************************************80
  !
  !! POLYGON_SAMPLE uniformly samples a polygon.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 August 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(int32) NV, the number of vertices.
  !
  !    Input, real(real64) V(2,NV), the vertices of the polygon, listed in
  !    counterclockwise order.
  !
  !    Input, integer(int32) N, the number of points to create.
  !
  !    Input/output, integer(int32) SEED, a seed for the random
  !    number generator.
  !
  !    Output, real(real64) S(2,N), the points.
  !

    integer(int32) n
    integer(int32) nv

    real(real64) area_cumulative(nv-2)
    real(real64) area_polygon
    real(real64) area_relative(nv-2)
    real(real64) area_triangle(nv-2)
    real(real64) area_percent
    integer(int32) i
    integer(int32) ip1
    integer(int32) j
    integer(int32) k
    real(real64) r(2)
    real(real64) r8_uniform_01
    integer(int32) seed
    real(real64) triangle_area
    integer(int32) triangles(3,nv-2)
    real(real64) s(2,n)
    real(real64) v(2,nv)
  !
  !  Triangulate the polygon.
  !
    call polygon_triangulate ( nv, v(1,1:nv), v(2,1:nv), triangles )
  !
  !  Determine the areas of each triangle.
  !
    do i = 1, nv - 2
      area_triangle(i) = triangle_area ( &
        v(1,triangles(1,i)), v(2,triangles(1,i)), &
        v(1,triangles(2,i)), v(2,triangles(2,i)), &
        v(1,triangles(3,i)), v(2,triangles(3,i)) )
    end do
  !
  !  Normalize the areas.
  !
    area_polygon = sum ( area_triangle(1:nv-2) )
    area_relative(1:nv-2) = area_triangle(1:nv-2) / area_polygon
  !
  !  Replace each area by the sum of itself and all previous ones.
  !
    area_cumulative(1) = area_relative(1)
    do i = 2, nv - 2
      area_cumulative(i) = area_relative(i) + area_cumulative(i-1)
    end do

    do j = 1, n
  !
  !  Choose triangle I at random, based on areas.
  !
      area_percent = r8_uniform_01 ( seed )

      do k = 1, nv - 2

        i = k

        if ( area_percent <= area_cumulative(k) ) then
          exit
        end if

      end do
  !
  !  Now choose a point at random in triangle I.
  !
      call r8vec_uniform_01 ( 2, seed, r )

      if ( 1.0e+00_real64 < sum ( r(1:2) ) ) then
        r(1:2) = 1.0e+00_real64 - r(1:2)
      end if

      s(1:2,j) = ( 1.0e+00_real64 - r(1) - r(2) ) * v(1:2,triangles(1,i)) &
                           + r(1)          * v(1:2,triangles(2,i)) &
                                  + r(2)   * v(1:2,triangles(3,i))
    end do
  end

  subroutine polygon_side_data ( n, side, area, radin, radout )

  !*****************************************************************************80
  !
  !! POLYGON_SIDE_DATA determines polygonal data from its side length.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    29 June 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of sides of the polygon.
  !    N must be at least 3.
  !
  !    Input, real(real64) SIDE, the length of one side of the polygon.
  !
  !    Output, real(real64) AREA, the area of the regular polygon.
  !
  !    Output, real(real64) RADIN, the inner radius of the polygon, that is,
  !    the radius of the largest circle that can be inscribed within
  !    the polygon.
  !
  !    Output, real(real64) RADOUT, the outer radius of the polygon, that is,
  !    the radius of the smallest circle that can be described about
  !    the polygon.
  !

    real(real64) angle
    real(real64) area
    integer(int32) n
    real(real64), parameter :: r8_pi = 3.141592653589793e+00_real64
    real(real64) radin
    real(real64) radout
    real(real64) side

    if ( n < 3 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLYGON_SIDE_DATA - Fatal error!'
      write ( *, '(a)' ) '  Input value of N must be at least 3'
      write ( *, '(a,i8)' ) '  but your input value was N = ', n
      stop 1
    end if

    angle = r8_pi / real ( n, real64)
    area = 0.25e+00_real64 * real ( n, real64) * side * side / tan ( angle )
    radin = 0.5e+00_real64 * side / tan ( angle )
    radout = 0.5e+00_real64 * side / sin ( angle )
  end

  subroutine polygon_triangulate ( n, x, y, triangles )

  !*****************************************************************************80
  !
  !! POLYGON_TRIANGULATE determines a triangulation of a polygon.
  !
  !  Discussion:
  !
  !    There are N-3 triangles in the triangulation.
  !
  !    For the first N-2 triangles, the first edge listed is always an
  !    internal diagonal.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 May 2014
  !
  !  Author:
  !
  !    Original C version by Joseph ORourke.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Joseph ORourke,
  !    Computational Geometry in C,
  !    Cambridge, 1998,
  !    ISBN: 0521649765,
  !    LC: QA448.e38_real64.
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of vertices.
  !
  !    Input, real(real64) X(N), Y(N), the coordinates of each vertex.
  !
  !    Output, integer(int32) TRIANGLES(3,N-2), the triangles of the 
  !    triangulation.
  !

    integer(int32) n

    real(real64) area
    logical diagonal
    logical ear(n)
    integer(int32) first
    integer(int32) i
    integer(int32) i0
    integer(int32) i1
    integer(int32) i2
    integer(int32) i3
    integer(int32) i4
    integer(int32) next(n)
    integer(int32) node
    integer(int32) node_m1
    integer(int32) prev(n)
    integer(int32) triangle_num
    integer(int32) triangles(3,n-2)
    real(real64) x(n)
    real(real64) y(n)
  !
  !  We must have at least 3 vertices.
  !
    if ( n < 3 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLYGON_TRIANGULATE - Fatal error!'
      write ( *, '(a)' ) '  N < 3.'
      stop 1
    end if
  !
  !  Consecutive vertices cannot be equal.
  !
    node_m1 = n
    do node = 1, n
      if ( x(node_m1) == x(node) .and. y(node_m1) == y(node) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'POLYGON_TRIANGULATE - Fatal error!'
        write ( *, '(a)' ) '  Two consecutive nodes are identical.'
        stop 1
      end if
      node_m1 = node
    end do
  !
  !  Area must be positive.
  !
    area = 0.0e+00_real64
    do node = 1, n - 2
      area = area + 0.5e+00_real64 * &
      ( &
          ( x(node+1) - x(node) ) * ( y(node+2) - y(node) ) &
        - ( x(node+2) - x(node) ) * ( y(node+1) - y(node) ) &
      )
    end do

    if ( area <= 0.0e+00_real64 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLYGON_TRIANGULATE - Fatal error!'
      write ( *, '(a)' ) '  Polygon has zero or negative area.'
      stop 1
    end if
  !
  !  PREV and NEXT point to the previous and next nodes.
  !
    i = 1
    prev(i) = n
    next(i) = i + 1

    do i = 2, n - 1
      prev(i) = i - 1
      next(i) = i + 1
    end do

    i = n
    prev(i) = i - 1
    next(i) = 1
  !
  !  EAR indicates whether the node and its immediate neighbors form an ear
  !  that can be sliced off immediately.
  !
    do i = 1, n
      ear(i) = diagonal ( prev(i), next(i), n, prev, next, x, y )
    end do

    triangle_num = 0

    i2 = 1

    do while ( triangle_num < n - 3 )
  !
  !  If I2 is an ear, gather information necessary to carry out
  !  the slicing operation and subsequent "healing".
  !
      if ( ear(i2) ) then

        i3 = next(i2)
        i4 = next(i3)
        i1 = prev(i2)
        i0 = prev(i1)
  !
  !  Make vertex I2 disappear.
  !
        next(i1) = i3
        prev(i3) = i1
  !
  !  Update the earity of I1 and I3, because I2 disappeared.
  !
        ear(i1) = diagonal ( i0, i3, n, prev, next, x, y )
        ear(i3) = diagonal ( i1, i4, n, prev, next, x, y )
  !
  !  Add the diagonal [I3, I1, I2] to the list.
  !
        triangle_num = triangle_num + 1
        triangles(1,triangle_num) = i3
        triangles(2,triangle_num) = i1
        triangles(3,triangle_num) = i2

      end if
  !
  !  Try the next vertex.
  !
      i2 = next(i2)

    end do
  !
  !  The last triangle is formed from the three remaining vertices.
  !
    i3 = next(i2)
    i1 = prev(i2)

    triangle_num = triangle_num + 1
    triangles(1,triangle_num) = i3
    triangles(2,triangle_num) = i1
    triangles(3,triangle_num) = i2
  end

  function r8_degrees ( radians )

  !*****************************************************************************80
  !
  !! R8_DEGREES converts an angle from radian to degree measure.
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
  !    Input, real(real64) RADIANS, the angle measurement in radians.
  !
  !    Output, real(real64) R8_DEGREES, the angle measurement in degrees.
  !

    real(real64) r8_degrees
    real(real64), parameter :: r8_pi = 3.1415926535897932384626434e+00_real64
    real(real64) radians

    r8_degrees = radians * 180.0e+00_real64 / r8_pi
  end

  function r8_uniform_01 ( seed )

  !*****************************************************************************80
  !
  !! R8_UNIFORM_01 returns a unit pseudorandom R8.
  !
  !  Discussion:
  !
  !    An R8 is a real(real64) value.
  !
  !    For now, the input quantity SEED is an integer variable.
  !
  !    This routine implements the recursion
  !
  !      seed = 16807 * seed mod ( 2^31 - 1 )
  !      r8_uniform_01 = seed / ( 2^31 - 1 )
  !
  !    The integer arithmetic never requires more than 32 bits,
  !    including a sign bit.
  !
  !    If the initial seed is 12345, then the first three computations are
  !
  !      Input     Output      R8_UNIFORM_01
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
  !    05 July 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Paul Bratley, Bennett Fox, Linus Schrage,
  !    A Guide to Simulation,
  !    Springer Verlag, pages 201-202, 1983.
  !
  !    Pierre L'Ecuyer,
  !    Random Number Generation,
  !    in Handbook of Simulation,
  !    edited by Jerry Banks,
  !    Wiley Interscience, page 95, 1998.
  !
  !    Bennett Fox,
  !    Algorithm 647:
  !    Implementation and Relative Efficiency of Quasirandom
  !    Sequence Generators,
  !    ACM Transactions on Mathematical Software,
  !    Volume 12, Number 4, pages 362-376, 1986.
  !
  !    Peter Lewis, Allen Goodman, James Miller
  !    A Pseudo-Random Number Generator for the System/360,
  !    IBM Systems Journal,
  !    Volume 8, pages 136-143, 1969.
  !
  !  Parameters:
  !
  !    Input/output, integer(int32) SEED, the "seed" value, which should
  !    NOT be 0. On output, SEED has been updated.
  !
  !    Output, real(real64) R8_UNIFORM_01, a new pseudorandom variate,
  !    strictly between 0 and 1.
  !

    integer(int32), parameter :: i4_huge = 2147483647
    integer(int32) k
    real(real64) r8_uniform_01
    integer(int32) seed

    if ( seed == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
      write ( *, '(a)' ) '  Input value of SEED = 0.'
      stop 1
    end if

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r8_uniform_01 = real ( seed, real64) * 4.656612875e-10_real64
  end

  subroutine r8mat_solve ( n, rhs_num, a, info )

  !*****************************************************************************80
  !
  !! R8MAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
  !
  !  Discussion:
  !
  !    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    06 August 2009
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the order of the matrix.
  !
  !    Input, integer(int32) RHS_NUM, the number of right hand sides.
  !    RHS_NUM must be at least 0.
  !
  !    Input/output, real(real64) A(N,N+RHS_NUM), contains in rows and
  !    columns 1 to N the coefficient matrix, and in columns N+1 through
  !    N+RHS_NUM, the right hand sides.  On output, the coefficient matrix
  !    area has been destroyed, while the right hand sides have
  !    been overwritten with the corresponding solutions.
  !
  !    Output, integer(int32) INFO, singularity flag.
  !    0, the matrix was not singular, the solutions were computed;
  !    J, factorization failed on step J, and the solutions could not
  !    be computed.
  !

    integer(int32) n
    integer(int32) rhs_num

    real(real64) a(n,n+rhs_num)
    real(real64) apivot
    real(real64) factor
    integer(int32) i
    integer(int32) info
    integer(int32) ipivot
    integer(int32) j
    real(real64) t(n+rhs_num)

    info = 0

    do j = 1, n
  !
  !  Choose a pivot row.
  !
      ipivot = j
      apivot = a(j,j)

      do i = j + 1, n
        if ( abs ( apivot ) < abs ( a(i,j) ) ) then
          apivot = a(i,j)
          ipivot = i
        end if
      end do

      if ( apivot == 0.0e+00_real64 ) then
        info = j
      end if
  !
  !  The pivot row moves into the J-th row.
  !
      if ( ipivot /= j ) then
        t(       1:n+rhs_num) = a(ipivot,1:n+rhs_num)
        a(ipivot,1:n+rhs_num) = a(j,     1:n+rhs_num)
        a(j,     1:n+rhs_num) = t(       1:n+rhs_num)
      end if
  !
  !  A(J,J) becomes 1.
  !
      a(j,j) = 1.0e+00_real64
      a(j,j+1:n+rhs_num) = a(j,j+1:n+rhs_num) / apivot
  !
  !  A(I,J) becomes 0.
  !
      do i = 1, n

        if ( i /= j ) then
          factor = a(i,j)
          a(i,j) = 0.0e+00_real64
          a(i,j+1:n+rhs_num) = a(i,j+1:n+rhs_num) - factor * a(j,j+1:n+rhs_num)
        end if

      end do

    end do
  end

  subroutine r8vec_uniform_01 ( n, seed, r )

  !*****************************************************************************80
  !
  !! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
  !
  !  Discussion:
  !
  !    An R8VEC is a vector of R8's.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 July 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Paul Bratley, Bennett Fox, Linus Schrage,
  !    A Guide to Simulation,
  !    Springer Verlag, pages 201-202, 1983.
  !
  !    Bennett Fox,
  !    Algorithm 647:
  !    Implementation and Relative Efficiency of Quasirandom
  !    Sequence Generators,
  !    ACM Transactions on Mathematical Software,
  !    Volume 12, Number 4, pages 362-376, 1986.
  !
  !    Peter Lewis, Allen Goodman, James Miller
  !    A Pseudo-Random Number Generator for the System/360,
  !    IBM Systems Journal,
  !    Volume 8, pages 136-143, 1969.
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of entries in the vector.
  !
  !    Input/output, integer(int32) SEED, the "seed" value, which
  !    should NOT be 0.  On output, SEED has been updated.
  !
  !    Output, real(real64) R(N), the vector of pseudorandom values.
  !

    integer(int32) n

    integer(int32) i
    integer(int32) k
    integer(int32) seed
    real(real64) r(n)

    if ( seed == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
      write ( *, '(a)' ) '  Input value of SEED = 0.'
      stop 1
    end if

    do i = 1, n

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + 2147483647
      end if

      r(i) = real ( seed, real64) * 4.656612875e-10_real64

    end do
  end

  subroutine segment_point_dist ( p1, p2, p, dist )

  !*****************************************************************************80
  !
  !! SEGMENT_POINT_DIST: distance ( line segment, point ).
  !
  !  Discussion:
  !
  !    A line segment is the finite portion of a line that lies between
  !    two points P1 and P2.
  !
  !    The nearest point will satisfy the condition
  !
  !      PN = (1-T) * P1 + T * P2.
  !
  !    T will always be between 0 and 1.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    03 May 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real(real64) P1(2), P2(2), the endpoints of the line segment.
  !
  !    Input, real(real64) P(2), the point whose nearest neighbor on the line
  !    segment is to be determined.
  !
  !    Output, real(real64) DIST, the distance from the point to the
  !    line segment.
  !

    real(real64) bot
    real(real64) dist
    real(real64) p(2)
    real(real64) p1(2)
    real(real64) p2(2)
    real(real64) pn(2)
    real(real64) t
  !
  !  If the line segment is actually a point, then the answer is easy.
  !
    if ( all ( p1(1:2) == p2(1:2) ) ) then

      t = 0.0e+00_real64

    else

      bot = sum ( ( p2(1:2) - p1(1:2) )**2 )

      t = sum ( ( p(1:2)  - p1(1:2) ) &
              * ( p2(1:2) - p1(1:2) ) ) / bot

      t = max ( t, 0.0e+00_real64 )
      t = min ( t, 1.0e+00_real64 )

    end if

    pn(1:2) = p1(1:2) + t * ( p2(1:2) - p1(1:2) )

    dist = sqrt ( sum ( ( p(1:2) - pn(1:2) )**2 ) )
  end

  subroutine segment_point_near ( p1, p2, p, pn, dist, t )

  !*****************************************************************************80
  !
  !! SEGMENT_POINT_NEAR: nearest point on line segment to point.
  !
  !  Discussion:
  !
  !    A line segment is the finite portion of a line that lies between
  !    two points P1 and P2.
  !
  !    The nearest point will satisfy the condition
  !
  !      PN = (1-T) * P1 + T * P2.
  !
  !    T will always be between 0 and 1.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    03 May 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real(real64) P1(2), P2(2), the endpoints of the line segment.
  !
  !    Input, real(real64) P(2), the point whose nearest neighbor
  !    on the line segment is to be determined.
  !
  !    Output, real(real64) PN(2), the point on the line segment which is
  !    nearest the point P.
  !
  !    Output, real(real64) DIST, the distance from the point to the 
  !    nearest point on the line segment.
  !
  !    Output, real(real64) T, the relative position of the point PN
  !    to the points P1 and P2.
  !

    real(real64) bot
    real(real64) dist
    real(real64) p(2)
    real(real64) p1(2)
    real(real64) p2(2)
    real(real64) pn(2)
    real(real64) t
  !
  !  If the line segment is actually a point, then the answer is easy.
  !
    if ( all ( p1(1:2) == p2(1:2) ) ) then

      t = 0.0e+00_real64

    else

      bot = sum ( ( p2(1:2) - p1(1:2) )**2 )

      t = sum ( ( p(1:2)  - p1(1:2) ) &
              * ( p2(1:2) - p1(1:2) ) ) / bot

      t = max ( t, 0.0e+00_real64 )
      t = min ( t, 1.0e+00_real64 )

    end if

    pn(1:2) = p1(1:2) + t * ( p2(1:2) - p1(1:2) )

    dist = sqrt ( sum ( ( p(1:2) - pn(1:2) )**2 ) )
  end

  function triangle_area ( xa, ya, xb, yb, xc, yc )

  !*****************************************************************************80
  !
  !! TRIANGLE_AREA computes the signed area of a triangle.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 May 2014
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real(real64) XA, YA, XB, YB, XC, YC, the coordinates of
  !    the vertices of the triangle, given in counterclockwise order.
  !
  !    Output, real(real64) TRIANGLE_AREA, the signed area of the triangle.
  !

    real(real64) triangle_area
    real(real64) value
    real(real64) xa
    real(real64) xb
    real(real64) xc
    real(real64) ya
    real(real64) yb
    real(real64) yc

    value = 0.5e+00_real64 * ( &
        ( xb - xa ) * ( yc - ya ) &
      - ( xc - xa ) * ( yb - ya ) )

    triangle_area = value
  end

  subroutine triangle_barycentric ( t, p, xsi )

  !*****************************************************************************80
  !
  !! TRIANGLE_BARYCENTRIC finds the barycentric coordinates of a point.
  !
  !  Discussion:
  !
  !    The barycentric coordinate of point P related to vertex A can be
  !    interpreted as the ratio of the area of the triangle with 
  !    vertex A replaced by vertex P to the area of the original 
  !    triangle.
  !
  !    This routine assumes that the triangle vertices are given in
  !    counter clockwise order.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    28 December 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real(real64) T(2,3), the triangle vertices.
  !    The vertices should be given in counter clockwise order.
  !
  !    Input, real(real64) P(2), the point to be checked.
  !
  !    Output, real(real64) XSI(3), the barycentric coordinates of P
  !    with respect to the triangle.
  !

    real(real64) a(2,3)
    integer(int32) info
    real(real64) p(2)
    real(real64) t(2,3)
    real(real64) xsi(3)
  !
  !  Set up the linear system
  !
  !    ( X2-X1  X3-X1 ) XSI(1)  = X-X1
  !    ( Y2-Y1  Y3-Y1 ) XSI(2)    Y-Y1
  !
  !  which is satisfied by the barycentric coordinates of P.
  !
    a(1,1) = t(1,2) - t(1,1)
    a(1,2) = t(1,3) - t(1,1)
    a(1,3) = p(1)   - t(1,1)

    a(2,1) = t(2,2) - t(2,1)
    a(2,2) = t(2,3) - t(2,1)
    a(2,3) = p(2)   - t(2,1)
  !
  !  Solve the linear system.
  !
    call r8mat_solve ( 2, 1, a, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TRIANGLE_BARYCENTRIC - Fatal error!'
      write ( *, '(a)' ) '  The linear system is singular.'
      write ( *, '(a)' ) '  The input data does not form a proper triangle.'
      stop 1
    end if

    xsi(1) = a(1,3)
    xsi(2) = a(2,3)
    xsi(3) = 1.0e+00_real64 - xsi(1) - xsi(2)
  end

  subroutine triangle_contains_point_1 ( t, p, inside )

  !*****************************************************************************80
  !
  !! TRIANGLE_CONTAINS_POINT_1 finds if a point is inside a triangle.
  !
  !  Discussion:
  !
  !    It is conventional to list the triangle vertices in counter clockwise
  !    order.  However, this routine does not require a particular order
  !    for the vertices.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    16 June 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real(real64) T(2,3), the triangle vertices.
  !
  !    Input, real(real64) P(2), the point to be checked.
  !
  !    Output, logical INSIDE, is TRUE if the point is 
  !    inside the triangle.
  !

    logical inside
    real(real64) p(2)
    real(real64) t(2,3)
    real(real64) xsi(3)

    call triangle_barycentric ( t, p, xsi )

    if ( any ( xsi(1:3) < 0.0e+00_real64 ) ) then
      inside = .false.
    else
      inside = .true.
    end if
  end

end module polygon_properties_mod
