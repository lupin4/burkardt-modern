!> polygon_triangulate — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module polygon_triangulate_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: angle_degree, between, collinear, diagonal, diagonalie, in_cone
  public :: intersect, intersect_prop, l4_xor, polygon_area, polygon_triangulate, triangle_area

contains

  pure function angle_degree ( x1, y1, x2, y2, x3, y3 ) &
        bind(C, name="angle_degree")

  !*****************************************************************************80
  !
  !! ANGLE_DEGREE returns the degree angle defined by three points.
  !
  !  Discussion:
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
  !    28 August 2016
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real X1, Y1, X2, Y2, X3, Y3, the coordinates of the points
  !    P1, P2, P3.
  !
  !    Output, real VALUE, the angle swept out by the rays, measured
  !    in degrees.  0 <= VALUE < 360.  If either ray has zero length,
  !    then VALUE is set to 0.
  !

    real(dp) :: angle_degree
    real(dp), parameter :: r8_pi = 3.141592653589793_dp
    real(dp) :: value
    real(dp) :: x
    real(dp), intent(in), value :: x1
    real(dp), intent(in), value :: x2
    real(dp), intent(in), value :: x3
    real(dp) :: y
    real(dp), intent(in), value :: y1
    real(dp), intent(in), value :: y2
    real(dp), intent(in), value :: y3

    x = ( x3 - x2 ) * ( x1 - x2 ) + ( y3 - y2 ) * ( y1 - y2 )

    y = ( x3 - x2 ) * ( y1 - y2 ) - ( y3 - y2 ) * ( x1 - x2 )

    if ( x == 0.0_dp .and. y == 0.0_dp ) then

      value = 0.0_dp

    else

      value = atan2 ( y, x )

      if ( value < 0.0_dp ) then
        value = value + 2.0_dp * r8_pi
      end if

      value = 180.0_dp * value / r8_pi

    end if

    angle_degree = value
  end function angle_degree

  function between ( xa, ya, xb, yb, xc, yc ) &
        bind(C, name="between")

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
  !    LC: QA448.e38_dp.
  !
  !  Parameters:
  !
  !    Input, real(dp) XA, YA, XB, YB, XC, YC, the coordinates of
  !    the vertices.
  !
  !    Output, logical BETWEEN, is TRUE if C is between A and B.
  !

    logical :: between
    logical :: value
    real(dp), intent(in), value :: xa
    real(dp), intent(in), value :: xb
    real(dp), intent(in), value :: xc
    real(dp) :: xmax
    real(dp) :: xmin
    real(dp), intent(in), value :: ya
    real(dp), intent(in), value :: yb
    real(dp), intent(in), value :: yc
    real(dp) :: ymax
    real(dp) :: ymin

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
  end function between

  function collinear ( xa, ya, xb, yb, xc, yc ) &
        bind(C, name="collinear")

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
  !    10 September 2016
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
  !    LC: QA448.e38_dp.
  !
  !  Parameters:
  !
  !    Input, real(dp) XA, YA, XB, YB, XC, YC, the coordinates of
  !    the vertices.
  !
  !    Output, logical COLLINEAR, is TRUE if the points are judged
  !    to be collinear.
  !

    real(dp) :: area
    logical :: collinear
    real(dp), parameter :: r8_eps = 2.220446049250313e-016_dp
    real(dp) :: side_ab_sq
    real(dp) :: side_bc_sq
    real(dp) :: side_ca_sq
    real(dp) :: side_max_sq
    logical :: value
    real(dp), intent(in), value :: xa
    real(dp), intent(in), value :: xb
    real(dp), intent(in), value :: xc
    real(dp), intent(in), value :: ya
    real(dp), intent(in), value :: yb
    real(dp), intent(in), value :: yc

    area = triangle_area ( xa, ya, xb, yb, xc, yc )

    side_ab_sq = ( xa - xb ) ** 2 + ( ya - yb ) ** 2
    side_bc_sq = ( xb - xc ) ** 2 + ( yb - yc ) ** 2
    side_ca_sq = ( xc - xa ) ** 2 + ( yc - ya ) ** 2

    side_max_sq = max ( side_ab_sq, max ( side_bc_sq, side_ca_sq ) )

    if ( side_max_sq <= r8_eps ) then
      value = .true.
    else if ( 2.0_dp * abs ( area ) <= r8_eps * side_max_sq ) then
      value = .true.
    else
      value = .false.
    end if

    collinear = value
  end function collinear

  function diagonal ( im1, ip1, n, prev_node, next_node, x, y ) &
        bind(C, name="diagonal")

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
  !    LC: QA448.e38_dp.
  !
  !  Parameters:
  !
  !    Input, integer(ip) IM1, IP1, the indices of two vertices.
  !
  !    Input, integer(ip) N, the number of vertices.
  !
  !    Input, integer(ip) PREV_NODE(N), the previous neighbor of
  !    each vertex.
  !
  !    Input, integer(ip) NEXT_NODE(N), the next neighbor of each vertex.
  !
  !    Input, real(dp) X(N), Y(N), the coordinates of each vertex.
  !
  !    Output, logical DIAGONAL, the value of the test.
  !

    integer(ip), intent(in), value :: n

    logical :: diagonal
    integer(ip), intent(in), value :: im1
    integer(ip), intent(in), value :: ip1
    integer(ip), intent(in) :: next_node(n)
    integer(ip), intent(in) :: prev_node(n)
    logical :: value1
    logical :: value2
    logical :: value3
    real(dp), intent(in) :: x(n)
    real(dp), intent(in) :: y(n)

    value1 = in_cone ( im1, ip1, n, prev_node, next_node, x, y )
    value2 = in_cone ( ip1, im1, n, prev_node, next_node, x, y )
    value3 = diagonalie ( im1, ip1, n, next_node, x, y )

    diagonal = ( value1 .and. value2 .and. value3 )
  end function diagonal

  function diagonalie ( im1, ip1, n, next_node, x, y ) &
        bind(C, name="diagonalie")

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
  !    LC: QA448.e38_dp.
  !
  !  Parameters:
  !
  !    Input, integer(ip) IM1, IP1, the indices of two vertices.
  !
  !    Input, integer(ip) N, the number of vertices.
  !
  !    Input, integer(ip) NEXT_NODE(N), the next neighbor of each vertex.
  !
  !    Input, real(dp) X(N), Y(N), the coordinates of each vertex.
  !
  !    Output, logical DIAGONALIE, the value of the test.
  !

    integer(ip), intent(in), value :: n

    logical :: diagonalie
    integer(ip) :: first
    integer(ip), intent(in), value :: im1
    integer(ip), intent(in), value :: ip1
    integer(ip) :: j
    integer(ip) :: jp1
    integer(ip), intent(in) :: next_node(n)
    logical :: value
    logical :: value2
    real(dp), intent(in) :: x(n)
    real(dp), intent(in) :: y(n)

    first = im1
    j = first
    jp1 = next_node(first)

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
      jp1 = next_node(j)

      if ( j == first ) then
        exit
      end if

    end do

    diagonalie = value
  end function diagonalie

  function in_cone ( im1, ip1, n, prev_node, next_node, x, y ) &
        bind(C, name="in_cone")

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
  !    LC: QA448.e38_dp.
  !
  !  Parameters:
  !
  !    Input, integer(ip) IM1, IP1, the indices of two vertices.
  !
  !    Input, integer(ip) N, the number of vertices.
  !
  !    Input, integer(ip) PREV_NODE(N), the previous neighbor of
  !    each vertex.
  !
  !    Input, integer(ip) NEXT_NODE(N), the next neighbor of each vertex.
  !
  !    Input, real(dp) X(N), Y(N), the coordinates of each vertex.
  !
  !    Output, logical IN_CONE, the value of the test.
  !

    integer(ip), intent(in), value :: n

    integer(ip) :: i
    integer(ip), intent(in), value :: im1
    integer(ip) :: im2
    logical :: in_cone
    integer(ip), intent(in), value :: ip1
    integer(ip), intent(in) :: next_node(n)
    integer(ip), intent(in) :: prev_node(n)
    real(dp) :: t1
    real(dp) :: t2
    real(dp) :: t3
    real(dp) :: t4
    real(dp) :: t5
    logical :: value
    real(dp), intent(in) :: x(n)
    real(dp), intent(in) :: y(n)

    im2 = prev_node(im1)
    i = next_node(im1)

    t1 = triangle_area ( x(im1), y(im1), x(i), y(i), x(im2), y(im2) )

    if ( 0.0_dp <= t1 ) then

      t2 = triangle_area ( x(im1), y(im1), x(ip1), y(ip1), x(im2), y(im2) )
      t3 = triangle_area ( x(ip1), y(ip1), x(im1), y(im1), x(i), y(i) )
      value = ( ( 0.0_dp < t2 ) .and. ( 0.0_dp < t3 ) )

    else

      t4 = triangle_area ( x(im1), y(im1), x(ip1), y(ip1), x(i), y(i) )
      t5 = triangle_area ( x(ip1), y(ip1), x(im1), y(im1), x(im2), y(im2) )
      value = .not. ( ( 0.0_dp <= t4 ) .and. ( 0.0_dp <= t5 ) )

    end if

    in_cone = value
  end function in_cone

  function intersect ( xa, ya, xb, yb, xc, yc, xd, yd ) &
        bind(C, name="intersect")

  !*****************************************************************************80
  !
  !! INTERSECT is true if lines VA:VB and VC:VD intersect.
  !
  !  Discussion:
  !
  !    Thanks to Gene Dial for correcting the call to intersect_prop(),
  !    08 September 2016.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    08 September 2016
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
  !    LC: QA448.e38_dp.
  !
  !  Parameters:
  !
  !    Input, real(dp) XA, YA, XB, YB, XC, YC, XD, YD, the X and Y
  !    coordinates of the four vertices.
  !
  !    Output, logical VALUE, the value of the test.
  !

    logical :: intersect
    logical :: value
    real(dp), intent(in), value :: xa
    real(dp), intent(in), value :: xb
    real(dp), intent(in), value :: xc
    real(dp), intent(in), value :: xd
    real(dp), intent(in), value :: ya
    real(dp), intent(in), value :: yb
    real(dp), intent(in), value :: yc
    real(dp), intent(in), value :: yd

    if ( intersect_prop ( xa, ya, xb, yb, xc, yc, xd, yd ) ) then
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
  end function intersect

  function intersect_prop ( xa, ya, xb, yb, xc, yc, xd, yd ) &
        bind(C, name="intersect_prop")

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
  !    LC: QA448.e38_dp.
  !
  !  Parameters:
  !
  !    Input, real(dp) XA, YA, XB, YB, XC, YC, XD, YD, the X and Y
  !    coordinates of the four vertices.
  !
  !    Output, logical INTERSECT_PROP, the result of the test.
  !

    logical :: intersect_prop
    real(dp) :: t1
    real(dp) :: t2
    real(dp) :: t3
    real(dp) :: t4
    logical :: value
    logical :: value1
    logical :: value2
    logical :: value3
    logical :: value4
    real(dp), intent(in), value :: xa
    real(dp), intent(in), value :: xb
    real(dp), intent(in), value :: xc
    real(dp), intent(in), value :: xd
    real(dp), intent(in), value :: ya
    real(dp), intent(in), value :: yb
    real(dp), intent(in), value :: yc
    real(dp), intent(in), value :: yd

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

      value1 = ( 0.0_dp < t1 )
      value2 = ( 0.0_dp < t2 )
      value3 = ( 0.0_dp < t3 )
      value4 = ( 0.0_dp < t4 )

      value = ( l4_xor ( value1, value2 ) ) .and. ( l4_xor ( value3, value4 ) )

    end if

    intersect_prop = value
  end function intersect_prop

  pure function l4_xor ( l1, l2 ) &
        bind(C, name="l4_xor")

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

    logical, intent(in), value :: l1
    logical, intent(in), value :: l2
    logical :: l4_xor
    logical :: value1
    logical :: value2

    value1 = (         l1   .and. ( .not. l2 ) )
    value2 = ( ( .not. l1 ) .and.         l2   )

    l4_xor = ( value1 .or. value2 )
  end function l4_xor

  pure function polygon_area ( n, x, y ) &
        bind(C, name="polygon_area")

  !*****************************************************************************80
  !
  !! POLYGON_AREA returns the area of a polygon.
  !
  !  Discussion:
  !
  !    The vertices should be listed in counter-clockwise order so that
  !    the area will be positive.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    10 September 2016
  !
  !  Author:
  !
  !    John Burkardt.
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of vertices.
  !
  !    Input, real(dp) X(N), Y(N), the vertex coordinates.
  !
  !    Output, real(dp) POLYGON_AREA, the area of the polygon.
  !

    integer(ip), intent(in), value :: n

    real(dp) :: area
    integer(ip) :: i
    integer(ip) :: im1
    real(dp) :: polygon_area
    real(dp), intent(in) :: x(n)
    real(dp), intent(in) :: y(n)

    area = 0.0_dp
    im1 = n

    do i = 1, n
      area = area + x(im1) * y(i) - x(i) * y(im1)
      im1 = i
    end do

    area = 0.5_dp * area

    polygon_area = area
  end function polygon_area

  subroutine polygon_triangulate ( n, x, y, triangles ) &
        bind(C, name="polygon_triangulate")

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
  !    Thanks to Gene Dial for pointing out a mistake in the area calculation,
  !    10 September 2016.
  !
  !    Gene Dial requested an angle tolerance of about 1 millionth radian or
  !    5.7E-05 degrees, 26 June 2018.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    26 June 2018
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
  !    LC: QA448.e38_dp.
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of vertices.
  !
  !    Input, real(dp) X(N), Y(N), the coordinates of each vertex.
  !
  !    Output, integer(ip) TRIANGLES(3,N-2), the triangles of the
  !    triangulation.
  !

    integer(ip), intent(in), value :: n

    real(dp) :: angle
    real(dp), parameter :: angle_tol = 5.7e-05_dp
    real(dp) :: area
    logical :: ear(n)
    integer(ip) :: i
    integer(ip) :: i0
    integer(ip) :: i1
    integer(ip) :: i2
    integer(ip) :: i3
    integer(ip) :: i4
    integer(ip) :: next_node(n)
    integer(ip) :: node
    integer(ip) :: node_m1
    integer(ip) :: node1
    integer(ip) :: node2
    integer(ip) :: node3
    integer(ip) :: prev_node(n)
    integer(ip) :: triangle_num
    integer(ip), intent(out) :: triangles(3,n-2)
    real(dp), intent(in) :: x(n)
    real(dp), intent(in) :: y(n)
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
  !  No node can be the vertex of an angle less than 1 degree
  !  in absolute value.
  !
    node1 = n

    do node2 = 1, n

      node3 = mod ( node2, n ) + 1

      angle = angle_degree ( &
        x(node1), y(node1), &
        x(node2), y(node2), &
        x(node3), y(node3) )

      if ( abs ( angle ) <= angle_tol ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'POLYGON_TRIANGULATE - Fatal error!'
        write ( *, '(a,g14.6)' ) '  Polygon has an angle smaller than ', angle_tol
        write ( *, '(a,i4)' ) '  occurring at node ', node2
        stop 1
      end if

      node1 = node2

    end do
  !
  !  Area must be positive.
  !
    area = polygon_area ( n, x, y )

    if ( area <= 0.0_dp ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLYGON_TRIANGULATE - Fatal error!'
      write ( *, '(a)' ) '  Polygon has zero or negative area.'
      stop 1
    end if
  !
  !  PREV and NEXT point to the previous and next nodes.
  !
    i = 1
    prev_node(i) = n
    next_node(i) = i + 1

    do i = 2, n - 1
      prev_node(i) = i - 1
      next_node(i) = i + 1
    end do

    i = n
    prev_node(i) = i - 1
    next_node(i) = 1
  !
  !  EAR indicates whether the node and its immediate neighbors form an ear
  !  that can be sliced off immediately.
  !
    do i = 1, n
      ear(i) = diagonal ( prev_node(i), next_node(i), n, prev_node, next_node, &
        x, y )
    end do

    triangle_num = 0

    i2 = 1

    do while ( triangle_num < n - 3 )
  !
  !  If I2 is an ear, gather information necessary to carry out
  !  the slicing operation and subsequent "healing".
  !
      if ( ear(i2) ) then

        i3 = next_node(i2)
        i4 = next_node(i3)
        i1 = prev_node(i2)
        i0 = prev_node(i1)
  !
  !  Make vertex I2 disappear.
  !
        next_node(i1) = i3
        prev_node(i3) = i1
  !
  !  Update the earity of I1 and I3, because I2 disappeared.
  !
        ear(i1) = diagonal ( i0, i3, n, prev_node, next_node, x, y )
        ear(i3) = diagonal ( i1, i4, n, prev_node, next_node, x, y )
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
      i2 = next_node(i2)

    end do
  !
  !  The last triangle is formed from the three remaining vertices.
  !
    i3 = next_node(i2)
    i1 = prev_node(i2)

    triangle_num = triangle_num + 1
    triangles(1,triangle_num) = i3
    triangles(2,triangle_num) = i1
    triangles(3,triangle_num) = i2
  end subroutine polygon_triangulate

  pure function triangle_area ( xa, ya, xb, yb, xc, yc ) &
        bind(C, name="triangle_area")

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
  !    Input, real(dp) XA, YA, XB, YB, XC, YC, the coordinates of
  !    the vertices of the triangle, given in counterclockwise order.
  !
  !    Output, real(dp) TRIANGLE_AREA, the signed area of the triangle.
  !

    real(dp) :: triangle_area
    real(dp) :: value
    real(dp), intent(in), value :: xa
    real(dp), intent(in), value :: xb
    real(dp), intent(in), value :: xc
    real(dp), intent(in), value :: ya
    real(dp), intent(in), value :: yb
    real(dp), intent(in), value :: yc

    value = 0.5_dp * ( &
        ( xb - xa ) * ( yc - ya ) &
      - ( xc - xa ) * ( yb - ya ) )

    triangle_area = value
  end function triangle_area

end module polygon_triangulate_mod
