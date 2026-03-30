!> geompack — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

function angle_rad_2d ( p1, p2, p3 )

!*****************************************************************************80
!
!! ANGLE_RAD_2D returns the angle swept out between two rays in 2D.
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

function diaedg ( x0, y0, x1, y1, x2, y2, x3, y3 )

!*****************************************************************************80
!
!! DIAEDG chooses a diagonal edge.
!
!  Discussion:
!
!    The routine determines whether 0--2 or 1--3 is the diagonal edge
!    that should be chosen, based on the circumcircle criterion, where
!    (X0,Y0), (X1,Y1), (X2,Y2), (X3,Y3) are the vertices of a simple
!    quadrilateral in counterclockwise order.
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, double precision X0, Y0, X1, Y1, X2, Y2, X3, Y3, the
!    coordinates of the vertices of a quadrilateral, given in
!    counter clockwise order.
!
!    Output, integer DIAEDG, chooses a diagonal:
!    +1, if diagonal edge 02 is chosen;
!    -1, if diagonal edge 13 is chosen;
!     0, if the four vertices are cocircular.
!
  implicit none

  double precision ca
  double precision cb
  integer diaedg
  double precision dx10
  double precision dx12
  double precision dx30
  double precision dx32
  double precision dy10
  double precision dy12
  double precision dy30
  double precision dy32
  double precision s
  double precision tol
  double precision tola
  double precision tolb
  double precision x0
  double precision x1
  double precision x2
  double precision x3
  double precision y0
  double precision y1
  double precision y2
  double precision y3

  tol = 100.0D+00 * epsilon ( tol )

  dx10 = x1 - x0
  dy10 = y1 - y0
  dx12 = x1 - x2
  dy12 = y1 - y2
  dx30 = x3 - x0
  dy30 = y3 - y0
  dx32 = x3 - x2
  dy32 = y3 - y2

  tola = tol * max ( abs ( dx10 ), abs ( dy10 ), abs ( dx30 ), abs ( dy30 ) )
  tolb = tol * max ( abs ( dx12 ), abs ( dy12 ), abs ( dx32 ), abs ( dy32 ) )

  ca = dx10 * dx30 + dy10 * dy30
  cb = dx12 * dx32 + dy12 * dy32

  if ( tola < ca .and. tolb < cb ) then

    diaedg = -1

  else if ( ca < -tola .and. cb < -tolb ) then

    diaedg = 1

  else

    tola = max ( tola, tolb )
    s = ( dx10 * dy30 - dx30 * dy10 ) * cb + ( dx32 * dy12 - dx12 * dy32 ) * ca

    if ( tola < s ) then
      diaedg = -1
    else if ( s < -tola ) then
      diaedg = 1
    else
      diaedg = 0
    end if

  end if
end

function lrline ( xu, yu, xv1, yv1, xv2, yv2, dv )

!*****************************************************************************80
!
!! LRLINE determines if a point is left of, right or, or on a directed line.
!
!  Discussion:
!
!    The directed line is parallel to, and at a signed distance DV from
!    a directed base line from (XV1,YV1) to (XV2,YV2).
!
!  Modified:
!
!    14 July 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, double precision XU, YU, the coordinates of the point whose
!    position relative to the directed line is to be determined.
!
!    Input, double precision XV1, YV1, XV2, YV2, the coordinates of two points
!    that determine the directed base line.
!
!    Input, double precision DV, the signed distance of the directed line
!    from the directed base line through the points (XV1,YV1) and (XV2,YV2).
!    DV is positive for a line to the left of the base line.
!
!    Output, integer LRLINE, the result:
!    +1, the point is to the right of the directed line;
!     0, the point is on the directed line;
!    -1, the point is to the left of the directed line.
!
  implicit none

  double precision dv
  double precision dx
  double precision dxu
  double precision dy
  double precision dyu
  integer lrline
  double precision t
  double precision tol
  double precision tolabs
  double precision xu
  double precision xv1
  double precision xv2
  double precision yu
  double precision yv1
  double precision yv2

  tol = 100.0D+00 * epsilon ( tol )

  dx = xv2 - xv1
  dy = yv2 - yv1
  dxu = xu - xv1
  dyu = yu - yv1

  tolabs = tol * max ( abs ( dx ), abs ( dy ), abs ( dxu ), &
    abs ( dyu ), abs ( dv ) )

  t = dy * dxu - dx * dyu + dv * sqrt ( dx * dx + dy * dy )

  if ( tolabs < t ) then
    lrline = 1
  else if ( -tolabs <= t ) then
    lrline = 0
  else
    lrline = -1
  end if
end

subroutine points_delaunay_naive_2d ( node_num, node_xy, maxtri, &
  triangle_num, triangle_node )

!*****************************************************************************80
!
!! POINTS_DELAUNAY_NAIVE_2D is a naive Delaunay triangulation scheme.
!
!  Discussion:
!
!    This routine is only suitable as a demonstration code for small
!    problems.  Its running time is of order NODE_NUM**4.  Much faster
!    algorithms are available.
!
!    Given a set of nodes in the plane, a triangulation is set of
!    triples of distinct nodes, forming triangles, so that every
!    point within the convex hull of the set of nodes is either
!    one of the nodes, or lies on an edge of one or more triangles,
!    or lies within exactly one triangle.
!
!    A Delaunay triangulation is a triangulation with additional
!    properties.
!
!    NODE_NUM must be at least 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 November 2000
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
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, double precision NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, integer MAXTRI, the maximum number of triangles.
!
!    Output, integer TRIANGLE_NUM, the number of triangles in
!    the triangulation.
!
!    Output, integer TRIANGLE_NODE(3,MAXTRI), the indices of the
!    triangle nodes.
!
  implicit none

  integer , parameter :: dim_num = 2
  integer maxtri
  integer node_num

  logical flag
  integer i
  integer j
  integer k
  integer m
  double precision node_xy(dim_num,node_num)
  integer triangle_node(3,maxtri)
  integer triangle_num
  double precision xn
  double precision yn
  double precision z(node_num)
  double precision zn

  triangle_num = 0

  if ( node_num < 3 ) then
  end if
!
!  Compute Z = X*X + Y*Y.
!
  z(1:node_num) = node_xy(1,1:node_num)**2 + node_xy(2,1:node_num)**2
!
!  For each triple (I,J,K):
!
  do i = 1, node_num - 2
    do j = i+1, node_num
      do k = i+1, node_num

        if ( j /= k ) then

          xn = ( node_xy(2,j) - node_xy(2,i) ) * ( z(k) - z(i) ) &
             - ( node_xy(2,k) - node_xy(2,i) ) * ( z(j) - z(i) )

          yn = ( node_xy(1,k) - node_xy(1,i) ) * ( z(j) - z(i) ) &
             - ( node_xy(1,j) - node_xy(1,i) ) * ( z(k) - z(i) )

          zn = ( node_xy(1,j) - node_xy(1,i) ) &
             * ( node_xy(2,k) - node_xy(2,i) ) &
             - ( node_xy(1,k) - node_xy(1,i) ) &
             * ( node_xy(2,j) - node_xy(2,i) )

          flag = ( zn < 0.0D+00 )

          if ( flag ) then
            do m = 1, node_num
              flag = flag .and. &
                ( ( node_xy(1,m) - node_xy(1,i) ) * xn &
                + ( node_xy(2,m) - node_xy(2,i) ) * yn &
                + ( z(m)   - z(i) )   * zn <= 0.0D+00 )
            end do
          end if

          if ( flag ) then
            if ( triangle_num < maxtri ) then
              triangle_num = triangle_num + 1
              triangle_node(1:3,triangle_num) = (/ i, j, k /)
            end if
          end if

        end if

      end do
    end do
  end do
end

subroutine points_hull_2d ( node_num, node_xy, hull_num, hull )

!*****************************************************************************80
!
!! POINTS_HULL_2D computes the convex hull of 2D points.
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
!    Output, integer HULL_NUM, the number of nodes that lie on
!    the convex hull.
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
       ( node_xy(1,i) == node_xy(1,q) .and. node_xy(2,i) < node_xy(2,q) ) ) then
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
      write ( *, '(a)' ) 'POINTS_HULL_2D - Fatal error!'
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
end

subroutine quad_convex_random ( seed, xy )

!*****************************************************************************80
!
!! QUAD_CONVEX_RANDOM returns a random convex quadrilateral.
!
!  Description:
!
!    The quadrilateral is constrained in that the vertices must all lie
!    with the unit square.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 June 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer SEED, a seed for the random number
!    generator.
!
!    Output, double precision XY(2,NODE_NUM), the coordinates of the
!    nodes of the quadrilateral, given in counterclockwise order.
!
  implicit none

  integer , parameter :: node_num = 4

  integer hull(node_num)
  integer hull_num
  integer j
  integer seed
  double precision xy(2,node_num)
  double precision xy_random(2,node_num)

  do
!
!  Generate 4 random points.
!
    call r8mat_uniform_01 ( 2, node_num, seed, xy_random )
!
!  Determine the convex hull.
!
    call points_hull_2d ( node_num, xy_random, hull_num, hull )
!
!  If HULL_NUM < NODE_NUM, then our convex hull is a triangle.
!  Try again.
!
    if ( hull_num == node_num ) then
      exit
    end if

  end do
!
!  Make an ordered copy of the random points.
!
  do j = 1, node_num
    xy(1:2,j) = xy_random(1:2,hull(j))
  end do
end

subroutine r8tris2 ( node_num, node_xy, triangle_num, triangle_node, &
  triangle_neighbor )

!*****************************************************************************80
!
!! R8TRIS2 constructs a Delaunay triangulation of 2D vertices.
!
!  Discussion:
!
!    The routine constructs the Delaunay triangulation of a set of 2D vertices
!    using an incremental approach and diagonal edge swaps.  Vertices are
!    first sorted in lexicographically increasing (X,Y) order, and
!    then are inserted one at a time from outside the convex hull.
!
!  Modified:
!
!    25 August 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of vertices.
!
!    Input/output, double precision NODE_XY(2,NODE_NUM), the coordinates
!    of the vertices.  On output, the vertices have been sorted into
!    dictionary order.
!
!    Output, integer TRIANGLE_NUM, the number of triangles in
!    the triangulation;  TRIANGLE_NUM is equal to 2*NODE_NUM - NB - 2, where
!    NB is the number of boundary vertices.
!
!    Output, integer TRIANGLE_NODE(3,TRIANGLE_NUM), the nodes
!    that make up each triangle.  The elements are indices of P.  The vertices
!    of the triangles are in counter clockwise order.
!
!    Output, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM),
!    the triangle neighbor list.  Positive elements are indices of TIL;
!    negative elements are used for links of a counter clockwise linked list
!    of boundary edges; LINK = -(3*I + J-1) where I, J = triangle, edge index;
!    TRIANGLE_NEIGHBOR(J,I) refers to the neighbor along edge from vertex J
!    to J+1 (mod 3).
!
  implicit none

  integer node_num

  double precision cmax
  integer e
  integer i
  integer ierr
  integer indx(node_num)
  integer j
  integer k
  integer l
  integer ledg
  integer lr
  integer lrline
  integer ltri
  integer m
  integer m1
  integer m2
  integer n
  double precision node_xy(2,node_num)
  integer redg
  integer rtri
  integer stack(node_num)
  integer t
  double precision tol
  integer top
  integer triangle_neighbor(3,node_num*2)
  integer triangle_num
  integer triangle_node(3,node_num*2)

  tol = 100.0D+00 * epsilon ( tol )

  ierr = 0
!
!  Sort the vertices by increasing (x,y).
!
  call r82vec_sort_heap_index_a ( node_num, node_xy, indx )

  call r82vec_permute ( node_num, node_xy, indx )
!
!  Make sure that the data points are "reasonably" distinct.
!
  m1 = 1

  do i = 2, node_num

    m = m1
    m1 = i

    k = 0

    do j = 1, 2

      cmax = max ( abs ( node_xy(j,m) ), abs ( node_xy(j,m1) ) )

      if ( tol * ( cmax + 1.0D+00 ) &
           < abs ( node_xy(j,m) - node_xy(j,m1) ) ) then
        k = j
        exit
      end if

    end do

    if ( k == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8TRIS2 - Fatal error!'
      write ( *, '(a,i8)' ) '  Fails for point number I = ', i
      write ( *, '(a,i8)' ) '  M = ', m
      write ( *, '(a,i8)' ) '  M1 = ', m1
      write ( *, '(a,2g14.6)' ) '  X,Y(M)  = ', node_xy(1:2,m)
      write ( *, '(a,2g14.6)' ) '  X,Y(M1) = ', node_xy(1:2,m1)
      ierr = 224
    end if

  end do
!
!  Starting from points M1 and M2, search for a third point M that
!  makes a "healthy" triangle (M1,M2,M)
!
  m1 = 1
  m2 = 2
  j = 3

  do

    if ( node_num < j ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8TRIS2 - Fatal error!'
      ierr = 225
    end if

    m = j

    lr = lrline ( node_xy(1,m), node_xy(2,m), node_xy(1,m1), &
      node_xy(2,m1), node_xy(1,m2), node_xy(2,m2), 0.0D+00 )

    if ( lr /= 0 ) then
      exit
    end if

    j = j + 1

  end do
!
!  Set up the triangle information for (M1,M2,M), and for any other
!  triangles you created because points were collinear with M1, M2.
!
  triangle_num = j - 2

  if ( lr == -1 ) then

    triangle_node(1,1) = m1
    triangle_node(2,1) = m2
    triangle_node(3,1) = m
    triangle_neighbor(3,1) = -3

    do i = 2, triangle_num

      m1 = m2
      m2 = i+1
      triangle_node(1,i) = m1
      triangle_node(2,i) = m2
      triangle_node(3,i) = m
      triangle_neighbor(1,i-1) = -3 * i
      triangle_neighbor(2,i-1) = i
      triangle_neighbor(3,i) = i - 1

    end do

    triangle_neighbor(1,triangle_num) = -3 * triangle_num - 1
    triangle_neighbor(2,triangle_num) = -5
    ledg = 2
    ltri = triangle_num

  else

    triangle_node(1,1) = m2
    triangle_node(2,1) = m1
    triangle_node(3,1) = m
    triangle_neighbor(1,1) = -4

    do i = 2, triangle_num
      m1 = m2
      m2 = i+1
      triangle_node(1,i) = m2
      triangle_node(2,i) = m1
      triangle_node(3,i) = m
      triangle_neighbor(3,i-1) = i
      triangle_neighbor(1,i) = -3 * i - 3
      triangle_neighbor(2,i) = i - 1
    end do

    triangle_neighbor(3,triangle_num) = -3 * triangle_num
    triangle_neighbor(2,1) = -3 * triangle_num - 2
    ledg = 2
    ltri = 1

  end if
!
!  Insert the vertices one at a time from outside the convex hull,
!  determine visible boundary edges, and apply diagonal edge swaps until
!  Delaunay triangulation of vertices (so far) is obtained.
!
  top = 0

  do i = j+1, node_num

    m = i
    m1 = triangle_node(ledg,ltri)

    if ( ledg <= 2 ) then
      m2 = triangle_node(ledg+1,ltri)
    else
      m2 = triangle_node(1,ltri)
    end if

    lr = lrline ( node_xy(1,m), node_xy(2,m), node_xy(1,m1), &
      node_xy(2,m1), node_xy(1,m2), node_xy(2,m2), 0.0D+00 )

    if ( 0 < lr ) then
      rtri = ltri
      redg = ledg
      ltri = 0
    else
      l = -triangle_neighbor(ledg,ltri)
      rtri = l / 3
      redg = mod(l,3) + 1
    end if

    call vbedg ( node_xy(1,m), node_xy(2,m), node_num, node_xy, triangle_num, &
      triangle_node, triangle_neighbor, ltri, ledg, rtri, redg )

    n = triangle_num + 1
    l = -triangle_neighbor(ledg,ltri)

    do

      t = l / 3
      e = mod ( l, 3 ) + 1
      l = -triangle_neighbor(e,t)
      m2 = triangle_node(e,t)

      if ( e <= 2 ) then
        m1 = triangle_node(e+1,t)
      else
        m1 = triangle_node(1,t)
      end if

      triangle_num = triangle_num + 1
      triangle_neighbor(e,t) = triangle_num
      triangle_node(1,triangle_num) = m1
      triangle_node(2,triangle_num) = m2
      triangle_node(3,triangle_num) = m
      triangle_neighbor(1,triangle_num) = t
      triangle_neighbor(2,triangle_num) = triangle_num - 1
      triangle_neighbor(3,triangle_num) = triangle_num + 1
      top = top + 1

      if ( node_num < top ) then
        ierr = 8
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8TRIS2 - Fatal error!'
        write ( *, '(a)' ) '  Stack overflow.'
      end if

      stack(top) = triangle_num

      if ( t == rtri .and. e == redg ) then
        exit
      end if

    end do

    triangle_neighbor(ledg,ltri) = -3 * n - 1
    triangle_neighbor(2,n) = -3 * triangle_num - 2
    triangle_neighbor(3,triangle_num) = -l
    ltri = n
    ledg = 2

    call swapec ( m, top, ltri, ledg, node_num, node_xy, triangle_num, &
      triangle_node, triangle_neighbor, stack, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8TRIS2 - Fatal error!'
      write ( *, '(a)' ) '  Error return from SWAPEC.'
    end if

  end do
!
!  Now account for the sorting that we did.
!
  do i = 1, 3
    do j = 1, triangle_num
      triangle_node(i,j) = indx ( triangle_node(i,j) )
    end do
  end do

  call perm_inverse ( node_num, indx )

  call r82vec_permute ( node_num, node_xy, indx )
end

subroutine swapec ( i, top, btri, bedg, node_num, node_xy, triangle_num, &
  triangle_node, triangle_neighbor, stack, ierr )

!*****************************************************************************80
!
!! SWAPEC swaps diagonal edges until all triangles are Delaunay.
!
!  Discussion:
!
!    The routine swaps diagonal edges in a 2D triangulation, based on
!    the empty circumcircle criterion, until all triangles are Delaunay,
!    given that I is the index of the new vertex added to the triangulation.
!
!  Modified:
!
!    14 July 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer I, the index of the new vertex.
!
!    Input/output, integer TOP, the index of the top of the stack.
!    On output, TOP is zero.
!
!    Input/output, integer BTRI, BEDG; on input, if positive, are
!    the triangle and edge indices of a boundary edge whose updated indices
!    must be recorded.  On output, these may be updated because of swaps.
!
!    Input, integer NODE_NUM, the number of points.
!
!    Input, double precision NODE_XY(2,NODE_NUM), the coordinates of
!    the points.
!
!    Input, integer TRIANGLE_NUM, the number of triangles.
!
!    Input/output, integer TRIANGLE_NODE(3,TRIANGLE_NUM), the 
!    triangle incidence list.  May be updated on output because of swaps.
!
!    Input/output, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM),
!    the triangle neighbor list; negative values are used for links of the
!    counter-clockwise linked list of boundary edges;  May be updated on output
!    because of swaps.
!    LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Workspace, integer STACK(MAXST); on input, entries 1 through 
!    TOP contain the indices of initial triangles (involving vertex I)
!    put in stack; the edges opposite I should be in interior;  entries
!    TOP+1 through MAXST are used as a stack.
!
!    Output, integer IERR is set to 8 for abnormal return.
!
  implicit none

  integer node_num
  integer triangle_num

  integer a
  integer b
  integer bedg
  integer btri
  integer c
  integer diaedg
  integer e
  integer ee
  integer em1
  integer ep1
  integer f
  integer fm1
  integer fp1
  integer i
  integer ierr
  integer i4_wrap
  integer l
  double precision node_xy(2,node_num)
  integer r
  integer s
  integer stack(node_num)
  integer swap
  integer t
  integer top
  integer triangle_neighbor(3,triangle_num)
  integer triangle_node(3,triangle_num)
  integer tt
  integer u
  double precision x
  double precision y
!
!  Determine whether triangles in stack are Delaunay, and swap
!  diagonal edge of convex quadrilateral if not.
!
  x = node_xy(1,i)
  y = node_xy(2,i)

  do

    if ( top <= 0 ) then
      exit
    end if

    t = stack(top)
    top = top - 1

    if ( triangle_node(1,t) == i ) then
      e = 2
      b = triangle_node(3,t)
    else if ( triangle_node(2,t) == i ) then
      e = 3
      b = triangle_node(1,t)
    else
      e = 1
      b = triangle_node(2,t)
    end if

    a = triangle_node(e,t)
    u = triangle_neighbor(e,t)

    if ( triangle_neighbor(1,u) == t ) then
      f = 1
      c = triangle_node(3,u)
    else if ( triangle_neighbor(2,u) == t ) then
      f = 2
      c = triangle_node(1,u)
    else
      f = 3
      c = triangle_node(2,u)
    end if

    swap = diaedg ( x, y, node_xy(1,a), node_xy(2,a), node_xy(1,c), &
      node_xy(2,c), node_xy(1,b), node_xy(2,b) )

    if ( swap == 1 ) then

      em1 = i4_wrap ( e - 1, 1, 3 )
      ep1 = i4_wrap ( e + 1, 1, 3 )
      fm1 = i4_wrap ( f - 1, 1, 3 )
      fp1 = i4_wrap ( f + 1, 1, 3 )

      triangle_node(ep1,t) = c
      triangle_node(fp1,u) = i
      r = triangle_neighbor(ep1,t)
      s = triangle_neighbor(fp1,u)
      triangle_neighbor(ep1,t) = u
      triangle_neighbor(fp1,u) = t
      triangle_neighbor(e,t) = s
      triangle_neighbor(f,u) = r

      if ( 0 < triangle_neighbor(fm1,u) ) then
        top = top + 1
        stack(top) = u
      end if

      if ( 0 < s ) then

        if ( triangle_neighbor(1,s) == u ) then
          triangle_neighbor(1,s) = t
        else if ( triangle_neighbor(2,s) == u ) then
          triangle_neighbor(2,s) = t
        else
          triangle_neighbor(3,s) = t
        end if

        top = top + 1

        if ( node_num < top ) then
          ierr = 8
        end if

        stack(top) = t

      else

        if ( u == btri .and. fp1 == bedg ) then
          btri = t
          bedg = e
        end if

        l = - ( 3 * t + e - 1 )
        tt = t
        ee = em1

        do while ( 0 < triangle_neighbor(ee,tt) )

          tt = triangle_neighbor(ee,tt)

          if ( triangle_node(1,tt) == a ) then
            ee = 3
          else if ( triangle_node(2,tt) == a ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        triangle_neighbor(ee,tt) = l

      end if

      if ( 0 < r ) then

        if ( triangle_neighbor(1,r) == t ) then
          triangle_neighbor(1,r) = u
        else if ( triangle_neighbor(2,r) == t ) then
          triangle_neighbor(2,r) = u
        else
          triangle_neighbor(3,r) = u
        end if

      else

        if ( t == btri .and. ep1 == bedg ) then
          btri = u
          bedg = f
        end if

        l = - ( 3 * u + f - 1 )
        tt = u
        ee = fm1

        do while ( 0 < triangle_neighbor(ee,tt) )

          tt = triangle_neighbor(ee,tt)

          if ( triangle_node(1,tt) == b ) then
            ee = 3
          else if ( triangle_node(2,tt) == b ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        triangle_neighbor(ee,tt) = l

      end if

    end if

  end do
end

subroutine triangle_circumcenter_2d ( t, center )

!*****************************************************************************80
!
!! TRIANGLE_CIRCUMCENTER_2D computes the circumcenter of a triangle in 2D.
!
!  Discussion:
!
!    The circumcenter of a triangle is the center of the circumcircle, the
!    circle that passes through the three vertices of the triangle.
!
!    The circumcircle contains the triangle, but it is not necessarily the
!    smallest triangle to do so.
!
!    If all angles of the triangle are no greater than 90 degrees, then
!    the center of the circumscribed circle will lie inside the triangle.
!    Otherwise, the center will lie outside the circle.
!
!    The circumcenter is the intersection of the perpendicular bisectors
!    of the sides of the triangle.
!
!    In geometry, the circumcenter of a triangle is often symbolized by "O".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision T(2,3), the triangle vertices.
!
!    Output, double precision CENTER(2), the circumcenter of the triangle.
!
  implicit none

  integer , parameter :: dim_num = 2

  double precision asq
  double precision bot
  double precision center(dim_num)
  double precision csq
  double precision t(dim_num,3)
  double precision top(dim_num)

  asq = ( t(1,2) - t(1,1) )**2 + ( t(2,2) - t(2,1) )**2
  csq = ( t(1,3) - t(1,1) )**2 + ( t(2,3) - t(2,1) )**2

  top(1) =  ( t(2,2) - t(2,1) ) * csq - ( t(2,3) - t(2,1) ) * asq
  top(2) =  ( t(1,2) - t(1,1) ) * csq - ( t(1,3) - t(1,1) ) * asq

  bot  =  ( t(2,2) - t(2,1) ) * ( t(1,3) - t(1,1) ) &
        - ( t(2,3) - t(2,1) ) * ( t(1,2) - t(1,1) )

  center(1:2) = t(1:2,1) + 0.5D+00 * top(1:2) / bot
end

subroutine vbedg ( x, y, node_num, node_xy, triangle_num, triangle_node, &
  triangle_neighbor, ltri, ledg, rtri, redg )

!*****************************************************************************80
!
!! VBEDG determines which boundary edges are visible to a point.
!
!  Discussion:
!
!    The point (X,Y) is assumed to be outside the convex hull of the
!    region covered by the 2D triangulation.
!
!  Modified:
!
!    25 August 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, double precision X, Y, the coordinates of a point outside the
!    convex hull of the current triangulation.
!
!    Input, integer NODE_NUM, the number of points.
!
!    Input, double precision NODE_XY(2,NODE_NUM), the coordinates of the
!    vertices.
!
!    Input, integer TRIANGLE_NUM, the number of triangles.
!
!    Input, integer TRIANGLE_NODE(3,TRIANGLE_NUM), the
!    triangle incidence list.
!
!    Input, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the
!    triangle neighbor list; negative values are used for links of a
!    counterclockwise linked list of boundary edges;
!      LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Input/output, integer LTRI, LEDG.  If LTRI /= 0 then these
!    values are assumed to be already computed and are not changed, else they
!    are updated.  On output, LTRI is the index of boundary triangle to the
!    left of the leftmost boundary triangle visible from (X,Y), and LEDG is
!    the boundary edge of triangle LTRI to the left of the leftmost boundary
!    edge visible from (X,Y).  1 <= LEDG <= 3.
!
!    Input/output, integer RTRI.  On input, the index of the
!    boundary triangle to begin the search at.  On output, the index of the
!    rightmost boundary triangle visible from (X,Y).
!
!    Input/output, integer REDG, the edge of triangle RTRI that
!    is visible from (X,Y).  1 <= REDG <= 3.
!
  implicit none

  integer , parameter :: dim_num = 2
  integer node_num
  integer triangle_num

  integer a
  integer b
  integer e
  integer i4_wrap
  integer l
  logical ldone
  integer ledg
  integer lr
  integer lrline
  integer ltri
  double precision node_xy(2,node_num)
  integer redg
  integer rtri
  integer t
  integer triangle_neighbor(3,triangle_num)
  integer triangle_node(3,triangle_num)
  double precision x
  double precision y
!
!  Find the rightmost visible boundary edge using links, then possibly
!  leftmost visible boundary edge using triangle neighbor information.
!
  if ( ltri == 0 ) then
    ldone = .false.
    ltri = rtri
    ledg = redg
  else
    ldone = .true.
  end if

  do

    l = -triangle_neighbor(redg,rtri)
    t = l / 3
    e = mod ( l, 3 ) + 1
    a = triangle_node(e,t)

    if ( e <= 2 ) then
      b = triangle_node(e+1,t)
    else
      b = triangle_node(1,t)
    end if

    lr = lrline ( x, y, node_xy(1,a), node_xy(2,a), node_xy(1,b), &
      node_xy(2,b), 0.0D+00 )

    if ( lr <= 0 ) then
      exit
    end if

    rtri = t
    redg = e

  end do

  if ( ldone ) then
  end if

  t = ltri
  e = ledg

  do

    b = triangle_node(e,t)
    e = i4_wrap ( e-1, 1, 3 )

    do while ( 0 < triangle_neighbor(e,t) )

      t = triangle_neighbor(e,t)

      if ( triangle_node(1,t) == b ) then
        e = 3
      else if ( triangle_node(2,t) == b ) then
        e = 1
      else
        e = 2
      end if

    end do

    a = triangle_node(e,t)

    lr = lrline ( x, y, node_xy(1,a), node_xy(2,a), node_xy(1,b), &
      node_xy(2,b), 0.0D+00 )

    if ( lr <= 0 ) then
      exit
    end if

  end do

  ltri = t
  ledg = e
end
