!> triangulation_quality — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module triangulation_quality_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: alpha_measure, arc_cosine, area_measure, bandwidth_mesh, mesh_base_one, q_measure

contains

  subroutine alpha_measure ( n, z, element_order, element_num, element_node, &
    alpha_min, alpha_ave, alpha_area ) &
        bind(C, name="alpha_measure")

  !*****************************************************************************80
  !
  !! ALPHA_MEASURE determines the triangulation quality measure ALPHA.
  !
  !  Discusion:
  !
  !    The ALPHA measure evaluates the uniformity of the shapes of the triangles
  !    defined by a triangulated pointset.
  !
  !    We compute the minimum angle among all the triangles in the triangulated
  !    dataset and divide by the maximum possible value (which, in degrees,
  !    is 60).  The best possible value is 1, and the worst 0.  A good
  !    triangulation should have an ALPHA score close to 1.
  !
  !    The code has been modified to 'allow' 6-node triangulations.
  !    However, no effort is made to actually process the midside nodes.
  !    Only information from the vertices is used.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 June 2009
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of points.
  !
  !    Input, real(dp) Z(2,N), the points.
  !
  !    Input, integer(ip) ELEMENT_ORDER, the order of the triangles.
  !
  !    Input, integer(ip) ELEMENT_NUM, the number of triangles.
  !
  !    Input, integer(ip) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM),
  !    the triangulation.
  !
  !    Output, real(dp) ALPHA_MIN, the minimum value of ALPHA over all
  !    triangles.
  !
  !    Output, real(dp) ALPHA_AVE, the value of ALPHA averaged over
  !    all triangles.
  !
  !    Output, real(dp) ALPHA_AREA, the value of ALPHA averaged over
  !    all triangles and weighted by area.
  !

    integer(ip), intent(in), value :: n                                  !! number of points
    integer(ip), intent(in), value :: element_num                        !! number of triangles
    integer(ip), intent(in), value :: element_order                      !! order of the triangles

    real(dp) :: a_angle
    integer(ip) :: a_index
    real(dp) :: a_x
    real(dp) :: a_y
    real(dp) :: ab_len
    real(dp) :: alpha
    real(dp), intent(out) :: alpha_area
    real(dp), intent(out) :: alpha_ave
    real(dp), intent(out) :: alpha_min
    real(dp) :: arc_cosine
    real(dp) :: area
    real(dp) :: area_total
    real(dp) :: b_angle
    integer(ip) :: b_index
    real(dp) :: b_x
    real(dp) :: b_y
    real(dp) :: bc_len
    real(dp) :: c_angle
    integer(ip) :: c_index
    real(dp) :: c_x
    real(dp) :: c_y
    real(dp) :: ca_len
    real(dp), parameter :: pi = 3.141592653589793_dp
    integer(ip) :: triangle
    integer(ip), intent(in) :: element_node(element_order,element_num)
    real(dp), intent(in) :: z(2,n)

    alpha_min = huge ( alpha )
    alpha_ave = 0.0_dp
    alpha_area = 0.0_dp
    area_total = 0.0_dp

    do triangle = 1, element_num

      a_index = element_node(1,triangle)
      b_index = element_node(2,triangle)
      c_index = element_node(3,triangle)

      a_x = z(1,a_index)
      a_y = z(2,a_index)
      b_x = z(1,b_index)
      b_y = z(2,b_index)
      c_x = z(1,c_index)
      c_y = z(2,c_index)

      area = 0.5_dp * abs ( a_x * ( b_y - c_y ) &
                           + b_x * ( c_y - a_y ) &
                           + c_x * ( a_y - b_y ) )

      ab_len = sqrt ( ( a_x - b_x )**2 + ( a_y - b_y )**2 )
      bc_len = sqrt ( ( b_x - c_x )**2 + ( b_y - c_y )**2 )
      ca_len = sqrt ( ( c_x - a_x )**2 + ( c_y - a_y )**2 )
  !
  !  Take care of a ridiculous special case.
  !
      if ( ab_len == 0.0_dp .and. &
           bc_len == 0.0_dp .and. &
           ca_len == 0.0_dp ) then

        a_angle = 2.0_dp * pi / 3.0_dp
        b_angle = 2.0_dp * pi / 3.0_dp
        c_angle = 2.0_dp * pi / 3.0_dp

      else

        if ( ca_len == 0.0_dp .or. ab_len == 0.0_dp ) then
          a_angle = pi
        else
          a_angle = arc_cosine ( ( ca_len**2 + ab_len**2 - bc_len**2 ) &
            / ( 2.0_dp * ca_len * ab_len ) )
        end if

        if ( ab_len == 0.0_dp .or. bc_len == 0.0_dp ) then
          b_angle = pi
        else
          b_angle = arc_cosine ( ( ab_len**2 + bc_len**2 - ca_len**2 ) &
            / ( 2.0_dp * ab_len * bc_len ) )
        end if

        if ( bc_len == 0.0_dp .or. ca_len == 0.0_dp ) then
          c_angle = pi
        else
          c_angle = arc_cosine ( ( bc_len**2 + ca_len**2 - ab_len**2 ) &
            / ( 2.0_dp * bc_len * ca_len ) )
        end if

      end if

      alpha_min = min ( alpha_min, a_angle )
      alpha_min = min ( alpha_min, b_angle )
      alpha_min = min ( alpha_min, c_angle )

      alpha_ave = alpha_ave + alpha_min

      alpha_area = alpha_area + area * alpha_min

      area_total = area_total + area

    end do

    alpha_ave = alpha_ave / real ( element_num, dp)
    alpha_area = alpha_area / area_total
  !
  !  Normalize angles from [0,pi/3] radians into qualities in [0,1].
  !
    alpha_min = alpha_min * 3.0_dp / pi
    alpha_ave = alpha_ave * 3.0_dp / pi
    alpha_area = alpha_area * 3.0_dp / pi
  end subroutine alpha_measure

  function arc_cosine ( c ) &
        bind(C, name="arc_cosine")

  !*****************************************************************************80
  !
  !! ARC_COSINE computes the arc cosine function, with argument truncation.
  !
  !  Discussion:
  !
  !    If you call your system ACOS routine with an input argument that is
  !    even slightly outside the range [-1.0, 1.0 ], you may get an unpleasant
  !    surprise (I did).
  !
  !    This routine simply truncates arguments outside the range.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    02 December 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real(dp) C, the argument.
  !
  !    Output, real(dp) ARC_COSINE, an angle whose cosine is C.
  !

    real(dp) :: arc_cosine
    real(dp), intent(in), value :: c                                     !! argument
    real(dp) :: c2

    c2 = c
    c2 = max ( c2, -1.0_dp )
    c2 = min ( c2, +1.0_dp )

    arc_cosine = acos ( c2 )
  end function arc_cosine

  subroutine area_measure ( n, z, element_order, element_num, element_node, &
    area_min, area_max, area_ratio, area_ave, area_std, area_negative, area_zero ) &
        bind(C, name="area_measure")

  !*****************************************************************************80
  !
  !! AREA_MEASURE determines the area ratio quality measure.
  !
  !  Discusion:
  !
  !    This measure computes the area of every triangle, and returns
  !    the ratio of the minimum to the maximum triangle.  A value of
  !    1 is "perfect", indicating that all triangles have the same area.
  !    A value of 0 is the worst possible result.
  !
  !    For these measurements, the absolute value of the area is considered,
  !    ignoring the triangle orientation.
  !
  !    However, the routine also returns a count of the number of triangles
  !    whose area is negative, or zero.
  !
  !    The code has been modified to 'allow' 6-node triangulations.
  !    However, no effort is made to actually process the midside nodes.
  !    Only information from the vertices is used.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    23 November 2011
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of points.
  !
  !    Input, real(dp) Z(2,N), the points.
  !
  !    Input, integer(ip) ELEMENT_ORDER, the order of the triangles.
  !
  !    Input, integer(ip) ELEMENT_NUM, the number of triangles.
  !
  !    Input, integer(ip) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM),
  !    the triangulation.
  !
  !    Output, real(dp) AREA_MIN, AREA_MAX, the minimum and maximum
  !    areas.
  !
  !    Output, real(dp) AREA_RATIO, the ratio of the minimum to the
  !    maximum area.
  !
  !    Output, real(dp) AREA_AVE, the average area.
  !
  !    Output, real(dp) AREA_STD, the standard deviation of the areas.
  !
  !    Output, integer(ip) AREA_NEGATIVE, the number of triangles with
  !    negative area.  This suggests an orientation error.
  !
  !    Output, integer(ip) AREA_ZERO, the number of triangles with zero
  !    area.
  !

    integer(ip), intent(in), value :: n                                  !! number of points
    integer(ip), intent(in), value :: element_num                        !! number of triangles
    integer(ip), intent(in), value :: element_order                      !! order of the triangles

    real(dp) :: area
    real(dp), intent(out) :: area_ave
    real(dp), intent(out) :: area_max
    real(dp), intent(out) :: area_min
    integer(ip), intent(out) :: area_negative
    real(dp), intent(out) :: area_ratio
    real(dp), intent(out) :: area_std
    integer(ip), intent(out) :: area_zero
    integer(ip) :: triangle
    integer(ip), intent(in) :: element_node(element_order,element_num)
    real(dp) :: x1
    real(dp) :: x2
    real(dp) :: x3
    real(dp) :: y1
    real(dp) :: y2
    real(dp) :: y3
    real(dp), intent(in) :: z(2,n)

    area_max = 0.0_dp
    area_min = huge ( area_min )
    area_ave = 0.0_dp

    area_negative = 0
    area_zero = 0

    do triangle = 1, element_num

      x1 = z(1,element_node(1,triangle))
      y1 = z(2,element_node(1,triangle))
      x2 = z(1,element_node(2,triangle))
      y2 = z(2,element_node(2,triangle))
      x3 = z(1,element_node(3,triangle))
      y3 = z(2,element_node(3,triangle))

      area = 0.5_dp * ( x1 * ( y2 - y3 ) &
                       + x2 * ( y3 - y1 ) &
                       + x3 * ( y1 - y2 ) )

      if ( area == 0.0_dp ) then
        area_zero = area_zero + 1
      end if

      if ( area < 0.0_dp ) then
        area_negative = area_negative + 1
      end if

      area_min = min ( area_min, abs ( area ) )
      area_max = max ( area_max, abs ( area ) )
      area_ave = area_ave + abs ( area )

    end do

    area_ave = area_ave / real ( element_num, dp)

    area_std = 0.0_dp
    do triangle = 1, element_num

      x1 = z(1,element_node(1,triangle))
      y1 = z(2,element_node(1,triangle))
      x2 = z(1,element_node(2,triangle))
      y2 = z(2,element_node(2,triangle))
      x3 = z(1,element_node(3,triangle))
      y3 = z(2,element_node(3,triangle))

      area = 0.5_dp * abs ( x1 * ( y2 - y3 ) &
                           + x2 * ( y3 - y1 ) &
                           + x3 * ( y1 - y2 ) )

      area_std = area_std + ( area - area_ave )**2
    end do
    area_std = sqrt ( area_std / real ( element_num, dp) )

    if ( 0.0_dp < area_max ) then
      area_ratio = area_min / area_max
    else
      area_ratio = 0.0_dp
    end if
  end subroutine area_measure

  subroutine bandwidth_mesh ( element_order, element_num, element_node, &
    ml, mu, m ) &
        bind(C, name="bandwidth_mesh")

  !*****************************************************************************80
  !
  !! BANDWIDTH_MESH: bandwidth of finite element mesh.
  !
  !  Discussion:
  !
  !    The quantity computed here is the "geometric" bandwidth determined
  !    by the finite element mesh alone.
  !
  !    If a single finite element variable is associated with each node
  !    of the mesh, and if the nodes and variables are numbered in the
  !    same way, then the geometric bandwidth is the same as the bandwidth
  !    of a typical finite element matrix.
  !
  !    The bandwidth M is defined in terms of the lower and upper bandwidths:
  !
  !      M = ML + 1 + MU
  !
  !    where
  !
  !      ML = maximum distance from any diagonal entry to a nonzero
  !      entry in the same row, but earlier column,
  !
  !      MU = maximum distance from any diagonal entry to a nonzero
  !      entry in the same row, but later column.
  !
  !    Because the finite element node adjacency relationship is symmetric,
  !    we are guaranteed that ML = MU.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    02 September 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) ELEMENT_ORDER, the order of the elements.
  !
  !    Input, integer(ip) ELEMENT_NUM, the number of elements.
  !
  !    Input, integer(ip) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
  !    ELEMENT_NODE(I,J) is the global index of local node I in element J.
  !
  !    Output, integer(ip) ML, MU, the lower and upper bandwidths
  !    of the matrix.
  !
  !    Output, integer(ip) M, the bandwidth of the matrix.
  !

    integer(ip), intent(in), value :: element_num                        !! number of elements
    integer(ip), intent(in), value :: element_order                      !! order of the elements

    integer(ip) :: element
    integer(ip), intent(in) :: element_node(element_order,element_num)
    integer(ip) :: global_i
    integer(ip) :: global_j
    integer(ip) :: local_i
    integer(ip) :: local_j
    integer(ip), intent(out) :: m
    integer(ip), intent(out) :: ml
    integer(ip), intent(out) :: mu

    ml = 0
    mu = 0

    do element = 1, element_num

      do local_i = 1, element_order
        global_i = element_node(local_i,element)

        do local_j = 1, element_order
          global_j = element_node(local_j,element)

          mu = max ( mu, global_j - global_i )
          ml = max ( ml, global_i - global_j )

        end do
      end do
    end do

    m = ml + 1 + mu
  end subroutine bandwidth_mesh

  subroutine mesh_base_one ( node_num, element_order, element_num, element_node ) &
        bind(C, name="mesh_base_one")

  !*****************************************************************************80
  !
  !! MESH_BASE_ONE ensures that the element definition is one-based.
  !
  !  Discussion:
  !
  !    The ELEMENT_NODE array contains nodes indices that form elements.
  !    The convention for node indexing might start at 0 or at 1.
  !    Since a FORTRAN90 program will naturally assume a 1-based indexing, it is
  !    necessary to check a given element definition and, if it is actually
  !    0-based, to convert it.
  !
  !    This function attempts to detect 9-based node indexing and correct it.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    02 October 2009
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, int NODE_NUM, the number of nodes.
  !
  !    Input, int ELEMENT_ORDER, the order of the elements.
  !
  !    Input, int ELEMENT_NUM, the number of elements.
  !
  !    Input/output, int ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), the element
  !    definitions.
  !

    integer(ip), intent(in), value :: element_num                        !! number of elements
    integer(ip), intent(in), value :: element_order                      !! order of the elements

    integer(ip) :: element
    integer(ip), intent(inout) :: element_node(element_order,element_num)
    integer(ip) :: node
    integer(ip) :: node_max
    integer(ip) :: node_min
    integer(ip), intent(in), value :: node_num                           !! number of nodes
    integer(ip) :: order

    node_min = node_num + 1
    node_max = -1

    node_min = minval ( element_node(1:element_order,1:element_num) )
    node_max = maxval ( element_node(1:element_order,1:element_num) )

    if ( node_min == 0 .and. node_max == node_num - 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )'MESH_BASE_ONE:'
      write ( *, '(a)' )'  The element indexing appears to be 0-based!'
      write ( *, '(a)' )'  This will be converted to 1-based.'
      element_node(1:element_order,1:element_num) = &
        element_node(1:element_order,1:element_num) + 1
    else if ( node_min == 1 .and. node_max == node_num  ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )'MESH_BASE_ONE:'
      write ( *, '(a)' )'  The element indexing appears to be 1-based!'
      write ( *, '(a)' )'  No conversion is necessary.'
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MESH_BASE_ONE - Warning!'
      write ( *, '(a)' ) '  The element indexing is not of a recognized type.'
      write ( *, '(a,i8)' ) '  NODE_MIN = ', node_min
      write ( *, '(a,i8)' ) '  NODE_MAX = ', node_max
      write ( *, '(a,i8)' ) '  NODE_NUM = ', node_num
    end if
  end subroutine mesh_base_one

  subroutine q_measure ( n, z, element_order, element_num, element_node, &
    q_min, q_max, q_ave, q_area ) &
        bind(C, name="q_measure")

  !*****************************************************************************80
  !
  !! Q_MEASURE determines the triangulated pointset quality measure Q.
  !
  !  Discussion:
  !
  !    The Q measure evaluates the uniformity of the shapes of the triangles
  !    defined by a triangulated pointset.
  !
  !    For a single triangle T, the value of Q(T) is defined as follows:
  !
  !      TAU_IN = radius of the inscribed circle,
  !      TAU_OUT = radius of the circumscribed circle,
  !
  !      Q(T) = 2 * TAU_IN / TAU_OUT
  !        = ( B + C - A ) * ( C + A - B ) * ( A + B - C ) / ( A * B * C )
  !
  !    where A, B and C are the lengths of the sides of the triangle T.
  !
  !    The Q measure computes the value of Q(T) for every triangle T in the
  !    triangulation, and then computes the minimum of this
  !    set of values:
  !
  !      Q_MEASURE = min ( all T in triangulation ) Q(T)
  !
  !    In an ideally regular mesh, all triangles would have the same
  !    equilateral shape, for which Q = 1.  A good mesh would have
  !    0.5 < Q.
  !
  !    Given the 2D coordinates of a set of N nodes, stored as Z(1:2,1:N),
  !    a triangulation is a list of TRIANGLE_NUM triples of node indices that form
  !    triangles.  Generally, a maximal triangulation is expected, namely,
  !    a triangulation whose image is a planar graph, but for which the
  !    addition of any new triangle would mean the graph was no longer planar.
  !    A Delaunay triangulation is a maximal triangulation which maximizes
  !    the minimum angle that occurs in any triangle.
  !
  !    The code has been modified to 'allow' 6-node triangulations.
  !    However, no effort is made to actually process the midside nodes.
  !    Only information from the vertices is used.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 June 2009
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Max Gunzburger and John Burkardt,
  !    Uniformity Measures for Point Samples in Hypercubes.
  !
  !    Per-Olof Persson and Gilbert Strang,
  !    A Simple Mesh Generator in MATLAB,
  !    SIAM Review,
  !    Volume 46, Number 2, pages 329-345, June 2004.
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of points.
  !
  !    Input, real(dp) Z(2,N), the points.
  !
  !    Input, integer(ip) ELEMENT_ORDER, the order of the triangles.
  !
  !    Input, integer(ip) ELEMENT_NUM, the number of triangles.
  !
  !    Input, integer(ip) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM),
  !    the triangulation.
  !
  !    Output, real(dp) Q_MIN, Q_MAX, the minimum and maximum values
  !    of Q over all triangles.
  !
  !    Output, real(dp) Q_AVE, the average value of Q.
  !
  !    Output, real(dp) Q_AREA, the average value of Q, weighted by
  !    the area of each triangle.
  !

    integer(ip), intent(in), value :: n                                  !! number of points
    integer(ip), intent(in), value :: element_num                        !! number of triangles
    integer(ip), intent(in), value :: element_order                      !! order of the triangles

    integer(ip) :: a_index
    real(dp) :: ab_length
    real(dp) :: area
    real(dp) :: area_total
    integer(ip) :: b_index
    real(dp) :: bc_length
    integer(ip) :: c_index
    real(dp) :: ca_length
    real(dp) :: q
    real(dp), intent(out) :: q_area
    real(dp), intent(out) :: q_ave
    real(dp), intent(out) :: q_max
    real(dp), intent(out) :: q_min
    integer(ip) :: triangle
    integer(ip), intent(in) :: element_node(element_order,element_num)
    real(dp) :: x1
    real(dp) :: x2
    real(dp) :: x3
    real(dp) :: y1
    real(dp) :: y2
    real(dp) :: y3
    real(dp), intent(in) :: z(2,n)

    q_min =   huge ( q_min )
    q_max = - huge ( q_max )
    q_ave = 0.0_dp
    q_area = 0.0_dp
    area_total = 0.0_dp

    do triangle = 1, element_num

      a_index = element_node(1,triangle)
      b_index = element_node(2,triangle)
      c_index = element_node(3,triangle)

      ab_length = sqrt ( &
          ( z(1,a_index) - z(1,b_index) )**2 &
        + ( z(2,a_index) - z(2,b_index) )**2 )

      bc_length = sqrt ( &
          ( z(1,b_index) - z(1,c_index) )**2 &
        + ( z(2,b_index) - z(2,c_index) )**2 )

      ca_length = sqrt ( &
          ( z(1,c_index) - z(1,a_index) )**2 &
        + ( z(2,c_index) - z(2,a_index) )**2 )

      q = ( bc_length + ca_length - ab_length ) &
        * ( ca_length + ab_length - bc_length ) &
        * ( ab_length + bc_length - ca_length ) &
        / ( ab_length * bc_length * ca_length )

      x1 = z(1,element_node(1,triangle))
      y1 = z(2,element_node(1,triangle))
      x2 = z(1,element_node(2,triangle))
      y2 = z(2,element_node(2,triangle))
      x3 = z(1,element_node(3,triangle))
      y3 = z(2,element_node(3,triangle))

      area = 0.5_dp * abs ( x1 * ( y2 - y3 ) &
                           + x2 * ( y3 - y1 ) &
                           + x3 * ( y1 - y2 ) )

      q_min = min ( q_min, q )
      q_max = max ( q_max, q )
      q_ave = q_ave + q
      q_area = q_area + q * area

      area_total = area_total + area

    end do

    q_ave = q_ave / real ( element_num, dp)

    if ( 0.0_dp < area_total ) then
      q_area = q_area / area_total
    else
      q_area = 0.0_dp
    end if
  end subroutine q_measure

end module triangulation_quality_mod
