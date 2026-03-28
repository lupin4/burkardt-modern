!> triangulation_quality — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module triangulation_quality_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  implicit none
  private

  public :: alpha_measure, arc_cosine, area_measure, bandwidth_mesh, mesh_base_one, q_measure

contains

  subroutine alpha_measure ( n, z, element_order, element_num, element_node, &
    alpha_min, alpha_ave, alpha_area )

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
  !    Input, integer(int32) N, the number of points.
  !
  !    Input, real(real64) Z(2,N), the points.
  !
  !    Input, integer(int32) TRIANGLE_ORDER, the order of the triangles.
  !
  !    Input, integer(int32) TRIANGLE_NUM, the number of triangles.
  !
  !    Input, integer(int32) TRIANGLE_NODE(TRIANGLE_ORDER,TRIANGLE_NUM), 
  !    the triangulation.
  !
  !    Output, real(real64) ALPHA_MIN, the minimum value of ALPHA over all
  !    triangles.
  !
  !    Output, real(real64) ALPHA_AVE, the value of ALPHA averaged over
  !    all triangles.
  !
  !    Output, real(real64) ALPHA_AREA, the value of ALPHA averaged over
  !    all triangles and weighted by area.
  !

    integer(int32) n
    integer(int32) element_num
    integer(int32) element_order

    real(real64) a_angle
    integer(int32) a_index
    real(real64) a_x
    real(real64) a_y
    real(real64) ab_len
    real(real64) alpha
    real(real64) alpha_area
    real(real64) alpha_ave
    real(real64) alpha_min
    real(real64) arc_cosine
    real(real64) area
    real(real64) area_total
    real(real64) b_angle
    integer(int32) b_index
    real(real64) b_x
    real(real64) b_y
    real(real64) bc_len
    real(real64) c_angle
    integer(int32) c_index
    real(real64) c_x
    real(real64) c_y
    real(real64) ca_len
    real(real64), parameter :: pi = 3.141592653589793e+00_real64
    integer(int32) triangle
    integer(int32) element_node(element_order,element_num)
    real(real64) z(2,n)

    alpha_min = huge ( alpha )
    alpha_ave = 0.0e+00_real64
    alpha_area = 0.0e+00_real64
    area_total = 0.0e+00_real64

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

      area = 0.5e+00_real64 * abs ( a_x * ( b_y - c_y ) &
                           + b_x * ( c_y - a_y ) &
                           + c_x * ( a_y - b_y ) )

      ab_len = sqrt ( ( a_x - b_x )**2 + ( a_y - b_y )**2 )
      bc_len = sqrt ( ( b_x - c_x )**2 + ( b_y - c_y )**2 )
      ca_len = sqrt ( ( c_x - a_x )**2 + ( c_y - a_y )**2 )
  !
  !  Take care of a ridiculous special case.
  !
      if ( ab_len == 0.0e+00_real64 .and. &
           bc_len == 0.0e+00_real64 .and. &
           ca_len == 0.0e+00_real64 ) then

        a_angle = 2.0e+00_real64 * pi / 3.0e+00_real64
        b_angle = 2.0e+00_real64 * pi / 3.0e+00_real64
        c_angle = 2.0e+00_real64 * pi / 3.0e+00_real64

      else

        if ( ca_len == 0.0e+00_real64 .or. ab_len == 0.0e+00_real64 ) then
          a_angle = pi
        else
          a_angle = arc_cosine ( ( ca_len**2 + ab_len**2 - bc_len**2 ) &
            / ( 2.0e+00_real64 * ca_len * ab_len ) )
        end if

        if ( ab_len == 0.0e+00_real64 .or. bc_len == 0.0e+00_real64 ) then
          b_angle = pi
        else
          b_angle = arc_cosine ( ( ab_len**2 + bc_len**2 - ca_len**2 ) &
            / ( 2.0e+00_real64 * ab_len * bc_len ) )
        end if

        if ( bc_len == 0.0e+00_real64 .or. ca_len == 0.0e+00_real64 ) then
          c_angle = pi
        else
          c_angle = arc_cosine ( ( bc_len**2 + ca_len**2 - ab_len**2 ) &
            / ( 2.0e+00_real64 * bc_len * ca_len ) )
        end if

      end if

      alpha_min = min ( alpha_min, a_angle )
      alpha_min = min ( alpha_min, b_angle )
      alpha_min = min ( alpha_min, c_angle )

      alpha_ave = alpha_ave + alpha_min

      alpha_area = alpha_area + area * alpha_min

      area_total = area_total + area

    end do

    alpha_ave = alpha_ave / real ( element_num, real64)
    alpha_area = alpha_area / area_total
  !
  !  Normalize angles from [0,pi/3] radians into qualities in [0,1].
  !
    alpha_min = alpha_min * 3.0e+00_real64 / pi
    alpha_ave = alpha_ave * 3.0e+00_real64 / pi
    alpha_area = alpha_area * 3.0e+00_real64 / pi
  end

  function arc_cosine ( c )

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
  !    Input, real(real64) C, the argument.
  !
  !    Output, real(real64) ARC_COSINE, an angle whose cosine is C.
  !

    real(real64) arc_cosine
    real(real64) c
    real(real64) c2

    c2 = c
    c2 = max ( c2, -1.0e+00_real64 )
    c2 = min ( c2, +1.0e+00_real64 )

    arc_cosine = acos ( c2 )
  end

  subroutine area_measure ( n, z, element_order, element_num, element_node, &
    area_min, area_max, area_ratio, area_ave, area_std, area_negative, area_zero )

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
  !    Input, integer(int32) N, the number of points.
  !
  !    Input, real(real64) Z(2,N), the points.
  !
  !    Input, integer(int32) TRIANGLE_ORDER, the order of the triangles.
  !
  !    Input, integer(int32) TRIANGLE_NUM, the number of triangles.
  !
  !    Input, integer(int32) TRIANGLE_NODE(TRIANGLE_ORDER,TRIANGLE_NUM), 
  !    the triangulation.
  !
  !    Output, real(real64) AREA_MIN, AREA_MAX, the minimum and maximum 
  !    areas.
  !
  !    Output, real(real64) AREA_RATIO, the ratio of the minimum to the
  !    maximum area.
  !
  !    Output, real(real64) AREA_AVE, the average area.
  !
  !    Output, real(real64) AREA_STD, the standard deviation of the areas.
  !
  !    Output, integer(int32) AREA_NEGATIVE, the number of triangles with
  !    negative area.  This suggests an orientation error.
  !
  !    Output, integer(int32) AREA_ZERO, the number of triangles with zero
  !    area.
  !

    integer(int32) n
    integer(int32) element_num
    integer(int32) element_order

    real(real64) area
    real(real64) area_ave
    real(real64) area_max
    real(real64) area_min
    integer(int32) area_negative
    real(real64) area_ratio
    real(real64) area_std
    integer(int32) area_zero
    integer(int32) triangle
    integer(int32) element_node(element_order,element_num)
    real(real64) x1
    real(real64) x2
    real(real64) x3
    real(real64) y1
    real(real64) y2
    real(real64) y3
    real(real64) z(2,n)

    area_max = 0.0e+00_real64
    area_min = huge ( area_min )
    area_ave = 0.0

    area_negative = 0
    area_zero = 0

    do triangle = 1, element_num

      x1 = z(1,element_node(1,triangle))
      y1 = z(2,element_node(1,triangle))
      x2 = z(1,element_node(2,triangle))
      y2 = z(2,element_node(2,triangle))
      x3 = z(1,element_node(3,triangle))
      y3 = z(2,element_node(3,triangle))

      area = 0.5e+00_real64 * ( x1 * ( y2 - y3 ) &
                       + x2 * ( y3 - y1 ) &
                       + x3 * ( y1 - y2 ) )

      if ( area == 0.0e+00_real64 ) then
        area_zero = area_zero + 1
      end if

      if ( area < 0.0e+00_real64 ) then
        area_negative = area_negative + 1
      end if

      area_min = min ( area_min, abs ( area ) )
      area_max = max ( area_max, abs ( area ) )
      area_ave = area_ave + abs ( area )

    end do

    area_ave = area_ave / real ( element_num, real64)

    area_std = 0.0e+00_real64
    do triangle = 1, element_num

      x1 = z(1,element_node(1,triangle))
      y1 = z(2,element_node(1,triangle))
      x2 = z(1,element_node(2,triangle))
      y2 = z(2,element_node(2,triangle))
      x3 = z(1,element_node(3,triangle))
      y3 = z(2,element_node(3,triangle))

      area = 0.5e+00_real64 * abs ( x1 * ( y2 - y3 ) &
                           + x2 * ( y3 - y1 ) &
                           + x3 * ( y1 - y2 ) )

      area_std = area_std + ( area - area_ave )**2
    end do
    area_std = sqrt ( area_std / real ( element_num, real64) )

    if ( 0.0e+00_real64 < area_max ) then
      area_ratio = area_min / area_max
    else
      area_ratio = 0.0e+00_real64
    end if
  end

  subroutine bandwidth_mesh ( element_order, element_num, element_node, &
    ml, mu, m )

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
  !    Input, integer(int32) ELEMENT_ORDER, the order of the elements.
  !
  !    Input, integer(int32) ELEMENT_NUM, the number of elements.
  !
  !    Input, integer(int32) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
  !    ELEMENT_NODE(I,J) is the global index of local node I in element J.
  !
  !    Output, integer(int32) ML, MU, the lower and upper bandwidths 
  !    of the matrix.
  !
  !    Output, integer(int32) M, the bandwidth of the matrix.
  !

    integer(int32) element_num
    integer(int32) element_order

    integer(int32) element
    integer(int32) element_node(element_order,element_num)
    integer(int32) global_i
    integer(int32) global_j
    integer(int32) local_i
    integer(int32) local_j
    integer(int32) m
    integer(int32) ml
    integer(int32) mu

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
  end

  subroutine mesh_base_one ( node_num, element_order, element_num, element_node )

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

    integer(int32) element_num
    integer(int32) element_order

    integer(int32) element
    integer(int32) element_node(element_order,element_num)
    integer(int32) node
    integer(int32) node_max
    integer(int32) node_min
    integer(int32) node_num
    integer(int32) order

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
  end

  subroutine q_measure ( n, z, element_order, element_num, element_node, &
    q_min, q_max, q_ave, q_area )

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
  !    Input, integer(int32) N, the number of points.
  !
  !    Input, real(real64) Z(2,N), the points.
  !
  !    Input, integer(int32) TRIANGLE_ORDER, the order of the triangles.
  !
  !    Input, integer(int32) TRIANGLE_NUM, the number of triangles.
  !
  !    Input, integer(int32) TRIANGLE_NODE(TRIANGLE_ORDER,TRIANGLE_NUM), 
  !    the triangulation.
  !
  !    Output, real(real64) Q_MIN, Q_MAX, the minimum and maximum values
  !    of Q over all triangles.
  !
  !    Output, real(real64) Q_AVE, the average value of Q.
  !
  !    Output, real(real64) Q_AREA, the average value of Q, weighted by
  !    the area of each triangle.
  !

    integer(int32) n
    integer(int32) element_num
    integer(int32) element_order

    integer(int32) a_index
    real(real64) ab_length
    real(real64) area
    real(real64) area_total
    integer(int32) b_index
    real(real64) bc_length
    integer(int32) c_index
    real(real64) ca_length
    real(real64) q
    real(real64) q_area
    real(real64) q_ave
    real(real64) q_max
    real(real64) q_min
    integer(int32) triangle
    integer(int32) element_node(element_order,element_num)
    real(real64) x1
    real(real64) x2
    real(real64) x3
    real(real64) y1
    real(real64) y2
    real(real64) y3
    real(real64) z(2,n)

    q_min =   huge ( q_min )
    q_max = - huge ( q_max )
    q_ave = 0.0e+00_real64
    q_area = 0.0e+00_real64
    area_total = 0.0e+00_real64

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

      area = 0.5e+00_real64 * abs ( x1 * ( y2 - y3 ) &
                           + x2 * ( y3 - y1 ) &
                           + x3 * ( y1 - y2 ) )

      q_min = min ( q_min, q )
      q_max = max ( q_max, q )
      q_ave = q_ave + q
      q_area = q_area + q * area

      area_total = area_total + area

    end do

    q_ave = q_ave / real ( element_num, real64)

    if ( 0.0e+00_real64 < area_total ) then
      q_area = q_area / area_total
    else
      q_area = 0.0e+00_real64
    end if
  end

end module triangulation_quality_mod
