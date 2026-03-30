!> quad_mesh — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

subroutine adj_set_q4_mesh ( node_num, element_num, element_node, &
  element_neighbor, adj_num, adj_col, adj )

!*****************************************************************************80
!
!! ADJ_SET_Q4_MESH sets adjacencies in a Q4 mesh.
!
!  Discussion:
!
!    This routine is called to set the adjacency values, after the
!    appropriate amount of memory has been set aside for storage.
!
!    The mesh is assumed to involve 4-node quadrilaterals.
!
!    Two nodes are "adjacent" if they are both nodes in some element.
!    Also, a node is considered to be adjacent to itself.
!
!    This routine can be used to create the compressed column storage
!    for a linear element finite element discretization of 
!    Poisson's equation in two dimensions.
!
!  Diagram:
!
!         side 3
!       4-------3
!    s  |       |  s
!    i  |       |  i
!    d  |       |  d
!    e  |       |  e
!       |       |
!    4  |       |  2
!       |       |
!       1-------2
!
!         side 1
!
!    The local node numbering
!
!
!   20-21-22-23-24
!    |  |  |  |  |
!    |  |  |  |  |
!   15-16-17-18-19
!    |  |  |  |  |
!    |  |  |  |  |
!   10-11-12-13-14
!    |  |  |  |  |
!    |  |  |  |  |
!    5--6--7--8--9
!    |  |  |  |  |
!    |  |  |  |  |
!    0--1--2--3--4
!
!    A sample grid.
!      
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(4,ELEMENT_NUM), lists the nodes
!    that make up each element in counterclockwise order.
!
!    Input, integer ELEMENT_NEIGHBOR(4,ELEMENT_NUM), for each 
!    side of an element, lists the neighboring element, or -1 if there is
!    no neighbor.
!
!    Input, integer ADJ_NUM, the number of adjacencies.
!
!    Input, integer ADJ_COL(NODE_NUM+1).  Information about 
!    column J is stored in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
!
!    Output, integer ADJ(ADJ_NUM), the adjacency information.
!
  implicit none

  integer adj_num
  integer node_num
  integer element_num
  integer , parameter :: element_order = 4

  integer adj(adj_num)
  integer adj_col(node_num+1)
  integer adj_copy(node_num)
  integer k1
  integer k2
  integer n1
  integer n2
  integer n3
  integer n4
  integer node
  integer number
  integer element
  integer element2
  integer element_neighbor(4,element_num)
  integer element_node(element_order,element_num)

  adj(1:adj_num) = -1
  adj_copy(1:node_num) = adj_col(1:node_num)
!
!  Set every node to be adjacent to itself.
!
  do node = 1, node_num
    adj(adj_copy(node)) = node
    adj_copy(node) = adj_copy(node) + 1
  end do
!
!  Examine each element.
!
  do element = 1, element_num

    n1 = element_node(1,element)
    n2 = element_node(2,element)
    n3 = element_node(3,element)
    n4 = element_node(4,element)
!
!  Add edges (1,3) and (2,4).  There is no need to check for redundancy,
!  since this is the only case when these nodes can share an element.
!
    adj(adj_copy(n1)) = n3
    adj_copy(n1) = adj_copy(n1) + 1
    adj(adj_copy(n3)) = n1
    adj_copy(n3) = adj_copy(n3) + 1

    adj(adj_copy(n2)) = n4
    adj_copy(n2) = adj_copy(n2) + 1
    adj(adj_copy(n4)) = n2
    adj_copy(n4) = adj_copy(n4) + 1
!
!  Add edge (1,2) if this is the first occurrence,
!  that is, if the edge (1,2) is on a boundary (ELEMENT2 <= 0)
!  or if this element is the first of the pair in which the edge
!  occurs (ELEMENT < ELEMENT2).
!
    element2 = element_neighbor(1,element)

    if ( element2 < 0 .or. element < element2 ) then
      adj(adj_copy(n1)) = n2
      adj_copy(n1) = adj_copy(n1) + 1
      adj(adj_copy(n2)) = n1
      adj_copy(n2) = adj_copy(n2) + 1
    end if
!
!  Add edge (2,3).
!
    element2 = element_neighbor(2,element)

    if ( element2 < 0 .or. element < element2 ) then
      adj(adj_copy(n2)) = n3
      adj_copy(n2) = adj_copy(n2) + 1
      adj(adj_copy(n3)) = n2
      adj_copy(n3) = adj_copy(n3) + 1
    end if
!
!  Add edge (3,4).
!
    element2 = element_neighbor(3,element)

    if ( element2 < 0 .or. element < element2 ) then
      adj(adj_copy(n4)) = n3
      adj_copy(n4) = adj_copy(n4) + 1
      adj(adj_copy(n3)) = n4
      adj_copy(n3) = adj_copy(n3) + 1
    end if
!
!  Add edge (4,1).
!
    element2 = element_neighbor(4,element)

    if ( element2 < 0 .or. element < element2 ) then
      adj(adj_copy(n1)) = n4
      adj_copy(n1) = adj_copy(n1) + 1
      adj(adj_copy(n4)) = n1
      adj_copy(n4) = adj_copy(n4) + 1
    end if
      
  end do
!
!  Ascending sort the entries for each node.
!
  do node = 1, node_num
    k1 = adj_col(node)
    k2 = adj_col(node+1)-1
    number = k2 + 1 - k1
    call i4vec_sort_heap_a ( number, adj(k1:k2) )
  end do
end

subroutine adj_size_q4_mesh ( node_num, element_num, element_node, &
  element_neighbor, adj_num, adj_col )

!*****************************************************************************80
!
!! ADJ_SIZE_Q4_MESH counts adjacencies in a Q4 mesh.
!
!  Discussion:
!
!    This routine is called to count the adjacencies, so that the
!    appropriate amount of memory can be set aside for storage when
!    the adjacency structure is created.
!
!    The mesh is assumed to involve 4-node quadrilaterals.
!
!    Two nodes are "adjacent" if they are both nodes in some quadrilateral.
!    Also, a node is considered to be adjacent to itself.
!
!  Diagram:
!
!         side 3
!       4-------3
!    s  |       |  s
!    i  |       |  i
!    d  |       |  d
!    e  |       |  e
!       |       |
!    4  |       |  2
!       |       |
!       1-------2
!
!         side 1
!
!    The local node numbering
!
!
!   20-21-22-23-24
!    |  |  |  |  |
!    |  |  |  |  |
!   15-16-17-18-19
!    |  |  |  |  |
!    |  |  |  |  |
!   10-11-12-13-14
!    |  |  |  |  |
!    |  |  |  |  |
!    5--6--7--8--9
!    |  |  |  |  |
!    |  |  |  |  |
!    0--1--2--3--4
!
!    A sample grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(4,ELEMENT_NUM), lists the 
!    nodes that make up each element, in counterclockwise order.
!
!    Input, integer ELEMENT_NEIGHBOR(4,ELEMENT_NUM), for each 
!    side of a element, lists the neighboring elment, or -1 if there is
!    no neighbor.
!
!    Output, integer ADJ_NUM, the number of adjacencies.
!
!    Output, integer ADJ_COL(NODE_NUM+1), Information about 
!    column J is stored in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
!
  implicit none

  integer element_num
  integer , parameter :: element_order = 4
  integer node_num

  integer adj_col(node_num+1)
  integer adj_num
  integer element
  integer element_neighbor(4,element_num)
  integer element_node(element_order,element_num)
  integer element2
  integer i
  integer n1
  integer n2
  integer n3
  integer n4
  integer node

  adj_num = 0
!
!  Set every node to be adjacent to itself.
!
  adj_col(1:node_num) = 1
!
!  Examine each element.
!
  do element = 1, element_num

    n1 = element_node(1,element)
    n2 = element_node(2,element)
    n3 = element_node(3,element)
    n4 = element_node(4,element)
!
!  Add edge (1,3).
!
    adj_col(n1) = adj_col(n1) + 1
    adj_col(n3) = adj_col(n3) + 1
!
!  Add edge (2,4).
!
    adj_col(n2) = adj_col(n2) + 1
    adj_col(n4) = adj_col(n4) + 1
!
!  Add edge (1,2) if this is the first occurrence,
!  that is, if the edge (1,2) is on a boundary (ELEMENT2 <= 0)
!  or if this element is the first of the pair in which the edge
!  occurs (ELEMENT < ELEMENT2).
!
    element2 = element_neighbor(1,element)

    if ( element2 < 0 .or. element < element2 ) then
      adj_col(n1) = adj_col(n1) + 1
      adj_col(n2) = adj_col(n2) + 1
    end if
!
!  Add edge (2,3).
!
    element2 = element_neighbor(2,element)

    if ( element2 < 0 .or. element < element2 ) then
      adj_col(n2) = adj_col(n2) + 1
      adj_col(n3) = adj_col(n3) + 1
    end if
!
!  Add edge (3,4).
!
    element2 = element_neighbor(3,element)

    if ( element2 < 0 .or. element < element2 ) then
      adj_col(n3) = adj_col(n3) + 1
      adj_col(n4) = adj_col(n4) + 1
    end if
!
!  Add edge (4,1).
!
    element2 = element_neighbor(4,element)

    if ( element2 < 0 .or. element < element2 ) then
      adj_col(n4) = adj_col(n4) + 1
      adj_col(n1) = adj_col(n1) + 1
    end if

  end do
!
!  We used ADJ_COL to count the number of entries in each column.
!  Convert it to pointers into the ADJ array.
!
  do node = node_num + 1, 2, -1
    adj_col(node) = adj_col(node-1)
  end do

  adj_col(1) = 1
  do node = 2, node_num + 1
    adj_col(node) = adj_col(node) + adj_col(node-1)
  end do
!
!  Finally, record the total number of adjacencies.
!
  adj_num = adj_col(node_num+1) - 1
end

subroutine area_q4_mesh ( node_num, element_num, node_xy, element_node, &
  element_area, mesh_area )

!*****************************************************************************80
!
!! AREA_Q4_MESH computes areas of elements in a Q4 mesh.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, double precision NODE_XY(2,NODE_NUM), the node coordinates.
!
!    Input, integer ELEMENT_NODE(4,ELEMENT_NUM), lists the 
!    nodes that make up each element, in counterclockwise order.
!
!    Output, double precision ELEMENT_AREA(ELEMENT_NUM), the element areas.
!
!    Output, double precision MESH_AREA, the mesh area.
!
  implicit none

  integer element_num
  integer node_num

  integer element
  double precision element_area(element_num)
  integer element_node(4,element_num)
  double precision mesh_area
  integer node
  double precision node_xy(2,node_num)
  double precision q4(2,4)

  do element = 1, element_num
    do node = 1, 4
      q4(1:2,node) = node_xy(1:2,element_node(node,element))
    end do
    call area_quad ( q4, element_area(element) )
  end do

  mesh_area = sum ( element_area(1:element_num) )
end

subroutine area_quad ( quad_xy, area )

!*****************************************************************************80
!
!! AREA_QUAD returns the area of a quadrilateral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision QUAD_XY(2,4), the coordinates of the nodes.
!
!    Output, double precision AREA, the area of the quadrilateral.
!
  implicit none

  double precision area
  double precision area1
  double precision area2
  double precision quad_xy(2,4)
  double precision t1(2,3)
  double precision t2(2,3)

  t1(1:2,1) = quad_xy(1:2,1)
  t1(1:2,2) = quad_xy(1:2,2)
  t1(1:2,3) = quad_xy(1:2,3)

  call triangle_area ( t1, area1 )

  t2(1:2,1) = quad_xy(1:2,3)
  t2(1:2,2) = quad_xy(1:2,4)
  t2(1:2,3) = quad_xy(1:2,1)

  call triangle_area ( t2, area2 )

  if ( area1 < 0.0D+00 .or. area2 < 0.0D+00 ) then

    t1(1:2,1) = quad_xy(1:2,2)
    t1(1:2,2) = quad_xy(1:2,3)
    t1(1:2,3) = quad_xy(1:2,4)

    call triangle_area ( t1, area1 )

    t2(1:2,1) = quad_xy(1:2,4)
    t2(1:2,2) = quad_xy(1:2,1)
    t2(1:2,3) = quad_xy(1:2,2)

    call triangle_area ( t2, area2 )

    if ( area1 < 0.0D+00 .or. area2 < 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'AREA_QUAD - Fatal error!'
      write ( *, '(a)' ) '  The quadrilateral nodes seem to be listed in'
      write ( *, '(a)' ) '  the wrong order, or the quadrilateral is'
      write ( *, '(a)' ) '  degenerate.'
      stop
    end if

  end if

  area = area1 + area2
end

subroutine boundary_edge_count_q4_mesh ( element_num, element_node, &
  boundary_edge_num )

!*****************************************************************************80
!
!! BOUNDARY_EDGE_COUNT_Q4_MESH counts the boundary edges.
!
!  Discussion:
!
!    This routine is given a Q4 mesh, an abstract list of sets of 4 nodes.
!    It is assumed that the nodes in each Q4 are listed
!    in a counterclockwise order, although the routine should work 
!    if the nodes are consistently listed in a clockwise order as well.
!
!    It is assumed that each edge of the mesh is either 
!    * an INTERIOR edge, which is listed twice, once with positive
!      orientation and once with negative orientation, or;
!    * a BOUNDARY edge, which will occur only once.
!
!    This routine should work even if the region has holes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(4,ELEMENT_NUM), the nodes 
!    that make up the elements.  These should be listed in counterclockwise 
!    order.
!
!    Output, integer BOUNDARY_EDGE_NUM, the number of boundary 
!    edges.
!
  implicit none

  integer element_num
  integer , parameter :: element_order = 4

  integer boundary_edge_num
  integer e1(4*element_num)
  integer e2(4*element_num)
  integer edge(2,4*element_num)
  integer interior_edge_num
  integer m
  integer n
  integer element_node(element_order,element_num)
  integer unique_num

  m = 2
  n = 4 * element_num
!
!  Set up the edge array.
!
  edge(1:2,              1:  element_num) = element_node(1:2,1:element_num)
  edge(1:2,  element_num+1:2*element_num) = element_node(2:3,1:element_num)
  edge(1:2,2*element_num+1:3*element_num) = element_node(3:4,1:element_num)

  edge(1  ,3*element_num+1:4*element_num) = element_node(4,  1:element_num)
  edge(2  ,3*element_num+1:4*element_num) = element_node(1,  1:element_num)
!
!  In each column, force the smaller entry to appear first.
!
  e1(1:n) = minval ( edge(1:2,1:n), dim = 1 )
  e2(1:n) = maxval ( edge(1:2,1:n), dim = 1 )

  edge(1,1:n) = e1(1:n)
  edge(2,1:n) = e2(1:n)
!
!  Ascending sort the column array.
!
  call i4col_sort_a ( m, n, edge )
!
!  Get the number of unique columns in EDGE.
!
  call i4col_sorted_unique_count ( m, n, edge, unique_num )

  interior_edge_num = 4 * element_num - unique_num

  boundary_edge_num = 4 * element_num - 2 * interior_edge_num
end

subroutine boundary_edge_count_euler_q4_mesh ( node_num, element_num, &
  hole_num, boundary_num )

!*****************************************************************************80
!
!! BOUNDARY_EDGE_COUNT_EULER_Q4_MESH counts boundary edges by Euler's formula.
!
!  Discussion:
!
!    We assume we are given information about a quadrilateral mesh
!    of a set of nodes in the plane.
!
!    Given the number of nodes, elements and holes, we are going to apply
!    Euler's formula to determine the number of edges that lie on the
!    boundary of the set of nodes.
!
!    The number of faces, including the infinite face and internal holes, 
!    is ELEMENT_NUM + HOLE_NUM + 1.
!
!    Let BOUNDARY_NUM denote the number of edges on the boundary.
!    Each of the ELEMENT_NUM quadrilaterals uses four edges.  Every edge
!    occurs in two different elements, so the number of edges must be
!    ( 4 * ELEMENT_NUM + BOUNDARY_NUM ) / 2.
!
!    The number of nodes used in the mesh is NODE_NUM.
!
!    Euler's formula asserts that, for a simple connected figure in the
!    plane with no edge crossings, NODE_NUM nodes, EDGE_NUM edges and
!    FACE_NUM faces:
!
!      NODE_NUM - EDGE_NUM + FACE_NUM = 2
!
!    In our context, this becomes
!
!      NODE_NUM - ( 4 * ELEMENT_NUM + BOUNDARY_NUM ) / 2 
!      + ELEMENT_NUM + HOLE_NUM + 1 = 2
!
!    or
!
!      BOUNDARY_NUM = 2 * NODE_NUM + 2 * HOLE_NUM - 2 * ELEMENT_NUM - 2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Marc de Berg, Marc Krevald, Mark Overmars, Otfried Schwarzkopf,
!    Computational Geometry, Section 9.1,
!    Springer, 2000.
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer HOLE_NUM, the number of internal holes.
!
!    Output, integer BOUNDARY_NUM, the number of edges that 
!    lie on the boundary of the mesh.
!
  implicit none

  integer boundary_num
  integer element_num
  integer hole_num
  integer node_num

  boundary_num = 2 * node_num + 2 * hole_num - 2 * element_num - 2
end

subroutine example1_q4_mesh ( node_num, element_num, node_xy, element_node, &
  element_neighbor )

!*****************************************************************************80
!
!! EXAMPLE1_Q4_MESH sets up example #1 Q4 mesh.
!
!  Discussion:
!
!    The appropriate values of NODE_NUM and ELEMENT_NUM can be found by
!    calling EXAMPLE1_Q4_MESH_SIZE first.
!
!   24---25---26---27---28
!    | 14 | 15 | 16 | 17 |
!   18---19---20---21---22---23
!    | 10 | -2 | 11 | 12 | 13 |
!   12---13---14---15---16---17
!    |  5 |  6 |  7 |  8 |  9 |
!    6----7----8----9---10---11
!    |  1 |  2 |  3 |  4 |
!    1----2----3----4----5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.  
!
!    Input, integer ELEMENT_NUM, the number of elements.  
!
!    Output, double precision NODE_XY(2,NODE_NUM), the coordinates of the
!    nodes.
!
!    Output, integer ELEMENT_NODE(4,ELEMENT_NUM), the nodes
!    that make up the elements.
!
!    Output, integer ELEMENT_NEIGHBOR(4,ELEMENT_NUM), the 
!    element neighbors on each side.  Negative values indicate edges that 
!    lie on the exterior.
!
  implicit none

  integer , parameter :: dim_num = 2
  integer element_num
  integer , parameter :: element_order = 4
  integer node_num

  double precision node_xy(dim_num,node_num)
  integer element_node(element_order,element_num)
  integer element_neighbor(4,element_num)

  node_xy = reshape ( (/ &
       0.0D+00, 0.0D+00, &
       1.0D+00, 0.0D+00, &
       2.0D+00, 0.0D+00, &
       3.0D+00, 0.0D+00, &
       4.0D+00, 0.0D+00, &
       0.0D+00, 1.0D+00, &
       1.0D+00, 1.0D+00, &
       2.0D+00, 1.0D+00, &
       3.0D+00, 1.0D+00, &
       4.0D+00, 1.0D+00, &
       5.0D+00, 1.0D+00, &
       0.0D+00, 2.0D+00, &
       1.0D+00, 2.0D+00, &
       2.0D+00, 2.0D+00, &
       3.0D+00, 2.0D+00, &
       4.0D+00, 2.0D+00, &
       5.0D+00, 2.0D+00, &
       0.0D+00, 3.0D+00, &
       1.0D+00, 3.0D+00, &
       2.0D+00, 3.0D+00, &
       3.0D+00, 3.0D+00, &
       4.0D+00, 3.0D+00, &
       5.0D+00, 3.0D+00, &
       0.0D+00, 4.0D+00, &
       1.0D+00, 4.0D+00, &
       2.0D+00, 4.0D+00, &
       3.0D+00, 4.0D+00, &
       4.0D+00, 4.0D+00 /), (/ dim_num, node_num /) )

  element_node(1:element_order,1:element_num ) = reshape ( (/ &
     1,  2,  7,  6, &
     2,  3,  8,  7, &
     3,  4,  9,  8, &
     4,  5, 10,  9, &
     6,  7, 13, 12, &
     7,  8, 14, 13, &
     8,  9, 15, 14, &
     9, 10, 16, 15, &
    10, 11, 17, 16, &
    12, 13, 19, 18, &
    14, 15, 21, 20, &
    15, 16, 22, 21, &
    16, 17, 23, 22, &
    18, 19, 25, 24, &
    19, 20, 26, 25, &
    20, 21, 27, 26, &
    21, 22, 28, 27 /), (/ element_order, element_num /) )

  element_neighbor(1:4,1:element_num) = reshape ( (/ &
       -1,  2,  5, -1, &
       -1,  3,  6,  1, &
       -1,  4,  7,  2, &
       -1, -1,  8,  3, &
        1,  6, 10, -1, &
        2,  7, -2,  5, &
        3,  8, 11,  6, &
        4,  9, 12,  7, &
       -1, -1, 13,  8, &
        5, -2, 14, -1, &
        7, 12, 16, -2, &
        8, 13, 17, 11, &
        9, -1, -1, 12, &
       10, 15, -1, -1, &
       -2, 16, -1, 14, &
       11, 17, -1, 15, &
       12, -1, -1, 16 /), (/ 4, element_num /) )
end

subroutine example1_q4_mesh_size ( node_num, element_num, hole_num )

!*****************************************************************************80
!
!! EXAMPLE1_Q4_MESH_SIZE sets sizes for example #1 Q4 mesh.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer NODE_NUM, the number of nodes.  
!
!    Output, integer ELEMENT_NUM, the number of elements. 
!
!    Output, integer HOLE_NUM, the number of holes.
!
  implicit none

  integer element_num
  integer hole_num
  integer node_num

  element_num = 17
  hole_num = 1
  node_num = 28
end

subroutine example2_q4_mesh ( node_num, element_num, node_xy, element_node, &
  element_neighbor )

!*****************************************************************************80
!
!! EXAMPLE2_Q4_MESH sets up example #2 Q4 mesh.
!
!  Discussion:
!
!    The region is a semicircle.  This example includes degenerate elements
!    (the first layer of elements is touching the origin, and so has a side
!    of length zero).  The elements are not parallelograms.  And the elements
!    vary in size.
!
!    Because of the treatment of node 1, algorithms for counting boundary 
!    edges may become "confused".
!
!    The appropriate values of NODE_NUM and ELEMENT_NUM can be found by
!    calling EXAMPLE1_Q4_MESH_SIZE first.
!
!   29---30---31---32---33---34---35---36---37
!    | 25 | 26 | 27 | 28 | 29 | 30 | 31 | 32 |
!   20---21---22---23---24---25---26---27---28
!    | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 |
!   11---12---13---14---15---16---17---18---19
!    |  9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |
!    2----3----4----5----6----7----8----9---10
!    |  1 |  2 |  3 |  4 |  5 |  6 |  7 |  8 |
!    1----1----1----1----1----1----1----1----1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.  
!
!    Input, integer ELEMENT_NUM, the number of elements.  
!
!    Output, double precision NODE_XY(2,NODE_NUM), the coordinates of the
!    nodes.
!
!    Output, integer ELEMENT_NODE(4,ELEMENT_NUM), the nodes
!    that make up the elements.
!
!    Output, integer ELEMENT_NEIGHBOR(4,ELEMENT_NUM), the 
!    element neighbors on each side.  Negative values indicate edges that 
!    lie on the exterior.
!
  implicit none

  integer , parameter :: dim_num = 2
  integer element_num
  integer , parameter :: element_order = 4
  integer node_num

  double precision a
  integer col
  integer element
  integer k
  integer element_node(element_order,element_num)
  integer element_neighbor(4,element_num)
  double precision node_xy(dim_num,node_num)
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision r
  integer row

  k = 1
  node_xy(1,k) = 0.0D+00
  node_xy(2,k) = 0.0D+00

  do row = 1, 4
    r = real ( row)
    do col = 0, 8
      a = real ( 8 - col) * pi / 8.0D+00
      k = k + 1
      node_xy(1,k) = r * cos ( a )
      node_xy(2,k) = r * sin ( a )
    end do
  end do

  element = 0
  do row = 0, 3
    do col = 0, 7
      element = element + 1
      if ( row == 0 ) then
        element_node(1,element) = 1
        element_node(2,element) = 1
        element_node(3,element) = col + 3
        element_node(4,element) = col + 2
      else
        element_node(1,element) = element_node(4,element-8)
        element_node(2,element) = element_node(3,element-8)
        element_node(3,element) = element_node(2,element) + 9
        element_node(4,element) = element_node(1,element) + 9
      end if
    end do
  end do

  element = 0
  do row = 0, 3
    do col = 0, 7
      element = element + 1
      if ( row == 0 ) then
        element_neighbor(1,element) = -1
      else
        element_neighbor(1,element) = element - 8
      end if
      if ( col == 7 ) then
        element_neighbor(2,element) = -1
      else
        element_neighbor(2,element) = element + 1
      end if
      if ( row == 3 ) then
        element_neighbor(3,element) = - 1
      else 
        element_neighbor(3,element) = element + 8
      end if
      if ( col == 0 ) then
        element_neighbor(4,element) = - 1
      else
        element_neighbor(4,element) = element - 1
      end if
    end do
  end do
end

subroutine example2_q4_mesh_size ( node_num, element_num, hole_num )

!*****************************************************************************80
!
!! EXAMPLE2_Q4_MESH_SIZE sets sizes for example #2 Q4 mesh.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer NODE_NUM, the number of nodes.  
!
!    Output, integer ELEMENT_NUM, the number of elements. 
!
!    Output, integer HOLE_NUM, the number of holes.
!
  implicit none

  integer element_num
  integer hole_num
  integer node_num

  element_num = 32
  hole_num = 0
  node_num = 37
end

subroutine neighbor_elements_q4_mesh ( element_num, element_node, &
  element_neighbor )

!*****************************************************************************80
!
!! NEIGHBOR_ELEMENTS_Q4_MESH determines element neighbors in a Q4 mesh.
!
!  Discussion:
!
!    A quadrilateral mesh of a set of nodes can be completely described by
!    the coordinates of the nodes, and the list of nodes that make up
!    each element.  However, in some cases, it is necessary to know
!    element adjacency information, that is, which element, if any,
!    is adjacent to a given element on a particular side.
!
!    This routine creates a data structure recording this information.
!
!    The primary amount of work occurs in sorting a list of 4 * ELEMENT_NUM
!    data items.
!
!    Note that COL is a work array allocated dynamically inside this
!    routine.  It is possible, for very large values of ELEMENT_NUM,
!    that the necessary amount of memory will not be accessible, and the
!    routine will fail.  This is a limitation of the implementation of
!    dynamic arrays in FORTRAN90.  One way to get around this would be
!    to require the user to declare ROW in the calling routine
!    as an allocatable array, get the necessary memory explicitly with
!    an ALLOCATE statement, and then pass ROW into this routine.
!
!    Of course, the point of dynamic arrays was to make it easy to
!    hide these sorts of temporary work arrays from the poor user!
!
!    This routine was revised to store the edge data in a column
!    array rather than a row array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(4,ELEMENT_NUM), the nodes 
!    that make up each element.
!
!    Output, integer ELEMENT_NEIGHBOR(4,ELEMENT_NUM), lists the
!    neighboring element on each side of a given element, or -1 if there is
!    no neighbor.
!
  implicit none

  integer element_num
  integer , parameter :: element_order = 4

  integer , allocatable :: col(:,:)
  integer element_neighbor(4,element_num)
  integer element_node(element_order,element_num)
  integer i
  integer icol
  integer j
  integer k
  integer l
  integer q
  integer q1
  integer q2
  integer side1
  integer side2

  allocate ( col (4,4*element_num) )
!
!  Step 1.
!  From the list of nodes for element Q, of the form: (I,J,K,L)
!  construct the four neighbor relations:
!
!    (I,J,1,Q) or (J,I,1,Q),
!    (J,K,2,Q) or (K,J,2,Q),
!    (K,L,3,Q) or (L,K,3,Q)
!    (K,I,4,Q) or (I,K,4,Q)
!
!  where we choose (I,J,1,Q) if I < J, or else (J,I,1,Q)
!
  do q = 1, element_num

    i = element_node(1,q)
    j = element_node(2,q)
    k = element_node(3,q)
    l = element_node(4,q)

    if ( i < j ) then
      col(1:4,4*(q-1)+1) = (/ i, j, 1, q /)
    else
      col(1:4,4*(q-1)+1) = (/ j, i, 1, q /)
    end if

    if ( j < k ) then
      col(1:4,4*(q-1)+2) = (/ j, k, 2, q /)
    else
      col(1:4,4*(q-1)+2) = (/ k, j, 2, q /)
    end if

    if ( k < l ) then
      col(1:4,4*(q-1)+3) = (/ k, l, 3, q /)
    else
      col(1:4,4*(q-1)+3) = (/ l, k, 3, q /)
    end if

    if ( l < i ) then
      col(1:4,4*(q-1)+4) = (/ l, i, 4, q /)
    else
      col(1:4,4*(q-1)+4) = (/ i, l, 4, q /)
    end if

  end do
!
!  Step 2. Perform an ascending dictionary sort on the neighbor relations.
!  We only intend to sort on rows 1 and 2; the routine we call here
!  sorts on rows 1 through 4 but that won't hurt us.
!
!  What we need is to find cases where two elements share an edge.
!  Say they share an edge defined by the nodes I and J.  Then there are
!  two columns of COL that start out ( I, J, ?, ? ).  By sorting COL,
!  we make sure that these two columns occur consecutively.  That will
!  make it easy to notice that the elements are neighbors.
!
  call i4col_sort_a ( 4, 4*element_num, col )
!
!  Step 3. Neighboring elements show up as consecutive columns with
!  identical first two entries.  Whenever you spot this happening,
!  make the appropriate entries in ELEMENT_NEIGHBOR.
!
  element_neighbor(1:4,1:element_num) = -1

  icol = 1

  do

    if ( 4 * element_num <= icol ) then
      exit
    end if

    if ( col(1,icol) /= col(1,icol+1) .or. col(2,icol) /= col(2,icol+1) ) then
      icol = icol + 1
      cycle
    end if

    side1 = col(3,icol)
    q1    = col(4,icol)
    side2 = col(3,icol+1)
    q2    = col(4,icol+1)

    element_neighbor(side1,q1) = q2
    element_neighbor(side2,q2) = q1

    icol = icol + 2

  end do

  deallocate ( col )
end

subroutine node_order_q4_mesh ( element_num, element_node, node_num, &
  node_order )

!*****************************************************************************80
!
!! NODE_ORDER_Q4_MESH determines the order of nodes in a Q4 mesh.
!
!  Discussion:
!
!    The order of a node is the number of elements that use that node
!    as a vertex.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(4,ELEMENT_NUM), 
!    the nodes that make up the elements.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Output, integer NODE_ORDER(NODE_NUM), the order of each node.
!
  implicit none

  integer element_num
  integer , parameter :: element_order = 4
  integer node_num

  integer element
  integer element_node(element_order,element_num)
  integer i
  integer node
  integer node_order(node_num)

  node_order(1:node_num) = 0

  do element = 1, element_num
    do i = 1, element_order
      node = element_node(i,element)
      if ( node < 1 .or. node_num < node ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'NODE_ORDER_Q4_MESH - Fatal error!'
        write ( *, '(a)' ) '  Illegal entry in ELEMENT_NODE.'
        stop
      else
        node_order(node) = node_order(node) + 1
      end if
    end do
  end do
end

subroutine plot_q4_mesh ( node_num, element_num, node_xy, element_node, &
  node_show, element_show, output_filename )

!*****************************************************************************80
!
!! PLOT_Q4_MESH plots a Q4 mesh.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2009
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
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(4,ELEMENT_NUM), the nodes
!    that form the elements.
!
!    Input, integer NODE_SHOW,
!    0, do not show nodes;
!    1, show nodes;
!    2, show nodes and label them.
!
!    Input, integer ELEMENT_SHOW,
!    0, do not show elements;
!    1, show elements;
!    2, show elements and label them.
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the name of the output file.
!
  implicit none

  integer node_num
  integer element_num
  integer , parameter :: element_order = 4

  double precision ave_x
  double precision ave_y
  character ( len = 40 ) date_time
  integer :: circle_size
  integer delta
  integer e
  integer element
  integer element_node(element_order,element_num)
  integer element_show
  integer i
  integer i4_wrap
  integer ios
  integer node
  integer node_show
  double precision node_xy(2,node_num)
  character ( len = * )  output_filename
  integer output_unit
  character ( len = 40 ) string
  double precision x_max
  double precision x_min
  integer x_ps
  integer :: x_ps_max = 576
  integer :: x_ps_max_clip = 594
  integer :: x_ps_min = 36
  integer :: x_ps_min_clip = 18
  double precision x_scale
  double precision y_max
  double precision y_min
  integer y_ps
  integer :: y_ps_max = 666
  integer :: y_ps_max_clip = 684
  integer :: y_ps_min = 126
  integer :: y_ps_min_clip = 108
  double precision y_scale
!
!  We need to do some figuring here, so that we can determine
!  the range of the data, and hence the height and width
!  of the piece of paper.
!
  x_max = maxval ( node_xy(1,1:node_num) )
  x_min = minval ( node_xy(1,1:node_num) )
  x_scale = x_max - x_min

  x_max = x_max + 0.05D+00 * x_scale
  x_min = x_min - 0.05D+00 * x_scale
  x_scale = x_max - x_min

  y_max = maxval ( node_xy(2,1:node_num) )
  y_min = minval ( node_xy(2,1:node_num) )
  y_scale = y_max - y_min

  y_max = y_max + 0.05D+00 * y_scale
  y_min = y_min - 0.05D+00 * y_scale
  y_scale = y_max - y_min

  if ( x_scale < y_scale ) then

    delta = nint ( real ( x_ps_max - x_ps_min) &
      * ( y_scale - x_scale ) / ( 2.0D+00 * y_scale ) )

    x_ps_max = x_ps_max - delta
    x_ps_min = x_ps_min + delta

    x_ps_max_clip = x_ps_max_clip - delta
    x_ps_min_clip = x_ps_min_clip + delta

    x_scale = y_scale

  else if ( y_scale < x_scale ) then

    delta = nint ( real ( y_ps_max - y_ps_min) &
      * ( x_scale - y_scale ) / ( 2.0D+00 * x_scale ) )

    y_ps_max      = y_ps_max - delta
    y_ps_min      = y_ps_min + delta

    y_ps_max_clip = y_ps_max_clip - delta
    y_ps_min_clip = y_ps_min_clip + delta

    y_scale = x_scale

  end if

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLOT_Q4_MESH - Fatal error!'
    write ( *, '(a)' ) '  Can not open output file.'
  end if

  write ( output_unit, '(a)' ) '%!PS-Adobe-3.0 EPSF-3.0'
  write ( output_unit, '(a)' ) '%%Creator: plot_q4_mesh.f90'
  write ( output_unit, '(a)' ) '%%Title: ' // trim ( output_filename )
  write ( output_unit, '(a)' ) '%%Pages: 1'
  write ( output_unit, '(a,i3,2x,i3,2x,i3,2x,i3)' ) '%%BoundingBox: ', &
    x_ps_min, y_ps_min, x_ps_max, y_ps_max
  write ( output_unit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( output_unit, '(a)' ) '%%LanguageLevel: 1'
  write ( output_unit, '(a)' ) '%%EndComments'
  write ( output_unit, '(a)' ) '%%BeginProlog'
  write ( output_unit, '(a)' ) '/inch {72 mul} def'
  write ( output_unit, '(a)' ) '%%EndProlog'
  write ( output_unit, '(a)' ) '%%Page: 1 1'
  write ( output_unit, '(a)' ) 'save'
  write ( output_unit, '(a)' ) '%'
  write ( output_unit, '(a)' ) '%  Set the RGB line color to very light gray.'
  write ( output_unit, '(a)' ) '%'
  write ( output_unit, '(a)' ) '0.900  0.900  0.900 setrgbcolor'
  write ( output_unit, '(a)' ) '%'
  write ( output_unit, '(a)' ) '%  Draw a gray border around the page.'
  write ( output_unit, '(a)' ) '%'
  write ( output_unit, '(a)' ) 'newpath'
  write ( output_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_min, ' moveto'
  write ( output_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_max, y_ps_min, ' lineto'
  write ( output_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_max, y_ps_max, ' lineto'
  write ( output_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_max, ' lineto'
  write ( output_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_min, ' lineto'
  write ( output_unit, '(a)' ) 'stroke'
  write ( output_unit, '(a)' ) '%'
  write ( output_unit, '(a)' ) '%  Set the RGB color to black.'
  write ( output_unit, '(a)' ) '%'
  write ( output_unit, '(a)' ) '0.000  0.000  0.000 setrgbcolor'
  write ( output_unit, '(a)' ) '%'
  write ( output_unit, '(a)' ) '%  Set the font and its size.'
  write ( output_unit, '(a)' ) '%'
  write ( output_unit, '(a)' ) '/Times-Roman findfont'
  write ( output_unit, '(a)' ) '0.50 inch scalefont'
  write ( output_unit, '(a)' ) 'setfont'
  write ( output_unit, '(a)' ) '%'
  write ( output_unit, '(a)' ) '%  Print a title.'
  write ( output_unit, '(a)' ) '%'
  write ( output_unit, '(a)' ) '%  210  702  moveto'
  write ( output_unit, '(a)' ) '%  (Q4 Mesh)  show'
  write ( output_unit, '(a)' ) '%'
  write ( output_unit, '(a)' ) '%  Define a clipping polygon.'
  write ( output_unit, '(a)' ) '%'
  write ( output_unit, '(a)' ) 'newpath'
  write ( output_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_min_clip, ' moveto'
  write ( output_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_max_clip, y_ps_min_clip, ' lineto'
  write ( output_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_max_clip, y_ps_max_clip, ' lineto'
  write ( output_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_max_clip, ' lineto'
  write ( output_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_min_clip, ' lineto'
  write ( output_unit, '(a)' ) 'clip newpath'
!
!  Draw the nodes.
!
  if ( node_num <= 200 ) then
    circle_size = 5
  else if ( node_num <= 500 ) then
    circle_size = 4
  else if ( node_num <= 1000 ) then
    circle_size = 3
  else if ( node_num <= 5000 ) then
    circle_size = 2
  else
    circle_size = 1
  end if

  if ( 1 <= node_show ) then
    write ( output_unit, '(a)' ) '%'
    write ( output_unit, '(a)' ) '%  Draw filled dots at the nodes.'
    write ( output_unit, '(a)' ) '%'
    write ( output_unit, '(a)' ) '%  Set the RGB color to blue.'
    write ( output_unit, '(a)' ) '%'
    write ( output_unit, '(a)' ) '0.000  0.150  0.750 setrgbcolor'
    write ( output_unit, '(a)' ) '%'

    do node = 1, node_num

      x_ps = int ( &
        ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min)   &
        + (         node_xy(1,node) - x_min ) * real ( x_ps_max) ) &
        / ( x_max                   - x_min ) )

      y_ps = int ( &
        ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min)   &
        + (         node_xy(2,node) - y_min ) * real ( y_ps_max) ) &
        / ( y_max                   - y_min ) )

      write ( output_unit, '(a,i4,2x,i4,2x,i4,2x,a)' ) 'newpath ', x_ps, y_ps, &
        circle_size, '0 360 arc closepath fill'

    end do

  end if
!
!  Label the nodes.
!
  if ( 2 <= node_show ) then

    write ( output_unit, '(a)' ) '%'
    write ( output_unit, '(a)' ) '%  Label the nodes:'
    write ( output_unit, '(a)' ) '%'
    write ( output_unit, '(a)' ) '%  Set the RGB color to darker blue.'
    write ( output_unit, '(a)' ) '%'
    write ( output_unit, '(a)' ) '0.000  0.250  0.850 setrgbcolor'
    write ( output_unit, '(a)' ) '/Times-Roman findfont'
    write ( output_unit, '(a)' ) '0.20 inch scalefont'
    write ( output_unit, '(a)' ) 'setfont'
    write ( output_unit, '(a)' ) '%'

    do node = 1, node_num

      x_ps = int ( &
        ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min)   &
        + (       + node_xy(1,node) - x_min ) * real ( x_ps_max) ) &
        / ( x_max                   - x_min ) )

      y_ps = int ( &
        ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min)   &
        + (         node_xy(2,node) - y_min ) * real ( y_ps_max) ) &
        / ( y_max                   - y_min ) )

      write ( string, '(i4)' ) node
      string = adjustl ( string )

      write ( output_unit, '(i4,2x,i4,a)' ) x_ps, y_ps+5, &
        ' moveto (' // trim ( string ) // ') show'

    end do

  end if
!
!  Draw the elements.
!
  if ( 1 <= element_show ) then
    write ( output_unit, '(a)' ) '%'
    write ( output_unit, '(a)' ) '%  Set the RGB color to red.'
    write ( output_unit, '(a)' ) '%'
    write ( output_unit, '(a)' ) '0.900  0.200  0.100 setrgbcolor'
    write ( output_unit, '(a)' ) '%'
    write ( output_unit, '(a)' ) '%  Draw the elements.'
    write ( output_unit, '(a)' ) '%'

    do element = 1, element_num

      write ( output_unit, '(a)' ) 'newpath'

      do i = 1, element_order + 1

        e = i4_wrap ( i, 1, element_order )

        node = element_node(e,element)

        x_ps = int ( &
          ( ( x_max - node_xy(1,node)         ) &
            * real ( x_ps_min)   &
          + (         node_xy(1,node) - x_min ) &
            * real ( x_ps_max) ) &
          / ( x_max                   - x_min ) )

        y_ps = int ( &
          ( ( y_max - node_xy(2,node)         ) &
            * real ( y_ps_min)   &
          + (         node_xy(2,node) - y_min ) &
            * real ( y_ps_max) ) &
          / ( y_max                   - y_min ) )

        if ( i == 1 ) then
          write ( output_unit, '(i3,2x,i3,2x,a)' ) x_ps, y_ps, ' moveto'
        else
          write ( output_unit, '(i3,2x,i3,2x,a)' ) x_ps, y_ps, ' lineto'
        end if

      end do

      write ( output_unit, '(a)' ) 'stroke'

    end do

  end if
!
!  Label the elements.
!
  if ( 2 <= element_show ) then

    write ( output_unit, '(a)' ) '%'
    write ( output_unit, '(a)' ) '%  Label the elements:'
    write ( output_unit, '(a)' ) '%'
    write ( output_unit, '(a)' ) '%  Set the RGB color to darker red.'
    write ( output_unit, '(a)' ) '%'
    write ( output_unit, '(a)' ) '0.950  0.250  0.150 setrgbcolor'
    write ( output_unit, '(a)' ) '/Times-Roman findfont'
    write ( output_unit, '(a)' ) '0.20 inch scalefont'
    write ( output_unit, '(a)' ) 'setfont'
    write ( output_unit, '(a)' ) '%'

    do element = 1, element_num

      ave_x = 0.0D+00
      ave_y = 0.0D+00

      do i = 1, element_order

        node = element_node(i,element)

        ave_x = ave_x + node_xy(1,node)
        ave_y = ave_y + node_xy(2,node)

      end do

      ave_x = ave_x / real ( element_order)
      ave_y = ave_y / real ( element_order)

      x_ps = int ( &
        ( ( x_max - ave_x         ) * real ( x_ps_min)   &
        + (       + ave_x - x_min ) * real ( x_ps_max) ) &
        / ( x_max         - x_min ) )

      y_ps = int ( &
        ( ( y_max - ave_y         ) * real ( y_ps_min)   &
        + (         ave_y - y_min ) * real ( y_ps_max) ) &
        / ( y_max         - y_min ) )

      write ( string, '(i4)' ) element
      string = adjustl ( string )

      write ( output_unit, '(i4,2x,i4,a)' ) x_ps, y_ps, ' moveto (' &
        // trim ( string ) // ') show'

    end do

  end if

  write ( output_unit, '(a)' ) '%'
  write ( output_unit, '(a)' ) 'restore  showpage'
  write ( output_unit, '(a)' ) '%'
  write ( output_unit, '(a)' ) '%  End of page.'
  write ( output_unit, '(a)' ) '%'
  write ( output_unit, '(a)' ) '%%Trailer'
  write ( output_unit, '(a)' ) '%%EOF'
  close ( unit = output_unit )
end

subroutine reference_to_physical_q4 ( q4, n, rs, xy )

!*****************************************************************************80
!
!! REFERENCE_TO_PHYSICAL_Q4 maps Q4 reference points to physical points.
!
!  Discussion:
!
!    XY(R,S) = XY(0,0) * (1-R) * (1-S)
!            + XY(1,0) *    R  * (1-S)
!            + XY(1,1) *    R  *    S
!            + XY(0,1) * (1-R) *    S
!
!  Reference Element Q4:
!
!    |
!    1  4-----3
!    |  |     |
!    |  |     |
!    S  |     |
!    |  |     |
!    |  |     |
!    0  1-----2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision Q4(2,4), the coordinates of the vertices.
!    The vertices are assumed to be the images of the reference vertices
!    (0,0), (1,0), (1,1) and (0,1) respectively.
!
!    Input, integer N, the number of points to transform.
!
!    Input, double precision RS(2,N), (R,S) points in the reference element.
!
!    Output, double precision XY(2,N), (X,Y) points in the physical element.
!
  implicit none

  integer n

  double precision psi(4,n)
  double precision q4(2,4)
  double precision rs(2,n)
  double precision xy(2,n)

  psi(1,1:n) = ( 1.0D+00 - rs(1,1:n) ) * ( 1.0D+00 - rs(2,1:n) )
  psi(2,1:n) =             rs(1,1:n)   * ( 1.0D+00 - rs(2,1:n) )
  psi(3,1:n) =             rs(1,1:n)   *             rs(2,1:n)
  psi(4,1:n) = ( 1.0D+00 - rs(1,1:n) ) *             rs(2,1:n)

  xy(1:2,1:n) = matmul ( q4(1:2,1:4), psi(1:4,1:n) )
end

subroutine sample_q4_mesh ( node_num, node_xy, element_num, element_node, &
  sample_num, seed, sample_xy, sample_element )

!*****************************************************************************80
!
!! SAMPLE_Q4_MESH returns random points in a Q4 mesh.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 March 2009
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
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(4,ELEMENT_NUM), the nodes
!    that form the elements.
!
!    Input, integer SAMPLE_NUM, the number of points to sample.
!
!    Input/output, integer SEED, a seed for the random 
!     number generator.
!
!    Output, double precision SAMPLE_XY(2,SAMPLE_NUM), the sample points.
!
!    Output, integer SAMPLE_ELEMENT(SAMPLE_NUM), the elements from
!    which each point was drawn.
!
  implicit none

  integer element_num
  integer node_num
  integer sample_num

  double precision area
  double precision area_cum(0:element_num)
  double precision area_total
  integer element
  integer element_node(4,element_num)
  integer i1
  integer i2
  integer i3
  integer i4
  integer left
  double precision node_xy(2,node_num)
  double precision quad_xy(2,4)
  double precision r
  double precision r8_uniform_01
  integer right
  integer sample
  integer sample_element(sample_num)
  double precision sample_xy(2,sample_num)
  integer seed
!
!  Compute the areas of the quadrilaterals.
!
  area_cum(0) = 0.0D+00

  do element = 1, element_num

    i1 = element_node(1,element)
    i2 = element_node(2,element)
    i3 = element_node(3,element)
    i4 = element_node(4,element)

    quad_xy(1:2,1) = node_xy(1:2,i1)
    quad_xy(1:2,2) = node_xy(1:2,i2)
    quad_xy(1:2,3) = node_xy(1:2,i3)
    quad_xy(1:2,4) = node_xy(1:2,i4)

    call area_quad ( quad_xy, area )

    area_cum(element) = area_cum(element-1) + area

  end do

  area_total = area_cum(element_num)

  area_cum(0:element_num) = area_cum(0:element_num) / area_total
!
!  A random value R indicates the corresponding quadrilateral whose
!  cumulative relative area first includes the number R.
!
  do sample = 1, sample_num

    r = r8_uniform_01 ( seed )

    call r8vec_bracket ( element_num + 1, area_cum, r, left, right )

    element = right - 1

    i1 = element_node(1,element)
    i2 = element_node(2,element)
    i3 = element_node(3,element)
    i4 = element_node(4,element)

    quad_xy(1:2,1) = node_xy(1:2,i1)
    quad_xy(1:2,2) = node_xy(1:2,i2)
    quad_xy(1:2,3) = node_xy(1:2,i3)
    quad_xy(1:2,4) = node_xy(1:2,i4)
    
    call sample_quad ( quad_xy, 1, seed, sample_xy(1:2,sample) )

    sample_element(sample) = element

  end do
end

subroutine sample_quad ( quad_xy, n, seed, xy )

!*****************************************************************************80
!
!! SAMPLE_QUAD returns random points in a quadrilateral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision QUAD_XY(2,4), the coordinates of the nodes.
!
!    Input, integer N, the number of points to sample.
!
!    Input/output, integer SEED, a seed for the random 
!     number generator.
!
!    Output, double precision XY(2,N), the sample points.
!
  implicit none

  integer n

  double precision area1
  double precision area2
  double precision area_total
  integer i
  double precision quad_xy(2,4)
  double precision r
  double precision r8_uniform_01
  integer seed
  double precision t1(2,3)
  double precision t2(2,3)
  double precision xy(2,n)

  t1(1:2,1) = quad_xy(1:2,1)
  t1(1:2,2) = quad_xy(1:2,2)
  t1(1:2,3) = quad_xy(1:2,3)

  call triangle_area ( t1, area1 )

  t2(1:2,1) = quad_xy(1:2,3)
  t2(1:2,2) = quad_xy(1:2,4)
  t2(1:2,3) = quad_xy(1:2,1)

  call triangle_area ( t2, area2 )

  if ( area1 < 0.0D+00 .or. area2 < 0.0D+00 ) then

    t1(1:2,1) = quad_xy(1:2,2)
    t1(1:2,2) = quad_xy(1:2,3)
    t1(1:2,3) = quad_xy(1:2,4)

    call triangle_area ( t1, area1 )

    t2(1:2,1) = quad_xy(1:2,4)
    t2(1:2,2) = quad_xy(1:2,1)
    t2(1:2,3) = quad_xy(1:2,2)

    call triangle_area ( t2, area2 )

    if ( area1 < 0.0D+00 .or. area2 < 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SAMPLE_QUAD - Fatal error!'
      write ( *, '(a)' ) '  The quadrilateral nodes seem to be listed in'
      write ( *, '(a)' ) '  the wrong order, or the quadrilateral is'
      write ( *, '(a)' ) '  degenerate.'
      stop
    end if

  end if

  area_total = area1 + area2
!
!  Choose a triangle at random, weighted by the areas.
!  Then choose a point in that triangle.
!
  do i = 1, n

    r = r8_uniform_01 ( seed )

    if ( r * area_total < area1 ) then
      call triangle_sample ( t1, 1, seed, xy(1:2,i) )
    else
      call triangle_sample ( t2, 1, seed, xy(1:2,i) )
    end if

  end do
end

subroutine triangle_area ( t, area )

!*****************************************************************************80
!
!! TRIANGLE_AREA computes the area of a triangle.
!
!  Discussion:
!
!    If the triangle's vertices are given in counter clockwise order,
!    the area will be positive.  If the triangle's vertices are given
!    in clockwise order, the area will be negative!
!
!    An earlier version of this routine always returned the absolute
!    value of the computed area.  I am convinced now that that is
!    a less useful result!  For instance, by returning the signed 
!    area of a triangle, it is possible to easily compute the area 
!    of a nonconvex polygon as the sum of the (possibly negative) 
!    areas of triangles formed by node 1 and successive pairs of vertices.
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
!    Input, double precision T(2,3), the triangle vertices.
!
!    Output, double precision AREA, the area of the triangle.
!
  implicit none

  integer , parameter :: dim_num = 2

  double precision area
  double precision t(dim_num,3)

  area = 0.5D+00 * ( &
      t(1,1) * ( t(2,2) - t(2,3) ) &
    + t(1,2) * ( t(2,3) - t(2,1) ) &
    + t(1,3) * ( t(2,1) - t(2,2) ) )
end

subroutine triangle_sample ( t, n, seed, p )

!*****************************************************************************80
!
!! TRIANGLE_SAMPLE returns random points in a triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision T(2,3), the triangle vertices.
!
!    Input, integer N, the number of points to generate.
!
!    Input/output, integer SEED, a seed for the random number 
!    generator.
!
!    Output, double precision P(2,N), random points in the triangle.
!
  implicit none

  integer , parameter :: dim_num = 2
  integer n

  double precision alpha(n)
  integer dim
  double precision p(dim_num,n)
  double precision p12(dim_num,n)
  double precision p13(dim_num,n)
  integer seed
  double precision t(dim_num,3)
!
!  For comparison between F90, C++ and MATLAB codes, call R8VEC_UNIFORM_01.
!
  call r8vec_uniform_01 ( n, seed, alpha )
!
!  For faster execution, call RANDOM_NUMBER.
!
! call random_number ( harvest = alpha(1:n) )
!
!  Interpret R as a percentage of the triangle's area.
!
!  Imagine a line L, parallel to side 1, so that the area between
!  vertex 1 and line L is R percent of the full triangle's area.
!
!  The line L will intersect sides 2 and 3 at a fraction
!  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
!
  alpha(1:n) = sqrt ( alpha(1:n) )
!
!  Determine the coordinates of the points on sides 2 and 3 intersected
!  by line L.
!
  do dim = 1, dim_num

    p12(dim,1:n) = ( 1.0D+00 - alpha(1:n) ) * t(dim,1) &
                             + alpha(1:n)   * t(dim,2)

    p13(dim,1:n) = ( 1.0D+00 - alpha(1:n) ) * t(dim,1) &
                             + alpha(1:n)   * t(dim,3)

  end do
!
!  Now choose, uniformly at random, a point on the line L.
!
!  For comparison between F90, C++ and MATLAB codes, call R8VEC_UNIFORM_01.
!
  call r8vec_uniform_01 ( n, seed, alpha )
!
!  For faster execution, call RANDOM_NUMBER.
!
! call random_number ( harvest = alpha(1:n) )

  do dim = 1, dim_num

    p(dim,1:n) = ( 1.0D+00 - alpha(1:n) ) * p12(dim,1:n) &
                           + alpha(1:n)   * p13(dim,1:n)

  end do
end
