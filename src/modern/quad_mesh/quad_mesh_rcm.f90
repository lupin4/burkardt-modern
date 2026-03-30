!> quad_mesh_rcm — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

function adj_bandwidth ( node_num, adj_num, adj_row, adj )

!*****************************************************************************80
!
!! ADJ_BANDWIDTH computes the bandwidth of an adjacency matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Author:
!
!    Original FORTRAN77 version by Alan George, Joseph Liu.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about
!    row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Output, integer ADJ_BANDWIDTH, the bandwidth of the adjacency
!    matrix.
!
  implicit none

  integer adj_num
  integer node_num

  integer adj(adj_num)
  integer adj_bandwidth
  integer adj_row(node_num+1)
  integer band_hi
  integer band_lo
  integer col
  integer i
  integer j

  band_lo = 0
  band_hi = 0

  do i = 1, node_num

    do j = adj_row(i), adj_row(i+1) - 1
      col = adj(j)
      band_lo = max ( band_lo, i - col )
      band_hi = max ( band_hi, col - i )
    end do

  end do

  adj_bandwidth = band_lo + 1 + band_hi
end

function adj_perm_bandwidth ( node_num, adj_num, adj_row, adj, perm, perm_inv )

!*****************************************************************************80
!
!! ADJ_PERM_BANDWIDTH computes the bandwidth of a permuted adjacency matrix.
!
!  Discussion:
!
!    The matrix is defined by the adjacency information and a permutation.
!
!    The routine also computes the bandwidth and the size of the envelope.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about
!    row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input, integer PERM(NODE_NUM), PERM_INV(NODE_NUM), the
!    permutation and inverse permutation.
!
!    Output, integer ADJ_PERM_BANDWIDTH, the bandwidth of the
!    permuted adjacency matrix.
!
  implicit none

  integer adj_num
  integer node_num

  integer adj(adj_num)
  integer adj_perm_bandwidth
  integer adj_row(node_num+1)
  integer band_hi
  integer band_lo
  integer col
  integer i
  integer j
  integer perm(node_num)
  integer perm_inv(node_num)

  band_lo = 0
  band_hi = 0

  do i = 1, node_num

    do j = adj_row(perm(i)), adj_row(perm(i)+1) - 1
      col = perm_inv(adj(j))
      band_lo = max ( band_lo, i - col )
      band_hi = max ( band_hi, col - i )
    end do

  end do

  adj_perm_bandwidth = band_lo + 1 + band_hi
end

subroutine adj_set_q4_mesh ( node_num, element_num, element_node, &
  element_neighbor, adj_num, adj_row, adj )

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
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about
!    column J is stored in entries ADJ_ROW(J) through ADJ_ROW(J+1)-1 of ADJ.
!
!    Output, integer ADJ(ADJ_NUM), the adjacency information.
!
  implicit none

  integer adj_num
  integer node_num
  integer element_num
  integer , parameter :: element_order = 4

  integer adj(adj_num)
  integer adj_row(node_num+1)
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
  adj_copy(1:node_num) = adj_row(1:node_num)
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
    k1 = adj_row(node)
    k2 = adj_row(node+1)-1
    number = k2 + 1 - k1
    call i4vec_sort_heap_a ( number, adj(k1:k2) )
  end do
end

subroutine adj_size_q4_mesh ( node_num, element_num, element_node, &
  element_neighbor, adj_num, adj_row )

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
!    Output, integer ADJ_ROW(NODE_NUM+1), Information about
!    column J is stored in entries ADJ_ROW(J) through ADJ_ROW(J+1)-1 of ADJ.
!
  implicit none

  integer element_num
  integer , parameter :: element_order = 4
  integer node_num

  integer adj_row(node_num+1)
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
  adj_row(1:node_num) = 1
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
    adj_row(n1) = adj_row(n1) + 1
    adj_row(n3) = adj_row(n3) + 1
!
!  Add edge (2,4).
!
    adj_row(n2) = adj_row(n2) + 1
    adj_row(n4) = adj_row(n4) + 1
!
!  Add edge (1,2) if this is the first occurrence,
!  that is, if the edge (1,2) is on a boundary (ELEMENT2 <= 0)
!  or if this element is the first of the pair in which the edge
!  occurs (ELEMENT < ELEMENT2).
!
    element2 = element_neighbor(1,element)

    if ( element2 < 0 .or. element < element2 ) then
      adj_row(n1) = adj_row(n1) + 1
      adj_row(n2) = adj_row(n2) + 1
    end if
!
!  Add edge (2,3).
!
    element2 = element_neighbor(2,element)

    if ( element2 < 0 .or. element < element2 ) then
      adj_row(n2) = adj_row(n2) + 1
      adj_row(n3) = adj_row(n3) + 1
    end if
!
!  Add edge (3,4).
!
    element2 = element_neighbor(3,element)

    if ( element2 < 0 .or. element < element2 ) then
      adj_row(n3) = adj_row(n3) + 1
      adj_row(n4) = adj_row(n4) + 1
    end if
!
!  Add edge (4,1).
!
    element2 = element_neighbor(4,element)

    if ( element2 < 0 .or. element < element2 ) then
      adj_row(n4) = adj_row(n4) + 1
      adj_row(n1) = adj_row(n1) + 1
    end if

  end do
!
!  We used ADJ_ROW to count the number of entries in each column.
!  Convert it to pointers into the ADJ array.
!
  do node = node_num + 1, 2, -1
    adj_row(node) = adj_row(node-1)
  end do

  adj_row(1) = 1
  do node = 2, node_num + 1
    adj_row(node) = adj_row(node) + adj_row(node-1)
  end do
!
!  Finally, record the total number of adjacencies.
!
  adj_num = adj_row(node_num+1) - 1
end

subroutine degree ( root, adj_num, adj_row, adj, mask, deg, iccsze, ls, &
  node_num )

!*****************************************************************************80
!
!! DEGREE computes the degrees of the nodes in the connected component.
!
!  Discussion:
!
!    The connected component is specified by MASK and ROOT.
!    Nodes for which MASK is zero are ignored.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2003
!
!  Author:
!
!    Original FORTRAN77 version by Alan George, Joseph Liu.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ROOT, the node that defines the connected
!    component.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about
!    row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input, integer MASK(NODE_NUM), is nonzero for those nodes
!    which are to be considered.
!
!    Output, integer DEG(NODE_NUM), contains, for each  node in
!    the connected component, its degree.
!
!    Output, integer ICCSIZE, the number of nodes in the
!    connected component.
!
!    Output, integer LS(NODE_NUM), stores in entries 1 through
!    ICCSIZE the nodes in the connected component, starting with ROOT, and
!    proceeding by levels.
!
!    Input, integer NODE_NUM, the number of nodes.
!
  implicit none

  integer adj_num
  integer node_num

  integer adj(adj_num)
  integer adj_row(node_num+1)
  integer deg(node_num)
  integer i
  integer iccsze
  integer ideg
  integer j
  integer jstop
  integer jstrt
  integer lbegin
  integer ls(node_num)
  integer lvlend
  integer lvsize
  integer mask(node_num)
  integer nbr
  integer node
  integer root
!
!  The sign of ADJ_ROW(I) is used to indicate if node I has been considered.
!
  ls(1) = root
  adj_row(root) = -adj_row(root)
  lvlend = 0
  iccsze = 1
!
!  LBEGIN is the pointer to the beginning of the current level, and
!  LVLEND points to the end of this level.
!
  do

    lbegin = lvlend + 1
    lvlend = iccsze
!
!  Find the degrees of nodes in the current level,
!  and at the same time, generate the next level.
!
    do i = lbegin, lvlend

      node = ls(i)
      jstrt = -adj_row(node)
      jstop = abs ( adj_row(node+1) ) - 1
      ideg = 0

      do j = jstrt, jstop

        nbr = adj(j)

        if ( mask(nbr) /= 0 ) then

          ideg = ideg + 1

          if ( 0 <= adj_row(nbr) ) then
            adj_row(nbr) = -adj_row(nbr)
            iccsze = iccsze + 1
            ls(iccsze) = nbr
          end if

        end if

      end do

      deg(node) = ideg

    end do
!
!  Compute the current level width.
!
    lvsize = iccsze - lvlend
!
!  If the current level width is nonzero, generate another level.
!
    if ( lvsize == 0 ) then
      exit
    end if

  end do
!
!  Reset ADJ_ROW to its correct sign and return.
!
  do i = 1, iccsze
    node = ls(i)
    adj_row(node) = -adj_row(node)
  end do
end

subroutine genrcm ( node_num, adj_num, adj_row, adj, perm )

!*****************************************************************************80
!
!! GENRCM finds the reverse Cuthill-Mckee ordering for a general graph.
!
!  Discussion:
!
!    For each connected component in the graph, the routine obtains
!    an ordering by calling RCM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 January 2003
!
!  Author:
!
!    Original FORTRAN77 version by Alan George, Joseph Liu.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about
!    row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Output, integer PERM(NODE_NUM), the RCM ordering.
!
!  Local Parameters:
!
!    Local, integer LEVEL_ROW(NODE_NUM+1), the index vector for a level
!    structure.  The level structure is stored in the currently unused
!    spaces in the permutation vector PERM.
!
!    Local, integer MASK(NODE_NUM), marks variables that have been numbered.
!
  implicit none

  integer adj_num
  integer node_num

  integer adj(adj_num)
  integer adj_row(node_num+1)
  integer i
  integer iccsze
  integer mask(node_num)
  integer level_num
  integer level_row(node_num+1)
  integer num
  integer perm(node_num)
  integer root

  mask(1:node_num) = 1

  num = 1

  do i = 1, node_num
!
!  For each masked connected component...
!
    if ( mask(i) /= 0 ) then

      root = i
!
!  Find a pseudo-peripheral node ROOT.  The level structure found by
!  ROOT_FIND is stored starting at PERM(NUM).
!
      call root_find ( root, adj_num, adj_row, adj, mask, level_num, &
        level_row, perm(num), node_num )

      write ( *, * ) '  ROOT = ',root
!
!  RCM orders the component using ROOT as the starting node.
!
      call rcm ( root, adj_num, adj_row, adj, mask, perm(num), iccsze, &
        node_num )

      num = num + iccsze
!
!  We can stop once every node is in one of the connected components.
!
      if ( node_num < num ) then
      end if

    end if

  end do
end

subroutine level_set ( root, adj_num, adj_row, adj, mask, level_num, &
  level_row, level, node_num )

!*****************************************************************************80
!
!! LEVEL_SET generates the connected level structure rooted at a given node.
!
!  Discussion:
!
!    Only nodes for which MASK is nonzero will be considered.
!
!    The root node chosen by the user is assigned level 1, and masked.
!    All (unmasked) nodes reachable from a node in level 1 are
!    assigned level 2 and masked.  The process continues until there
!    are no unmasked nodes adjacent to any node in the current level.
!    The number of levels may vary between 2 and NODE_NUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 October 2003
!
!  Author:
!
!    Original FORTRAN77 version by Alan George, Joseph Liu.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ROOT, the node at which the level structure
!    is to be rooted.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about
!    row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input/output, integer MASK(NODE_NUM).  On input, only nodes
!    with nonzero MASK are to be processed.  On output, those nodes which were
!    included in the level set have MASK set to 1.
!
!    Output, integer LEVEL_NUM, the number of levels in the level
!    structure.  ROOT is in level 1.  The neighbors of ROOT
!    are in level 2, and so on.
!
!    Output, integer LEVEL_ROW(NODE_NUM+1), LEVEL(NODE_NUM),
!    the rooted level structure.
!
!    Input, integer NODE_NUM, the number of nodes.
!
  implicit none

  integer adj_num
  integer node_num

  integer adj(adj_num)
  integer adj_row(node_num+1)
  integer i
  integer iccsze
  integer j
  integer jstop
  integer jstrt
  integer lbegin
  integer level_num
  integer level_row(node_num+1)
  integer level(node_num)
  integer lvlend
  integer lvsize
  integer mask(node_num)
  integer nbr
  integer node
  integer root

  mask(root) = 0
  level(1) = root
  level_num = 0
  lvlend = 0
  iccsze = 1
!
!  LBEGIN is the pointer to the beginning of the current level, and
!  LVLEND points to the end of this level.
!
  do

    lbegin = lvlend + 1
    lvlend = iccsze
    level_num = level_num + 1
    level_row(level_num) = lbegin
!
!  Generate the next level by finding all the masked neighbors of nodes
!  in the current level.
!
    do i = lbegin, lvlend

      node = level(i)
      jstrt = adj_row(node)
      jstop = adj_row(node+1) - 1

      do j = jstrt, jstop

        nbr = adj(j)

        if ( mask(nbr) /= 0 ) then
          iccsze = iccsze + 1
          level(iccsze) = nbr
          mask(nbr) = 0
        end if

      end do

    end do
!
!  Compute the current level width (the number of nodes encountered.)
!  If it is positive, generate the next level.
!
    lvsize = iccsze - lvlend

    if ( lvsize <= 0 ) then
      exit
    end if

  end do

  level_row(level_num+1) = lvlend + 1
!
!  Reset MASK to 1 for the nodes in the level structure.
!
  mask(level(1:iccsze)) = 1
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

subroutine r8col_permute ( m, n, p, a )

!*****************************************************************************80
!
!! R8COL_PERMUTE permutes an R8COL in place.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    The same logic can be used to permute an array of objects of any
!    arithmetic type, or an array of objects of any complexity.  The only
!    temporary storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!
!    Input:
!
!      M = 2
!      N = 5
!      P = (   2,    4,    5,    1,    3 )
!      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
!          (11.0, 22.0, 33.0, 44.0, 55.0 )
!
!    Output:
!
!      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
!             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the dimension of objects.
!
!    Input, integer N, the number of objects.
!
!    Input, integer P(N), the permutation.  P(I) = J means
!    that the I-th element of the output array should be the J-th
!    element of the input array.
!
!    Input/output, double precision A(M,N), the array to be permuted.
!
  implicit none

  integer m
  integer n

  double precision a(m,n)
  double precision a_temp(m)
  integer , parameter :: base = 1
  integer ierror
  integer iget
  integer iput
  integer istart
  integer p(n)

  call perm_check ( n, p, base, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8COL_PERMUTE - Fatal error!'
    write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
    stop
  end if
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = - p(istart)
      cycle

    else

      a_temp(1:m) = a(1:m,istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = - p(iput)

        if ( iget < 1 .or. n < iget ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8COL_PERMUTE - Fatal error!'
          write ( *, '(a)' ) '  A permutation index is out of range.'
          write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
          stop
        end if

        if ( iget == istart ) then
          a(1:m,iput) = a_temp(1:m)
          exit
        end if

        a(1:m,iput) = a(1:m,iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = - p(1:n)
end

subroutine rcm ( root, adj_num, adj_row, adj, mask, perm, iccsze, node_num )

!*****************************************************************************80
!
!! RCM renumbers a connected component by the reverse Cuthill McKee algorithm.
!
!  Discussion:
!
!    The connected component is specified by a node ROOT and a mask.
!    The numbering starts at the root node.
!
!    An outline of the algorithm is as follows:
!
!    X(1) = ROOT.
!
!    for ( I = 1 to N-1)
!      Find all unlabeled neighbors of X(I),
!      assign them the next available labels, in order of increasing degree.
!
!    When done, reverse the ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2007
!
!  Author:
!
!    Original FORTRAN77 version by Alan George, Joseph Liu.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ROOT, the node that defines the connected
!    component.  It is used as the starting point for the RCM ordering.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about
!    row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input/output, integer MASK(NODE_NUM), a mask for the nodes.
!    Only those nodes with nonzero input mask values are considered by the
!    routine.  The nodes numbered by RCM will have their mask values
!    set to zero.
!
!    Output, integer PERM(NODE_NUM), the RCM ordering.
!
!    Output, integer ICCSZE, the size of the connected component
!    that has been numbered.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!  Local Parameters:
!
!    Workspace, integer DEG(NODE_NUM), a temporary vector used to hold
!    the degree of the nodes in the section graph specified by mask and root.
!
  implicit none

  integer adj_num
  integer node_num

  integer adj(adj_num)
  integer adj_row(node_num+1)
  integer deg(node_num)
  integer fnbr
  integer i
  integer iccsze
  integer j
  integer jstop
  integer jstrt
  integer k
  integer l
  integer lbegin
  integer lnbr
  integer lperm
  integer lvlend
  integer mask(node_num)
  integer nbr
  integer node
  integer perm(node_num)
  integer root
!
!  Find the degrees of the nodes in the component specified by MASK and ROOT.
!
  call degree ( root, adj_num, adj_row, adj, mask, deg, iccsze, perm, node_num )

  mask(root) = 0

  if ( iccsze <= 1 ) then
  end if

  lvlend = 0
  lnbr = 1
!
!  LBEGIN and LVLEND point to the beginning and
!  the end of the current level respectively.
!
  do while ( lvlend < lnbr )

    lbegin = lvlend + 1
    lvlend = lnbr

    do i = lbegin, lvlend
!
!  For each node in the current level...
!
      node = perm(i)
      jstrt = adj_row(node)
      jstop = adj_row(node+1) - 1
!
!  Find the unnumbered neighbors of NODE.
!
!  FNBR and LNBR point to the first and last neighbors
!  of the current node in PERM.
!
      fnbr = lnbr + 1

      do j = jstrt, jstop

        nbr = adj(j)

        if ( mask(nbr) /= 0 ) then
          lnbr = lnbr + 1
          mask(nbr) = 0
          perm(lnbr) = nbr
        end if

      end do
!
!  If no neighbors, skip to next node in this level.
!
      if ( lnbr <= fnbr ) then
        cycle
      end if
!
!  Sort the neighbors of NODE in increasing order by degree.
!  Linear insertion is used.
!
      k = fnbr

      do while ( k < lnbr )

        l = k
        k = k + 1
        nbr = perm(k)

        do while ( fnbr < l )

          lperm = perm(l)

          if ( deg(lperm) <= deg(nbr) ) then
            exit
          end if

          perm(l+1) = lperm
          l = l - 1

        end do

        perm(l+1) = nbr

      end do

    end do

  end do
!
!  We now have the Cuthill-McKee ordering.  Reverse it.
!
  call i4vec_reverse ( iccsze, perm )
end

subroutine root_find ( root, adj_num, adj_row, adj, mask, level_num, &
  level_row, level, node_num )

!*****************************************************************************80
!
!! ROOT_FIND finds a pseudo-peripheral node.
!
!  Discussion:
!
!    The diameter of a graph is the maximum distance (number of edges)
!    between any two nodes of the graph.
!
!    The eccentricity of a node is the maximum distance between that
!    node and any other node of the graph.
!
!    A peripheral node is a node whose eccentricity equals the
!    diameter of the graph.
!
!    A pseudo-peripheral node is an approximation to a peripheral node;
!    it may be a peripheral node, but all we know is that we tried our
!    best.
!
!    The routine is given a graph, and seeks pseudo-peripheral nodes,
!    using a modified version of the scheme of Gibbs, Poole and
!    Stockmeyer.  It determines such a node for the section subgraph
!    specified by MASK and ROOT.
!
!    The routine also determines the level structure associated with
!    the given pseudo-peripheral node; that is, how far each node
!    is from the pseudo-peripheral node.  The level structure is
!    returned as a list of nodes LS, and pointers to the beginning
!    of the list of nodes that are at a distance of 0, 1, 2, ...,
!    NODE_NUM-1 from the pseudo-peripheral node.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 October 2003
!
!  Author:
!
!    Original FORTRAN77 version by Alan George, Joseph Liu.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!    Norman Gibbs, William Poole, Paul Stockmeyer,
!    An Algorithm for Reducing the Bandwidth and Profile of a Sparse Matrix,
!    SIAM Journal on Numerical Analysis,
!    Volume 13, pages 236-250, 1976.
!
!    Norman Gibbs,
!    Algorithm 509: A Hybrid Profile Reduction Algorithm,
!    ACM Transactions on Mathematical Software,
!    Volume 2, pages 378-387, 1976.
!
!  Parameters:
!
!    Input/output, integer ROOT.  On input, ROOT is a node in the
!    the component of the graph for which a pseudo-peripheral node is
!    sought.  On output, ROOT is the pseudo-peripheral node obtained.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about
!    row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input, integer MASK(NODE_NUM), specifies a section subgraph.
!    Nodes for which MASK is zero are ignored by FNROOT.
!
!    Output, integer LEVEL_NUM, is the number of levels in the
!    level structure rooted at the node ROOT.
!
!    Output, integer LEVEL_ROW(NODE_NUM+1), LEVEL(NODE_NUM), the
!    level structure array pair containing the level structure found.
!
!    Input, integer NODE_NUM, the number of nodes.
!
  implicit none

  integer adj_num
  integer node_num

  integer adj(adj_num)
  integer adj_row(node_num+1)
  integer iccsze
  integer j
  integer jstrt
  integer k
  integer kstop
  integer kstrt
  integer level(node_num)
  integer level_num
  integer level_num2
  integer level_row(node_num+1)
  integer mask(node_num)
  integer mindeg
  integer nabor
  integer ndeg
  integer node
  integer root
!
!  Determine the level structure rooted at ROOT.
!
  call level_set ( root, adj_num, adj_row, adj, mask, level_num, &
    level_row, level, node_num )
!
!  Count the number of nodes in this level structure.
!
  iccsze = level_row(level_num+1) - 1
!
!  Extreme case:
!    A complete graph has a level set of only a single level.
!    Every node is equally good (or bad).
!
  if ( level_num == 1 ) then
  end if
!
!  Extreme case:
!    A "line graph" 0--0--0--0--0 has every node in its only level.
!    By chance, we've stumbled on the ideal root.
!
  if ( level_num == iccsze ) then
  end if
!
!  Pick any node from the last level that has minimum degree
!  as the starting point to generate a new level set.
!
  do

    mindeg = iccsze

    jstrt = level_row(level_num)
    root = level(jstrt)

    if ( jstrt < iccsze ) then

      do j = jstrt, iccsze

        node = level(j)
        ndeg = 0
        kstrt = adj_row(node)
        kstop = adj_row(node+1) - 1

        do k = kstrt, kstop
          nabor = adj(k)
          if ( 0 < mask(nabor) ) then
            ndeg = ndeg + 1
          end if
        end do

        if ( ndeg < mindeg ) then
          root = node
          mindeg = ndeg
        end if

      end do

    end if
!
!  Generate the rooted level structure associated with this node.
!
    call level_set ( root, adj_num, adj_row, adj, mask, level_num2, &
      level_row, level, node_num )
!
!  If the number of levels did not increase, accept the new ROOT.
!
    if ( level_num2 <= level_num ) then
      exit
    end if

    level_num = level_num2
!
!  In the unlikely case that ROOT is one endpoint of a line graph,
!  we can exit now.
!
    if ( iccsze <= level_num ) then
      exit
    end if

  end do
end
