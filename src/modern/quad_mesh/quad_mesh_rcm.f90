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

subroutine bandwidth ( element_order, element_num, element_node, ml, mu, m )

!*****************************************************************************80
!
!! BANDWIDTH determines the bandwidth associated with a finite element mesh.
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
!    Input, integer ELEMENT_ORDER, the order of the elements.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
!    ELEMENT_NODE(I,J) is the global index of local node I in element J.
!
!    Output, integer ML, MU, the lower and upper bandwidths
!    of the matrix.
!
!    Output, integer M, the bandwidth of the matrix.
!
  implicit none

  integer element_num
  integer element_order

  integer element
  integer element_node(element_order,element_num)
  integer global_i
  integer global_j
  integer local_i
  integer local_j
  integer m
  integer ml
  integer mu

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

subroutine i4col_compare ( m, n, a, i, j, isgn )

!*****************************************************************************80
!
!! I4COL_COMPARE compares columns I and J of an I4COL.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      ISGN = -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(M,N), an array of N columns of vectors
!    of length M.
!
!    Input, integer I, J, the columns to be compared.
!    I and J must be between 1 and N.
!
!    Output, integer ISGN, the results of the comparison:
!    -1, column I < column J,
!     0, column I = column J,
!    +1, column J < column I.
!
  implicit none

  integer m
  integer n

  integer a(m,n)
  integer i
  integer isgn
  integer j
  integer k
!
!  Check.
!
  if ( i < 1 .or. n < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index I is out of bounds.'
    stop
  end if

  if ( j < 1 .or. n < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index J is out of bounds.'
    stop
  end if

  isgn = 0

  if ( i == j ) then
  end if

  k = 1

  do while ( k <= m )

    if ( a(k,i) < a(k,j) ) then
      isgn = -1
      return
    else if ( a(k,j) < a(k,i) ) then
      isgn = +1
    end if

    k = k + 1

  end do
end

subroutine i4col_sort_a ( m, n, a )

!*****************************************************************************80
!
!! I4COL_SORT_A ascending sorts an I4COL.
!
!  Discussion:
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows of A, and the length of
!    a vector of data.
!
!    Input, integer N, the number of columns of A.
!
!    Input/output, integer A(M,N).
!    On input, the array of N columns of M-vectors.
!    On output, the columns of A have been sorted in ascending
!    lexicographic order.
!
  implicit none

  integer m
  integer n

  integer a(m,n)
  integer i
  integer indx
  integer isgn
  integer j

  if ( m <= 0 ) then
  end if

  if ( n <= 1 ) then
  end if
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

      call i4col_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4col_compare ( m, n, a, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do
end

subroutine i4col_swap ( m, n, a, i, j )

!*****************************************************************************80
!
!! I4COL_SWAP swaps columns I and J of an I4COL.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      A = (
!        1  4  3  2
!        5  8  7  6
!        9 12 11 10 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns
!    in the array.
!
!    Input/output, integer A(M,N), an array of N columns
!    of length M.
!
!    Input, integer I, J, the columns to be swapped.
!
  implicit none

  integer m
  integer n

  integer a(m,n)
  integer col(m)
  integer i
  integer j

  if ( i < 1 .or. n < i .or. j < 1 .or. n < j ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_SWAP - Fatal error!'
    write ( *, '(a)' ) '  I or J is out of bounds.'
    write ( *, '(a,i8)' ) '  I =    ', i
    write ( *, '(a,i8)' ) '  J =    ', j
    write ( *, '(a,i8)' ) '  N =    ', n
    stop

  end if

  if ( i == j ) then
  end if

  col(1:m) = a(1:m,i)
  a(1:m,i) = a(1:m,j)
  a(1:m,j) = col(1:m)
end

subroutine i4vec_heap_d ( n, a )

!*****************************************************************************80
!
!! I4VEC_HEAP_D reorders an I4VEC into an descending heap.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    A descending heap is an array A with the property that, for every index J,
!    A(J) >= A(2*J) and A(J) >= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer N, the size of the input array.
!
!    Input/output, integer A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  implicit none

  integer n

  integer a(n)
  integer i
  integer ifree
  integer key
  integer m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n/2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = a(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( n < m ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the larger of the two values,
!  and update M if necessary.
!
        if ( a(m) < a(m+1) ) then
          m = m + 1
        end if

      end if
!
!  If the large descendant is larger than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( a(m) <= key ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot IFREE.
!
    a(ifree) = key

  end do
end

subroutine i4vec_reverse ( n, a )

!*****************************************************************************80
!
!! I4VEC_REVERSE reverses the elements of an I4VEC.
!
!  Example:
!
!    Input:
!
!      N = 5,
!      A = ( 11, 12, 13, 14, 15 ).
!
!    Output:
!
!      A = ( 15, 14, 13, 12, 11 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, integer A(N), the array to be reversed.
!
  implicit none

  integer n

  integer a(n)

  a(1:n) = a(n:1:-1)
end

subroutine i4vec_sort_heap_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_HEAP_A ascending sorts an I4VEC using heap sort.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, integer A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer n

  integer a(n)
  integer n1
  integer t

  if ( n <= 1 ) then
  end if
!
!  1: Put A into descending heap form.
!
  call i4vec_heap_d ( n, a )
!
!  2: Sort A.
!
!  The largest object in the heap is in A(1).
!  Move it to position A(N).
!
  t    = a(1)
  a(1) = a(n)
  a(n) = t
!
!  Consider the diminished heap of size N1.
!
  do n1 = n - 1, 2, -1
!
!  Restore the heap structure of A(1) through A(N1).
!
    call i4vec_heap_d ( n1, a )
!
!  Take the largest object from A(1) and move it to A(N1).
!
    t     = a(1)
    a(1)  = a(n1)
    a(n1) = t

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

subroutine perm_check ( n, p, base, ierror )

!*****************************************************************************80
!
!! PERM_CHECK checks that a vector represents a permutation.
!
!  Discussion:
!
!    The routine verifies that each of the integers from BASE to
!    to BASE+N-1 occurs among the N entries of the permutation.
!
!    Set the input quantity BASE to 0, if P is a 0-based permutation,
!    or to 1 if P is a 1-based permutation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries.
!
!    Input, integer P(N), the array to check.
!
!    Input, integer BASE, the index base.
!
!    Output, integer IERROR, error flag.
!    0, the array represents a permutation.
!    nonzero, the array does not represent a permutation.  The smallest
!    missing value is equal to IERROR.
!
  implicit none

  integer n

  integer base
  integer find
  integer ierror
  integer p(n)
  integer seek

  ierror = 0

  do seek = base, base + n - 1

    ierror = 1

    do find = 1, n
      if ( p(find) == seek ) then
        ierror = 0
        exit
      end if
    end do

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PERM_CHECK - Fatal error!'
      write ( *, '(a)' ) '  The input array does not represent'
      write ( *, '(a)' ) '  a proper permutation.'
      stop
    end if

  end do
end

subroutine perm_inverse3 ( n, perm, perm_inv )

!*****************************************************************************80
!
!! PERM_INVERSE3 produces the inverse of a given permutation.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items permuted.
!
!    Input, integer PERM(N), a permutation.
!
!    Output, integer PERM_INV(N), the inverse permutation.
!
  implicit none

  integer n

  integer i
  integer perm(n)
  integer perm_inv(n)

  do i = 1, n
    perm_inv(perm(i)) = i
  end do
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

subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer N, the number of items to be sorted.
!
!    Input/output, integer INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ISGN, results of comparison of elements
!    I and J. (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer i
  integer , save :: i_save = 0
  integer indx
  integer isgn
  integer j
  integer , save :: j_save = 0
  integer , save :: k = 0
  integer , save :: k1 = 0
  integer n
  integer , save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
  end if
end
