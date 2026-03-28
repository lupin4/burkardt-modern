!> tet_mesh_boundary — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

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
!    12 June 2005
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
  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a,i8,a)' ) '  Column index I = ', i, ' is less than 1.'
    stop
  end if

  if ( n < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a,i8,a)' ) '  N = ', n, ' is less than column index I = ', i
    stop
  end if

  if ( j < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a,i8,a)' ) '  Column index J = ', j, ' is less than 1.'
    stop
  end if

  if ( n < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a,i8,a)' ) '  N = ', n, ' is less than column index J = ', j
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
!! I4COL_SORT_A ascending sorts an I4COL of columns.
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

subroutine i4col_swap ( m, n, a, j1, j2 )

!*****************************************************************************80
!
!! I4COL_SWAP swaps columns J1 and J2 of a integer array of column data.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, J1 = 2, J2 = 4
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
!    Input/output, integer A(M,N), an array of N columns of 
!    length M.
!
!    Input, integer J1, J2, the columns to be swapped.
!
  implicit none

  integer m
  integer n

  integer a(m,n)
  integer col(m)
  integer j1
  integer j2

  if ( j1 < 1 .or. n < j1 .or. j2 < 1 .or. n < j2 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_SWAP - Fatal error!'
    write ( *, '(a)' ) '  J1 or J2 is out of bounds.'
    write ( *, '(a,i8)' ) '  J1 =    ', j1
    write ( *, '(a,i8)' ) '  J2 =    ', j2
    write ( *, '(a,i8)' ) '  N =     ', n
    stop

  end if

  if ( j1 == j2 ) then
  end if

  col(1:m)  = a(1:m,j1)
  a(1:m,j1) = a(1:m,j2)
  a(1:m,j2) = col(1:m)
end

subroutine i4i4i4_sort_a ( i1, i2, i3, j1, j2, j3 )

!*****************************************************************************80
!
!! I4I4I4_SORT_A ascending sorts a triple of I4's.
!
!  Discussion:
!
!    The program allows the reasonable call:
!
!      call i4i4i4_sort_a ( i1, i2, i3, i1, i2, i3 )
!
!    and this will return the reasonable result.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I1, I2, I3, the values to sort.
!
!    Output, integer J1, J2, J3, the sorted values.
!
  implicit none

  integer i1
  integer i2
  integer i3
  integer j1
  integer j2
  integer j3
  integer k1
  integer k2
  integer k3
!
!  Copy arguments, so that the user can make "reasonable" calls like:
!
!    call i4i4i4_sort_a ( i1, i2, i3, i1, i2, i3 )
!
  k1 = i1
  k2 = i2
  k3 = i3

  j1 = min ( min ( k1, k2 ), min ( k2, k3 ) )
  j2 = min ( max ( k1, k2 ), &
       min ( max ( k2, k3 ), max ( k3, k1 ) ) )
  j3 = max ( max ( k1, k2 ), max ( k2, k3 ) )
end

subroutine i4vec_cum ( n, a, a_cum )

!*****************************************************************************80
!
!! I4VEC_CUM computes the cumulutive sum of the entries of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Example:
!
!    Input:
!
!      A = (/ 1, 2, 3, 4 /)
!
!    Output:
!
!      A_CUM = (/ 1, 3, 6, 10 /)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, integer A(N), the vector to be summed.
!
!    Output, integer A_CUM(N), the cumulative sum of the
!    entries of A.
!
  implicit none

  integer n

  integer a(n)
  integer a_cum(1:n)
  integer i

  a_cum(1) = a(1)

  do i = 2, n
    a_cum(i) = a_cum(i-1) + a(i)
  end do
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
!    29 September 2009
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
  implicit none

  integer element_num
  integer element_order

  integer element
  integer element_node(element_order,element_num)
  integer node
  integer node_max
  integer node_min
  integer node_num
  integer order

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
    write ( *, '(a)' )'MESH_BASE_ONE - Warning!'
    write ( *, '(a)' )' The element indexing is not of a recognized type.'
  end if
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
!    A Nijenhuis and H Wilf,
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
!    Input, integer ISGN, results of comparison of elements I 
!    and J.
!    (Used only when the previous call returned INDX less than 0).
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

subroutine tet_mesh_boundary_count ( element_order, element_num, element_node, &
  node_num, boundary_node_num, boundary_element_num, boundary_node_mask )

!*****************************************************************************80
!
!! TET_MESH_BOUNDARY_COUNT counts boundary faces and nodes in a tet mesh.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 December 2010
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
!    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), the
!    nodes that make up each element.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Output, integer BOUNDARY_NODE_NUM, the number of 
!    boundary nodes.
!
!    Output, integer BOUNDARY_ELEMENT_NUM, the number of 
!    boundary faces.
!
!    Output, integer BOUNDARY_NODE_MASK(NODE_NUM), is 0 for
!    interior nodes, 1 for boundary nodes.
!
  implicit none

  integer element_num
  integer element_order
  integer node_num

  integer a
  integer b
  integer boundary_element_num
  integer boundary_node_mask(node_num)
  integer boundary_node_num
  integer c
  integer element
  integer element_node(element_order,element_num)
  integer element1
  integer element2
  integer f
  integer face
  integer face1
  integer face2
  integer , allocatable, dimension ( :, : ) :: faces
  integer i
  integer j
  integer k
  integer l

  if ( element_order /= 4 .and. element_order /= 10 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TET_MESH_BOUNDARY_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Unexpected element order.'
    stop
  end if

  allocate ( faces(1:5,1:4*element_num) )
!
!  Step 1.
!  From the list of nodes forming tetrahedron T, of the form: 
!
!    (I,J,K,L)
!
!  or
!
!    (I,J,K,L,I+J,I+K,I+L,J+K,J+L,K+L),
!    (1,2,3,4, 5,  6,  7,  8,  9, 10 ),
!
!  construct the four face relations:
!
!    F = 1: (J,K,L,F,T)
!    F = 2: (I,K,L,F,T)
!    F = 3: (I,J,L,F,T)
!    F = 4: (I,J,K,F,T)
!
!  If T is actually order 10, we can retrieve the indices of the midside
!  nodes from the values of F and T.  In that case, the 4 faces are:
!
!    F = 1: 2, 3, 4, 8, 10, 9
!    F = 2: 1, 3, 4, 6, 10, 7
!    F = 3: 1, 2, 4, 5,  9, 7
!    F = 4: 1, 2, 3, 5,  8, 6
!
!  In order to make matching easier, we reorder each triple of nodes
!  into ascending order.
!
  do element = 1, element_num

    i = element_node(1,element)
    j = element_node(2,element)
    k = element_node(3,element)
    l = element_node(4,element)

    call i4i4i4_sort_a ( j, k, l, a, b, c )

    faces(1:5,4*(element-1)+1) = (/ a, b, c, 1, element /)

    call i4i4i4_sort_a ( i, k, l, a, b, c )

    faces(1:5,4*(element-1)+2) = (/ a, b, c, 2, element /)

    call i4i4i4_sort_a ( i, j, l, a, b, c )

    faces(1:5,4*(element-1)+3) = (/ a, b, c, 3, element /)

    call i4i4i4_sort_a ( i, j, k, a, b, c )

    faces(1:5,4*(element-1)+4) = (/ a, b, c, 4, element /)

  end do
!
!  Step 2. Perform an ascending dictionary sort on the neighbor relations.
!  We only intend to sort on rows 1:3; the routine we call here
!  sorts on rows 1 through 5 but that won't hurt us.
!
!  What we need is to find cases where two tetrahedrons share a face.
!  By sorting the columns of the FACES array, we will put shared faces
!  next to each other.
!
  call i4col_sort_a ( 5, 4*element_num, faces )
!
!  Step 3. Neighboring faces show up as consecutive columns with
!  identical first three entries.  Count columns which don't have
!  a following column that matches the first three entries.
!
  boundary_element_num = 0
  face = 1

  boundary_node_mask(1:node_num) = 0

  do while ( face <= 4 * element_num )

    if ( face < 4 * element_num ) then

      if ( all ( faces(1:3,face) == faces(1:3,face+1) ) ) then
        face = face + 2
        cycle
      end if

    end if

    boundary_element_num = boundary_element_num + 1
!
!  The vertices of the triangle are boundary nodes.
!
    boundary_node_mask(faces(1:3,face)) = 1
!
!  For quadratic tetrahedrons, we need to add three more side nodes.
!  We retrieve the face index by F = FACES(4,*).
!  We can determine the local indices from the value of F.
!  We can determine the global indices from ELEMENT_NODE.
!
    if ( element_order == 10 ) then

      f = faces(4,face)
      element = faces(5,face)

      if ( f == 1 ) then
        boundary_node_mask(element_node(8,element)) = 1
        boundary_node_mask(element_node(10,element)) = 1
        boundary_node_mask(element_node(9,element)) = 1
      else if ( f == 2 ) then
        boundary_node_mask(element_node(6,element)) = 1
        boundary_node_mask(element_node(10,element)) = 1
        boundary_node_mask(element_node(7,element)) = 1
      else if ( f == 3 ) then
        boundary_node_mask(element_node(5,element)) = 1
        boundary_node_mask(element_node(9,element)) = 1
        boundary_node_mask(element_node(7,element)) = 1
      else if ( f == 4 ) then
        boundary_node_mask(element_node(5,element)) = 1
        boundary_node_mask(element_node(8,element)) = 1
        boundary_node_mask(element_node(6,element)) = 1
      end if

    end if

    face = face + 1

  end do

  boundary_node_num = sum ( boundary_node_mask(1:node_num) )

  deallocate ( faces )
end

subroutine tet_mesh_boundary_set ( element_order, element_num, element_node, &
  boundary_element_order, boundary_element_num, boundary_element_node )

!*****************************************************************************80
!
!! TET_MESH_BOUNDARY_SET sets the boundary faces in a tet mesh.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 2010
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
!    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), the
!    nodes that make up each element.
!
!    Input, integer BOUNDARY_ELEMENT_ORDER, the order of the
!    boundary faces.
!
!    Input, integer BOUNDARY_ELEMENT_NUM, the number of 
!    boundary faces.
!
!    Output, integer !    BOUNDARY_ELEMENT_NODE(BOUNDARY_ELEMENT_ORDER,BOUNDARY_ELEMENT_NUM),
!    the nodes that make up each boundary face.
!
  implicit none

  integer boundary_element_num
  integer boundary_element_order
  integer element_num
  integer element_order

  integer a
  integer b
  integer boundary_element
  integer boundary_element_node(boundary_element_order,&
    boundary_element_num)
  integer c
  integer element
  integer element_node(element_order,element_num)
  integer element1
  integer element2
  integer f
  integer face
  integer face1
  integer face2
  integer , allocatable, dimension ( :, : ) :: faces
  integer i
  integer j
  integer k
  integer l

  if ( element_order /= 4 .and. element_order /= 10 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TET_MESH_BOUNDARY_SET - Fatal error!'
    write ( *, '(a)' ) '  Unexpected element order.'
    stop
  end if

  allocate ( faces(1:5,1:4*element_num) )
!
!  Step 1.
!  From the list of nodes forming tetrahedron T, of the form: 
!
!    (I,J,K,L)
!
!  or
!
!    (I,J,K,L,I+J,I+K,I+L,J+K,J+L,K+L),
!    (1,2,3,4, 5,  6,  7,  8,  9, 10 ),
!
!  construct the four face relations:
!
!    F = 1: (J,K,L,F,T)
!    F = 2: (I,K,L,F,T)
!    F = 3: (I,J,L,F,T)
!    F = 4: (I,J,K,F,T)
!
!  If T is actually order 10, we can retrieve the indices of the midside
!  nodes from the values of F and T.  In that case, the 4 faces are:
!
!    F = 1: 2, 3, 4, 8, 10, 9
!    F = 2: 1, 3, 4, 6, 10, 7
!    F = 3: 1, 2, 4, 5,  9, 7
!    F = 4: 1, 2, 3, 5,  8, 6
!
!  In order to make matching easier, we reorder each triple of nodes
!  into ascending order.
!
  do element = 1, element_num

    i = element_node(1,element)
    j = element_node(2,element)
    k = element_node(3,element)
    l = element_node(4,element)

    call i4i4i4_sort_a ( j, k, l, a, b, c )

    faces(1:5,4*(element-1)+1) = (/ a, b, c, 1, element /)

    call i4i4i4_sort_a ( i, k, l, a, b, c )

    faces(1:5,4*(element-1)+2) = (/ a, b, c, 2, element /)

    call i4i4i4_sort_a ( i, j, l, a, b, c )

    faces(1:5,4*(element-1)+3) = (/ a, b, c, 3, element /)

    call i4i4i4_sort_a ( i, j, k, a, b, c )

    faces(1:5,4*(element-1)+4) = (/ a, b, c, 4, element /)

  end do
!
!  Step 2. Perform an ascending dictionary sort on the neighbor relations.
!  We only intend to sort on rows 1:3; the routine we call here
!  sorts on rows 1 through 5 but that won't hurt us.
!
!  What we need is to find cases where two tetrahedrons share a face.
!  By sorting the columns of the FACES array, we will put shared faces
!  next to each other.
!
  call i4col_sort_a ( 5, 4*element_num, faces )
!
!  Step 3. Neighboring faces show up as consecutive columns with
!  identical first three entries.  Count columns which don't have
!  a following column that matches the first three entries.
!
  boundary_element = 0
  face = 1

  do while ( face <= 4 * element_num )

    if ( face < 4 * element_num ) then

      if ( all ( faces(1:3,face) == faces(1:3,face+1) ) ) then
        face = face + 2
        cycle
      end if

    end if

    boundary_element = boundary_element + 1

    f = faces(4,face)
    element = faces(5,face)

    if ( f == 1 ) then
      boundary_element_node(1,boundary_element) = element_node(2,element)
      boundary_element_node(2,boundary_element) = element_node(3,element)
      boundary_element_node(3,boundary_element) = element_node(4,element)
    else if ( f == 2 ) then
      boundary_element_node(1,boundary_element) = element_node(1,element)
      boundary_element_node(2,boundary_element) = element_node(3,element)
      boundary_element_node(3,boundary_element) = element_node(4,element)
    else if ( f == 3 ) then
      boundary_element_node(1,boundary_element) = element_node(1,element)
      boundary_element_node(2,boundary_element) = element_node(2,element)
      boundary_element_node(3,boundary_element) = element_node(4,element)
    else if ( f == 4 ) then
      boundary_element_node(1,boundary_element) = element_node(1,element)
      boundary_element_node(2,boundary_element) = element_node(2,element)
      boundary_element_node(3,boundary_element) = element_node(3,element)
    end if
!
!  For quadratic tetrahedrons, we need to add three more side nodes.
!
    if ( element_order == 10 ) then

      if ( f == 1 ) then
        boundary_element_node(4,boundary_element) = element_node(8,element)
        boundary_element_node(5,boundary_element) = element_node(10,element)
        boundary_element_node(6,boundary_element) = element_node(9,element)
      else if ( f == 2 ) then
        boundary_element_node(4,boundary_element) = element_node(6,element)
        boundary_element_node(5,boundary_element) = element_node(10,element)
        boundary_element_node(6,boundary_element) = element_node(7,element)
      else if ( f == 3 ) then
        boundary_element_node(4,boundary_element) = element_node(5,element)
        boundary_element_node(5,boundary_element) = element_node(9,element)
        boundary_element_node(6,boundary_element) = element_node(7,element)
      else if ( f == 4 ) then
        boundary_element_node(4,boundary_element) = element_node(5,element)
        boundary_element_node(5,boundary_element) = element_node(8,element)
        boundary_element_node(6,boundary_element) = element_node(6,element)
      end if

    end if

    face = face + 1

  end do

  deallocate ( faces )
end
