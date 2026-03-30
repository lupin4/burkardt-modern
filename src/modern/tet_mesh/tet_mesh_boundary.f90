!> tet_mesh_boundary — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

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
