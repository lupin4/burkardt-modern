!> bezier_surface — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

subroutine bezier_patch_evaluate ( node_num, node_xyz, rectangle_num, &
  rectangle_node, patch, point_num, point_uv, point_xyz )

!*****************************************************************************80
!
!! BEZIER_PATCH_EVALUATE evaluates a Bezier patch.
!
!  Discussion:
!
!    Given a Bezier surface defined as a collection of patches,
!    it is desired to evaluate the surface, in a single patch,
!    for a variety of (U,V) parameter values.
!
!    The user specifies the patch index, and the (U,V) coordinates of
!    several points.  The routine returns the (X,Y,Z) coordinates of the
!    points.  The special parameter values (0,0), (1,0), (0,1) and (1,1) should
!    return the coordinates of the SW, SE, NW and NE corners of the patch.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer NODE_XYZ(3,NODE_NUM), the node coordinates.
!
!    Input, integer RECTANGLE_NUM, the number of rectangles.
!
!    Input, integer RECTANGLE_NODE(16,RECTANGLE_NUM), the nodes that make up
!    each rectangle.
!
!    Input, integer PATCH, the index of the Bezier patch.
!
!    Input, integer POINT_NUM, the number of points at which evaluation
!    is desired.
!
!    Input, double precision POINT_UV(2,POINT_NUM), the (U,V) parameter
!    coordinates of the points.
!
!    Output, double precision POINT_XYZ(3,POINT_NUM), the (X,Y,Z)
!    physical coordinates of the points.
!
  implicit none

  integer node_num
  integer point_num
  integer rectangle_num

  integer i
  double precision node_xyz(3,node_num)
  integer patch
  double precision patch_xmat(4,4)
  double precision patch_ymat(4,4)
  double precision patch_zmat(4,4)
  integer point
  double precision point_uv(2,point_num)
  double precision point_xyz(3,point_num)
  integer rectangle_node(16,rectangle_num)
  double precision u
  double precision uvec(4)
  double precision v
  double precision vvec(4)

  patch_xmat(1:4,1:4) = reshape ( node_xyz(1,rectangle_node(1:16,patch)), &
    (/ 4, 4 /) )
  patch_ymat(1:4,1:4) = reshape ( node_xyz(2,rectangle_node(1:16,patch)), &
    (/ 4, 4 /) )
  patch_zmat(1:4,1:4) = reshape ( node_xyz(3,rectangle_node(1:16,patch)), &
    (/ 4, 4 /) )

  do point = 1, point_num

    u = point_uv(1,point)
    v = point_uv(2,point)

    uvec(1:4) = (/                  ( 1.0D+00 - u )**3, &
                   3.0D+00 * u    * ( 1.0D+00 - u )**2, &
                   3.0D+00 * u**2 * ( 1.0D+00 - u ),    &
                             u**3 /)

    vvec(1:4) = (/                  ( 1.0D+00 - v )**3, &
                   3.0D+00 * v    * ( 1.0D+00 - v )**2, &
                   3.0D+00 * v**2 * ( 1.0D+00 - v ),    &
                             v**3 /)

    point_xyz(1,point) = dot_product ( uvec(1:4), &
      matmul ( patch_xmat(1:4,1:4), vvec(1:4) ) )

    point_xyz(2,point) = dot_product ( uvec(1:4), &
      matmul ( patch_ymat(1:4,1:4), vvec(1:4) ) )

    point_xyz(3,point) = dot_product ( uvec(1:4), &
      matmul ( patch_zmat(1:4,1:4), vvec(1:4) ) )

  end do
end

subroutine bezier_surface_neighbors ( rectangle_num, &
  rectangle_node, rectangle_neighbor )

!*****************************************************************************80
!
!! BEZIER_SURFACE_NEIGHBORS determines Bezier rectangle neighbors.
!
!  Discussion:
!
!    A (bicubic) Bezier surface is constructed of patches.  Each
!    patch is an (X,Y,Z) image of the unit rectangle in (U,V) space.
!    Two patches may be said to be "neighbors" if their images share
!    a side.
!
!    It is perfectly possible for each Bezier patch to have NO neighbors.
!    It is perfectly possible for each Bezier patch to have many neighbors
!    sharing a single side.
!
!    However, in the most common case, and the one handled here, we
!    may assume that the (X,Y,Z) patches fit together in a way similar
!    to the way quilt patches form a surface, so that most rectangles
!    have four neighbors, one on each side, while boundary rectangles
!    might have fewer neighbors.
!
!    This routine creates a data structure recording the neighbor information.
!
!    The primary amount of work occurs in sorting a list of
!    4 * RECTANGLE_NUM data items representing the sides of each patch.
!
!    The nodes of a single rectangle are assumed to be numbered as follows:
!
!    V
!    |   13-14-15-16
!    |    |  |  |  |
!    |    9-10-11-12
!    |    |  |  |  |
!    |    5--6--7--8
!    |    |  |  |  |
!    |    1--2--3--4
!    |
!    +--------------->U
!
!    We assume that a neighbor patch will agree in its (X,Y,Z) values
!    for all four nodes along a single side (although the nodes might
!    be listed in reverse order.)
!
!    We assume that there is at most one neighbor patch on each side.
!    (However, even in a simple case like the Utah teapot, this is not
!    true, because at the top and bottom, many patches have one side
!    that degenerates to a single point, and several such patches
!    meet at that point!)
!
!    We choose to number the sides of each patch as follows:
!
!
!    V      SIDE 3
!    |
!    |   13-14-15-16
!    | S  |  |  |  |
!    | I  9-10-11-12
!    | D  |  |  |  | SIDE 2
!    | E  5--6--7--8
!    |    |  |  |  |
!    | 4  1--2--3--4
!    |       SIDE 1
!    +--------------->U
!
!    And these indices for the sides correspond to the first index
!    in the RECTANGLE_NEIGHBOR array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer RECTANGLE_NUM, the number of rectangles.
!
!    Input, integer RECTANGLE_NODE(16,RECTANGLE_NUM), the nodes that make up
!    each rectangle.
!
!    Output, integer RECTANGLE_NEIGHBOR(4,RECTANGLE_NUM), the indices of
!    the rectanges that are direct neighbors of a given rectangle.
!    RECTANGLE_NEIGHBOR(1,I) is the index of the rectangle which touches
!    side 1, defined by nodes 1 and 4, and so on.  RECTANGLE_NEIGHBOR(1,I)
!    is negative if there is no neighbor on that side.
!
  implicit none

  integer rectangle_num

  integer i1
  integer i2
  integer i3
  integer i4
  integer irow
  integer rectangle
  integer rectangle1
  integer rectangle2
  integer row(4*rectangle_num,4)
  integer rectangle_node(16,rectangle_num)
  integer rectangle_neighbor(4,rectangle_num)
  integer side1
  integer side2
!
!  Step 1.
!  From the list of vertices for rectangle T,
!  construct records that describe the four side segments, by listing
!  the first and last nodes of each segment.
!  To make matching easier later, we sort each pair of nodes.
!
  do rectangle = 1, rectangle_num

    i1 = rectangle_node(1,rectangle)
    i2 = rectangle_node(4,rectangle)
    i3 = rectangle_node(13,rectangle)
    i4 = rectangle_node(16,rectangle)

    if ( i1 < i2 ) then
      row(4*(rectangle-1)+1,1:4) = (/ i1, i2, 1, rectangle /)
    else
      row(4*(rectangle-1)+1,1:4) = (/ i2, i1, 1, rectangle /)
    end if

    if ( i2 < i4 ) then
      row(4*(rectangle-1)+2,1:4) = (/ i2, i4, 2, rectangle /)
    else
      row(4*(rectangle-1)+2,1:4) = (/ i4, i2, 2, rectangle /)
    end if

    if ( i4 < i3 ) then
      row(4*(rectangle-1)+3,1:4) = (/ i4, i3, 3, rectangle /)
    else
      row(4*(rectangle-1)+3,1:4) = (/ i3, i4, 3, rectangle /)
    end if

    if ( i3 < i1 ) then
      row(4*(rectangle-1)+4,1:4) = (/ i3, i1, 4, rectangle /)
    else
      row(4*(rectangle-1)+4,1:4) = (/ i1, i3, 4, rectangle /)
    end if

  end do
!
!  Step 2. Perform an ascending dictionary sort on the neighbor relations.
!  We only intend to sort on columns 1 and 2; the routine we call here
!  sorts on columns 1 through 3 but that won't hurt us.
!
!  What we need is to find cases where two rectangles share an edge.
!  Say they share an edge defined by the nodes I and J.  Then there are
!  two rows of ROW that start out ( I, J, ?, ? ).  By sorting ROW,
!  we make sure that these two rows occur consecutively.  That will
!  make it easy to notice that the rectangles are neighbors.
!
  call i4row_sort_a ( 4*rectangle_num, 4, row )
!
!  Step 3. Neighboring rectangles show up as consecutive rows with
!  identical first two entries.  Whenever you spot this happening,
!  make the appropriate entries in RECTANGLE_NEIGHBOR.
!
  rectangle_neighbor(1:4,1:rectangle_num) = -1

  irow = 1

  do

    if ( 4 * rectangle_num <= irow ) then
      exit
    end if

    if ( row(irow,1) /= row(irow+1,1) .or. row(irow,2) /= row(irow+1,2) ) then
      irow = irow + 1
      cycle
    end if

    side1 = row(irow,3)
    rectangle1 = row(irow,4)
    side2 = row(irow+1,3)
    rectangle2 = row(irow+1,4)

    rectangle_neighbor(side1,rectangle1) = rectangle2
    rectangle_neighbor(side2,rectangle2) = rectangle1

    irow = irow + 2

  end do
end

subroutine bezier_surface_node_print ( node_num, node_xyz )

!*****************************************************************************80
!
!! BEZIER_SURFACE_NODE_PRINT prints nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, double precision NODE_XYZ(3,NODE_NUM), the coordinates of the
!    nodes.
!
  implicit none

  integer , parameter :: dim_num = 3
  integer node_num

  double precision node_xyz(dim_num,node_num)

  call r8mat_transpose_print ( dim_num, node_num, node_xyz, &
    '  Bezier Surface Nodes:' )
end

subroutine bezier_surface_node_read ( node_file_name, node_num, node_xyz )

!*****************************************************************************80
!
!! BEZIER_SURFACE_NODE_READ reads nodes from a node file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) NODE_FILE_NAME, the name of the node file.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Output, double precision NODE_XYZ(3,NODE_NUM), the coordinates of the
!    nodes.
!
  implicit none

  integer , parameter :: dim_num = 3
  integer node_num

  character ( len = * ) node_file_name
  double precision node_xyz(dim_num,node_num)

  call r8mat_data_read ( node_file_name, dim_num, node_num, node_xyz )
end

subroutine bezier_surface_node_size ( node_file_name, node_num )

!*****************************************************************************80
!
!! BEZIER_SURFACE_NODE_SIZE counts nodes in a node file.
!
!  Discussion:
!
!    This version of the routine simply counts the number of lines
!    in the file (although it ignores comment lines beginning with "#").
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) NODE_FILE_NAME, the name of the node file.
!
!    Output, integer NODE_NUM, the number of nodes.
!
  implicit none

  character ( len = * ) node_file_name
  integer node_num

  call file_row_count ( node_file_name, node_num )
end

subroutine bezier_surface_node_write ( node_file_name, node_num, node_xyz )

!*****************************************************************************80
!
!! BEZIER_SURFACE_NODE_WRITE writes nodes to a node file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) NODE_FILE_NAME, the name of the node file.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, double precision NODE_XYZ(3,NODE_NUM), the coordinates of the
!    nodes.
!
  implicit none

  integer , parameter :: dim_num = 3
  integer node_num

  logical, parameter :: header = .true.
  character ( len = * ) node_file_name
  double precision node_xyz(dim_num,node_num)

  call r8mat_write ( node_file_name, dim_num, node_num, node_xyz, header )
end

subroutine bezier_surface_rectangle_print ( rectangle_num, rectangle_node )

!*****************************************************************************80
!
!! BEZIER_SURFACE_RECTANGLE_PRINT prints rectangles.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer RECTANGLE_NUM, the number of rectangles.
!
!    Input, integer RECTANGLE_NODE(16,RECTANGLE_NUM), the nodes
!    that make up each rectangle.
!
  implicit none

  integer , parameter :: dim_num = 16
  integer rectangle_num

  integer rectangle_node(dim_num,rectangle_num)

  call i4mat_transpose_print ( dim_num, rectangle_num, rectangle_node, &
    '  Bezier Rectangles:' )
end

subroutine bezier_surface_rectangle_read ( rectangle_file_name, &
  rectangle_num, rectangle_node )

!*****************************************************************************80
!
!! BEZIER_SURFACE_RECTANGLE_READ reads rectangles from a rectangle file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) RECTANGLE_FILE_NAME, the name of the
!    rectangle file.
!
!    Input, integer RECTANGLE_NUM, the number of rectangles.
!
!    Output, double precision RECTANGLE_NODE(16,RECTANGLE_NUM),
!    the nodes that make up each rectangle.
!
  implicit none

  integer , parameter :: dim_num = 16
  integer rectangle_num

  character ( len = * ) rectangle_file_name
  integer rectangle_node(dim_num,rectangle_num)

  call i4mat_data_read ( rectangle_file_name, dim_num, rectangle_num, &
    rectangle_node )
end

subroutine bezier_surface_rectangle_size ( rectangle_file_name, rectangle_num )

!*****************************************************************************80
!
!! BEZIER_SURFACE_RECTANGLE_SIZE counts rectangles in a rectangle file.
!
!  Discussion:
!
!    This version of the routine simply counts the number of lines
!    in the file (although it ignores comment lines beginning with "#").
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) RECTANGLE_FILE_NAME, the name of the
!    rectangle file.
!
!    Output, integer NODE_NUM, the number of rectangles.
!
  implicit none

  character ( len = * ) rectangle_file_name
  integer rectangle_num

  call file_row_count ( rectangle_file_name, rectangle_num )
end

subroutine bezier_surface_rectangle_write ( rectangle_file_name, &
  rectangle_num, rectangle_node )

!*****************************************************************************80
!
!! BEZIER_SURFACE_RECTANGLE_WRITE writes rectangle to a rectangle file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) RECTANGLE_FILE_NAME, the name of the
!    rectangle file.
!
!    Input, integer RECTANGLE_NUM, the number of rectangle.
!
!    Input, integer RECTANGLE_NODE(16,RECTANGLE_NUM), the nodes that
!    make up each rectangle.
!
  implicit none

  integer , parameter :: dim_num = 16
  integer rectangle_num

  logical, parameter :: header = .true.
  character ( len = * ) rectangle_file_name
  integer rectangle_node(dim_num,rectangle_num)

  call i4mat_write ( rectangle_file_name, dim_num, rectangle_num, &
    rectangle_node, header )
end

subroutine i4row_compare ( m, n, a, i, j, isgn )

!*****************************************************************************80
!
!! I4ROW_COMPARE compares two rows of an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of integer values, regarded
!    as an array of M rows of length N.
!
!  Example:
!
!    Input:
!
!    M = 3, N = 4, I = 2, J = 3
!
!    A = (
!      1  2  3  4
!      5  6  7  8
!      9 10 11 12 )
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
!    14 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(M,N), an array of M rows of vectors of length N.
!
!    Input, integer I, J, the rows to be compared.
!    I and J must be between 1 and M.
!
!    Output, integer ISGN, the results of the comparison:
!    -1, row I < row J,
!     0, row I = row J,
!    +1, row J < row I.
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
!  Check that I and J are legal.
!
  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index I is less than 1.'
    write ( *, '(a,i8)' ) '  I = ', i
    stop
  else if ( m < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index I is out of bounds.'
    write ( *, '(a,i8)' ) '  I = ', i
    write ( *, '(a,i8)' ) '  Maximum legal value is M = ', m
    stop
  end if

  if ( j < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index J is less than 1.'
    write ( *, '(a,i8)' ) '  J = ', j
    stop
  else if ( m < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index J is out of bounds.'
    write ( *, '(a,i8)' ) '  J = ', j
    write ( *, '(a,i8)' ) '  Maximum legal value is M = ', m
    stop
  end if

  isgn = 0

  if ( i == j ) then
  end if

  k = 1

  do while ( k <= n )

    if ( a(i,k) < a(j,k) ) then
      isgn = -1
      return
    else if ( a(j,k) < a(i,k) ) then
      isgn = +1
    end if

    k = k + 1

  end do
end

subroutine i4row_sort_a ( m, n, a )

!*****************************************************************************80
!
!! I4ROW_SORT_A ascending sorts the rows of an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of integer values, regarded
!    as an array of M rows of length N.
!
!    In lexicographic order, the statement "X < Y", applied to two
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, X is less than Y if, at the first index where they
!    differ, the X value is less than the Y value.
!
!  Example:
!
!    Input:
!
!      M = 5, N = 3
!
!      A =
!        3  2  1
!        2  4  3
!        3  1  8
!        2  4  2
!        1  9  9
!
!    Output:
!
!      A =
!        1  9  9
!        2  4  2
!        2  4  3
!        3  1  8
!        3  2  1
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
!    Input, integer M, the number of rows of A.
!
!    Input, integer N, the number of columns of A.
!
!    Input/output, integer A(M,N).
!    On input, the array of M rows of N-vectors.
!    On output, the rows of A have been sorted in ascending
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

  if ( m <= 1 ) then
  end if

  if ( n <= 0 ) then
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

    call sort_heap_external ( m, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call i4row_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4row_compare ( m, n, a, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do
end

subroutine i4row_swap ( m, n, a, i1, i2 )

!*****************************************************************************80
!
!! I4ROW_SWAP swaps two rows of an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of integer values, regarded
!    as an array of M rows of length N.
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
!    Input, integer M, N, the number of rows and columns.
!
!    Input/output, integer A(M,N), an array of data.
!
!    Input, integer I1, I2, the two rows to swap.
!
  implicit none

  integer m
  integer n

  integer a(m,n)
  integer i1
  integer i2
  integer row(n)
!
!  Check.
!
  if ( i1 < 1 .or. m < i1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_SWAP - Fatal error!'
    write ( *, '(a)' ) '  I1 is out of range.'
    stop
  end if

  if ( i2 < 1 .or. m < i2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_SWAP - Fatal error!'
    write ( *, '(a)' ) '  I2 is out of range.'
    stop
  end if

  if ( i1 == i2 ) then
  end if

  row(1:n)  = a(i1,1:n)
  a(i1,1:n) = a(i2,1:n)
  a(i2,1:n) = row(1:n)
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
!    FORTRAN77 original version by Albert Nijenhuis, Herbert Wilf.
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
!    Input, integer ISGN, results of comparison of elements I and J.
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
