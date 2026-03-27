!> tet_mesh_quality � Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module tet_mesh_quality_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: i4vec_histogram, mesh_base_one, r8_swap, r8mat_det_4d, r8mat_solve, r8vec_cross_3d
  public :: r8vec_length, r8vec_mean, r8vec_variance, tet_mesh_node_order, tet_mesh_quality1, tet_mesh_quality2
  public :: tet_mesh_quality3, tet_mesh_quality4, tet_mesh_quality5, tetrahedron_circumsphere_3d, tetrahedron_edge_length_3d, tetrahedron_insphere_3d
  public :: tetrahedron_quality1_3d, tetrahedron_quality2_3d, tetrahedron_quality3_3d, tetrahedron_quality4_3d, tetrahedron_volume_3d

contains

  pure subroutine i4vec_histogram ( n, a, histo_num, histo_gram ) &
        bind(C, name="i4vec_histogram")

  !*****************************************************************************80
  !
  !! I4VEC_HISTOGRAM computes a histogram of the elements of an I4VEC.
  !
  !  Discussion:
  !
  !    It is assumed that the entries in the vector A are nonnegative.
  !    Only values between 0 and HISTO_NUM will be histogrammed.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    29 August 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of elements of A.
  !
  !    Input, integer(ip) A(N), the array to examine.
  !
  !    Input, integer(ip) HISTO_NUM, the maximum value for which a
  !    histogram entry will be computed.
  !
  !    Output, integer(ip) HISTO_GRAM(0:HISTO_NUM), contains the number of
  !    entries of A with the values of 0 through HISTO_NUM.
  !

    integer(ip), intent(in), value :: histo_num
    integer(ip), intent(in), value :: n

    integer(ip), intent(in) :: a(n)
    integer(ip), intent(out) :: histo_gram(0:histo_num)
    integer(ip) :: i

    histo_gram(0:histo_num) = 0

    do i = 1, n

      if ( 0 <= a(i) .and. a(i) <= histo_num ) then
        histo_gram(a(i)) = histo_gram(a(i)) + 1
      end if

    end do
  end subroutine i4vec_histogram

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

    integer(ip), intent(in), value :: element_num
    integer(ip), intent(in), value :: element_order

    integer(ip) :: element
    integer(ip), intent(inout) :: element_node(element_order,element_num)
    integer(ip) :: node
    integer(ip) :: node_max
    integer(ip) :: node_min
    integer(ip), intent(in), value :: node_num
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
      write ( *, '(a)' )'MESH_BASE_ZERO:'
      write ( *, '(a)' )'  The element indexing appears to be 1-based!'
      write ( *, '(a)' )'  No conversion is necessary.'
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )'MESH_BASE_ZERO - Warning!'
      write ( *, '(a)' )' The element indexing is not of a recognized type.'
    end if
  end subroutine mesh_base_one

  pure subroutine r8_swap ( x, y ) &
        bind(C, name="r8_swap")

  !*****************************************************************************80
  !
  !! R8_SWAP switches two R8's.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    01 May 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, real(dp) X, Y.  On output, the values of X and
  !    Y have been interchanged.
  !

    real(dp), intent(inout) :: x
    real(dp), intent(inout) :: y
    real(dp) :: z

    z = x
    x = y
    y = z
  end subroutine r8_swap

  pure function r8mat_det_4d ( a ) &
        bind(C, name="r8mat_det_4d")

  !*****************************************************************************80
  !
  !! R8MAT_DET_4D computes the determinant of a 4 by 4 matrix.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    16 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real(dp) A(4,4), the matrix whose determinant is desired.
  !
  !    Output, real(dp) R8MAT_DET_4D, the determinant of the matrix.
  !

    real(dp), intent(in) :: a(4,4)
    real(dp) :: r8mat_det_4d

    r8mat_det_4d = &
        a(1,1) * ( &
          a(2,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
        - a(2,3) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
        + a(2,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) ) &
      - a(1,2) * ( &
          a(2,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
        - a(2,3) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
        + a(2,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) ) &
      + a(1,3) * ( &
          a(2,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
        - a(2,2) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
        + a(2,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) ) &
      - a(1,4) * ( &
          a(2,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
        - a(2,2) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
        + a(2,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) )
  end function r8mat_det_4d

  subroutine r8mat_solve ( n, rhs_num, a, info ) &
        bind(C, name="r8mat_solve")

  !*****************************************************************************80
  !
  !! R8MAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    29 August 2003
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the order of the matrix.
  !
  !    Input, integer(ip) rhs_num, the number of right hand sides.  rhs_num
  !    must be at least 0.
  !
  !    Input/output, real(dp) A(N,N+rhs_num), contains in rows and
  !    columns 1 to N the coefficient matrix, and in columns N+1 through
  !    N+rhs_num, the right hand sides.  On output, the coefficient matrix
  !    area has been destroyed, while the right hand sides have
  !    been overwritten with the corresponding solutions.
  !
  !    Output, integer(ip) INFO, singularity flag.
  !    0, the matrix was not singular, the solutions were computed;
  !    J, factorization failed on step J, and the solutions could not
  !    be computed.
  !

    integer(ip), intent(in), value :: n
    integer(ip), intent(in), value :: rhs_num

    real(dp), intent(inout) :: a(n,n+rhs_num)
    real(dp) :: apivot
    real(dp) :: factor
    integer(ip) :: i
    integer(ip), intent(out) :: info
    integer(ip) :: ipivot
    integer(ip) :: j

    info = 0

    do j = 1, n
  !
  !  Choose a pivot row.
  !
      ipivot = j
      apivot = a(j,j)

      do i = j+1, n
        if ( abs ( apivot ) < abs ( a(i,j) ) ) then
          apivot = a(i,j)
          ipivot = i
        end if
      end do

      if ( apivot == 0.0_dp ) then
        info = j
      end if
  !
  !  Interchange.
  !
      do i = 1, n + rhs_num
        call r8_swap ( a(ipivot,i), a(j,i) )
      end do
  !
  !  A(J,J) becomes 1.
  !
      a(j,j) = 1.0_dp
      a(j,j+1:n+rhs_num) = a(j,j+1:n+rhs_num) / apivot
  !
  !  A(I,J) becomes 0.
  !
      do i = 1, n

        if ( i /= j ) then

          factor = a(i,j)
          a(i,j) = 0.0_dp
          a(i,j+1:n+rhs_num) = a(i,j+1:n+rhs_num) - factor * a(j,j+1:n+rhs_num)

        end if

      end do

    end do
  end subroutine r8mat_solve

  pure subroutine r8vec_cross_3d ( v1, v2, v3 ) &
        bind(C, name="r8vec_cross_3d")

  !*****************************************************************************80
  !
  !! R8VEC_CROSS_3D computes the cross product of two vectors in 3D.
  !
  !  Discussion:
  !
  !    The cross product in 3D can be regarded as the determinant of the
  !    symbolic matrix:
  !
  !          |  i  j  k |
  !      det | x1 y1 z1 |
  !          | x2 y2 z2 |
  !
  !      = ( y1 * z2 - z1 * y2 ) * i
  !      + ( z1 * x2 - x1 * z2 ) * j
  !      + ( x1 * y2 - y1 * x2 ) * k
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 August 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real(dp) V1(3), V2(3), the two vectors.
  !
  !    Output, real(dp) V3(3), the cross product vector.
  !

    integer(ip), parameter :: dim_num = 3

    real(dp), intent(in) :: v1(dim_num)
    real(dp), intent(in) :: v2(dim_num)
    real(dp), intent(out) :: v3(dim_num)

    v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
    v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
    v3(3) = v1(1) * v2(2) - v1(2) * v2(1)
  end subroutine r8vec_cross_3d

  pure function r8vec_length ( dim_num, x ) &
        bind(C, name="r8vec_length")

  !*****************************************************************************80
  !
  !! R8VEC_LENGTH returns the Euclidean length of an R8VEC.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    08 August 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) DIM_NUM, the spatial dimension.
  !
  !    Input, real(dp) X(DIM_NUM), the vector.
  !
  !    Output, real(dp) R8VEC_LENGTH, the Euclidean length of the vector.
  !

    integer(ip), intent(in), value :: dim_num

    real(dp) :: r8vec_length
    real(dp), intent(in) :: x(dim_num)

    r8vec_length = sqrt ( sum ( ( x(1:dim_num) )**2 ) )
  end function r8vec_length

  pure subroutine r8vec_mean ( n, a, mean ) &
        bind(C, name="r8vec_mean")

  !*****************************************************************************80
  !
  !! R8VEC_MEAN returns the mean of an R8VEC.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    02 February 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of entries in the vector.
  !
  !    Input, real(dp) A(N), the vector whose mean is desired.
  !
  !    Output, real(dp) MEAN, the mean of the vector entries.
  !

    integer(ip), intent(in), value :: n

    real(dp), intent(in) :: a(n)
    real(dp), intent(out) :: mean

    mean = sum ( a(1:n) ) / real ( n, dp)
  end subroutine r8vec_mean

  pure subroutine r8vec_variance ( n, a, variance ) &
        bind(C, name="r8vec_variance")

  !*****************************************************************************80
  !
  !! R8VEC_VARIANCE returns the variance of an R8VEC.
  !
  !  Discussion:
  !
  !    The variance of a vector X of length N is defined as
  !
  !      mean ( X(1:n) ) = sum ( X(1:n) ) / n
  !
  !      var ( X(1:n) ) = sum ( ( X(1:n) - mean )**2 ) / ( n - 1 )
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 February 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of entries in the vector.
  !    N should be at least 2.
  !
  !    Input, real(dp) A(N), the vector.
  !
  !    Output, real(dp) VARIANCE, the variance of the vector.
  !

    integer(ip), intent(in), value :: n

    real(dp), intent(in) :: a(n)
    real(dp) :: mean
    real(dp), intent(out) :: variance

    if ( n < 2 ) then

      variance = 0.0_dp

    else

      mean = sum ( a(1:n) ) / real ( n, dp)

      variance = sum ( ( a(1:n) - mean )**2 )

      variance = variance / real ( n - 1, dp)

    end if
  end subroutine r8vec_variance

  subroutine tet_mesh_node_order ( tetra_order, tetra_num, tetra_node, &
    node_num, node_order ) &
        bind(C, name="tet_mesh_node_order")

  !*****************************************************************************80
  !
  !! TET_MESH_NODE_ORDER: determine the order of nodes in a tetra mesh.
  !
  !  Discussion:
  !
  !    The order of a node is the number of tetrahedrons that use that node
  !    as a vertex.
  !
  !    Tetrahedrons of order 4 or 10 are allowed as input.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 October 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) TETRA_ORDER, the order of the mesh, either 4 or 10.
  !
  !    Input, integer(ip) TETRA_NUM, the number of tetrahedrons.
  !
  !    Input, integer(ip) TETRA_NODE(TETRA_ORDER,TETRA_NUM), the nodes
  !    that make up the tetrahedrons.
  !
  !    Input, integer(ip) NODE_NUM, the number of nodes.
  !
  !    Output, integer(ip) NODE_ORDER(NODE_NUM), the order of each node.
  !

    integer(ip), intent(in), value :: node_num
    integer(ip), intent(in), value :: tetra_num
    integer(ip), intent(in), value :: tetra_order

    integer(ip) :: i
    integer(ip) :: node
    integer(ip), intent(out) :: node_order(node_num)
    integer(ip) :: tetra
    integer(ip), intent(in) :: tetra_node(tetra_order,tetra_num)

    node_order(1:node_num) = 0

    do tetra = 1, tetra_num
      do i = 1, tetra_order
        node = tetra_node(i,tetra)
        if ( node < 1 .or. node_num < node ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TET_MESH_NODE_ORDER - Fatal error!'
          write ( *, '(a)' ) '  Illegal entry in TETRA_NODE.'
          stop
        else
          node_order(node) = node_order(node) + 1
        end if
      end do
    end do
  end subroutine tet_mesh_node_order

  subroutine tet_mesh_quality1 ( node_num, node_xyz, tetra_order, tetra_num, &
    tetra_node, tetra_quality ) &
        bind(C, name="tet_mesh_quality1")

  !*****************************************************************************80
  !
  !! TET_MESH_QUALITY1 returns the quality of each tet in a mesh.
  !
  !  Discussion:
  !
  !    The overall tet mesh quality measure is the minimum of the corresponding
  !    tetrahedron quality measure, over all tetrahedrons in the tet mesh.
  !
  !    Although tetrahedrons of order 10 are allowed as input,
  !    only the first 4 nodes (presumed to be the vertices) are used
  !    in the calculation.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 November 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) NODE_NUM, the number of nodes.
  !
  !    Input, real(dp) NODE_XYZ(3,NODE_NUM), the nodes.
  !
  !    Input, integer(ip) TETRA_ORDER, the order of the mesh, either 4 or 10.
  !
  !    Input, integer(ip) TETRA_NUM, the number of tetrahedrons.
  !
  !    Input, integer(ip) TETRA_NODE(TETRA_ORDER,TETRA_NUM), the nodes
  !    that make up the tetrahedrons.
  !
  !    Output, real(dp) TETRA_QUALITY(TETRA_NUM), the quality
  !    measure for each tetrahedron.
  !

    integer(ip), parameter :: dim_num = 3
    integer(ip), intent(in), value :: node_num
    integer(ip), intent(in), value :: tetra_num
    integer(ip), intent(in), value :: tetra_order

    real(dp), intent(in) :: node_xyz(dim_num,node_num)
    integer(ip) :: tetra
    integer(ip), intent(in) :: tetra_node(tetra_order,tetra_num)
    real(dp), intent(out) :: tetra_quality(tetra_num)
    real(dp) :: tetra_xyz(dim_num,4)

    do tetra = 1, tetra_num

      tetra_xyz(1:dim_num,1:4) = node_xyz(1:dim_num,tetra_node(1:4,tetra))

      call tetrahedron_quality1_3d ( tetra_xyz, tetra_quality(tetra) )

    end do
  end subroutine tet_mesh_quality1

  subroutine tet_mesh_quality2 ( node_num, node_xyz, tetra_order, tetra_num, &
    tetra_node, tetra_quality ) &
        bind(C, name="tet_mesh_quality2")

  !*****************************************************************************80
  !
  !! TET_MESH_QUALITY2 returns the quality of each tet in a mesh.
  !
  !  Discussion:
  !
  !    The overall tet mesh quality measure is the minimum of the
  !    corresponding tetrahedron quality measure, over all tetrahedrons in the
  !    tet mesh.
  !
  !    Although tetrahedrons of order 10 are allowed as input,
  !    only the first 4 nodes (presumed to be the vertices) are used
  !    in the calculation.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 November 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) NODE_NUM, the number of nodes.
  !
  !    Input, real(dp) NODE_XYZ(3,NODE_NUM), the nodes.
  !
  !    Input, integer(ip) TETRA_ORDER, the order of the mesh, either 4 or 10.
  !
  !    Input, integer(ip) TETRA_NUM, the number of tetrahedrons.
  !
  !    Input, integer(ip) TETRA_NODE(TETRA_ORDER,TETRA_NUM), the nodes
  !    that make up the tetrahedrons.
  !
  !    Output, real(dp) TETRA_QUALITY(TETRA_NUM), the quality
  !    measure for each tetrahedron.
  !

    integer(ip), parameter :: dim_num = 3
    integer(ip), intent(in), value :: node_num
    integer(ip), intent(in), value :: tetra_num
    integer(ip), intent(in), value :: tetra_order

    real(dp), intent(in) :: node_xyz(dim_num,node_num)
    integer(ip) :: tetra
    integer(ip), intent(in) :: tetra_node(tetra_order,tetra_num)
    real(dp), intent(out) :: tetra_quality(tetra_num)
    real(dp) :: tetra_xyz(dim_num,4)

    do tetra = 1, tetra_num

      tetra_xyz(1:dim_num,1:4) = node_xyz(1:dim_num,tetra_node(1:4,tetra))

      call tetrahedron_quality2_3d ( tetra_xyz, tetra_quality(tetra) )

    end do
  end subroutine tet_mesh_quality2

  subroutine tet_mesh_quality3 ( node_num, node_xyz, tetra_order, tetra_num, &
    tetra_node, tetra_quality ) &
        bind(C, name="tet_mesh_quality3")

  !*****************************************************************************80
  !
  !! TET_MESH_QUALITY3 returns the quality of each tet in a mesh.
  !
  !  Discussion:
  !
  !    The overall tet mesh quality measure is the minimum of the
  !    corresponding tetrahedron quality measure, over all tetrahedrons in the
  !    tet mesh.
  !
  !    Although tetrahedrons of order 10 are allowed as input,
  !    only the first 4 nodes (presumed to be the vertices) are used
  !    in the calculation.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 November 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) NODE_NUM, the number of nodes.
  !
  !    Input, real(dp) NODE_XYZ(3,NODE_NUM), the nodes.
  !
  !    Input, integer(ip) TETRA_ORDER, the order of the mesh, either 4 or 10.
  !
  !    Input, integer(ip) TETRA_NUM, the number of tetrahedrons.
  !
  !    Input, integer(ip) TETRA_NODE(TETRA_ORDER,TETRA_NUM), the nodes
  !    that make up the tetrahedrons.
  !
  !    Output, real(dp) TETRA_QUALITY(TETRA_NUM), the quality
  !    measure for each tetrahedron.
  !

    integer(ip), parameter :: dim_num = 3
    integer(ip), intent(in), value :: node_num
    integer(ip), intent(in), value :: tetra_num
    integer(ip), intent(in), value :: tetra_order

    real(dp), intent(in) :: node_xyz(dim_num,node_num)
    integer(ip) :: tetra
    integer(ip), intent(in) :: tetra_node(tetra_order,tetra_num)
    real(dp), intent(out) :: tetra_quality(tetra_num)
    real(dp) :: tetra_xyz(dim_num,4)

    do tetra = 1, tetra_num

      tetra_xyz(1:dim_num,1:4) = node_xyz(1:dim_num,tetra_node(1:4,tetra))

      call tetrahedron_quality3_3d ( tetra_xyz, tetra_quality(tetra) )

    end do
  end subroutine tet_mesh_quality3

  subroutine tet_mesh_quality4 ( node_num, node_xyz, tetra_order, tetra_num, &
    tetra_node, tetra_quality ) &
        bind(C, name="tet_mesh_quality4")

  !*****************************************************************************80
  !
  !! TET_MESH_QUALITY4 returns the quality of each tet in a mesh.
  !
  !  Discussion:
  !
  !    The overall tet mesh quality measure is the minimum of the
  !    corresponding tetrahedron quality measure, over all tetrahedrons in the
  !    tet mesh.
  !
  !    Although tetrahedrons of order 10 are allowed as input,
  !    only the first 4 nodes (presumed to be the vertices) are used
  !    in the calculation.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 November 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) NODE_NUM, the number of nodes.
  !
  !    Input, real(dp) NODE_XYZ(3,NODE_NUM), the nodes.
  !
  !    Input, integer(ip) TETRA_ORDER, the order of the mesh, either 4 or 10.
  !
  !    Input, integer(ip) TETRA_NUM, the number of tetrahedrons.
  !
  !    Input, integer(ip) TETRA_NODE(TETRA_ORDER,TETRA_NUM), the nodes
  !    that make up the tetrahedrons.
  !
  !    Output, real(dp) TETRA_QUALITY(TETRA_NUM), the quality
  !    measure for each tetrahedron.
  !

    integer(ip), parameter :: dim_num = 3
    integer(ip), intent(in), value :: node_num
    integer(ip), intent(in), value :: tetra_num
    integer(ip), intent(in), value :: tetra_order

    real(dp), intent(in) :: node_xyz(dim_num,node_num)
    integer(ip) :: tetra
    integer(ip), intent(in) :: tetra_node(tetra_order,tetra_num)
    real(dp), intent(out) :: tetra_quality(tetra_num)
    real(dp) :: tetra_xyz(dim_num,4)

    do tetra = 1, tetra_num

      tetra_xyz(1:dim_num,1:4) = node_xyz(1:dim_num,tetra_node(1:4,tetra))

      call tetrahedron_quality4_3d ( tetra_xyz, tetra_quality(tetra) )

    end do
  end subroutine tet_mesh_quality4

  subroutine tet_mesh_quality5 ( node_num, node_xyz, tetra_order, tetra_num, &
    tetra_node, tetra_quality ) &
        bind(C, name="tet_mesh_quality5")

  !*****************************************************************************80
  !
  !! TET_MESH_QUALITY5 returns the quality of each tet in a mesh.
  !
  !  Discussion:
  !
  !    The overall tet mesh quality measure is the ratio of the minimum
  !    tetrahedron volume to the maximum tetrahedron volume.
  !
  !    Although tetrahedrons of order 10 are allowed as input,
  !    only the first 4 nodes (presumed to be the vertices) are used
  !    in the calculation.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 November 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) NODE_NUM, the number of nodes.
  !
  !    Input, real(dp) NODE_XYZ(3,NODE_NUM), the nodes.
  !
  !    Input, integer(ip) TETRA_ORDER, the order of the mesh, either 4 or 10.
  !
  !    Input, integer(ip) TETRA_NUM, the number of tetrahedrons.
  !
  !    Input, integer(ip) TETRA_NODE(TETRA_ORDER,TETRA_NUM), the nodes
  !    that make up the tetrahedrons.
  !
  !    Output, real(dp) TETRA_QUALITY(TETRA_NUM), the quality
  !    measure for each tetrahedron.
  !

    integer(ip), parameter :: dim_num = 3
    integer(ip), intent(in), value :: node_num
    integer(ip), intent(in), value :: tetra_num
    integer(ip), intent(in), value :: tetra_order

    real(dp), intent(in) :: node_xyz(dim_num,node_num)
    integer(ip) :: tetra
    integer(ip), intent(in) :: tetra_node(tetra_order,tetra_num)
    real(dp), intent(out) :: tetra_quality(tetra_num)
    real(dp) :: tetra_xyz(dim_num,4)
    real(dp) :: volume_max

    do tetra = 1, tetra_num

      tetra_xyz(1:dim_num,1:4) = node_xyz(1:dim_num,tetra_node(1:4,tetra))

      call tetrahedron_volume_3d ( tetra_xyz, tetra_quality(tetra) )

    end do

    volume_max = maxval ( tetra_quality(1:tetra_num) )

    if ( 0.0_dp < volume_max ) then
      tetra_quality(1:tetra_num) = tetra_quality(1:tetra_num) / volume_max
    end if
  end subroutine tet_mesh_quality5

  subroutine tetrahedron_circumsphere_3d ( tetra, r, pc ) &
        bind(C, name="tetrahedron_circumsphere_3d")

  !*****************************************************************************80
  !
  !! TETRAHEDRON_CIRCUMSPHERE_3D computes the circumsphere of a tetrahedron in 3D.
  !
  !  Discussion:
  !
  !    The circumsphere, or circumscribed sphere, of a tetrahedron is the
  !    sphere that passes through the four vertices.  The circumsphere is
  !    not necessarily the smallest sphere that contains the tetrahedron.
  !
  !    Surprisingly, the diameter of the sphere can be found by solving
  !    a 3 by 3 linear system.  This is because the vectors P2 - P1,
  !    P3 - P1 and P4 - P1 are secants of the sphere, and each forms a
  !    right triangle with the diameter through P1.  Hence, the dot product of
  !    P2 - P1 with that diameter is equal to the square of the length
  !    of P2 - P1, and similarly for P3 - P1 and P4 - P1.  This determines
  !    the diameter vector originating at P1, and hence the radius and
  !    center.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    10 August 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Adrian Bowyer and John Woodwark,
  !    A Programmer's Geometry,
  !    Butterworths, 1983.
  !
  !  Parameters:
  !
  !    Input, real(dp) TETRA(3,4) the tetrahedron vertices.
  !
  !    Output, real(dp) R, PC(3), the center of the
  !    circumscribed sphere, and its radius.  If the linear system is
  !    singular, then R = -1, PC(1:3) = 0.
  !

    integer(ip), parameter :: dim_num = 3
    integer(ip), parameter :: rhs_num = 1

    real(dp) :: a(dim_num,dim_num+rhs_num)
    integer(ip) :: i
    integer(ip) :: info
    integer(ip) :: j
    real(dp), intent(out) :: pc(dim_num)
    real(dp), intent(out) :: r
    real(dp), intent(in) :: tetra(dim_num,4)
  !
  !  Set up the linear system.
  !
    a(1:dim_num,1:3) = transpose ( tetra(1:dim_num,2:4) )

    do j = 1, dim_num
      a(1:dim_num,j) = a(1:dim_num,j) - tetra(j,1)
    end do

    do i = 1, 3
      a(i,4) = sum ( a(i,1:3)**2 )
    end do
  !
  !  Solve the linear system.
  !
    call r8mat_solve ( dim_num, rhs_num, a, info )
  !
  !  If the system was singular, return a consolation prize.
  !
    if ( info /= 0 ) then
      r = -1.0_dp
      pc(1:dim_num) = 0.0_dp
    end if
  !
  !  Compute the radius and center.
  !
    r = 0.5_dp * sqrt ( sum ( a(1:dim_num,4)**2 ) )

    pc(1:dim_num) = tetra(1:dim_num,1) + 0.5_dp * a(1:dim_num,4)
  end subroutine tetrahedron_circumsphere_3d

  pure subroutine tetrahedron_edge_length_3d ( tetra, edge_length ) &
        bind(C, name="tetrahedron_edge_length_3d")

  !*****************************************************************************80
  !
  !! TETRAHEDRON_EDGE_LENGTH_3D returns edge lengths of a tetrahedron in 3D.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    09 August 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real(dp) TETRA(3,4), the tetrahedron vertices.
  !
  !    Output, real(dp) EDGE_LENGTH(6), the length of the edges.
  !

    integer(ip), parameter :: dim_num = 3

    real(dp) :: r8vec_length
    real(dp), intent(out) :: edge_length(6)
    integer(ip) :: j1
    integer(ip) :: j2
    integer(ip) :: k
    real(dp), intent(in) :: tetra(dim_num,4)

    k = 0
    do j1 = 1, 3
      do j2 = j1+1, 4
        k = k + 1
        edge_length(k) = r8vec_length ( dim_num, &
          tetra(1:dim_num,j2) - tetra(1:dim_num,j1) )
      end do
    end do
  end subroutine tetrahedron_edge_length_3d

  subroutine tetrahedron_insphere_3d ( tetra, r, pc ) &
        bind(C, name="tetrahedron_insphere_3d")

  !*****************************************************************************80
  !
  !! TETRAHEDRON_INSPHERE_3D finds the insphere of a tetrahedron in 3D.
  !
  !  Discussion:
  !
  !    The insphere of a tetrahedron is the inscribed sphere, which touches
  !    each face of the tetrahedron at a single point.
  !
  !    The points of contact are the centroids of the triangular faces
  !    of the tetrahedron.  Therefore, the point of contact for a face
  !    can be computed as the average of the vertices of that face.
  !
  !    The sphere can then be determined as the unique sphere through
  !    the four given centroids.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    08 August 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    David Eberly,
  !    Centers of a Simplex,
  !    http://www.geometrictools.com
  !
  !  Parameters:
  !
  !    Input, real(dp) TETRA(3,4), the vertices of the tetrahedron.
  !
  !    Output, real(dp) R, PC(3), the radius and the center
  !    of the sphere.
  !

    integer(ip), parameter :: dim_num = 3

    real(dp) :: b(4,4)
    real(dp) :: r8mat_det_4d
    real(dp) :: r8vec_length
    real(dp) :: gamma
    real(dp) :: l123
    real(dp) :: l124
    real(dp) :: l134
    real(dp) :: l234
    real(dp) :: n123(1:dim_num)
    real(dp) :: n124(1:dim_num)
    real(dp) :: n134(1:dim_num)
    real(dp) :: n234(1:dim_num)
    real(dp), intent(out) :: pc(1:dim_num)
    real(dp), intent(out) :: r
    real(dp), intent(in) :: tetra(1:dim_num,4)
    real(dp) :: v21(1:dim_num)
    real(dp) :: v31(1:dim_num)
    real(dp) :: v41(1:dim_num)
    real(dp) :: v32(1:dim_num)
    real(dp) :: v42(1:dim_num)
    real(dp) :: v43(1:dim_num)

    v21(1:dim_num) = tetra(1:dim_num,2) - tetra(1:dim_num,1)
    v31(1:dim_num) = tetra(1:dim_num,3) - tetra(1:dim_num,1)
    v41(1:dim_num) = tetra(1:dim_num,4) - tetra(1:dim_num,1)
    v32(1:dim_num) = tetra(1:dim_num,3) - tetra(1:dim_num,2)
    v42(1:dim_num) = tetra(1:dim_num,4) - tetra(1:dim_num,2)
    v43(1:dim_num) = tetra(1:dim_num,4) - tetra(1:dim_num,3)

    call r8vec_cross_3d ( v21, v31, n123 )
    call r8vec_cross_3d ( v41, v21, n124 )
    call r8vec_cross_3d ( v31, v41, n134 )
    call r8vec_cross_3d ( v42, v32, n234 )

    l123 = r8vec_length ( dim_num, n123 )
    l124 = r8vec_length ( dim_num, n124 )
    l134 = r8vec_length ( dim_num, n134 )
    l234 = r8vec_length ( dim_num, n234 )

    pc(1:dim_num) = ( l234 * tetra(1:dim_num,1)   &
                    + l134 * tetra(1:dim_num,2)   &
                    + l124 * tetra(1:dim_num,3)   &
                    + l123 * tetra(1:dim_num,4) ) &
                  / ( l234 + l134 + l124 + l123 )

    b(1:dim_num,1:4) = tetra(1:dim_num,1:4)
    b(4,1:4) = 1.0_dp

    gamma = abs ( r8mat_det_4d ( b ) )

  ! gamma = abs ( &
  !     ( tetra(1,2) * tetra(2,3) * tetra(3,4) &
  !     - tetra(1,3) * tetra(2,4) * tetra(3,2) &
  !     + tetra(1,4) * tetra(2,2) * tetra(3,3) ) &
  !   - ( tetra(1,1) * tetra(2,3) * tetra(3,4) &
  !     - tetra(1,3) * tetra(2,4) * tetra(3,1) &
  !     + tetra(1,4) * tetra(2,1) * tetra(3,3) ) &
  !   + ( tetra(1,1) * tetra(2,2) * tetra(3,4) &
  !     - tetra(1,2) * tetra(2,4) * tetra(3,1) &
  !     + tetra(1,4) * tetra(2,1) * tetra(3,2) ) &
  !   - ( tetra(1,1) * tetra(2,2) * tetra(3,3) &
  !     - tetra(1,2) * tetra(2,3) * tetra(3,1) &
  !     + tetra(1,3) * tetra(2,1) * tetra(3,2) ) )

    r = gamma / ( l234 + l134 + l124 + l123 )
  end subroutine tetrahedron_insphere_3d

  subroutine tetrahedron_quality1_3d ( tetra, quality ) &
        bind(C, name="tetrahedron_quality1_3d")

  !*****************************************************************************80
  !
  !! TETRAHEDRON_QUALITY1_3D: "quality" of a tetrahedron in 3D.
  !
  !  Discussion:
  !
  !    The quality of a tetrahedron is 3.0 times the ratio of the radius of
  !    the inscribed sphere divided by that of the circumscribed sphere.
  !
  !    An equilateral tetrahredron achieves the maximum possible quality of 1.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    09 August 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real(dp) TETRA(3,4), the tetrahedron vertices.
  !
  !    Output, real(dp) QUALITY, the quality of the tetrahedron.
  !

    integer(ip), parameter :: dim_num = 3

    real(dp) :: pc(dim_num)
    real(dp), intent(out) :: quality
    real(dp) :: r_in
    real(dp) :: r_out
    real(dp), intent(in) :: tetra(dim_num,4)

    call tetrahedron_circumsphere_3d ( tetra, r_out, pc )

    call tetrahedron_insphere_3d ( tetra, r_in, pc )

    quality = 3.0_dp * r_in / r_out
  end subroutine tetrahedron_quality1_3d

  subroutine tetrahedron_quality2_3d ( tetra, quality2 ) &
        bind(C, name="tetrahedron_quality2_3d")

  !*****************************************************************************80
  !
  !! TETRAHEDRON_QUALITY2_3D: "quality" of a tetrahedron in 3D.
  !
  !  Discussion:
  !
  !    The quality measure #2 of a tetrahedron is:
  !
  !      QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
  !
  !    where
  !
  !      RIN = radius of the inscribed sphere;
  !      LMAX = length of longest side of the tetrahedron.
  !
  !    An equilateral tetrahredron achieves the maximum possible quality of 1.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    16 August 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Qiang Du and Desheng Wang,
  !    The Optimal Centroidal Voronoi Tesselations and the Gersho's
  !    Conjecture in the Three-Dimensional Space,
  !    Computers and Mathematics with Applications,
  !    Volume 49, 2005, pages 1355-1373.
  !
  !  Parameters:
  !
  !    Input, real(dp) TETRA(3,4), the tetrahedron vertices.
  !
  !    Output, real(dp) QUALITY2, the quality of the tetrahedron.
  !

    integer(ip), parameter :: dim_num = 3

    real(dp) :: edge_length(6)
    real(dp) :: l_max
    real(dp) :: pc(dim_num)
    real(dp), intent(out) :: quality2
    real(dp) :: r_in
    real(dp), intent(in) :: tetra(dim_num,4)

    call tetrahedron_edge_length_3d ( tetra, edge_length )

    l_max = maxval ( edge_length(1:6) )

    call tetrahedron_insphere_3d ( tetra, r_in, pc )

    quality2 = 2.0_dp * sqrt ( 6.0_dp ) * r_in / l_max
  end subroutine tetrahedron_quality2_3d

  pure subroutine tetrahedron_quality3_3d ( tetra, quality3 ) &
        bind(C, name="tetrahedron_quality3_3d")

  !******************************************************************************
  !
  !! TETRAHEDRON_QUALITY3_3D computes the mean ratio of a tetrahedron.
  !
  !  Discussion:
  !
  !    This routine computes QUALITY3, the eigenvalue or mean ratio of
  !    a tetrahedron.
  !
  !      QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of squares of edge lengths).
  !
  !    This value may be used as a shape quality measure for the tetrahedron.
  !
  !    For an equilateral tetrahedron, the value of this quality measure
  !    will be 1.  For any other tetrahedron, the value will be between
  !    0 and 1.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 August 2005
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
  !    Input, real(dp) TETRA(3,4), the vertices of the tetrahedron.
  !
  !    Output, real(dp) QUALITY3, the mean ratio of the tetrahedron.
  !

    integer(ip), parameter :: dim_num = 3

    real(dp) :: ab(dim_num)
    real(dp) :: ac(dim_num)
    real(dp) :: ad(dim_num)
    real(dp) :: bc(dim_num)
    real(dp) :: bd(dim_num)
    real(dp) :: cd(dim_num)
    real(dp) :: denom
    real(dp) :: lab
    real(dp) :: lac
    real(dp) :: lad
    real(dp) :: lbc
    real(dp) :: lbd
    real(dp) :: lcd
    real(dp), intent(out) :: quality3
    real(dp), intent(in) :: tetra(dim_num,4)
    real(dp) :: volume
  !
  !  Compute the vectors representing the sides of the tetrahedron.
  !
    ab(1:3) = tetra(1:dim_num,2) - tetra(1:dim_num,1)
    ac(1:3) = tetra(1:dim_num,3) - tetra(1:dim_num,1)
    ad(1:3) = tetra(1:dim_num,4) - tetra(1:dim_num,1)
    bc(1:3) = tetra(1:dim_num,3) - tetra(1:dim_num,2)
    bd(1:3) = tetra(1:dim_num,4) - tetra(1:dim_num,2)
    cd(1:3) = tetra(1:dim_num,4) - tetra(1:dim_num,3)
  !
  !  Compute the squares of the lengths of the sides.
  !
    lab = sum ( ab(1:dim_num)**2 )
    lac = sum ( ac(1:dim_num)**2 )
    lad = sum ( ad(1:dim_num)**2 )
    lbc = sum ( bc(1:dim_num)**2 )
    lbd = sum ( bd(1:dim_num)**2 )
    lcd = sum ( cd(1:dim_num)**2 )
  !
  !  Compute the volume.
  !
    volume = abs ( &
        ab(1) * ( ac(2) * ad(3) - ac(3) * ad(2) ) &
      + ab(2) * ( ac(3) * ad(1) - ac(1) * ad(3) ) &
      + ab(3) * ( ac(1) * ad(2) - ac(2) * ad(1) ) ) / 6.0_dp

    denom = lab + lac + lad + lbc + lbd + lcd

    if ( denom == 0.0_dp ) then
      quality3 = 0.0_dp
    else
      quality3 = 12.0_dp * ( 3.0_dp * volume )**( 2.0_dp / 3.0_dp ) / denom
    end if
  end subroutine tetrahedron_quality3_3d

  pure subroutine tetrahedron_quality4_3d ( tetra, quality4 ) &
        bind(C, name="tetrahedron_quality4_3d")

  !******************************************************************************
  !
  !! TETRAHEDRON_QUALITY4_3D computes the minimum solid angle of a tetrahedron.
  !
  !  Discussion:
  !
  !    This routine computes a quality measure for a tetrahedron, based
  !    on the sine of half the minimum of the four solid angles.
  !
  !    The quality measure for an equilateral tetrahedron should be 1,
  !    since the solid angles of such a tetrahedron are each equal to pi.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 August 2005
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
  !    Input, real(dp) TETRA(3,4), the vertices of the tetrahedron.
  !
  !    Output, real(dp) QUALITY4, the value of the quality measure.
  !

    integer(ip), parameter :: dim_num = 3

    real(dp) :: a(dim_num)
    real(dp) :: ab(dim_num)
    real(dp) :: ac(dim_num)
    real(dp) :: ad(dim_num)
    real(dp) :: b(dim_num)
    real(dp) :: bc(dim_num)
    real(dp) :: bd(dim_num)
    real(dp) :: c(dim_num)
    real(dp) :: cd(dim_num)
    real(dp) :: d(dim_num)
    real(dp) :: denom
    real(dp) :: l1
    real(dp) :: l2
    real(dp) :: l3
    real(dp) :: lab
    real(dp) :: lac
    real(dp) :: lad
    real(dp) :: lbc
    real(dp) :: lbd
    real(dp) :: lcd
    real(dp), intent(out) :: quality4
    real(dp), intent(in) :: tetra(dim_num,4)
    real(dp) :: volume
  !
  !  Compute the vectors that represent the sides.
  !
    ab(1:dim_num) = tetra(1:dim_num,2) - tetra(1:dim_num,1)
    ac(1:dim_num) = tetra(1:dim_num,3) - tetra(1:dim_num,1)
    ad(1:dim_num) = tetra(1:dim_num,4) - tetra(1:dim_num,1)
    bc(1:dim_num) = tetra(1:dim_num,3) - tetra(1:dim_num,2)
    bd(1:dim_num) = tetra(1:dim_num,4) - tetra(1:dim_num,2)
    cd(1:dim_num) = tetra(1:dim_num,4) - tetra(1:dim_num,3)
  !
  !  Compute the lengths of the sides.
  !
    lab = sqrt ( sum ( ab(1:dim_num)**2 ) )
    lac = sqrt ( sum ( ac(1:dim_num)**2 ) )
    lad = sqrt ( sum ( ad(1:dim_num)**2 ) )
    lbc = sqrt ( sum ( bc(1:dim_num)**2 ) )
    lbd = sqrt ( sum ( bd(1:dim_num)**2 ) )
    lcd = sqrt ( sum ( cd(1:dim_num)**2 ) )
  !
  !  Compute the volume
  !
    volume = abs ( &
        ab(1) * ( ac(2) * ad(3) - ac(3) * ad(2) ) &
      + ab(2) * ( ac(3) * ad(1) - ac(1) * ad(3) ) &
      + ab(3) * ( ac(1) * ad(2) - ac(2) * ad(1) ) ) / 6.0_dp

    quality4 = 1.0_dp

    l1 = lab + lac
    l2 = lab + lad
    l3 = lac + lad

    denom = ( l1 + lbc ) * ( l1 - lbc ) &
          * ( l2 + lbd ) * ( l2 - lbd ) &
          * ( l3 + lcd ) * ( l3 - lcd )

    if ( denom <= 0.0_dp ) then
      quality4 = 0.0_dp
    else
      quality4 = min ( quality4, 12.0_dp * volume / sqrt ( denom ) )
    end if

    l1 = lab + lbc
    l2 = lab + lbd
    l3 = lbc + lbd

    denom = ( l1 + lac ) * ( l1 - lac ) &
          * ( l2 + lad ) * ( l2 - lad ) &
          * ( l3 + lcd ) * ( l3 - lcd )

    if ( denom <= 0.0_dp ) then
      quality4 = 0.0_dp
    else
      quality4 = min ( quality4, 12.0_dp * volume / sqrt ( denom ) )
    end if

    l1 = lac + lbc
    l2 = lac + lcd
    l3 = lbc + lcd

    denom = ( l1 + lab ) * ( l1 - lab ) &
          * ( l2 + lad ) * ( l2 - lad ) &
          * ( l3 + lbd ) * ( l3 - lbd )

    if ( denom <= 0.0_dp ) then
      quality4 = 0.0_dp
    else
      quality4 = min ( quality4, 12.0_dp * volume / sqrt ( denom ) )
    end if

    l1 = lad + lbd
    l2 = lad + lcd
    l3 = lbd + lcd

    denom = ( l1 + lab ) * ( l1 - lab ) &
          * ( l2 + lac ) * ( l2 - lac ) &
          * ( l3 + lbc ) * ( l3 - lbc )

    if ( denom <= 0.0_dp ) then
      quality4 = 0.0_dp
    else
      quality4 = min ( quality4, 12.0_dp * volume / sqrt ( denom ) )
    end if

    quality4 = quality4 * 1.5_dp * sqrt ( 6.0_dp )
  end subroutine tetrahedron_quality4_3d

  pure subroutine tetrahedron_volume_3d ( tetra, volume ) &
        bind(C, name="tetrahedron_volume_3d")

  !*****************************************************************************80
  !
  !! TETRAHEDRON_VOLUME_3D computes the volume of a tetrahedron in 3D.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    30 December 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real(dp) TETRA(3,4), the vertices of the tetrahedron.
  !
  !    Output, real(dp) VOLUME, the volume of the tetrahedron.
  !

    integer(ip), parameter :: dim_num = 3

    real(dp) :: a(4,4)
    real(dp) :: r8mat_det_4d
    real(dp), intent(in) :: tetra(dim_num,4)
    real(dp), intent(out) :: volume

    a(1:dim_num,1:4) = tetra(1:dim_num,1:4)
    a(4,1:4) = 1.0_dp

    volume = abs ( r8mat_det_4d ( a ) ) / 6.0_dp
  end subroutine tetrahedron_volume_3d

end module tet_mesh_quality_mod
