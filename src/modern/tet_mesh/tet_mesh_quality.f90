!> tet_mesh_quality — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module tet_mesh_quality_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  implicit none
  private

  public :: i4vec_histogram, mesh_base_one, r8_swap, r8mat_det_4d, r8mat_solve, r8vec_cross_3d
  public :: r8vec_length, r8vec_mean, r8vec_variance, tet_mesh_node_order, tet_mesh_quality1, tet_mesh_quality2
  public :: tet_mesh_quality3, tet_mesh_quality4, tet_mesh_quality5, tetrahedron_circumsphere_3d, tetrahedron_edge_length_3d, tetrahedron_insphere_3d
  public :: tetrahedron_quality1_3d, tetrahedron_quality2_3d, tetrahedron_quality3_3d, tetrahedron_quality4_3d, tetrahedron_volume_3d

contains

  subroutine i4vec_histogram ( n, a, histo_num, histo_gram )

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
  !    Input, integer(int32) N, the number of elements of A.
  !
  !    Input, integer(int32) A(N), the array to examine.
  !
  !    Input, integer(int32) HISTO_NUM, the maximum value for which a
  !    histogram entry will be computed.
  !
  !    Output, integer(int32) HISTO_GRAM(0:HISTO_NUM), contains the number of
  !    entries of A with the values of 0 through HISTO_NUM.
  !

    integer(int32) histo_num
    integer(int32) n

    integer(int32) a(n)
    integer(int32) histo_gram(0:histo_num)
    integer(int32) i

    histo_gram(0:histo_num) = 0

    do i = 1, n

      if ( 0 <= a(i) .and. a(i) <= histo_num ) then
        histo_gram(a(i)) = histo_gram(a(i)) + 1
      end if

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
      write ( *, '(a)' )'MESH_BASE_ZERO:'
      write ( *, '(a)' )'  The element indexing appears to be 1-based!'
      write ( *, '(a)' )'  No conversion is necessary.'
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )'MESH_BASE_ZERO - Warning!'
      write ( *, '(a)' )' The element indexing is not of a recognized type.'
    end if
  end

  subroutine r8_swap ( x, y )

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
  !    Input/output, real(real64) X, Y.  On output, the values of X and
  !    Y have been interchanged.
  !

    real(real64) x
    real(real64) y
    real(real64) z

    z = x
    x = y
    y = z
  end

  function r8mat_det_4d ( a )

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
  !    Input, real(real64) A(4,4), the matrix whose determinant is desired.
  !
  !    Output, real(real64) R8MAT_DET_4D, the determinant of the matrix.
  !

    real(real64) a(4,4)
    real(real64) r8mat_det_4d

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
  end

  subroutine r8mat_solve ( n, rhs_num, a, info )

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
  !    Input, integer(int32) N, the order of the matrix.
  !
  !    Input, integer(int32) rhs_num, the number of right hand sides.  rhs_num
  !    must be at least 0.
  !
  !    Input/output, real(real64) A(N,N+rhs_num), contains in rows and
  !    columns 1 to N the coefficient matrix, and in columns N+1 through
  !    N+rhs_num, the right hand sides.  On output, the coefficient matrix
  !    area has been destroyed, while the right hand sides have
  !    been overwritten with the corresponding solutions.
  !
  !    Output, integer(int32) INFO, singularity flag.
  !    0, the matrix was not singular, the solutions were computed;
  !    J, factorization failed on step J, and the solutions could not
  !    be computed.
  !

    integer(int32) n
    integer(int32) rhs_num

    real(real64) a(n,n+rhs_num)
    real(real64) apivot
    real(real64) factor
    integer(int32) i
    integer(int32) info
    integer(int32) ipivot
    integer(int32) j

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

      if ( apivot == 0.0e+00_real64 ) then
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
      a(j,j) = 1.0e+00_real64
      a(j,j+1:n+rhs_num) = a(j,j+1:n+rhs_num) / apivot
  !
  !  A(I,J) becomes 0.
  !
      do i = 1, n

        if ( i /= j ) then

          factor = a(i,j)
          a(i,j) = 0.0e+00_real64
          a(i,j+1:n+rhs_num) = a(i,j+1:n+rhs_num) - factor * a(j,j+1:n+rhs_num)

        end if

      end do

    end do
  end

  subroutine r8vec_cross_3d ( v1, v2, v3 )

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
  !    Input, real(real64) V1(3), V2(3), the two vectors.
  !
  !    Output, real(real64) V3(3), the cross product vector.
  !

    integer(int32), parameter :: dim_num = 3

    real(real64) v1(dim_num)
    real(real64) v2(dim_num)
    real(real64) v3(dim_num)

    v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
    v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
    v3(3) = v1(1) * v2(2) - v1(2) * v2(1)
  end

  function r8vec_length ( dim_num, x )

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
  !    Input, integer(int32) DIM_NUM, the spatial dimension.
  !
  !    Input, real(real64) X(DIM_NUM), the vector.
  !
  !    Output, real(real64) R8VEC_LENGTH, the Euclidean length of the vector.
  !

    integer(int32) dim_num

    real(real64) r8vec_length
    real(real64) x(dim_num)

    r8vec_length = sqrt ( sum ( ( x(1:dim_num) )**2 ) )
  end

  subroutine r8vec_mean ( n, a, mean )

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
  !    Input, integer(int32) N, the number of entries in the vector.
  !
  !    Input, real(real64) A(N), the vector whose mean is desired.
  !
  !    Output, real(real64) MEAN, the mean of the vector entries.
  !

    integer(int32) n

    real(real64) a(n)
    real(real64) mean

    mean = sum ( a(1:n) ) / real ( n, real64)
  end

  subroutine r8vec_variance ( n, a, variance )

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
  !    Input, integer(int32) N, the number of entries in the vector.
  !    N should be at least 2.
  !
  !    Input, real(real64) A(N), the vector.
  !
  !    Output, real(real64) VARIANCE, the variance of the vector.
  !

    integer(int32) n

    real(real64) a(n)
    real(real64) mean
    real(real64) variance

    if ( n < 2 ) then

      variance = 0.0e+00_real64

    else

      mean = sum ( a(1:n) ) / real ( n, real64)

      variance = sum ( ( a(1:n) - mean )**2 )

      variance = variance / real ( n - 1, real64)

    end if
  end

  subroutine tet_mesh_node_order ( tetra_order, tetra_num, tetra_node, &
    node_num, node_order )

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
  !    Input, integer(int32) TETRA_ORDER, the order of the mesh, either 4 or 10.
  !
  !    Input, integer(int32) TETRA_NUM, the number of tetrahedrons.
  !
  !    Input, integer(int32) TETRA_NODE(TETRA_ORDER,TETRA_NUM), the nodes
  !    that make up the tetrahedrons.
  !
  !    Input, integer(int32) NODE_NUM, the number of nodes.
  !
  !    Output, integer(int32) NODE_ORDER(NODE_NUM), the order of each node.
  !

    integer(int32) node_num
    integer(int32) tetra_num
    integer(int32) tetra_order

    integer(int32) i
    integer(int32) node
    integer(int32) node_order(node_num)
    integer(int32) tetra
    integer(int32) tetra_node(tetra_order,tetra_num)

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
  end

  subroutine tet_mesh_quality1 ( node_num, node_xyz, tetra_order, tetra_num, &
    tetra_node, tetra_quality )

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
  !    Input, integer(int32) NODE_NUM, the number of nodes.
  !
  !    Input, real(real64) NODE_XYZ(3,NODE_NUM), the nodes.
  !
  !    Input, integer(int32) TETRA_ORDER, the order of the mesh, either 4 or 10.
  !
  !    Input, integer(int32) TETRA_NUM, the number of tetrahedrons.
  !
  !    Input, integer(int32) TETRA_NODE(TETRA_ORDER,TETRA_NUM), the nodes
  !    that make up the tetrahedrons.
  !
  !    Output, real(real64) TETRA_QUALITY(TETRA_NUM), the quality
  !    measure for each tetrahedron.
  !

    integer(int32), parameter :: dim_num = 3
    integer(int32) node_num
    integer(int32) tetra_num
    integer(int32) tetra_order

    real(real64) node_xyz(dim_num,node_num)
    integer(int32) tetra
    integer(int32) tetra_node(tetra_order,tetra_num)
    real(real64) tetra_quality(tetra_num)
    real(real64) tetra_xyz(dim_num,4)

    do tetra = 1, tetra_num

      tetra_xyz(1:dim_num,1:4) = node_xyz(1:dim_num,tetra_node(1:4,tetra))

      call tetrahedron_quality1_3d ( tetra_xyz, tetra_quality(tetra) )

    end do
  end

  subroutine tet_mesh_quality2 ( node_num, node_xyz, tetra_order, tetra_num, &
    tetra_node, tetra_quality )

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
  !    Input, integer(int32) NODE_NUM, the number of nodes.
  !
  !    Input, real(real64) NODE_XYZ(3,NODE_NUM), the nodes.
  !
  !    Input, integer(int32) TETRA_ORDER, the order of the mesh, either 4 or 10.
  !
  !    Input, integer(int32) TETRA_NUM, the number of tetrahedrons.
  !
  !    Input, integer(int32) TETRA_NODE(TETRA_ORDER,TETRA_NUM), the nodes
  !    that make up the tetrahedrons.
  !
  !    Output, real(real64) TETRA_QUALITY(TETRA_NUM), the quality
  !    measure for each tetrahedron.
  !

    integer(int32), parameter :: dim_num = 3
    integer(int32) node_num
    integer(int32) tetra_num
    integer(int32) tetra_order

    real(real64) node_xyz(dim_num,node_num)
    integer(int32) tetra
    integer(int32) tetra_node(tetra_order,tetra_num)
    real(real64) tetra_quality(tetra_num)
    real(real64) tetra_xyz(dim_num,4)

    do tetra = 1, tetra_num

      tetra_xyz(1:dim_num,1:4) = node_xyz(1:dim_num,tetra_node(1:4,tetra))

      call tetrahedron_quality2_3d ( tetra_xyz, tetra_quality(tetra) )

    end do
  end

  subroutine tet_mesh_quality3 ( node_num, node_xyz, tetra_order, tetra_num, &
    tetra_node, tetra_quality )

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
  !    Input, integer(int32) NODE_NUM, the number of nodes.
  !
  !    Input, real(real64) NODE_XYZ(3,NODE_NUM), the nodes.
  !
  !    Input, integer(int32) TETRA_ORDER, the order of the mesh, either 4 or 10.
  !
  !    Input, integer(int32) TETRA_NUM, the number of tetrahedrons.
  !
  !    Input, integer(int32) TETRA_NODE(TETRA_ORDER,TETRA_NUM), the nodes
  !    that make up the tetrahedrons.
  !
  !    Output, real(real64) TETRA_QUALITY(TETRA_NUM), the quality
  !    measure for each tetrahedron.
  !

    integer(int32), parameter :: dim_num = 3
    integer(int32) node_num
    integer(int32) tetra_num
    integer(int32) tetra_order

    real(real64) node_xyz(dim_num,node_num)
    integer(int32) tetra
    integer(int32) tetra_node(tetra_order,tetra_num)
    real(real64) tetra_quality(tetra_num)
    real(real64) tetra_xyz(dim_num,4)

    do tetra = 1, tetra_num

      tetra_xyz(1:dim_num,1:4) = node_xyz(1:dim_num,tetra_node(1:4,tetra))

      call tetrahedron_quality3_3d ( tetra_xyz, tetra_quality(tetra) )

    end do
  end

  subroutine tet_mesh_quality4 ( node_num, node_xyz, tetra_order, tetra_num, &
    tetra_node, tetra_quality )

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
  !    Input, integer(int32) NODE_NUM, the number of nodes.
  !
  !    Input, real(real64) NODE_XYZ(3,NODE_NUM), the nodes.
  !
  !    Input, integer(int32) TETRA_ORDER, the order of the mesh, either 4 or 10.
  !
  !    Input, integer(int32) TETRA_NUM, the number of tetrahedrons.
  !
  !    Input, integer(int32) TETRA_NODE(TETRA_ORDER,TETRA_NUM), the nodes
  !    that make up the tetrahedrons.
  !
  !    Output, real(real64) TETRA_QUALITY(TETRA_NUM), the quality
  !    measure for each tetrahedron.
  !

    integer(int32), parameter :: dim_num = 3
    integer(int32) node_num
    integer(int32) tetra_num
    integer(int32) tetra_order

    real(real64) node_xyz(dim_num,node_num)
    integer(int32) tetra
    integer(int32) tetra_node(tetra_order,tetra_num)
    real(real64) tetra_quality(tetra_num)
    real(real64) tetra_xyz(dim_num,4)

    do tetra = 1, tetra_num

      tetra_xyz(1:dim_num,1:4) = node_xyz(1:dim_num,tetra_node(1:4,tetra))

      call tetrahedron_quality4_3d ( tetra_xyz, tetra_quality(tetra) )

    end do
  end

  subroutine tet_mesh_quality5 ( node_num, node_xyz, tetra_order, tetra_num, &
    tetra_node, tetra_quality )

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
  !    Input, integer(int32) NODE_NUM, the number of nodes.
  !
  !    Input, real(real64) NODE_XYZ(3,NODE_NUM), the nodes.
  !
  !    Input, integer(int32) TETRA_ORDER, the order of the mesh, either 4 or 10.
  !
  !    Input, integer(int32) TETRA_NUM, the number of tetrahedrons.
  !
  !    Input, integer(int32) TETRA_NODE(TETRA_ORDER,TETRA_NUM), the nodes
  !    that make up the tetrahedrons.
  !
  !    Output, real(real64) TETRA_QUALITY(TETRA_NUM), the quality
  !    measure for each tetrahedron.
  !

    integer(int32), parameter :: dim_num = 3
    integer(int32) node_num
    integer(int32) tetra_num
    integer(int32) tetra_order

    real(real64) node_xyz(dim_num,node_num)
    integer(int32) tetra
    integer(int32) tetra_node(tetra_order,tetra_num)
    real(real64) tetra_quality(tetra_num)
    real(real64) tetra_xyz(dim_num,4)
    real(real64) volume_max

    do tetra = 1, tetra_num

      tetra_xyz(1:dim_num,1:4) = node_xyz(1:dim_num,tetra_node(1:4,tetra))

      call tetrahedron_volume_3d ( tetra_xyz, tetra_quality(tetra) )

    end do

    volume_max = maxval ( tetra_quality(1:tetra_num) )

    if ( 0.0e+00_real64 < volume_max ) then
      tetra_quality(1:tetra_num) = tetra_quality(1:tetra_num) / volume_max
    end if
  end

  subroutine tetrahedron_circumsphere_3d ( tetra, r, pc )

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
  !    Input, real(real64) TETRA(3,4) the tetrahedron vertices.
  !
  !    Output, real(real64) R, PC(3), the center of the
  !    circumscribed sphere, and its radius.  If the linear system is
  !    singular, then R = -1, PC(1:3) = 0.
  !

    integer(int32), parameter :: dim_num = 3
    integer(int32), parameter :: rhs_num = 1

    real(real64) a(dim_num,dim_num+rhs_num)
    integer(int32) i
    integer(int32) info
    integer(int32) j
    real(real64) pc(dim_num)
    real(real64) r
    real(real64) tetra(dim_num,4)
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
      r = -1.0e+00_real64
      pc(1:dim_num) = 0.0e+00_real64
    end if
  !
  !  Compute the radius and center.
  !
    r = 0.5e+00_real64 * sqrt ( sum ( a(1:dim_num,4)**2 ) )

    pc(1:dim_num) = tetra(1:dim_num,1) + 0.5e+00_real64 * a(1:dim_num,4)
  end

  subroutine tetrahedron_edge_length_3d ( tetra, edge_length )

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
  !    Input, real(real64) TETRA(3,4), the tetrahedron vertices.
  !
  !    Output, real(real64) EDGE_LENGTH(6), the length of the edges.
  !

    integer(int32), parameter :: dim_num = 3

    real(real64) r8vec_length
    real(real64) edge_length(6)
    integer(int32) j1
    integer(int32) j2
    integer(int32) k
    real(real64) tetra(dim_num,4)

    k = 0
    do j1 = 1, 3
      do j2 = j1+1, 4
        k = k + 1
        edge_length(k) = r8vec_length ( dim_num, &
          tetra(1:dim_num,j2) - tetra(1:dim_num,j1) )
      end do
    end do
  end

  subroutine tetrahedron_insphere_3d ( tetra, r, pc )

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
  !    Input, real(real64) TETRA(3,4), the vertices of the tetrahedron.
  !
  !    Output, real(real64) R, PC(3), the radius and the center
  !    of the sphere.
  !

    integer(int32), parameter :: dim_num = 3

    real(real64) b(4,4)
    real(real64) r8mat_det_4d
    real(real64) r8vec_length
    real(real64) gamma
    real(real64) l123
    real(real64) l124
    real(real64) l134
    real(real64) l234
    real(real64) n123(1:dim_num)
    real(real64) n124(1:dim_num)
    real(real64) n134(1:dim_num)
    real(real64) n234(1:dim_num)
    real(real64) pc(1:dim_num)
    real(real64) r
    real(real64) tetra(1:dim_num,4)
    real(real64) v21(1:dim_num)
    real(real64) v31(1:dim_num)
    real(real64) v41(1:dim_num)
    real(real64) v32(1:dim_num)
    real(real64) v42(1:dim_num)
    real(real64) v43(1:dim_num)

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
    b(4,1:4) = 1.0e+00_real64

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
  end

  subroutine tetrahedron_quality1_3d ( tetra, quality )

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
  !    Input, real(real64) TETRA(3,4), the tetrahedron vertices.
  !
  !    Output, real(real64) QUALITY, the quality of the tetrahedron.
  !

    integer(int32), parameter :: dim_num = 3

    real(real64) pc(dim_num)
    real(real64) quality
    real(real64) r_in
    real(real64) r_out
    real(real64) tetra(dim_num,4)

    call tetrahedron_circumsphere_3d ( tetra, r_out, pc )

    call tetrahedron_insphere_3d ( tetra, r_in, pc )

    quality = 3.0e+00_real64 * r_in / r_out
  end

  subroutine tetrahedron_quality2_3d ( tetra, quality2 )

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
  !    Input, real(real64) TETRA(3,4), the tetrahedron vertices.
  !
  !    Output, real(real64) QUALITY2, the quality of the tetrahedron.
  !

    integer(int32), parameter :: dim_num = 3

    real(real64) edge_length(6)
    real(real64) l_max
    real(real64) pc(dim_num)
    real(real64) quality2
    real(real64) r_in
    real(real64) tetra(dim_num,4)

    call tetrahedron_edge_length_3d ( tetra, edge_length )

    l_max = maxval ( edge_length(1:6) )

    call tetrahedron_insphere_3d ( tetra, r_in, pc )

    quality2 = 2.0e+00_real64 * sqrt ( 6.0e+00_real64 ) * r_in / l_max
  end

  subroutine tetrahedron_quality3_3d ( tetra, quality3 )

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
  !    Input, real(real64) TETRA(3,4), the vertices of the tetrahedron.
  !
  !    Output, real(real64) QUALITY3, the mean ratio of the tetrahedron.
  !

    integer(int32), parameter :: dim_num = 3

    real(real64) ab(dim_num)
    real(real64) ac(dim_num)
    real(real64) ad(dim_num)
    real(real64) bc(dim_num)
    real(real64) bd(dim_num)
    real(real64) cd(dim_num)
    real(real64) denom
    real(real64) lab
    real(real64) lac
    real(real64) lad
    real(real64) lbc
    real(real64) lbd
    real(real64) lcd
    real(real64) quality3
    real(real64) tetra(dim_num,4)
    real(real64) volume
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
      + ab(3) * ( ac(1) * ad(2) - ac(2) * ad(1) ) ) / 6.0e+00_real64

    denom = lab + lac + lad + lbc + lbd + lcd

    if ( denom == 0.0e+00_real64 ) then
      quality3 = 0.0e+00_real64
    else
      quality3 = 12.0e+00_real64 * ( 3.0e+00_real64 * volume )**( 2.0e+00_real64 / 3.0e+00_real64 ) / denom
    end if
  end

  subroutine tetrahedron_quality4_3d ( tetra, quality4 )

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
  !    Input, real(real64) TETRA(3,4), the vertices of the tetrahedron.
  !
  !    Output, real(real64) QUALITY4, the value of the quality measure.
  !

    integer(int32), parameter :: dim_num = 3

    real(real64) a(dim_num)
    real(real64) ab(dim_num)
    real(real64) ac(dim_num)
    real(real64) ad(dim_num)
    real(real64) b(dim_num)
    real(real64) bc(dim_num)
    real(real64) bd(dim_num)
    real(real64) c(dim_num)
    real(real64) cd(dim_num)
    real(real64) d(dim_num)
    real(real64) denom
    real(real64) l1
    real(real64) l2
    real(real64) l3
    real(real64) lab
    real(real64) lac
    real(real64) lad
    real(real64) lbc
    real(real64) lbd
    real(real64) lcd
    real(real64) quality4
    real(real64) tetra(dim_num,4)
    real(real64) volume
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
      + ab(3) * ( ac(1) * ad(2) - ac(2) * ad(1) ) ) / 6.0e+00_real64

    quality4 = 1.0e+00_real64

    l1 = lab + lac
    l2 = lab + lad
    l3 = lac + lad

    denom = ( l1 + lbc ) * ( l1 - lbc ) &
          * ( l2 + lbd ) * ( l2 - lbd ) &
          * ( l3 + lcd ) * ( l3 - lcd )

    if ( denom <= 0.0e+00_real64 ) then
      quality4 = 0.0e+00_real64
    else
      quality4 = min ( quality4, 12.0e+00_real64 * volume / sqrt ( denom ) )
    end if

    l1 = lab + lbc
    l2 = lab + lbd
    l3 = lbc + lbd

    denom = ( l1 + lac ) * ( l1 - lac ) &
          * ( l2 + lad ) * ( l2 - lad ) &
          * ( l3 + lcd ) * ( l3 - lcd )

    if ( denom <= 0.0e+00_real64 ) then
      quality4 = 0.0e+00_real64
    else
      quality4 = min ( quality4, 12.0e+00_real64 * volume / sqrt ( denom ) )
    end if

    l1 = lac + lbc
    l2 = lac + lcd
    l3 = lbc + lcd

    denom = ( l1 + lab ) * ( l1 - lab ) &
          * ( l2 + lad ) * ( l2 - lad ) &
          * ( l3 + lbd ) * ( l3 - lbd )

    if ( denom <= 0.0e+00_real64 ) then
      quality4 = 0.0e+00_real64
    else
      quality4 = min ( quality4, 12.0e+00_real64 * volume / sqrt ( denom ) )
    end if

    l1 = lad + lbd
    l2 = lad + lcd
    l3 = lbd + lcd

    denom = ( l1 + lab ) * ( l1 - lab ) &
          * ( l2 + lac ) * ( l2 - lac ) &
          * ( l3 + lbc ) * ( l3 - lbc )

    if ( denom <= 0.0e+00_real64 ) then
      quality4 = 0.0e+00_real64
    else
      quality4 = min ( quality4, 12.0e+00_real64 * volume / sqrt ( denom ) )
    end if

    quality4 = quality4 * 1.5e+00_real64 * sqrt ( 6.0e+00_real64 )
  end

  subroutine tetrahedron_volume_3d ( tetra, volume )

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
  !    Input, real(real64) TETRA(3,4), the vertices of the tetrahedron.
  !
  !    Output, real(real64) VOLUME, the volume of the tetrahedron.
  !

    integer(int32), parameter :: dim_num = 3

    real(real64) a(4,4)
    real(real64) r8mat_det_4d
    real(real64) tetra(dim_num,4)
    real(real64) volume

    a(1:dim_num,1:4) = tetra(1:dim_num,1:4)
    a(4,1:4) = 1.0e+00_real64

    volume = abs ( r8mat_det_4d ( a ) ) / 6.0e+00_real64
  end

end module tet_mesh_quality_mod
