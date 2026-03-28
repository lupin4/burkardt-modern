!> tri_surface_io — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module tri_surface_io_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  implicit none
  private

  public :: tri_surface_print, tri_surface_read, tri_surface_size, tri_surface_size_print, tri_surface_write

contains

  subroutine tri_surface_print ( node_file_name, triangle_file_name, dim_num, &
    node_num, order_num, triangle_num, node_xyz, triangle_node )

  !*****************************************************************************80
  !
  !! TRI_SURFACE_PRINT prints graphics information from TRI_SURFACE files.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 September 2008
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) NODE_FILE_NAME, the name of the node file.
  !
  !    Input, character ( len = * ) TRIANGLE_FILE_NAME, the name of the
  !    triangle file.
  !
  !    Input, integer(int32) DIM_NUM, the spatial dimension.
  !
  !    Input, integer(int32) NODE_NUM, the number of points.
  !
  !    Input, integer(int32) ORDER_NUM, the order of the triangles.
  !
  !    Input, integer(int32) TRIANGLE_NUM, the number of triangles.
  !
  !    Input, real(real64) NODE_XYZ(DIM_NUM,NODE_NUM), the node coordinates.
  !
  !    Input, integer(int32) TRIANGLE_NODE(ORDER_NUM,TRIANGLE_NUM),
  !    the nodes that form the triangles.
  !

    integer(int32) dim_num
    integer(int32) node_num
    integer(int32) order_num
    integer(int32) triangle_num

    character ( len = *  ) node_file_name
    real(real64) node_xyz(dim_num,node_num)
    character ( len = *  ) triangle_file_name
    integer(int32) triangle_node(order_num,triangle_num)

    call r8mat_transpose_print ( dim_num, node_num, node_xyz, &
      '  Node coordinates' )

    call i4mat_transpose_print ( order_num, triangle_num, triangle_node, &
      '  Triangle nodes' )
  end

  subroutine tri_surface_read ( node_file_name, triangle_file_name, dim_num, &
    node_num, order_num, triangle_num, node_xyz, triangle_node )

  !*****************************************************************************80
  !
  !! TRI_SURFACE_READ reads graphics information from a pair of TRI_SURFACE files.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 September 2008
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) NODE_FILE_NAME, the name of the node file.
  !
  !    Input, character ( len = * ) TRIANGLE_FILE_NAME, the name of the 
  !    triangle file.
  !
  !    Input, integer(int32) DIM_NUM, the spatial dimension.
  !
  !    Input, integer(int32) NODE_NUM, the number of points.
  !
  !    Input, integer(int32) ORDER_NUM, the order of the triangles.
  !
  !    Input, integer(int32) TRIANGLE_NUM, the number of triangles.
  !
  !    Output, real(real64) NODE_XYZ(DIM_NUM,NODE_NUM), the node coordinates.
  !
  !    Output, integer(int32) TRIANGLE_NODE(ORDER_NUM,TRIANGLE_NUM),
  !    the nodes that form the triangles.
  !

    integer(int32) dim_num
    integer(int32) node_num
    integer(int32) order_num
    integer(int32) triangle_num

    character ( len = * ) node_file_name
    real(real64) node_xyz(dim_num,node_num)
    character ( len = * ) triangle_file_name
    integer(int32) triangle_node(order_num,triangle_num)

    call r8mat_data_read ( node_file_name, dim_num, node_num, node_xyz )

    call i4mat_data_read ( triangle_file_name, order_num, triangle_num, &
      triangle_node )
  end

  subroutine tri_surface_size ( node_file_name, triangle_file_name, dim_num, &
    node_num, order_num, triangle_num )

  !*****************************************************************************80
  !
  !! TRI_SURFACE_SIZE determines the size of a TRI_SURFACE object.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 September 2008
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) NODE_FILE_NAME, the name of the node file.
  !
  !    Input, character ( len = * ) TRIANGLE_FILE_NAME, the name of the 
  !    triangle file.
  !
  !    Output, integer(int32) DIM_NUM, the spatial dimension.
  !
  !    Output, integer(int32) NODE_NUM, the number of nodes.
  !
  !    Output, integer(int32) ORDER_NUM, the order of the triangles.
  !
  !    Output, integer(int32) TRIANGLE_NUM, the number of triangles.
  !

    integer(int32) dim_num
    character ( len = * ) node_file_name
    integer(int32) node_num
    integer(int32) order_num
    character ( len = * ) triangle_file_name
    integer(int32) triangle_num

    dim_num = -1
    node_num = -1
    triangle_num = -1

    call r8mat_header_read ( node_file_name, dim_num, node_num )

    if ( dim_num < 2 .or. 3 < dim_num ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TRI_SURFACE_SIZE - Warning!'
      write ( *, '(a,i8)' ) '  The spatial dimension DIM_NUM = ', dim_num
      write ( *, '(a)' ) '  This seems an unlikely value.'
    end if

    call i4mat_header_read ( triangle_file_name, order_num, triangle_num )

    if ( order_num /= 3 .and. order_num /= 6 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TRI_SURFACE_SIZE - Fatal error!'
      write ( *, '(a,i8)' ) '  The "order" of the triangles seems to be ', &
        order_num
      write ( *, '(a)' ) '  Only the values 3 and 6 are acceptable.'
      stop
    end if
  end

  subroutine tri_surface_size_print ( node_file_name, triangle_file_name, &
    dim_num, node_num, order_num, triangle_num )

  !*****************************************************************************80
  !
  !! TRI_SURFACE_SIZE_PRINT prints sizes associated with a TRI_SURFACE file.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 September 2008
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) NODE_FILE_NAME, the name of the node file.
  !
  !    Input, character ( len = * ) TRIANGLE_FILE_NAME, the name of the 
  !    triangle file.
  !
  !    Input, integer(int32) DIM_NUM, the number of spatial dimensions.
  !
  !    Input, integer(int32) NODE_NUM, the number of nodes.
  !
  !    Input, integer(int32) ORDER_NUM, the order of the triangles.
  !
  !    Input, integer(int32) TRIANGLE_NUM, the number of triangles.
  !

    integer(int32) dim_num
    character ( len = * ) node_file_name
    integer(int32) node_num
    integer(int32) order_num
    character ( len = * ) triangle_file_name
    integer(int32) triangle_num

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRI_SURFACE_SIZE_PRINT:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Node file     "' // trim ( node_file_name ) // '".'
    write ( *, '(a)' ) '  Triangle file "' // trim ( triangle_file_name ) // '".'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Spatial dimension  = ', dim_num
    write ( *, '(a,i8)' ) '  Nodes              = ', node_num
    write ( *, '(a,i8)' ) '  Triangle order     = ', order_num
    write ( *, '(a,i8)' ) '  Triangles          = ', triangle_num
  end

  subroutine tri_surface_write ( node_file_name, triangle_file_name, dim_num, &
    node_num, order_num, triangle_num, node_xyz, triangle_node )

  !*****************************************************************************80
  !
  !! TRI_SURFACE_WRITE writes graphics information to a pair of TRI_SURFACE files.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 September 2008
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) NODE_FILE_NAME, the name of the node file.
  !
  !    Input, character ( len = * ) TRIANGLE_FILE_NAME, the name of the 
  !    triangle file.
  !
  !    Input, integer(int32) DIM_NUM, the spatial dimension.
  !
  !    Input, integer(int32) NODE_NUM, the number of points.
  !
  !    Input, integer(int32) ORDER_NUM, the order of the triangles.
  !
  !    Input, integer(int32) TRIANGLE_NUM, the number of triangles.
  !
  !    Input, real(real64) NODE_XYZ(DIM_NUM,NODE_NUM), the node coordinates.
  !
  !    Input, integer(int32) TRIANGLE_NODE(ORDER_NUM,TRIANGLE_NUM),
  !    the nodes that form the triangles.
  !

    integer(int32) dim_num
    integer(int32) node_num
    integer(int32) order_num
    integer(int32) triangle_num

    character ( len = * ) node_file_name
    real(real64) node_xyz(dim_num,node_num)
    character ( len = * ) triangle_file_name
    integer(int32) triangle_node(order_num,triangle_num)

    call r8mat_write ( node_file_name, dim_num, node_num, node_xyz )

    call i4mat_write ( triangle_file_name, order_num, triangle_num, &
      triangle_node )
  end

end module tri_surface_io_mod
