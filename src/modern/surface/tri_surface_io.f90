!> tri_surface_io — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module tri_surface_io_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: tri_surface_print, tri_surface_read, tri_surface_size, tri_surface_size_print, tri_surface_write

contains

  subroutine tri_surface_print ( node_file_name, triangle_file_name, dim_num, &
    node_num, order_num, triangle_num, node_xyz, triangle_node ) &
        bind(C, name="tri_surface_print")

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
  !    Input, integer(ip) DIM_NUM, the spatial dimension.
  !
  !    Input, integer(ip) NODE_NUM, the number of points.
  !
  !    Input, integer(ip) ORDER_NUM, the order of the triangles.
  !
  !    Input, integer(ip) TRIANGLE_NUM, the number of triangles.
  !
  !    Input, real(dp) NODE_XYZ(DIM_NUM,NODE_NUM), the node coordinates.
  !
  !    Input, integer(ip) TRIANGLE_NODE(ORDER_NUM,TRIANGLE_NUM),
  !    the nodes that form the triangles.
  !

    integer(ip), intent(in), value :: dim_num
    integer(ip), intent(in), value :: node_num
    integer(ip), intent(in), value :: order_num
    integer(ip), intent(in), value :: triangle_num

    character ( len = *  ), intent(in) :: node_file_name
    real(dp), intent(in) :: node_xyz(dim_num,node_num)
    character ( len = *  ), intent(in) :: triangle_file_name
    integer(ip), intent(in) :: triangle_node(order_num,triangle_num)

    call r8mat_transpose_print ( dim_num, node_num, node_xyz, &
      '  Node coordinates' )

    call i4mat_transpose_print ( order_num, triangle_num, triangle_node, &
      '  Triangle nodes' )
  end subroutine tri_surface_print

  subroutine tri_surface_read ( node_file_name, triangle_file_name, dim_num, &
    node_num, order_num, triangle_num, node_xyz, triangle_node ) &
        bind(C, name="tri_surface_read")

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
  !    Input, integer(ip) DIM_NUM, the spatial dimension.
  !
  !    Input, integer(ip) NODE_NUM, the number of points.
  !
  !    Input, integer(ip) ORDER_NUM, the order of the triangles.
  !
  !    Input, integer(ip) TRIANGLE_NUM, the number of triangles.
  !
  !    Output, real(dp) NODE_XYZ(DIM_NUM,NODE_NUM), the node coordinates.
  !
  !    Output, integer(ip) TRIANGLE_NODE(ORDER_NUM,TRIANGLE_NUM),
  !    the nodes that form the triangles.
  !

    integer(ip), intent(in), value :: dim_num
    integer(ip), intent(in), value :: node_num
    integer(ip), intent(in), value :: order_num
    integer(ip), intent(in), value :: triangle_num

    character ( len = * ), intent(in) :: node_file_name
    real(dp), intent(out) :: node_xyz(dim_num,node_num)
    character ( len = * ), intent(in) :: triangle_file_name
    integer(ip), intent(out) :: triangle_node(order_num,triangle_num)

    call r8mat_data_read ( node_file_name, dim_num, node_num, node_xyz )

    call i4mat_data_read ( triangle_file_name, order_num, triangle_num, &
      triangle_node )
  end subroutine tri_surface_read

  subroutine tri_surface_size ( node_file_name, triangle_file_name, dim_num, &
    node_num, order_num, triangle_num ) &
        bind(C, name="tri_surface_size")

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
  !    Output, integer(ip) DIM_NUM, the spatial dimension.
  !
  !    Output, integer(ip) NODE_NUM, the number of nodes.
  !
  !    Output, integer(ip) ORDER_NUM, the order of the triangles.
  !
  !    Output, integer(ip) TRIANGLE_NUM, the number of triangles.
  !

    integer(ip), intent(out) :: dim_num
    character ( len = * ), intent(in) :: node_file_name
    integer(ip), intent(out) :: node_num
    integer(ip), intent(out) :: order_num
    character ( len = * ), intent(in) :: triangle_file_name
    integer(ip), intent(out) :: triangle_num

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
  end subroutine tri_surface_size

  subroutine tri_surface_size_print ( node_file_name, triangle_file_name, &
    dim_num, node_num, order_num, triangle_num ) &
        bind(C, name="tri_surface_size_print")

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
  !    Input, integer(ip) DIM_NUM, the number of spatial dimensions.
  !
  !    Input, integer(ip) NODE_NUM, the number of nodes.
  !
  !    Input, integer(ip) ORDER_NUM, the order of the triangles.
  !
  !    Input, integer(ip) TRIANGLE_NUM, the number of triangles.
  !

    integer(ip), intent(in), value :: dim_num
    character ( len = * ), intent(in) :: node_file_name
    integer(ip), intent(in), value :: node_num
    integer(ip), intent(in), value :: order_num
    character ( len = * ), intent(in) :: triangle_file_name
    integer(ip), intent(in), value :: triangle_num

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
  end subroutine tri_surface_size_print

  subroutine tri_surface_write ( node_file_name, triangle_file_name, dim_num, &
    node_num, order_num, triangle_num, node_xyz, triangle_node ) &
        bind(C, name="tri_surface_write")

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
  !    Input, integer(ip) DIM_NUM, the spatial dimension.
  !
  !    Input, integer(ip) NODE_NUM, the number of points.
  !
  !    Input, integer(ip) ORDER_NUM, the order of the triangles.
  !
  !    Input, integer(ip) TRIANGLE_NUM, the number of triangles.
  !
  !    Input, real(dp) NODE_XYZ(DIM_NUM,NODE_NUM), the node coordinates.
  !
  !    Input, integer(ip) TRIANGLE_NODE(ORDER_NUM,TRIANGLE_NUM),
  !    the nodes that form the triangles.
  !

    integer(ip), intent(in), value :: dim_num
    integer(ip), intent(in), value :: node_num
    integer(ip), intent(in), value :: order_num
    integer(ip), intent(in), value :: triangle_num

    character ( len = * ), intent(in) :: node_file_name
    real(dp), intent(in) :: node_xyz(dim_num,node_num)
    character ( len = * ), intent(in) :: triangle_file_name
    integer(ip), intent(in) :: triangle_node(order_num,triangle_num)

    call r8mat_write ( node_file_name, dim_num, node_num, node_xyz )

    call i4mat_write ( triangle_file_name, order_num, triangle_num, &
      triangle_node )
  end subroutine tri_surface_write

end module tri_surface_io_mod
