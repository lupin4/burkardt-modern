!> pyramid_grid -- Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module pyramid_grid_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: pyramid_grid_size, pyramid_unit_grid, pyramid_unit_grid_plot, pyramid_unit_vertices, r8_print

contains

  pure function pyramid_grid_size ( n ) &
        bind(C, name="pyramid_grid_size")

  !*****************************************************************************80
  !
  !! PYRAMID_GRID_SIZE sizes a pyramid grid.
  !
  !  Discussion:
  !
  !    0:  x
  !
  !    1:  x  x
  !        x  x
  !
  !    2:  x  x  x
  !        x  x  x
  !        x  x  x
  !
  !    3:  x  x  x  x
  !        x  x  x  x
  !        x  x  x  x
  !        x  x  x  x
  !
  !    N  Size
  !
  !    0     1
  !    1     5 = 1 + 4
  !    2    14 = 1 + 4 + 9
  !    3    30 = 1 + 4 + 9 + 16
  !    4    55 = 1 + 4 + 9 + 16 + 25
  !    5    91 = 1 + 4 + 9 + 16 + 25 + 36
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 August 2014
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of subintervals.
  !
  !    Output, integer(ip) PYRAMID_GRID_SIZE, the number of
  !    nodes in the grid of size N.
  !

    integer(ip), intent(in), value :: n
    integer(ip) :: np1
    integer(ip) :: pyramid_grid_size
    integer(ip) :: value

    np1 = n + 1

    value = ( np1 * ( np1 + 1 ) * ( 2 * np1 + 1 ) ) / 6

    pyramid_grid_size = value
  end function pyramid_grid_size

  pure subroutine pyramid_unit_grid ( n, ng, pg ) &
        bind(C, name="pyramid_unit_grid")

  !*****************************************************************************80
  !
  !! PYRAMID_UNIT_GRID computes grid points in the unit pyramid.
  !
  !  Discussion:
  !
  !    The unit pyramid has base (-1,-1,0), (+1,1,0), (+1,+1,0), (-1,+1,0)
  !    and vertex (0,0,1).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 August 2014
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of subintervals.
  !
  !    Input, integer(ip) NG, the number of nodes to generate,
  !    as determined by pyramid_grid_size().
  !
  !    Output, real(dp) PG(3,NG), the grid point coordinates.
  !

    integer(ip), intent(in), value :: ng
    integer(ip) :: g
    integer(ip) :: hi
    integer(ip) :: i
    integer(ip) :: j
    integer(ip) :: k
    integer(ip) :: lo
    integer(ip), intent(in), value :: n
    real(dp), intent(out) :: pg(3,ng)

    g = 0

    do k = n, 0, -1
      hi = n - k
      lo = - hi
      do j = lo, hi, 2
        do i = lo, hi, 2
          g = g + 1
          pg(1,g) = real ( i, dp) / real ( n, dp)
          pg(2,g) = real ( j, dp) / real ( n, dp)
          pg(3,g) = real ( k, dp) / real ( n, dp)
        end do
      end do
    end do
  end subroutine pyramid_unit_grid

  subroutine pyramid_unit_grid_plot ( n, ng, pg, header ) &
        bind(C, name="pyramid_unit_grid_plot")

  !*****************************************************************************80
  !
  !! PYRAMID_UNIT_GRID_PLOT sets up a GNUPLOT plot of a unit pyramid grid.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 August 2014
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of subintervals.
  !
  !    Input, integer(ip) NG, the number of nodes to generate,
  !    as determined by pyramid_grid_size().
  !
  !    Input, real(dp) PG(3,NG), the grid point coordinates.
  !
  !    Input, character ( len = * ) HEADER, the header for the files.
  !

    integer(ip), intent(in), value :: ng
    character ( len = 255 ) :: command_filename
    integer(ip) :: command_unit
    character ( len = * ), intent(in) :: header
    integer(ip) :: i
    integer(ip) :: j
    integer(ip) :: l
    integer(ip), intent(in), value :: n
    character ( len = 255 ) :: node_filename
    integer(ip) :: node_unit
    real(dp), intent(in) :: pg(3,ng)
    character ( len = 255 ) :: plot_filename
    real(dp) :: v1(3)
    real(dp) :: v2(3)
    real(dp) :: v3(3)
    real(dp) :: v4(3)
    real(dp) :: v5(3)
    character ( len = 255 ) :: vertex_filename
    integer(ip) :: vertex_unit
  !
  !  Create the vertex file.
  !
    call pyramid_unit_vertices ( v1, v2, v3, v4, v5 )

    call get_unit ( vertex_unit )
    vertex_filename = trim ( header ) // '_vertices.txt'
    open ( unit = vertex_unit, file = vertex_filename, &
      status = 'replace' )

    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v2(1), v2(2), v2(3)
    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v3(1), v3(2), v3(3)
    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v4(1), v4(2), v4(3)
    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v5(1), v5(2), v5(3)
    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v2(1), v2(2), v2(3)
    write ( vertex_unit, '(a)' ) ''

    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v1(1), v1(2), v1(3)
    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v2(1), v2(2), v2(3)
    write ( vertex_unit, '(a)' ) ''

    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v1(1), v1(2), v1(3)
    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v3(1), v3(2), v3(3)
    write ( vertex_unit, '(a)' ) ''

    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v1(1), v1(2), v1(3)
    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v4(1), v4(2), v4(3)
    write ( vertex_unit, '(a)' ) ''

    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v1(1), v1(2), v1(3)
    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v5(1), v5(2), v5(3)

    close ( unit = vertex_unit )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Created vertex file "' // &
      trim ( vertex_filename ) // '".'
  !
  !  Create the node file.
  !
    call get_unit ( node_unit )
    node_filename = trim ( header ) // '_nodes.txt'
    open ( unit = node_unit, file = node_filename, &
      status = 'replace' )
    do j = 1, ng
      write ( node_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) pg(1:3,j)
    end do
    close ( unit = node_unit )
    write ( *, '(a)' ) '  Created node file "' // &
      trim ( node_filename ) // '".'
  !
  !  Create the command file.
  !
    call get_unit ( command_unit )
    command_filename = trim ( header ) // '_commands.txt'
    open ( unit = command_unit, file = command_filename, &
      status = 'replace' )
    write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
    write ( command_unit, '(a)' ) '#'
    write ( command_unit, '(a)' ) '# Usage:'
    write ( command_unit, '(a)' ) '#  gnuplot < ' // &
      trim ( command_filename )
    write ( command_unit, '(a)' ) '#'
    write ( command_unit, '(a)' ) 'set term png'
    plot_filename = trim ( header ) // '.png'
    write ( command_unit, '(a)' ) 'set output "' // &
      trim ( plot_filename ) // '"'
    write ( command_unit, '(a)' ) 'set xlabel "<--- X --->"'
    write ( command_unit, '(a)' ) 'set ylabel "<--- Y --->"'
    write ( command_unit, '(a)' ) 'set zlabel "<--- Z --->"'
    write ( command_unit, '(a)' ) &
      'set title "' // trim ( header ) // '"'
    write ( command_unit, '(a)' ) 'set grid'
    write ( command_unit, '(a)' ) 'set key off'
    write ( command_unit, '(a)' ) 'set view equal xyz'
    write ( command_unit, '(a)' ) 'set view 80, 40'
    write ( command_unit, '(a)' ) 'set style data lines'
    write ( command_unit, '(a)' ) 'set timestamp'
    write ( command_unit, '(a)' ) 'splot "' // &
      trim ( vertex_filename ) // &
      '" with lines lw 3, \'
    write ( command_unit, '(a)' ) '     "' // &
      trim ( node_filename ) // '" with points pt 7 lt 0'
    close ( unit = command_unit )

    write ( *, '(a)' ) &
      '  Created command file "' // trim ( command_filename ) // '".'
  end subroutine pyramid_unit_grid_plot

  pure subroutine pyramid_unit_vertices ( v1, v2, v3, v4, v5 ) &
        bind(C, name="pyramid_unit_vertices")

  !*****************************************************************************80
  !
  !! PYRAMID_UNIT_VERTICES returns the vertices of the unit pyramid.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 August 2014
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Output, real(dp) V1(3), V2(3), V3(3), V4(3), V5(3), the vertices.
  !

    real(dp), intent(out) :: v1(3)
    real(dp), intent(out) :: v2(3)
    real(dp), intent(out) :: v3(3)
    real(dp), intent(out) :: v4(3)
    real(dp), intent(out) :: v5(3)

    v1(1:3) = (/  0.0_dp,  0.0_dp, +1.0_dp /)
    v2(1:3) = (/ -1.0_dp, -1.0_dp,  0.0_dp /)
    v3(1:3) = (/ +1.0_dp, -1.0_dp,  0.0_dp /)
    v4(1:3) = (/ +1.0_dp, +1.0_dp,  0.0_dp /)
    v5(1:3) = (/ -1.0_dp, +1.0_dp,  0.0_dp /)
  end subroutine pyramid_unit_vertices

  subroutine r8_print ( r, title ) &
        bind(C, name="r8_print")

  !*****************************************************************************80
  !
  !! R8_PRINT prints an R8.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 August 2014
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real(dp) R, the value.
  !
  !    Input, character ( len = * ) TITLE, a title.
  !

    real(dp), intent(in), value :: r
    character ( len = * ), intent(in) :: title

    write ( *, '(a,2x,g14.6)' ) trim ( title ), r
  end subroutine r8_print

end module pyramid_grid_mod
