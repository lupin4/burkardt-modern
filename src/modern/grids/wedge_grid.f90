!> wedge_grid -- Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module wedge_grid_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: wedge_grid, wedge_grid_size, wedge_grid_plot, wedge_vertices

contains

  pure subroutine wedge_grid ( n, ng, g ) &
        bind(C, name="wedge_grid")

  !*****************************************************************************80
  !
  !! WEDGE_GRID computes grid points in the unit wedge in 3D.
  !
  !  Discussion:
  !
  !    The interior of the unit wedge in 3D is defined by the constraints:
  !      0 <= X
  !      0 <= Y
  !           X + Y <= 1
  !     -1 <= Z <= +1
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 August 2014
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of subintervals.
  !    0 <= N.
  !
  !    Input, integer(ip) NG, the number of grid points.
  !    This can be computed by WEDGE_GRID_SIZE, or else determined by
  !    NG =(N+1)*((N+1)*(N+2))/2.
  !
  !    Output, real(dp) G(3,NG), the coordinates
  !    of the grid points.
  !

    integer(ip), intent(in), value :: n
    integer(ip), intent(in), value :: ng
    integer(ip) :: i
    real(dp) :: ir
    integer(ip) :: j
    real(dp) :: jr
    integer(ip) :: k
    real(dp) :: kr
    real(dp) :: nr
    integer(ip) :: p
    real(dp), intent(out) :: g(3,ng)

    if ( n == 0 ) then
      g(1,1) = 0.5_dp
      g(2,1) = 0.5_dp
      g(3,1) = 0.0_dp
    end if

    p = 0
    nr = real ( n, dp)

    do k = 0, n
      kr = real ( 2 * k - n, dp) / nr
      do j = 0, n
        jr = real ( j, dp) / nr
        do i = 0, n - j
          ir = real ( i, dp) / nr
          p = p + 1
          g(1,p) = ir
          g(2,p) = jr
          g(3,p) = kr
        end do
      end do
    end do
  end subroutine wedge_grid

  pure subroutine wedge_grid_size ( n, ng ) &
        bind(C, name="wedge_grid_size")

  !*****************************************************************************80
  !
  !! WEDGE_GRID_SIZE counts the points in a grid of the unit wedge in 3D.
  !
  !  Discussion:
  !
  !    The interior of the unit wedge in 3D is defined by the constraints:
  !      0 <= X
  !      0 <= Y
  !           X + Y <= 1
  !     -1 <= Z <= +1
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 August 2014
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of subintervals.
  !    0 <= N.
  !
  !    Output, integer(ip) NG, the number of grid points.
  !

    integer(ip), intent(in), value :: n
    integer(ip), intent(out) :: ng

    ng = ( n + 1 ) * ( ( n + 1 ) * ( n + 2 ) ) / 2
  end subroutine wedge_grid_size

  subroutine wedge_grid_plot ( n, ng, g, header ) &
        bind(C, name="wedge_grid_plot")

  !*****************************************************************************80
  !
  !! WEDGE_GRID_PLOT sets up a GNUPLOT plot of a unit wedge grid.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 August 2014
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of subintervals.
  !
  !    Input, integer(ip) NG, the number of nodes.
  !
  !    Input, real(dp) G(3,NG), the grid point coordinates.
  !
  !    Input, character ( len = * ) HEADER, the header for the files.
  !

    integer(ip), intent(in), value :: ng
    character ( len = 255 ) :: command_filename
    integer(ip) :: command_unit
    real(dp), intent(in) :: g(3,ng)
    character ( len = * ), intent(in) :: header
    integer(ip) :: j
    integer(ip), intent(in), value :: n
    character ( len = 255 ) :: node_filename
    integer(ip) :: node_unit
    character ( len = 255 ) :: plot_filename
    real(dp) :: v1(3)
    real(dp) :: v2(3)
    real(dp) :: v3(3)
    real(dp) :: v4(3)
    real(dp) :: v5(3)
    real(dp) :: v6(3)
    character ( len = 255 ) :: vertex_filename
    integer(ip) :: vertex_unit
  !
  !  Create the vertex file.
  !
    call wedge_vertices ( v1, v2, v3, v4, v5, v6 )

    call get_unit ( vertex_unit )
    vertex_filename = trim ( header ) // '_vertices.txt'
    open ( unit = vertex_unit, file = vertex_filename, &
      status = 'replace' )

    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v1(1), v1(2), v1(3)
    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v2(1), v2(2), v2(3)
    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v3(1), v3(2), v3(3)
    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v1(1), v1(2), v1(3)
    write ( vertex_unit, '(a)' ) ''

    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v4(1), v4(2), v4(3)
    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v5(1), v5(2), v5(3)
    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v6(1), v6(2), v6(3)
    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v4(1), v4(2), v4(3)
    write ( vertex_unit, '(a)' ) ''

    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v1(1), v1(2), v1(3)
    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v4(1), v4(2), v4(3)
    write ( vertex_unit, '(a)' ) ''

    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v2(1), v2(2), v2(3)
    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v5(1), v5(2), v5(3)
    write ( vertex_unit, '(a)' ) ''

    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v3(1), v3(2), v3(3)
    write ( vertex_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) &
      v6(1), v6(2), v6(3)
    write ( vertex_unit, '(a)' ) ''

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
      write ( node_unit, '(g14.6,2x,g14.6,2x,g14.6)' ) g(1:3,j)
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
    write ( command_unit, '(a)' ) '#set view equal xyz'
    write ( command_unit, '(a)' ) 'set view 80, 85'
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
  end subroutine wedge_grid_plot

  pure subroutine wedge_vertices ( v1, v2, v3, v4, v5, v6 ) &
        bind(C, name="wedge_vertices")

  !*****************************************************************************80
  !
  !! WEDGE_VERTICES returns the vertices of the unit wege.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 August 2014
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Output, real(dp) V1(3), V2(3), V3(3), V4(3), V5(3), V6(3),
  !    the vertices.
  !

    real(dp), intent(out) :: v1(3)
    real(dp), intent(out) :: v2(3)
    real(dp), intent(out) :: v3(3)
    real(dp), intent(out) :: v4(3)
    real(dp), intent(out) :: v5(3)
    real(dp), intent(out) :: v6(3)

    v1(1:3) = (/  0.0_dp,  0.0_dp, -1.0_dp /)
    v2(1:3) = (/  1.0_dp,  0.0_dp, -1.0_dp /)
    v3(1:3) = (/  0.0_dp,  1.0_dp, -1.0_dp /)
    v4(1:3) = (/  0.0_dp,  0.0_dp, +1.0_dp /)
    v5(1:3) = (/  1.0_dp,  0.0_dp, +1.0_dp /)
    v6(1:3) = (/  0.0_dp,  1.0_dp, +1.0_dp /)
  end subroutine wedge_vertices

end module wedge_grid_mod
