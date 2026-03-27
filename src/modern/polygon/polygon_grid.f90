!> polygon_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module polygon_grid_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: polygon_grid_count, polygon_grid_display, polygon_grid_points

contains

  pure subroutine polygon_grid_count ( n, nv, ng ) &
        bind(C, name="polygon_grid_count")

  !*****************************************************************************80
  !
  !! POLYGON_GRID_COUNT counts the grid points inside a polygon.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    09 May 2015
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of subintervals on a side.
  !
  !    Input, integer(ip) NV, the number of vertices.
  !    3 <= NV.
  !
  !    Output, integer(ip) NG, the number of grid points.
  !

    integer(ip), intent(in), value :: n
    integer(ip), intent(in), value :: nv
    integer(ip), intent(out) :: ng

    ng = 1 + nv * ( n * ( n + 1 ) ) / 2
  end subroutine polygon_grid_count

  subroutine polygon_grid_display ( n, nv, v, ng, xg, prefix ) &
        bind(C, name="polygon_grid_display")

  !*****************************************************************************80
  !
  !! POLYGON_GRID_DISPLAY displays grid points inside a polygon.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    10 May 2015
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of subintervals.
  !
  !    Input, integer(ip) NV, the number of vertices in the polygon.
  !
  !    Input, real(dp) V(2,NV), the coordinates of the vertices.
  !
  !    Input, integer(ip) NG, the number of grid points.
  !
  !    Input, real(dp) XG(2,NG), the grid points.
  !
  !    Input, character ( len = * ) PREFIX, a string used to name the files.
  !

    integer(ip), intent(in), value :: ng
    integer(ip), intent(in), value :: nv

    integer(ip) :: command_unit
    character ( len = 255 ) :: command_filename
    integer(ip) :: grid_unit
    character ( len = 255 ) :: grid_filename
    integer(ip) :: j
    integer(ip), intent(in), value :: n
    character ( len = 255 ) :: plot_filename
    character ( len = * ), intent(in) :: prefix
    real(dp), intent(in) :: v(2,nv)
    real(dp) :: vc(2)
    integer(ip) :: vertex_unit
    character ( len = 255 ) :: vertex_filename
    real(dp), intent(in) :: xg(2,ng)
  !
  !  Write the vertex file.
  !
    vc(1) = sum ( v(1,1:nv) ) / real ( nv, dp)
    vc(2) = sum ( v(2,1:nv) ) / real ( nv, dp)

    call get_unit ( vertex_unit )
    vertex_filename = trim ( prefix ) // '_vertex.txt'
    open ( unit = vertex_unit, file = vertex_filename, status = 'replace' )
    do j = 1, nv
      write ( vertex_unit, '(2x,g14.6,2x,g14.6)' ) v(1,j), v(2,j)
    end do
    write ( vertex_unit, '(2x,g14.6,2x,g14.6)' ) v(1,1), v(2,1)
    do j = 1, nv
      write ( vertex_unit, '(a)' ) ''
      write ( vertex_unit, '(2x,g14.6,2x,g14.6)' ) v(1,j), v(2,j)
      write ( vertex_unit, '(2x,g14.6,2x,g14.6)' ) vc(1), vc(2)
    end do
    close ( unit = vertex_unit )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) &
      '  Created vertex file "' // trim ( vertex_filename ) // '".'
  !
  !  Write the gridpoint file.
  !
    call get_unit ( grid_unit )
    grid_filename = trim ( prefix ) // '_grid.txt'
    open ( unit = grid_unit, file = grid_filename, status = 'replace' )
    do j = 1, ng
      write ( grid_unit, '(2x,g14.6,2x,g14.6)' ) xg(1,j), xg(2,j)
    end do
    close ( unit = grid_unit )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Created grid file "' // trim ( grid_filename ) // '".'
  !
  !  Write the command file.
  !
    plot_filename = trim ( prefix ) // '.png'
    call get_unit ( command_unit )
    command_filename = trim ( prefix ) // '_commands.txt'
    open ( unit = command_unit, file = command_filename, status = 'replace' )
    write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
    write ( command_unit, '(a)' ) '#'
    write ( command_unit, '(a)' ) '# Usage:'
    write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
    write ( command_unit, '(a)' ) '#'
    write ( command_unit, '(a)' ) 'set term png'
    write ( command_unit, '(a)' ) 'set output "' // trim ( plot_filename ) // '"'
    write ( command_unit, '(a)' ) 'set xlabel "<--- X --->"'
    write ( command_unit, '(a)' ) 'set ylabel "<--- Y --->"'
    write ( command_unit, '(a)' ) 'set title "' // trim ( prefix ) // '"'
    write ( command_unit, '(a)' ) 'set grid'
    write ( command_unit, '(a)' ) 'set key off'
    write ( command_unit, '(a)' ) 'set size ratio -1'
    write ( command_unit, '(a)' ) 'set style data lines'

    write ( command_unit, '(a)' ) 'plot "' // trim ( grid_filename ) // &
      '" using 1:2 with points lt 3 pt 3,\'
    write ( command_unit, '(a)' ) '    "' // trim ( vertex_filename ) // &
      '" using 1:2 lw 3 linecolor rgb "black"'
    write ( command_unit, '(a)' ) 'quit'
    close ( unit = command_unit )

    write ( *, '(a)' ) &
      '  Created command file "' // trim ( command_filename ) // '".'
  end subroutine polygon_grid_display

  pure subroutine polygon_grid_points ( n, nv, v, ng, xg ) &
        bind(C, name="polygon_grid_points")

  !*****************************************************************************80
  !
  !! POLYGON_GRID_POINTS computes points on a polygonal grid.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    09 May 2015
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of subintervals.
  !
  !    Input, integer(ip) NV, the number of vertices in the polygon.
  !
  !    Input, real(dp) V(2,NV), the coordinates of the vertices.
  !
  !    Input, integer(ip) NG, the number of grid points.
  !
  !    Output, real(dp) XG(2,NG), the coordinates of the grid points.
  !

    integer(ip), intent(in), value :: ng
    integer(ip), intent(in), value :: nv

    integer(ip) :: i
    integer(ip) :: j
    integer(ip) :: k
    integer(ip) :: l
    integer(ip) :: lp1
    integer(ip), intent(in), value :: n
    integer(ip) :: p
    real(dp), intent(in) :: v(2,nv)
    real(dp) :: vc(2)
    real(dp), intent(out) :: xg(2,ng)
  !
  !  Determine the centroid, and use it as the first grid point.
  !
    p = 1
    vc(1) = sum ( v(1,1:nv) ) / real ( nv, dp)
    vc(2) = sum ( v(2,1:nv) ) / real ( nv, dp)
    xg(1:2,p) = vc(1:2)
  !
  !  Consider each triangle formed by two consecutive vertices and the centroid,
  !  but skip the first line of points.
  !
    do l = 1, nv
      lp1 = mod ( l, nv ) + 1
      do i = 1, n
        do j = 0, n - i
          k = n - i - j
          p = p + 1
          xg(1:2,p) = ( real ( i, dp) * v(1:2,l)   &
                      + real ( j, dp) * v(1:2,lp1) &
                      + real ( k, dp) * vc(1:2) )  &
                      / real ( n, dp)
        end do
      end do
    end do
  end subroutine polygon_grid_points

end module polygon_grid_mod
