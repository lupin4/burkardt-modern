!> polygon_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module polygon_grid_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  implicit none
  private

  public :: polygon_grid_count, polygon_grid_display, polygon_grid_points

contains

  subroutine polygon_grid_count ( n, nv, ng )

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
  !    Input, integer(int32) N, the number of subintervals on a side.
  !
  !    Input, integer(int32) NV, the number of vertices.
  !    3 <= NV.
  !
  !    Output, integer(int32) NG, the number of grid points.
  !

    integer(int32) n
    integer(int32) ng
    integer(int32) nv

    ng = 1 + nv * ( n * ( n + 1 ) ) / 2
  end

  subroutine polygon_grid_display ( n, nv, v, ng, xg, prefix )

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
  !    Input, integer(int32) N, the number of subintervals.
  !
  !    Input, integer(int32) NV, the number of vertices in the polygon.
  !
  !    Input, real(real64) V(2,NV), the coordinates of the vertices.
  !
  !    Input, integer(int32) NG, the number of grid points.
  !
  !    Input, real(real64) XG(2,NG), the grid points.
  !
  !    Input, character ( len = * ) PREFIX, a string used to name the files.
  !

    integer(int32) ng
    integer(int32) nv

    integer(int32) command_unit
    character ( len = 255 ) command_filename
    integer(int32) grid_unit
    character ( len = 255 ) grid_filename
    integer(int32) j
    integer(int32) n
    character ( len = 255 ) plot_filename
    character ( len = * ) prefix
    real(real64) v(2,nv)
    real(real64) vc(2)
    integer(int32) vertex_unit
    character ( len = 255 ) vertex_filename
    real(real64) xg(2,ng)
  !
  !  Write the vertex file.
  !
    vc(1) = sum ( v(1,1:nv) ) / real ( nv, real64)
    vc(2) = sum ( v(2,1:nv) ) / real ( nv, real64)

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
  end

  subroutine polygon_grid_points ( n, nv, v, ng, xg )

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
  !    Input, integer(int32) N, the number of subintervals.
  !
  !    Input, integer(int32) NV, the number of vertices in the polygon.
  !
  !    Input, real(real64) V(2,NV), the coordinates of the vertices.
  !
  !    Input, integer(int32) NG, the number of grid points.
  !
  !    Output, real(real64) XG(2,NG), the coordinates of the grid points.
  !

    integer(int32) ng
    integer(int32) nv

    integer(int32) i
    integer(int32) j
    integer(int32) k
    integer(int32) l
    integer(int32) lp1
    integer(int32) n
    integer(int32) p
    real(real64) v(2,nv)
    real(real64) vc(2)
    real(real64) xg(2,ng)
  !
  !  Determine the centroid, and use it as the first grid point.
  !
    p = 1
    vc(1) = sum ( v(1,1:nv) ) / real ( nv, real64)
    vc(2) = sum ( v(2,1:nv) ) / real ( nv, real64)
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
          xg(1:2,p) = ( real ( i, real64) * v(1:2,l)   &
                      + real ( j, real64) * v(1:2,lp1) & 
                      + real ( k, real64) * vc(1:2) )  &
                      / real ( n, real64)
        end do
      end do
    end do
  end

end module polygon_grid_mod
