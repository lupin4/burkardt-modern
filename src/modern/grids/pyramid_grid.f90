!> pyramid_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

function pyramid_grid_size ( n )

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
!    Input, integer N, the number of subintervals.
!
!    Output, integer PYRAMID_GRID_SIZE, the number of
!    nodes in the grid of size N.
!
  implicit none

  integer n
  integer np1
  integer pyramid_grid_size
  integer value

  np1 = n + 1

  value = ( np1 * ( np1 + 1 ) * ( 2 * np1 + 1 ) ) / 6

  pyramid_grid_size = value
end

subroutine pyramid_unit_grid ( n, ng, pg )

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
!    Input, integer N, the number of subintervals.
!
!    Input, integer NG, the number of nodes to generate,
!    as determined by pyramid_grid_size().
!
!    Output, double precision PG(3,NG), the grid point coordinates.
!
  implicit none

  integer ng

  integer g
  integer hi
  integer i
  integer j
  integer k
  integer lo
  integer n
  double precision pg(3,ng)

  g = 0

  do k = n, 0, -1
    hi = n - k
    lo = - hi
    do j = lo, hi, 2
      do i = lo, hi, 2
        g = g + 1
        pg(1,g) = real ( i) / real ( n)
        pg(2,g) = real ( j) / real ( n)
        pg(3,g) = real ( k) / real ( n)
      end do
    end do
  end do
end

subroutine pyramid_unit_grid_plot ( n, ng, pg, header )

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
!    Input, integer N, the number of subintervals.
!
!    Input, integer NG, the number of nodes to generate,
!    as determined by pyramid_grid_size().
!
!    Input, double precision PG(3,NG), the grid point coordinates.
!
!    Input, character ( len = * ) HEADER, the header for the files.
!
  implicit none

  integer ng

  character ( len = 255 ) command_filename
  integer command_unit
  character ( len = * ) header
  integer i
  integer j
  integer l
  integer n
  character ( len = 255 ) node_filename
  integer node_unit
  double precision pg(3,ng)
  character ( len = 255 ) plot_filename
  double precision v1(3)
  double precision v2(3)
  double precision v3(3)
  double precision v4(3)
  double precision v5(3)
  character ( len = 255 ) vertex_filename
  integer vertex_unit
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
end

subroutine pyramid_unit_vertices ( v1, v2, v3, v4, v5 )

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
!    Output, double precision V1(3), V2(3), V3(3), V4(3), V5(3), the vertices.
!
  implicit none

  double precision v1(3)
  double precision v2(3)
  double precision v3(3)
  double precision v4(3)
  double precision v5(3)

  v1(1:3) = (/  0.0D+00,  0.0D+00, +1.0D+00 /)
  v2(1:3) = (/ -1.0D+00, -1.0D+00,  0.0D+00 /)
  v3(1:3) = (/ +1.0D+00, -1.0D+00,  0.0D+00 /)
  v4(1:3) = (/ +1.0D+00, +1.0D+00,  0.0D+00 /)
  v5(1:3) = (/ -1.0D+00, +1.0D+00,  0.0D+00 /)
end

subroutine r8_print ( r, title )

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
!    Input, double precision R, the value.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  double precision r
  character ( len = * ) title

  write ( *, '(a,2x,g14.6)' ) trim ( title ), r
end
