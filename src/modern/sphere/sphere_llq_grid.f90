!> sphere_llq_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

subroutine sphere_llq_grid_display ( ng, xg, line_num, line_data, prefix )

!*****************************************************************************80
!
!! SPHERE_LLQ_GRID_DISPLAY displays an LLQ grid on a sphere.
!
!  Discussion:
!
!    A SPHERE LLQ grid imposes a grid of quadrilaterals on a sphere,
!    using latitude and longitude lines.
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
!    Input, integer NG, the number of points.
!
!    Input, double precision XG(3,NG), the points.
!
!    Input, integer(8) LINE_NUM, the number of grid lines.
!
!    Input, integer(8) LINE_DATA(2,LINE_NUM), contains pairs of 
!    point indices for line segments that make up the grid.
!
!    Input, character ( len = * ) PREFIX, a prefix for the filenames.
!
  implicit none

  integer line_num
  integer ng

  character ( len = 255 ) command_filename
  integer command_unit
  integer j
  integer j1
  integer j2
  integer l
  integer line_data(2,line_num)
  character ( len = 255 ) line_filename
  integer line_unit
  character ( len = 255 ) node_filename
  integer node_unit
  character ( len = 255 ) plot_filename
  character ( len = * ) prefix
  double precision xg(3,ng)
!
!  Create graphics data files.
!
  call get_unit ( node_unit )
  node_filename = trim ( prefix ) // '_nodes.txt'
  open ( unit = node_unit, file = node_filename, status = 'replace' )
  do j = 1, ng
    write ( node_unit, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) xg(1:3,j)
  end do
  close ( unit = node_unit )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Created node file "' // trim ( node_filename ) // '".'

  call get_unit ( line_unit )
  line_filename = trim ( prefix ) // '_lines.txt'
  open ( unit = line_unit, file = line_filename, status = 'replace' )
  do l = 1, line_num
    if ( 1 < l ) then
      write ( line_unit, '(a)' ) ''
      write ( line_unit, '(a)' ) ''
    end if
    j1 = line_data(1,l)
    j2 = line_data(2,l)
    write ( line_unit, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) xg(1:3,j1)
    write ( line_unit, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) xg(1:3,j2)
  end do
  close ( unit = line_unit )
  write ( *, '(a)' ) '  Created line file "' // trim ( line_filename ) // '".'
!
!  Create graphics command file.
!
  call get_unit ( command_unit )
  command_filename = trim ( prefix ) // '_commands.txt'
  open ( unit = command_unit, file = command_filename, status = 'replace' )
  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  plot_filename = trim ( prefix ) // '.png'
  write ( command_unit, '(a)' ) 'set output "' // trim ( plot_filename ) // '"'
  write ( command_unit, '(a)' ) 'set xlabel "<--- X --->"'
  write ( command_unit, '(a)' ) 'set ylabel "<--- Y --->"'
  write ( command_unit, '(a)' ) 'set zlabel "<--- Z --->"'
  write ( command_unit, '(a)' ) 'set title "' // trim ( prefix ) // '"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set key off'
  write ( command_unit, '(a)' ) 'set style data points'
  write ( command_unit, '(a)' ) 'set timestamp'
  write ( command_unit, '(a)' ) 'set view equal xyz'
  write ( command_unit, '(a)' ) 'splot "' // &
    trim ( line_filename ) // &
    '" with lines lw 3, \'
  write ( command_unit, '(a)' ) '     "' // &
    trim ( node_filename ) // '" with points pt 7 lt 0'
  write ( command_unit, '(a)' ) 'quit'
  close ( unit = command_unit )

  write ( *, '(a)' ) &
    '  Created command file "' // trim ( command_filename ) // '".'
end

subroutine sphere_llq_grid_lines ( lat_num, long_num, line_num, line )

!*****************************************************************************80
!
!! SPHERE_LLQ_GRID_LINES: latitude/longitude quadrilateral grid lines.
!
!  Discussion:
!
!    A SPHERE LLQ grid imposes a grid of quadrilaterals on a sphere,
!    using latitude and longitude lines.
!
!    The point numbering system is the same used in SPHERE_LLQ_POINTS,
!    and that routine may be used to compute the coordinates of the points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LAT_NUM, LONG_NUM, the number of latitude and
!    longitude lines to draw.  The latitudes do not include the North and South
!    poles, which will be included automatically, so LAT_NUM = 5, for instance,
!    will result in points along 7 lines of latitude.
!
!    Input, integer LINE_NUM, the number of grid lines.
!
!    Output, integer LINE(2,LINE_NUM), contains pairs of point 
!    indices for line segments that make up the grid.
!
  implicit none

  integer line_num

  integer i
  integer j
  integer lat_num
  integer l
  integer line(2,line_num)
  integer long_num
  integer new
  integer newcol
  integer old

  l = 0
!
!  "Vertical" lines.
!
  do j = 0, long_num - 1

    old = 1
    new = j + 2

    l = l + 1
    line(1:2,l) = (/ old, new /)

    do i = 1, lat_num - 1

      old = new
      new = old + long_num

      l = l + 1
      line(1:2,l) = (/ old, new /)

    end do

    old = new

    l = l + 1
    line(1:2,l) = (/ old, 1 + lat_num * long_num + 1 /)

  end do
!
!  "Horizontal" lines.
!
  do i = 1, lat_num

    new = 1 + ( i - 1 ) * long_num + 1

    do j = 0, long_num - 2
      old = new
      new = old + 1
      l = l + 1
      line(1:2,l) = (/ old, new /)
    end do

    old = new
    new = 1 + ( i - 1 ) * long_num + 1
    l = l + 1
    line(1:2,l) = (/ old, new /)

  end do
end

subroutine sphere_llq_grid_line_count ( lat_num, long_num, line_num )

!*****************************************************************************80
!
!! SPHERE_LLQ_GRID_LINE_COUNT counts latitude/longitude quad grid lines.
!
!  Discussion:
!
!    The number returned is the number of pairs of points to be connected.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LAT_NUM, LONG_NUM, the number of latitude and
!    longitude lines to draw.  The latitudes do not include the North and South
!    poles, which will be included automatically, so LAT_NUM = 5, for instance,
!    will result in points along 7 lines of latitude.
!
!    Output, integer LINE_NUM, the number of grid lines.
!
  implicit none

  integer lat_num
  integer line_num
  integer long_num

  line_num = long_num * ( lat_num + 1 ) &
           + lat_num * long_num
end

subroutine sphere_llq_grid_points ( r, pc, lat_num, long_num, point_num, p )

!*****************************************************************************80
!
!! SPHERE_LLQ_GRID_POINTS produces points for a latitude/longitude quad grid.
!
!  Discussion:
!
!    A SPHERE LLQ grid imposes a grid of quadrilaterals on a sphere,
!    using latitude and longitude lines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision R, the radius of the sphere.
!
!    Input, double precision PC(3), the center of the sphere.
!
!    Input, integer LAT_NUM, LONG_NUM, the number of latitude 
!    and longitude lines to draw.  The latitudes do not include the North and 
!    South poles, which will be included automatically, so LAT_NUM = 5, for 
!    instance, will result in points along 7 lines of latitude.
!
!    Input, integer POINT_NUM, the number of points.
!
!    Output, double precision P(3,POINT_NUM), the grid points.
!
  implicit none

  integer lat_num
  integer long_num
  integer point_num

  integer lat
  integer long
  integer n
  double precision p(3,point_num)
  double precision pc(3)
  double precision phi
  double precision r
  double precision , parameter :: r8_pi = 3.141592653589793D+00
  double precision theta

  n = 0
!
!  The north pole.
!
  theta = 0.0D+00
  phi = 0.0D+00
  n = n + 1
  p(1,n) = pc(1) + r * sin ( phi ) * cos ( theta )
  p(2,n) = pc(2) + r * sin ( phi ) * sin ( theta )
  p(3,n) = pc(3) + r * cos ( phi )
!
!  Do each intermediate ring of latitude.
!
  do lat = 1, lat_num

    phi = real ( lat) * r8_pi &
        / real ( lat_num + 1)
!
!  Along that ring of latitude, compute points at various longitudes.
!
    do long = 0, long_num - 1

      theta = real ( long) * 2.0D+00 * r8_pi &
            / real ( long_num)

      n = n + 1
      p(1,n) = pc(1) + r * sin ( phi ) * cos ( theta )
      p(2,n) = pc(2) + r * sin ( phi ) * sin ( theta )
      p(3,n) = pc(3) + r * cos ( phi )

    end do
  end do
!
!  The south pole.
!
  theta = 0.0D+00
  phi = r8_pi
  n = n + 1
  p(1,n) = pc(1) + r * sin ( phi ) * cos ( theta )
  p(2,n) = pc(2) + r * sin ( phi ) * sin ( theta )
  p(3,n) = pc(3) + r * cos ( phi )
end

subroutine sphere_llq_grid_point_count ( lat_num, long_num, point_num )

!*****************************************************************************80
!
!! SPHERE_LLQ_GRID_POINT_COUNT counts points for a latitude/longitude grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LAT_NUM, LONG_NUM, the number of latitude 
!    and longitude lines to draw.  The latitudes do not include the North and 
!    South poles, which will be included automatically, so LAT_NUM = 5, for 
!    instance, will result in points along 7 lines of latitude.
!
!    Output, integer POINT_NUM, the number of grid points.
!
  implicit none

  integer lat_num
  integer long_num
  integer point_num

  point_num = 2 + lat_num * long_num
end
