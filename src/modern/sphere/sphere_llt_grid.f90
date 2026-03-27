!> sphere_llt_grid -- Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module sphere_llt_grid_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: sphere_llt_grid_display, sphere_llt_grid_line_count, sphere_llt_grid_lines, sphere_llt_grid_point_count, sphere_llt_grid_points

contains

  subroutine sphere_llt_grid_display ( ng, xg, line_num, line_data, prefix ) &
        bind(C, name="sphere_llt_grid_display")

  !*****************************************************************************80
  !
  !! SPHERE_LLT_GRID_DISPLAY displays an LLT grid on a sphere.
  !
  !  Discussion:
  !
  !    A SPHERE LLT grid imposes a grid of triangles on a sphere,
  !    using latitude and longitude lines.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 May 2015
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) NG, the number of points.
  !
  !    Input, real(dp) XG(3,NG), the points.
  !
  !    Input, integer(ip) LINE_NUM, the number of grid lines.
  !
  !    Input, integer(ip) LINE_DATA(2,LINE_NUM), contains pairs of
  !    point indices for line segments that make up the grid.
  !
  !    Input, character ( len = * ) PREFIX, a prefix for the filenames.
  !

    integer(ip), intent(in), value :: line_num
    integer(ip), intent(in), value :: ng
    real(dp), intent(in) :: xg(3,ng)
    integer(ip), intent(in) :: line_data(2,line_num)
    character ( len = * ), intent(in) :: prefix

    character ( len = 255 ) :: command_filename
    integer(ip) :: command_unit
    integer(ip) :: j
    integer(ip) :: j1
    integer(ip) :: j2
    integer(ip) :: l
    character ( len = 255 ) :: line_filename
    integer(ip) :: line_unit
    character ( len = 255 ) :: node_filename
    integer(ip) :: node_unit
    character ( len = 255 ) :: plot_filename
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
  end subroutine sphere_llt_grid_display

  pure subroutine sphere_llt_grid_line_count ( lat_num, long_num, line_num ) &
        bind(C, name="sphere_llt_grid_line_count")

  !*****************************************************************************80
  !
  !! SPHERE_LLT_GRID_LINE_COUNT counts latitude/longitude triangle grid lines.
  !
  !  Discussion:
  !
  !    An LLT grid is a grid of triangles bounded by latitude and longitude
  !    lines over the surface of a sphere in 3D.
  !
  !    The number returned is the number of pairs of points to be connected.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 May 2015
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) LAT_NUM, LONG_NUM, the number of latitude and
  !    longitude lines to draw.  The latitudes do not include the North and South
  !    poles, which will be included automatically, so LAT_NUM = 5, for instance,
  !    will result in points along 7 lines of latitude.
  !
  !    Output, integer(ip) LINE_NUM, the number of grid lines.
  !

    integer(ip), intent(in), value :: lat_num
    integer(ip), intent(in), value :: long_num
    integer(ip), intent(out) :: line_num

    line_num = long_num * ( lat_num + 1 ) &
             + long_num *   lat_num &
             + long_num * ( lat_num - 1 )
  end subroutine sphere_llt_grid_line_count

  subroutine sphere_llt_grid_lines ( lat_num, long_num, line_num, line ) &
        bind(C, name="sphere_llt_grid_lines")

  !*****************************************************************************80
  !
  !! SPHERE_LLT_GRID_LINES: latitude/longitude triangle grid lines.
  !
  !  Discussion:
  !
  !    A SPHERE LLT grid imposes a grid of triangles on a sphere,
  !    using latitude and longitude lines.
  !
  !    The point numbering system is the same used in SPHERE_LLT_POINTS,
  !    and that routine may be used to compute the coordinates of the points.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 May 2015
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) LAT_NUM, LONG_NUM, the number of latitude and
  !    longitude lines to draw.  The latitudes do not include the North and South
  !    poles, which will be included automatically, so LAT_NUM = 5, for instance,
  !    will result in points along 7 lines of latitude.
  !
  !    Input, integer(ip) LINE_NUM, the number of grid lines.
  !
  !    Output, integer(ip) LINE(2,LINE_NUM), contains pairs of point
  !    indices for line segments that make up the grid.
  !

    integer(ip), intent(in), value :: lat_num
    integer(ip), intent(in), value :: long_num
    integer(ip), intent(in), value :: line_num
    integer(ip), intent(out) :: line(2,line_num)

    integer(ip) :: i
    integer(ip) :: j
    integer(ip) :: l
    integer(ip) :: new
    integer(ip) :: newcol
    integer(ip) :: old

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
  !
  !  "Diagonal" lines.
  !
    do j = 0, long_num - 1

      old = 1
      new = j + 2
      newcol = j

      do i = 1, lat_num - 1

        old = new
        new = old + long_num + 1
        newcol = newcol + 1
        if ( long_num - 1 < newcol ) then
          newcol = 0
          new = new - long_num
        end if

        l = l + 1
        line(1:2,l) = (/ old, new /)

      end do

    end do
  end subroutine sphere_llt_grid_lines

  pure subroutine sphere_llt_grid_point_count ( lat_num, long_num, point_num ) &
        bind(C, name="sphere_llt_grid_point_count")

  !*****************************************************************************80
  !
  !! SPHERE_LLT_GRID_POINT_COUNT counts points for a latitude/longitude grid.
  !
  !  Discussion:
  !
  !    An LLT grid is a grid of triangles defined by latitude and longitude
  !    lines over the surface of a sphere in 3D.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 May 2015
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) LAT_NUM, LONG_NUM, the number of latitude
  !    and longitude lines to draw.  The latitudes do not include the North and
  !    South poles, which will be included automatically, so LAT_NUM = 5, for
  !    instance, will result in points along 7 lines of latitude.
  !
  !    Output, integer(ip) POINT_NUM, the number of grid points.
  !

    integer(ip), intent(in), value :: lat_num
    integer(ip), intent(in), value :: long_num
    integer(ip), intent(out) :: point_num

    point_num = 2 + lat_num * long_num
  end subroutine sphere_llt_grid_point_count

  subroutine sphere_llt_grid_points ( r, pc, lat_num, long_num, point_num, p ) &
        bind(C, name="sphere_llt_grid_points")

  !*****************************************************************************80
  !
  !! SPHERE_LLT_GRID_POINTS: points for a latitude/longitude triangle grid.
  !
  !  Discussion:
  !
  !    A SPHERE LLT grid imposes a grid of triangles on a sphere,
  !    using latitude and longitude lines.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 May 2015
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real(dp) R, the radius of the sphere.
  !
  !    Input, real(dp) PC(3), the center of the sphere.
  !
  !    Input, integer(ip) LAT_NUM, LONG_NUM, the number of latitude
  !    and longitude lines to draw.  The latitudes do not include the North and
  !    South poles, which will be included automatically, so LAT_NUM = 5, for
  !    instance, will result in points along 7 lines of latitude.
  !
  !    Input, integer(ip) POINT_NUM, the number of points.
  !
  !    Output, real(dp) P(3,POINT_NUM), the grid points.
  !

    real(dp), intent(in), value :: r
    real(dp), intent(in) :: pc(3)
    integer(ip), intent(in), value :: lat_num
    integer(ip), intent(in), value :: long_num
    integer(ip), intent(in), value :: point_num
    real(dp), intent(out) :: p(3,point_num)

    integer(ip) :: lat
    integer(ip) :: lon
    integer(ip) :: n
    real(dp) :: phi
    real(dp), parameter :: r8_pi = 3.141592653589793_dp
    real(dp) :: theta

    n = 0
  !
  !  The north pole.
  !
    theta = 0.0_dp
    phi = 0.0_dp
    n = n + 1
    p(1,n) = pc(1) + r * sin ( phi ) * cos ( theta )
    p(2,n) = pc(2) + r * sin ( phi ) * sin ( theta )
    p(3,n) = pc(3) + r * cos ( phi )
  !
  !  Do each intermediate ring of latitude.
  !
    do lat = 1, lat_num

      phi = real ( lat, dp) * r8_pi &
          / real ( lat_num + 1, dp)
  !
  !  Along that ring of latitude, compute points at various longitudes.
  !
      do lon = 0, long_num - 1

        theta = real ( lon, dp) * 2.0_dp * r8_pi &
              / real ( long_num, dp)

        n = n + 1
        p(1,n) = pc(1) + r * sin ( phi ) * cos ( theta )
        p(2,n) = pc(2) + r * sin ( phi ) * sin ( theta )
        p(3,n) = pc(3) + r * cos ( phi )

      end do
    end do
  !
  !  The south pole.
  !
    theta = 0.0_dp
    phi = r8_pi
    n = n + 1
    p(1,n) = pc(1) + r * sin ( phi ) * cos ( theta )
    p(2,n) = pc(2) + r * sin ( phi ) * sin ( theta )
    p(3,n) = pc(3) + r * cos ( phi )
  end subroutine sphere_llt_grid_points

end module sphere_llt_grid_mod
