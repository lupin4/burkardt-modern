!> sphere_fibonacci_grid -- Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module sphere_fibonacci_grid_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: sphere_fibonacci_grid_display, sphere_fibonacci_grid_points

contains

  subroutine sphere_fibonacci_grid_display ( ng, xg, prefix ) &
        bind(C, name="sphere_fibonacci_grid_display")

  !*****************************************************************************80
  !
  !! SPHERE_FIBONACCI_GRID_DISPLAY displays sphere points on a Fibonacci spiral.
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
  !    Input, integer(ip) NG, the number of points.
  !
  !    Input, real(dp) XG(3,NG), the Fibonacci spiral points.
  !
  !    Input, character ( len = * ) PREFIX, a prefix for the filenames.
  !

    integer(ip), intent(in), value :: ng
    real(dp), intent(in) :: xg(3,ng)
    character ( len = * ), intent(in) :: prefix

    character ( len = 255 ) :: command_filename
    integer(ip) :: command_unit
    character ( len = 255 ) :: data_filename
    integer(ip) :: data_unit
    integer(ip) :: j
    character ( len = 255 ) :: plot_filename
  !
  !  Create graphics data files.
  !
    call get_unit ( data_unit )
    data_filename = trim ( prefix ) // '_data.txt'
    open ( unit = data_unit, file = data_filename, status = 'replace' )
    do j = 1, ng
      write ( data_unit, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) xg(1:3,j)
    end do
    close ( unit = data_unit )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Created data file "' // trim ( data_filename ) // '".'
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
    write ( command_unit, '(a)' ) 'splot "' // trim ( data_filename ) // '"'
    write ( command_unit, '(a)' ) 'quit'
    close ( unit = command_unit )

    write ( *, '(a)' ) &
      '  Created command file "' // trim ( command_filename ) // '".'
  end subroutine sphere_fibonacci_grid_display

  subroutine sphere_fibonacci_grid_points ( ng, xyz ) &
        bind(C, name="sphere_fibonacci_grid_points")

  !*****************************************************************************80
  !
  !! SPHERE_FIBONACCI_GRID_POINTS computes sphere points on a Fibonacci spiral.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 October 2013
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Richard Swinbank, James Purser,
  !    Fibonacci grids: A novel approach to global modelling,
  !    Quarterly Journal of the Royal Meteorological Society,
  !    Volume 132, Number 619, July 2006 Part B, pages 1769-1793.
  !
  !  Parameters:
  !
  !    Input, integer(ip) NG, the number of points.
  !
  !    Output, real(dp) XYZ(3,NG), the Fibonacci spiral points.
  !

    integer(ip), intent(in), value :: ng
    real(dp), intent(out) :: xyz(3,ng)

    real(dp) :: cphi
    integer(ip) :: j
    real(dp) :: i_r8
    real(dp) :: ng_r8
    real(dp) :: r8_phi
    real(dp), parameter :: r8_pi = 3.141592653589793_dp
    real(dp) :: sphi
    real(dp) :: theta

    r8_phi = ( 1.0_dp + sqrt ( 5.0_dp ) ) / 2.0_dp
    ng_r8 = real ( ng, dp)

    do j = 1, ng
      i_r8 = real ( - ng - 1 + 2 * j, dp)
      theta = 2.0_dp * r8_pi * i_r8 / r8_phi
      sphi = i_r8 / ng_r8
      cphi = sqrt ( ( ng_r8 + i_r8 ) * ( ng_r8 - i_r8 ) ) / ng_r8
      xyz(1,j) = cphi * sin ( theta )
      xyz(2,j) = cphi * cos ( theta )
      xyz(3,j) = sphi
    end do
  end subroutine sphere_fibonacci_grid_points

end module sphere_fibonacci_grid_mod
