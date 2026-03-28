!> sphere_fibonacci_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module sphere_fibonacci_grid_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  implicit none
  private

  public :: sphere_fibonacci_grid_display, sphere_fibonacci_grid_points

contains

  subroutine sphere_fibonacci_grid_display ( ng, xg, prefix )

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
  !    Input, integer(int32) NG, the number of points.
  !
  !    Input, real(real64) XG(3,NG), the Fibonacci spiral points.
  !
  !    Input, character ( len = * ) PREFIX, a prefix for the filenames.
  !

    integer(int32) ng

    character ( len = 255 ) command_filename
    integer(int32) command_unit
    character ( len = 255 ) data_filename
    integer(int32) data_unit
    integer(int32) j
    character ( len = 255 ) plot_filename
    character ( len = * ) prefix
    real(real64) xg(3,ng)
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
  end

  subroutine sphere_fibonacci_grid_points ( ng, xyz )

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
  !    Input, integer(int32) NG, the number of points.
  !
  !    Output, real(real64) XYZ(3,NG), the Fibonacci spiral points.
  !

    integer(int32) ng

    real(real64) cphi
    integer(int32) i
    real(real64) i_r8
    integer(int32) j
    real(real64) ng_r8
    real(real64) r8_phi
    real(real64), parameter :: r8_pi = 3.141592653589793e+00_real64
    real(real64) sphi
    real(real64) theta
    real(real64) xyz(3,ng)

    r8_phi = ( 1.0e+00_real64 + sqrt ( 5.0e+00_real64 ) ) / 2.0e+00_real64
    ng_r8 = real ( ng, real64)

    do j = 1, ng
      i_r8 = real ( - ng - 1 + 2 * j, real64)
      theta = 2.0e+00_real64 * r8_pi * i_r8 / r8_phi
      sphi = i_r8 / ng_r8
      cphi = sqrt ( ( ng_r8 + i_r8 ) * ( ng_r8 - i_r8 ) ) / ng_r8
      xyz(1,j) = cphi * sin ( theta )
      xyz(2,j) = cphi * cos ( theta )
      xyz(3,j) = sphi
    end do
  end

end module sphere_fibonacci_grid_mod
