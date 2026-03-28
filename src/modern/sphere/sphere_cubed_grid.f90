!> sphere_cubed_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module sphere_cubed_grid_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  implicit none
  private

  public :: sphere_cubed_grid_ijk_to_xyz, sphere_cubed_grid_line_count, sphere_cubed_grid_lines, sphere_cubed_grid_lines_display, sphere_cubed_grid_point_count, sphere_cubed_grid_points
  public :: sphere_cubed_grid_points_display, sphere_cubed_grid_points_face

contains

  subroutine sphere_cubed_grid_ijk_to_xyz ( n, i, j, k, xyz )

  !*****************************************************************************80
  !
  !! SPHERE_CUBED_GRID_IJK_TO_XYZ: cubed sphere IJK to XYZ coordinates.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 May 2015
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of sections into which each 
  !    face of the cube is to be divided.
  !
  !    Input, integer(int32) I, J, K, indices between 0 and N.  Normally,
  !    at least one of the indices should have the value 0 or N.
  !
  !    Output, real(real64) XYZ(3), coordinates of the point.
  !

    integer(int32) i
    integer(int32) j
    integer(int32) k
    integer(int32) n
    real(real64), parameter :: r8_pi = 3.141592653589793e+00_real64
    real(real64) xc
    real(real64) xyz(3)
    real(real64) xyzn
    real(real64) yc
    real(real64) zc

    if ( i == 0 ) then
      xc = -1.0e+00_real64
    else if ( i == n ) then
      xc = +1.0e+00_real64
    else
      xc = tan ( real ( 2 * i - n, real64) * 0.25e+00_real64 * r8_pi &
        / real ( n, real64) )
    end if

    if ( j == 0 ) then
      yc = -1.0e+00_real64
    else if ( j == n ) then
      yc = +1.0e+00_real64
    else
      yc = tan ( real ( 2 * j - n, real64) * 0.25e+00_real64 * r8_pi &
        / real ( n, real64) )
    end if

    if ( k == 0 ) then
      zc = -1.0e+00_real64
    else if ( k == n ) then
      zc = +1.0e+00_real64
    else
      zc = tan ( real ( 2 * k - n, real64) * 0.25e+00_real64 * r8_pi &
        / real ( n, real64) )
    end if

    xyzn = sqrt ( xc ** 2 + yc ** 2 + zc ** 2 )

    xyz(1) = xc / xyzn
    xyz(2) = yc / xyzn
    xyz(3) = zc / xyzn
  end

  subroutine sphere_cubed_grid_line_count ( n, line_num )

  !*****************************************************************************80
  !
  !! SPHERE_CUBED_GRID_LINE_COUNT counts lines on a cubed sphere grid.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 May 2015
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of sections into which each 
  !    face of the cube is to be divided.
  !
  !    Output, integer(int32) LINE_NUM, the number of lines.
  !

    integer(int32) line_num
    integer(int32) n

    line_num = 0
  !
  !  If N = 1, the corners form 12 lines.
  !
    if ( n == 1 ) then
      line_num = 12
      return
  !
  !  If 1 < N, each of 8 corners connects to three neighboring edges.
  !
    else
      line_num = line_num + 8 * 3
    end if
  !
  !  If 2 < N, then each of the 12 edges includes lines.
  !
    if ( 2 < n ) then
      line_num = line_num + 12 * ( n - 2 )
    end if
  !
  !  Lines that belong to one of the six faces.
  !
    if ( 1 < n ) then
      line_num = line_num + 6 * 2 * n * ( n - 1 )
    end if
  end

  subroutine sphere_cubed_grid_lines ( n, line_num, line_data )

  !*****************************************************************************80
  !
  !! SPHERE_CUBED_GRID_LINES computes the lines on a cubed sphere grid.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 May 2015
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of sections into which each face 
  !    of the cube is to be divided.
  !
  !    Input, integer(int32) LINE_NUM, the number of lines.
  !
  !    Output, real(real64) LINE_DATA(3,2,LINE_NUM), for each line I, the 
  !    X/Y/Z coordinates of the start and end of a line segment on the grid.
  !

    integer(int32) line_num

    integer(int32) i
    integer(int32) j
    integer(int32) l
    real(real64) line_data(3,2,line_num)
    integer(int32) n

    l = 0
  !
  !  If N = 1, the corners form 12 lines.
  !
    if ( n == 1 ) then

      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, 0, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, 0, line_data(1:3,2,l) )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, 0, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n, 0, line_data(1:3,2,l) )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n, 0, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, 0, line_data(1:3,2,l) )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, 0, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, 0, line_data(1:3,2,l) )

      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, n, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, n, line_data(1:3,2,l) )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, n, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n, n, line_data(1:3,2,l) )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n, n, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, n, line_data(1:3,2,l) )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, n, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, n, line_data(1:3,2,l) )

      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, 0, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, n, line_data(1:3,2,l) )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, 0, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, n, line_data(1:3,2,l) )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n, 0, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n, n, line_data(1:3,2,l) )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, 0, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, n, line_data(1:3,2,l) )
      return
  !
  !  If 1 < N, each of 8 corners connects to three neighboring edges.
  !
    else

      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, 0, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, 1, 0, 0, line_data(1:3,2,l) )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, 0, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 1, 0, line_data(1:3,2,l) )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, 0, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, 1, line_data(1:3,2,l) )

      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n,   0, 0, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, n-1, 0, 0, line_data(1:3,2,l) )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, 0, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 1, 0, line_data(1:3,2,l) )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, 0, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, 1, line_data(1:3,2,l) )

      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n,   n, 0, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, n-1, n, 0, line_data(1:3,2,l) )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n,   0, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n-1, 0, line_data(1:3,2,l) )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n, 0, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n, 1, line_data(1:3,2,l) )

      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, 0, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, 1, n, 0, line_data(1:3,2,l) )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n,   0, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n-1, 0, line_data(1:3,2,l) )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, 0, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, 1, line_data(1:3,2,l) )

      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, n, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, 1, 0, n, line_data(1:3,2,l) )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, n, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 1, n, line_data(1:3,2,l) )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, n,   line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, n-1, line_data(1:3,2,l) )

      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n,   0, n, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, n-1, 0, n, line_data(1:3,2,l) )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, n, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 1, n, line_data(1:3,2,l) )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, n,   line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, n-1, line_data(1:3,2,l) )

      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n,   n, n, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, n-1, n, n, line_data(1:3,2,l) )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n,   n, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n-1, n, line_data(1:3,2,l) )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n, n,   line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n, n-1, line_data(1:3,2,l) )

      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, n, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, 1, n, n, line_data(1:3,2,l) )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n,   n, line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n-1, n, line_data(1:3,2,l) )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, n,   line_data(1:3,1,l) )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, n-1, line_data(1:3,2,l) )

    end if
  !
  !  If 2 < N, then each of the 12 edges includes lines.
  !
    if ( 2 < n ) then

      do i = 1, n - 2
        l = l + 1
        call sphere_cubed_grid_ijk_to_xyz ( n, i,   0, 0, line_data(1:3,1,l) )
        call sphere_cubed_grid_ijk_to_xyz ( n, i+1, 0, 0, line_data(1:3,2,l) )
      end do
      do i = 1, n - 2
        l = l + 1
        call sphere_cubed_grid_ijk_to_xyz ( n, n,   i, 0, line_data(1:3,1,l) )
        call sphere_cubed_grid_ijk_to_xyz ( n, n, i+1, 0, line_data(1:3,2,l) )
      end do
      do i = 1, n - 2
        l = l + 1
        call sphere_cubed_grid_ijk_to_xyz ( n, n-i,   n, 0, line_data(1:3,1,l) )
        call sphere_cubed_grid_ijk_to_xyz ( n, n-i-1, n, 0, line_data(1:3,2,l) )
      end do
      do i = 1, n - 2
        l = l + 1
        call sphere_cubed_grid_ijk_to_xyz ( n, 0, n-i,   0, line_data(1:3,1,l) )
        call sphere_cubed_grid_ijk_to_xyz ( n, 0, n-i-1, 0, line_data(1:3,2,l) )
      end do

      do i = 1, n - 2
        l = l + 1
        call sphere_cubed_grid_ijk_to_xyz ( n, i,   0, n, line_data(1:3,1,l) )
        call sphere_cubed_grid_ijk_to_xyz ( n, i+1, 0, n, line_data(1:3,2,l) )
      end do
      do i = 1, n - 2
        l = l + 1
        call sphere_cubed_grid_ijk_to_xyz ( n, n,   i, n, line_data(1:3,1,l) )
        call sphere_cubed_grid_ijk_to_xyz ( n, n, i+1, n, line_data(1:3,2,l) )
      end do
      do i = 1, n - 2
        l = l + 1
        call sphere_cubed_grid_ijk_to_xyz ( n, n-i,   n, n, line_data(1:3,1,l) )
        call sphere_cubed_grid_ijk_to_xyz ( n, n-i-1, n, n, line_data(1:3,2,l) )
      end do
      do i = 1, n - 2
        l = l + 1
        call sphere_cubed_grid_ijk_to_xyz ( n, 0, n-i,   n, line_data(1:3,1,l) )
        call sphere_cubed_grid_ijk_to_xyz ( n, 0, n-i-1, n, line_data(1:3,2,l) )
      end do

      do i = 1, n - 2
        l = l + 1
        call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, i,   line_data(1:3,1,l) )
        call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, i+1, line_data(1:3,2,l) )
      end do
      do i = 1, n - 2
        l = l + 1
        call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, i,   line_data(1:3,1,l) )
        call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, i+1, line_data(1:3,2,l) )
      end do
      do i = 1, n - 2
        l = l + 1
        call sphere_cubed_grid_ijk_to_xyz ( n, n, n, i,   line_data(1:3,1,l) )
        call sphere_cubed_grid_ijk_to_xyz ( n, n, n, i+1, line_data(1:3,2,l) )
      end do
      do i = 1, n - 2
        l = l + 1
        call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, i,   line_data(1:3,1,l) )
        call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, i+1, line_data(1:3,2,l) )
      end do

    end if
  !
  !  Lines that belong to one of the six faces.
  !
    if ( 1 < n ) then
  !
  !  000 : nn0
  !
      do i = 1, n - 1
        do j = 0, n - 1
          l = l + 1
          call sphere_cubed_grid_ijk_to_xyz ( n, i, j,   0, line_data(1:3,1,l) )
          call sphere_cubed_grid_ijk_to_xyz ( n, i, j+1, 0, line_data(1:3,2,l) )
        end do
      end do
      do j = 1, n - 1
        do i = 0, n - 1
          l = l + 1
          call sphere_cubed_grid_ijk_to_xyz ( n, i,   j, 0, line_data(1:3,1,l) )
          call sphere_cubed_grid_ijk_to_xyz ( n, i+1, j, 0, line_data(1:3,2,l) )
        end do
      end do
  !
  !  00n : nnn
  !
      do i = 1, n - 1
        do j = 0, n - 1
          l = l + 1
          call sphere_cubed_grid_ijk_to_xyz ( n, i, j,   n, line_data(1:3,1,l) )
          call sphere_cubed_grid_ijk_to_xyz ( n, i, j+1, n, line_data(1:3,2,l) )
        end do
      end do
      do j = 1, n - 1
        do i = 0, n - 1
          l = l + 1
          call sphere_cubed_grid_ijk_to_xyz ( n, i,   j, n, line_data(1:3,1,l) )
          call sphere_cubed_grid_ijk_to_xyz ( n, i+1, j, n, line_data(1:3,2,l) )
        end do
      end do
  !
  !  000:n0n
  !
      do i = 1, n - 1
        do j = 0, n - 1
          l = l + 1
          call sphere_cubed_grid_ijk_to_xyz ( n, i, 0, j,   line_data(1:3,1,l)   )
          call sphere_cubed_grid_ijk_to_xyz ( n, i, 0, j+1, line_data(1:3,2,l) )
        end do
      end do
      do j = 1, n - 1
        do i = 0, n - 1
          l = l + 1
          call sphere_cubed_grid_ijk_to_xyz ( n, i,   0, j, line_data(1:3,1,l) )
          call sphere_cubed_grid_ijk_to_xyz ( n, i+1, 0, j, line_data(1:3,2,l) )
        end do
      end do
  !
  !  0n0:nnn
  !
      do i = 1, n - 1
        do j = 0, n - 1
          l = l + 1
          call sphere_cubed_grid_ijk_to_xyz ( n, i, n, j,   line_data(1:3,1,l)   )
          call sphere_cubed_grid_ijk_to_xyz ( n, i, n, j+1, line_data(1:3,2,l) )
        end do
      end do
      do j = 1, n - 1
        do i = 0, n - 1
          l = l + 1
          call sphere_cubed_grid_ijk_to_xyz ( n, i,   n, j, line_data(1:3,1,l) )
          call sphere_cubed_grid_ijk_to_xyz ( n, i+1, n, j, line_data(1:3,2,l) )
        end do
      end do
  !
  !  000:0nn
  !
      do i = 1, n - 1
        do j = 0, n - 1
          l = l + 1
          call sphere_cubed_grid_ijk_to_xyz ( n, 0, i, j,   line_data(1:3,1,l)   )
          call sphere_cubed_grid_ijk_to_xyz ( n, 0, i, j+1, line_data(1:3,2,l) )
        end do
      end do
      do j = 1, n - 1
        do i = 0, n - 1
          l = l + 1
          call sphere_cubed_grid_ijk_to_xyz ( n, 0, i,   j, line_data(1:3,1,l) )
          call sphere_cubed_grid_ijk_to_xyz ( n, 0, i+1, j, line_data(1:3,2,l) )
        end do
      end do
  !
  !  n00:nnn
  !
      do i = 1, n - 1
        do j = 0, n - 1
          l = l + 1
          call sphere_cubed_grid_ijk_to_xyz ( n, n, i, j,   line_data(1:3,1,l)   )
          call sphere_cubed_grid_ijk_to_xyz ( n, n, i, j+1, line_data(1:3,2,l) )
        end do
      end do
      do j = 1, n - 1
        do i = 0, n - 1
          l = l + 1
          call sphere_cubed_grid_ijk_to_xyz ( n, n, i,   j, line_data(1:3,1,l) )
          call sphere_cubed_grid_ijk_to_xyz ( n, n, i+1, j, line_data(1:3,2,l) )
        end do
      end do

    end if

    if ( l /= line_num ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'SPHERE_CUBED_GRID_LINES - Fatal error!'
      write ( *, '(a,i6)' ) '  LINE_NUM = ', line_num
      write ( *, '(a,i6)' ) '  L = ', l
      stop 1
    end if
  end

  subroutine sphere_cubed_grid_lines_display ( line_num, line_data, prefix )

  !*****************************************************************************80
  !
  !! SPHERE_CUBED_GRID_LINES_DISPLAY displays a cubed grid on a sphere.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 May 2015
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
  !    Input, integer(int32) LINE_NUM, the number of grid lines.
  !
  !    Input, real(real64) LINE_DATA(3,2,LINE_NUM), contains pairs of 
  !    point indices for line segments that make up the grid.
  !
  !    Input, character ( len = * ) PREFIX, a prefix for the filenames.
  !

    integer(int32) line_num

    character ( len = 255 ) command_filename
    integer(int32) command_unit
    integer(int32) j
    integer(int32) l
    real(real64) line_data(3,2,line_num)
    character ( len = 255 ) line_filename
    integer(int32) line_unit
    character ( len = 255 ) plot_filename
    character ( len = * ) prefix
  !
  !  Create graphics data files.
  !
    call get_unit ( line_unit )
    line_filename = trim ( prefix ) // '_lines.txt'
    open ( unit = line_unit, file = line_filename, status = 'replace' )
    do l = 1, line_num
      if ( 1 < l ) then
        write ( line_unit, '(a)' ) ''
        write ( line_unit, '(a)' ) ''
      end if
      write ( line_unit, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) line_data(1:3,1,l)
      write ( line_unit, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) line_data(1:3,2,l)
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
      '" with lines lw 3'
    write ( command_unit, '(a)' ) 'quit'
    close ( unit = command_unit )

    write ( *, '(a)' ) &
      '  Created command file "' // trim ( command_filename ) // '".'
  end

  subroutine sphere_cubed_grid_point_count ( n, ns )

  !*****************************************************************************80
  !
  !! SPHERE_CUBED_POINT_COUNT counts the points on a cubed sphere grid.
  !
  !  Discussion:
  !
  !    For a value of N = 3, for instance, each of the 6 cube faces will
  !    be divided into 3 sections, so that a single cube face will have
  !    (3+1)x(3+1) points:
  !
  !      X---X---X---X
  !      | 1 | 4 | 7 |
  !      X---X---X---X
  !      | 2 | 5 | 8 |
  !      X---X---X---X
  !      | 3 | 6 | 9 |
  !      X---X---X---X
  !
  !    The number of points is simply (N+1)^3 - (N-1)^3.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 May 2015
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of sections into which 
  !    each face of the cube is to be divided.
  !
  !    Output, integer(int32) NS, the number of points.
  !

    integer(int32) n
    integer(int32) ns

    ns = ( n + 1 ) ** 3 - ( n - 1 ) ** 3
  end

  subroutine sphere_cubed_grid_points ( n, ns, xyz )

  !*****************************************************************************80
  !
  !! SPHERE_CUBED_GRID_POINTS computes the points on a cubed sphere grid.
  !
  !  Discussion:
  !
  !    For a value of N = 3, for instance, each of the 6 cube faces will
  !    be divided into 3 sections, so that a single cube face will have
  !    (3+1)x(3+1) points:
  !
  !      X---X---X---X
  !      | 1 | 4 | 7 |
  !      X---X---X---X
  !      | 2 | 5 | 8 |
  !      X---X---X---X
  !      | 3 | 6 | 9 |
  !      X---X---X---X
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 May 2015
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of sections into which each 
  !    face of the cube is to be divided.
  !
  !    Input, integer(int32) NS, the number of points.
  !
  !    Output, real(real64) XYZ(3,NS), distinct points on the unit sphere
  !    generated by a cubed sphere grid.
  !

    integer(int32) ns

    integer(int32) n
    integer(int32) ns2
    real(real64) xyz(3,ns)

    ns2 = 0
  !
  !  Bottom full.
  !
    call sphere_cubed_grid_points_face ( n, 0, 0, 0, n, n, 0, ns2, xyz )
  !
  !  To avoid repetition, draw the middles as grids of n-2 x n-1 points.
  !
    call sphere_cubed_grid_points_face ( n, 0, 0, 1, 0,   n-1, n-1, ns2, xyz )
    call sphere_cubed_grid_points_face ( n, 0, n, 1, n-1, n,   n-1, ns2, xyz )
    call sphere_cubed_grid_points_face ( n, n, 1, 1, n,   n,   n-1, ns2, xyz )
    call sphere_cubed_grid_points_face ( n, 1, 0, 1, n,   0,   n-1, ns2, xyz )
  !
  !  Top full.
  !
    call sphere_cubed_grid_points_face ( n, 0, 0, n, n, n, n, ns2, xyz )

    if ( ns2 /= ns ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPHERE_CUBED_GRID_POINTS - Fatal error!'
      write ( *, '(a,i8,a)' ) '  Expected to generated NS = ', ns, ' points.'
      write ( *, '(a,i8,a)' ) '  Generated ', ns2, ' points.'
      stop
    end if
  end

  subroutine sphere_cubed_grid_points_display ( ng, xg, prefix )

  !*****************************************************************************80
  !
  !! SPHERE_CUBED_GRID_POINTS_DISPLAY displays an LLQ grid on a sphere.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 May 2015
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

    integer(int32) line_num
    integer(int32) ng

    character ( len = 255 ) command_filename
    integer(int32) command_unit
    integer(int32) j
    character ( len = 255 ) node_filename
    integer(int32) node_unit
    character ( len = 255 ) plot_filename
    character ( len = * ) prefix
    real(real64) xg(3,ng)
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
      trim ( node_filename ) // '" with points pt 7 lt 0'
    write ( command_unit, '(a)' ) 'quit'
    close ( unit = command_unit )

    write ( *, '(a)' ) &
      '  Created command file "' // trim ( command_filename ) // '".'
  end

  subroutine sphere_cubed_grid_points_face ( n, i1, j1, k1, i2, j2, k2, ns, xyz )

  !*****************************************************************************80
  !
  !! SPHERE_CUBED_GRID_POINTS_FACE: points on one face of a cubed sphere grid.
  !
  !  Discussion:
  !
  !    This routine starts with NS = 0, and is called repeatedly to
  !    add points for another face.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 May 2015
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(int32) N, the number of sections into which each face 
  !    of the cube is to be divided.
  !
  !    Input, integer(int32) I1, J1, K1, I2, J2, K2, the logical indices, 
  !    between 0 and N, of two corners of the face grid.  It is guaranteed that 
  !    I1 <= I2, J1 <= J2, and K1 <= K2.  
  !
  !    Input/output, integer(int32) NS, the number of points.
  !
  !    Input/output, real XYZ(3,NS), distinct points on the unit sphere
  !    generated by a cubed sphere grid.
  !

    integer(int32) i
    integer(int32) i1
    integer(int32) i2
    integer(int32) j
    integer(int32) j1
    integer(int32) j2
    integer(int32) k
    integer(int32) k1
    integer(int32) k2
    integer(int32) n
    integer(int32) ns
    real(real64), parameter :: r8_pi = 3.141592653589793e+00_real64
    real(real64) xyz(3,*)
    real(real64) xyzn
    real(real64) xc
    real(real64) yc
    real(real64) zc

    do i = i1, i2

      if ( i1 < i2 ) then
        xc = tan ( real ( 2 * i - n, real64) * 0.25e+00_real64 * r8_pi &
          / real ( n, real64) )
      else if ( i1 == 0 ) then
        xc = -1.0e+00_real64
      else if ( i1 == n ) then
        xc = +1.0e+00_real64
      else
        xc = 0.0e+00_real64
      end if

      do j = j1, j2

        if ( j1 < j2 ) then
          yc = tan ( real ( 2 * j - n, real64) * 0.25e+00_real64 * r8_pi &
            / real ( n, real64) )
        else if ( j1 == 0 ) then
          yc = -1.0e+00_real64
        else if ( j1 == n ) then
          yc = +1.0e+00_real64
        else
          yc = 0.0e+00_real64
        end if

        do k = k1, k2

          if ( k1 < k2 ) then
            zc = tan ( real ( 2 * k - n, real64) * 0.25e+00_real64 * r8_pi &
              / real ( n, real64) )
          else if ( k1 == 0 ) then
            zc = -1.0e+00_real64
          else if ( k1 == n ) then
            zc = +1.0e+00_real64
          else
            zc = 0.0e+00_real64
          end if

          xyzn = sqrt ( xc ** 2 + yc ** 2 + zc ** 2 )

          ns = ns + 1
          xyz(1,ns) = xc / xyzn
          xyz(2,ns) = yc / xyzn
          xyz(3,ns) = zc / xyzn

        end do
      end do
    end do
  end

end module sphere_cubed_grid_mod
