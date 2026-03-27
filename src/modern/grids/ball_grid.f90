!> ball_grid -- Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module ball_grid_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: ball_grid, ball_grid_count

contains

  subroutine ball_grid ( n, r, c, ng, bg ) &
        bind(C, name="ball_grid")

  !*****************************************************************************80
  !
  !! BALL_GRID computes grid points inside a ball.
  !
  !  Discussion:
  !
  !    The grid is defined by specifying the radius and center of the ball,
  !    and the number of subintervals N into which the horizontal radius
  !    should be divided.  Thus, a value of N = 2 will result in 5 points
  !    along that horizontal line.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    10 November 2011
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of subintervals.
  !
  !    Input, real(dp) R, the radius of the ball.
  !
  !    Input, real(dp) C(3), the coordinates of the center of the ball.
  !
  !    Input, integer(ip) NG, the number of grid points, as determined by
  !    BALL_GRID_COUNT.
  !
  !    Output, real(dp) BG(3,NG), the grid points inside the ball.
  !

    integer(ip), intent(in), value :: ng
    real(dp), intent(out) :: bg(3,ng)
    real(dp), intent(in) :: c(3)
    integer(ip) :: i
    integer(ip) :: j
    integer(ip) :: k
    integer(ip), intent(in), value :: n
    integer(ip) :: p
    real(dp), intent(in), value :: r
    real(dp) :: x
    real(dp) :: y
    real(dp) :: z

    p = 0

    do i = 0, n

      x = c(1) + r * real ( 2 * i, dp) / real ( 2 * n + 1, dp)

      do j = 0, n

        y = c(2) + r * real ( 2 * j, dp) / real ( 2 * n + 1, dp)

        do k = 0, n

          z = c(3) + r * real ( 2 * k, dp) / real ( 2 * n + 1, dp)

          if ( r * r < ( x - c(1) )**2 &
                     + ( y - c(2) )**2 &
                     + ( z - c(3) )**2 ) then
            exit
          end if

          p = p + 1
          bg(1,p) = x
          bg(2,p) = y
          bg(3,p) = z

          if ( i > 0 ) then
            p = p + 1
            bg(1,p) = 2.0_dp * c(1) - x
            bg(2,p) = y
            bg(3,p) = z
          end if

          if ( j > 0 ) then
            p = p + 1
            bg(1,p) = x
            bg(2,p) = 2.0_dp * c(2) - y
            bg(3,p) = z
          end if

          if ( k > 0 ) then
            p = p + 1
            bg(1,p) = x
            bg(2,p) = y
            bg(3,p) = 2.0_dp * c(3) - z
          end if

          if ( i > 0 .and. j > 0 ) then
            p = p + 1
            bg(1,p) = 2.0_dp * c(1) - x
            bg(2,p) = 2.0_dp * c(2) - y
            bg(3,p) = z
          end if

          if ( i > 0 .and. k > 0 ) then
            p = p + 1
            bg(1,p) = 2.0_dp * c(1) - x
            bg(2,p) = y
            bg(3,p) = 2.0_dp * c(3) - z
          end if

          if ( j > 0 .and. k > 0 ) then
            p = p + 1
            bg(1,p) = x
            bg(2,p) = 2.0_dp * c(2) - y
            bg(3,p) = 2.0_dp * c(3) - z
          end if

          if ( i > 0 .and. j > 0 .and. k > 0 ) then
            p = p + 1
            bg(1,p) = 2.0_dp * c(1) - x
            bg(2,p) = 2.0_dp * c(2) - y
            bg(3,p) = 2.0_dp * c(3) - z
          end if

        end do
      end do
    end do
  end subroutine ball_grid

  subroutine ball_grid_count ( n, r, c, ng ) &
        bind(C, name="ball_grid_count")

  !*****************************************************************************80
  !
  !! BALL_GRID computes grid points inside a ball.
  !
  !  Discussion:
  !
  !    The grid is defined by specifying the radius and center of the ball,
  !    and the number of subintervals N into which the horizontal radius
  !    should be divided.  Thus, a value of N = 2 will result in 5 points
  !    along that horizontal line.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    10 November 2011
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of subintervals.
  !
  !    Input, real(dp) R, the radius of the ball.
  !
  !    Input, real(dp) C(3), the coordinates of the center of the ball.
  !
  !    Output, integer(ip) NG, the number of grid points inside the ball.
  !

    real(dp), intent(in) :: c(3)
    integer(ip) :: i
    integer(ip) :: j
    integer(ip) :: k
    integer(ip), intent(in), value :: n
    integer(ip), intent(out) :: ng
    real(dp), intent(in), value :: r
    real(dp) :: x
    real(dp) :: y
    real(dp) :: z

    ng = 0

    do i = 0, n

      x = c(1) + r * real ( 2 * i, dp) / real ( 2 * n + 1, dp)

      do j = 0, n

        y = c(2) + r * real ( 2 * j, dp) / real ( 2 * n + 1, dp)

        do k = 0, n

          z = c(3) + r * real ( 2 * k, dp) / real ( 2 * n + 1, dp)

          if ( r * r < ( x - c(1) )**2 &
                     + ( y - c(2) )**2 &
                     + ( z - c(3) )**2 ) then
            exit
          end if

          ng = ng + 1

          if ( i > 0 ) then
            ng = ng + 1
          end if

          if ( j > 0 ) then
            ng = ng + 1
          end if

          if ( k > 0 ) then
            ng = ng + 1
          end if

          if ( i > 0 .and. j > 0 ) then
            ng = ng + 1
          end if

          if ( i > 0 .and. k > 0 ) then
            ng = ng + 1
          end if

          if ( j > 0 .and. k > 0 ) then
            ng = ng + 1
          end if

          if ( i > 0 .and. j > 0 .and. k > 0 ) then
            ng = ng + 1
          end if

        end do
      end do
    end do
  end subroutine ball_grid_count

end module ball_grid_mod
