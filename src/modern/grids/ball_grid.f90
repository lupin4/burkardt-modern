!> ball_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module ball_grid_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  implicit none
  private

  public :: ball_grid, ball_grid_count

contains

  subroutine ball_grid ( n, r, c, ng, bg )

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
  !    Input, integer(int32) N, the number of subintervals.
  !
  !    Input, real(real64) R, the radius of the ball.
  !
  !    Input, real(real64) C(3), the coordinates of the center of the ball.
  !
  !    Input, integer(int32) NG, the number of grid points, as determined by
  !    BALL_GRID_COUNT.
  !
  !    Output, real(real64) BG(3,NG), the grid points inside the ball.
  !

    integer(int32) ng

    real(real64) bg(3,ng)
    real(real64) c(3)
    integer(int32) i
    integer(int32) j
    integer(int32) k
    integer(int32) n
    integer(int32) p
    real(real64) r
    real(real64) x
    real(real64) y
    real(real64) z

    p = 0

    do i = 0, n

      x = c(1) + r * real ( 2 * i, real64) / real ( 2 * n + 1, real64)

      do j = 0, n

        y = c(2) + r * real ( 2 * j, real64) / real ( 2 * n + 1, real64)

        do k = 0, n

          z = c(3) + r * real ( 2 * k, real64) / real ( 2 * n + 1, real64)

          if ( r * r < ( x - c(1) )**2 &
                     + ( y - c(2) )**2 &
                     + ( z - c(3) )**2 ) then
            exit
          end if

          p = p + 1
          bg(1,p) = x
          bg(2,p) = y
          bg(3,p) = z

          if ( 0 < i ) then
            p = p + 1
            bg(1,p) = 2.0e+00_real64 * c(1) - x
            bg(2,p) = y
            bg(3,p) = z
          end if

          if ( 0 < j ) then
            p = p + 1
            bg(1,p) = x
            bg(2,p) = 2.0e+00_real64 * c(2) - y
            bg(3,p) = z
          end if

          if ( 0 < k ) then
            p = p + 1
            bg(1,p) = x
            bg(2,p) = y
            bg(3,p) = 2.0e+00_real64 * c(3) - z
          end if

          if ( 0 < i .and. 0 < j ) then
            p = p + 1
            bg(1,p) = 2.0e+00_real64 * c(1) - x
            bg(2,p) = 2.0e+00_real64 * c(2) - y
            bg(3,p) = z
          end if

          if ( 0 < i .and. 0 < k ) then
            p = p + 1
            bg(1,p) = 2.0e+00_real64 * c(1) - x
            bg(2,p) = y
            bg(3,p) = 2.0e+00_real64 * c(3) - z
          end if

          if ( 0 < j .and. 0 < k ) then
            p = p + 1
            bg(1,p) = x
            bg(2,p) = 2.0e+00_real64 * c(2) - y
            bg(3,p) = 2.0e+00_real64 * c(3) - z
          end if

          if ( 0 < i .and. 0 < j .and. 0 < k ) then
            p = p + 1
            bg(1,p) = 2.0e+00_real64 * c(1) - x
            bg(2,p) = 2.0e+00_real64 * c(2) - y
            bg(3,p) = 2.0e+00_real64 * c(3) - z
          end if

        end do
      end do
    end do
  end

  subroutine ball_grid_count ( n, r, c, ng )

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
  !    Input, integer(int32) N, the number of subintervals.
  !
  !    Input, real(real64) R, the radius of the ball.
  !
  !    Input, real(real64) C(3), the coordinates of the center of the ball.
  !
  !    Output, integer(int32) NG, the number of grid points inside the ball.
  !

    real(real64) c(3)
    integer(int32) i
    integer(int32) j
    integer(int32) k
    integer(int32) n
    integer(int32) ng
    real(real64) r
    real(real64) x
    real(real64) y
    real(real64) z

    ng = 0

    do i = 0, n

      x = c(1) + r * real ( 2 * i, real64) / real ( 2 * n + 1, real64)

      do j = 0, n

        y = c(2) + r * real ( 2 * j, real64) / real ( 2 * n + 1, real64)

        do k = 0, n

          z = c(3) + r * real ( 2 * k, real64) / real ( 2 * n + 1, real64)

          if ( r * r < ( x - c(1) )**2 &
                     + ( y - c(2) )**2 &
                     + ( z - c(3) )**2 ) then
            exit
          end if

          ng = ng + 1

          if ( 0 < i ) then
            ng = ng + 1
          end if

          if ( 0 < j ) then
            ng = ng + 1
          end if

          if ( 0 < k ) then
            ng = ng + 1
          end if

          if ( 0 < i .and. 0 < j ) then
            ng = ng + 1
          end if

          if ( 0 < i .and. 0 < k ) then
            ng = ng + 1
          end if

          if ( 0 < j .and. 0 < k ) then
            ng = ng + 1
          end if

          if ( 0 < i .and. 0 < j .and. 0 < k ) then
            ng = ng + 1
          end if

        end do
      end do
    end do
  end

end module ball_grid_mod
