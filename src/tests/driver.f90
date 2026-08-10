program regression_driver
!! burkardt-modern test driver.
!!
!! Serves TWO purposes at once:
!!
!!   1. DIFFERENTIAL. run_tests.sh builds this same source against the ORIGINAL
!!      objects and against the MODERN objects and diffs stdout. Any numeric
!!      drift introduced by the modernization fails the run.
!!
!!   2. ANALYTIC. Every value below is also checked against a closed form that
!!      is independent of the implementation (pi*r^2, 4/3 pi r^3, gcd, 5!, the
!!      3-4-5 triangle, ...). A pure diff can only prove the two trees AGREE --
!!      it cannot catch a bug they SHARE. The analytic layer can, and sets a
!!      non-zero exit status when a closed form is violated.
!!
!! Written in Fortran 90 style deliberately: this file must compile clean under
!! both -std=legacy (original tree) and -std=f2018 (modern tree).
!!
!! Copyright The Fantastic Planet - By David Clabaugh
  implicit none

  real ( kind = 8 ), parameter :: PI = 3.141592653589793D+00
  real ( kind = 8 ), parameter :: TOL = 1.0D-12

  ! --- externals from core/geometry.f90 ---
  real ( kind = 8 ) r8_asin, sphere01_area_nd, sphere01_volume_nd
  real ( kind = 8 ) torus_area_3d, hexagon_area_2d, hexagon01_area_2d
  real ( kind = 8 ) radians_to_degrees, degrees_to_radians, r8vec_norm
  integer ( kind = 4 ) i4_gcd
  external triangle_area_2d, triangle_area_3d, parallelogram_area_2d
  external tetrahedron_volume_3d, polygon_area_2d, line_exp_point_dist_2d
  external triangle_centroid_2d, circle_imp_point_dist_2d

  real ( kind = 8 ) :: t2(2,3), t3(3,3), tet(3,4), quad(2,4), sq(2,4)
  real ( kind = 8 ) :: p1(2), p2(2), pq(2), pc(2), cen(2)
  real ( kind = 8 ) :: x, area, vol, dist, r
  integer ( kind = 4 ) :: g, n, nfail

  nfail = 0

  write (*,'(a)') '--- burkardt-modern regression + analytic driver ---'

! ---------------------------------------------------------------- areas
  t2 = reshape([0.0D+00,0.0D+00, 4.0D+00,0.0D+00, 0.0D+00,3.0D+00],[2,3])
  call triangle_area_2d ( t2, area )
  call chk ( 'triangle_area_2d      ', area, 6.0D+00, nfail )

  t3 = reshape([0.0D+00,0.0D+00,0.0D+00, 1.0D+00,0.0D+00,0.0D+00, &
                0.0D+00,1.0D+00,0.0D+00],[3,3])
  call triangle_area_3d ( t3, area )
  call chk ( 'triangle_area_3d      ', area, 0.5D+00, nfail )

  quad = reshape([0.0D+00,0.0D+00, 1.0D+00,0.0D+00, 1.0D+00,1.0D+00, &
                  0.0D+00,1.0D+00],[2,4])
  call parallelogram_area_2d ( quad, area )
  call chk ( 'parallelogram_area_2d ', area, 1.0D+00, nfail )

  sq = reshape([0.0D+00,0.0D+00, 1.0D+00,0.0D+00, 1.0D+00,1.0D+00, &
                0.0D+00,1.0D+00],[2,4])
  call polygon_area_2d ( 4, sq, area )
  call chk ( 'polygon_area_2d       ', area, 1.0D+00, nfail )

  tet = reshape([0.0D+00,0.0D+00,0.0D+00, 1.0D+00,0.0D+00,0.0D+00, &
                 0.0D+00,1.0D+00,0.0D+00, 0.0D+00,0.0D+00,1.0D+00],[3,4])
  call tetrahedron_volume_3d ( tet, vol )
  call chk ( 'tetrahedron_volume_3d ', vol, 1.0D+00/6.0D+00, nfail )

  call chk ( 'hexagon01_area_2d     ', hexagon01_area_2d ( ), &
             3.0D+00*sqrt(3.0D+00)/2.0D+00, nfail )
  call chk ( 'hexagon_area_2d(2)    ', hexagon_area_2d ( 2.0D+00 ), &
             4.0D+00*3.0D+00*sqrt(3.0D+00)/2.0D+00, nfail )

! ------------------------------------------------------- unit hypersphere
! surface area of S^(n-1) = 2 pi^(n/2) / Gamma(n/2); volume = that / n
  call chk ( 'sphere01_area_nd(2)   ', sphere01_area_nd(2), 2.0D+00*PI, nfail )
  call chk ( 'sphere01_area_nd(3)   ', sphere01_area_nd(3), 4.0D+00*PI, nfail )
  call chk ( 'sphere01_area_nd(4)   ', sphere01_area_nd(4), 2.0D+00*PI*PI, nfail )
  call chk ( 'sphere01_area_nd(5)   ', sphere01_area_nd(5), 8.0D+00*PI*PI/3.0D+00, nfail )

  call chk ( 'sphere01_volume_nd(2) ', sphere01_volume_nd(2), PI, nfail )
  call chk ( 'sphere01_volume_nd(3) ', sphere01_volume_nd(3), 4.0D+00*PI/3.0D+00, nfail )
  call chk ( 'sphere01_volume_nd(4) ', sphere01_volume_nd(4), PI*PI/2.0D+00, nfail )
  call chk ( 'sphere01_volume_nd(5) ', sphere01_volume_nd(5), 8.0D+00*PI*PI/15.0D+00, nfail )

! area/volume identity: V(n) = A(n) / n   (unit ball)
  do n = 2, 5
    call chk ( 'ball identity V=A/n   ', sphere01_volume_nd(n), &
               sphere01_area_nd(n)/real(n,kind=8), nfail )
  end do

  call chk ( 'torus_area_3d(3,1)    ', torus_area_3d(3.0D+00,1.0D+00), &
             4.0D+00*PI*PI*3.0D+00, nfail )

! ------------------------------------------------------------ scalar math
  call chk ( 'r8_asin(0.5)          ', r8_asin(0.5D+00), PI/6.0D+00, nfail )
  call chk ( 'r8_asin(clamp>1)      ', r8_asin(1.0000001D+00), PI/2.0D+00, nfail )
  call chk ( 'r8_asin(clamp<-1)     ', r8_asin(-1.0000001D+00), -PI/2.0D+00, nfail )

  g = i4_gcd ( 1071, 462 )
  call chki ( 'i4_gcd(1071,462)      ', g, 21, nfail )
  call chki ( 'i4_gcd(17,5)          ', i4_gcd(17,5), 1, nfail )

  call chk ( 'degrees_to_radians    ', degrees_to_radians(180.0D+00), PI, nfail )
  call chk ( 'radians_to_degrees    ', radians_to_degrees(PI), 180.0D+00, nfail )
! round trip must be the identity
  call chk ( 'deg/rad round trip    ', radians_to_degrees(degrees_to_radians(37.5D+00)), &
             37.5D+00, nfail )

  call chk ( 'r8vec_norm(3-4-5)     ', r8vec_norm(3,[3.0D+00,4.0D+00,0.0D+00]), &
             5.0D+00, nfail )

! -------------------------------------------------------------- distances
  p1 = [0.0D+00, 0.0D+00]
  p2 = [1.0D+00, 0.0D+00]
  pq = [0.5D+00, 2.0D+00]
  call line_exp_point_dist_2d ( p1, p2, pq, dist )
  call chk ( 'line_exp_point_dist_2d', dist, 2.0D+00, nfail )

  r  = 1.0D+00
  pc = [0.0D+00, 0.0D+00]
  pq = [3.0D+00, 0.0D+00]
  call circle_imp_point_dist_2d ( r, pc, pq, dist )
  call chk ( 'circle_imp_point_dist ', dist, 2.0D+00, nfail )

! -------------------------------------------------------------- centroid
  call triangle_centroid_2d ( t2, cen )
  call chk ( 'triangle_centroid x   ', cen(1), 4.0D+00/3.0D+00, nfail )
  call chk ( 'triangle_centroid y   ', cen(2), 1.0D+00, nfail )

! ------------------------------------------------------------------ done
  write (*,'(a,i0,a)') '--- analytic failures: ', nfail, ' ---'
  if ( nfail > 0 ) then
    write (*,'(a)') 'ANALYTIC FAIL'
    stop 1
  end if

contains

  subroutine chk ( name, got, want, nf )
    character(len=*), intent(in) :: name
    real ( kind = 8 ), intent(in) :: got, want
    integer ( kind = 4 ), intent(inout) :: nf
    real ( kind = 8 ) :: scale
    scale = max ( 1.0D+00, abs ( want ) )
    ! the value is printed unconditionally so the ORIGINAL-vs-MODERN diff sees it
    write (*,'(a,a,f22.12)') name, ': ', got
    if ( abs ( got - want ) > TOL * scale ) then
      write (*,'(a,a,f22.12)') name, ':   EXPECTED ', want
      nf = nf + 1
    end if
  end subroutine chk

  subroutine chki ( name, got, want, nf )
    character(len=*), intent(in) :: name
    integer ( kind = 4 ), intent(in) :: got, want
    integer ( kind = 4 ), intent(inout) :: nf
    write (*,'(a,a,i12)') name, ': ', got
    if ( got /= want ) then
      write (*,'(a,a,i12)') name, ':   EXPECTED ', want
      nf = nf + 1
    end if
  end subroutine chki

end program regression_driver
