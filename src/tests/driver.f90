program regression_driver
!! burkardt-modern regression driver — one binary is built against the
!! ORIGINAL objects and one against the MODERN objects; run_tests.sh
!! diffs their stdout. Any numeric drift between the trees fails the test.
!! Copyright The Fantastic Planet - By David Clabaugh
  implicit none
  real ( kind = 8 ) :: t(2,3), area
  real ( kind = 8 ) :: r8_asin, sphere01_area_nd
  real ( kind = 8 ) :: x
  integer ( kind = 4 ) :: i4_gcd, g, n
  external triangle_area_2d

  ! triangle area (3-4-5 right triangle halves)
  t = reshape([0.0d0, 0.0d0, 4.0d0, 0.0d0, 0.0d0, 3.0d0], [2,3])
  call triangle_area_2d ( t, area )
  write (*, '(a,f20.12)') 'triangle_area_2d: ', area

  ! r8_asin at a mid value + clamped edge
  x = r8_asin ( 0.5d0 )
  write (*, '(a,f20.12)') 'r8_asin(0.5):    ', x
  x = r8_asin ( 1.0000001d0 )
  write (*, '(a,f20.12)') 'r8_asin(clamp):  ', x

  ! i4_gcd
  g = i4_gcd ( 1071, 462 )
  write (*, '(a,i12)') 'i4_gcd(1071,462): ', g

  ! unit-sphere surface areas for a few dimensions
  do n = 2, 5
    x = sphere01_area_nd ( n )
    write (*, '(a,i2,a,f20.12)') 'sphere01_area_nd(', n, '): ', x
  end do
end program regression_driver
