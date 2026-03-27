!> hypersphere_properties -- Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module hypersphere_properties_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: cartesian_to_hypersphere, hypersphere_01_area, hypersphere_01_area_values, hypersphere_01_interior_uniform, hypersphere_01_surface_uniform, hypersphere_01_volume
  public :: hypersphere_01_volume_values, hypersphere_area, hypersphere_stereograph, hypersphere_stereograph_inverse, hypersphere_surface_uniform, hypersphere_to_cartesian
  public :: hypersphere_volume, r8_uniform_01, r8mat_norm_fro_affine, r8mat_normal_01, r8mat_uniform_01, r8vec_normal_01
  public :: r8vec_uniform_01, sphere_stereograph, sphere_stereograph_inverse

contains

  pure subroutine cartesian_to_hypersphere ( m, n, c, x, r, theta ) &
        bind(C, name="cartesian_to_hypersphere")

  !*****************************************************************************80
  !
  !! CARTESIAN_TO_HYPERSPHERE: Cartesian to hypersphere coordinate transform.
  !
  !  Discussion:
  !
  !    We allow the trivial case M = 1; in that case alone, the value R
  !    must be assumed to have a sign.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 May 2013
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) M, the spatial dimension.
  !    1 <= M.
  !
  !    Input, integer(ip) N, the number of points to transform.
  !
  !    Input, real(dp) C(M), the center of the hypersphere.
  !
  !    Input, real(dp) X(M,N), the Cartesian coordinates of the points.
  !
  !    Output, real(dp) R(N), the radius of the points on the
  !    hypersphere.  Except for the trivial case M = 1, R is assumed nonnegative.
  !
  !    Output, real(dp) THETA(M-1,N), the coordinate angles of the
  !    points, measured in radians.
  !

    integer(ip), intent(in), value :: m
    integer(ip), intent(in), value :: n
    real(dp), intent(in) :: c(m)
    real(dp), intent(in) :: x(m,n)
    real(dp), intent(out) :: r(n)
    real(dp), intent(out) :: theta(m-1,n)
    integer(ip) :: i
    integer(ip) :: i1
    integer(ip) :: j
    real(dp) :: top(n)
    real(dp) :: x2(m,n)
  !
  !  Handle special case of M = 1.
  !
    if ( m == 1 ) then
      r(1:n) = x(1,1:n) - c(1)
    end if
  !
  !  Subtract the center.
  !
    do i = 1, m
      x2(i,1:n) = x(i,1:n) - c(i)
    end do
  !
  !  Compute R.
  !
    do j = 1, n
      r(j) = sqrt ( sum ( x2(1:m,j)**2 ) )
    end do
  !
  !  Compute M-2 components of THETA.
  !
    theta(1:m-1,n) = 0.0_dp

    do i = 2, m - 1
      do i1 = 1, i - 1
        theta(i1,1:n) = theta(i1,1:n) + x2(i,1:n) ** 2
      end do
    end do

    do i = 1, m - 2
      theta(i,1:n) = theta(i,1:n) + x2(m,1:n) ** 2
    end do

    do i = 1, m - 2
      theta(i,1:n) = atan2 ( sqrt ( theta(i,1:n) ), x2(i,1:n) )
    end do
  !
  !  Compute last component of THETA.
  !
    top(1:n) = sqrt ( x2(m,1:n) ** 2 + x2(m-1,1:n) ** 2 ) + x2(m-1,1:n)

    do j = 1, n
      theta(m-1,j) = 2.0_dp * atan2 ( x2(m,j), top(j) )
    end do
  end subroutine cartesian_to_hypersphere

  pure function hypersphere_01_area ( m ) &
        bind(C, name="hypersphere_01_area")

  !*****************************************************************************80
  !
  !! HYPERSPHERE_01_AREA computes the surface area of a hypersphere.
  !
  !  Discussion:
  !
  !    The unit hypersphere satisfies:
  !
  !      sum ( 1 <= I <= M ) X(I) * X(I) = 1
  !
  !    Results include:
  !
  !     M   Area
  !
  !     2    2        * PI
  !     3    4        * PI
  !     4  ( 2 /   1) * PI^2
  !     5  ( 8 /   3) * PI^2
  !     6  ( 1 /   1) * PI^3
  !     7  (16 /  15) * PI^3
  !     8  ( 1 /   3) * PI^4
  !     9  (32 / 105) * PI^4
  !    10  ( 1 /  12) * PI^5
  !
  !    For the unit hypersphere, Area(M) = M * Volume(M)
  !
  !    Sphere_Unit_Area ( M ) = 2 * PI^(M/2) / Gamma ( M / 2 )
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 September 2003
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) M, the dimension of the space.
  !
  !    Output, real(dp) HYPERSPHERE_01_AREA, the area.
  !

    real(dp) :: hypersphere_01_area
    integer(ip), intent(in), value :: m
    real(dp) :: area
    integer(ip) :: i
    integer(ip) :: m2
    real(dp), parameter :: r8_pi = 3.141592653589793_dp

    if ( mod ( m, 2 ) == 0 ) then
      m2 = m / 2
      area = 2.0_dp * ( r8_pi ) ** m2
      do i = 1, m2 - 1
        area = area / real ( i, dp)
      end do
    else
      m2 = ( m - 1 ) / 2
      area = ( r8_pi ** m2 ) * ( 2.0_dp ** m )
      do i = m2 + 1, 2 * m2
        area = area / real ( i, dp)
      end do
    end if

    hypersphere_01_area = area
  end function hypersphere_01_area

  subroutine hypersphere_01_area_values ( n_data, n, area ) &
        bind(C, name="hypersphere_01_area_values")

  !*****************************************************************************80
  !
  !! HYPERSPHERE_01_AREA_VALUES returns some areas of the unit hypersphere.
  !
  !  Discussion:
  !
  !    The formula for the surface area of the unit hypersphere is:
  !
  !      Sphere_Unit_Area ( N ) = 2 * pi^(N/2) / Gamma ( N / 2 )
  !
  !    Some values of the function include:
  !
  !       N   Area
  !
  !       2    2        * PI
  !       3  ( 4 /    ) * PI
  !       4  ( 2 /   1) * PI^2
  !       5  ( 8 /   3) * PI^2
  !       6  ( 1 /   1) * PI^3
  !       7  (16 /  15) * PI^3
  !       8  ( 1 /   3) * PI^4
  !       9  (32 / 105) * PI^4
  !      10  ( 1 /  12) * PI^5
  !
  !    For the unit hypersphere, Area(N) = N * Volume(N)
  !
  !    In Mathematica, the function can be evaluated by:
  !
  !      2 * Pi^(n/2) / Gamma[n/2]
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    20 August 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Stephen Wolfram,
  !    The Mathematica Book,
  !    Fourth Edition,
  !    Cambridge University Press, 1999,
  !    ISBN: 0-521-64314-7,
  !    LC: QA76.95.W65.
  !
  !  Parameters:
  !
  !    Input/output, integer(ip) N_DATA.  The user sets N_DATA to 0
  !    before the first call.  On each call, the routine increments N_DATA by 1,
  !    and returns the corresponding data; when there is no more data, the
  !    output value of N_DATA will be 0 again.
  !
  !    Output, integer(ip) N, the spatial dimension.
  !
  !    Output, real(dp) AREA, the area
  !    in that dimension.
  !

    integer(ip), parameter :: n_max = 20
    integer(ip), intent(inout) :: n_data
    integer(ip), intent(out) :: n
    real(dp), intent(out) :: area
    real(dp), save, dimension ( n_max ) :: area_vec = (/ &
      0.2000000000000000e+01_dp, &
      0.6283185307179586e+01_dp, &
      0.1256637061435917e+02_dp, &
      0.1973920880217872e+02_dp, &
      0.2631894506957162e+02_dp, &
      0.3100627668029982e+02_dp, &
      0.3307336179231981e+02_dp, &
      0.3246969701133415e+02_dp, &
      0.2968658012464836e+02_dp, &
      0.2550164039877345e+02_dp, &
      0.2072514267328890e+02_dp, &
      0.1602315322625507e+02_dp, &
      0.1183817381218268e+02_dp, &
      0.8389703410491089e+01_dp, &
      0.5721649212349567e+01_dp, &
      0.3765290085742291e+01_dp, &
      0.2396678817591364e+01_dp, &
      0.1478625959000308e+01_dp, &
      0.8858104195716824e+00_dp, &
      0.5161378278002812e+00_dp /)
    integer(ip), save, dimension ( n_max ) :: n_vec = (/ &
       1, &
       2, &
       3, &
       4, &
       5, &
       6, &
       7, &
       8, &
       9, &
      10, &
      11, &
      12, &
      13, &
      14, &
      15, &
      16, &
      17, &
      18, &
      19, &
      20 /)

    if ( n_data < 0 ) then
      n_data = 0
    end if

    n_data = n_data + 1

    if ( n_max < n_data ) then
      n_data = 0
      n = 0
      area = 0.0_dp
    else
      n = n_vec(n_data)
      area = area_vec(n_data)
    end if
  end subroutine hypersphere_01_area_values

  subroutine hypersphere_01_interior_uniform ( m, n, seed, x ) &
        bind(C, name="hypersphere_01_interior_uniform")

  !*****************************************************************************80
  !
  !! HYPERSPHERE_01_INTERIOR_UNIFORM: uniform points inside unit hypersphere.
  !
  !  Discussion:
  !
  !    The hypersphere has center 0 and radius 1.
  !
  !    This routine is valid for any spatial dimension.
  !
  !    We first generate a point ON the hypersphere, and then distribute it
  !    IN the hypersphere.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 August 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Russell Cheng,
  !    Random Variate Generation,
  !    in Handbook of Simulation,
  !    edited by Jerry Banks,
  !    Wiley, 1998, pages 168.
  !
  !    Reuven Rubinstein,
  !    Monte Carlo Optimization, Simulation, and Sensitivity
  !    of Queueing Networks,
  !    Krieger, 1992,
  !    ISBN: 0894647644,
  !    LC: QA298.R79.
  !
  !  Parameters:
  !
  !    Input, integer(ip) M, the dimension of the space.
  !
  !    Input, integer(ip) N, the number of points.
  !
  !    Input/output, integer(ip) SEED, a seed for the random
  !    number generator.
  !
  !    Output, real(dp) X(M,N), the points.
  !

    integer(ip), intent(in), value :: m
    integer(ip), intent(in), value :: n
    integer(ip), intent(inout) :: seed
    real(dp), intent(out) :: x(m,n)
    real(dp) :: exponent
    integer(ip) :: j
    real(dp) :: norm
    real(dp) :: r
    real(dp) :: r8_uniform_01

    exponent = 1.0_dp / real ( m, dp)

    do j = 1, n
  !
  !  Fill a vector with normally distributed values.
  !
      call r8vec_normal_01 ( m, seed, x(1:m,j) )
  !
  !  Compute the length of the vector.
  !
      norm = sqrt ( sum ( x(1:m,j)**2 ) )
  !
  !  Normalize the vector.
  !
      x(1:m,j) = x(1:m,j) / norm
  !
  !  Now compute a value to map the point ON the hypersphere INTO the hypersphere.
  !
      r = r8_uniform_01 ( seed )

      x(1:m,j) = r ** exponent * x(1:m,j)

    end do
  end subroutine hypersphere_01_interior_uniform

  subroutine hypersphere_01_surface_uniform ( m, n, seed, x ) &
        bind(C, name="hypersphere_01_surface_uniform")

  !*****************************************************************************80
  !
  !! HYPERSPHERE_01_SURFACE_UNIFORM: uniform points on unit hypersphere surface.
  !
  !  Discussion:
  !
  !    The hypersphere has center 0 and radius 1.
  !
  !    This procedure is valid for any spatial dimension.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    12 November 2010
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Russell Cheng,
  !    Random Variate Generation,
  !    in Handbook of Simulation,
  !    edited by Jerry Banks,
  !    Wiley, 1998, pages 168.
  !
  !    George Marsaglia,
  !    Choosing a point from the surface of a sphere,
  !    Annals of Mathematical Statistics,
  !    Volume 43, Number 2, April 1972, pages 645-646.
  !
  !    Reuven Rubinstein,
  !    Monte Carlo Optimization, Simulation, and Sensitivity
  !    of Queueing Networks,
  !    Krieger, 1992,
  !    ISBN: 0894647644,
  !    LC: QA298.R79.
  !
  !  Parameters:
  !
  !    Input, integer(ip) M, the dimension of the space.
  !
  !    Input, integer(ip) N, the number of points.
  !
  !    Input/output, integer(ip) SEED, a seed for the random
  !    number generator.
  !
  !    Output, real(dp) X(M,N), the points.
  !

    integer(ip), intent(in), value :: m
    integer(ip), intent(in), value :: n
    integer(ip), intent(inout) :: seed
    real(dp), intent(out) :: x(m,n)
    integer(ip) :: j
    real(dp) :: norm
  !
  !  Fill a matrix with normally distributed values.
  !
    call r8mat_normal_01 ( m, n, seed, x )
  !
  !  Normalize each column.
  !
    do j = 1, n
  !
  !  Compute the length of the vector.
  !
      norm = sqrt ( sum ( x(1:m,j)**2 ) )
  !
  !  Normalize the vector.
  !
      x(1:m,j) = x(1:m,j) / norm

    end do
  end subroutine hypersphere_01_surface_uniform

  pure function hypersphere_01_volume ( m ) &
        bind(C, name="hypersphere_01_volume")

  !*****************************************************************************80
  !
  !! HYPERSPHERE_01_VOLUME computes the volume of a unit hypersphere.
  !
  !  Discussion:
  !
  !    The unit hypersphere satisfies:
  !
  !      sum ( 1 <= I <= M ) X(I) * X(I) = 1
  !
  !    Results include:
  !
  !     M    Volume
  !
  !     1    2
  !     2    1        * PI
  !     3  ( 4 /   3) * PI
  !     4  ( 1 /   2) * PI^2
  !     5  ( 8 /  15) * PI^2
  !     6  ( 1 /   6) * PI^3
  !     7  (16 / 105) * PI^3
  !     8  ( 1 /  24) * PI^4
  !     9  (32 / 945) * PI^4
  !    10  ( 1 / 120) * PI^5
  !
  !    For the unit hypersphere, Volume(M) = 2 * PI * Volume(M-2)/ M
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 September 2003
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) M, the spatial dimension.
  !
  !    Output, real(dp) HYPERSPHERE_01_VOLUME, the volume.
  !

    real(dp) :: hypersphere_01_volume
    integer(ip), intent(in), value :: m
    integer(ip) :: i
    integer(ip) :: m2
    real(dp), parameter :: r8_pi = 3.141592653589793_dp
    real(dp) :: volume

    if ( mod ( m, 2 ) == 0 ) then
      m2 = m / 2
      volume = r8_pi ** m2
      do i = 1, m2
        volume = volume / real ( i, dp)
      end do
    else
      m2 = ( m - 1 ) / 2
      volume = ( r8_pi ** m2 ) * ( 2.0_dp ** m )
      do i = m2 + 1, 2 * m2 + 1
        volume = volume / real ( i, dp)
      end do
    end if

    hypersphere_01_volume = volume
  end function hypersphere_01_volume

  subroutine hypersphere_01_volume_values ( n_data, n, volume ) &
        bind(C, name="hypersphere_01_volume_values")

  !*****************************************************************************80
  !
  !! HYPERSPHERE_01_VOLUME_VALUES returns some volumes of the unit hypersphere.
  !
  !  Discussion:
  !
  !    The formula for the volume of the unit hypersphere is
  !
  !      Volume(N) = 2 * pi^(N/2) / ( N * Gamma ( N / 2 ) )
  !
  !    This function satisfies the relationships:
  !
  !      Volume(N) = 2 * pi * Volume(N-2) / N
  !      Volume(N) = Area(N) / N
  !
  !    Some values of the function include:
  !
  !       N  Volume
  !
  !       1    1
  !       2    1        * PI
  !       3  ( 4 /   3) * PI
  !       4  ( 1 /   2) * PI^2
  !       5  ( 8 /  15) * PI^2
  !       6  ( 1 /   6) * PI^3
  !       7  (16 / 105) * PI^3
  !       8  ( 1 /  24) * PI^4
  !       9  (32 / 945) * PI^4
  !      10  ( 1 / 120) * PI^5
  !
  !    In Mathematica, the function can be evaluated by:
  !
  !      2 * Pi^(n/2) / ( n * Gamma[n/2] )
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 August 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Stephen Wolfram,
  !    The Mathematica Book,
  !    Fourth Edition,
  !    Cambridge University Press, 1999,
  !    ISBN: 0-521-64314-7,
  !    LC: QA76.95.W65.
  !
  !  Parameters:
  !
  !    Input/output, integer(ip) N_DATA.  The user sets N_DATA to 0
  !    before the first call.  On each call, the routine increments N_DATA by 1,
  !    and returns the corresponding data; when there is no more data, the
  !    output value of N_DATA will be 0 again.
  !
  !    Output, integer(ip) N, the spatial dimension.
  !
  !    Output, real(dp) VOLUME, the volume.
  !

    integer(ip), parameter :: n_max = 20
    integer(ip), intent(inout) :: n_data
    integer(ip), intent(out) :: n
    real(dp), intent(out) :: volume
    integer(ip), save, dimension ( n_max ) :: n_vec = (/ &
       1,  2, &
       3,  4, &
       5,  6, &
       7,  8, &
       9, 10, &
      11, 12, &
      13, 14, &
      15, 16, &
      17, 18, &
      19, 20 /)
    real(dp), save, dimension ( n_max ) :: volume_vec = (/ &
      0.2000000000000000e+01_dp, &
      0.3141592653589793e+01_dp, &
      0.4188790204786391e+01_dp, &
      0.4934802200544679e+01_dp, &
      0.5263789013914325e+01_dp, &
      0.5167712780049970e+01_dp, &
      0.4724765970331401e+01_dp, &
      0.4058712126416768e+01_dp, &
      0.3298508902738707e+01_dp, &
      0.2550164039877345e+01_dp, &
      0.1884103879389900e+01_dp, &
      0.1335262768854589e+01_dp, &
      0.9106287547832831e+00_dp, &
      0.5992645293207921e+00_dp, &
      0.3814432808233045e+00_dp, &
      0.2353306303588932e+00_dp, &
      0.1409811069171390e+00_dp, &
      0.8214588661112823e-01_dp, &
      0.4662160103008855e-01_dp, &
      0.2580689139001406e-01_dp /)

    if ( n_data < 0 ) then
      n_data = 0
    end if

    n_data = n_data + 1

    if ( n_max < n_data ) then
      n_data = 0
      n = 0
      volume = 0.0_dp
    else
      n = n_vec(n_data)
      volume = volume_vec(n_data)
    end if
  end subroutine hypersphere_01_volume_values

  pure function hypersphere_area ( m, r ) &
        bind(C, name="hypersphere_area")

  !*****************************************************************************80
  !
  !! HYPERSPHERE_AREA computes the surface area of a hypersphere.
  !
  !  Discussion:
  !
  !    A hypersphere satisfies the equation:
  !
  !      sum ( ( P(1:M) - C(1:M) )^2 ) = R^2
  !
  !    M   Area
  !
  !    2      2       * PI   * R
  !    3      4       * PI   * R^2
  !    4      2       * PI^2 * R^3
  !    5      (8/3)   * PI^2 * R^4
  !    6                PI^3 * R^5
  !    7      (16/15) * PI^3 * R^6
  !
  !    Sphere_Area ( M, R ) = 2 * PI^(M/2) * R^(M-1) / Gamma ( M / 2 )
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 September 2003
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) M, the dimension of the space.
  !
  !    Input, real(dp) R, the radius.
  !
  !    Output, real(dp) HYPERSPHERE_AREA, the area.
  !

    real(dp) :: hypersphere_area
    integer(ip), intent(in), value :: m
    real(dp), intent(in), value :: r

    hypersphere_area = r ** ( m - 1  ) * hypersphere_01_area ( m )
  end function hypersphere_area

  pure subroutine hypersphere_stereograph ( m, n, x, x2 ) &
        bind(C, name="hypersphere_stereograph")

  !*****************************************************************************80
  !
  !! HYPERSPHERE_STEREOGRAPH: stereographic mapping of points on a hypersphere.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 May 2013
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) M, the spatial dimension.
  !    M must be at least 2.
  !
  !    Input, integer(ip) N, the number of points.
  !
  !    Input, real(dp) X(M,N), the points to be mapped.
  !
  !    Output, real(dp) X2(M-1,N), the stereographically mapped points.
  !

    integer(ip), intent(in), value :: m
    integer(ip), intent(in), value :: n
    real(dp), intent(in) :: x(m,n)
    real(dp), intent(out) :: x2(m-1,n)
    integer(ip) :: i

    x2(1:m-1,1:n) = x(1:m-1,1:n)

    do i = 1, m - 1
      x2(i,1:n) = x2(i,1:n) / ( 1.0_dp - x(m,1:n) )
    end do
  end subroutine hypersphere_stereograph

  pure subroutine hypersphere_stereograph_inverse ( m, n, x2, x ) &
        bind(C, name="hypersphere_stereograph_inverse")

  !*****************************************************************************80
  !
  !! HYPERSPHERE_STEREOGRAPH_INVERSE inverts a stereographic map.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 May 2013
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) M, the spatial dimension.
  !    M must be at least 2.
  !
  !    Input, integer(ip) N, the number of points.
  !
  !    Input, real(dp) X2(M-1,N), points in the plane.
  !
  !    Output, real(dp) X(M,N), points mapped back to the hypersphere.
  !

    integer(ip), intent(in), value :: m
    integer(ip), intent(in), value :: n
    real(dp), intent(in) :: x2(m-1,n)
    real(dp), intent(out) :: x(m,n)
    real(dp) :: d(n)
    integer(ip) :: i
    integer(ip) :: j

    x(1:m-1,1:n) = 2.0_dp * x2(1:m-1,1:n)

    do j = 1, n
      d(j) = sum ( x2(1:m-1,j) ** 2 )
    end do

    x(m,1:n) = d(1:n) - 1.0_dp

    do i = 1, m
      x(i,1:n) = x(i,1:n) / ( d(1:n) + 1.0_dp )
    end do
  end subroutine hypersphere_stereograph_inverse

  subroutine hypersphere_surface_uniform ( m, n, r, c, seed, x ) &
        bind(C, name="hypersphere_surface_uniform")

  !*****************************************************************************80
  !
  !! HYPERSPHERE_SURFACE_UNIFORM: uniform hypersphere surface samples
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    23 September 2018
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Russell Cheng,
  !    Random Variate Generation,
  !    in Handbook of Simulation,
  !    edited by Jerry Banks,
  !    Wiley, 1998, pages 168.
  !
  !    George Marsaglia,
  !    Choosing a point from the surface of a sphere,
  !    Annals of Mathematical Statistics,
  !    Volume 43, Number 2, April 1972, pages 645-646.
  !
  !    Reuven Rubinstein,
  !    Monte Carlo Optimization, Simulation, and Sensitivity
  !    of Queueing Networks,
  !    Wiley, 1986, page 234.
  !
  !  Parameters:
  !
  !    Input, integer(ip) M, the dimension of the space.
  !
  !    Input, integer(ip) N, the number of points.
  !
  !    Input, real(dp) R, the radius.
  !
  !    Input, real(dp) C(M), the center.
  !
  !    Input/output, integer(ip) SEED, a seed for the random number
  !    generator.
  !
  !    Output, real(dp) X(M,N), the points.
  !

    integer(ip), intent(in), value :: m
    integer(ip), intent(in), value :: n
    real(dp), intent(in), value :: r
    real(dp), intent(in) :: c(m)
    integer(ip), intent(inout) :: seed
    real(dp), intent(out) :: x(m,n)
    integer(ip) :: i

    call hypersphere_01_surface_uniform ( m, n, seed, x )
  !
  !  Scale by the radius.
  !
    x(1:m,1:n) = r * x(1:m,1:n)
  !
  !  Shift to the center.
  !
    do i = 1, m
      x(i,1:n) = x(i,1:n) + c(i)
    end do
  end subroutine hypersphere_surface_uniform

  pure subroutine hypersphere_to_cartesian ( m, n, c, r, theta, x ) &
        bind(C, name="hypersphere_to_cartesian")

  !*****************************************************************************80
  !
  !! HYPERSPHERE_TO_CARTESIAN: hypersphere to Cartesian coordinate transform.
  !
  !  Discussion:
  !
  !    We allow the trivial case M = 1; in that case alone, the value R
  !    must be assumed to have a sign.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 May 2013
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer M, the spatial dimension.
  !    1 <= M.
  !
  !    Input, integer N, the number of points to transform.
  !
  !    Input, real C(M), the center of the hypersphere.
  !
  !    Input, real R(N), the radius of the points on the hypersphere.
  !    Except for the trivial case M = 1, R is assumed nonnegative.
  !
  !    Input, real THETA(M-1,N), the coordinate angles of the points,
  !    measured in radians.
  !
  !    Output, real X(M,N), the Cartesian coordinates of the points.
  !

    integer(ip), intent(in), value :: m
    integer(ip), intent(in), value :: n
    real(dp), intent(in) :: c(m)
    real(dp), intent(in) :: r(n)
    real(dp), intent(in) :: theta(m-1,n)
    real(dp), intent(out) :: x(m,n)
    integer(ip) :: i
    integer(ip) :: i1
    integer(ip) :: i2

    if ( m == 1 ) then

      x(1,1:n) = r(1:n)

    else

      do i = 1, m
        x(i,1:n) = r(1:n)
      end do

      do i1 = 1, m - 1
        x(i1,1:n) = x(i1,1:n) * cos ( theta(i1,1:n) )
        do i2 = i1 + 1, m
          x(i2,1:n) = x(i2,1:n) * sin ( theta(i1,1:n) )
        end do
      end do
    end if
  !
  !  Add the center.
  !
    do i = 1, m
      x(i,1:n) = x(i,1:n) + c(i)
    end do
  end subroutine hypersphere_to_cartesian

  pure function hypersphere_volume ( m, r ) &
        bind(C, name="hypersphere_volume")

  !*****************************************************************************80
  !
  !! HYPERSPHERE_VOLUME computes the volume of a hypersphere.
  !
  !  Discussion:
  !
  !    A hypersphere satisfies the equation:
  !
  !      sum ( ( X(1:N) - PC(1:N) )^2 ) = R^2
  !
  !    where R is the radius and PC is the center.
  !
  !    Results include:
  !
  !    M     Volume
  !    -     -----------------------
  !    2                PI   * R^2
  !    3     (4/3)    * PI   * R^3
  !    4     (1/2)    * PI^2 * R^4
  !    5     (8/15)   * PI^2 * R^5
  !    6     (1/6)    * PI^3 * R^6
  !    7     (16/105) * PI^3 * R^7
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 September 2003
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) M, the dimension of the space.
  !
  !    Input, real(dp) R, the radius.
  !
  !    Output, real(dp) HYPERSPHERE_VOLUME, the volume.
  !

    real(dp) :: hypersphere_volume
    integer(ip), intent(in), value :: m
    real(dp), intent(in), value :: r

    hypersphere_volume = ( r ** m ) * hypersphere_01_volume ( m )
  end function hypersphere_volume

  function r8_uniform_01 ( seed ) &
        bind(C, name="r8_uniform_01")

  !*****************************************************************************80
  !
  !! R8_UNIFORM_01 returns a unit pseudorandom R8.
  !
  !  Discussion:
  !
  !    An R8 is a real(dp) value.
  !
  !    For now, the input quantity SEED is an integer variable.
  !
  !    This routine implements the recursion
  !
  !      seed = 16807 * seed mod ( 2^31 - 1 )
  !      r8_uniform_01 = seed / ( 2^31 - 1 )
  !
  !    The integer arithmetic never requires more than 32 bits,
  !    including a sign bit.
  !
  !    If the initial seed is 12345, then the first three computations are
  !
  !      Input     Output      R8_UNIFORM_01
  !      SEED      SEED
  !
  !         12345   207482415  0.096616
  !     207482415  1790989824  0.833995
  !    1790989824  2035175616  0.947702
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 July 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Paul Bratley, Bennett Fox, Linus Schrage,
  !    A Guide to Simulation,
  !    Springer Verlag, pages 201-202, 1983.
  !
  !    Pierre L'Ecuyer,
  !    Random Number Generation,
  !    in Handbook of Simulation,
  !    edited by Jerry Banks,
  !    Wiley Interscience, page 95, 1998.
  !
  !    Bennett Fox,
  !    Algorithm 647:
  !    Implementation and Relative Efficiency of Quasirandom
  !    Sequence Generators,
  !    ACM Transactions on Mathematical Software,
  !    Volume 12, Number 4, pages 362-376, 1986.
  !
  !    Peter Lewis, Allen Goodman, James Miller
  !    A Pseudo-Random Number Generator for the System/360,
  !    IBM Systems Journal,
  !    Volume 8, pages 136-143, 1969.
  !
  !  Parameters:
  !
  !    Input/output, integer(ip) SEED, the "seed" value, which should
  !    NOT be 0. On output, SEED has been updated.
  !
  !    Output, real(dp) R8_UNIFORM_01, a new pseudorandom variate,
  !    strictly between 0 and 1.
  !

    integer(ip), parameter :: i4_huge = 2147483647
    integer(ip) :: k
    real(dp) :: r8_uniform_01
    integer(ip), intent(inout) :: seed

    if ( seed == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
      write ( *, '(a)' ) '  Input value of SEED = 0.'
      stop 1
    end if

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r8_uniform_01 = real ( seed, dp) * 4.656612875e-10_dp
  end function r8_uniform_01

  pure function r8mat_norm_fro_affine ( m, n, a1, a2 ) &
        bind(C, name="r8mat_norm_fro_affine")

  !*****************************************************************************80
  !
  !! R8MAT_NORM_FRO_AFFINE returns the Frobenius norm of an R8MAT difference.
  !
  !  Discussion:
  !
  !    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
  !
  !    The Frobenius norm is defined as
  !
  !      R8MAT_NORM_FRO = sqrt (
  !        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J) * A(I,J) )
  !
  !    The matrix Frobenius norm is not derived from a vector norm, but
  !    is compatible with the vector L2 norm, so that:
  !
  !      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    24 March 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) M, the number of rows.
  !
  !    Input, integer(ip) N, the number of columns.
  !
  !    Input, real(dp) A1(M,N), A2(M,N), the matrices for whose
  !    difference the Frobenius norm is desired.
  !
  !    Output, real(dp) R8MAT_NORM_FRO_AFFINE, the Frobenius
  !    norm of A1 - A2.
  !

    integer(ip), intent(in), value :: m
    integer(ip), intent(in), value :: n
    real(dp), intent(in) :: a1(m,n)
    real(dp), intent(in) :: a2(m,n)
    real(dp) :: r8mat_norm_fro_affine
    real(dp) :: t
    integer(ip) :: i
    integer(ip) :: j

    t = 0.0_dp
    do j = 1, n
      do i = 1, m
        t = t + ( a1(i,j) - a2(i,j) )**2
      end do
    end do
    t = sqrt ( t )

    r8mat_norm_fro_affine = t
  end function r8mat_norm_fro_affine

  subroutine r8mat_normal_01 ( m, n, seed, r ) &
        bind(C, name="r8mat_normal_01")

  !*****************************************************************************80
  !
  !! R8MAT_NORMAL_01 returns a unit pseudonormal R8MAT.
  !
  !  Discussion:
  !
  !    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    12 November 2010
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Paul Bratley, Bennett Fox, Linus Schrage,
  !    A Guide to Simulation,
  !    Second Edition,
  !    Springer, 1987,
  !    ISBN: 0387964673,
  !    LC: QA76.9.C65.B73.
  !
  !    Bennett Fox,
  !    Algorithm 647:
  !    Implementation and Relative Efficiency of Quasirandom
  !    Sequence Generators,
  !    ACM Transactions on Mathematical Software,
  !    Volume 12, Number 4, December 1986, pages 362-376.
  !
  !    Pierre L'Ecuyer,
  !    Random Number Generation,
  !    in Handbook of Simulation,
  !    edited by Jerry Banks,
  !    Wiley, 1998,
  !    ISBN: 0471134031,
  !    LC: T57.62.H37.
  !
  !    Peter Lewis, Allen Goodman, James Miller,
  !    A Pseudo-Random Number Generator for the System/360,
  !    IBM Systems Journal,
  !    Volume 8, 1969, pages 136-143.
  !
  !  Parameters:
  !
  !    Input, integer(ip) M, N, the number of rows and columns
  !    in the array.
  !
  !    Input/output, integer(ip) SEED, the "seed" value, which
  !    should NOT be 0.  On output, SEED has been updated.
  !
  !    Output, real(dp) R(M,N), the array of pseudonormal values.
  !

    integer(ip), intent(in), value :: m
    integer(ip), intent(in), value :: n
    integer(ip), intent(inout) :: seed
    real(dp), intent(out) :: r(m,n)

    call r8vec_normal_01 ( m * n, seed, r )
  end subroutine r8mat_normal_01

  pure subroutine r8mat_uniform_01 ( m, n, seed, r ) &
        bind(C, name="r8mat_uniform_01")

  !*****************************************************************************80
  !
  !! R8MAT_UNIFORM_01 fills an R8MAT with unit pseudorandom numbers.
  !
  !  Discussion:
  !
  !    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    11 August 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Paul Bratley, Bennett Fox, Linus Schrage,
  !    A Guide to Simulation,
  !    Springer Verlag, pages 201-202, 1983.
  !
  !    Bennett Fox,
  !    Algorithm 647:
  !    Implementation and Relative Efficiency of Quasirandom
  !    Sequence Generators,
  !    ACM Transactions on Mathematical Software,
  !    Volume 12, Number 4, pages 362-376, 1986.
  !
  !    Peter Lewis, Allen Goodman, James Miller,
  !    A Pseudo-Random Number Generator for the System/360,
  !    IBM Systems Journal,
  !    Volume 8, pages 136-143, 1969.
  !
  !  Parameters:
  !
  !    Input, integer(ip) M, N, the number of rows and columns in
  !    the array.
  !
  !    Input/output, integer(ip) SEED, the "seed" value, which
  !    should NOT be 0.  On output, SEED has been updated.
  !
  !    Output, real(dp) R(M,N), the array of pseudorandom values.
  !

    integer(ip), intent(in), value :: m
    integer(ip), intent(in), value :: n
    integer(ip), intent(inout) :: seed
    real(dp), intent(out) :: r(m,n)
    integer(ip) :: i
    integer(ip), parameter :: i4_huge = 2147483647
    integer(ip) :: j
    integer(ip) :: k

    do j = 1, n
      do i = 1, m

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed < 0 ) then
          seed = seed + i4_huge
        end if

        r(i,j) = real ( seed, dp) * 4.656612875e-10_dp

      end do
    end do
  end subroutine r8mat_uniform_01

  subroutine r8vec_normal_01 ( n, seed, x ) &
        bind(C, name="r8vec_normal_01")

  !*****************************************************************************80
  !
  !! R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
  !
  !  Discussion:
  !
  !    An R8VEC is a vector of R8's.
  !
  !    The standard normal probability distribution function (PDF) has
  !    mean 0 and standard deviation 1.
  !
  !    This routine can generate a vector of values on one call.  It
  !    has the feature that it should provide the same results
  !    in the same order no matter how we break up the task.
  !
  !    Before calling this routine, the user may call RANDOM_SEED
  !    in order to set the seed of the random number generator.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    06 August 2013
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of values desired.
  !
  !    Input/output, integer(ip) SEED, a seed for the random
  !    number generator.
  !
  !    Output, real(dp) X(N), a sample of the standard normal PDF.
  !
  !  Local parameters:
  !
  !    Local, real(dp) R(N+1), is used to store some uniform
  !    random values.  Its dimension is N+1, but really it is only needed
  !    to be the smallest even number greater than or equal to N.
  !
  !    Local, integer X_LO_INDEX, X_HI_INDEX, records the range of entries of
  !    X that we need to compute.
  !

    integer(ip), intent(in), value :: n
    integer(ip), intent(inout) :: seed
    real(dp), intent(out) :: x(n)
    integer(ip) :: m
    real(dp) :: r(n+1)
    real(dp), parameter :: r8_pi = 3.141592653589793_dp
    real(dp) :: r8_uniform_01
    integer(ip) :: x_hi_index
    integer(ip) :: x_lo_index
  !
  !  Record the range of X we need to fill in.
  !
    x_lo_index = 1
    x_hi_index = n
  !
  !  Maybe we don't need any more values.
  !
    if ( x_hi_index - x_lo_index + 1 == 1 ) then

      r(1) = r8_uniform_01 ( seed )

      if ( r(1) == 0.0_dp ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_NORMAL_01 - Fatal error!'
        write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
        stop 1
      end if

      r(2) = r8_uniform_01 ( seed )

      x(x_hi_index) = &
               sqrt ( -2.0_dp * log ( r(1) ) ) * cos ( 2.0_dp * r8_pi * r(2) )
  !
  !  If we require an even number of values, that's easy.
  !
    else if ( mod ( x_hi_index - x_lo_index + 1, 2 ) == 0 ) then

      m = ( x_hi_index - x_lo_index + 1 ) / 2

      call r8vec_uniform_01 ( 2*m, seed, r )

      x(x_lo_index:x_hi_index-1:2) = &
        sqrt ( -2.0_dp * log ( r(1:2*m-1:2) ) ) &
        * cos ( 2.0_dp * r8_pi * r(2:2*m:2) )

      x(x_lo_index+1:x_hi_index:2) = &
        sqrt ( -2.0_dp * log ( r(1:2*m-1:2) ) ) &
        * sin ( 2.0_dp * r8_pi * r(2:2*m:2) )
  !
  !  If we require an odd number of values, we generate an even number,
  !  and handle the last pair specially, storing one in X(N), and
  !  saving the other for later.
  !
    else

      x_hi_index = x_hi_index - 1

      m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

      call r8vec_uniform_01 ( 2*m, seed, r )

      x(x_lo_index:x_hi_index-1:2) = &
        sqrt ( -2.0_dp * log ( r(1:2*m-3:2) ) ) &
        * cos ( 2.0_dp * r8_pi * r(2:2*m-2:2) )

      x(x_lo_index+1:x_hi_index:2) = &
        sqrt ( -2.0_dp * log ( r(1:2*m-3:2) ) ) &
        * sin ( 2.0_dp * r8_pi * r(2:2*m-2:2) )

      x(n) = sqrt ( -2.0_dp * log ( r(2*m-1) ) ) &
        * cos ( 2.0_dp * r8_pi * r(2*m) )

    end if
  end subroutine r8vec_normal_01

  subroutine r8vec_uniform_01 ( n, seed, r ) &
        bind(C, name="r8vec_uniform_01")

  !*****************************************************************************80
  !
  !! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
  !
  !  Discussion:
  !
  !    An R8VEC is a vector of R8's.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 July 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Paul Bratley, Bennett Fox, Linus Schrage,
  !    A Guide to Simulation,
  !    Springer Verlag, pages 201-202, 1983.
  !
  !    Bennett Fox,
  !    Algorithm 647:
  !    Implementation and Relative Efficiency of Quasirandom
  !    Sequence Generators,
  !    ACM Transactions on Mathematical Software,
  !    Volume 12, Number 4, pages 362-376, 1986.
  !
  !    Peter Lewis, Allen Goodman, James Miller
  !    A Pseudo-Random Number Generator for the System/360,
  !    IBM Systems Journal,
  !    Volume 8, pages 136-143, 1969.
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of entries in the vector.
  !
  !    Input/output, integer(ip) SEED, the "seed" value, which
  !    should NOT be 0.  On output, SEED has been updated.
  !
  !    Output, real(dp) R(N), the vector of pseudorandom values.
  !

    integer(ip), intent(in), value :: n
    integer(ip), intent(inout) :: seed
    real(dp), intent(out) :: r(n)
    integer(ip) :: i
    integer(ip) :: k

    if ( seed == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
      write ( *, '(a)' ) '  Input value of SEED = 0.'
      stop 1
    end if

    do i = 1, n

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + 2147483647
      end if

      r(i) = real ( seed, dp) * 4.656612875e-10_dp

    end do
  end subroutine r8vec_uniform_01

  pure subroutine sphere_stereograph ( m, n, p, q ) &
        bind(C, name="sphere_stereograph")

  !*****************************************************************************80
  !
  !! SPHERE_STEREOGRAPH computes the stereographic image of points on a sphere.
  !
  !  Discussion:
  !
  !    We start with a sphere of radius 1 and center (0,0,0).
  !
  !    The north pole N = (0,0,1) is the point of tangency to the sphere
  !    of a plane, and the south pole S = (0,0,-1) is the focus for the
  !    stereographic projection.
  !
  !    For any point P on the sphere, the stereographic projection Q of the
  !    point is defined by drawing the line from S through P, and computing
  !    Q as the intersection of this line with the plane.
  !
  !    Actually, we allow the spatial dimension M to be arbitrary.  Values
  !    of M make sense starting with 2.  The north and south poles are
  !    selected as the points (0,0,...,+1) and (0,0,...,-1).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    26 September 2012
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    C F Marcus,
  !    The stereographic projection in vector notation,
  !    Mathematics Magazine,
  !    Volume 39, Number 2, March 1966, pages 100-102.
  !
  !  Parameters:
  !
  !    Input, integer(ip) M, the spatial dimension.
  !
  !    Input, integer(ip) N, the number of points.
  !
  !    Input, real(dp) P(M,N), a set of points on the unit sphere.
  !
  !    Output, real(dp) Q(M,N), the coordinates of the
  !    image points.
  !

    integer(ip), intent(in), value :: m
    integer(ip), intent(in), value :: n
    real(dp), intent(in) :: p(m,n)
    real(dp), intent(out) :: q(m,n)
    integer(ip) :: i
    integer(ip) :: j

    do j = 1, n
      do i = 1, m - 1
        q(i,j) = 2.0_dp * p(i,j) / ( 1.0_dp + p(m,j) )
      end do
      q(m,j) = 1.0_dp
    end do
  end subroutine sphere_stereograph

  pure subroutine sphere_stereograph_inverse ( m, n, q, p ) &
        bind(C, name="sphere_stereograph_inverse")

  !*****************************************************************************80
  !
  !! SPHERE_STEREOGRAPH_INVERSE computes stereographic preimages of points.
  !
  !  Discussion:
  !
  !    We start with a sphere of radius 1 and center (0,0,0).
  !
  !    The north pole N = (0,0,1) is the point of tangency to the sphere
  !    of a plane, and the south pole S = (0,0,-1) is the focus for the
  !    stereographic projection.
  !
  !    For any point Q on the plane, the stereographic inverse projection
  !    P of the point is defined by drawing the line from S through Q, and
  !    computing P as the intersection of this line with the sphere.
  !
  !    Actually, we allow the spatial dimension M to be arbitrary.  Values
  !    of M make sense starting with 2.  The north and south poles are
  !    selected as the points (0,0,...,+1) and (0,0,...,-1).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    26 September 2012
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    C F Marcus,
  !    The stereographic projection in vector notation,
  !    Mathematics Magazine,
  !    Volume 39, Number 2, March 1966, pages 100-102.
  !
  !  Parameters:
  !
  !    Input, integer(ip) M, the spatial dimension.
  !
  !    Input, integer(ip) N, the number of points.
  !
  !    Input, real(dp) Q(M,N), the points, which are presumed to lie
  !    on the plane Z = 1.
  !
  !    Output, real(dp) P(M,N), the stereographic
  !    inverse projections of the points.
  !

    integer(ip), intent(in), value :: m
    integer(ip), intent(in), value :: n
    real(dp), intent(in) :: q(m,n)
    real(dp), intent(out) :: p(m,n)
    integer(ip) :: j
    real(dp) :: qn

    do j = 1, n

      qn = sum ( q(1:m-1,j)**2 )

      p(1:m-1,j) = 4.0_dp * q(1:m-1,j) / ( 4.0_dp + qn )

      p(m,j) = ( 4.0_dp - qn ) / ( 4.0_dp + qn )

    end do
  end subroutine sphere_stereograph_inverse

end module hypersphere_properties_mod
