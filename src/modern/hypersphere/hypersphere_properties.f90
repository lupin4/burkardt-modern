!> hypersphere_properties — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

subroutine cartesian_to_hypersphere ( m, n, c, x, r, theta )

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
!    Input, integer M, the spatial dimension.
!    1 <= M.
!
!    Input, integer N, the number of points to transform.
!
!    Input, double precision C(M), the center of the hypersphere.
!
!    Input, double precision X(M,N), the Cartesian coordinates of the points.
!
!    Output, double precision R(N), the radius of the points on the 
!    hypersphere.  Except for the trivial case M = 1, R is assumed nonnegative.
!
!    Output, double precision THETA(M-1,N), the coordinate angles of the 
!    points, measured in radians.
!
  implicit none

  integer m
  integer n

  double precision c(m)
  integer i
  integer i1
  integer j
  double precision r(n)
  double precision theta(m-1,n)
  double precision top(n)
  double precision x(m,n)
  double precision x2(m,n)
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
  theta(1:m-1,n) = 0.0D+00

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
    theta(m-1,j) = 2.0 * atan2 ( x2(m,j), top(j) )
  end do
end

function hypersphere_01_area ( m )

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
!    Input, integer M, the dimension of the space.
!
!    Output, double precision HYPERSPHERE_01_AREA, the area.
!
  implicit none

  double precision area
  double precision hypersphere_01_area
  integer i
  integer m
  integer m2
  double precision , parameter :: r8_pi = 3.141592653589793D+00

  if ( mod ( m, 2 ) == 0 ) then
    m2 = m / 2
    area = 2.0D+00 * ( r8_pi ) ** m2
    do i = 1, m2 - 1
      area = area / real ( i)
    end do
  else
    m2 = ( m - 1 ) / 2
    area = ( r8_pi ** m2 ) * ( 2.0D+00 ** m )
    do i = m2 + 1, 2 * m2
      area = area / real ( i)
    end do
  end if

  hypersphere_01_area = area
end

subroutine hypersphere_01_area_values ( n_data, n, area )

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
!    Input/output, integer N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer N, the spatial dimension.
!
!    Output, double precision AREA, the area
!    in that dimension.
!
  implicit none

  integer , parameter :: n_max = 20

  double precision area
  double precision , save, dimension ( n_max ) :: area_vec = (/ &
    0.2000000000000000D+01, &
    0.6283185307179586D+01, &
    0.1256637061435917D+02, &
    0.1973920880217872D+02, &
    0.2631894506957162D+02, &
    0.3100627668029982D+02, &
    0.3307336179231981D+02, &
    0.3246969701133415D+02, &
    0.2968658012464836D+02, &
    0.2550164039877345D+02, &
    0.2072514267328890D+02, &
    0.1602315322625507D+02, &
    0.1183817381218268D+02, &
    0.8389703410491089D+01, &
    0.5721649212349567D+01, &
    0.3765290085742291D+01, &
    0.2396678817591364D+01, &
    0.1478625959000308D+01, &
    0.8858104195716824D+00, &
    0.5161378278002812D+00 /)
  integer n_data
  integer n
  integer , save, dimension ( n_max ) :: n_vec = (/ &
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
    area = 0.0D+00
  else
    n = n_vec(n_data)
    area = area_vec(n_data)
  end if
end

subroutine hypersphere_01_interior_uniform ( m, n, seed, x )

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
!    Input, integer M, the dimension of the space.
!
!    Input, integer N, the number of points.
!
!    Input/output, integer SEED, a seed for the random
!    number generator.
!
!    Output, double precision X(M,N), the points.
!
  implicit none

  integer m
  integer n

  double precision exponent
  integer j
  double precision norm
  double precision r
  double precision r8_uniform_01
  integer seed
  double precision x(m,n)

  exponent = 1.0D+00 / real ( m)

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
end

subroutine hypersphere_01_surface_uniform ( m, n, seed, x )

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
!    Input, integer M, the dimension of the space.
!
!    Input, integer N, the number of points.
!
!    Input/output, integer SEED, a seed for the random
!    number generator.
!
!    Output, double precision X(M,N), the points.
!
  implicit none

  integer m
  integer n

  integer j
  double precision norm
  integer seed
  double precision x(m,n)
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
end

function hypersphere_01_volume ( m )

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
!    Input, integer M, the spatial dimension.
!
!    Output, double precision HYPERSPHERE_01_VOLUME, the volume.
!
  implicit none

  double precision hypersphere_01_volume
  integer i
  integer m
  integer m2
  double precision , parameter :: r8_pi = 3.141592653589793D+00
  double precision volume

  if ( mod ( m, 2 ) == 0 ) then
    m2 = m / 2
    volume = r8_pi ** m2
    do i = 1, m2
      volume = volume / real ( i)
    end do
  else
    m2 = ( m - 1 ) / 2
    volume = ( r8_pi ** m2 ) * ( 2.0D+00 ** m )
    do i = m2 + 1, 2 * m2 + 1
      volume = volume / real ( i)
    end do
  end if

  hypersphere_01_volume = volume
end

subroutine hypersphere_01_volume_values ( n_data, n, volume )

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
!    Input/output, integer N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer N, the spatial dimension.
!
!    Output, double precision VOLUME, the volume.
!
  implicit none

  integer , parameter :: n_max = 20

  integer n_data
  integer n
  integer , save, dimension ( n_max ) :: n_vec = (/ &
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
  double precision volume
  double precision , save, dimension ( n_max ) :: volume_vec = (/ &
    0.2000000000000000D+01, &
    0.3141592653589793D+01, &
    0.4188790204786391D+01, &
    0.4934802200544679D+01, &
    0.5263789013914325D+01, &
    0.5167712780049970D+01, &
    0.4724765970331401D+01, &
    0.4058712126416768D+01, &
    0.3298508902738707D+01, &
    0.2550164039877345D+01, &
    0.1884103879389900D+01, &
    0.1335262768854589D+01, &
    0.9106287547832831D+00, &
    0.5992645293207921D+00, &
    0.3814432808233045D+00, &
    0.2353306303588932D+00, &
    0.1409811069171390D+00, &
    0.8214588661112823D-01, &
    0.4662160103008855D-01, &
    0.2580689139001406D-01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    volume = 0.0D+00
  else
    n = n_vec(n_data)
    volume = volume_vec(n_data)
  end if
end

function hypersphere_area ( m, r )

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
!    Input, integer M, the dimension of the space.
!
!    Input, double precision R, the radius.
!
!    Output, double precision HYPERSPHERE_AREA, the area.
!
  implicit none

  double precision hypersphere_01_area
  double precision hypersphere_area
  integer m
  double precision r

  hypersphere_area = r ** ( m - 1  ) * hypersphere_01_area ( m )
end

subroutine hypersphere_stereograph ( m, n, x, x2 )

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
!    Input, integer M, the spatial dimension.
!    M must be at least 2.
!
!    Input, integer N, the number of points.
!
!    Input, double precision X(M,N), the points to be mapped.
!
!    Output, double precision X2(M-1,N), the stereographically mapped points.
!
  implicit none

  integer m
  integer n

  integer i
  double precision x(m,n)
  double precision x2(m-1,n)

  x2(1:m-1,1:n) = x(1:m-1,1:n)

  do i = 1, m - 1
    x2(i,1:n) = x2(i,1:n) / ( 1.0D+00 - x(m,1:n) )
  end do
end

subroutine hypersphere_stereograph_inverse ( m, n, x2, x )

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
!    Input, integer M, the spatial dimension.
!    M must be at least 2.
!
!    Input, integer N, the number of points.
!
!    Input, double precision X2(M-1,N), points in the plane.
!
!    Input, double precision X(M,N), points mapped back to the hypersphere.
!
  implicit none

  integer m
  integer n

  double precision d(n)
  integer i
  integer j
  double precision x(m,n)
  double precision x2(m-1,n)

  x(1:m-1,1:n) = 2.0D+00 * x2(1:m-1,1:n)

  do j = 1, n
    d(j) = sum ( x2(1:m-1,j) ** 2 )
  end do

  x(m,1:n) = d(1:n) - 1.0D+00
  
  do i = 1, m
    x(i,1:n) = x(i,1:n) / ( d(1:n) + 1.0D+00 )
  end do
end

subroutine hypersphere_surface_uniform ( m, n, r, c, seed, x )

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
!    Input, integer M, the dimension of the space.
!
!    Input, integer N, the number of points.
!
!    Input, double precision R, the radius.
!
!    Input, double precision C(M), the center.
!
!    Input/output, integer SEED, a seed for the random number 
!    generator.
!
!    Output, double precision X(M,N), the points.
!
  implicit none

  integer m
  integer n

  double precision c(m)
  integer i
  double precision r
  integer seed
  double precision x(m,n)

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
end

subroutine hypersphere_to_cartesian ( m, n, c, r, theta, x )

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
  implicit none

  integer m
  integer n

  double precision c(m)
  integer i
  integer i1
  integer i2
  double precision r(n)
  double precision theta(m-1,n)
  double precision x(m,n)

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
end

function hypersphere_volume ( m, r )

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
!    Input, integer M, the dimension of the space.
!
!    Input, double precision R, the radius.
!
!    Output, double precision HYPERSPHERE_VOLUME, the volume.
!
  implicit none

  double precision hypersphere_01_volume
  double precision hypersphere_volume
  integer m
  double precision r

  hypersphere_volume = ( r ** m ) * hypersphere_01_volume ( m )
end

subroutine sphere_stereograph ( m, n, p, q )

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
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, double precision P(M,N), a set of points on the unit sphere.
!
!    Output, double precision Q(M,N), the coordinates of the
!    image points.
!
  implicit none

  integer m
  integer n

  integer i
  integer j
  double precision p(m,n)
  double precision q(m,n)

  do j = 1, n
    do i = 1, m - 1
      q(i,j) = 2.0D+00 * p(i,j) / ( 1.0D+00 + p(m,j) )
    end do
    q(m,j) = 1.0D+00
  end do
end

subroutine sphere_stereograph_inverse ( m, n, q, p )

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
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, double precision Q(M,N), the points, which are presumed to lie
!    on the plane Z = 1.
!
!    Output, double precision P(M,N), the stereographic
!    inverse projections of the points.
!
  implicit none

  integer m
  integer n

  integer j
  double precision p(m,n)
  double precision q(m,n)
  double precision qn

  do j = 1, n

    qn = sum ( q(1:m-1,j)**2 )

    p(1:m-1,j) = 4.0D+00 * q(1:m-1,j) / ( 4.0D+00 + qn )

    p(m,j) = ( 4.0D+00 - qn ) / ( 4.0D+00 + qn )

  end do
end
