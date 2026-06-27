!> sphere_quad — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

function arc_cosine ( c )

!*****************************************************************************80
!
!! ARC_COSINE computes the arc cosine function, with argument truncation.
!
!  Discussion:
!
!    If you call your system ACOS routine with an input argument that is
!    even slightly outside the range [-1.0, 1.0 ], you may get an unpleasant 
!    surprise (I did).
!
!    This routine simply truncates arguments outside the range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision C, the argument.
!
!    Output, double precision ARC_COSINE, an angle whose cosine is C.
!
  implicit none

  double precision arc_cosine
  double precision c
  double precision c2

  c2 = c
  c2 = max ( c2, - 1.0D+00 )
  c2 = min ( c2, + 1.0D+00 )

  arc_cosine = acos ( c2 )
end

function arc_sine ( s )

!*****************************************************************************80
!
!! ARC_SINE computes the arc sine function, with argument truncation.
!
!  Discussion:
!
!    If you call your system ASIN routine with an input argument that is
!    even slightly outside the range [-1.0, 1.0 ], you may get an unpleasant 
!    surprise (I did).
!
!    This routine simply truncates arguments outside the range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision S, the argument.
!
!    Output, double precision ARC_SINE, an angle whose sine is S.
!
  implicit none

  double precision arc_sine
  double precision s
  double precision s2

  s2 = s
  s2 = max ( s2, - 1.0D+00 )
  s2 = min ( s2, + 1.0D+00 )

  arc_sine = asin ( s2 )
end

function atan4 ( y, x )

!*****************************************************************************80
!
!! ATAN4 computes the inverse tangent of the ratio Y / X.
!
!  Discussion:
!
!    ATAN4 returns an angle whose tangent is ( Y / X ), a job which
!    the built in functions ATAN and ATAN2 already do.
!
!    However:
!
!    * ATAN4 always returns a positive angle, between 0 and 2 PI,
!      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
!      and [-PI,+PI] respectively;
!
!    * ATAN4 accounts for the signs of X and Y, (as does ATAN2).  The ATAN
!     function by contrast always returns an angle in the first or fourth
!     quadrants.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision Y, X, two quantities which represent the 
!    tangent of an angle.  If Y is not zero, then the tangent is (Y/X).
!
!    Output, double precision ATAN4, an angle between 0 and 2 * PI, 
!    whose tangent is (Y/X), and which lies in the appropriate quadrant so 
!    that the signs of its cosine and sine match those of X and Y.
!
  implicit none

  double precision abs_x
  double precision abs_y
  double precision atan4
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision theta
  double precision theta_0
  double precision x
  double precision y
!
!  Special cases:
!
  if ( x == 0.0D+00 ) then

    if ( 0.0D+00 < y ) then
      theta = pi / 2.0D+00
    else if ( y < 0.0D+00 ) then
      theta = 3.0D+00 * pi / 2.0D+00
    else if ( y == 0.0D+00 ) then
      theta = 0.0D+00
    end if

  else if ( y == 0.0D+00 ) then

    if ( 0.0D+00 < x ) then
      theta = 0.0D+00
    else if ( x < 0.0D+00 ) then
      theta = PI
    end if
!
!  We assume that ATAN2 is correct when both arguments are positive.
!
  else

    abs_y = abs ( y )
    abs_x = abs ( x )

    theta_0 = atan2 ( abs_y, abs_x )

    if ( 0.0D+00 < x .and. 0.0D+00 < y ) then
      theta = theta_0
    else if ( x < 0.0D+00 .and. 0.0D+00 < y ) then
      theta = pi - theta_0
    else if ( x < 0.0D+00 .and. y < 0.0D+00 ) then
      theta = pi + theta_0
    else if ( 0.0D+00 < x .and. y < 0.0D+00 ) then
      theta = 2.0D+00 * pi - theta_0
    end if

  end if

  atan4 = theta
end

subroutine icos_shape ( point_num, edge_num, face_num, face_order_max, &
  point_coord, edge_point, face_order, face_point )

!*****************************************************************************80
!
!! ICOS_SHAPE describes an icosahedron.
!
!  Discussion:
!
!    The input data required for this routine can be retrieved from ICOS_SIZE.
!
!    The vertices lie on the unit sphere.
!
!    The dual of an icosahedron is a dodecahedron.
!
!    The data has been rearranged from a previous assignment.  
!    The STRIPACK program refuses to triangulate data if the first
!    three nodes are "collinear" on the sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer POINT_NUM, the number of points (12).
!
!    Input, integer EDGE_NUM, the number of edges (30).
!
!    Input, integer FACE_NUM, the number of faces (20).
!
!    Input, integer FACE_ORDER_MAX, the maximum number of 
!    vertices per face (3).
!
!    Output, double precision POINT_COORD(3,POINT_NUM), the points.
!
!    Output, integer EDGE_POINT(2,EDGE_NUM), the points that 
!    make up each edge, listed in ascending order of their indexes.
!
!    Output, integer FACE_ORDER(FACE_NUM), the number of vertices
!    per face.
!
!    Output, integer FACE_POINT(FACE_ORDER_MAX,FACE_NUM); 
!    FACE_POINT(I,J) is the index of the I-th point in the J-th face.  The
!    points are listed in the counter clockwise direction defined
!    by the outward normal at the face.  The nodes of each face are ordered 
!    so that the lowest index occurs first.  The faces are then sorted by
!    nodes.
!
  implicit none

  integer edge_num
  integer , parameter :: edge_order = 2
  integer face_num
  integer face_order_max
  integer point_num

  double precision a
  double precision b
  integer edge_point(edge_order,edge_num)
  integer face_order(face_num)
  integer face_point(face_order_max,face_num)
  double precision phi
  double precision point_coord(3,point_num)
  double precision z
!
!  Set the point coordinates.
!
  phi = 0.5D+00 * ( sqrt ( 5.0D+00 ) + 1.0D+00 )

  a = phi / sqrt ( 1.0D+00 + phi * phi )
  b = 1.0D+00 / sqrt ( 1.0D+00 + phi * phi )
  z = 0.0D+00
!
!  A*A + B*B + Z*Z = 1.
!
  point_coord(1:3,1:point_num) = reshape ( (/ &
      a,  b,  z, &
      a, -b,  z, &
      b,  z,  a, &
      b,  z, -a, &
      z,  a,  b, &
      z,  a, -b, &
      z, -a,  b, &
      z, -a, -b, &
     -b,  z,  a, &
     -b,  z, -a, &
     -a,  b,  z, &
     -a, -b,  z /), (/ 3, point_num /) )
!
!  Set the edges.
!
  edge_point(1:edge_order,1:edge_num) = reshape ( (/ &
     1,  2, &
     1,  3, &
     1,  4, &
     1,  5, &
     1,  6, &
     2,  3, &
     2,  4, &
     2,  7, &
     2,  8, &
     3,  5, &
     3,  7, &
     3,  9, &
     4,  6, &
     4,  8, &
     4, 10, &
     5,  6, &
     5,  9, &
     5, 11, &
     6, 10, &
     6, 11, &
     7,  8, &
     7,  9, &
     7, 12, &
     8, 10, &
     8, 12, &
     9, 11, &
     9, 12, &
    10, 11, &
    10, 12, &
    11, 12 /), (/ edge_order, edge_num /) )
!
!  Set the face orders.
!
  face_order(1:face_num) = (/ &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3 /)
!
!  Set the faces.
!
  face_point(1:face_order_max,1:face_num) = reshape ( (/ &
     1,  2,  4, &
     1,  3,  2, &
     1,  4,  6, &
     1,  5,  3, &
     1,  6,  5, &
     2,  3,  7, &
     2,  7,  8, &
     2,  8,  4, &
     3,  5,  9, &
     3,  9,  7, &
     4,  8, 10, &
     4, 10,  6, &
     5,  6, 11, &
     5, 11,  9, &
     6, 10, 11, &
     7,  9, 12, &
     7, 12,  8, &
     8, 12, 10, &
     9, 11, 12, &
    10, 12, 11 /), (/ face_order_max, face_num /) )
end

subroutine icos_size ( point_num, edge_num, face_num, face_order_max )

!*****************************************************************************80
!
!! ICOS_SIZE gives "sizes" for an icosahedron in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer POINT_NUM, the number of points.
!
!    Output, integer EDGE_NUM, the number of edges.
!
!    Output, integer FACE_NUM, the number of faces.
!
!    Output, integer FACE_ORDER_MAX, the maximum order of any face.
!
  implicit none

  integer edge_num
  integer face_num
  integer face_order_max
  integer point_num

  point_num = 12
  edge_num = 30
  face_num = 20
  face_order_max = 3
end

function r8_gamma ( x )

!*****************************************************************************80
!
!! R8_GAMMA evaluates Gamma(X) for a real argument.
!
!  Discussion:
!
!    This routine calculates the gamma function for a real argument X.
!
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the gamma
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for 12 <= X are from reference 2.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    An Overview of Software Development for Special Functions,
!    in Numerical Analysis Dundee, 1975,
!    edited by GA Watson,
!    Lecture Notes in Mathematics 506,
!    Springer, 1976.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, double precision X, the argument of the function.
!
!    Output, double precision R8_GAMMA, the value of the function.
!
  implicit none

  double precision , dimension ( 7 ) :: c = (/ &
   -1.910444077728D-03, &
    8.4171387781295D-04, &
   -5.952379913043012D-04, &
    7.93650793500350248D-04, &
   -2.777777777777681622553D-03, &
    8.333333333333333331554247D-02, &
    5.7083835261D-03 /)
  double precision , parameter :: eps = 2.22D-16
  double precision fact
  integer i
  integer n
  double precision , dimension ( 8 ) :: p = (/ &
    -1.71618513886549492533811D+00, &
     2.47656508055759199108314D+01, &
    -3.79804256470945635097577D+02, &
     6.29331155312818442661052D+02, &
     8.66966202790413211295064D+02, &
    -3.14512729688483675254357D+04, &
    -3.61444134186911729807069D+04, &
     6.64561438202405440627855D+04 /)
  logical parity
  double precision , parameter :: pi = 3.1415926535897932384626434D+00
  double precision , dimension ( 8 ) :: q = (/ &
    -3.08402300119738975254353D+01, &
     3.15350626979604161529144D+02, &
    -1.01515636749021914166146D+03, &
    -3.10777167157231109440444D+03, &
     2.25381184209801510330112D+04, &
     4.75584627752788110767815D+03, &
    -1.34659959864969306392456D+05, &
    -1.15132259675553483497211D+05 /)
  double precision r8_gamma
  double precision res
  double precision , parameter :: sqrtpi = 0.9189385332046727417803297D+00
  double precision sum
  double precision x
  double precision , parameter :: xbig = 171.624D+00
  double precision xden
  double precision , parameter :: xinf = 1.0D+30
  double precision , parameter :: xminin = 2.23D-308
  double precision xnum
  double precision y
  double precision y1
  double precision ysq
  double precision z

  parity = .false.
  fact = 1.0D+00
  n = 0
  y = x
!
!  Argument is negative.
!
  if ( y <= 0.0D+00 ) then

    y = - x
    y1 = aint ( y )
    res = y - y1

    if ( res /= 0.0D+00 ) then

      if ( y1 /= aint ( y1 * 0.5D+00 ) * 2.0D+00 ) then
        parity = .true.
      end if

      fact = - pi / sin ( pi * res )
      y = y + 1.0D+00

    else

      res = xinf
      r8_gamma = res
    end if

  end if
!
!  Argument is positive.
!
  if ( y < eps ) then
!
!  Argument < EPS.
!
    if ( xminin <= y ) then
      res = 1.0D+00 / y
    else
      res = xinf
      r8_gamma = res
    end if

  else if ( y < 12.0D+00 ) then

    y1 = y
!
!  0.0 < argument < 1.0.
!
    if ( y < 1.0D+00 ) then

      z = y
      y = y + 1.0D+00
!
!  1.0 < argument < 12.0.
!  Reduce argument if necessary.
!
    else

      n = int ( y ) - 1
      y = y - real ( n)
      z = y - 1.0D+00

    end if
!
!  Evaluate approximation for 1.0 < argument < 2.0.
!
    xnum = 0.0D+00
    xden = 1.0D+00
    do i = 1, 8
      xnum = ( xnum + p(i) ) * z
      xden = xden * z + q(i)
    end do

    res = xnum / xden + 1.0D+00
!
!  Adjust result for case  0.0 < argument < 1.0.
!
    if ( y1 < y ) then

      res = res / y1
!
!  Adjust result for case 2.0 < argument < 12.0.
!
    else if ( y < y1 ) then

      do i = 1, n
        res = res * y
        y = y + 1.0D+00
      end do

    end if

  else
!
!  Evaluate for 12.0 <= argument.
!
    if ( y <= xbig ) then

      ysq = y * y
      sum = c(7)
      do i = 1, 6
        sum = sum / ysq + c(i)
      end do
      sum = sum / y - y + sqrtpi
      sum = sum + ( y - 0.5D+00 ) * log ( y )
      res = exp ( sum )

    else

      res = xinf
      r8_gamma = res
    end if

  end if
!
!  Final adjustments and return.
!
  if ( parity ) then
    res = - res
  end if

  if ( fact /= 1.0D+00 ) then
    res = fact / res
  end if

  r8_gamma = res
end

function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a double precision value.
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
!    Input/output, integer SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer , parameter :: i4_huge = 2147483647
  integer k
  double precision r8_uniform_01
  integer seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed) * 4.656612875D-10
end

function r8vec_norm ( n, a )

!*****************************************************************************80
!
!! R8VEC_NORM returns the L2 norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The vector L2 norm is defined as:
!
!      R8VEC_NORM = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in A.
!
!    Input, double precision A(N), the vector whose L2 norm is desired.
!
!    Output, double precision R8VEC_NORM, the L2 norm of A.
!
  implicit none

  integer n

  double precision a(n)
  double precision r8vec_norm

  r8vec_norm = sqrt ( sum ( a(1:n)**2 ) )
end

subroutine r8vec_polarize ( n, a, p, a_normal, a_parallel )

!*****************************************************************************80
!
!! R8VEC_POLARIZE decomposes an R8VEC into normal and parallel components.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The (nonzero) vector P defines a direction.
!
!    The vector A can be written as the sum
!
!      A = A_normal + A_parallel
!
!    where A_parallel is a linear multiple of P, and A_normal
!    is perpendicular to P.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, double precision A(N), the vector to be polarized.
!
!    Input, double precision P(N), the polarizing direction.
!
!    Output, double precision A_NORMAL(N), A_PARALLEL(N), the normal
!    and parallel components of A.
!
  implicit none

  integer n

  double precision a(n)
  double precision a_dot_p
  double precision a_normal(n)
  double precision a_parallel(n)
  double precision p(n)
  double precision p_norm

  p_norm = sqrt ( sum ( p(1:n)**2 ) )

  if ( p_norm == 0.0D+00 ) then
    a_normal(1:n) = a(1:n)
    a_parallel(1:n) = 0.0D+00
  end if

  a_dot_p = dot_product ( a(1:n), p(1:n) ) / p_norm

  a_parallel(1:n) = a_dot_p * p(1:n) / p_norm

  a_normal(1:n) = a(1:n) - a_parallel(1:n)
end

subroutine s_cat ( s1, s2, s3 )

!*****************************************************************************80
!
!! S_CAT concatenates two strings to make a third string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, the "prefix" string.
!
!    Input, character ( len = * ) S2, the "postfix" string.
!
!    Output, character ( len = * ) S3, the string made by
!    concatenating S1 and S2, ignoring any trailing blanks.
!
  implicit none

  character ( len = * ) s1
  character ( len = * ) s2
  character ( len = * ) s3

  if ( s1 == ' ' .and. s2 == ' ' ) then
    s3 = ' '
  else if ( s1 == ' ' ) then
    s3 = s2
  else if ( s2 == ' ' ) then
    s3 = s1
  else
    s3 = trim ( s1 ) // trim ( s2 )
  end if
end

subroutine sphere01_distance_xyz ( xyz1, xyz2, dist )

!*****************************************************************************80
!
!! SPHERE01_DISTANCE_XYZ computes great circle distances on a unit sphere.
!
!  Discussion:
!
!    XYZ coordinates are used.
!
!    We assume the points XYZ1 and XYZ2 lie on the unit sphere.
!
!    This computation is a special form of the Vincenty formula.
!    It should be less sensitive to errors associated with very small 
!    or very large angular separations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    "Great-circle distance",
!    Wikipedia.
!
!  Parameters:
!
!    Input, double precision XYZ1(3), the coordinates of the first point.
!
!    Input, double precision XYZ2(3), the coordinates of the second point.
!
!    Output, double precision DIST, the great circle distance between
!    the points.
!
  implicit none

  double precision arc_sine
  double precision atan4
  double precision bot
  double precision dist
  double precision lat1
  double precision lat2
  double precision lon1
  double precision lon2
  double precision top
  double precision xyz1(3)
  double precision xyz2(3)

  lat1 = arc_sine ( xyz1(3) )
  lon1 = atan4 ( xyz1(2), xyz1(1) )

  lat2 = arc_sine ( xyz2(3) )
  lon2 = atan4 ( xyz2(2), xyz2(1) )

  top = ( cos ( lat2 ) * sin ( lon1 - lon2 ) )**2 &
      + ( cos ( lat1 ) * sin ( lat2 ) &
      -   sin ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 ) )**2

  top = sqrt ( top )

  bot = sin ( lat1 ) * sin ( lat2 ) &
      + cos ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 )

  dist = atan2 ( top, bot )
end

subroutine sphere01_monomial_integral ( e, integral )

!*****************************************************************************80
!
!! SPHERE01_MONOMIAL_INTEGRAL returns monomial integrals on the unit sphere.
!
!  Discussion:
!
!    The integration region is 
!
!      X^2 + Y^2 + Z^2 = 1.
!
!    The monomial is F(X,Y,Z) = X^E(1) * Y^E(2) * Z^E(3).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Academic Press, 1984, page 263.
!
!  Parameters:
!
!    Input, integer E(3), the exponents of X, Y and Z in the 
!    monomial.  Each exponent must be nonnegative.
!
!    Output, double precision INTEGRAL, the integral.
!
  implicit none

  integer e(3)
  integer i
  double precision integral
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision r8_gamma

  if ( any ( e(1:3) < 0 ) ) then
    integral = - huge ( integral )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPHERE01_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  All exponents must be nonnegative.'
    write ( *, '(a,i8)' ) '  E(1) = ', e(1)
    write ( *, '(a,i8)' ) '  E(2) = ', e(2)
    write ( *, '(a,i8)' ) '  E(3) = ', e(3)
    stop
  end if

  if ( all ( e(1:3) == 0 ) ) then

    integral = 2.0D+00 * sqrt ( pi**3 ) / r8_gamma ( 1.5D+00 )

  else if ( any ( mod ( e(1:3), 2 ) == 1 ) ) then

    integral = 0.0D+00

  else

    integral = 2.0D+00

    do i = 1, 3
      integral = integral * r8_gamma ( 0.5D+00 * real ( e(i) + 1) )
    end do

    integral = integral &
      / r8_gamma ( 0.5D+00 * ( real ( sum ( e(1:3) + 1 )) ) )

  end if
end

subroutine sphere01_quad_icos1c ( factor, fun, node_num, result )

!*****************************************************************************80
!
!! SPHERE01_QUAD_ICOS1C: centroid rule, subdivide then project.
!
!  Discussion:
!
!    This function estimates an integral over the surface of the unit sphere.
!
!    This function sets up an icosahedral grid, and subdivides each
!    edge of the icosahedron into FACTOR subedges.  These edges define a grid
!    within each triangular icosahedral face.  The centroids of these
!    triangles can be determined.  All of these calculations are done,
!    essentially, on the FLAT faces of the icosahedron.  Only then are
!    the triangle vertices and centroids projected to the sphere.  
!
!    The resulting grid of spherical triangles and projected centroids
!    is used to apply a centroid quadrature rule over the surface of
!    the unit sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer FACTOR, the subdivision factor, which must
!    be at least 1.
!
!    Input, external :: FUN, evaluates the integrand, of the form:
!      subroutine fun ( n, x, v )
!      integer n
!      double precision v(n)
!      double precision x(3,n)
!
!    Output, integer NODE_NUM, the number of evaluation points.
!
!    Output, double precision RESULT, the estimated integral.
!
  implicit none

  integer node_num

  integer a
  double precision a_xyz(3)
  double precision a2_xyz(3)
  double precision area
  double precision area_total
  integer b
  double precision b_xyz(3)
  double precision b2_xyz(3)
  integer c
  double precision c_xyz(3)
  double precision c2_xyz(3)
  integer edge_num
  integer , allocatable, dimension ( :, : ) :: edge_point
  integer f1
  integer f2
  integer f3
  integer face
  integer face_num
  integer , allocatable, dimension ( : ) :: face_order
  integer , allocatable, dimension ( :, : ) :: face_point
  integer face_order_max
  integer factor
  external             fun
  double precision node_xyz(3)
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision , allocatable, dimension ( :, : ) :: point_coord
  integer point_num
  double precision result
  double precision v
!
!  Size the icosahedron.
!
  call icos_size ( point_num, edge_num, face_num, face_order_max )
!
!  Set the icosahedron.
!
  allocate ( point_coord(1:3,1:point_num) )
  allocate ( edge_point(1:2,1:edge_num) )
  allocate ( face_order(1:face_num) )
  allocate ( face_point(1:face_order_max,1:face_num) )

  call icos_shape ( point_num, edge_num, face_num, face_order_max, &
    point_coord, edge_point, face_order, face_point )
!
!  Initialize the integral data.
!
  result = 0.0D+00
  area_total = 0.0D+00
  node_num = 0
!
!  Pick a face of the icosahedron, and identify its vertices as A, B, C.
!
  do face = 1, face_num

    a = face_point(1,face)
    b = face_point(2,face)
    c = face_point(3,face)

    a_xyz(1:3) = point_coord(1:3,a)
    b_xyz(1:3) = point_coord(1:3,b)
    c_xyz(1:3) = point_coord(1:3,c)
!
!  Some subtriangles will have the same direction as the face.
!  Generate each in turn, by determining the barycentric coordinates
!  of the centroid (F1,F2,F3), from which we can also work out the barycentric
!  coordinates of the vertices of the subtriangle.
!
    do f3 = 1, 3 * factor - 2, 3
      do f2 = 1, 3 * factor - f3 - 1, 3

        f1 = 3 * factor - f3 - f2

        call sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1, f2, f3, &
          node_xyz )

        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 + 2, f2 - 1, f3 - 1, a2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 - 1, f2 + 2, f3 - 1, b2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 - 1, f2 - 1, f3 + 2, c2_xyz )

        call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, area )

        call fun ( 1, node_xyz, v )    

        node_num = node_num + 1
        result = result + area * v
        area_total = area_total + area

      end do
    end do
!
!  The other subtriangles have the opposite direction from the face.
!  Generate each in turn, by determining the barycentric coordinates
!  of the centroid (F1,F2,F3), from which we can also work out the barycentric
!  coordinates of the vertices of the subtriangle.
!
    do f3 = 2, 3 * factor - 4, 3
      do f2 = 2, 3 * factor - f3 - 2, 3

        f1 = 3 * factor - f3 - f2

        call sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1, f2, f3, &
          node_xyz )

        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 - 2, f2 + 1, f3 + 1, a2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 + 1, f2 - 2, f3 + 1, b2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 + 1, f2 + 1, f3 - 2, c2_xyz )

        call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, area )

        call fun ( 1, node_xyz, v )  

        node_num = node_num + 1  
        result = result + area * v
        area_total = area_total + area

      end do
    end do

  end do
!
!  Discard allocated memory.
!
  deallocate ( edge_point )
  deallocate ( face_order )
  deallocate ( face_point )
  deallocate ( point_coord )
end

subroutine sphere01_quad_icos1m ( factor, fun, node_num, result )

!*****************************************************************************80
!
!! SPHERE01_QUAD_ICOS1M: midside rule, subdivide then project.
!
!  Discussion:
!
!    This function estimates an integral over the surface of the unit sphere.
!
!    This function sets up an icosahedral grid, and subdivides each
!    edge of the icosahedron into FACTOR subedges.  These edges define a grid
!    within each triangular icosahedral face.  The midsides of these
!    triangles can be determined.  All of these calculations are done,
!    essentially, on the FLAT faces of the icosahedron.  Only then are
!    the triangle vertices and midsides projected to the sphere.  
!
!    The resulting grid of spherical triangles and projected midsides
!    is used to apply a midside quadrature rule over the surface of
!    the unit sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer FACTOR, the subdivision factor, which must
!    be at least 1.
!
!    Input, external :: FUN, evaluates the integrand, of the form:
!      subroutine fun ( n, x, v )
!      integer n
!      double precision v(n)
!      double precision x(3,n)
!
!    Output, integer NODE_NUM, the number of evaluation points.
!
!    Output, double precision RESULT, the estimated integral.
!
  implicit none

  integer node_num

  integer a
  double precision a_xyz(3)
  double precision a2_xyz(3)
  double precision a3_xyz(3)
  double precision area
  double precision area_total
  integer b
  double precision b_xyz(3)
  double precision b2_xyz(3)
  double precision b3_xyz(3)
  integer c
  double precision c_xyz(3)
  double precision c2_xyz(3)
  double precision c3_xyz(3)
  integer edge
  integer edge_num
  integer , allocatable, dimension ( :, : ) :: edge_point
  integer f
  integer f1
  integer f2
  integer f3
  integer face
  integer face_num
  integer , allocatable, dimension ( : ) :: face_order
  integer , allocatable, dimension ( :, : ) :: face_point
  integer face_order_max
  integer factor
  external fun
  integer j
  integer node
  double precision node_norm
  double precision node_xyz(3)
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision , allocatable, dimension ( :, : ) :: point_coord
  integer point_num
  double precision result
  double precision va
  double precision vb
  double precision vc
!
!  Size the icosahedron.
!
  call icos_size ( point_num, edge_num, face_num, face_order_max )
!
!  Set the icosahedron.
!
  allocate ( point_coord(1:3,1:point_num) )
  allocate ( edge_point(1:2,1:edge_num) )
  allocate ( face_order(1:face_num) )
  allocate ( face_point(1:face_order_max,1:face_num) )

  call icos_shape ( point_num, edge_num, face_num, face_order_max, &
    point_coord, edge_point, face_order, face_point )
!
!  Initialize the integral data.
!
  result = 0.0D+00
  node_num = 0
  area_total = 0.0D+00
!
!  Consider each face.
!
  do face = 1, face_num

    a = face_point(1,face)
    b = face_point(2,face)
    c = face_point(3,face)

    a_xyz(1:3) = point_coord(1:3,a)
    b_xyz(1:3) = point_coord(1:3,b)
    c_xyz(1:3) = point_coord(1:3,c)
!
!  Deal with subtriangles that have same orientation as face.
!
    do f1 = 0, factor - 1
      do f2 = 0, factor - f1 - 1
        f3 = factor - f1 - f2

        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 + 1, f2,     f3 - 1, a2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1,     f2 + 1, f3 - 1, b2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1,     f2,     f3,     c2_xyz )

        call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, area )

        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, 2 * f1 + 1, 2 * f2 + 1, 2 * f3 - 2, a3_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, 2 * f1,     2 * f2 + 1, 2 * f3 - 1, b3_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, 2 * f1 + 1, 2 * f2,     2 * f3 - 1, c3_xyz )

        node_num = node_num + 3
        call fun ( 1, a3_xyz, va )
        call fun ( 1, b3_xyz, vb )   
        call fun ( 1, c3_xyz, vc )   
        result = result + area * ( va + vb + vc ) / 3.0D+00
        area_total = area_total + area

      end do
    end do
!
!  Deal with subtriangles that have opposite orientation as face.
!
    do f3 = 0, factor - 2
      do f2 = 1, factor - f3 - 1
        f1 = factor - f2 - f3

        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 - 1, f2,     f3 + 1, a2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1,     f2 - 1, f3 + 1, b2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1,     f2,     f3,     c2_xyz )

        call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, area )

        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, 2 * f1 - 1, 2 * f2 - 1, 2 * f3 + 2, a3_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, 2 * f1,     2 * f2 - 1, 2 * f3 + 1, b3_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, 2 * f1 - 1, 2 * f2,     2 * f3 + 1, c3_xyz )

        node_num = node_num + 3
        call fun ( 1, a3_xyz, va )
        call fun ( 1, b3_xyz, vb )   
        call fun ( 1, c3_xyz, vc )   
        result = result + area * ( va + vb + vc ) / 3.0D+00
        area_total = area_total + area

      end do
    end do
  end do
!
!  Discard allocated memory.
!
  deallocate ( edge_point )
  deallocate ( face_order )
  deallocate ( face_point )
  deallocate ( point_coord )
end

subroutine sphere01_quad_icos1v ( factor, fun, node_num, result )

!*****************************************************************************80
!
!! SPHERE01_QUAD_ICOS1V: vertex rule, subdivide then project.
!
!  Discussion:
!
!    This function estimates an integral over the surface of the unit sphere.
!
!    This function sets up an icosahedral grid, and subdivides each
!    edge of the icosahedron into FACTOR subedges.  These edges define a grid
!    within each triangular icosahedral face.  The vertices of these
!    triangles can be determined.  All of these calculations are done,
!    essentially, on the FLAT faces of the icosahedron.  Only then are
!    the triangle vertices projected to the sphere.  
!
!    The resulting grid of spherical triangles is used to apply a vertex
!    quadrature rule over the surface of the unit sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer FACTOR, the subdivision factor, which must
!    be at least 1.
!
!    Input, external :: FUN, evaluates the integrand, of the form:
!      subroutine fun ( n, x, v )
!      integer n
!      double precision v(n)
!      double precision x(3,n)
!
!    Output, integer NODE_NUM, the number of evaluation points.
!
!    Output, double precision RESULT, the estimated integral.
!
  implicit none

  integer node_num

  integer a
  double precision a_xyz(3)
  double precision a2_xyz(3)
  double precision area
  double precision area_total
  integer b
  double precision b_xyz(3)
  double precision b2_xyz(3)
  integer c
  double precision c_xyz(3)
  double precision c2_xyz(3)
  integer edge
  integer edge_num
  integer , allocatable, dimension ( :, : ) :: edge_point
  integer f
  integer f1
  integer f2
  integer f3
  integer face
  integer face_num
  integer , allocatable, dimension ( : ) :: face_order
  integer , allocatable, dimension ( :, : ) :: face_point
  integer face_order_max
  integer factor
  external fun
  integer j
  integer node
  double precision node_norm
  double precision node_xyz(3)
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision , allocatable, dimension ( :, : ) :: point_coord
  integer point_num
  double precision result
  double precision va
  double precision vb
  double precision vc
!
!  Size the icosahedron.
!
  call icos_size ( point_num, edge_num, face_num, face_order_max )
!
!  Set the icosahedron.
!
  allocate ( point_coord(1:3,1:point_num) )
  allocate ( edge_point(1:2,1:edge_num) )
  allocate ( face_order(1:face_num) )
  allocate ( face_point(1:face_order_max,1:face_num) )

  call icos_shape ( point_num, edge_num, face_num, face_order_max, &
    point_coord, edge_point, face_order, face_point )
!
!  Initialize the integral data.
!
  result = 0.0D+00
  node_num = 0
  area_total = 0.0D+00
!
!  Consider each face.
!
  do face = 1, face_num

    a = face_point(1,face)
    b = face_point(2,face)
    c = face_point(3,face)

    a_xyz(1:3) = point_coord(1:3,a)
    b_xyz(1:3) = point_coord(1:3,b)
    c_xyz(1:3) = point_coord(1:3,c)
!
!  Deal with subtriangles that have same orientation as face.
!
    do f1 = 0, factor - 1
      do f2 = 0, factor - f1 - 1
        f3 = factor - f1 - f2

        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 + 1, f2,     f3 - 1, a2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1,     f2 + 1, f3 - 1, b2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1,     f2,     f3,     c2_xyz )

        call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, area )

        node_num = node_num + 3
        call fun ( 1, a2_xyz, va )
        call fun ( 1, b2_xyz, vb )   
        call fun ( 1, c2_xyz, vc )   
        result = result + area * ( va + vb + vc ) / 3.0D+00
        area_total = area_total + area

      end do
    end do
!
!  Deal with subtriangles that have opposite orientation as face.
!
    do f3 = 0, factor - 2
      do f2 = 1, factor - f3 - 1
        f1 = factor - f2 - f3

        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 - 1, f2,     f3 + 1, a2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1,     f2 - 1, f3 + 1, b2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1,     f2,     f3,     c2_xyz )

        call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, area )

        node_num = node_num + 3
        call fun ( 1, a2_xyz, va )
        call fun ( 1, b2_xyz, vb )   
        call fun ( 1, c2_xyz, vc )   
        result = result + area * ( va + vb + vc ) / 3.0D+00
        area_total = area_total + area

      end do
    end do
  end do
!
!  Discard allocated memory.
!
  deallocate ( edge_point )
  deallocate ( face_order )
  deallocate ( face_point )
  deallocate ( point_coord )
end

subroutine sphere01_quad_icos2v ( factor, fun, node_num, result )

!*****************************************************************************80
!
!! SPHERE01_QUAD_ICOS2V: vertex rule, subdivide then project.
!
!  Discussion:
!
!    This function estimates an integral over the surface of the unit sphere.
!
!    This function sets up an icosahedral grid, and subdivides each
!    edge of the icosahedron into FACTOR subedges.  These edges define a grid
!    within each triangular icosahedral face.  The vertices of these
!    triangles can be determined.  All of these calculations are done,
!    essentially, on the FLAT faces of the icosahedron.  Only then are
!    the triangle vertices projected to the sphere.  
!
!    The resulting grid of spherical triangles is used to apply a vertex
!    quadrature rule over the surface of the unit sphere.
!
!    This is a revision of SPHERE01_QUAD_ICOS2V that attempted to use a more
!    sophisticated scheme to map points from the planar triangle to the surface
!    of the unit sphere.  Very little improvement to the estimated integral
!    was observed, so development of this scheme has been set aside for now.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer FACTOR, the subdivision factor, which must
!    be at least 1.
!
!    Input, external :: FUN, evaluates the integrand, of the form:
!      subroutine fun ( n, x, v )
!      integer n
!      double precision v(n)
!      double precision x(3,n)
!
!    Output, integer NODE_NUM, the number of evaluation points.
!
!    Output, double precision RESULT, the estimated integral.
!
  implicit none

  integer node_num

  integer a
  double precision a_xyz(3)
  double precision a2_xyz(3)
  double precision area
  double precision area_total
  integer b
  double precision b_xyz(3)
  double precision b2_xyz(3)
  integer c
  double precision c_xyz(3)
  double precision c2_xyz(3)
  integer edge
  integer edge_num
  integer , allocatable, dimension ( :, : ) :: edge_point
  integer f
  integer f1
  integer f2
  integer f3
  integer face
  integer face_num
  integer , allocatable, dimension ( : ) :: face_order
  integer , allocatable, dimension ( :, : ) :: face_point
  integer face_order_max
  integer factor
  external fun
  integer j
  integer node
  double precision node_norm
  double precision node_xyz(3)
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision , allocatable, dimension ( :, : ) :: point_coord
  integer point_num
  double precision result
  double precision va
  double precision vb
  double precision vc
!
!  Size the icosahedron.
!
  call icos_size ( point_num, edge_num, face_num, face_order_max )
!
!  Set the icosahedron.
!
  allocate ( point_coord(1:3,1:point_num) )
  allocate ( edge_point(1:2,1:edge_num) )
  allocate ( face_order(1:face_num) )
  allocate ( face_point(1:face_order_max,1:face_num) )

  call icos_shape ( point_num, edge_num, face_num, face_order_max, &
    point_coord, edge_point, face_order, face_point )
!
!  Initialize the integral data.
!
  result = 0.0D+00
  node_num = 0
  area_total = 0.0D+00
!
!  Consider each face.
!
  do face = 1, face_num

    a = face_point(1,face)
    b = face_point(2,face)
    c = face_point(3,face)

    a_xyz(1:3) = point_coord(1:3,a)
    b_xyz(1:3) = point_coord(1:3,b)
    c_xyz(1:3) = point_coord(1:3,c)
!
!  Deal with subtriangles that have same orientation as face.
!
    do f1 = 0, factor - 1
      do f2 = 0, factor - f1 - 1
        f3 = factor - f1 - f2

        call sphere01_triangle_project2 ( &
          a_xyz, b_xyz, c_xyz, f1 + 1, f2,     f3 - 1, a2_xyz )
        call sphere01_triangle_project2 ( &
          a_xyz, b_xyz, c_xyz, f1,     f2 + 1, f3 - 1, b2_xyz )
        call sphere01_triangle_project2 ( &
          a_xyz, b_xyz, c_xyz, f1,     f2,     f3,     c2_xyz )

        call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, area )

        node_num = node_num + 3
        call fun ( 1, a2_xyz, va )
        call fun ( 1, b2_xyz, vb )   
        call fun ( 1, c2_xyz, vc )   
        result = result + area * ( va + vb + vc ) / 3.0D+00
        area_total = area_total + area

      end do
    end do
!
!  Deal with subtriangles that have opposite orientation as face.
!
    do f3 = 0, factor - 2
      do f2 = 1, factor - f3 - 1
        f1 = factor - f2 - f3

        call sphere01_triangle_project2 ( &
          a_xyz, b_xyz, c_xyz, f1 - 1, f2,     f3 + 1, a2_xyz )
        call sphere01_triangle_project2 ( &
          a_xyz, b_xyz, c_xyz, f1,     f2 - 1, f3 + 1, b2_xyz )
        call sphere01_triangle_project2 ( &
          a_xyz, b_xyz, c_xyz, f1,     f2,     f3,     c2_xyz )

        call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, area )

        node_num = node_num + 3
        call fun ( 1, a2_xyz, va )
        call fun ( 1, b2_xyz, vb )   
        call fun ( 1, c2_xyz, vc )   
        result = result + area * ( va + vb + vc ) / 3.0D+00
        area_total = area_total + area

      end do
    end do
  end do
!
!  Discard allocated memory.
!
  deallocate ( edge_point )
  deallocate ( face_order )
  deallocate ( face_point )
  deallocate ( point_coord )
end

subroutine sphere01_quad_llc ( f, h, n, result )

!*****************************************************************************80
!
!! SPHERE01_QUAD_LLC: Longitude/Latitude grid with centroid rule.
!
!  Discussion:
!
!    The sphere is broken up into spherical triangles, whose sides
!    do not exceed the length H.  Then a centroid rule is used on
!    each spherical triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external :: F, evaluates the integrand, of the form:
!      subroutine f ( n, x, v )
!      integer n
!      double precision v(n)
!      double precision x(3,n)
!
!    Input, double precision H, the maximum length of a side of the spherical
!    quadrilaterals.
!
!    Output, integer N, the number of points used.
!
!    Output, double precision RESULT, the approximate integral.
!
  implicit none

  double precision area
  external f
  double precision h
  integer i
  integer j
  integer n
  double precision phi
  integer phi_num
  double precision phi1
  double precision phi2
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision result
  double precision sector_area
  double precision sphere_area
  double precision theta
  integer theta_num
  double precision theta1
  double precision theta2
  double precision v(1)
  double precision x(3)
  double precision x1(3)
  double precision x11(3)
  double precision x12(3)
  double precision x2(3)
  double precision x21(3)
  double precision x22(3)
!
!  Choose PHI and THETA counts that make short sides.
!
  phi_num = int ( pi / h )

  if ( h * real ( phi_num) < pi ) then
    phi_num = phi_num + 1
  end if

  theta_num = int ( 2.0D+00 * pi / h )

  if ( h * real ( theta_num) < pi ) then
    theta_num = theta_num + 1
  end if

  n = 0
  result = 0.0D+00
!
!  Only one THETA (and hence, only one PHI.)
!
  if ( theta_num == 1 ) then

    sphere_area = 4.0D+00 * pi

    theta = 0.0D+00
    phi = pi / 2.0D+00
    call tp_to_xyz ( theta, phi, x )

    call f ( 1, x, v )
    n = n + 1
    result = sphere_area * v(1)
!
!  Several THETA's, one PHI.
!
  else if ( phi_num == 1 ) then

    sphere_area = 4.0D+00 * pi
    sector_area = sphere_area / real ( theta_num)

    result = 0.0D+00

    do j = 1, theta_num

      theta = real ( ( j - 1 ) * 2) * pi &
            / real ( theta_num)
      phi = pi / 2.0D+00
      call tp_to_xyz ( theta, phi, x )
      call f ( 1, x, v )
      n = n + 1
      result = result + sector_area * v(1)

    end do
!
!  At least two PHI's.
!
  else

    result = 0.0D+00
!
!  Picture in top row, with V1 = north pole:
!
!        V1
!       /  \
!      /    \
!    V12----V22
!
    phi1 = 0.0D+00
    phi2 = pi / real ( phi_num)

    do j = 1, theta_num

      theta1 = real ( j - 1) * 2.0D+00 * pi &
             / real ( theta_num)
      theta2 = real ( j    ) * 2.0D+00 * pi &
             / real ( theta_num)

      call tp_to_xyz ( theta1, phi1, x1 )
      call tp_to_xyz ( theta1, phi2, x12 )
      call tp_to_xyz ( theta2, phi2, x22 )

      call sphere01_triangle_vertices_to_area ( x1, x12, x22, area )
      call sphere01_triangle_vertices_to_centroid ( x1, x12, x22, x )
      call f ( 1, x, v )
      n = n + 1
      result = result + area * v(1)

    end do
!
!  Picture in all intermediate rows:
!
!    V11--V21
!     | \  |
!     |  \ |
!    V12--V22
!
    do i = 2, phi_num-1

      phi1 = real ( i - 1) * pi / real ( phi_num)
      phi2 = real ( i) * pi / real ( phi_num)

      do j = 1, theta_num

        theta1 = real ( j - 1) * 2.0D+00 * pi &
               / real ( theta_num)
        theta2 = real ( j) * 2.0D+00 * pi &
               / real ( theta_num)

        call tp_to_xyz ( theta1, phi1, x11 )
        call tp_to_xyz ( theta2, phi1, x21 )
        call tp_to_xyz ( theta1, phi2, x12 )
        call tp_to_xyz ( theta2, phi2, x22 )

        call sphere01_triangle_vertices_to_area ( x11, x12, x22, area )
        call sphere01_triangle_vertices_to_centroid ( x11, x12, x22, x )
        call f ( 1, x, v )
        n = n + 1
        result = result + area * v(1)

        call sphere01_triangle_vertices_to_area ( x22, x21, x11, area )
        call sphere01_triangle_vertices_to_centroid ( x22, x21, x11, x )
        call f ( 1, x, v )
        n = n + 1
        result = result + area * v(1)

      end do

    end do
!
!  Picture in last row, with V2 = south pole:
!
!    V11----V21
!      \    /
!       \  /
!        V2
!
    phi1 = real ( phi_num - 1) * pi &
         / real ( phi_num)
    phi2 =                                  pi

    do j = 1, theta_num

      theta1 = real ( j - 1) * 2.0D+00 * pi &
             / real ( theta_num)
      theta2 = real ( j) * 2.0D+00 * pi &
             / real ( theta_num)

      call tp_to_xyz ( theta1, phi1, x11 )
      call tp_to_xyz ( theta2, phi1, x21 )
      call tp_to_xyz ( theta2, phi2, x2 )

      call sphere01_triangle_vertices_to_area ( x11, x2, x21, area )
      call sphere01_triangle_vertices_to_centroid ( x11, x2, x21, x )
      call f ( 1, x, v )
      n = n + 1
      result = result + area * v(1)

    end do

  end if
end

subroutine sphere01_quad_llm ( f, h, n, result )

!*****************************************************************************80
!
!! SPHERE01_QUAD_LLM: longitude/latitude grid plus midside rule.
!
!  Discussion:
!
!    The sphere is broken up into spherical triangles, whose sides
!    do not exceed the length H.  Then the function is evaluated
!    at the midsides, and the average is multiplied by the area.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external :: F, evaluates the integrand, of the form:
!      subroutine f ( n, x, v )
!      integer n
!      double precision v(n)
!      double precision x(3,n)
!
!    Input, double precision H, the maximum length of a side of the spherical
!    quadrilaterals.
!
!    Output, integer N, the number of points used.
!
!    Output, double precision RESULT, the approximate integral.
!
  implicit none

  double precision area
  external f
  double precision h
  integer i
  integer j
  double precision m1(3)
  double precision m2(3)
  double precision m3(3)
  integer n
  double precision phi
  integer phi_num
  double precision phi1
  double precision phi2
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision result
  double precision sector_area
  double precision sphere_area
  double precision theta
  integer theta_num
  double precision theta1
  double precision theta2
  double precision v(1)
  double precision x(3)
  double precision x1(3)
  double precision x11(3)
  double precision x12(3)
  double precision x2(3)
  double precision x21(3)
  double precision x22(3)
!
!  Choose PHI and THETA counts that make short sides.
!
  phi_num = int ( pi / h )

  if ( h * real ( phi_num) < pi ) then
    phi_num = phi_num + 1
  end if

  theta_num = int ( 2.0D+00 * pi / h )

  if ( h * real ( theta_num) < pi ) then
    theta_num = theta_num + 1
  end if

  n = 0
  result = 0.0D+00
!
!  Only one THETA (and hence, only one PHI.)
!
  if ( theta_num == 1 ) then

    sphere_area = 4.0D+00 * pi

    theta = 0.0D+00
    phi = pi / 2.0D+00
    call tp_to_xyz ( theta, phi, x )
    call f ( 1, x, v )
    n = n + 1
    result = sphere_area * v(1)
!
!  Several THETA's, one PHI.
!
  else if ( phi_num == 1 ) then

    sphere_area = 4.0D+00 * pi
    sector_area = sphere_area / real ( theta_num)

    result = 0.0D+00

    do j = 1, theta_num

      theta = real ( ( j - 1 ) * 2) * pi &
            / real ( theta_num)
      phi = pi / 2.0D+00
      call tp_to_xyz ( theta, phi, x )
      call f ( 1, x, v )
      n = n + 1
      result = result + sector_area * v(1)

    end do
!
!  At least two PHI's.
!
  else

    result = 0.0D+00
!
!  Picture:
!
!        V1
!       /  \
!      /    \
!    V12----V22
!
    phi1 = 0.0D+00
    phi2 = pi / real ( phi_num)

    do j = 1, theta_num

      theta1 = real ( j - 1) * 2.0D+00 * pi &
             / real ( theta_num)
      theta2 = real ( j) * 2.0D+00 * pi &
             / real ( theta_num)

      call tp_to_xyz ( theta1, phi1, x1 )
      call tp_to_xyz ( theta1, phi2, x12 )
      call tp_to_xyz ( theta2, phi2, x22 )

      call sphere01_triangle_vertices_to_area ( x1, x12, x22, area )

      call sphere01_triangle_vertices_to_midpoints ( x1, x12, x22, &
        m1, m2, m3 )

      call f ( 1, m1, v )
      n = n + 1
      result = result + area * v(1) / 3.0D+00
      call f ( 1, m2, v )
      n = n + 1
      result = result + area * v(1) / 3.0D+00
      call f ( 1, m3, v )
      n = n + 1
      result = result + area * v(1) / 3.0D+00

    end do
!
!  Picture:
!
!    V11--V21
!     | \  |
!     |  \ |
!    V12--V22
!
    do i = 2, phi_num-1

      phi1 = real ( i - 1) * pi &
           / real ( phi_num)
      phi2 = real ( i) * pi &
           / real ( phi_num)

      do j = 1, theta_num

        theta1 = real ( j - 1) * 2.0D+00 * pi &
               / real ( theta_num)
        theta2 = real ( j) * 2.0D+00 * pi &
               / real ( theta_num)

        call tp_to_xyz ( theta1, phi1, x11 )
        call tp_to_xyz ( theta2, phi1, x21 )
        call tp_to_xyz ( theta1, phi2, x12 )
        call tp_to_xyz ( theta2, phi2, x22 )

        call sphere01_triangle_vertices_to_area ( x11, x12, x22, area )

        call sphere01_triangle_vertices_to_midpoints ( x11, x12, x22, &
          m1, m2, m3 )

        call f ( 1, m1, v )
        n = n + 1
        result = result + area * v(1) / 3.0D+00
        call f ( 1, m2, v )
        n = n + 1
        result = result + area * v(1) / 3.0D+00
        call f ( 1, m3, v )
        n = n + 1
        result = result + area * v(1) / 3.0D+00

        call sphere01_triangle_vertices_to_area ( x22, x21, x11, area )

        call sphere01_triangle_vertices_to_midpoints ( x22, x21, x11, &
          m1, m2, m3 )

        call f ( 1, m1, v )
        n = n + 1
        result = result + area * v(1) / 3.0D+00
        call f ( 1, m2, v )
        n = n + 1
        result = result + area * v(1) / 3.0D+00
        call f ( 1, m3, v )
        n = n + 1
        result = result + area * v(1) / 3.0D+00

      end do

    end do
!
!  Picture:
!
!    V11----V21
!      \    /
!       \  /
!        V2
!
    phi1 = real ( phi_num - 1) * pi &
         / real ( phi_num)
    phi2 =                                  pi

    do j = 1, theta_num

      theta1 = real ( j - 1) * 2.0D+00 * pi &
             / real ( theta_num)
      theta2 = real ( j) * 2.0D+00 * pi &
             / real ( theta_num)

      call tp_to_xyz ( theta1, phi1, x11 )
      call tp_to_xyz ( theta2, phi1, x21 )
      call tp_to_xyz ( theta2, phi2, x2 )

      call sphere01_triangle_vertices_to_area ( x11, x2, x21, area )

      call sphere01_triangle_vertices_to_midpoints ( x11, x2, x21, &
        m1, m2, m3 )

      call f ( 1, m1, v )
      n = n + 1
      result = result + area * v(1) / 3.0D+00
      call f ( 1, m2, v )
      n = n + 1
      result = result + area * v(1) / 3.0D+00
      call f ( 1, m3, v )
      n = n + 1
      result = result + area * v(1) / 3.0D+00

    end do

  end if
end

subroutine sphere01_quad_llv ( f, h, n, result )

!*****************************************************************************80
!
!! SPHERE01_QUAD_LLV: longitude/latitude grid with vertex rule.
!
!  Discussion:
!
!    The sphere is broken up into spherical triangles, whose sides
!    do not exceed the length H.  Then the function is evaluated
!    at the vertices, and the average is multiplied by the area.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external :: F, evaluates the integrand, of the form:
!      subroutine f ( n, x, v )
!      integer n
!      double precision v(n)
!      double precision x(3,n)
!
!    Input, double precision H, the maximum length of a side of the spherical
!    quadrilaterals.
!
!    Output, integer N, the number of points used.
!
!    Output, double precision RESULT, the approximate integral.
!
  implicit none

  double precision area
  external f
  double precision h
  integer i
  integer j
  integer n
  double precision phi
  integer phi_num
  double precision phi1
  double precision phi2
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision result
  double precision sector_area
  double precision sphere_area
  double precision theta
  integer theta_num
  double precision theta1
  double precision theta2
  double precision v(1)
  double precision x(3)
  double precision x1(3)
  double precision x11(3)
  double precision x12(3)
  double precision x2(3)
  double precision x21(3)
  double precision x22(3)
!
!  Choose PHI and THETA counts that make short sides.
!
  phi_num = int ( pi / h )

  if ( h * real ( phi_num) < pi ) then
    phi_num = phi_num + 1
  end if

  theta_num = int ( 2.0D+00 * pi / h )

  if ( h * real ( theta_num) < pi ) then
    theta_num = theta_num + 1
  end if

  n = 0
  result = 0.0D+00
!
!  Only one THETA (and hence, only one PHI.)
!
  if ( theta_num == 1 ) then

    sphere_area = 4.0D+00 * pi

    theta = 0.0D+00
    phi = pi / 2.0D+00
    call tp_to_xyz ( theta, phi, x )
    call f ( 1, x, v )
    result = sphere_area * v(1)
!
!  Several THETA's, one PHI.
!
  else if ( phi_num == 1 ) then

    sphere_area = 4.0D+00 * pi
    sector_area = sphere_area / real ( theta_num)

    result = 0.0D+00

    do j = 1, theta_num

      theta = real ( ( j - 1 ) * 2) * pi &
        / real ( theta_num)
      phi = pi / 2.0D+00
      call tp_to_xyz ( theta, phi, x )
      call f ( 1, x, v )
      n = n + 1
      result = result + sector_area * v(1)

    end do
!
!  At least two PHI's.
!
  else

    result = 0.0D+00
!
!  Picture:
!
!        V1
!       /  \
!      /    \
!    V12----V22
!
    phi1 = 0.0D+00
    phi2 = pi / real ( phi_num)

    do j = 1, theta_num

      theta1 = real ( j - 1) * 2.0D+00 * pi &
             / real ( theta_num)
      theta2 = real ( j) * 2.0D+00 * pi &
             / real ( theta_num)

      call tp_to_xyz ( theta1, phi1, x1 )
      call tp_to_xyz ( theta1, phi2, x12 )
      call tp_to_xyz ( theta2, phi2, x22 )

      call sphere01_triangle_vertices_to_area ( x1, x12, x22, area )

      call f ( 1, x1, v )
      n = n + 1
      result = result + area * v(1) / 3.0D+00
      call f ( 1, x12, v )
      n = n + 1
      result = result + area * v(1) / 3.0D+00
      call f ( 1, x22, v )
      n = n + 1
      result = result + area * v(1) / 3.0D+00

    end do
!
!  Picture:
!
!    V11--V21
!     | \  |
!     |  \ |
!    V12--V22
!
    do i = 2, phi_num-1

      phi1 = real ( i - 1) * pi &
           / real ( phi_num)
      phi2 = real ( i) * pi &
           / real ( phi_num)

      do j = 1, theta_num

        theta1 = real ( j - 1) * 2.0D+00 * pi &
               / real ( theta_num)
        theta2 = real ( j) * 2.0D+00 * pi &
               / real ( theta_num)

        call tp_to_xyz ( theta1, phi1, x11 )
        call tp_to_xyz ( theta2, phi1, x21 )
        call tp_to_xyz ( theta1, phi2, x12 )
        call tp_to_xyz ( theta2, phi2, x22 )

        call sphere01_triangle_vertices_to_area ( x11, x12, x22, area )

        call f ( 1, x11, v )
        n = n + 1
        result = result + area * v(1) / 3.0D+00
        call f ( 1, x12, v )
        n = n + 1
        result = result + area * v(1) / 3.0D+00
        call f ( 1, x22, v )
        n = n + 1
        result = result + area * v(1) / 3.0D+00

        call sphere01_triangle_vertices_to_area ( x22, x21, x11, area )

        call f ( 1, x22, v )
        n = n + 1
        result = result + area * v(1) / 3.0D+00
        call f ( 1, x21, v )
        n = n + 1
        result = result + area * v(1) / 3.0D+00
        call f ( 1, x11, v )
        n = n + 1
        result = result + area * v(1) / 3.0D+00

      end do

    end do
!
!  Picture:
!
!    V11----V21
!      \    /
!       \  /
!        V2
!
    phi1 = real ( phi_num - 1) * pi &
         / real ( phi_num)
    phi2 =                                  pi

    do j = 1, theta_num

      theta1 = real ( j - 1) * 2.0D+00 * pi &
             / real ( theta_num)
      theta2 = real ( j) * 2.0D+00 * pi &
             / real ( theta_num)

      call tp_to_xyz ( theta1, phi1, x11 )
      call tp_to_xyz ( theta2, phi1, x21 )
      call tp_to_xyz ( theta2, phi2, x2 )

      call sphere01_triangle_vertices_to_area ( x11, x2, x21, area )

      call f ( 1, x11, v )
      n = n + 1
      result = result + area * v(1) / 3.0D+00
      call f ( 1, x2, v )
      n = n + 1
      result = result + area * v(1) / 3.0D+00
      call f ( 1, x21, v )
      n = n + 1
      result = result + area * v(1) / 3.0D+00

    end do

  end if
end

subroutine sphere01_quad_mc ( f, h, seed, n, result )

!*****************************************************************************80
!
!! SPHERE01_QUAD_MC uses the Monte Carlo rule for sphere quadrature.
!
!  Discussion:
!
!    A number of points N are chosen at random on the sphere, with N
!    being determined so that, if the points were laid out on a regular
!    grid, the average spacing would be no more than H.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external :: F, evaluates the integrand, of the form:
!      subroutine f ( n, x, v )
!      integer n
!      double precision v(n)
!      double precision x(3,n)
!
!    Input, double precision H, the maximum length of a side of the spherical
!    quadrilaterals.
!
!    Input/output, integer SEED, a seed for the random
!    number generator.
!
!    Input, integer N, the number of points used.
!
!    Output, double precision RESULT, the approximate integral.
!
  implicit none

  integer n

  external f
  double precision h
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision result
  integer seed
  double precision sphere_area
  double precision v(n)
  double precision x(3,n)

  sphere_area = 4.0D+00 * pi

  call sphere01_sample_3d ( n, seed, x )

  call f ( n, x, v )

  result = sphere_area * sum ( v(1:n) ) / real ( n)
end

subroutine sphere01_quad_mc_size ( h, n )

!*****************************************************************************80
!
!! SPHERE01_QUAD_MC_SIZE sizes a Monte Carlo rule for sphere quadrature.
!
!  Discussion:
!
!    A number of points N are chosen at random on the sphere, with N
!    being determined so that, if the points were laid out on a regular
!    grid, the average spacing would be no more than H.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision H, the maximum length of a side of the spherical
!    quadrilaterals.
!
!    Output, integer N, the number of points to use.
!
  implicit none

  double precision h
  integer n
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision sphere_area
!
!  The sphere's area is 4 * PI.
!  Choose N so that we divide this area into N subareas of PI * H * H.
!
  sphere_area = 4.0D+00 * pi

  n = int ( sphere_area / h**2 )
  n = max ( n, 1 )
end

subroutine sphere01_sample_3d ( n, seed, x )

!*****************************************************************************80
!
!! SPHERE01_SAMPLE_3D picks random points on a sphere in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of samples.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, double precision X(3,N), the sample points.
!
  implicit none

  integer n

  integer j
  double precision phi
  double precision , parameter :: pi = 3.141592653589793D+00
  integer seed
  double precision theta
  double precision vdot
  double precision x(3,n)

  do j = 1, n
!
!  Pick a uniformly random VDOT, which must be between -1 and 1.
!  This represents the dot product of the random vector with the Z unit vector.
!
!  Note: this works because the surface area of the sphere between
!  Z and Z + dZ is independent of Z.  So choosing Z uniformly chooses
!  a patch of area uniformly.
!
    call random_number ( harvest = vdot )
    vdot = 2.0D+00 * vdot - 1.0D+00

    phi = acos ( vdot )
!
!  Pick a uniformly random rotation between 0 and 2 Pi around the
!  axis of the Z vector.
!
    call random_number ( harvest = theta )
    theta = 2.0D+00 * pi * theta

    x(1,j) = cos ( theta ) * sin ( phi )
    x(2,j) = sin ( theta ) * sin ( phi )
    x(3,j) =                 cos ( phi )

  end do
end

subroutine sphere01_triangle_angles_to_area ( a, b, c, area )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_ANGLES_TO_AREA computes the area of a spherical triangle.
!
!  Discussion:
!
!    A unit sphere in 3D satisfies the equation:
!
!      X^2 + Y^2 + Z^2 = 1
!
!    A spherical triangle is specified by three points on the surface
!    of the sphere.
!
!    The area formula is known as Girard's formula.
!
!    The area of a spherical triangle is:
!
!      AREA = ( A + B + C - PI )
!
!    where A, B and C are the (surface) angles of the triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision A, B, C, the angles of the triangle.
!
!    Output, double precision AREA, the area of the sphere.
!
  implicit none

  double precision area
  double precision a
  double precision b
  double precision c
  double precision , parameter :: pi = 3.141592653589793D+00
!
!  Apply Girard's formula.
!
  area = a + b + c - pi
end

subroutine sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1, f2, f3, &
  node_xyz )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_PROJECT projects from plane to spherical triangle.
!
!  Discussion:
!
!    We assume that points A, B and C lie on the unit sphere, and they
!    thus define a spherical triangle.
!
!    They also, of course, define a planar triangle.
!
!    Let (F1,F2,F3) be the barycentric coordinates of a point in this 
!    planar triangle.
!
!    This function determines the coordinates of the point in the planar
!    triangle identified by the barycentric coordinates, and returns the
!    coordinates of the projection of that point onto the unit sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision A_XYZ(3), B_XYZ(3), C_XYZ(3), the coordinates
!    of the points A, B, and C.
!
!    Input, integer F1, F2, F3, the barycentric coordinates
!    of a point in the triangle ABC.  Normally, these coordinates would
!    be real numbers, and would sum to 1.  For convenience, we allow these
!    to be integers which must be divided by F1+F2+F3.
!
!    Output, double precision NODE_XYZ(3), the coordinates of the 
!    point on the unit sphere which is the projection of the point on the plane
!    whose barycentric coordinates with respect to A, B, and C is
!    (F1,F2,F3)/(F1+F2+F3).
!
  implicit none

  double precision a_xyz(3)
  double precision b_xyz(3)
  double precision c_xyz(3)
  integer f1
  integer f2
  integer f3
  double precision node_norm
  double precision node_xyz(3)
  double precision r8vec_norm

  node_xyz(1:3) = &
    ( real ( f1) * a_xyz(1:3)   &
    + real (      f2) * b_xyz(1:3)   &
    + real (           f3) * c_xyz(1:3) ) &
    / real ( f1 + f2 + f3)

  node_norm = r8vec_norm ( 3, node_xyz(1:3) )

  node_xyz(1:3) = node_xyz(1:3) / node_norm
end

subroutine sphere01_triangle_project2 ( a_xyz, b_xyz, c_xyz, f1, f2, f3, &
  node_xyz )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_PROJECT2 projects from plane to spherical triangle.
!
!  Discussion:
!
!    We assume that points A, B and C lie on the unit sphere, and they
!    thus define a spherical triangle.
!
!    They also, of course, define a planar triangle.
!
!    Let (F1,F2,F3) be the barycentric coordinates of a point in this 
!    planar triangle.
!
!    This function determines the coordinates of the point in the planar
!    triangle identified by the barycentric coordinates, and returns the
!    coordinates of the projection of that point onto the unit sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision A_XYZ(3), B_XYZ(3), C_XYZ(3), the coordinates
!    of the points A, B, and C.
!
!    Input, integer F1, F2, F3, the barycentric coordinates
!    of a point in the triangle ABC.  Normally, these coordinates would
!    be real numbers, and would sum to 1.  For convenience, we allow these
!    to be integers which must be divided by F1+F2+F3.
!
!    Output, double precision NODE_XYZ(3), the coordinates of the 
!    point on the unit sphere which is the projection of the point on the 
!    plane whose barycentric coordinates with respect to A, B, and C is
!    (F1,F2,F3)/(F1+F2+F3).
!
  implicit none

  double precision a_xyz(3)
  double precision ab(3)
  double precision ac(3)
  double precision acn(3)
  double precision acp(3)
  double precision angle
  double precision b_xyz(3)
  double precision bn(3)
  double precision bp(3)
  double precision c_xyz(3)
  double precision cn(3)
  double precision cp(3)
  integer f1
  integer f2
  integer f3
  double precision node_xyz(3)
  double precision r8vec_norm
  double precision theta_ab
  double precision theta_ac
  double precision theta_bc
!
!  This check avoids 0/0 calculations later.
!
  if ( f2 == 0 .and. f3 == 0 ) then
    node_xyz(1:3) = a_xyz(1:3)
    return
  else if ( f1 == 0 .and. f3 == 0 ) then
    node_xyz(1:3) = b_xyz(1:3)
    return
  else if ( f1 == 0 .and. f2 == 0 ) then
    node_xyz(1:3) = c_xyz(1:3)
  end if
!
!  Determine the angular distances (A,B) and (A,C).
!
  call sphere01_distance_xyz ( a_xyz, b_xyz, theta_ab )

  call sphere01_distance_xyz ( a_xyz, c_xyz, theta_ac )
!
!  Polarize B = BP + BN
!  Normalize BN, 
!  Same for C.
!
  call r8vec_polarize ( 3, b_xyz, a_xyz, bn, bp )
  bn(1:3) = bn(1:3) / r8vec_norm ( 3, bn )

  call r8vec_polarize ( 3, c_xyz, a_xyz, cn, cp )
  cn(1:3) = cn(1:3) / r8vec_norm ( 3, cn )
!
!  Determine AB and AC that use cos ( ( F2 + F3 ) / ( F1 + F2 + F3 ) ) of A
!  and cos ( F1 / ( F1 + F2 + F3 ) ) of B or C.
!
  angle = ( real ( f2 + f3) * theta_ab ) &
    / real ( f1 + f2 + f3)
  ab(1:3) = cos ( angle ) * a_xyz(1:3) + sin ( angle ) * bn(1:3)

  angle = ( real ( f2 + f3) * theta_ac ) &
    / real ( f1 + f2 + f3)
  ac(1:3) = cos ( angle ) * a_xyz(1:3) + sin ( angle ) * cn(1:3)
!
!  Determine the angular distance between AB and AC.
!
  call sphere01_distance_xyz ( ab(1:3), ac(1:3), theta_bc )
!
!  Polarize AC = ACP + ACN, normalize ACN.
!
  call r8vec_polarize ( 3, ac, ab, acn, acp )
  acn(1:3) = acn(1:3) / r8vec_norm ( 3, acn )
!
!  The interval between AB and AC is marked by F2+F3+1 vertices 0 through F2+F3.
!
  angle = ( real ( f3) * theta_bc ) / real ( f2 + f3)

  node_xyz(1:3) = cos ( angle ) * ab(1:3) + sin ( angle ) * acn(1:3)
end

subroutine sphere01_triangle_sample ( n, v1, v2, v3, seed, x )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_SAMPLE: sample points from triangle on unit sphere.
!
!  Discussion:
!
!    The sphere has center 0 and radius 1.
!
!    A spherical triangle on the surface of the unit sphere contains those 
!    points with radius R = 1, bounded by the vertices V1, V2, V3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Arvo,
!    Stratified sampling of spherical triangles,
!    Computer Graphics Proceedings, Annual Conference Series, 
!    ACM SIGGRAPH '95, pages 437-438, 1995.
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, double precision V1(3), V2(3), V3(3), the XYZ coordinates of
!    the vertices of the spherical triangle.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, double precision X(3,N), the XYZ coordinates of the 
!    sample points.
!
  implicit none

  integer n

  double precision a
  double precision alpha
  double precision area
  double precision area_hat
  double precision b
  double precision beta
  double precision c
  double precision gamma
  integer j
  double precision q
  double precision r8_uniform_01
  double precision r8vec_norm
  double precision s
  integer seed
  double precision t
  double precision u
  double precision v
  double precision v1(3)
  double precision v2(3)
  double precision v3(3)
  double precision v31(3)
  double precision v4(3)
  double precision v42(3)
  double precision w
  double precision x(3,n)
  double precision xsi1
  double precision xsi2
  double precision z

  call sphere01_triangle_vertices_to_sides ( v1, v2, v3, a, b, c )

  call sphere01_triangle_sides_to_angles ( a, b, c, alpha, beta, gamma )

  call sphere01_triangle_angles_to_area ( alpha, beta, gamma, area )

  do j = 1, n
!
!  Select the new area.
!
    xsi1 = r8_uniform_01 ( seed )
    area_hat = xsi1 * area
!
!  Compute the sine and cosine of the angle phi.
!
    s = sin ( area_hat - alpha )
    t = cos ( area_hat - alpha )
!
!  Compute the pair that determines beta_hat.
!
    u = t - cos ( alpha )
    v = s + sin ( alpha ) * cos ( c )
!
!  Q is the cosine of the new edge length b_hat.
!
    q = ( ( v * t - u * s ) * cos ( alpha ) - v ) &
      / ( ( v * s + u * t ) * sin ( alpha ) )
!
!  We very occasionally get a Q value out of bounds.
!
    q = max ( q, - 1.0D+00 )
    q = min ( q, + 1.0D+00 )
!
!  V31 = normalized ( V3 - ( V3 dot V1 ) * V1 )
!
    w = dot_product ( v3, v1 )
    v31(1:3) = v3(1:3) - w * v1(1:3)
    v31(1:3) = v31(1:3) / r8vec_norm ( 3, v31(1:3) )
!
!  V4 is the third vertex of the subtriangle V1, V2, V4.
!
    v4(1:3) = q * v1(1:3) + sqrt ( 1.0D+00 - q * q ) * v31(1:3)
!
!  Select cos theta, which will sample along the edge from V2 to V4.
!
    xsi2 = r8_uniform_01 ( seed )
    z = 1.0D+00 - xsi2 * ( 1.0D+00 - dot_product ( v4, v2 ) )
!
!  V42 = normalized ( V4 - ( V4 dot V2 ) * V2 )
!
    w = dot_product ( v4, v2 )
    v42(1:3) = v4(1:3) - w * v2(1:3)
    v42(1:3) = v42(1:3) / r8vec_norm ( 3, v42(1:3) )
!
!  Construct the point.
!
    x(1:3,j) = z * v2(1:3) + sqrt ( 1.0D+00 - z * z ) * v42(1:3)

  end do
end

subroutine sphere01_triangle_sides_to_angles ( as, bs, cs, a, b, c )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_SIDES_TO_ANGLES computes spherical triangle angles.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision AS, BS, CS, the (geodesic) length of the 
!    sides of the triangle.
!
!    Output, double precision A, B, C, the spherical angles of the triangle.
!    Angle A is opposite the side of length AS, and so on.
!
  implicit none

  double precision a
  double precision as
  double precision asu
  double precision b
  double precision bs
  double precision bsu
  double precision c
  double precision cs
  double precision csu
  double precision ssu
  double precision tan_a2
  double precision tan_b2
  double precision tan_c2

  asu = as
  bsu = bs
  csu = cs
  ssu = ( asu + bsu + csu ) / 2.0D+00

  tan_a2 = sqrt ( ( sin ( ssu - bsu ) * sin ( ssu - csu ) ) / &
                  ( sin ( ssu ) * sin ( ssu - asu )     ) )

  a = 2.0D+00 * atan ( tan_a2 )

  tan_b2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - csu ) ) / &
                  ( sin ( ssu ) * sin ( ssu - bsu )     ) )

  b = 2.0D+00 * atan ( tan_b2 )

  tan_c2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - bsu ) ) / &
                  ( sin ( ssu ) * sin ( ssu - csu )     ) )

  c = 2.0D+00 * atan ( tan_c2 )
end

subroutine sphere01_triangle_vertices_to_area ( v1, v2, v3, area )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_VERTICES_TO_AREA computes the area of a spherical triangle.
!
!  Discussion:
!
!    A sphere in 3D satisfies the equation:
!
!      X^2 + Y^2 + Z^2 = 1
!
!    A spherical triangle is specified by three points on the surface
!    of the sphere.
!
!    The area formula is known as Girard's formula.
!
!    The area of a spherical triangle is:
!
!      AREA = ( A + B + C - PI )
!
!    where A, B and C are the (surface) angles of the triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision V1(3), V2(3), V3(3), the vertices of the triangle.
!
!    Output, double precision AREA, the area of the sphere.
!
  implicit none

  double precision area
  double precision a
  double precision as
  double precision b
  double precision bs
  double precision c
  double precision cs
  double precision v1(3)
  double precision v2(3)
  double precision v3(3)
!
!  Compute the lengths of the sides of the spherical triangle.
!
  call sphere01_triangle_vertices_to_sides ( v1, v2, v3, as, bs, cs )
!
!  Get the spherical angles.
!
  call sphere01_triangle_sides_to_angles ( as, bs, cs, a, b, c )
!
!  Get the area.
!
  call sphere01_triangle_angles_to_area ( a, b, c, area )
end

subroutine sphere01_triangle_vertices_to_centroid ( v1, v2, v3, vs )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_VERTICES_TO_CENTROID gets a spherical triangle "centroid".
!
!  Discussion:
!
!    A sphere in 3D satisfies the equation:
!
!      X^2 + Y^2 + Z^2 = 1
!
!    A spherical triangle is specified by three points on the sphere.
!
!    The (true) centroid of a spherical triangle is the point
!
!      VT = (XT,YT,ZT) = Integral ( X, Y, Z ) dArea / Integral 1 dArea
!
!    Note that the true centroid does NOT, in general, lie on the sphere.
!
!    The "flat" centroid VF is the centroid of the planar triangle defined by
!    the vertices of the spherical triangle.
!
!    The "spherical" centroid VS of a spherical triangle is computed by
!    the intersection of the geodesic bisectors of the triangle angles.
!    The spherical centroid lies on the sphere.
!
!    VF, VT and VS lie on a line through the center of the sphere.  We can
!    easily calculate VF by averaging the vertices, and from this determine
!    VS by normalizing.
!
!    (Of course, we still will not have actually computed VT, which lies
!    somewhere between VF and VS!)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision V1(3), V2(3), V3(3), the vertices of the triangle.
!
!    Output, double precision VS(3), the coordinates of the "spherical
!    centroid" of the spherical triangle.
!
  implicit none

  double precision norm
  double precision v1(3)
  double precision v2(3)
  double precision v3(3)
  double precision vs(3)

  vs(1:3) = ( v1(1:3) + v2(1:3) + v3(1:3) ) / 3.0D+00

  norm = sqrt ( sum ( vs(1:3)**2 ) )

  vs(1:3) = vs(1:3) / norm
end

subroutine sphere01_triangle_vertices_to_midpoints ( v1, v2, v3, m1, m2, m3 )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_VERTICES_TO_MIDPOINTS: midsides of a spherical triangle.
!
!  Discussion:
!
!    The points are assumed to lie on the unit sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision V1(3), V2(3), V3(3), the vertices of the triangle.
!
!    Output, double precision M1(3), M2(3), M3(3), the coordinates of 
!    the midpoints of the sides of the spherical triangle.
!
  implicit none

  double precision m1(3)
  double precision m2(3)
  double precision m3(3)
  double precision norm
  double precision v1(3)
  double precision v2(3)
  double precision v3(3)

  m1(1:3) = ( v1(1:3) + v2(1:3) ) / 2.0D+00
  norm = sqrt ( sum ( m1(1:3)**2 ) )
  m1(1:3) = m1(1:3) / norm

  m2(1:3) = ( v2(1:3) + v3(1:3) ) / 2.0D+00
  norm = sqrt ( sum ( m2(1:3)**2 ) )
  m2(1:3) = m2(1:3) / norm

  m3(1:3) = ( v3(1:3) + v1(1:3) ) / 2.0D+00
  norm = sqrt ( sum ( m3(1:3)**2 ) )
  m3(1:3) = m3(1:3) / norm
end

subroutine sphere01_triangle_vertices_to_sides ( v1, v2, v3, as, bs, cs )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_VERTICES_TO_SIDES computes spherical triangle sides.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision V1(3), V2(3), V3(3), the vertices of the spherical
!    triangle.
!
!    Output, double precision AS, BS, CS, the (geodesic) length of the 
!    sides of the triangle.
!
  implicit none

  double precision as
  double precision bs
  double precision cs
  double precision v1(3)
  double precision v2(3)
  double precision v3(3)

  as = acos ( dot_product ( v2(1:3), v3(1:3) ) )
  bs = acos ( dot_product ( v3(1:3), v1(1:3) ) )
  cs = acos ( dot_product ( v1(1:3), v2(1:3) ) )
end

subroutine tp_to_xyz ( theta, phi, v )

!*****************************************************************************80
!
!! TP_TO_XYZ converts spherical TP coordinates to XYZ coordinates.
!
!  Discussion:
!
!    The sphere is assumed to have radius 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision THETA, PHI, the spherical coordinates 
!    of a point.
!
!    Output, double precision V(3), the XYZ coordinates.
!
  implicit none

  double precision phi
  double precision theta
  double precision v(3)

  v(1) = cos ( theta ) * sin ( phi )
  v(2) = sin ( theta ) * sin ( phi )
  v(3) =                 cos ( phi )
end
