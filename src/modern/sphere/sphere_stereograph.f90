!> sphere_stereograph — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

subroutine plane_normal_basis_3d ( pp, normal, pq, pr )

!*****************************************************************************80
!
!! PLANE_NORMAL_BASIS_3D finds two perpendicular vectors in a plane in 3D.
!
!  Discussion:
!
!    The normal form of a plane in 3D is:
!
!      PP is a point on the plane,
!      N is a normal vector to the plane.
!
!    The two vectors to be computed, PQ and PR, can be regarded as
!    the basis of a Cartesian coordinate system for points in the plane.
!    Any point in the plane can be described in terms of the "origin" 
!    point PP plus a weighted sum of the two vectors PQ and PR:
!
!      P = PP + a * PQ + b * PR.
!
!    The vectors PQ and PR have unit length, and are perpendicular to N
!    and to each other.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision PP(3), a point on the plane.  (Actually,
!    we never need to know these values to do the calculation!)
!
!    Input, double precision NORMAL(3), a normal vector N to the plane.  The
!    vector must not have zero length, but it is not necessary for N
!    to have unit length.
!
!    Output, double precision PQ(3), a vector of unit length,
!    perpendicular to the vector N and the vector PR.
!
!    Output, double precision PR(3), a vector of unit length,
!    perpendicular to the vector N and the vector PQ.
!
  implicit none

  integer , parameter :: dim_num = 3

  double precision r8vec_norm
  double precision normal(dim_num)
  double precision normal_norm
  double precision pp(dim_num)
  double precision pq(dim_num)
  double precision pr(dim_num)
  double precision pr_norm
!
!  Compute the length of NORMAL.
!
  normal_norm = r8vec_norm ( dim_num, normal )

  if ( normal_norm == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_NORMAL_BASIS_3D - Fatal error!'
    write ( *, '(a)' ) '  The normal vector is 0.'
    stop
  end if
!
!  Find a vector PQ that is normal to NORMAL and has unit length.
!
  call r8vec_any_normal ( dim_num, normal, pq )
!
!  Now just take the cross product NORMAL x PQ to get the PR vector.
!
  call r8vec_cross_product_3d ( normal, pq, pr )

  pr_norm = r8vec_norm ( dim_num, pr )

  pr(1:dim_num) = pr(1:dim_num) / pr_norm
end

function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
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
!    31 May 2007
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
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer k
  double precision r8_uniform_01
  integer seed

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed) * 4.656612875D-10
end

function r8mat_norm_fro_affine ( m, n, a1, a2 )

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
!    Input, integer M, the number of rows.
!
!    Input, integer N, the number of columns.
!
!    Input, double precision A1(M,N), A2(M,N), the matrices for whose 
!    difference the Frobenius norm is desired.
!
!    Output, double precision R8MAT_NORM_FRO_AFFINE, the Frobenius 
!    norm of A1 - A2.
!
  implicit none

  integer m
  integer n

  double precision a1(m,n)
  double precision a2(m,n)
  double precision r8mat_norm_fro_affine

  r8mat_norm_fro_affine = sqrt ( sum ( ( a1(1:m,1:n) - a2(1:m,1:n) )**2 ) )
end

subroutine r8mat_normal_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_NORMAL_01 returns a unit pseudonormal R8MAT.
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
!    Input, integer M, N, the number of rows and columns
!    in the array.
!
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, double precision R(M,N), the array of pseudonormal values.
!
  implicit none

  integer m
  integer n

  integer seed
  double precision r(m,n)

  call r8vec_normal_01 ( m * n, seed, r )
end

subroutine r8mat_uniform_01 ( m, n, seed, r )

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
!    Input, integer M, N, the number of rows and columns in
!    the array.
!
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, double precision R(M,N), the array of pseudorandom values.
!
  implicit none

  integer m
  integer n

  integer i
  integer , parameter :: i4_huge = 2147483647
  integer j
  integer k
  integer seed
  double precision r(m,n)

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r(i,j) = real ( seed) * 4.656612875D-10

    end do
  end do
end

subroutine r8vec_any_normal ( dim_num, v1, v2 )

!*****************************************************************************80
!
!! R8VEC_ANY_NORMAL returns some normal vector to V1.
!
!  Discussion:
!
!    If DIM_NUM < 2, then no normal vector can be returned.
!
!    If V1 is the zero vector, then any unit vector will do.
!
!    No doubt, there are better, more robust algorithms.  But I will take
!    just about ANY reasonable unit vector that is normal to V1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, double precision V1(DIM_NUM), the vector.
!
!    Output, double precision V2(DIM_NUM), a vector that is
!    normal to V2, and has unit Euclidean length.
!
  implicit none

  integer dim_num

  integer i
  integer j
  integer k
  double precision r8vec_norm
  double precision v1(dim_num)
  double precision v2(dim_num)
  double precision vj
  double precision vk

  if ( dim_num < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_ANY_NORMAL - Fatal error!'
    write ( *, '(a)' ) '  Called with DIM_NUM < 2.'
    stop
  end if

  if ( r8vec_norm ( dim_num, v1 ) == 0.0D+00 ) then
    v2(1) = 1.0D+00
    v2(2:dim_num) = 0.0D+00
  end if
!
!  Seek the largest entry in V1, VJ = V1(J), and the
!  second largest, VK = V1(K).
!
!  Since V1 does not have zero norm, we are guaranteed that
!  VJ, at least, is not zero.
!
  j = - 1
  vj = 0.0D+00

  k = - 1
  vk = 0.0D+00

  do i = 1, dim_num

    if ( abs ( vk ) < abs ( v1(i) ) .or. k < 1 ) then

      if ( abs ( vj ) < abs ( v1(i) ) .or. j < 1 ) then
        k = j
        vk = vj
        j = i
        vj = v1(i)
      else
        k = i
        vk = v1(i)
      end if

    end if

  end do
!
!  Setting V2 to zero, except that V2(J) = -VK, and V2(K) = VJ,
!  will just about do the trick.
!
  v2(1:dim_num) = 0.0D+00

  v2(j) = - vk / sqrt ( vk * vk + vj * vj )
  v2(k) =   vj / sqrt ( vk * vk + vj * vj )
end

subroutine r8vec_cross_product_3d ( v1, v2, v3 )

!*****************************************************************************80
!
!! R8VEC_CROSS_PRODUCT_3D computes the cross product of two R8VEC's in 3D.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The cross product in 3D can be regarded as the determinant of the
!    symbolic matrix:
!
!          |  i  j  k |
!      det | x1 y1 z1 |
!          | x2 y2 z2 |
!
!      = ( y1 * z2 - z1 * y2 ) * i
!      + ( z1 * x2 - x1 * z2 ) * j
!      + ( x1 * y2 - y1 * x2 ) * k
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision V1(3), V2(3), the two vectors.
!
!    Output, double precision V3(3), the cross product vector.
!
  implicit none

  double precision v1(3)
  double precision v2(3)
  double precision v3(3)

  v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
  v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
  v3(3) = v1(1) * v2(2) - v1(2) * v2(1)
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

function r8vec_norm_affine ( n, v0, v1 )

!*****************************************************************************80
!
!! R8VEC_NORM_AFFINE returns the affine norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The affine vector L2 norm is defined as:
!
!      R8VEC_NORM_AFFINE(V0,V1)
!        = sqrt ( sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the vectors.
!
!    Input, double precision V0(N), the base vector.
!
!    Input, double precision V1(N), the vector whose affine norm is desired.
!
!    Output, double precision R8VEC_NORM_AFFINE, the L2 norm of V1-V0.
!
  implicit none

  integer n

  double precision r8vec_norm_affine
  double precision v0(n)
  double precision v1(n)

  r8vec_norm_affine = sqrt ( sum ( ( v0(1:n) - v1(1:n) )**2 ) )
end

subroutine r8vec_normal_01 ( n, seed, x )

!*****************************************************************************80
!
!! R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
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
!    The Box-Muller method is used, which is efficient, but
!    generates an even number of values each time.  On any call
!    to this routine, an even number of new values are generated.
!    Depending on the situation, one value may be left over.
!    In that case, it is saved for the next call.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of values desired.  If N is
!    negative, then the code will flush its internal memory; in particular,
!    if there is a saved value to be used on the next call, it is
!    instead discarded.  This is useful if the user has reset the
!    random number seed, for instance.
!
!    Input/output, integer SEED, a seed for the random
!    number generator.
!
!    Output, double precision X(N), a sample of the standard normal PDF.
!
!  Local parameters:
!
!    Local, integer MADE, records the number of values that have
!    been computed.  On input with negative N, this value overwrites
!    the return value of N, so the user can get an accounting of
!    how much work has been done.
!
!    Local, double precision R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer SAVED, is 0 or 1 depending on whether there
!    is a single saved value left over from the previous call.
!
!    Local, integer X_LO_INDEX, X_HI_INDEX, records the range
!    of entries of X that we need to compute.  This starts off as 1:N, but
!    is adjusted if we have a saved value that can be immediately stored
!    in X(1), and so on.
!
!    Local, double precision Y, the value saved from the previous call, if
!    SAVED is 1.
!
  implicit none

  integer n

  integer m
  integer , save :: made = 0
  double precision , parameter :: pi = 3.141592653589793D+00
  double precision r(n+1)
  double precision r8_uniform_01
  integer , save :: saved = 0
  integer seed
  double precision x(n)
  integer x_hi_index
  integer x_lo_index
  double precision , save :: y = 0.0D+00
!
!  I'd like to allow the user to reset the internal data.
!  But this won't work properly if we have a saved value Y.
!  I'm making a crock option that allows the user to signal
!  explicitly that any internal memory should be flushed,
!  by passing in a negative value for N.
!
  if ( n < 0 ) then
    n = made
    made = 0
    saved = 0
    y = 0.0D+00
    return
  else if ( n == 0 ) then
  end if
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  Use up the old value, if we have it.
!
  if ( saved == 1 ) then
    x(1) = y
    saved = 0
    x_lo_index = 2
  end if
!
!  Maybe we don't need any more values.
!
  if ( x_hi_index - x_lo_index + 1 == 0 ) then
!
!  If we need just one new value, do that here to avoid null arrays.
!
  else if ( x_hi_index - x_lo_index + 1 == 1 ) then

    r(1) = r8_uniform_01 ( seed )

    if ( r(1) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_NORMAL_01 - Fatal error!'
      write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
      stop
    end if

    r(2) = r8_uniform_01 ( seed )

    x(x_hi_index) = &
             sqrt ( - 2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * pi * r(2) )
    y =      sqrt ( - 2.0D+00 * log ( r(1) ) ) * sin ( 2.0D+00 * pi * r(2) )

    saved = 1

    made = made + 2
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index, 2 ) == 1 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m:2) )

    made = made + x_hi_index - x_lo_index + 1
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
      sqrt ( - 2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(n) = sqrt ( - 2.0D+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * pi * r(2*m) )

    y = sqrt ( - 2.0D+00 * log ( r(2*m-1) ) ) &
      * sin ( 2.0D+00 * pi * r(2*m) )

    saved = 1

    made = made + x_hi_index - x_lo_index + 2

  end if
end

subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
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
!    Input, integer N, the number of entries in the vector.
!
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, double precision R(N), the vector of pseudorandom values.
!
  implicit none

  integer n

  integer i
  integer k
  integer seed
  double precision r(n)

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed) * 4.656612875D-10

  end do
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

subroutine sphere_stereograph2 ( m, n, p, focus, center, q )

!*****************************************************************************80
!
!! SPHERE_STEREOGRAPH2 computes the stereographic image of points on a sphere.
!
!  Discussion:
!
!    We start with a sphere of center C.
!
!    F is a point on the sphere which is the focus of the mapping,
!    and the antipodal point 2*C-F is the point of tangency
!    to the sphere of a plane.
!
!    For any point P on the sphere, the stereographic projection Q of the
!    point is defined by drawing the line from F through P, and computing
!    Q as the intersection of this line with the plane.
!
!    The spatial dimension M is arbitrary, but should be at least 2.
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
!    Input, integer the number of points.
!
!    Input, double precision P[M*N], a set of points on the unit sphere.
!
!    Input, double precision FOCUS[M], the coordinates of the focus point.
!
!    Input, double precision CENTER[M], the coordinates of the center of 
!    the sphere.
!
!    Output, double precision Q[M*N], the coordinates of the
!    image points,
!
  implicit none

  integer m
  integer n

  double precision center(m)
  double precision cf_dot_pf
  double precision cf_normsq
  double precision focus(m)
  integer i
  integer j
  double precision p(m,n)
  double precision q(m,n)
  double precision s

  do j = 1, n 
    cf_normsq = 0.0D+00
    cf_dot_pf = 0.0D+00
    do i = 1, m
      cf_normsq = cf_normsq + ( center(i) - focus(i) ) ** 2
      cf_dot_pf = cf_dot_pf + ( center(i) - focus(i) ) * ( p(i,j) - focus(i) )
    end do
    s = 2.0D+00 * cf_normsq / cf_dot_pf
    do i = 1, m
      q(i,j) = s * p(i,j) + ( 1.0D+00 - s ) * focus(i)
    end do
  end do
end

subroutine sphere_stereograph2_inverse ( m, n, q, focus, center, p )

!*****************************************************************************80
!
!! SPHERE_STEREOGRAPH2_INVERSE computes stereographic preimages of points.
!
!  Discussion:
!
!    We start with a sphere of center C.
!
!    F is a point on the sphere which is the focus of the mapping,
!    and the antipodal point 2*C-F is the point of tangency
!    to the sphere of a plane.
!
!    For any point Q on the plane, the stereographic inverse projection
!    P of the point is defined by drawing the line from F through Q, and
!    computing P as the intersection of this line with the sphere.
!
!    The spatial dimension M is arbitrary, but should be at least 2.
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
!    on the plane.
!
!    Input, double precision FOCUS(M), the coordinates of the focus point.
!
!    Input, double precision CENTER(M), the coordinates of the center 
!    of the sphere.
!
!    Output, double precision P(M,N), the stereographic
!    inverse projections of the points.
!
  implicit none

  integer m
  integer n

  double precision center(m)
  double precision cf_dot_qf
  double precision focus(m)
  integer i
  integer j
  double precision p(m,n)
  double precision q(m,n)
  double precision qf_normsq
  double precision s

  do j = 1, n

    cf_dot_qf = 0.0D+00
    qf_normsq = 0.0D+00
    do i = 1, m
      cf_dot_qf = cf_dot_qf + ( center(i) - focus(i) ) * ( q(i,j) - focus(i) )
      qf_normsq = qf_normsq + ( q(i,j) - focus(i) ) ** 2
    end do

    s = 2.0D+00 * cf_dot_qf / qf_normsq
    do i = 1, m
      p(i,j) = s * q(i,j) + ( 1.0D+00 - s ) * focus(i)
    end do

  end do
end

subroutine uniform_on_sphere01_map ( dim_num, n, seed, x )

!*****************************************************************************80
!
!! UNIFORM_ON_SPHERE01_MAP maps uniform points onto the unit sphere.
!
!  Discussion:
!
!    The sphere has center 0 and radius 1.
!
!    This procedure is valid for any spatial dimension DIM_NUM.
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
!    Reuven Rubinstein,
!    Monte Carlo Optimization, Simulation, and Sensitivity
!    of Queueing Networks,
!    Krieger, 1992,
!    ISBN: 0894647644,
!    LC: QA298.R79.
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Input, integer N, the number of points.
!
!    Input/output, integer SEED, a seed for the random
!    number generator.
!
!    Output, double precision X(DIM_NUM,N), the points.
!
  implicit none

  integer dim_num
  integer n

  integer j
  double precision norm
  integer seed
  double precision x(dim_num,n)
!
!  Fill a matrix with normally distributed values.
!
  call r8mat_normal_01 ( dim_num, n, seed, x )
!
!  Normalize each column.
!
  do j = 1, n
!
!  Compute the length of the vector.
!
    norm = sqrt ( sum ( x(1:dim_num,j)**2 ) )
!
!  Normalize the vector.
!
    x(1:dim_num,j) = x(1:dim_num,j) / norm

  end do
end
