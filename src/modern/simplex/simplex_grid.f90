!> simplex_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

subroutine comp_next_grlex ( kc, xc )

!*****************************************************************************80
!
!! COMP_NEXT_GRLEX returns the next composition in grlex order.
!
!  Discussion:
!
!    Example:
!
!    KC = 3
!
!    #   XC(1  XC(2) XC(3)  Degree
!      +------------------------
!    1 |  0     0     0        0
!      |
!    2 |  0     0     1        1
!    3 |  0     1     0        1
!    4 |  1     0     0        1
!      |
!    5 |  0     0     2        2
!    6 |  0     1     1        2
!    7 |  0     2     0        2
!    8 |  1     0     1        2
!    9 |  1     1     0        2
!   10 |  2     0     0        2
!      |
!   11 |  0     0     3        3
!   12 |  0     1     2        3
!   13 |  0     2     1        3
!   14 |  0     3     0        3
!   15 |  1     0     2        3
!   16 |  1     1     1        3
!   17 |  1     2     0        3
!   18 |  2     0     1        3
!   19 |  2     1     0        3
!   20 |  3     0     0        3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer KC, the number of parts of the composition.
!    1 <= KC.
!
!    Input/output, integer XC(KC), the current composition.
!    Each entry of XC must be nonnegative.
!    On return, XC has been replaced by the next composition in the
!    grlex order.
!
  implicit none

  integer kc

  integer i
  integer im1
  integer j
  integer t
  integer xc(kc)
!
!  Ensure that 1 <= KC.
!
  if ( kc < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COMP_NEXT_GRLEX - Fatal error!'
    write ( *, '(a)' ) '  KC < 1'
    stop 1
  end if
!
!  Ensure that 0 <= XC(I).
!
  do i = 1, kc
    if ( xc(i) < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COMP_NEXT_GRLEX - Fatal error!'
      write ( *, '(a)' ) '  XC(I) < 0'
      stop 1
    end if
  end do
!
!  Find I, the index of the rightmost nonzero entry of X.
!
  i = 0
  do j = kc, 1, -1
    if ( 0 < xc(j) ) then
      i = j
      exit
    end if
  end do
!
!  set T = X(I)
!  set XC(I) to zero,
!  increase XC(I-1) by 1,
!  increment XC(KC) by T-1.
!
  if ( i == 0 ) then
    xc(kc) = 1
    return
  else if ( i == 1 ) then
    t = xc(1) + 1
    im1 = kc
  else if ( 1 < i ) then
    t = xc(i)
    im1 = i - 1
  end if

  xc(i) = 0
  xc(im1) = xc(im1) + 1
  xc(kc) = xc(kc) + t - 1
end

subroutine comp_random ( n, k, seed, a )

!*****************************************************************************80
!
!! COMP_RANDOM selects a random composition of the integer N into K parts.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2003
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Second Edition,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer N, the integer to be decomposed.
!
!    Input, integer K, the number of parts in the composition.
!
!    Input/output, integer SEED, a seed for the random number
!    generator.
!
!    Output, integer A(K), the parts of the composition.
!
  implicit none

  integer k

  integer a(k)
  integer i
  integer l
  integer m
  integer n
  integer seed

  call ksub_random ( n + k - 1, k - 1, seed, a )

  a(k) = n + k
  l = 0

  do i = 1, k
    m = a(i)
    a(i) = a(i) - l - 1
    l = m
  end do
end

function i4_uniform_ab ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
!
!  Discussion:
!
!    An I4 is an integer value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 October 2012
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
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer A, B, the limits of the interval.
!
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer I4_UNIFORM_AB, a number between A and B.
!
  implicit none

  integer a
  integer b
  integer , parameter :: i4_huge = 2147483647
  integer i4_uniform_ab
  integer k
  real r
  integer seed
  integer value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM_AB - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r = real ( seed) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b )) - 0.5E+00 ) & 
    +             r   * ( real ( max ( a, b )) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r)

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform_ab = value
end

subroutine ksub_random ( n, k, seed, a )

!*****************************************************************************80
!
!! KSUB_RANDOM selects a random subset of size K from a set of size N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2003
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Second Edition,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer N, the size of the set from which subsets
!    are drawn.
!
!    Input, integer K, number of elements in desired subsets.
!    K must be between 0 and N.
!
!    Input/output, integer SEED, a seed for the random number
!    generator.
!
!    Output, integer A(K).  A(I) is the I-th element of the
!    output set.  The elements of A are in order.
!
  implicit none

  integer k

  integer a(k)
  integer i
  integer i4_uniform_ab
  integer ids
  integer ihi
  integer ip
  integer ir
  integer is
  integer ix
  integer l
  integer ll
  integer m
  integer m0
  integer n
  integer seed

  if ( k < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  K = ', k
    write ( *, '(a)' ) '  but 0 <= K is required!'
    stop
  else if ( n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a,i8)' ) '  K = ', k
    write ( *, '(a)' ) '  K <= N is required!'
    stop
  end if

  if ( k == 0 ) then
  end if

  do i = 1, k
    a(i) = ( ( i - 1 ) * n ) / k
  end do

  do i = 1, k

    do

      ix = i4_uniform_ab ( 1, n, seed )

      l = 1 + ( ix * k - 1 ) / n

      if ( a(l) < ix ) then
        exit
      end if

    end do

    a(l) = a(l) + 1

  end do

  ip = 0
  is = k

  do i = 1, k

    m = a(i)
    a(i) = 0

    if ( m /= ( ( i - 1 ) * n ) / k ) then
      ip = ip + 1
      a(ip) = m
    end if

  end do

  ihi = ip

  do i = 1, ihi
    ip = ihi + 1 - i
    l = 1 + ( a(ip) * k - 1 ) / n
    ids = a(ip) - ( ( l - 1 ) * n ) / k
    a(ip) = 0
    a(is) = l
    is = is - ids
  end do

  do ll = 1, k

    l = k + 1 - ll

    if ( a(l) /= 0 ) then
      ir = l
      m0 = 1 + ( ( a(l) - 1 ) * n ) / k
      m = ( a(l) * n ) / k - m0 + 1
    end if

    ix = i4_uniform_ab ( m0, m0 + m - 1, seed )

    i = l + 1

    do while ( i <= ir )

      if ( ix < a(i) ) then
        exit
      end if

      ix = ix + 1
      a(i-1) = a(i)
      i = i + 1

    end do

    a(i-1) = ix
    m = m - 1

  end do
end

subroutine simplex_grid_index_all ( m, n, ng, grid )

!*****************************************************************************80
!
!! SIMPLEX_GRID_INDEX_ALL returns all the simplex grid indices.
!
!  Discussion:
!
!    The number of grid indices can be determined by calling 
!      ng = simplex_grid_size ( m, n )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of subintervals.
!
!    Input, integer NG, the number of values in the grid.
!
!    Output, integer GRID(M+1,NG), the current, and then the next,
!    grid index.
!
  implicit none

  integer m
  integer ng

  integer g(m+1)
  integer grid(m+1,ng)
  integer i
  integer k
  integer n

  do i = 1, m
    g(i) = 0
  end do
  g(m+1) = n

  k = 1
  grid(1:m+1,k) = g(1:m+1)

  do while ( k < ng )
    call comp_next_grlex ( m + 1, g )
    k = k + 1
    grid(1:m+1,k) = g(1:m+1)
  end do
end

subroutine simplex_grid_index_next ( m, n, g )

!*****************************************************************************80
!
!! SIMPLEX_GRID_INDEX_NEXT returns the next simplex grid index.
!
!  Discussion:
!
!    The vector G has dimension M+1.  The first M entries may be regarded
!    as grid coordinates.  These coordinates must have a sum between 0 and N.
!    The M+1 entry contains the remainder, that is N minus the sum of the
!    first M coordinates.
!
!    Each time the function is called, it is given a current grid index, and
!    computes the next one.  The very first index is all zero except for a 
!    final value of N, and the very last index has all zero except for an'
!    intial value of N.
!
!    For example, here are the coordinates in order for M = 3, N = 3:
!
!     0  0 0 0 3
!     1  0 0 1 2
!     2  0 0 2 1
!     3  0 0 3 0
!     4  0 1 0 2
!     5  0 1 1 1
!     6  0 1 2 0
!     7  0 2 0 1
!     8  0 2 1 0
!     9  0 3 0 0
!    10  1 0 0 2
!    11  1 0 1 1
!    12  1 0 2 0
!    13  1 1 0 1
!    14  1 1 1 0
!    15  1 2 0 0
!    16  2 0 0 1
!    17  2 0 1 0
!    18  2 1 0 0
!    19  3 0 0 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of subintervals.
!
!    Input/output, integer G(M+1), the current, and then the next,
!    grid index.
!
  implicit none

  integer m

  integer g(m+1)
  integer n

  call comp_next_grlex ( m + 1, g )
end

subroutine simplex_grid_index_sample ( m, n, seed, g )

!*****************************************************************************80
!
!! SIMPLEX_GRID_INDEX_SAMPLE returns a random simplex grid index.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of subintervals in
!    each dimension.
!
!    Input, integer SEED, a seed for the random number generator.
!
!    Output, integer G(M+1), a randomly selected index in the 
!    simplex grid.
!
!    Output, integer SEED, the updated random number seed.
!
  implicit none

  integer m

  integer g(m+1)
  integer n
  integer seed

  call comp_random ( n, m + 1, seed, g )
end

subroutine simplex_grid_index_to_point ( m, n, ng, g, v, x )

!*****************************************************************************80
!
!! SIMPLEX_GRID_INDEX_TO_POINT returns  points corresponding to simplex indices.
!
!  Discussion:
!
!    The M-dimensional simplex is defined by M+1 vertices.
!
!    Given a regular grid that uses N subintervals along the edge between
!    each pair of vertices, a simplex grid index G is a set of M+1 values
!    each between 0 and N, and summing to N. 
!
!    This function determines the coordinates X of the point corresponding
!    to the index G.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of subintervals.
!
!    Input, integer NG, the number of grid indices to be converted.
!
!    Input, integer G(M+1,NG), the grid indices of 1 
!    or more points.
!
!    Input, double precision V(M,M+1), the coordinates of the vertices 
!    of the simplex.
!
!    Output, double precision X(M,NG), the coordinates of one or more points.
!
  implicit none

  integer m
  integer n
  integer ng

  integer g(m+1,ng)
  integer i
  integer j
  integer k
  double precision v(m,m+1)
  double precision x(m,ng)
  
  do j = 1, ng
    do i = 1, m
      x(i,j) = 0.0D+00
      do k = 1, m + 1
        x(i,j) = x(i,j) + v(i,k) * real ( g(k,j))
      end do
      x(i,j) = x(i,j) / real ( n)
    end do
  end do
end

subroutine simplex_grid_size ( m, n, ng )

!*****************************************************************************80
!
!! SIMPLEX_GRID_SIZE counts the grid points inside a simplex.
!
!  Discussion:
!
!    The size of a grid with parameters M, N is C(M+N,N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of subintervals.
!
!    Output, integer NG, the number of grid points.
!
  implicit none

  integer i
  integer m
  integer n
  integer ng

  ng = 1

  do i = 1, m
    ng = ( ng * ( n + i ) ) / i
  end do
end
