!> cube_grid — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

subroutine cube_grid ( n, ns, a, b, c, x )

!*****************************************************************************80
!
!! CUBE_GRID: grid points over the interior of a cube in 3D.
!
!  Discussion:
!
!    In 3D, a logically rectangular grid is to be created.
!    In the I-th dimension, the grid will use S(I) points.
!    The total number of grid points is 
!      N = product ( 1 <= I <= 3 ) S(I)
!
!    Over the interval [A(i),B(i)], we have 5 choices for grid centering:
!      1: 0,   1/3, 2/3, 1
!      2: 1/5, 2/5, 3/5, 4/5
!      3: 0,   1/4, 2/4, 3/4
!      4: 1/4, 2/4, 3/4, 1
!      5: 1/8, 3/8, 5/8, 7/8
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 August 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!    N = product ( 1 <= I <= 3 ) NS(I).
!
!    Input, integer NS(3), the number of points along 
!    each dimension.
!
!    Input, double precision A(3), B(3), the endpoints for each dimension.
!
!    Input, integer C(3), the grid centering for each dimension.
!    1 <= C(*) <= 5.
!
!    Output, double precision X(3,N) = X(3,S(1)*S(2)*S(3)), the points.
!
  implicit none

  integer , parameter :: m = 3
  integer n

  double precision a(m)
  double precision b(m)
  integer c(m)
  integer i
  integer j
  integer ns(m)
  integer s
  double precision x(m,n)
  double precision , allocatable :: xs(:)
!
!  Create the 1D grids in each dimension.
!
  do i = 1, m

    s = ns(i)

    allocate ( xs(1:s) )

    do j = 1, s

      if ( c(i) == 1 ) then

        if ( s == 1 ) then
          xs(j) = 0.5D+00 * ( a(i) + b(i) )
        else
          xs(j) = (   real ( s - j) * a(i)   &
                    + real (     j - 1) * b(i) ) & 
                    / real ( s     - 1)
        end if
      else if ( c(i) == 2 ) then
        xs(j) = (   real ( s - j + 1) * a(i)   &
                  + real (     j) * b(i) ) & 
                  / real ( s     + 1)
      else if ( c(i) == 3 ) then
        xs(j) = (   real ( s - j + 1) * a(i)   &
                  + real (     j - 1) * b(i) ) & 
                  / real ( s)
      else if ( c(i) == 4 ) then
        xs(j) = (   real ( s - j) * a(i)   &
                  + real (     j) * b(i) ) & 
                  / real ( s)
      else if ( c(i) == 5 ) then
        xs(j) = (   real ( 2 * s - 2 * j + 1) * a(i)   &
                  + real (         2 * j - 1) * b(i) ) & 
                  / real ( 2 * s)
      end if

    end do

    call r8vec_direct_product ( i, s, xs, m, n, x )

    deallocate ( xs )

  end do
end

subroutine r8vec_direct_product ( factor_index, factor_order, factor_value, &
  factor_num, point_num, x )

!*****************************************************************************80
!
!! R8VEC_DIRECT_PRODUCT creates a direct product of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    To explain what is going on here, suppose we had to construct
!    a multidimensional quadrature rule as the product of K rules
!    for 1D quadrature.
!
!    The product rule will be represented as a list of points and weights.
!
!    The J-th item in the product rule will be associated with
!      item J1 of 1D rule 1,
!      item J2 of 1D rule 2,
!      ...,
!      item JK of 1D rule K.
!
!    In particular,
!      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
!    and
!      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
!
!    So we can construct the quadrature rule if we can properly
!    distribute the information in the 1D quadrature rules.
!
!    This routine carries out that task for the abscissas X.
!
!    Another way to do this would be to compute, one by one, the
!    set of all possible indices (J1,J2,...,JK), and then index
!    the appropriate information.  An advantage of the method shown
!    here is that you can process the K-th set of information and
!    then discard it.
!
!  Example:
!
!    Rule 1:
!      Order = 4
!      X(1:4) = ( 1, 2, 3, 4 )
!
!    Rule 2:
!      Order = 3
!      X(1:3) = ( 10, 20, 30 )
!
!    Rule 3:
!      Order = 2
!      X(1:2) = ( 100, 200 )
!
!    Product Rule:
!      Order = 24
!      X(1:24) =
!        ( 1, 10, 100 )
!        ( 2, 10, 100 )
!        ( 3, 10, 100 )
!        ( 4, 10, 100 )
!        ( 1, 20, 100 )
!        ( 2, 20, 100 )
!        ( 3, 20, 100 )
!        ( 4, 20, 100 )
!        ( 1, 30, 100 )
!        ( 2, 30, 100 )
!        ( 3, 30, 100 )
!        ( 4, 30, 100 )
!        ( 1, 10, 200 )
!        ( 2, 10, 200 )
!        ( 3, 10, 200 )
!        ( 4, 10, 200 )
!        ( 1, 20, 200 )
!        ( 2, 20, 200 )
!        ( 3, 20, 200 )
!        ( 4, 20, 200 )
!        ( 1, 30, 200 )
!        ( 2, 30, 200 )
!        ( 3, 30, 200 )
!        ( 4, 30, 200 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer FACTOR_INDEX, the index of the factor being
!    processed.  The first factor processed must be factor 1!
!
!    Input, integer FACTOR_ORDER, the order of the factor.
!
!    Input, double precision FACTOR_VALUE(FACTOR_ORDER), the factor values
!    for factor FACTOR_INDEX.
!
!    Input, integer FACTOR_NUM, the number of factors.
!
!    Input, integer POINT_NUM, the number of elements in the
!    direct product.
!
!    Input/output, double precision X(FACTOR_NUM,POINT_NUM), the elements of
!    the direct product, which are built up gradually.
!
!  Local Parameters:
!
!    Local, integer START, the first location of a block of 
!    values to set.
!
!    Local, integer CONTIG, the number of consecutive values 
!    to set.
!
!    Local, integer SKIP, the distance from the current value 
!    of START to the next location of a block of values to set.
!
!    Local, integer REP, the number of blocks of values to set.
!
  implicit none

  integer factor_num
  integer factor_order
  integer point_num

  integer , save :: contig
  integer factor_index
  double precision factor_value(factor_order)
  integer j
  integer k
  integer , save :: rep
  integer , save :: skip
  integer start
  double precision x(factor_num,point_num)

  if ( factor_index == 1 ) then
    contig = 1
    skip = 1
    rep = point_num
    x(1:factor_num,1:point_num) = 0.0D+00
  end if

  rep = rep / factor_order
  skip = skip * factor_order

  do j = 1, factor_order

    start = 1 + ( j - 1 ) * contig

    do k = 1, rep
      x(factor_index,start:start+contig-1) = factor_value(j)
      start = start + skip
    end do

  end do

  contig = contig * factor_order
end
