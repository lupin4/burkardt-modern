!> hypercube_grid -- Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module hypercube_grid_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: hypercube_grid, i4vec_product, r8vec_direct_product

contains

  subroutine hypercube_grid ( m, n, ns, a, b, c, x ) &
        bind(C, name="hypercube_grid")

  !*****************************************************************************80
  !
  !! HYPERCUBE_GRID: grid points over the interior of a hypercube in M dimensions.
  !
  !  Discussion:
  !
  !    In M dimensional space, a logically rectangular grid is to be created.
  !    In the I-th dimension, the grid will use S(I) points.
  !    The total number of grid points is
  !      N = product ( 1 <= I <= M ) S(I)
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
  !    29 August 2014
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) M, the spatial dimension.
  !
  !    Input, integer(ip) N, the number of points.
  !    N = product ( 1 <= I <= M ) NS(I).
  !
  !    Input, integer(ip) NS(M), the number of points along
  !    each dimension.
  !
  !    Input, real(dp) A(M), B(M), the endpoints for each dimension.
  !
  !    Input, integer(ip) C(M), the grid centering for each dimension.
  !    1 <= C(*) <= 5.
  !
  !    Output, real(dp) X(M,N) = X(M*S(1),S(2),...,S(M)), the points.
  !

    integer(ip), intent(in), value :: m
    integer(ip), intent(in), value :: n
    real(dp), intent(in) :: a(m)
    real(dp), intent(in) :: b(m)
    integer(ip), intent(in) :: c(m)
    integer(ip) :: i
    integer(ip) :: j
    integer(ip), intent(in) :: ns(m)
    integer(ip) :: s
    real(dp), intent(out) :: x(m,n)
    real(dp), allocatable :: xs(:)
  !
  !  Create the 1D grids in each dimension.
  !
    do i = 1, m

      s = ns(i)

      allocate ( xs(1:s) )

      do j = 1, s

        if ( c(i) == 1 ) then

          if ( s == 1 ) then
            xs(j) = 0.5_dp * ( a(i) + b(i) )
          else
            xs(j) = (   real ( s - j, dp) * a(i)   &
                      + real (     j - 1, dp) * b(i) ) &
                      / real ( s     - 1, dp)
          end if
        else if ( c(i) == 2 ) then
          xs(j) = (   real ( s - j + 1, dp) * a(i)   &
                    + real (     j, dp) * b(i) ) &
                    / real ( s     + 1, dp)
        else if ( c(i) == 3 ) then
          xs(j) = (   real ( s - j + 1, dp) * a(i)   &
                    + real (     j - 1, dp) * b(i) ) &
                    / real ( s, dp)
        else if ( c(i) == 4 ) then
          xs(j) = (   real ( s - j, dp) * a(i)   &
                    + real (     j, dp) * b(i) ) &
                    / real ( s, dp)
        else if ( c(i) == 5 ) then
          xs(j) = (   real ( 2 * s - 2 * j + 1, dp) * a(i)   &
                    + real (         2 * j - 1, dp) * b(i) ) &
                    / real ( 2 * s, dp)
        end if

      end do

      call r8vec_direct_product ( i, s, xs, m, n, x )

      deallocate ( xs )

    end do
  end subroutine hypercube_grid

  pure function i4vec_product ( n, a ) &
        bind(C, name="i4vec_product")

  !*****************************************************************************80
  !
  !! I4VEC_PRODUCT returns the product of the entries of an I4VEC.
  !
  !  Discussion:
  !
  !    An I4VEC is a vector of I4's.
  !
  !    In FORTRAN90, this facility is offered by the built in
  !    PRODUCT function:
  !
  !      I4VEC_PRODUCT ( N, A ) = PRODUCT ( A(1:N) )
  !
  !    In MATLAB, this facility is offered by the built in
  !    PROD function:
  !
  !      I4VEC_PRODUCT ( N, A ) = PROD ( A(1:N) )
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    22 September 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of entries in the array.
  !
  !    Input, integer(ip) A(N), the array.
  !
  !    Output, integer(ip) I4VEC_PRODUCT, the product of the entries.
  !

    integer(ip), intent(in), value :: n
    integer(ip), intent(in) :: a(n)
    integer(ip) :: i4vec_product

    i4vec_product = product ( a(1:n) )
  end function i4vec_product

  subroutine r8vec_direct_product ( factor_index, factor_order, factor_value, &
    factor_num, point_num, x ) &
        bind(C, name="r8vec_direct_product")

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
  !    Input, integer(ip) FACTOR_INDEX, the index of the factor being
  !    processed.  The first factor processed must be factor 1!
  !
  !    Input, integer(ip) FACTOR_ORDER, the order of the factor.
  !
  !    Input, real(dp) FACTOR_VALUE(FACTOR_ORDER), the factor values
  !    for factor FACTOR_INDEX.
  !
  !    Input, integer(ip) FACTOR_NUM, the number of factors.
  !
  !    Input, integer(ip) POINT_NUM, the number of elements in the
  !    direct product.
  !
  !    Input/output, real(dp) X(FACTOR_NUM,POINT_NUM), the elements of
  !    the direct product, which are built up gradually.
  !
  !  Local Parameters:
  !
  !    Local, integer(ip) START, the first location of a block of
  !    values to set.
  !
  !    Local, integer(ip) CONTIG, the number of consecutive values
  !    to set.
  !
  !    Local, integer(ip) SKIP, the distance from the current value
  !    of START to the next location of a block of values to set.
  !
  !    Local, integer(ip) REP, the number of blocks of values to set.
  !

    integer(ip), intent(in), value :: factor_num
    integer(ip), intent(in), value :: factor_order
    integer(ip), intent(in), value :: point_num
    integer(ip), save :: contig
    integer(ip), intent(in), value :: factor_index
    real(dp), intent(in) :: factor_value(factor_order)
    integer(ip) :: j
    integer(ip) :: k
    integer(ip), save :: rep
    integer(ip), save :: skip
    integer(ip) :: start
    real(dp), intent(inout) :: x(factor_num,point_num)

    if ( factor_index == 1 ) then
      contig = 1
      skip = 1
      rep = point_num
      x(1:factor_num,1:point_num) = 0.0_dp
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
  end subroutine r8vec_direct_product

end module hypercube_grid_mod
