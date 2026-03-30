!> felippa — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

subroutine comp_next ( n, k, a, more, h, t )

!*****************************************************************************80
!
!! COMP_NEXT computes the compositions of the integer N into K parts.
!
!  Discussion:
!
!    A composition of the integer N into K parts is an ordered sequence
!    of K nonnegative integers which sum to N.  The compositions (1,2,1)
!    and (1,1,2) are considered to be distinct.
!
!    The routine computes one composition on each call until there are no more.
!    For instance, one composition of 6 into 3 parts is
!    3+2+1, another would be 6+0+0.
!
!    On the first call to this routine, set MORE = FALSE.  The routine
!    will compute the first element in the sequence of compositions, and
!    return it, as well as setting MORE = TRUE.  If more compositions
!    are desired, call again, and again.  Each time, the routine will
!    return with a new composition.
!
!    However, when the LAST composition in the sequence is computed
!    and returned, the routine will reset MORE to FALSE, signaling that
!    the end of the sequence has been reached.
!
!    This routine originally used a SAVE statement to maintain the
!    variables H and T.  I have decided that it is safer
!    to pass these variables as arguments, even though the user should
!    never alter them.  This allows this routine to safely shuffle
!    between several ongoing calculations.
!
!
!    There are 28 compositions of 6 into three parts.  This routine will
!    produce those compositions in the following order:
!
!     I         A
!     -     ---------
!     1     6   0   0
!     2     5   1   0
!     3     4   2   0
!     4     3   3   0
!     5     2   4   0
!     6     1   5   0
!     7     0   6   0
!     8     5   0   1
!     9     4   1   1
!    10     3   2   1
!    11     2   3   1
!    12     1   4   1
!    13     0   5   1
!    14     4   0   2
!    15     3   1   2
!    16     2   2   2
!    17     1   3   2
!    18     0   4   2
!    19     3   0   3
!    20     2   1   3
!    21     1   2   3
!    22     0   3   3
!    23     2   0   4
!    24     1   1   4
!    25     0   2   4
!    26     1   0   5
!    27     0   1   5
!    28     0   0   6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 July 2008
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt
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
!    Input, integer N, the integer whose compositions are desired.
!
!    Input, integer K, the number of parts in the composition.
!
!    Input/output, integer A(K), the parts of the composition.
!
!    Input/output, logical MORE, set by the user to start the
!    computation, and by the routine to terminate it.
!
!    Input/output, integer H, T, two internal parameters needed
!    for the computation.  The user should allocate space for these in the
!    calling program, include them in the calling sequence, but never alter
!    them!
!
  implicit none

  integer k

  integer a(k)
  integer h
  logical more
  integer n
  integer t
!
!  The first computation.
!
  if ( .not. more ) then

    t = n
    h = 0
    a(1) = n
    a(2:k) = 0
!
!  The next computation.
!
  else

    if ( 1 < t ) then
      h = 0
    end if

    h = h + 1
    t = a(h)
    a(h) = 0
    a(1) = t - 1
    a(h+1) = a(h+1) + 1

  end if
!
!  This is the last element of the sequence if all the
!  items are in the last slot.
!
  more = ( a(k) /= n )
end

subroutine hexa_unit_monomial ( expon, value )

!*****************************************************************************80
!
!! HEXA_UNIT_MONOMIAL integrates a monomial over the unit hexahedron.
!
!  Discussion:
!
!    This routine integrates a monomial of the form
!
!      product ( 1 <= dim <= 3 ) x(dim)^expon(dim)
!
!    The combination 0^0 should be treated as 1.
!
!    The integration region is:
!
!    - 1.0 <= X <= 1.0
!    - 1.0 <= Y <= 1.0
!    - 1.0 <= Z <= 1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer EXPON(3), the exponents.
!
!    Output, double precision VALUE, the integral of the monomial.
!
  implicit none

  integer expon(3)
  integer i
  double precision value

  value = 1.0D+00

  do i = 1, 3

    if ( mod ( expon(i), 2 ) == 1 ) then
      value = 0.0D+00
    else if ( expon(i) == -1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HEXA_UNIT_MONOMIAL - Fatal error!'
      write ( *, '(a)' ) '  Exponent of -1 encountered.'
      stop 1
    else
      value = value * 2.0D+00 / real ( expon(i) + 1)
    end if

  end do
end

subroutine hexa_unit_monomial_test ( degree_max )

!*****************************************************************************80
!
!! HEXA_UNIT_MONOMIAL_TEST tests HEXA_UNIT_MONOMIAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DEGREE_MAX, the maximum total degree of the
!    monomials to check.
!
  implicit none

  integer alpha
  integer beta
  integer degree_max
  integer expon(3)
  integer gamma
  double precision hexa_unit_volume
  double precision value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HEXA_UNIT_MONOMIAL_TEST'
  write ( *, '(a)' ) '  For the unit hexahedron,'
  write ( *, '(a)' ) '  HEXA_UNIT_MONOMIAL returns the exact value of the'
  write ( *, '(a)' ) '  integral of X^ALPHA Y^BETA Z^GAMMA'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Volume = ', hexa_unit_volume ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     ALPHA      BETA     GAMMA      INTEGRAL'
  write ( *, '(a)' ) ' '

  do alpha = 0, degree_max
    expon(1) = alpha
    do beta = 0, degree_max - alpha
      expon(2) = beta
      do gamma = 0, degree_max - alpha - beta
        expon(3) = gamma
        call hexa_unit_monomial ( expon, value )
        write ( *, '(2x,i8,2x,i8,2x,i8,2x,g14.6)' ) expon(1:3), value
      end do
    end do
  end do
end

subroutine hexa_unit_quad_test ( degree_max )

!*****************************************************************************80
!
!! HEXA_UNIT_QUAD_TEST tests the rules for the unit hexahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DEGREE_MAX, the maximum total degree of the
!    monomials to check.
!
  implicit none

  integer , parameter :: dim_num = 3

  integer degree_max
  integer expon(dim_num)
  integer h
  double precision hexa_unit_volume
  integer k
  logical more
  integer order
  integer order_1d(dim_num)
  double precision quad
  integer t
  double precision , allocatable :: v(:)
  double precision , allocatable :: w(:)
  double precision , allocatable :: xyz(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HEXA_UNIT_QUAD_TEST'
  write ( *, '(a)' ) '  For the unit hexahedron,'
  write ( *, '(a)' ) '  we approximate monomial integrals with'
  write ( *, '(a)' ) &
    '  HEXA_UNIT_RULE, which returns N1 by N2 by N3 point rules.'

  more = .false.

  do

    call subcomp_next ( degree_max, dim_num, expon, more, h, t )

    if ( any ( mod ( expon(1:dim_num), 2 ) == 1 ) ) then
      if ( .not. more ) then
        exit
      else
        cycle
      end if
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,2x,i2,2x,i2,2x,i2)' ) &
      '  Monomial exponents: ', expon(1:dim_num)
    write ( *, '(a)' ) ' '

    do k = 1, 5

      order_1d(1:dim_num) = k
      order = product ( order_1d(1:dim_num) )
      allocate ( v(1:order) )
      allocate ( w(1:order) )
      allocate ( xyz(1:dim_num,1:order) )
      call hexa_unit_rule ( order_1d, w, xyz )
      call monomial_value ( dim_num, order, expon, xyz, v )
      quad = hexa_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
      write ( *, '(2x,i6,2x,i6,2x,i6,2x,g14.6)' ) order_1d(1:dim_num), quad
      deallocate ( v )
      deallocate ( w )
      deallocate ( xyz )

    end do
!
!  Try a rule of mixed orders.
!
    order_1d(1) = 3
    order_1d(2) = 5
    order_1d(3) = 2
    order = product ( order_1d(1:dim_num) )
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call hexa_unit_rule ( order_1d, w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = hexa_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,i6,2x,i6,2x,g14.6)' ) order_1d(1:dim_num), quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    write ( *, '(a)' ) ' '
    call hexa_unit_monomial ( expon, quad )
    write ( *, '(2x,a,2x,6x,2x,6x,2x,g14.6)' ) ' Exact', quad

    if ( .not. more ) then
      exit
    end if

  end do
end

subroutine hexa_unit_rule ( order_1d, w, xyz )

!*****************************************************************************80
!
!! HEXA_UNIT_RULE returns a quadrature rule for the unit hexahedron.
!
!  Discussion:
!
!    The integration region is:
!
!    - 1.0 <= X <= 1.0
!    - 1.0 <= Y <= 1.0
!    - 1.0 <= Z <= 1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Input, integer ORDER_1D(3), the order of the rule in each
!    dimension.  1 <= ORDER_1D(I) <= 5.
!
!    Output, double precision W(ORDER_1D(1)*ORDER_1D(2)*ORDER_1D(3)), 
!    the weights.
!
!    Output, double precision XYZ(3,ORDER_1D(1)*ORDER_1D(2)*ORDER_1D(3)), 
!    the abscissas.
!
  implicit none

  integer , parameter :: dim_num = 3

  integer dim
  integer order
  integer order_1d(dim_num)
  double precision w(order_1d(1)*order_1d(2)*order_1d(3))
  double precision , allocatable :: w_1d(:)
  double precision , allocatable :: x_1d(:)
  double precision xyz(3,order_1d(1)*order_1d(2)*order_1d(3))

  order = product ( order_1d(1:dim_num) )

  do dim = 1, dim_num

    allocate ( w_1d(order_1d(dim)) )
    allocate ( x_1d(order_1d(dim)) )

    if ( order_1d(dim) == 1 ) then
      call line_unit_o01 ( w_1d, x_1d )
    else if ( order_1d(dim) == 2 ) then
      call line_unit_o02 ( w_1d, x_1d )
    else if ( order_1d(dim) == 3 ) then
      call line_unit_o03 ( w_1d, x_1d )
    else if ( order_1d(dim) == 4 ) then
      call line_unit_o04 ( w_1d, x_1d )
    else if ( order_1d(dim) == 5 ) then
      call line_unit_o05 ( w_1d, x_1d )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HEXA_UNIT_RULE - Fatal error!'
      write ( *, '(a)' ) '  Illegal value of ORDER_1D(*).'
      stop 1
    end if

    call r8vec_direct_product ( dim, order_1d(dim), x_1d, &
      dim_num, order, xyz )

    call r8vec_direct_product2 ( dim, order_1d(dim), w_1d, &
      dim_num, order, w )

    deallocate ( w_1d )
    deallocate ( x_1d )

  end do
end

function hexa_unit_volume ( )

!*****************************************************************************80
!
!! HEXA_UNIT_VOLUME: volume of a unit hexahedron.
!
!  Discussion:
!
!    The integration region is:
!
!    - 1.0 <= X <= 1.0
!    - 1.0 <= Y <= 1.0
!    - 1.0 <= Z <= 1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, double precision HEXA_UNIT_VOLUME, the volume.
!
  implicit none

  double precision hexa_unit_volume

  hexa_unit_volume = 8.0D+00
end

subroutine line_unit_monomial ( alpha, value )

!*****************************************************************************80
!
!! LINE_UNIT_MONOMIAL: monomial integral in a unit line.
!
!  Discussion:
!
!    This function returns the integral of X^ALPHA over the unit line.
!
!    The integration region is:
!
!    - 1.0 <= X <= 1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ALPHA, the exponent of X.
!    ALPHA must not be -1.
!
!    Output, double precision value, the integral of the monomial.
!
  implicit none

  integer alpha
  double precision value

  if ( alpha == - 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE_UNIT_MONOMIAL - Fatal error!'
    write ( *, '(a)' ) '  ALPHA = -1 is not a legal input.'
    stop 1
  else if ( mod ( alpha, 2 ) == 1 ) then
    value = 0.0D+00
  else
    value = 2.0D+00 / real ( alpha + 1)
  end if
end

subroutine line_unit_monomial_test ( degree_max )

!*****************************************************************************80
!
!! LINE_UNIT_MONOMIAL_TEST tests LINE_UNIT_MONOMIAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DEGREE_MAX, the maximum total degree of the
!    monomials to check.
!
  implicit none

  integer alpha
  integer degree_max
  double precision line_unit_volume
  double precision value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LINE_UNIT_MONOMIAL_TEST'
  write ( *, '(a)' ) '  For the unit line,'
  write ( *, '(a)' ) '  LINE_UNIT_MONOMIAL returns the exact value of the'
  write ( *, '(a)' ) '  integral of X^ALPHA'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Volume = ', line_unit_volume ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     ALPHA      INTEGRAL'
  write ( *, '(a)' ) ' '

  do alpha = 0, degree_max
    call line_unit_monomial ( alpha, value )
    write ( *, '(2x,i8,2x,g14.6)' ) alpha, value
  end do
end

subroutine line_unit_o01 ( w, x )

!*****************************************************************************80
!
!! LINE_UNIT_O01 returns a 1 point quadrature rule for the unit line.
!
!  Discussion:
!
!    The integration region is:
!
!    - 1.0 <= X <= 1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(1), the weights.
!
!    Output, double precision X(1), the abscissas.
!
  implicit none

  integer , parameter :: order = 1

  double precision line_unit_volume
  double precision w(order)
  double precision :: w_save(1) = (/ &
    2.0D+00 /)
  double precision x(order)
  double precision :: x_save(1) = (/ &
    0.0D+00 /)

  w(1:order) = w_save(1:order) / line_unit_volume ( )
  x(1:order) = x_save(1:order)
end

subroutine line_unit_o02 ( w, x )

!*****************************************************************************80
!
!! LINE_UNIT_O02 returns a 2 point quadrature rule for the unit line.
!
!  Discussion:
!
!    The integration region is:
!
!    - 1.0 <= X <= 1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(2), the weights.
!
!    Output, double precision X(2), the abscissas.
!
  implicit none

  integer , parameter :: order = 2

  double precision line_unit_volume
  double precision w(order)
  double precision :: w_save(2) = (/ &
    1.0000000000000000000D+00, &
    1.0000000000000000000D+00 /)
  double precision x(order)
  double precision :: x_save(2) = (/ &
    -0.57735026918962576451D+00, &
     0.57735026918962576451D+00 /)

  w(1:order) = w_save(1:order) / line_unit_volume ( )
  x(1:order) = x_save(1:order)
end

subroutine line_unit_o03 ( w, x )

!*****************************************************************************80
!
!! LINE_UNIT_O03 returns a 3 point quadrature rule for the unit line.
!
!  Discussion:
!
!    The integration region is:
!
!    - 1.0 <= X <= 1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(3), the weights.
!
!    Output, double precision X(3), the abscissas.
!
  implicit none

  integer , parameter :: order = 3

  double precision line_unit_volume
  double precision w(order)
  double precision :: w_save(3) = (/ &
    0.55555555555555555556D+00, &
    0.88888888888888888889D+00, &
    0.55555555555555555556D+00 /)
  double precision x(order)
  double precision :: x_save(3) = (/ &
    -0.77459666924148337704D+00, &
     0.00000000000000000000D+00, &
     0.77459666924148337704D+00 /)

  w(1:order) = w_save(1:order) / line_unit_volume ( )
  x(1:order) = x_save(1:order)
end

subroutine line_unit_o04 ( w, x )

!*****************************************************************************80
!
!! LINE_UNIT_O04 returns a 4 point quadrature rule for the unit line.
!
!  Discussion:
!
!    The integration region is:
!
!    - 1.0 <= X <= 1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(4), the weights.
!
!    Output, double precision X(4), the abscissas.
!
  implicit none

  integer , parameter :: order = 4

  double precision line_unit_volume
  double precision w(order)
  double precision :: w_save(4) = (/ &
    0.34785484513745385737D+00, &
    0.65214515486254614263D+00, &
    0.65214515486254614263D+00, &
    0.34785484513745385737D+00 /)
  double precision x(order)
  double precision :: x_save(4) = (/ &
    -0.86113631159405257522D+00, &
    -0.33998104358485626480D+00, &
     0.33998104358485626480D+00, &
     0.86113631159405257522D+00 /)

  w(1:order) = w_save(1:order) / line_unit_volume ( )
  x(1:order) = x_save(1:order)
end

subroutine line_unit_o05 ( w, x )

!*****************************************************************************80
!
!! LINE_UNIT_O05 returns a 5 point quadrature rule for the unit line.
!
!  Discussion:
!
!    The integration region is:
!
!    - 1.0 <= X <= 1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(5), the weights.
!
!    Output, double precision X(5), the abscissas.
!
  implicit none

  integer , parameter :: order = 5

  double precision line_unit_volume
  double precision w(order)
  double precision :: w_save(5) = (/ &
    0.23692688505618908751D+00, &
    0.47862867049936646804D+00, &
    0.56888888888888888889D+00, &
    0.47862867049936646804D+00, &
    0.23692688505618908751D+00 /)
  double precision x(order)
  double precision :: x_save(5) = (/ &
    -0.90617984593866399280D+00, &
    -0.53846931010568309104D+00, &
     0.00000000000000000000D+00, &
     0.53846931010568309104D+00, &
     0.90617984593866399280D+00 /)

  w(1:order) = w_save(1:order) / line_unit_volume ( )
  x(1:order) = x_save(1:order)
end

subroutine line_unit_quad_test ( degree_max )

!*****************************************************************************80
!
!! LINE_UNIT_QUAD_TEST tests the rules for the unit line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DEGREE_MAX, the maximum total degree of the
!    monomials to check.
!
  implicit none

  integer , parameter :: dim_num = 1

  integer degree_max
  integer expon(dim_num)
  integer h
  double precision line_unit_volume
  logical more
  integer order
  double precision quad
  integer t
  double precision , allocatable :: v(:)
  double precision , allocatable :: w(:)
  double precision , allocatable :: x(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LINE_UNIT_QUAD_TEST'
  write ( *, '(a)' ) '  For the unit line,'
  write ( *, '(a)' ) '  we approximate monomial integrals with:'
  write ( *, '(a)' ) '  LINE_UNIT_O01, a 1 point rule.'
  write ( *, '(a)' ) '  LINE_UNIT_O02, a 2 point rule.'
  write ( *, '(a)' ) '  LINE_UNIT_O03, a 3 point rule.'
  write ( *, '(a)' ) '  LINE_UNIT_O04, a 4 point rule.'
  write ( *, '(a)' ) '  LINE_UNIT_O05, a 5 point rule.'

  more = .false.

  do

    call subcomp_next ( degree_max, dim_num, expon, more, h, t )

    if ( any ( mod ( expon(1:dim_num), 2 ) == 1 ) ) then
      if ( .not. more ) then
        exit
      else
        cycle
      end if
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,2x,i2)' ) '  Monomial exponents: ', expon(1:dim_num)
    write ( *, '(a)' ) ' '

    order = 1
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( x(1:dim_num,1:order) )
    call line_unit_o01 ( w, x )
    call monomial_value ( dim_num, order, expon, x, v )
    quad = line_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    order = 2
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( x(1:dim_num,1:order) )
    call line_unit_o02 ( w, x )
    call monomial_value ( dim_num, order, expon, x, v )
    quad = line_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    order = 3
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( x(1:dim_num,1:order) )
    call line_unit_o03 ( w, x )
    call monomial_value ( dim_num, order, expon, x, v )
    quad = line_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    order = 4
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( x(1:dim_num,1:order) )
    call line_unit_o04 ( w, x )
    call monomial_value ( dim_num, order, expon, x, v )
    quad = line_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    order = 5
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( x(1:dim_num,1:order) )
    call line_unit_o05 ( w, x )
    call monomial_value ( dim_num, order, expon, x, v )
    quad = line_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    write ( *, '(a)' ) ' '
    call line_unit_monomial ( expon(1), quad )
    write ( *, '(2x,a,2x,g14.6)' ) ' Exact', quad

    if ( .not. more ) then
      exit
    end if

  end do
end

function line_unit_volume ( )

!*****************************************************************************80
!
!! LINE_UNIT_VOLUME: volume of a unit line.
!
!  Discussion:
!
!    The integration region is:
!
!    - 1.0 <= X <= 1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, double precision LINE_UNIT_VOLUME, the volume.
!
  implicit none

  double precision line_unit_volume

  line_unit_volume = 2.0D+00
end

subroutine monomial_value ( dim_num, point_num, expon, x, v )

!*****************************************************************************80
!
!! MONOMIAL_VALUE evaluates a monomial.
!
!  Discussion:
!
!    F(X) = product ( 1 <= DIM <= DIM_NUM ) X(I)^EXPON(I)
!
!    with the convention that 0^0 = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, integer EXPON(DIM_NUM), the exponents.
!
!    Input, double precision X(DIM_NUM,POINT_NUM), the evaluation points.
!
!    Output, double precision V(POINT_NUM), the monomial values.
!
  implicit none

  integer dim_num
  integer point_num

  integer dim
  integer expon(dim_num)
  double precision v(point_num)
  double precision x(dim_num,point_num)

  v(1:point_num) = 1.0D+00

  do dim = 1, dim_num
    if ( expon(dim) /= 0.0D+00 ) then
      v(1:point_num) = v(1:point_num) * x(dim,1:point_num)**expon(dim)
    end if
  end do
end

subroutine pyra_unit_monomial ( expon, value )

!*****************************************************************************80
!
!! PYRA_UNIT_MONOMIAL: monomial integral in a unit pyramid.
!
!  Discussion:
!
!    This routine returns the integral of
!
!      product ( 1 <= I <= 3 ) X(I)^EXPON(I)
!
!    over the unit pyramid.
!
!    The integration region is:
!
!    - ( 1 - Z ) <= X <= 1 - Z
!    - ( 1 - Z ) <= Y <= 1 - Z
!              0 <= Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer EXPON(3), the exponents.
!
!    Output, double precision VALUE, the integral of the monomial.
!
  implicit none

  integer expon(3)
  integer i
  integer i_hi
  double precision r8_choose
  double precision r8_mop
  double precision value

  value = 0.0D+00

  if ( mod ( expon(1), 2 ) == 0 .and. mod ( expon(2), 2 ) == 0 ) then

    i_hi = 2 + expon(1) + expon(2)

    do i = 0, i_hi
      value = value + r8_mop ( i ) * r8_choose ( i_hi, i ) &
      / real ( i + expon(3) + 1)
    end do

    value = value &
          * 2.0D+00 / real ( expon(1) + 1) &
          * 2.0D+00 / real ( expon(2) + 1)

  end if
end

subroutine pyra_unit_monomial_test ( degree_max )

!*****************************************************************************80
!
!! PYRA_UNIT_MONOMIAL_TEST tests PYRA_UNIT_MONOMIAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DEGREE_MAX, the maximum total degree of the
!    monomials to check.
!
  implicit none

  integer alpha
  integer beta
  integer degree_max
  integer expon(3)
  integer gamma
  double precision pyra_unit_volume
  double precision value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PYRA_UNIT_MONOMIAL_TEST'
  write ( *, '(a)' ) '  For the unit pyramid,'
  write ( *, '(a)' ) '  PYRA_UNIT_MONOMIAL returns the exact value of the'
  write ( *, '(a)' ) '  integral of X^ALPHA Y^BETA Z^GAMMA'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Volume = ', pyra_unit_volume ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     ALPHA      BETA     GAMMA      INTEGRAL'
  write ( *, '(a)' ) ' '

  do alpha = 0, degree_max
    expon(1) = alpha
    do beta = 0, degree_max - alpha
      expon(2) = beta
      do gamma = 0, degree_max - alpha - beta
        expon(3) = gamma
        call pyra_unit_monomial ( expon, value )
        write ( *, '(2x,i8,2x,i8,2x,i8,2x,g14.6)' ) expon(1:3), value
      end do
    end do
  end do
end

subroutine pyra_unit_o01 ( w, xyz )

!*****************************************************************************80
!
!! PYRA_UNIT_O01 returns a 1 point quadrature rule for the unit pyramid.
!
!  Discussion:
!
!    The integration region is:
!
!    - ( 1 - Z ) <= X <= 1 - Z
!    - ( 1 - Z ) <= Y <= 1 - Z
!              0 <= Z <= 1.
!
!    When Z is zero, the integration region is a square lying in the (X,Y)
!    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the
!    radius of the square diminishes, and when Z reaches 1, the square has
!    contracted to the single point (0,0,1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(1), the weights.
!
!    Output, double precision XYZ(3,1), the abscissas.
!
  implicit none

  integer , parameter :: order = 1

  double precision w(order)
  double precision :: w_save(1) = (/ &
    1.0D+00 /)
  double precision xyz(3,order)
  double precision :: xyz_save(3,1) = reshape ( (/ &
    0.0D+00, 0.0D+00, 0.25D+00 /), &
  (/ 3, 1 /) )

  w(1:order) = w_save(1:order)
  xyz(1:3,1:order) = xyz_save(1:3,1:order)
end

subroutine pyra_unit_o05 ( w, xyz )

!*****************************************************************************80
!
!! PYRA_UNIT_O05 returns a 5 point quadrature rule for the unit pyramid.
!
!  Discussion:
!
!    The integration region is:
!
!    - ( 1 - Z ) <= X <= 1 - Z
!    - ( 1 - Z ) <= Y <= 1 - Z
!              0 <= Z <= 1.
!
!    When Z is zero, the integration region is a square lying in the (X,Y)
!    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the
!    radius of the square diminishes, and when Z reaches 1, the square has
!    contracted to the single point (0,0,1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(5), the weights.
!
!    Output, double precision XYZ(3,5), the abscissas.
!
  implicit none

  integer , parameter :: order = 5

  double precision w(order)
  double precision :: w_save(5) = (/ &
   0.21093750000000000000D+00, &
   0.21093750000000000000D+00, &
   0.21093750000000000000D+00, &
   0.21093750000000000000D+00, &
   0.15625000000000000000D+00 /)
  double precision xyz(3,order)
  double precision :: xyz_save(3,5) = reshape ( (/ &
  -0.48686449556014765641D+00, &
  -0.48686449556014765641D+00, &
   0.16666666666666666667D+00, &
   0.48686449556014765641D+00, &
  -0.48686449556014765641D+00, &
   0.16666666666666666667D+00, &
   0.48686449556014765641D+00, &
   0.48686449556014765641D+00, &
   0.16666666666666666667D+00, &
  -0.48686449556014765641D+00, &
   0.48686449556014765641D+00, &
   0.16666666666666666667D+00, &
   0.00000000000000000000D+00, &
   0.00000000000000000000D+00, &
   0.70000000000000000000D+00 /), &
  (/ 3, 5 /) )

  w(1:order) = w_save(1:order)
  xyz(1:3,1:order) = xyz_save(1:3,1:order)
end

subroutine pyra_unit_o06 ( w, xyz )

!*****************************************************************************80
!
!! PYRA_UNIT_O06 returns a 6 point quadrature rule for the unit pyramid.
!
!  Discussion:
!
!    The integration region is:
!
!    - ( 1 - Z ) <= X <= 1 - Z
!    - ( 1 - Z ) <= Y <= 1 - Z
!              0 <= Z <= 1.
!
!    When Z is zero, the integration region is a square lying in the (X,Y)
!    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the
!    radius of the square diminishes, and when Z reaches 1, the square has
!    contracted to the single point (0,0,1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(6), the weights.
!
!    Output, double precision XYZ(3,6), the abscissas.
!
  implicit none

  integer , parameter :: order = 6

  double precision w(order)
  double precision :: w_save(6) = (/ &
   0.21000000000000000000D+00, &
   0.21000000000000000000D+00, &
   0.21000000000000000000D+00, &
   0.21000000000000000000D+00, &
   0.06000000000000000000D+00, &
   0.10000000000000000000D+00 /)
  double precision xyz(3,order)
  double precision :: xyz_save(3,6) = reshape ( (/ &
  -0.48795003647426658968D+00, &
  -0.48795003647426658968D+00, &
   0.16666666666666666667D+00, &
   0.48795003647426658968D+00, &
  -0.48795003647426658968D+00, &
   0.16666666666666666667D+00, &
   0.48795003647426658968D+00, &
   0.48795003647426658968D+00, &
   0.16666666666666666667D+00, &
  -0.48795003647426658968D+00, &
   0.48795003647426658968D+00, &
   0.16666666666666666667D+00, &
   0.00000000000000000000D+00, &
   0.00000000000000000000D+00, &
   0.58333333333333333333D+00, &
   0.00000000000000000000D+00, &
   0.00000000000000000000D+00, &
   0.75000000000000000000D+00 /), &
  (/ 3, 6 /) )

  w(1:order) = w_save(1:order)
  xyz(1:3,1:order) = xyz_save(1:3,1:order)
end

subroutine pyra_unit_o08 ( w, xyz )

!*****************************************************************************80
!
!! PYRA_UNIT_O08 returns an 8 point quadrature rule for the unit pyramid.
!
!  Discussion:
!
!    The integration region is:
!
!    - ( 1 - Z ) <= X <= 1 - Z
!    - ( 1 - Z ) <= Y <= 1 - Z
!              0 <= Z <= 1.
!
!    When Z is zero, the integration region is a square lying in the (X,Y)
!    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the
!    radius of the square diminishes, and when Z reaches 1, the square has
!    contracted to the single point (0,0,1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(8), the weights.
!
!    Output, double precision XYZ(3,8), the abscissas.
!
  implicit none

  integer , parameter :: order = 8

  double precision w(order)
  double precision :: w_save(8) = (/ &
   0.075589411559869072938D+00, &
   0.075589411559869072938D+00, &
   0.075589411559869072938D+00, &
   0.075589411559869072938D+00, &
   0.17441058844013092706D+00, &
   0.17441058844013092706D+00, &
   0.17441058844013092706D+00, &
   0.17441058844013092706D+00 /)
  double precision xyz(3,order)
  double precision :: xyz_save(3,8) = reshape ( (/ &
  -0.26318405556971359557D+00, &
  -0.26318405556971359557D+00, &
   0.54415184401122528880D+00, &
   0.26318405556971359557D+00, &
  -0.26318405556971359557D+00, &
   0.54415184401122528880D+00, &
   0.26318405556971359557D+00, &
   0.26318405556971359557D+00, &
   0.54415184401122528880D+00, &
  -0.26318405556971359557D+00, &
   0.26318405556971359557D+00, &
   0.54415184401122528880D+00, &
  -0.50661630334978742377D+00, &
  -0.50661630334978742377D+00, &
   0.12251482265544137787D+00, &
   0.50661630334978742377D+00, &
  -0.50661630334978742377D+00, &
   0.12251482265544137787D+00, &
   0.50661630334978742377D+00, &
   0.50661630334978742377D+00, &
   0.12251482265544137787D+00, &
  -0.50661630334978742377D+00, &
   0.50661630334978742377D+00, &
   0.12251482265544137787D+00 /), &
  (/ 3, 8 /) )

  w(1:order) = w_save(1:order)
  xyz(1:3,1:order) = xyz_save(1:3,1:order)
end

subroutine pyra_unit_o08b ( w, xyz )

!*****************************************************************************80
!
!! PYRA_UNIT_O08B returns an 8 point quadrature rule for the unit pyramid.
!
!  Discussion:
!
!    The integration region is:
!
!    - ( 1 - Z ) <= X <= 1 - Z
!    - ( 1 - Z ) <= Y <= 1 - Z
!              0 <= Z <= 1.
!
!    When Z is zero, the integration region is a square lying in the (X,Y)
!    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the
!    radius of the square diminishes, and when Z reaches 1, the square has
!    contracted to the single point (0,0,1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(8), the weights.
!
!    Output, double precision XYZ(3,8), the abscissas.
!
  implicit none

  integer , parameter :: order = 1

  double precision w(order)
  double precision :: w_save(8) = (/ &
   0.16438287736328777572D+00, &
   0.16438287736328777572D+00, &
   0.16438287736328777572D+00, &
   0.16438287736328777572D+00, &
   0.085617122636712224276D+00, &
   0.085617122636712224276D+00, &
   0.085617122636712224276D+00, &
   0.085617122636712224276D+00 /)
  double precision xyz(3,order)
  double precision :: xyz_save(3,8) = reshape ( (/ &
  -0.51197009372656270107D+00, &
  -0.51197009372656270107D+00, &
   0.11024490204163285720D+00, &
   0.51197009372656270107D+00, &
  -0.51197009372656270107D+00, &
   0.11024490204163285720D+00, &
   0.51197009372656270107D+00, &
   0.51197009372656270107D+00, &
   0.11024490204163285720D+00, &
  -0.51197009372656270107D+00, &
   0.51197009372656270107D+00, &
   0.11024490204163285720D+00, &
  -0.28415447557052037456D+00, &
  -0.28415447557052037456D+00, &
   0.518326526529795714229D+00, &
   0.28415447557052037456D+00, &
  -0.28415447557052037456D+00, &
   0.518326526529795714229D+00, &
   0.28415447557052037456D+00, &
   0.28415447557052037456D+00, &
   0.518326526529795714229D+00, &
  -0.28415447557052037456D+00, &
   0.28415447557052037456D+00, &
   0.518326526529795714229D+00 /), &
  (/ 3, 8 /) )

  w(1:order) = w_save(1:order)
  xyz(1:3,1:order) = xyz_save(1:3,1:order)
end

subroutine pyra_unit_o09 ( w, xyz )

!*****************************************************************************80
!
!! PYRA_UNIT_O09 returns a 9 point quadrature rule for the unit pyramid.
!
!  Discussion:
!
!    The integration region is:
!
!    - ( 1 - Z ) <= X <= 1 - Z
!    - ( 1 - Z ) <= Y <= 1 - Z
!              0 <= Z <= 1.
!
!    When Z is zero, the integration region is a square lying in the (X,Y)
!    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the
!    radius of the square diminishes, and when Z reaches 1, the square has
!    contracted to the single point (0,0,1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(9), the weights.
!
!    Output, double precision XYZ(3,9), the abscissas.
!
  implicit none

  integer , parameter :: order = 9

  double precision w(order)
  double precision :: w_save(9) = (/ &
   0.13073389672275944791D+00, &
   0.13073389672275944791D+00, &
   0.13073389672275944791D+00, &
   0.13073389672275944791D+00, &
   0.10989110327724055209D+00, &
   0.10989110327724055209D+00, &
   0.10989110327724055209D+00, &
   0.10989110327724055209D+00, &
   0.03750000000000000000D+00 /)
  double precision xyz(3,order)
  double precision :: xyz_save(3,9) = reshape ( (/ &
  -0.52966422253852215131D+00, &
  -0.52966422253852215131D+00, &
   0.08176876558246862335D+00, &
   0.52966422253852215131D+00, &
  -0.52966422253852215131D+00, &
   0.08176876558246862335D+00, &
   0.52966422253852215131D+00, &
   0.52966422253852215131D+00, &
   0.08176876558246862335D+00, &
  -0.52966422253852215131D+00, &
   0.52966422253852215131D+00, &
   0.08176876558246862335D+00, &
  -0.34819753825720418039D+00, &
  -0.34819753825720418039D+00, &
   0.400374091560388519511D+00, &
   0.34819753825720418039D+00, &
  -0.34819753825720418039D+00, &
   0.400374091560388519511D+00, &
   0.34819753825720418039D+00, &
   0.34819753825720418039D+00, &
   0.400374091560388519511D+00, &
  -0.34819753825720418039D+00, &
   0.34819753825720418039D+00, &
   0.400374091560388519511D+00, &
   0.00000000000000000000D+00, &
   0.00000000000000000000D+00, &
   0.83333333333333333333D+00 /), &
  (/ 3, 9 /) )

  w(1:order) = w_save(1:order)
  xyz(1:3,1:order) = xyz_save(1:3,1:order)
end

subroutine pyra_unit_o13 ( w, xyz )

!*****************************************************************************80
!
!! PYRA_UNIT_O13 returns a 13 point quadrature rule for the unit pyramid.
!
!  Discussion:
!
!    The integration region is:
!
!    - ( 1 - Z ) <= X <= 1 - Z
!    - ( 1 - Z ) <= Y <= 1 - Z
!              0 <= Z <= 1.
!
!    When Z is zero, the integration region is a square lying in the (X,Y)
!    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the
!    radius of the square diminishes, and when Z reaches 1, the square has
!    contracted to the single point (0,0,1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(13), the weights.
!
!    Output, double precision XYZ(3,13), the abscissas.
!
  implicit none

  integer , parameter :: order = 13

  double precision w(order)
  double precision :: w_save(13) = (/ &
   0.063061594202898550725D+00, &
   0.063061594202898550725D+00, &
   0.063061594202898550725D+00, &
   0.063061594202898550725D+00, &
   0.042101946815575556199D+00, &
   0.042101946815575556199D+00, &
   0.042101946815575556199D+00, &
   0.042101946815575556199D+00, &
   0.13172030707666776585D+00, &
   0.13172030707666776585D+00, &
   0.13172030707666776585D+00, &
   0.13172030707666776585D+00, &
   0.05246460761943250889D+00 /)
  double precision xyz(3,order)
  double precision :: xyz_save(3,13) = reshape ( (/ &
  -0.38510399211870384331D+00, &
  -0.38510399211870384331D+00, &
  0.428571428571428571429D+00, &
   0.38510399211870384331D+00, &
  -0.38510399211870384331D+00, &
  0.428571428571428571429D+00, &
   0.38510399211870384331D+00, &
   0.38510399211870384331D+00, &
  0.428571428571428571429D+00, &
  -0.38510399211870384331D+00, &
   0.38510399211870384331D+00, &
  0.428571428571428571429D+00, &
  -0.40345831960728204766D+00, &
   0.00000000000000000000D+00, &
  0.33928571428571428571D+00,  &
   0.40345831960728204766D+00, &
   0.00000000000000000000D+00, &
  0.33928571428571428571D+00,  &
   0.00000000000000000000D+00, &
  -0.40345831960728204766D+00, &
  0.33928571428571428571D+00,  &
   0.00000000000000000000D+00, &
   0.40345831960728204766D+00, &
  0.33928571428571428571D+00,  &
  -0.53157877436961973359D+00, &
  -0.53157877436961973359D+00, &
  0.08496732026143790850D+00,  &
   0.53157877436961973359D+00, &
  -0.53157877436961973359D+00, &
  0.08496732026143790850D+00,  &
   0.53157877436961973359D+00, &
   0.53157877436961973359D+00, &
  0.08496732026143790850D+00,  &
  -0.53157877436961973359D+00, &
   0.53157877436961973359D+00, &
  0.08496732026143790850D+00,  &
   0.00000000000000000000D+00, &
   0.00000000000000000000D+00, &
  0.76219701803768503595D+00 /), &
  (/ 3, 13 /) )

  w(1:order) = w_save(1:order)
  xyz(1:3,1:order) = xyz_save(1:3,1:order)
end

subroutine pyra_unit_o18 ( w, xyz )

!*****************************************************************************80
!
!! PYRA_UNIT_O18 returns an 18 point quadrature rule for the unit pyramid.
!
!  Discussion:
!
!    The integration region is:
!
!    - ( 1 - Z ) <= X <= 1 - Z
!    - ( 1 - Z ) <= Y <= 1 - Z
!              0 <= Z <= 1.
!
!    When Z is zero, the integration region is a square lying in the (X,Y)
!    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the
!    radius of the square diminishes, and when Z reaches 1, the square has
!    contracted to the single point (0,0,1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(18), the weights.
!
!    Output, double precision XYZ(3,18), the abscissas.
!
  implicit none

  integer , parameter :: order = 18

  double precision w(order)
  double precision :: w_save(18) = (/ &
   0.023330065296255886709D+00, &
   0.037328104474009418735D+00, &
   0.023330065296255886709D+00, &
   0.037328104474009418735D+00, &
   0.059724967158415069975D+00, &
   0.037328104474009418735D+00, &
   0.023330065296255886709D+00, &
   0.037328104474009418735D+00, &
   0.023330065296255886709D+00, &
   0.05383042853090460712D+00, &
   0.08612868564944737139D+00, &
   0.05383042853090460712D+00, &
   0.08612868564944737139D+00, &
   0.13780589703911579422D+00, &
   0.08612868564944737139D+00, &
   0.05383042853090460712D+00, &
   0.08612868564944737139D+00, &
   0.05383042853090460712D+00 /)
  double precision xyz(3,order)
  double precision :: xyz_save(3,18) = reshape ( (/ &
  -0.35309846330877704481D+00, &
  -0.35309846330877704481D+00, &
  0.544151844011225288800D+00, &
   0.00000000000000000000D+00, &
  -0.35309846330877704481D+00, &
  0.544151844011225288800D+00, &
   0.35309846330877704481D+00, &
  -0.35309846330877704481D+00, &
  0.544151844011225288800D+00, &
  -0.35309846330877704481D+00, &
   0.00000000000000000000D+00, &
  0.544151844011225288800D+00, &
   0.00000000000000000000D+00, &
   0.00000000000000000000D+00, &
  0.544151844011225288800D+00, &
   0.35309846330877704481D+00, &
   0.00000000000000000000D+00, &
  0.544151844011225288800D+00, &
  -0.35309846330877704481D+00, &
   0.35309846330877704481D+00, &
  0.544151844011225288800D+00, &
   0.00000000000000000000D+00, &
   0.35309846330877704481D+00, &
  0.544151844011225288800D+00, &
   0.35309846330877704481D+00, &
   0.35309846330877704481D+00, &
  0.544151844011225288800D+00, &
  -0.67969709567986745790D+00, &
  -0.67969709567986745790D+00, &
  0.12251482265544137787D+00, &
   0.00000000000000000000D+00, &
  -0.67969709567986745790D+00, &
  0.12251482265544137787D+00, &
   0.67969709567986745790D+00, &
  -0.67969709567986745790D+00, &
  0.12251482265544137787D+00, &
  -0.67969709567986745790D+00, &
   0.00000000000000000000D+00, &
  0.12251482265544137787D+00, &
   0.00000000000000000000D+00, &
   0.00000000000000000000D+00, &
  0.12251482265544137787D+00, &
   0.67969709567986745790D+00, &
   0.00000000000000000000D+00, &
  0.12251482265544137787D+00, &
  -0.67969709567986745790D+00, &
   0.67969709567986745790D+00, &
  0.12251482265544137787D+00, &
   0.00000000000000000000D+00, &
   0.67969709567986745790D+00, &
  0.12251482265544137787D+00, &
   0.67969709567986745790D+00, &
   0.67969709567986745790D+00, &
  0.12251482265544137787D+00 /), &
  (/ 3, 18 /) )

  w(1:order) = w_save(1:order)
  xyz(1:3,1:order) = xyz_save(1:3,1:order)
end

subroutine pyra_unit_o27 ( w, xyz )

!*****************************************************************************80
!
!! PYRA_UNIT_O27 returns a 27 point quadrature rule for the unit pyramid.
!
!  Discussion:
!
!    The integration region is:
!
!    - ( 1 - Z ) <= X <= 1 - Z
!    - ( 1 - Z ) <= Y <= 1 - Z
!              0 <= Z <= 1.
!
!    When Z is zero, the integration region is a square lying in the (X,Y)
!    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the
!    radius of the square diminishes, and when Z reaches 1, the square has
!    contracted to the single point (0,0,1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(27), the weights.
!
!    Output, double precision XYZ(3,27), the abscissas.
!
  implicit none

  integer , parameter :: order = 27

  double precision w(order)
  double precision :: w_save(27) = (/ &
   0.036374157653908938268D+00, &
   0.05819865224625430123D+00, &
   0.036374157653908938268D+00, &
   0.05819865224625430123D+00, &
   0.09311784359400688197D+00, &
   0.05819865224625430123D+00, &
   0.036374157653908938268D+00, &
   0.05819865224625430123D+00, &
   0.036374157653908938268D+00, &
   0.033853303069413431019D+00, &
   0.054165284911061489631D+00, &
   0.033853303069413431019D+00, &
   0.054165284911061489631D+00, &
   0.08666445585769838341D+00, &
   0.054165284911061489631D+00, &
   0.033853303069413431019D+00, &
   0.054165284911061489631D+00, &
   0.033853303069413431019D+00, &
   0.006933033103838124540D+00, &
   0.011092852966140999264D+00, &
   0.006933033103838124540D+00, &
   0.011092852966140999264D+00, &
   0.017748564745825598822D+00, &
   0.011092852966140999264D+00, &
   0.006933033103838124540D+00, &
   0.011092852966140999264D+00, &
   0.006933033103838124540D+00 /)
  double precision xyz(3,order)
  double precision :: xyz_save(3,27) = reshape ( (/ &
  -0.7180557413198889387D+00, &
   -0.7180557413198889387D+00, &
   0.07299402407314973216D+00, &
   0.00000000000000000000D+00, &
  -0.7180557413198889387D+00, &
   0.07299402407314973216D+00, &
   0.7180557413198889387D+00, &
   -0.7180557413198889387D+00, &
   0.07299402407314973216D+00, &
  -0.7180557413198889387D+00, &
    0.00000000000000000000D+00, &
  0.07299402407314973216D+00, &
   0.00000000000000000000D+00, &
   0.00000000000000000000D+00, &
  0.07299402407314973216D+00, &
   0.7180557413198889387D+00, &
    0.00000000000000000000D+00, &
  0.07299402407314973216D+00, &
  -0.7180557413198889387D+00, &
    0.7180557413198889387D+00, &
   0.07299402407314973216D+00, &
   0.00000000000000000000D+00, &
   0.7180557413198889387D+00, &
   0.07299402407314973216D+00, &
   0.7180557413198889387D+00, &
    0.7180557413198889387D+00, &
   0.07299402407314973216D+00, &
  -0.50580870785392503961D+00, &
  -0.50580870785392503961D+00, &
  0.34700376603835188472D+00, &
   0.00000000000000000000D+00, &
  -0.50580870785392503961D+00, &
  0.34700376603835188472D+00, &
   0.50580870785392503961D+00, &
  -0.50580870785392503961D+00, &
  0.34700376603835188472D+00, &
  -0.50580870785392503961D+00, &
   0.00000000000000000000D+00, &
  0.34700376603835188472D+00, &
   0.00000000000000000000D+00, &
   0.00000000000000000000D+00, &
  0.34700376603835188472D+00, &
   0.50580870785392503961D+00, &
   0.00000000000000000000D+00, &
  0.34700376603835188472D+00, &
  -0.50580870785392503961D+00, &
   0.50580870785392503961D+00, &
  0.34700376603835188472D+00, &
   0.00000000000000000000D+00, &
   0.50580870785392503961D+00, &
  0.34700376603835188472D+00, &
   0.50580870785392503961D+00, &
   0.50580870785392503961D+00, &
  0.34700376603835188472D+00, &
  -0.22850430565396735360D+00, &
  -0.22850430565396735360D+00, &
  0.70500220988849838312D+00, &
   0.00000000000000000000D+00, &
  -0.22850430565396735360D+00, &
  0.70500220988849838312D+00, &
   0.22850430565396735360D+00, &
  -0.22850430565396735360D+00, &
  0.70500220988849838312D+00, &
  -0.22850430565396735360D+00, &
   0.00000000000000000000D+00, &
  0.70500220988849838312D+00, &
   0.00000000000000000000D+00, &
   0.00000000000000000000D+00, &
  0.70500220988849838312D+00, &
   0.22850430565396735360D+00, &
   0.00000000000000000000D+00, &
  0.70500220988849838312D+00, &
  -0.22850430565396735360D+00, &
   0.22850430565396735360D+00, &
  0.70500220988849838312D+00, &
   0.00000000000000000000D+00, &
   0.22850430565396735360D+00, &
  0.70500220988849838312D+00, &
   0.22850430565396735360D+00, &
   0.22850430565396735360D+00, &
  0.70500220988849838312D+00 /), &
  (/ 3, 27 /) )

  w(1:order) = w_save(1:order)
  xyz(1:3,1:order) = xyz_save(1:3,1:order)
end

subroutine pyra_unit_o48 ( w, xyz )

!*****************************************************************************80
!
!! PYRA_UNIT_O48 returns a 48 point quadrature rule for the unit pyramid.
!
!  Discussion:
!
!    The integration region is:
!
!    - ( 1 - Z ) <= X <= 1 - Z
!    - ( 1 - Z ) <= Y <= 1 - Z
!              0 <= Z <= 1.
!
!    When Z is zero, the integration region is a square lying in the (X,Y)
!    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the
!    radius of the square diminishes, and when Z reaches 1, the square has
!    contracted to the single point (0,0,1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Output, double precision W(48), the weights.
!
!    Output, double precision XYZ(3,48), the abscissas.
!
  implicit none

  integer , parameter :: order = 48

  double precision w(order)
  double precision :: w_save(48) = (/ &
  2.01241939442682455D-002, &
  2.01241939442682455D-002, &
  2.01241939442682455D-002, &
  2.01241939442682455D-002, &
  2.60351137043010779D-002, &
  2.60351137043010779D-002, &
  2.60351137043010779D-002, &
  2.60351137043010779D-002, &
  1.24557795239745531D-002, &
  1.24557795239745531D-002, &
  1.24557795239745531D-002, &
  1.24557795239745531D-002, &
  1.87873998794808156D-003, &
  1.87873998794808156D-003, &
  1.87873998794808156D-003, &
  1.87873998794808156D-003, &
  4.32957927807745280D-002, &
  4.32957927807745280D-002, &
  4.32957927807745280D-002, &
  4.32957927807745280D-002, &
  1.97463249834127288D-002, &
  1.97463249834127288D-002, &
  1.97463249834127288D-002, &
  1.97463249834127288D-002, &
  5.60127223523590526D-002, &
  5.60127223523590526D-002, &
  5.60127223523590526D-002, &
  5.60127223523590526D-002, &
  2.55462562927473852D-002, &
  2.55462562927473852D-002, &
  2.55462562927473852D-002, &
  2.55462562927473852D-002, &
  2.67977366291788643D-002, &
  2.67977366291788643D-002, &
  2.67977366291788643D-002, &
  2.67977366291788643D-002, &
  1.22218992265373354D-002, &
  1.22218992265373354D-002, &
  1.22218992265373354D-002, &
  1.22218992265373354D-002, &
  4.04197740453215038D-003, &
  4.04197740453215038D-003, &
  4.04197740453215038D-003, &
  4.04197740453215038D-003, &
  1.84346316995826843D-003, &
  1.84346316995826843D-003, &
  1.84346316995826843D-003, &
  1.84346316995826843D-003 /)
  double precision xyz(3,order)
  double precision :: xyz_save(3,48) = reshape ( (/ &
  0.88091731624450909D+00, &
   0.0000000000000000D+00, &
   4.85005494469969989D-02, &
 -0.88091731624450909D+00, &
   0.0000000000000000D+00, &
   4.85005494469969989D-02, &
   0.0000000000000000D+00, &
   0.88091731624450909D+00, &
  4.85005494469969989D-02, &
   0.0000000000000000D+00, &
  -0.88091731624450909D+00, &
  4.85005494469969989D-02, &
  0.70491874112648223D+00, &
   0.0000000000000000D+00, &
   0.23860073755186201D+00, &
 -0.70491874112648223D+00, &
   0.0000000000000000D+00, &
   0.23860073755186201D+00, &
   0.0000000000000000D+00, &
   0.70491874112648223D+00, &
  0.23860073755186201D+00, &
   0.0000000000000000D+00, &
  -0.70491874112648223D+00, &
  0.23860073755186201D+00, &
  0.44712732143189760D+00, &
   0.0000000000000000D+00, &
   0.51704729510436798D+00, &
 -0.44712732143189760D+00, &
   0.0000000000000000D+00, &
   0.51704729510436798D+00, &
   0.0000000000000000D+00, &
   0.44712732143189760D+00, &
  0.51704729510436798D+00, &
   0.0000000000000000D+00, &
  -0.44712732143189760D+00, &
  0.51704729510436798D+00, &
  0.18900486065123448D+00, &
   0.0000000000000000D+00, &
   0.79585141789677305D+00, &
 -0.18900486065123448D+00, &
   0.0000000000000000D+00, &
   0.79585141789677305D+00, &
   0.0000000000000000D+00, &
   0.18900486065123448D+00, &
  0.79585141789677305D+00, &
   0.0000000000000000D+00, &
  -0.18900486065123448D+00, &
  0.79585141789677305D+00, &
  0.36209733410322176D+00, &
   0.36209733410322176D+00, &
  4.85005494469969989D-02, &
 -0.36209733410322176D+00, &
   0.36209733410322176D+00, &
  4.85005494469969989D-02, &
 -0.36209733410322176D+00, &
  -0.36209733410322176D+00, &
  4.85005494469969989D-02, &
  0.36209733410322176D+00, &
  -0.36209733410322176D+00, &
  4.85005494469969989D-02, &
  0.76688932060387538D+00, &
   0.76688932060387538D+00, &
  4.85005494469969989D-02, &
 -0.76688932060387538D+00, &
   0.76688932060387538D+00, &
  4.85005494469969989D-02, &
 -0.76688932060387538D+00, &
  -0.76688932060387538D+00, &
  4.85005494469969989D-02, &
  0.76688932060387538D+00, &
  -0.76688932060387538D+00, &
  4.85005494469969989D-02, &
  0.28975386476618070D+00, &
   0.28975386476618070D+00, &
  0.23860073755186201D+00, &
 -0.28975386476618070D+00, &
   0.28975386476618070D+00, &
  0.23860073755186201D+00, &
 -0.28975386476618070D+00, &
  -0.28975386476618070D+00, &
  0.23860073755186201D+00, &
  0.28975386476618070D+00, &
  -0.28975386476618070D+00, &
  0.23860073755186201D+00, &
  0.61367241226233160D+00, &
   0.61367241226233160D+00, &
  0.23860073755186201D+00, &
 -0.61367241226233160D+00, &
   0.61367241226233160D+00, &
  0.23860073755186201D+00, &
 -0.61367241226233160D+00, &
  -0.61367241226233160D+00, &
  0.23860073755186201D+00, &
  0.61367241226233160D+00, &
  -0.61367241226233160D+00, &
  0.23860073755186201D+00, &
  0.18378979287798017D+00, &
   0.18378979287798017D+00, &
  0.51704729510436798D+00, &
 -0.18378979287798017D+00, &
   0.18378979287798017D+00, &
  0.51704729510436798D+00, &
 -0.18378979287798017D+00, &
  -0.18378979287798017D+00, &
  0.51704729510436798D+00, &
  0.18378979287798017D+00, &
  -0.18378979287798017D+00, &
  0.51704729510436798D+00, &
  0.38925011625173161D+00, &
   0.38925011625173161D+00, &
  0.51704729510436798D+00, &
 -0.38925011625173161D+00, &
   0.38925011625173161D+00, &
  0.51704729510436798D+00, &
 -0.38925011625173161D+00, &
  -0.38925011625173161D+00, &
  0.51704729510436798D+00, &
  0.38925011625173161D+00, &
  -0.38925011625173161D+00, &
  0.51704729510436798D+00, &
  7.76896479525748113D-02, &
   7.76896479525748113D-02, &
  0.79585141789677305D+00, &
 -7.76896479525748113D-02, &
   7.76896479525748113D-02, &
  0.79585141789677305D+00, &
 -7.76896479525748113D-02, &
  -7.76896479525748113D-02, &
  0.79585141789677305D+00, &
  7.76896479525748113D-02, &
  -7.76896479525748113D-02, &
  0.79585141789677305D+00, &
  0.16453962988669860D+00, &
   0.16453962988669860D+00, &
  0.79585141789677305D+00, &
 -0.16453962988669860D+00, &
   0.16453962988669860D+00, &
  0.79585141789677305D+00, &
 -0.16453962988669860D+00, &
  -0.16453962988669860D+00, &
  0.79585141789677305D+00, &
  0.16453962988669860D+00, &
  -0.16453962988669860D+00, &
  0.79585141789677305D+00 /), &
  (/ 3, 48 /) )

  w(1:order) = w_save(1:order)
  xyz(1:3,1:order) = xyz_save(1:3,1:order)
end

subroutine pyra_unit_quad_test ( degree_max )

!*****************************************************************************80
!
!! PYRA_UNIT_QUAD_TEST tests the rules for the unit pyramid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DEGREE_MAX, the maximum total degree of the
!    monomials to check.
!
  implicit none

  integer , parameter :: dim_num = 3

  integer degree_max
  integer expon(dim_num)
  integer h
  logical more
  integer order
  double precision quad
  integer t
  double precision pyra_unit_volume
  double precision , allocatable :: v(:)
  double precision , allocatable :: w(:)
  double precision , allocatable :: xyz(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PYRA_UNIT_QUAD_TEST'
  write ( *, '(a)' ) '  For the unit pyramid,'
  write ( *, '(a)' ) '  we approximate monomial integrals with:'
  write ( *, '(a)' ) '  PYRA_UNIT_O01,'
  write ( *, '(a)' ) '  PYRA_UNIT_O05,'
  write ( *, '(a)' ) '  PYRA_UNIT_O06,'
  write ( *, '(a)' ) '  PYRA_UNIT_O08,'
  write ( *, '(a)' ) '  PYRA_UNIT_O08b,'
  write ( *, '(a)' ) '  PYRA_UNIT_O09,'
  write ( *, '(a)' ) '  PYRA_UNIT_O13,'
  write ( *, '(a)' ) '  PYRA_UNIT_O18,'
  write ( *, '(a)' ) '  PYRA_UNIT_O27,'
  write ( *, '(a)' ) '  PYRA_UNIT_O48.'

  more = .false.

  do

    call subcomp_next ( degree_max, dim_num, expon, more, h, t )

    if ( mod ( expon(1), 2 ) == 1 .or. &
         mod ( expon(2), 2 ) == 1 ) then
      cycle
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,2x,i2,2x,i2,2x,i2)' ) &
      '  Monomial exponents: ', expon(1:dim_num)
    write ( *, '(a)' ) ' '

    order = 1
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call pyra_unit_o01 ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = pyra_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 5
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call pyra_unit_o05 ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = pyra_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 6
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call pyra_unit_o06 ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = pyra_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 8
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call pyra_unit_o08 ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = pyra_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 8
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call pyra_unit_o08b ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = pyra_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 9
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call pyra_unit_o09 ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = pyra_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 13
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call pyra_unit_o13 ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = pyra_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 18
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call pyra_unit_o18 ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = pyra_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 27
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call pyra_unit_o27 ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = pyra_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 48
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call pyra_unit_o48 ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = pyra_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    write ( *, '(a)' ) ' '
    call pyra_unit_monomial ( expon, quad )
    write ( *, '(2x,a,2x,g14.6)' ) ' Exact', quad

    if ( .not. more ) then
      exit
    end if

  end do
end

function pyra_unit_volume ( )

!*****************************************************************************80
!
!! PYRA_UNIT_VOLUME: volume of a unit pyramid with square base.
!
!  Discussion:
!
!    The volume of this unit pyramid is 4/3.
!
!    The integration region is:
!
!      - ( 1 - Z ) <= X <= 1 - Z
!      - ( 1 - Z ) <= Y <= 1 - Z
!                0 <= Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, double precision PYRA_UNIT_VOLUME, the volume.
!
  implicit none

  double precision pyra_unit_volume

  pyra_unit_volume = 4.0D+00 / 3.0D+00
end

subroutine quad_unit_monomial ( expon, value )

!*****************************************************************************80
!
!! QUAD_UNIT_MONOMIAL integrates a monomial over the unit quadrilateral.
!
!  Discussion:
!
!    This routine integrates a monomial of the form
!
!      product ( 1 <= dim <= 2 ) x(dim)^expon(dim)
!
!    where the exponents are nonnegative integers.  Note that
!    if the combination 0^0 is encountered, it should be treated
!    as 1.
!
!    The integration region is:
!
!    - 1.0 <= X <= 1.0
!    - 1.0 <= Y <= 1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer EXPON(2), the exponents.
!
!    Output, double precision VALUE, the integral of the monomial.
!
  implicit none

  integer expon(2)
  integer i
  double precision value

  value = 1.0D+00

  do i = 1, 2

    if ( mod ( expon(i), 2 ) == 1 ) then
      value = 0.0D+00
    else if ( expon(i) == -1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'QUAD_UNIT_MONOMIAL - Fatal error!'
      write ( *, '(a)' ) '  Exponent of -1 encountered.'
      stop 1
    else
      value = value * 2.0D+00 / real ( expon(i) + 1)
    end if

  end do
end

subroutine quad_unit_monomial_test ( degree_max )

!*****************************************************************************80
!
!! QUAD_UNIT_MONOMIAL_TEST tests QUAD_UNIT_MONOMIAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DEGREE_MAX, the maximum total degree of the
!    monomials to check.
!
  implicit none

  integer alpha
  integer beta
  integer degree_max
  integer expon(2)
  double precision quad_unit_volume
  double precision value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUAD_UNIT_MONOMIAL_TEST'
  write ( *, '(a)' ) '  For the unit quadrilateral,'
  write ( *, '(a)' ) '  QUAD_UNIT_MONOMIAL returns the exact value of the'
  write ( *, '(a)' ) '  integral of X^ALPHA Y^BETA'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Volume = ', quad_unit_volume ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     ALPHA      BETA      INTEGRAL'
  write ( *, '(a)' ) ' '

  do alpha = 0, degree_max
    expon(1) = alpha
    do beta = 0, degree_max - alpha
      expon(2) = beta
      call quad_unit_monomial ( expon, value )
      write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) expon(1:2), value
    end do
  end do
end

subroutine quad_unit_quad_test ( degree_max )

!*****************************************************************************80
!
!! QUAD_UNIT_QUAD_TEST tests the rules for the unit quadrilateral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DEGREE_MAX, the maximum total degree of the
!    monomials to check.
!
  implicit none

  integer , parameter :: dim_num = 2

  integer degree_max
  integer expon(dim_num)
  integer h
  integer k
  logical more
  integer order
  integer order_1d(dim_num)
  double precision quad
  double precision quad_unit_volume
  integer t
  double precision , allocatable :: v(:)
  double precision , allocatable :: w(:)
  double precision , allocatable :: xy(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUAD_UNIT_QUAD_TEST'
  write ( *, '(a)' ) '  For the unit quadrilateral,'
  write ( *, '(a)' ) '  we approximate monomial integrals with'
  write ( *, '(a)' ) '  QUAD_UNIT_RULE, which returns M by N point rules.'

  more = .false.

  do

    call subcomp_next ( degree_max, dim_num, expon, more, h, t )

    if ( any ( mod ( expon(1:dim_num), 2 ) == 1 ) ) then
      if ( .not. more ) then
        exit
      else
        cycle
      end if
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,2x,i2,2x,i2)' ) '  Monomial exponents: ', expon(1:dim_num)
    write ( *, '(a)' ) ' '

    do k = 1, 5

      order_1d(1:dim_num) = k
      order = product ( order_1d(1:dim_num) )
      allocate ( v(1:order) )
      allocate ( w(1:order) )
      allocate ( xy(1:dim_num,1:order) )
      call quad_unit_rule ( order_1d, w, xy )
      call monomial_value ( dim_num, order, expon, xy, v )
      quad = quad_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
      write ( *, '(2x,i6,2x,i6,2x,g14.6)' ) order_1d(1:dim_num), quad
      deallocate ( v )
      deallocate ( w )
      deallocate ( xy )

    end do
!
!  Try a rule of mixed orders.
!
    order_1d(1) = 3
    order_1d(2) = 5
    order = product ( order_1d(1:dim_num) )
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xy(1:dim_num,1:order) )
    call quad_unit_rule ( order_1d, w, xy )
    call monomial_value ( dim_num, order, expon, xy, v )
    quad = quad_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,i6,2x,g14.6)' ) order_1d(1:dim_num), quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xy )

    write ( *, '(a)' ) ' '
    call quad_unit_monomial ( expon, quad )
    write ( *, '(2x,a,2x,6x,2x,g14.6)' ) ' Exact', quad

    if ( .not. more ) then
      exit
    end if

  end do
end

subroutine quad_unit_rule ( order_1d, w, xy )

!*****************************************************************************80
!
!! QUAD_UNIT_RULE returns a quadrature rule for the unit quadrilateral.
!
!  Discussion:
!
!    The integration region is:
!
!    - 1.0 <= X <= 1.0
!    - 1.0 <= Y <= 1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Input, integer ORDER_1D(2), the order of the rule in
!    each dimension.  1 <= ORDER_1D(I) <= 5.
!
!    Output, double precision W(ORDER_1D(1)*ORDER_1D(2)), the weights.
!
!    Output, double precision XY(2,ORDER_1D(1)*ORDER_1D(2)), the abscissas.
!
  implicit none

  integer , parameter :: dim_num = 2

  integer dim
  integer order
  integer order_1d(2)
  double precision w(order_1D(1)*order_1D(2))
  double precision , allocatable :: w_1d(:)
  double precision , allocatable :: x_1d(:)
  double precision xy(2,order_1D(1)*order_1D(2))

  order = product ( order_1d(1:dim_num) )

  do dim = 1, dim_num

    allocate ( w_1d(order_1D(dim)) )
    allocate ( x_1d(order_1D(dim)) )

    if ( order_1D(dim) == 1 ) then
      call line_unit_o01 ( w_1d, x_1d )
    else if ( order_1D(dim) == 2 ) then
      call line_unit_o02 ( w_1d, x_1d )
    else if ( order_1D(dim) == 3 ) then
      call line_unit_o03 ( w_1d, x_1d )
    else if ( order_1D(dim) == 4 ) then
      call line_unit_o04 ( w_1d, x_1d )
    else if ( order_1D(dim) == 5 ) then
      call line_unit_o05 ( w_1d, x_1d )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'QUAD_UNIT_RULE - Fatal error!'
      write ( *, '(a)' ) '  Illegal value of ORDER_1D(*).'
      stop 1
    end if

    call r8vec_direct_product ( dim, order_1D(dim), x_1d, &
      dim_num, order, xy )

    call r8vec_direct_product2 ( dim, order_1D(dim), w_1d, &
      dim_num, order, w )

    deallocate ( w_1d )
    deallocate ( x_1d )

  end do
end

function quad_unit_volume ( )

!*****************************************************************************80
!
!! QUAD_UNIT_VOLUME: volume of a unit quadrilateral.
!
!  Discussion:
!
!    The integration region is:
!
!    - 1.0 <= X <= 1.0
!    - 1.0 <= Y <= 1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, double precision QUAD_UNIT_VOLUME, the volume.
!
  implicit none

  double precision quad_unit_volume

  quad_unit_volume = 4.0D+00
end

subroutine subcomp_next ( n, k, a, more, h, t )

!*****************************************************************************80
!
!! SUBCOMP_NEXT computes the next subcomposition of N into K parts.
!
!  Discussion:
!
!    A composition of the integer N into K parts is an ordered sequence
!    of K nonnegative integers which sum to a value of N.
!
!    A subcomposition of the integer N into K parts is a composition
!    of M into K parts, where 0 <= M <= N.
!
!    A subcomposition of the integer N into K parts is also a lattice
!    point in the simplex whose vertices are the origin, and the K direction
!    vectors N*E(I) for I = 1 to K.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the integer whose subcompositions
!    are desired.
!
!    Input, integer K, the number of parts in the subcomposition.
!
!    Input/output, integer A(K), the parts of the subcomposition.
!
!    Input/output, logical MORE, set by the user to start the computation,
!    and by the routine to terminate it.
!
!    Input/output, integer H, T, two internal parameters needed for the
!    computation.  The user should allocate space for these in the calling
!    program, include them in the calling sequence, but never alter them!
!
  implicit none

  integer k

  integer a(k)
  integer h
  logical more
  logical, save :: more2 = .false.
  integer n
  integer , save :: n2 = 0
  integer t
!
!  The first computation.
!
  if ( .not. more ) then

    n2 = 0
    a(1:k) = 0
    more2 = .false.
    h = 0
    t = 0

    more = .true.
!
!  Do the next element at the current value of N.
!
  else if ( more2 ) then

    call comp_next ( n2, k, a, more2, h, t )

  else

    more2 = .false.
    n2 = n2 + 1

    call comp_next ( n2, k, a, more2, h, t )

  end if
!
!  Termination occurs if MORE2 = FALSE and N2 = N.
!
  if ( .not. more2 .and. n2 == n ) then
    more = .false.
  end if
end

subroutine tetr_unit_monomial ( expon, value )

!*****************************************************************************80
!
!! TETR_UNIT_MONOMIAL integrates a monomial over the unit tetrahedron.
!
!  Discussion:
!
!    This routine integrates a monomial of the form
!
!      product ( 1 <= dim <= 3 ) x(dim)^expon(dim)
!
!    where the exponents are nonnegative integers.  Note that
!    if the combination 0^0 is encountered, it should be treated
!    as 1.
!
!    Integral ( over unit tetrahedron ) x^l y^m z^n dx dy =
!    l! * m! * n! / ( m + n + 3 )!
!
!    The integration region is:
!
!      0 <= X
!      0 <= Y
!      0 <= Z
!      0 <= X + Y + Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer EXPON(3), the exponents.
!
!    Output, double precision VALUE, the integral of the monomial.
!
  implicit none

  integer expon(3)
  integer i
  integer k
  double precision value
!
!  The first computation ends with VALUE = 1.0;
!
  value = 1.0D+00
!
!  The first loop simply calculates 1, so we short circuit it.
!
! k = 0
!
! do i = 1, expon(1)
!   k = k + 1
!   value = value * real ( i) / real ( k)
! end do

  k = expon(1)
  do i = 1, expon(2)
    k = k + 1
    value = value * real ( i) / real ( k)
  end do

  do i = 1, expon(3)
    k = k + 1
    value = value * real ( i) / real ( k)
  end do

  k = k + 1
  value = value / real ( k)

  k = k + 1
  value = value / real ( k)

  k = k + 1
  value = value / real ( k)
end

subroutine tetr_unit_monomial_test ( degree_max )

!*****************************************************************************80
!
!! TETR_UNIT_MONOMIAL_TEST tests TETR_UNIT_MONOMIAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DEGREE_MAX, the maximum total degree of the
!    monomials to check.
!
  implicit none

  integer alpha
  integer beta
  integer degree_max
  integer expon(3)
  integer gamma
  double precision tetr_unit_volume
  double precision value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TETR_UNIT_MONOMIAL_TEST'
  write ( *, '(a)' ) '  For the unit tetrahedron,'
  write ( *, '(a)' ) '  TETR_UNIT_MONOMIAL returns the exact value of the'
  write ( *, '(a)' ) '  integral of X^ALPHA Y^BETA Z^GAMMA'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Volume = ', tetr_unit_volume ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     ALPHA      BETA     GAMMA      INTEGRAL'
  write ( *, '(a)' ) ' '

  do alpha = 0, degree_max
    expon(1) = alpha
    do beta = 0, degree_max - alpha
      expon(2) = beta
      do gamma = 0, degree_max - alpha - beta
        expon(3) = gamma
        call tetr_unit_monomial ( expon, value )
        write ( *, '(2x,i8,2x,i8,2x,i8,2x,g14.6)' ) expon(1:3), value
      end do
    end do
  end do
end

subroutine tetr_unit_o01 ( w, xyz )

!*****************************************************************************80
!
!! TETR_UNIT_O01 returns a 1 point quadrature rule for the unit tetrahedron.
!
!  Discussion:
!
!    The integration region is:
!
!      0 <= X
!      0 <= Y
!      0 <= Z
!      X + Y + Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(1), the weights.
!
!    Output, double precision XYZ(3,1), the abscissas.
!
  implicit none

  integer , parameter :: order = 1

  double precision w(order)
  double precision :: w_save(1) = (/ &
    1.0000000000000000000D+00 /)
  double precision xyz(3,order)
  double precision :: xyz_save(3,1) = reshape ( (/ &
    0.25000000000000000000D+00,  0.25000000000000000000D+00,  &
    0.25000000000000000000D+00 /), &
  (/ 3, 1 /) )

  w(1:order) = w_save(1:order)
  xyz(1:3,1:order) = xyz_save(1:3,1:order)
end

subroutine tetr_unit_o04 ( w, xyz )

!*****************************************************************************80
!
!! TETR_UNIT_O04 returns a 4 point quadrature rule for the unit tetrahedron.
!
!  Discussion:
!
!    The integration region is:
!
!      0 <= X
!      0 <= Y
!      0 <= Z
!      X + Y + Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(4), the weights.
!
!    Output, double precision XYZ(3,4), the abscissas.
!
  implicit none

  integer , parameter :: order = 4

  double precision w(order)
  double precision :: w_save(4) = (/ &
    0.25000000000000000000D+00, &
    0.25000000000000000000D+00, &
    0.25000000000000000000D+00, &
    0.25000000000000000000D+00 /)
  double precision xyz(3,order)
  double precision :: xyz_save(3,4) = reshape ( (/ &
    0.58541019662496845446D+00, &
    0.13819660112501051518D+00, &
    0.13819660112501051518D+00, &
    0.13819660112501051518D+00, &
    0.58541019662496845446D+00, &
    0.13819660112501051518D+00, &
    0.13819660112501051518D+00, &
    0.13819660112501051518D+00, &
    0.58541019662496845446D+00, &
    0.13819660112501051518D+00, &
    0.13819660112501051518D+00, &
    0.13819660112501051518D+00 /), (/ 3, 4 /) )

  w(1:order) = w_save(1:order)
  xyz(1:3,1:order) = xyz_save(1:3,1:order)
end

subroutine tetr_unit_o08 ( w, xyz )

!*****************************************************************************80
!
!! TETR_UNIT_O08 returns an 8 point quadrature rule for the unit tetrahedron.
!
!  Discussion:
!
!    The integration region is:
!
!      0 <= X
!      0 <= Y
!      0 <= Z
!      X + Y + Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(8), the weights.
!
!    Output, double precision XYZ(3,8), the abscissas.
!
  implicit none

  integer , parameter :: order = 8

  double precision w(order)
  double precision :: w_save(8) = (/ &
    0.13852796651186214232D+00, &
    0.13852796651186214232D+00, &
    0.13852796651186214232D+00, &
    0.13852796651186214232D+00, &
    0.11147203348813785768D+00, &
    0.11147203348813785768D+00, &
    0.11147203348813785768D+00, &
    0.11147203348813785768D+00 /)
  double precision xyz(3,order)
  double precision :: xyz_save(3,8) = reshape ( (/ &
    0.015835909865720057993D+00, &
    0.32805469671142664734D+00,  &
    0.32805469671142664734D+00,  &
    0.32805469671142664734D+00,  &
    0.015835909865720057993D+00, &
    0.32805469671142664734D+00,  &
    0.32805469671142664734D+00,  &
    0.32805469671142664734D+00,  &
    0.015835909865720057993D+00, &
    0.32805469671142664734D+00,  &
    0.32805469671142664734D+00,  &
    0.32805469671142664734D+00,  &
    0.67914317820120795168D+00,  &
    0.10695227393293068277D+00,  &
    0.10695227393293068277D+00,  &
    0.10695227393293068277D+00,  &
    0.67914317820120795168D+00,  &
    0.10695227393293068277D+00,  &
    0.10695227393293068277D+00,  &
    0.10695227393293068277D+00,  &
    0.67914317820120795168D+00,  &
    0.10695227393293068277D+00,  &
    0.10695227393293068277D+00,  &
    0.10695227393293068277D+00 /), &
  (/ 3, 8 /) )

  w(1:order) = w_save(1:order)
  xyz(1:3,1:order) = xyz_save(1:3,1:order)
end

subroutine tetr_unit_o08b ( w, xyz )

!*****************************************************************************80
!
!! TETR_UNIT_O08B returns an 8 point quadrature rule for the unit tetrahedron.
!
!  Discussion:
!
!    The integration region is:
!
!      0 <= X
!      0 <= Y
!      0 <= Z
!      X + Y + Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(8), the weights.
!
!    Output, double precision XYZ(3,8), the abscissas.
!
  implicit none

  integer , parameter :: order = 8

  double precision w(order)
  double precision :: w_save(8) = (/ &
    0.025000000000000000000D+00, &
    0.025000000000000000000D+00, &
    0.025000000000000000000D+00, &
    0.025000000000000000000D+00, &
    0.22500000000000000000D+00, &
    0.22500000000000000000D+00, &
    0.22500000000000000000D+00, &
    0.22500000000000000000D+00 /)
  double precision xyz(3,order)
  double precision :: xyz_save(3,8) = reshape ( (/ &
    1.00000000000000000000D+00, &
    0.00000000000000000000D+00, &
    0.00000000000000000000D+00, &
    0.00000000000000000000D+00, &
    1.00000000000000000000D+00, &
    0.00000000000000000000D+00, &
    0.00000000000000000000D+00, &
    0.00000000000000000000D+00, &
    1.00000000000000000000D+00, &
    0.00000000000000000000D+00, &
    0.00000000000000000000D+00, &
    0.00000000000000000000D+00, &
    0.00000000000000000000D+00, &
    0.33333333333333333333D+00, &
    0.33333333333333333333D+00, &
    0.33333333333333333333D+00, &
    0.00000000000000000000D+00, &
    0.33333333333333333333D+00, &
    0.33333333333333333333D+00, &
    0.33333333333333333333D+00, &
    0.00000000000000000000D+00, &
    0.33333333333333333333D+00, &
    0.33333333333333333333D+00, &
    0.33333333333333333333D+00 /), &
  (/ 3, 8 /) )

  w(1:order) = w_save(1:order)
  xyz(1:3,1:order) = xyz_save(1:3,1:order)
end

subroutine tetr_unit_o14 ( w, xyz )

!*****************************************************************************80
!
!! TETR_UNIT_O14 returns a 14 point quadrature rule for the unit tetrahedron.
!
!  Discussion:
!
!    The integration region is:
!
!      0 <= X
!      0 <= Y
!      0 <= Z
!      X + Y + Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(ORDER), the weights.
!
!    Output, double precision XYZ(3,ORDER), the abscissas.
!
  implicit none

  integer , parameter :: order = 14

  double precision w(order)
  double precision :: w_save(14) = (/ &
    0.073493043116361949544D+00, &
    0.073493043116361949544D+00, &
    0.073493043116361949544D+00, &
    0.073493043116361949544D+00, &
    0.11268792571801585080D+00, &
    0.11268792571801585080D+00, &
    0.11268792571801585080D+00, &
    0.11268792571801585080D+00, &
    0.042546020777081466438D+00, &
    0.042546020777081466438D+00, &
    0.042546020777081466438D+00, &
    0.042546020777081466438D+00, &
    0.042546020777081466438D+00, &
    0.042546020777081466438D+00 /)
  double precision xyz(3,order)
  double precision :: xyz_save(3,14) = reshape ( (/ &
    0.72179424906732632079D+00, &
    0.092735250310891226402D+00, &
    0.092735250310891226402D+00, &
    0.092735250310891226402D+00, &
    0.72179424906732632079D+00, &
    0.092735250310891226402D+00, &
    0.092735250310891226402D+00, &
    0.092735250310891226402D+00, &
    0.72179424906732632079D+00, &
    0.092735250310891226402D+00, &
    0.092735250310891226402D+00, &
    0.092735250310891226402D+00, &
    0.067342242210098170608D+00, &
    0.31088591926330060980D+00, &
    0.31088591926330060980D+00, &
    0.31088591926330060980D+00, &
    0.067342242210098170608D+00, &
    0.31088591926330060980D+00, &
    0.31088591926330060980D+00, &
    0.31088591926330060980D+00, &
    0.067342242210098170608D+00, &
    0.31088591926330060980D+00, &
    0.31088591926330060980D+00, &
    0.31088591926330060980D+00, &
    0.045503704125649649492D+00, &
    0.045503704125649649492D+00, &
    0.45449629587435035051D+00, &
    0.045503704125649649492D+00, &
    0.45449629587435035051D+00, &
    0.045503704125649649492D+00, &
    0.045503704125649649492D+00, &
    0.45449629587435035051D+00, &
    0.45449629587435035051D+00, &
    0.45449629587435035051D+00, &
    0.045503704125649649492D+00, &
    0.045503704125649649492D+00, &
    0.45449629587435035051D+00, &
    0.045503704125649649492D+00, &
    0.45449629587435035051D+00, &
    0.45449629587435035051D+00, &
    0.45449629587435035051D+00, &
    0.045503704125649649492D+00 /), &
  (/ 3, 14 /) )

  w(1:order) = w_save(1:order)
  xyz(1:3,1:order) = xyz_save(1:3,1:order)
end

subroutine tetr_unit_o14b ( w, xyz )

!*****************************************************************************80
!
!! TETR_UNIT_O14B returns a 14 point quadrature rule for the unit tetrahedron.
!
!  Discussion:
!
!    The integration region is:
!
!      0 <= X
!      0 <= Y
!      0 <= Z
!      X + Y + Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(14), the weights.
!
!    Output, double precision XYZ(3,14), the abscissas.
!
  implicit none

  integer , parameter :: order = 14

  double precision w(order)
  double precision :: w_save(14) = (/ &
    0.13283874668559071814D+00, &
    0.13283874668559071814D+00, &
    0.13283874668559071814D+00, &
    0.13283874668559071814D+00, &
    0.088589824742980710434D+00, &
    0.088589824742980710434D+00, &
    0.088589824742980710434D+00, &
    0.088589824742980710434D+00, &
    0.019047619047619047619D+00, &
    0.019047619047619047619D+00, &
    0.019047619047619047619D+00, &
    0.019047619047619047619D+00, &
    0.019047619047619047619D+00, &
    0.019047619047619047619D+00  /)
  double precision xyz(3,order) 
  double precision :: xyz_save(3,14) = reshape ( (/ &
    0.056881379520423421748D+00, &
    0.31437287349319219275D+00, &
    0.31437287349319219275D+00, &
    0.31437287349319219275D+00, &
    0.056881379520423421748D+00, &
    0.31437287349319219275D+00, &
    0.31437287349319219275D+00, &
    0.31437287349319219275D+00, &
    0.056881379520423421748D+00, &
    0.31437287349319219275D+00, &
    0.31437287349319219275D+00, &
    0.31437287349319219275D+00, &
    0.69841970432438656092D+00, &
    0.10052676522520447969D+00, &
    0.10052676522520447969D+00, &
    0.10052676522520447969D+00, &
    0.69841970432438656092D+00, &
    0.10052676522520447969D+00, &
    0.10052676522520447969D+00, &
    0.10052676522520447969D+00, &
    0.69841970432438656092D+00, &
    0.10052676522520447969D+00, &
    0.10052676522520447969D+00, &
    0.10052676522520447969D+00, &
    0.50000000000000000000D+00, &
    0.50000000000000000000D+00, &
    0.00000000000000000000D+00, &
    0.50000000000000000000D+00, &
    0.00000000000000000000D+00, &
    0.50000000000000000000D+00, &
    0.50000000000000000000D+00, &
    0.00000000000000000000D+00, &
    0.00000000000000000000D+00, &
    0.00000000000000000000D+00, &
    0.50000000000000000000D+00, &
    0.50000000000000000000D+00, &
    0.00000000000000000000D+00, &
    0.50000000000000000000D+00, &
    0.00000000000000000000D+00, &
    0.00000000000000000000D+00, &
    0.00000000000000000000D+00, &
    0.50000000000000000000D+00 /), &
  (/ 3, 14 /) )

  w(1:order) = w_save(1:order)
  xyz(1:3,1:order) = xyz_save(1:3,1:order)
end

subroutine tetr_unit_o15 ( w, xyz )

!*****************************************************************************80
!
!! TETR_UNIT_O15 returns a 15 point quadrature rule for the unit tetrahedron.
!
!  Discussion:
!
!    The integration region is:
!
!      0 <= X
!      0 <= Y
!      0 <= Z
!      X + Y + Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(15), the weights.
!
!    Output, double precision XYZ(3,15), the abscissas.
!
  implicit none

  integer , parameter :: order = 15

  double precision w(order)
  double precision :: w_save(15) = (/ &
    0.071937083779018620010D+00, &
    0.071937083779018620010D+00, &
    0.071937083779018620010D+00, &
    0.071937083779018620010D+00, &
    0.069068207226272385281D+00, &
    0.069068207226272385281D+00, &
    0.069068207226272385281D+00, &
    0.069068207226272385281D+00, &
    0.052910052910052910053D+00, &
    0.052910052910052910053D+00, &
    0.052910052910052910053D+00, &
    0.052910052910052910053D+00, &
    0.052910052910052910053D+00, &
    0.052910052910052910053D+00, &
    0.11851851851851851852D+00 /)
  double precision xyz(3,order)
  double precision :: xyz_save(3,15) = reshape ( (/ &
    0.72408676584183090163D+00, &
  0.091971078052723032789D+00, &
  0.091971078052723032789D+00, &
    0.091971078052723032789D+00, &
  0.72408676584183090163D+00, &
  0.091971078052723032789D+00, &
    0.091971078052723032789D+00, &
  0.091971078052723032789D+00, &
  0.72408676584183090163D+00, &
    0.091971078052723032789D+00, &
  0.091971078052723032789D+00, &
  0.091971078052723032789D+00, &
    0.040619116511110274837D+00, &
  0.31979362782962990839D+00, &
  0.31979362782962990839D+00, &
    0.31979362782962990839D+00, &
  0.040619116511110274837D+00, &
  0.31979362782962990839D+00, &
    0.31979362782962990839D+00, &
  0.31979362782962990839D+00, &
  0.040619116511110274837D+00, &
    0.31979362782962990839D+00, &
  0.31979362782962990839D+00, &
  0.31979362782962990839D+00, &
    0.44364916731037084426D+00, &
  0.44364916731037084426D+00, &
  0.056350832689629155741D+00, &
    0.44364916731037084426D+00, &
  0.056350832689629155741D+00, &
  0.44364916731037084426D+00, &
    0.44364916731037084426D+00, &
  0.056350832689629155741D+00, &
  0.056350832689629155741D+00, &
    0.056350832689629155741D+00, &
  0.44364916731037084426D+00, &
  0.44364916731037084426D+00, &
    0.056350832689629155741D+00, &
  0.44364916731037084426D+00, &
  0.056350832689629155741D+00, &
    0.056350832689629155741D+00, &
  0.056350832689629155741D+00, &
  0.44364916731037084426D+00, &
    0.25000000000000000000D+00, &
  0.25000000000000000000D+00, &
  0.25000000000000000000D+00 /), &
  (/ 3, 15 /) )

  w(1:order) = w_save(1:order)
  xyz(1:3,1:order) = xyz_save(1:3,1:order)
end

subroutine tetr_unit_o15b ( w, xyz )

!*****************************************************************************80
!
!! TETR_UNIT_O15B returns a 15 point quadrature rule for the unit tetrahedron.
!
!  Discussion:
!
!    The integration region is:
!
!      0 <= X
!      0 <= Y
!      0 <= Z
!      X + Y + Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(15), the weights.
!
!    Output, double precision XYZ(3,15), the abscissas.
!
  implicit none

  integer , parameter :: order = 15

  double precision w(order)
  double precision :: w_save(15) = (/ &
    0.036160714285714285714D+00, &
    0.036160714285714285714D+00, &
    0.036160714285714285714D+00, &
    0.036160714285714285714D+00, &
    0.069871494516173816465D+00, &
    0.069871494516173816465D+00, &
    0.069871494516173816465D+00, &
    0.069871494516173816465D+00, &
    0.065694849368318756074D+00, &
    0.065694849368318756074D+00, &
    0.065694849368318756074D+00, &
    0.065694849368318756074D+00, &
    0.065694849368318756074D+00, &
    0.065694849368318756074D+00, &
    0.18170206858253505484D+00 /)
  double precision xyz(3,order)
  double precision :: xyz_save(3,15) = reshape ( (/ &
    0.00000000000000000000D+00, &
  0.33333333333333333333D+00, &
  0.33333333333333333333D+00, &
    0.33333333333333333333D+00, &
  0.00000000000000000000D+00, &
  0.33333333333333333333D+00, &
    0.33333333333333333333D+00, &
  0.33333333333333333333D+00, &
  0.00000000000000000000D+00, &
    0.33333333333333333333D+00, &
  0.33333333333333333333D+00, &
  0.33333333333333333333D+00, &
    0.72727272727272727273D+00, &
  0.090909090909090909091D+00, &
  0.090909090909090909091D+00, &
    0.090909090909090909091D+00, &
  0.72727272727272727273D+00, &
  0.090909090909090909091D+00, &
    0.090909090909090909091D+00, &
  0.090909090909090909091D+00, &
  0.72727272727272727273D+00, &
    0.090909090909090909091D+00, &
  0.090909090909090909091D+00, &
  0.090909090909090909091D+00, &
    0.43344984642633570176D+00, &
  0.43344984642633570176D+00, &
  0.066550153573664298240D+00, &
    0.43344984642633570176D+00, &
  0.066550153573664298240D+00, &
  0.43344984642633570176D+00, &
    0.43344984642633570176D+00, &
  0.066550153573664298240D+00, &
  0.066550153573664298240D+00, &
    0.066550153573664298240D+00, &
  0.43344984642633570176D+00, &
  0.43344984642633570176D+00, &
    0.066550153573664298240D+00, &
  0.43344984642633570176D+00, &
  0.066550153573664298240D+00, &
    0.066550153573664298240D+00, &
  0.066550153573664298240D+00, &
  0.43344984642633570176D+00, &
    0.25000000000000000000D+00, &
  0.25000000000000000000D+00, &
  0.250000000000000000D+00 /), &
  (/ 3, 15 /) )

  w(1:order) = w_save(1:order)
  xyz(1:3,1:order) = xyz_save(1:3,1:order)
end

subroutine tetr_unit_o24 ( w, xyz )

!*****************************************************************************80
!
!! TETR_UNIT_O24 returns a 24 point quadrature rule for the unit tetrahedron.
!
!  Discussion:
!
!    The integration region is:
!
!      0 <= X
!      0 <= Y
!      0 <= Z
!      X + Y + Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(24), the weights.
!
!    Output, double precision XYZ(3,24), the abscissas.
!
  implicit none

  integer , parameter :: order = 24

  double precision w(order)
  double precision :: w_save(24) = (/ &
    0.039922750257869636194D+00, &
    0.039922750257869636194D+00, &
    0.039922750257869636194D+00, &
    0.039922750257869636194D+00, &
    0.010077211055345822612D+00, &
    0.010077211055345822612D+00, &
    0.010077211055345822612D+00, &
    0.010077211055345822612D+00, &
    0.055357181543927398338D+00, &
    0.055357181543927398338D+00, &
    0.055357181543927398338D+00, &
    0.055357181543927398338D+00, &
    0.048214285714285714286D+00, &
    0.048214285714285714286D+00, &
    0.048214285714285714286D+00, &
    0.048214285714285714286D+00, &
    0.048214285714285714286D+00, &
    0.048214285714285714286D+00, &
    0.048214285714285714286D+00, &
    0.048214285714285714286D+00, &
    0.048214285714285714286D+00, &
    0.048214285714285714286D+00, &
    0.048214285714285714286D+00, &
    0.048214285714285714286D+00 /)
  double precision xyz(3,order)
  double precision :: xyz_save(3,24) = reshape ( (/ &
    0.35619138622025439121D+00, &
    0.21460287125991520293D+00, &
    0.21460287125991520293D+00, &
    0.21460287125991520293D+00, &
    0.35619138622025439121D+00, &
    0.21460287125991520293D+00, &
    0.21460287125991520293D+00, &
    0.21460287125991520293D+00, &
    0.35619138622025439121D+00, &
    0.21460287125991520293D+00, &
    0.21460287125991520293D+00, &
    0.21460287125991520293D+00, &
    0.87797812439616594065D+00, &
    0.040673958534611353116D+00, &
    0.040673958534611353116D+00, &
    0.040673958534611353116D+00, &
    0.87797812439616594065D+00, &
    0.040673958534611353116D+00, &
    0.040673958534611353116D+00, &
    0.040673958534611353116D+00, &
    0.87797812439616594065D+00, &
    0.040673958534611353116D+00, &
    0.040673958534611353116D+00, &
    0.040673958534611353116D+00, &
    0.032986329573173468968D+00, &
    0.32233789014227551034D+00, &
    0.32233789014227551034D+00, &
    0.32233789014227551034D+00, &
    0.032986329573173468968D+00, &
    0.32233789014227551034D+00, &
    0.32233789014227551034D+00, &
    0.32233789014227551034D+00, &
    0.032986329573173468968D+00, &
    0.32233789014227551034D+00, &
    0.32233789014227551034D+00, &
    0.32233789014227551034D+00, &
    0.60300566479164914137D+00, &
    0.26967233145831580803D+00, &
    0.063661001875017525299D+00, &
    0.60300566479164914137D+00, &
    0.063661001875017525299D+00, &
    0.26967233145831580803D+00, &
    0.60300566479164914137D+00, &
    0.063661001875017525299D+00, &
    0.063661001875017525299D+00, &
    0.063661001875017525299D+00, &
    0.60300566479164914137D+00, &
    0.26967233145831580803D+00, &
    0.063661001875017525299D+00, &
    0.60300566479164914137D+00, &
    0.063661001875017525299D+00, &
    0.063661001875017525299D+00, &
    0.063661001875017525299D+00, &
    0.60300566479164914137D+00, &
    0.26967233145831580803D+00, &
    0.60300566479164914137D+00, &
    0.063661001875017525299D+00, &
    0.26967233145831580803D+00, &
    0.063661001875017525299D+00, &
    0.60300566479164914137D+00, &
    0.26967233145831580803D+00, &
    0.063661001875017525299D+00, &
    0.063661001875017525299D+00, &
    0.063661001875017525299D+00, &
    0.26967233145831580803D+00, &
    0.60300566479164914137D+00, &
    0.063661001875017525299D+00, &
    0.26967233145831580803D+00, &
    0.063661001875017525299D+00, &
    0.063661001875017525299D+00, &
    0.063661001875017525299D+00, &
    0.26967233145831580803D+00 /), &
  (/ 3, 24 /) )

  w(1:order) = w_save(1:order)
  xyz(1:3,1:order) = xyz_save(1:3,1:order)
end

subroutine tetr_unit_quad_test ( degree_max )

!*****************************************************************************80
!
!! TETR_UNIT_QUAD_TEST tests the rules for the unit tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DEGREE_MAX, the maximum total degree of the
!    monomials to check.
!
  implicit none

  integer , parameter :: dim_num = 3

  integer degree_max
  integer expon(dim_num)
  integer h
  logical more
  integer order
  double precision quad
  integer t
  double precision tetr_unit_volume
  double precision , allocatable :: v(:)
  double precision , allocatable :: w(:)
  double precision , allocatable :: xyz(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TETR_UNIT_QUAD_TEST'
  write ( *, '(a)' ) '  For the unit tetrahedron,'
  write ( *, '(a)' ) '  we approximate monomial integrals with:'
  write ( *, '(a)' ) '  TETR_UNIT_O01,'
  write ( *, '(a)' ) '  TETR_UNIT_O04,'
  write ( *, '(a)' ) '  TETR_UNIT_O08,'
  write ( *, '(a)' ) '  TETR_UNIT_O08b,'
  write ( *, '(a)' ) '  TETR_UNIT_O14,'
  write ( *, '(a)' ) '  TETR_UNIT_O14b,'
  write ( *, '(a)' ) '  TETR_UNIT_O15,'
  write ( *, '(a)' ) '  TETR_UNIT_O15b,'
  write ( *, '(a)' ) '  TETR_UNIT_O24.'

  more = .false.

  do

    call subcomp_next ( degree_max, dim_num, expon, more, h, t )

    write ( *, '(a)' ) ' '
    write ( *, '(a,2x,i2,2x,i2,2x,i2)' ) &
      '  Monomial exponents: ', expon(1:dim_num)
    write ( *, '(a)' ) ' '

    order = 1
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call tetr_unit_o01 ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = tetr_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 4
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call tetr_unit_o04 ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = tetr_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 8
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call tetr_unit_o08 ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = tetr_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 8
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call tetr_unit_o08b ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = tetr_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 14
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call tetr_unit_o14 ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = tetr_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 14
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call tetr_unit_o14b ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = tetr_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 15
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call tetr_unit_o15 ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = tetr_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 15
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call tetr_unit_o15b ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = tetr_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 24
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call tetr_unit_o24 ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = tetr_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    write ( *, '(a)' ) ' '
    call tetr_unit_monomial ( expon, quad )
    write ( *, '(2x,a,2x,g14.6)' ) ' Exact', quad

    if ( .not. more ) then
      exit
    end if

  end do
end

function tetr_unit_volume ( )

!*****************************************************************************80
!
!! TETR_UNIT_VOLUME returns the volume of the unit tetrahedron.
!
!  Discussion:
!
!    The integration region is:
!
!      0 <= X,
!      0 <= Y,
!      0 <= Z,
!      X + Y + Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, double precision TETR_UNIT_VOLUME, the volume.
!
  implicit none

  double precision tetr_unit_volume

  tetr_unit_volume = 1.0D+00 / 6.0D+00
end

subroutine trig_unit_monomial ( expon, value )

!*****************************************************************************80
!
!! TRIG_UNIT_MONOMIAL integrates a monomial over the unit triangle.
!
!  Discussion:
!
!    This routine integrates a monomial of the form
!
!      product ( 1 <= dim <= 2 ) x(dim)^expon(dim)
!
!    where the exponents are nonnegative integers.  Note that
!    if the combination 0^0 is encountered, it should be treated
!    as 1.
!
!    Integral ( over unit triangle ) x^m y^n dx dy = m! * n! / ( m + n + 2 )!
!
!    The integration region is:
!
!      0 <= X
!      0 <= Y
!      X + Y <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer EXPON(2), the exponents.
!
!    Output, double precision VALUE, the integral of the monomial.
!
  implicit none

  integer expon(2)
  integer i
  integer k
  double precision value
!
!  The first computation ends with VALUE = 1.0;
!
  value = 1.0D+00

! k = 0
!
!  The first loop simply computes 1 so we short circuit it!
!
! do i = 1, expon(1)
!   k = k + 1
!   value = value * real ( i) / real ( k)
! end do

  k = expon(1)

  do i = 1, expon(2)
    k = k + 1
    value = value * real ( i) / real ( k)
  end do

  k = k + 1
  value = value / real ( k)

  k = k + 1
  value = value / real ( k)
end

subroutine trig_unit_monomial_test ( degree_max )

!*****************************************************************************80
!
!! TRIG_UNIT_MONOMIAL_TEST tests TRIG_UNIT_MONOMIAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DEGREE_MAX, the maximum total degree of the
!    monomials to check.
!
  implicit none

  integer alpha
  integer beta
  integer degree_max
  integer expon(2)
  double precision trig_unit_volume
  double precision value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIG_UNIT_MONOMIAL_TEST'
  write ( *, '(a)' ) '  For the unit triangle,'
  write ( *, '(a)' ) '  TRIG_UNIT_MONOMIAL returns the exact value of the'
  write ( *, '(a)' ) '  integral of X^ALPHA Y^BETA'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Volume = ', trig_unit_volume ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     ALPHA      BETA      INTEGRAL'
  write ( *, '(a)' ) ' '

  do alpha = 0, degree_max
    expon(1) = alpha
    do beta = 0, degree_max - alpha
      expon(2) = beta
      call trig_unit_monomial ( expon, value )
      write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) expon(1:2), value
    end do
  end do
end

subroutine trig_unit_o01 ( w, xy )

!*****************************************************************************80
!
!! TRIG_UNIT_O01 returns a 1 point quadrature rule for the unit triangle.
!
!  Discussion:
!
!    This rule is precise for monomials through degree 1.
!
!    The integration region is:
!
!      0 <= X
!      0 <= Y
!      X + Y <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(1), the weights.
!
!    Output, double precision XY(2,1), the abscissas.
!
  implicit none

  integer , parameter :: order = 1

  double precision w(order)
  double precision :: w_save(1) = (/ &
    1.0D+00 /)
  double precision xy(2,order)
  double precision :: xy_save(2,1) = reshape ( (/ &
    0.33333333333333333333D+00, &
    0.33333333333333333333D+00 /), &
  (/ 2, 1 /) )

  w(1:order) = w_save(1:order)
  xy(1:2,1:order) = xy_save(1:2,1:order)
end

subroutine trig_unit_o03 ( w, xy )

!*****************************************************************************80
!
!! TRIG_UNIT_O03 returns a 3 point quadrature rule for the unit triangle.
!
!  Discussion:
!
!    This rule is precise for monomials through degree 2.
!
!    The integration region is:
!
!      0 <= X
!      0 <= Y
!      X + Y <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(3), the weights.
!
!    Output, double precision XY(2,3), the abscissas.
!
  implicit none

  integer , parameter :: order = 3

  double precision w(order)
  double precision :: w_save(3) = (/ &
    0.33333333333333333333D+00, &
    0.33333333333333333333D+00, &
    0.33333333333333333333D+00 /)
  double precision xy(2,order)
  double precision :: xy_save(2,3) = reshape ( (/ &
    0.66666666666666666667D+00, &
    0.16666666666666666667D+00, &
    0.16666666666666666667D+00, &
    0.66666666666666666667D+00, &
    0.16666666666666666667D+00, &
    0.16666666666666666667D+00 /), &
  (/ 2, 3 /) )

  w(1:order) = w_save(1:order)
  xy(1:2,1:order) = xy_save(1:2,1:order)
end

subroutine trig_unit_o03b ( w, xy )

!*****************************************************************************80
!
!! TRIG_UNIT_O03B returns a 3 point quadrature rule for the unit triangle.
!
!  Discussion:
!
!    This rule is precise for monomials through degree 2.
!
!    The integration region is:
!
!      0 <= X
!      0 <= Y
!      X + Y <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(3), the weights.
!
!    Output, double precision XY(2,3), the abscissas.
!
  implicit none

  integer , parameter :: order = 3

  double precision w(order)
  double precision :: w_save(3) = (/ &
    0.33333333333333333333D+00, &
    0.33333333333333333333D+00, &
    0.33333333333333333333D+00 /)
  double precision xy(2,order)
  double precision :: xy_save(2,3) = reshape ( (/ &
    0.0D+00, &
    0.5D+00, &
    0.5D+00, &
    0.0D+00, &
    0.5D+00, &
    0.5D+00 /), &
  (/ 2, 3 /) )

  w(1:order) = w_save(1:order)
  xy(1:2,1:order) = xy_save(1:2,1:order)
end

subroutine trig_unit_o06 ( w, xy )

!*****************************************************************************80
!
!! TRIG_UNIT_O06 returns a 6 point quadrature rule for the unit triangle.
!
!  Discussion:
!
!    This rule is precise for monomials through degree 4.
!
!    The integration region is:
!
!      0 <= X
!      0 <= Y
!      X + Y <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(6), the weights.
!
!    Output, double precision XY(2,6), the abscissas.
!
  implicit none

  integer , parameter :: order = 6

  double precision w(order)
  double precision :: w_save(6) = (/ &
    0.22338158967801146570D+00, &
    0.22338158967801146570D+00, &
    0.22338158967801146570D+00, &
    0.10995174365532186764D+00, &
    0.10995174365532186764D+00, &
    0.10995174365532186764D+00 /)
  double precision xy(2,order)
  double precision :: xy_save(2,6) = reshape ( (/ &
    0.10810301816807022736D+00, &
    0.44594849091596488632D+00, &
    0.44594849091596488632D+00, &
    0.10810301816807022736D+00, &
    0.44594849091596488632D+00, &
    0.44594849091596488632D+00, &
    0.81684757298045851308D+00, &
    0.091576213509770743460D+00, &
    0.091576213509770743460D+00, &
    0.81684757298045851308D+00, &
    0.091576213509770743460D+00, &
    0.091576213509770743460D+00 /), &
  (/ 2, 6 /) )

  w(1:order) = w_save(1:order)
  xy(1:2,1:order) = xy_save(1:2,1:order)
end

subroutine trig_unit_o06b ( w, xy )

!*****************************************************************************80
!
!! TRIG_UNIT_O06B returns a 6 point quadrature rule for the unit triangle.
!
!  Discussion:
!
!    This rule is precise for monomials through degree 3.
!
!    The integration region is:
!
!      0 <= X
!      0 <= Y
!      X + Y <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(6), the weights.
!
!    Output, double precision XY(2,6), the abscissas.
!
  implicit none

  integer , parameter :: order = 6

  double precision w(order)
  double precision :: w_save(6) = (/ &
    0.30000000000000000000D+00, &
    0.30000000000000000000D+00, &
    0.30000000000000000000D+00, &
    0.033333333333333333333D+00, &
    0.033333333333333333333D+00, &
    0.033333333333333333333D+00 /)
  double precision xy(2,order)
  double precision :: xy_save(2,6) = reshape ( (/ &
    0.66666666666666666667D+00, &
    0.16666666666666666667D+00, &
    0.16666666666666666667D+00, &
    0.66666666666666666667D+00, &
    0.16666666666666666667D+00, &
    0.16666666666666666667D+00, &
    0.0D+00, &
    0.5D+00, &
    0.5D+00, &
    0.0D+00, &
    0.5D+00, &
    0.5D+00 /), &
  (/ 2, 6 /) )

  w(1:order) = w_save(1:order)
  xy(1:2,1:order) = xy_save(1:2,1:order)
end

subroutine trig_unit_o07 ( w, xy )

!*****************************************************************************80
!
!! TRIG_UNIT_O07 returns a 7 point quadrature rule for the unit triangle.
!
!  Discussion:
!
!    This rule is precise for monomials through degree 5.
!
!    The integration region is:
!
!      0 <= X
!      0 <= Y
!      X + Y <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(7), the weights.
!
!    Output, double precision XY(2,7), the abscissas.
!
  implicit none

  integer , parameter :: order = 7

  double precision w(order)
  double precision :: w_save(7) = (/ &
    0.12593918054482715260D+00, &
    0.12593918054482715260D+00, &
    0.12593918054482715260D+00, &
    0.13239415278850618074D+00, &
    0.13239415278850618074D+00, &
    0.13239415278850618074D+00, &
    0.22500000000000000000D+00 /)
  double precision xy(2,order)
  double precision :: xy_save(2,7) = reshape ( (/ &
    0.79742698535308732240D+00, &
    0.10128650732345633880D+00, &
    0.10128650732345633880D+00, &
    0.79742698535308732240D+00, &
    0.10128650732345633880D+00, &
    0.10128650732345633880D+00, &
    0.059715871789769820459D+00, &
    0.47014206410511508977D+00, &
    0.47014206410511508977D+00, &
    0.059715871789769820459D+00, &
    0.47014206410511508977D+00, &
    0.47014206410511508977D+00, &
    0.33333333333333333333D+00, &
    0.33333333333333333333D+00 /), &
  (/ 2, 7 /) )

  w(1:order) = w_save(1:order)
  xy(1:2,1:order) = xy_save(1:2,1:order)
end

subroutine trig_unit_o12 ( w, xy )

!*****************************************************************************80
!
!! TRIG_UNIT_O12 returns a 12 point quadrature rule for the unit triangle.
!
!  Discussion:
!
!    This rule is precise for monomials through degree 6.
!
!    The integration region is:
!
!      0 <= X
!      0 <= Y
!      X + Y <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Output, double precision W(12), the weights.
!
!    Output, double precision XY(2,12), the abscissas.
!
  implicit none

  integer , parameter :: order = 12

  double precision w(order)
  double precision :: w_save(12) = (/ &
     0.050844906370206816921D+00, &
     0.050844906370206816921D+00, &
     0.050844906370206816921D+00, &
     0.11678627572637936603D+00, &
     0.11678627572637936603D+00, &
     0.11678627572637936603D+00, &
     0.082851075618373575194D+00, &
     0.082851075618373575194D+00, &
     0.082851075618373575194D+00, &
     0.082851075618373575194D+00, &
     0.082851075618373575194D+00, &
     0.082851075618373575194D+00 /)
  double precision xy(2,order)
  double precision :: xy_save(2,12) = reshape ( (/ &
    0.87382197101699554332D+00, &
    0.063089014491502228340D+00, &
    0.063089014491502228340D+00, &
    0.87382197101699554332D+00, &
    0.063089014491502228340D+00, &
    0.063089014491502228340D+00, &
    0.50142650965817915742D+00, &
    0.24928674517091042129D+00, &
    0.24928674517091042129D+00, &
    0.50142650965817915742D+00, &
    0.24928674517091042129D+00, &
    0.24928674517091042129D+00, &
    0.053145049844816947353D+00, &
    0.31035245103378440542D+00, &
    0.31035245103378440542D+00, &
    0.053145049844816947353D+00, &
    0.053145049844816947353D+00, &
    0.63650249912139864723D+00, &
    0.31035245103378440542D+00, &
    0.63650249912139864723D+00, &
    0.63650249912139864723D+00, &
    0.053145049844816947353D+00, &
    0.63650249912139864723D+00, &
    0.31035245103378440542D+00 /), &
  (/ 2, 12 /) )

  w(1:order) = w_save(1:order)
  xy(1:2,1:order) = xy_save(1:2,1:order)
end

subroutine trig_unit_quad_test ( degree_max )

!*****************************************************************************80
!
!! TRIG_UNIT_QUAD_TEST tests the rules for the unit triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DEGREE_MAX, the maximum total degree of the
!    monomials to check.
!
  implicit none

  integer , parameter :: dim_num = 2

  integer degree_max
  integer expon(dim_num)
  integer h
  logical more
  integer order
  double precision quad
  integer t
  double precision trig_unit_volume
  double precision , allocatable :: v(:)
  double precision , allocatable :: w(:)
  double precision , allocatable :: xy(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIG_UNIT_QUAD_TEST'
  write ( *, '(a)' ) '  For the unit triangle,'
  write ( *, '(a)' ) '  we approximate monomial integrals with:'
  write ( *, '(a)' ) '  TRIG_UNIT_O01,'
  write ( *, '(a)' ) '  TRIG_UNIT_O03,'
  write ( *, '(a)' ) '  TRIG_UNIT_O03b,'
  write ( *, '(a)' ) '  TRIG_UNIT_O06,'
  write ( *, '(a)' ) '  TRIG_UNIT_O06b,'
  write ( *, '(a)' ) '  TRIG_UNIT_O07,'
  write ( *, '(a)' ) '  TRIG_UNIT_O012,'

  more = .false.

  do

    call subcomp_next ( degree_max, dim_num, expon, more, h, t )

    write ( *, '(a)' ) ' '
    write ( *, '(a,2x,i2,2x,i2)' ) '  Monomial exponents: ', expon(1:dim_num)
    write ( *, '(a)' ) ' '

    order = 1
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xy(1:dim_num,1:order) )
    call trig_unit_o01 ( w, xy )
    call monomial_value ( dim_num, order, expon, xy, v )
    quad = trig_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xy )

    order = 3
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xy(1:dim_num,1:order) )
    call trig_unit_o03 ( w, xy )
    call monomial_value ( dim_num, order, expon, xy, v )
    quad = trig_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xy )

    order = 3
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xy(1:dim_num,1:order) )
    call trig_unit_o03b ( w, xy )
    call monomial_value ( dim_num, order, expon, xy, v )
    quad = trig_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xy )

    order = 6
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xy(1:dim_num,1:order) )
    call trig_unit_o06 ( w, xy )
    call monomial_value ( dim_num, order, expon, xy, v )
    quad = trig_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xy )

    order = 6
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xy(1:dim_num,1:order) )
    call trig_unit_o06b ( w, xy )
    call monomial_value ( dim_num, order, expon, xy, v )
    quad = trig_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xy )

    order = 7
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xy(1:dim_num,1:order) )
    call trig_unit_o07 ( w, xy )
    call monomial_value ( dim_num, order, expon, xy, v )
    quad = trig_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xy )

    order = 12
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xy(1:dim_num,1:order) )
    call trig_unit_o12 ( w, xy )
    call monomial_value ( dim_num, order, expon, xy, v )
    quad = trig_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xy )

    write ( *, '(a)' ) ' '
    call trig_unit_monomial ( expon, quad )
    write ( *, '(2x,a,2x,g14.6)' ) ' Exact', quad

    if ( .not. more ) then
      exit
    end if

  end do
end

function trig_unit_volume ( )

!*****************************************************************************80
!
!! TRIG_UNIT_VOLUME: volume of a unit triangle.
!
!  Discussion:
!
!    The integration region is:
!
!      0 <= X
!      0 <= Y
!      X + Y <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, double precision TRIG_UNIT_VOLUME, the volume.
!
  implicit none

  double precision trig_unit_volume

  trig_unit_volume = 0.5D+00
end

subroutine wedg_unit_monomial ( expon, value )

!*****************************************************************************80
!
!! WEDG_UNIT_MONOMIAL: monomial integral in a unit wedge.
!
!  Discussion:
!
!    This routine returns the integral of
!
!      product ( 1 <= I <= 3 ) X(I)^EXPON(I)
!
!    over the unit wedge.
!
!    The integration region is:
!
!      0 <= X
!      0 <= Y
!      X + Y <= 1
!      -1 <= Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer EXPON(3), the exponents.
!
!    Output, double precision VALUE, the integral of the monomial.
!
  implicit none

  integer expon(3)
  integer i
  integer k
  double precision value
!
!  The first computation ends with VALUE = 1.0;
!
  value = 1.0D+00

! k = 0
!
!  The first loop simply computes 1 so we short circuit it!
!
! do i = 1, expon(1)
!   k = k + 1
!   value = value * real ( i) / real ( k)
! end do

  k = expon(1)

  do i = 1, expon(2)
    k = k + 1
    value = value * real ( i) / real ( k)
  end do

  k = k + 1
  value = value / real ( k)

  k = k + 1
  value = value / real ( k)
!
!  Now account for integration in Z.
!
  if ( expon(3) == - 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WEDG_UNIT_MONOMIAL - Fatal error!'
    write ( *, '(a)' ) '  EXPON(3) = -1 is not a legal input.'
    stop 1
  else if ( mod ( expon(3), 2 ) == 1 ) then
    value = 0.0D+00
  else
    value = value * 2.0D+00 / real ( expon(3) + 1)
  end if
end

subroutine wedg_unit_monomial_test ( degree_max )

!*****************************************************************************80
!
!! WEDG_UNIT_MONOMIAL_TEST tests WEDG_UNIT_MONOMIAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DEGREE_MAX, the maximum total degree of the
!    monomials to check.
!
  implicit none

  integer alpha
  integer beta
  integer degree_max
  integer expon(3)
  integer gamma
  double precision value
  double precision wedg_unit_volume

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WEDG_UNIT_MONOMIAL_TEST'
  write ( *, '(a)' ) '  For the unit wedge,'
  write ( *, '(a)' ) '  WEDG_UNIT_MONOMIAL returns the exact value of the'
  write ( *, '(a)' ) '  integral of X^ALPHA Y^BETA Z^GAMMA'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Volume = ', wedg_unit_volume ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     ALPHA      BETA     GAMMA      INTEGRAL'
  write ( *, '(a)' ) ' '

  do alpha = 0, degree_max
    expon(1) = alpha
    do beta = 0, degree_max - alpha
      expon(2) = beta
      do gamma = 0, degree_max - alpha - beta
        expon(3) = gamma
        call wedg_unit_monomial ( expon, value )
        write ( *, '(2x,i8,2x,i8,2x,i8,2x,g14.6)' ) expon(1:3), value
      end do
    end do
  end do
end

subroutine wedg_unit_quad_test ( degree_max )

!*****************************************************************************80
!
!! WEDG_UNIT_QUAD_TEST tests the rules for the unit wedge.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DEGREE_MAX, the maximum total degree of the
!    monomials to check.
!
  implicit none

  integer , parameter :: dim_num = 3
  integer , parameter :: test_num = 7

  integer degree_max
  integer expon(dim_num)
  integer h
  integer line_order
  integer :: line_order_array(test_num) = (/ &
    1, 2, 2, 3, 2, 3, 4 /)
  logical more
  integer order
  double precision quad
  integer t
  integer test
  integer trig_order
  integer trig_order_index
  integer :: trig_order_array(test_num) = (/ &
    1, 3, -3, 6, -6, 7, 12 /)
  double precision wedg_unit_volume
  double precision , allocatable :: v(:)
  double precision , allocatable :: w(:)
  double precision , allocatable :: xyz(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WEDG_UNIT_QUAD_TEST'
  write ( *, '(a)' ) '  For the unit wedge,'
  write ( *, '(a)' ) '  we approximate monomial integrals with WEDG_UNIT_RULE.'

  more = .false.

  do

    call subcomp_next ( degree_max, dim_num, expon, more, h, t )

    if ( mod ( expon(3), 2 ) == 1 ) then
      if ( .not. more ) then
        exit
      else
        cycle
      end if
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,2x,i2,2x,i2,2x,i2)' ) &
      '  Monomial exponents: ', expon(1:dim_num)
    write ( *, '(a)' ) ' '

    do test = 1, test_num

      line_order = line_order_array(test)
      trig_order = trig_order_array(test)

      order = line_order * abs ( trig_order )

      allocate ( v(1:order) )
      allocate ( w(1:order) )
      allocate ( xyz(1:dim_num,1:order) )
      call wedg_unit_rule ( line_order, trig_order, w, xyz )
      call monomial_value ( dim_num, order, expon, xyz, v )
      quad = wedg_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
      write ( *, '(2x,i6,2x,i6,2x,i6,2x,g14.6)' ) &
        trig_order, line_order, order, quad
      deallocate ( v )
      deallocate ( w )
      deallocate ( xyz )

    end do

    write ( *, '(a)' ) ' '
    call wedg_unit_monomial ( expon, quad )
    write ( *, '(2x,a,2x,6x,2x,6x,2x,g14.6)' ) ' Exact', quad

    if ( .not. more ) then
      exit
    end if

  end do
end

subroutine wedg_unit_rule ( line_order, trig_order, w, xyz )

!*****************************************************************************80
!
!! WEDG_UNIT_RULE returns a quadrature rule for the unit wedge.
!
!  Discussion:
!
!    It is usually sensible to take LINE_ORDER and TRIG_ORDER so that
!    the line and triangle rules are roughly the same precision.  For that
!    criterion, we recommend the following combinations:
!
!      TRIG_ORDER  LINE_ORDER  Precision
!      ----------  ----------  ---------
!          1           1       1 x 1
!          3           2       2 x 3
!         -3           2       2 x 3
!          6           3       4 x 5
!         -6           2       3 x 3
!          7           3       5 x 5
!         12           4       6 x 7
!
!    The integration region is:
!
!      0 <= X
!      0 <= Y
!      X + Y <= 1
!      -1 <= Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carlos Felippa,
!    A compendium of FEM integration formulas for symbolic work,
!    Engineering Computation,
!    Volume 21, Number 8, 2004, pages 867-890.
!
!  Parameters:
!
!    Input, integer LINE_ORDER, the index of the line rule.
!    The index of the rule is equal to the order of the rule.
!    1 <= LINE_ORDER <= 5.
!
!    Input, integer TRIG_ORDER, the indes of the triangle rule.
!    The index of the rule is 1, 3, -3, 6, -6, 7 or 12.
!
!    Output, double precision W(LINE_ORDER*abs(TRIG_ORDER)), the weights.
!
!    Output, double precision XYZ(3,LINE_ORDER*abs(TRIG_ORDER)), the abscissas.
!
  implicit none

  integer line_order
  integer trig_order

  integer i
  integer j
  integer k
  double precision line_w(line_order)
  double precision line_x(line_order)
  double precision trig_w(abs(trig_order))
  double precision trig_xy(2,abs(trig_order))
  double precision w(line_order*abs(trig_order))
  double precision xyz(3,line_order*abs(trig_order))

  if ( line_order == 1 ) then
    call line_unit_o01 ( line_w, line_x )
  else if ( line_order == 2 ) then
    call line_unit_o02 ( line_w, line_x )
  else if ( line_order == 3 ) then
    call line_unit_o03 ( line_w, line_x )
  else if ( line_order == 4 ) then
    call line_unit_o04 ( line_w, line_x )
  else if ( line_order == 5 ) then
    call line_unit_o05 ( line_w, line_x )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WEDG_UNIT_RULE - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of LINE_ORDER.'
    stop 1
  end if

  if ( trig_order == 1 ) then
    call trig_unit_o01 ( trig_w, trig_xy )
  else if ( trig_order == 3 ) then
    call trig_unit_o03 ( trig_w, trig_xy )
  else if ( trig_order == - 3 ) then
    call trig_unit_o03b ( trig_w, trig_xy )
  else if ( trig_order == 6 ) then
    call trig_unit_o06 ( trig_w, trig_xy )
  else if ( trig_order == - 6 ) then
    call trig_unit_o06b ( trig_w, trig_xy )
  else if ( trig_order == 7 ) then
    call trig_unit_o07 ( trig_w, trig_xy )
  else if ( trig_order == 12 ) then
    call trig_unit_o12 ( trig_w, trig_xy )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WEDG_UNIT_RULE - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of TRIG_ORDER.'
    stop 1
  end if

  k = 0
  do i = 1, line_order
    do j = 1, abs ( trig_order )
      k = k + 1
      w(k) = line_w(i) * trig_w(j)
      xyz(1:2,k) = trig_xy(1:2,j)
      xyz(3,k) = line_x(i)
    end do
  end do
end

function wedg_unit_volume ( )

!*****************************************************************************80
!
!! WEDG_UNIT_VOLUME: volume of a unit wedge.
!
!  Discussion:
!
!    The integration region is:
!
!      0 <= X
!      0 <= Y
!      X + Y <= 1
!      -1 <= Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, double precision WEDG_UNIT_VOLUME, the volume.
!
  implicit none

  double precision wedg_unit_volume

  wedg_unit_volume = 1.0D+00
end

subroutine wedg_unit_write_test ( )

!*****************************************************************************80
!
!! WEDG_UNIT_WRITE_TEST writes out some rules for the unit wedge.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer , parameter :: dim_num = 3
  integer , parameter :: rule_num = 7

  integer line_order
  integer :: line_order_array(rule_num) = (/ &
    1, 2, 2, 3, 2, 3, 4 /)
  integer order
  integer rule
  integer trig_order
  integer :: trig_order_array(rule_num) = (/ &
    1, 3, -3, 6, -6, 7, 12 /)
  double precision , allocatable :: w(:)
  character ( len = 255 ) w_filename
  double precision , allocatable :: x(:,:)
  character ( len = 255 ) x_filename

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WEDG_UNIT_WRITE_TEST'
  write ( *, '(a)' ) '  For the unit wedge,'
  write ( *, '(a)' ) '  write some rules to a file'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Rule  Trig    Line   Total  W_File X_File'
  write ( *, '(a)' ) '         Order   Order  Order'
  write ( *, '(a)' ) ' '

  do rule = 1, rule_num

    if ( rule == 1 ) then
      w_filename = 'wedge_felippa_1x1_w.txt'
      x_filename = 'wedge_felippa_1x1_x.txt'
    else if ( rule == 2 ) then
      w_filename = 'wedge_felippa_3x2_w.txt'
      x_filename = 'wedge_felippa_3x2_x.txt'
    else if ( rule == 3 ) then
      w_filename = 'wedge_felippa_3bx2_w.txt'
      x_filename = 'wedge_felippa_3bx2_x.txt'
    else if ( rule == 4 ) then
      w_filename = 'wedge_felippa_6x3_w.txt'
      x_filename = 'wedge_felippa_6x3_x.txt'
    else if ( rule == 5 ) then
      w_filename = 'wedge_felippa_6bx2_w.txt'
      x_filename = 'wedge_felippa_6bx2_x.txt'
    else if ( rule == 6 ) then
      w_filename = 'wedge_felippa_7x3_w.txt'
      x_filename = 'wedge_felippa_7x3_x.txt'
    else if ( rule == 7 ) then
      w_filename = 'wedge_felippa_12x4_w.txt'
      x_filename = 'wedge_felippa_12x4_x.txt'
    end if

    line_order = line_order_array(rule)
    trig_order = trig_order_array(rule)

    order = line_order * abs ( trig_order )

    allocate ( w(1:order) )
    allocate ( x(1:dim_num,1:order) )
    call wedg_unit_rule ( line_order, trig_order, w, x )
    call r8mat_write ( w_filename, 1, order, w )
    call r8mat_write ( x_filename, dim_num, order, x )
    write ( *, '(2x,i6,2x,i6,2x,i6,2x,i6,2x,a25,2x,a25)' ) &
      rule, trig_order, line_order, order, w_filename, x_filename

    deallocate ( w )
    deallocate ( x )

  end do
end
