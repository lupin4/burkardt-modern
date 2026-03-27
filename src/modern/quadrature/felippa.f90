!> felippa � Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module felippa_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: comp_next, hexa_unit_monomial, hexa_unit_monomial_test, hexa_unit_quad_test, hexa_unit_rule, hexa_unit_volume
  public :: line_unit_monomial, line_unit_monomial_test, line_unit_o01, line_unit_o02, line_unit_o03, line_unit_o04
  public :: line_unit_o05, line_unit_quad_test, line_unit_volume, monomial_value, pyra_unit_monomial, pyra_unit_monomial_test
  public :: pyra_unit_o01, pyra_unit_o05, pyra_unit_o06, pyra_unit_o08, pyra_unit_o08b, pyra_unit_o09
  public :: pyra_unit_o13, pyra_unit_o18, pyra_unit_o27, pyra_unit_o48, pyra_unit_quad_test, pyra_unit_volume
  public :: quad_unit_monomial, quad_unit_monomial_test, quad_unit_quad_test, quad_unit_rule, quad_unit_volume, r8_choose
  public :: r8_mop, r8mat_det_4d, r8vec_direct_product, r8vec_direct_product2, subcomp_next, tetr_unit_monomial
  public :: tetr_unit_monomial_test, tetr_unit_o01, tetr_unit_o04, tetr_unit_o08, tetr_unit_o08b, tetr_unit_o14
  public :: tetr_unit_o14b, tetr_unit_o15, tetr_unit_o15b, tetr_unit_o24, tetr_unit_quad_test, tetr_unit_volume
  public :: trig_unit_monomial, trig_unit_monomial_test, trig_unit_o01, trig_unit_o03, trig_unit_o03b, trig_unit_o06
  public :: trig_unit_o06b, trig_unit_o07, trig_unit_o12, trig_unit_quad_test, trig_unit_volume, wedg_unit_monomial
  public :: wedg_unit_monomial_test, wedg_unit_quad_test, wedg_unit_rule, wedg_unit_volume, wedg_unit_write_test

contains

  pure subroutine comp_next ( n, k, a, more, h, t ) &
        bind(C, name="comp_next")

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
  !    Input, integer(ip) N, the integer whose compositions are desired.
  !
  !    Input, integer(ip) K, the number of parts in the composition.
  !
  !    Input/output, integer(ip) A(K), the parts of the composition.
  !
  !    Input/output, logical MORE, set by the user to start the
  !    computation, and by the routine to terminate it.
  !
  !    Input/output, integer(ip)  H, T, two internal parameters needed
  !    for the computation.  The user should allocate space for these in the
  !    calling program, include them in the calling sequence, but never alter
  !    them!
  !

    integer(ip), intent(in), value :: k

    integer(ip), intent(inout) :: a(k)
    integer(ip), intent(inout) :: h
    logical, intent(inout) :: more
    integer(ip), intent(in), value :: n
    integer(ip), intent(inout) :: t
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
  end subroutine comp_next

  subroutine hexa_unit_monomial ( expon, value ) &
        bind(C, name="hexa_unit_monomial")

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
  !    Input, integer(ip) EXPON(3), the exponents.
  !
  !    Output, real(dp) VALUE, the integral of the monomial.
  !

    integer(ip), intent(in) :: expon(3)
    integer(ip) :: i
    real(dp), intent(out) :: value

    value = 1.0_dp

    do i = 1, 3

      if ( mod ( expon(i), 2 ) == 1 ) then
        value = 0.0_dp
      else if ( expon(i) == -1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HEXA_UNIT_MONOMIAL - Fatal error!'
        write ( *, '(a)' ) '  Exponent of -1 encountered.'
        stop 1
      else
        value = value * 2.0_dp / real ( expon(i) + 1, dp)
      end if

    end do
  end subroutine hexa_unit_monomial

  subroutine hexa_unit_monomial_test ( degree_max ) &
        bind(C, name="hexa_unit_monomial_test")

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
  !    Input, integer(ip) DEGREE_MAX, the maximum total degree of the
  !    monomials to check.
  !

    integer(ip) :: alpha
    integer(ip) :: beta
    integer(ip), intent(in), value :: degree_max
    integer(ip) :: expon(3)
    integer(ip) :: gamma
    real(dp) :: hexa_unit_volume
    real(dp) :: value

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
  end subroutine hexa_unit_monomial_test

  subroutine hexa_unit_quad_test ( degree_max ) &
        bind(C, name="hexa_unit_quad_test")

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
  !    Input, integer(ip) DEGREE_MAX, the maximum total degree of the
  !    monomials to check.
  !

    integer(ip), parameter :: dim_num = 3

    integer(ip), intent(in), value :: degree_max
    integer(ip) :: expon(dim_num)
    integer(ip) :: h
    real(dp) :: hexa_unit_volume
    integer(ip) :: k
    logical :: more
    integer(ip) :: order
    integer(ip) :: order_1d(dim_num)
    real(dp) :: quad
    integer(ip) :: t
    real(dp), allocatable :: v(:)
    real(dp), allocatable :: w(:)
    real(dp), allocatable :: xyz(:,:)

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
  end subroutine hexa_unit_quad_test

  subroutine hexa_unit_rule ( order_1d, w, xyz ) &
        bind(C, name="hexa_unit_rule")

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
  !    Input, integer(ip) ORDER_1D(3), the order of the rule in each
  !    dimension.  1 <= ORDER_1D(I) <= 5.
  !
  !    Output, real(dp) W(ORDER_1D(1)*ORDER_1D(2)*ORDER_1D(3)), 
  !    the weights.
  !
  !    Output, real(dp) XYZ(3,ORDER_1D(1)*ORDER_1D(2)*ORDER_1D(3)), 
  !    the abscissas.
  !

    integer(ip), parameter :: dim_num = 3

    integer(ip) :: dim
    integer(ip) :: order
    integer(ip), intent(in) :: order_1d(dim_num)
    real(dp), intent(out) :: w(order_1d(1)*order_1d(2)*order_1d(3))
    real(dp), allocatable :: w_1d(:)
    real(dp), allocatable :: x_1d(:)
    real(dp), intent(out) :: xyz(3,order_1d(1)*order_1d(2)*order_1d(3))

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
  end subroutine hexa_unit_rule

  pure function hexa_unit_volume ( ) &
        bind(C, name="hexa_unit_volume")

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
  !    Output, real(dp) HEXA_UNIT_VOLUME, the volume.
  !

    real(dp) :: hexa_unit_volume

    hexa_unit_volume = 8.0_dp
  end function hexa_unit_volume

  subroutine line_unit_monomial ( alpha, value ) &
        bind(C, name="line_unit_monomial")

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
  !    Input, integer(ip) ALPHA, the exponent of X.
  !    ALPHA must not be -1.
  !
  !    Output, real(dp) value, the integral of the monomial.
  !

    integer(ip), intent(in), value :: alpha
    real(dp), intent(out) :: value

    if ( alpha == - 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LINE_UNIT_MONOMIAL - Fatal error!'
      write ( *, '(a)' ) '  ALPHA = -1 is not a legal input.'
      stop 1
    else if ( mod ( alpha, 2 ) == 1 ) then
      value = 0.0_dp
    else
      value = 2.0_dp / real ( alpha + 1, dp)
    end if
  end subroutine line_unit_monomial

  subroutine line_unit_monomial_test ( degree_max ) &
        bind(C, name="line_unit_monomial_test")

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
  !    Input, integer(ip) DEGREE_MAX, the maximum total degree of the
  !    monomials to check.
  !

    integer(ip) :: alpha
    integer(ip), intent(in), value :: degree_max
    real(dp) :: line_unit_volume
    real(dp) :: value

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
  end subroutine line_unit_monomial_test

  pure subroutine line_unit_o01 ( w, x ) &
        bind(C, name="line_unit_o01")

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
  !    Output, real(dp) W(1), the weights.
  !
  !    Output, real(dp) X(1), the abscissas.
  !

    integer(ip), parameter :: order = 1

    real(dp) :: line_unit_volume
    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(1) = (/ &
      2.0_dp /)
    real(dp), intent(out) :: x(order)
    real(dp) :: x_save(1) = (/ &
      0.0_dp /)

    w(1:order) = w_save(1:order) / line_unit_volume ( )
    x(1:order) = x_save(1:order)
  end subroutine line_unit_o01

  pure subroutine line_unit_o02 ( w, x ) &
        bind(C, name="line_unit_o02")

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
  !    Output, real(dp) W(2), the weights.
  !
  !    Output, real(dp) X(2), the abscissas.
  !

    integer(ip), parameter :: order = 2

    real(dp) :: line_unit_volume
    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(2) = (/ &
      1.0000000000000000000_dp, &
      1.0000000000000000000_dp /)
    real(dp), intent(out) :: x(order)
    real(dp) :: x_save(2) = (/ &
      -0.57735026918962576451_dp, &
       0.57735026918962576451_dp /)

    w(1:order) = w_save(1:order) / line_unit_volume ( )
    x(1:order) = x_save(1:order)
  end subroutine line_unit_o02

  pure subroutine line_unit_o03 ( w, x ) &
        bind(C, name="line_unit_o03")

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
  !    Output, real(dp) W(3), the weights.
  !
  !    Output, real(dp) X(3), the abscissas.
  !

    integer(ip), parameter :: order = 3

    real(dp) :: line_unit_volume
    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(3) = (/ &
      0.55555555555555555556_dp, &
      0.88888888888888888889_dp, &
      0.55555555555555555556_dp /)
    real(dp), intent(out) :: x(order)
    real(dp) :: x_save(3) = (/ &
      -0.77459666924148337704_dp, &
       0.00000000000000000000_dp, &
       0.77459666924148337704_dp /)

    w(1:order) = w_save(1:order) / line_unit_volume ( )
    x(1:order) = x_save(1:order)
  end subroutine line_unit_o03

  pure subroutine line_unit_o04 ( w, x ) &
        bind(C, name="line_unit_o04")

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
  !    Output, real(dp) W(4), the weights.
  !
  !    Output, real(dp) X(4), the abscissas.
  !

    integer(ip), parameter :: order = 4

    real(dp) :: line_unit_volume
    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(4) = (/ &
      0.34785484513745385737_dp, &
      0.65214515486254614263_dp, &
      0.65214515486254614263_dp, &
      0.34785484513745385737_dp /)
    real(dp), intent(out) :: x(order)
    real(dp) :: x_save(4) = (/ &
      -0.86113631159405257522_dp, &
      -0.33998104358485626480_dp, &
       0.33998104358485626480_dp, &
       0.86113631159405257522_dp /)

    w(1:order) = w_save(1:order) / line_unit_volume ( )
    x(1:order) = x_save(1:order)
  end subroutine line_unit_o04

  pure subroutine line_unit_o05 ( w, x ) &
        bind(C, name="line_unit_o05")

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
  !    Output, real(dp) W(5), the weights.
  !
  !    Output, real(dp) X(5), the abscissas.
  !

    integer(ip), parameter :: order = 5

    real(dp) :: line_unit_volume
    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(5) = (/ &
      0.23692688505618908751_dp, &
      0.47862867049936646804_dp, &
      0.56888888888888888889_dp, &
      0.47862867049936646804_dp, &
      0.23692688505618908751_dp /)
    real(dp), intent(out) :: x(order)
    real(dp) :: x_save(5) = (/ &
      -0.90617984593866399280_dp, &
      -0.53846931010568309104_dp, &
       0.00000000000000000000_dp, &
       0.53846931010568309104_dp, &
       0.90617984593866399280_dp /)

    w(1:order) = w_save(1:order) / line_unit_volume ( )
    x(1:order) = x_save(1:order)
  end subroutine line_unit_o05

  subroutine line_unit_quad_test ( degree_max ) &
        bind(C, name="line_unit_quad_test")

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
  !    Input, integer(ip) DEGREE_MAX, the maximum total degree of the
  !    monomials to check.
  !

    integer(ip), parameter :: dim_num = 1

    integer(ip), intent(in), value :: degree_max
    integer(ip) :: expon(dim_num)
    integer(ip) :: h
    real(dp) :: line_unit_volume
    logical :: more
    integer(ip) :: order
    real(dp) :: quad
    integer(ip) :: t
    real(dp), allocatable :: v(:)
    real(dp), allocatable :: w(:)
    real(dp), allocatable :: x(:,:)

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
      call line_unit_monomial ( expon, quad )
      write ( *, '(2x,a,2x,g14.6)' ) ' Exact', quad

      if ( .not. more ) then
        exit
      end if

    end do
  end subroutine line_unit_quad_test

  pure function line_unit_volume ( ) &
        bind(C, name="line_unit_volume")

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
  !    Output, real(dp) LINE_UNIT_VOLUME, the volume.
  !

    real(dp) :: line_unit_volume

    line_unit_volume = 2.0_dp
  end function line_unit_volume

  pure subroutine monomial_value ( dim_num, point_num, expon, x, v ) &
        bind(C, name="monomial_value")

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
  !    Input, integer(ip) DIM_NUM, the spatial dimension.
  !
  !    Input, integer(ip) POINT_NUM, the number of points.
  !
  !    Input, integer(ip) EXPON(DIM_NUM), the exponents.
  !
  !    Input, real(dp) X(DIM_NUM,POINT_NUM), the evaluation points.
  !
  !    Output, real(dp) V(POINT_NUM), the monomial values.
  !

    integer(ip), intent(in), value :: dim_num
    integer(ip), intent(in), value :: point_num

    integer(ip) :: dim
    integer(ip), intent(in) :: expon(dim_num)
    real(dp), intent(out) :: v(point_num)
    real(dp), intent(in) :: x(dim_num,point_num)

    v(1:point_num) = 1.0_dp

    do dim = 1, dim_num
      if ( expon(dim) /= 0.0_dp ) then
        v(1:point_num) = v(1:point_num) * x(dim,1:point_num)**expon(dim)
      end if
    end do
  end subroutine monomial_value

  pure subroutine pyra_unit_monomial ( expon, value ) &
        bind(C, name="pyra_unit_monomial")

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
  !    Input, integer(ip) EXPON(3), the exponents.
  !
  !    Output, real(dp) VALUE, the integral of the monomial.
  !

    integer(ip), intent(in) :: expon(3)
    integer(ip) :: i
    integer(ip) :: i_hi
    real(dp) :: r8_choose
    real(dp) :: r8_mop
    real(dp), intent(out) :: value

    value = 0.0_dp

    if ( mod ( expon(1), 2 ) == 0 .and. mod ( expon(2), 2 ) == 0 ) then

      i_hi = 2 + expon(1) + expon(2)

      do i = 0, i_hi
        value = value + r8_mop ( i ) * r8_choose ( i_hi, i ) &
        / real ( i + expon(3) + 1, dp)
      end do

      value = value &
            * 2.0_dp / real ( expon(1) + 1, dp) &
            * 2.0_dp / real ( expon(2) + 1, dp)

    end if
  end subroutine pyra_unit_monomial

  subroutine pyra_unit_monomial_test ( degree_max ) &
        bind(C, name="pyra_unit_monomial_test")

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
  !    Input, integer(ip) DEGREE_MAX, the maximum total degree of the
  !    monomials to check.
  !

    integer(ip) :: alpha
    integer(ip) :: beta
    integer(ip), intent(in), value :: degree_max
    integer(ip) :: expon(3)
    integer(ip) :: gamma
    real(dp) :: pyra_unit_volume
    real(dp) :: value

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
  end subroutine pyra_unit_monomial_test

  pure subroutine pyra_unit_o01 ( w, xyz ) &
        bind(C, name="pyra_unit_o01")

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
  !    Output, real(dp) W(1), the weights.
  !
  !    Output, real(dp) XYZ(3,1), the abscissas.
  !

    integer(ip), parameter :: order = 1

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(1) = (/ &
      1.0_dp /)
    real(dp), intent(out) :: xyz(3,order)
    real(dp) :: xyz_save(3,1) = reshape ( (/ &
      0.0_dp, 0.0_dp, 0.25_dp /), &
    (/ 3, 1 /) )

    w(1:order) = w_save(1:order)
    xyz(1:3,1:order) = xyz_save(1:3,1:order)
  end subroutine pyra_unit_o01

  pure subroutine pyra_unit_o05 ( w, xyz ) &
        bind(C, name="pyra_unit_o05")

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
  !    Output, real(dp) W(5), the weights.
  !
  !    Output, real(dp) XYZ(3,5), the abscissas.
  !

    integer(ip), parameter :: order = 5

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(5) = (/ &
     0.21093750000000000000_dp, &
     0.21093750000000000000_dp, &
     0.21093750000000000000_dp, &
     0.21093750000000000000_dp, &
     0.15625000000000000000_dp /)
    real(dp), intent(out) :: xyz(3,order)
    real(dp) :: xyz_save(3,5) = reshape ( (/ &
    -0.48686449556014765641_dp, &
    -0.48686449556014765641_dp, &
     0.16666666666666666667_dp, &
     0.48686449556014765641_dp, &
    -0.48686449556014765641_dp, &
     0.16666666666666666667_dp, &
     0.48686449556014765641_dp, &
     0.48686449556014765641_dp, &
     0.16666666666666666667_dp, &
    -0.48686449556014765641_dp, &
     0.48686449556014765641_dp, &
     0.16666666666666666667_dp, &
     0.00000000000000000000_dp, &
     0.00000000000000000000_dp, &
     0.70000000000000000000_dp /), &
    (/ 3, 5 /) )

    w(1:order) = w_save(1:order)
    xyz(1:3,1:order) = xyz_save(1:3,1:order)
  end subroutine pyra_unit_o05

  pure subroutine pyra_unit_o06 ( w, xyz ) &
        bind(C, name="pyra_unit_o06")

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
  !    Output, real(dp) W(6), the weights.
  !
  !    Output, real(dp) XYZ(3,6), the abscissas.
  !

    integer(ip), parameter :: order = 6

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(6) = (/ &
     0.21000000000000000000_dp, &
     0.21000000000000000000_dp, &
     0.21000000000000000000_dp, &
     0.21000000000000000000_dp, &
     0.06000000000000000000_dp, &
     0.10000000000000000000_dp /)
    real(dp), intent(out) :: xyz(3,order)
    real(dp) :: xyz_save(3,6) = reshape ( (/ &
    -0.48795003647426658968_dp, &
    -0.48795003647426658968_dp, &
     0.16666666666666666667_dp, &
     0.48795003647426658968_dp, &
    -0.48795003647426658968_dp, &
     0.16666666666666666667_dp, &
     0.48795003647426658968_dp, &
     0.48795003647426658968_dp, &
     0.16666666666666666667_dp, &
    -0.48795003647426658968_dp, &
     0.48795003647426658968_dp, &
     0.16666666666666666667_dp, &
     0.00000000000000000000_dp, &
     0.00000000000000000000_dp, &
     0.58333333333333333333_dp, &
     0.00000000000000000000_dp, &
     0.00000000000000000000_dp, &
     0.75000000000000000000_dp /), &
    (/ 3, 6 /) )

    w(1:order) = w_save(1:order)
    xyz(1:3,1:order) = xyz_save(1:3,1:order)
  end subroutine pyra_unit_o06

  pure subroutine pyra_unit_o08 ( w, xyz ) &
        bind(C, name="pyra_unit_o08")

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
  !    Output, real(dp) W(8), the weights.
  !
  !    Output, real(dp) XYZ(3,8), the abscissas.
  !

    integer(ip), parameter :: order = 8

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(8) = (/ &
     0.075589411559869072938_dp, &
     0.075589411559869072938_dp, &
     0.075589411559869072938_dp, &
     0.075589411559869072938_dp, &
     0.17441058844013092706_dp, &
     0.17441058844013092706_dp, &
     0.17441058844013092706_dp, &
     0.17441058844013092706_dp /)
    real(dp), intent(out) :: xyz(3,order)
    real(dp) :: xyz_save(3,8) = reshape ( (/ &
    -0.26318405556971359557_dp, &
    -0.26318405556971359557_dp, &
     0.54415184401122528880_dp, &
     0.26318405556971359557_dp, &
    -0.26318405556971359557_dp, &
     0.54415184401122528880_dp, &
     0.26318405556971359557_dp, &
     0.26318405556971359557_dp, &
     0.54415184401122528880_dp, &
    -0.26318405556971359557_dp, &
     0.26318405556971359557_dp, &
     0.54415184401122528880_dp, &
    -0.50661630334978742377_dp, &
    -0.50661630334978742377_dp, &
     0.12251482265544137787_dp, &
     0.50661630334978742377_dp, &
    -0.50661630334978742377_dp, &
     0.12251482265544137787_dp, &
     0.50661630334978742377_dp, &
     0.50661630334978742377_dp, &
     0.12251482265544137787_dp, &
    -0.50661630334978742377_dp, &
     0.50661630334978742377_dp, &
     0.12251482265544137787_dp /), &
    (/ 3, 8 /) )

    w(1:order) = w_save(1:order)
    xyz(1:3,1:order) = xyz_save(1:3,1:order)
  end subroutine pyra_unit_o08

  pure subroutine pyra_unit_o08b ( w, xyz ) &
        bind(C, name="pyra_unit_o08b")

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
  !    Output, real(dp) W(8), the weights.
  !
  !    Output, real(dp) XYZ(3,8), the abscissas.
  !

    integer(ip), parameter :: order = 1

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(8) = (/ &
     0.16438287736328777572_dp, &
     0.16438287736328777572_dp, &
     0.16438287736328777572_dp, &
     0.16438287736328777572_dp, &
     0.085617122636712224276_dp, &
     0.085617122636712224276_dp, &
     0.085617122636712224276_dp, &
     0.085617122636712224276_dp /)
    real(dp), intent(out) :: xyz(3,order)
    real(dp) :: xyz_save(3,8) = reshape ( (/ &
    -0.51197009372656270107_dp, &
    -0.51197009372656270107_dp, &
     0.11024490204163285720_dp, &
     0.51197009372656270107_dp, &
    -0.51197009372656270107_dp, &
     0.11024490204163285720_dp, &
     0.51197009372656270107_dp, &
     0.51197009372656270107_dp, &
     0.11024490204163285720_dp, &
    -0.51197009372656270107_dp, &
     0.51197009372656270107_dp, &
     0.11024490204163285720_dp, &
    -0.28415447557052037456_dp, &
    -0.28415447557052037456_dp, &
     0.518326526529795714229_dp, &
     0.28415447557052037456_dp, &
    -0.28415447557052037456_dp, &
     0.518326526529795714229_dp, &
     0.28415447557052037456_dp, &
     0.28415447557052037456_dp, &
     0.518326526529795714229_dp, &
    -0.28415447557052037456_dp, &
     0.28415447557052037456_dp, &
     0.518326526529795714229_dp /), &
    (/ 3, 8 /) )

    w(1:order) = w_save(1:order)
    xyz(1:3,1:order) = xyz_save(1:3,1:order)
  end subroutine pyra_unit_o08b

  pure subroutine pyra_unit_o09 ( w, xyz ) &
        bind(C, name="pyra_unit_o09")

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
  !    Output, real(dp) W(9), the weights.
  !
  !    Output, real(dp) XYZ(3,9), the abscissas.
  !

    integer(ip), parameter :: order = 9

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(9) = (/ &
     0.13073389672275944791_dp, &
     0.13073389672275944791_dp, &
     0.13073389672275944791_dp, &
     0.13073389672275944791_dp, &
     0.10989110327724055209_dp, &
     0.10989110327724055209_dp, &
     0.10989110327724055209_dp, &
     0.10989110327724055209_dp, &
     0.03750000000000000000_dp /)
    real(dp), intent(out) :: xyz(3,order)
    real(dp) :: xyz_save(3,9) = reshape ( (/ &
    -0.52966422253852215131_dp, &
    -0.52966422253852215131_dp, &
     0.08176876558246862335_dp, &
     0.52966422253852215131_dp, &
    -0.52966422253852215131_dp, &
     0.08176876558246862335_dp, &
     0.52966422253852215131_dp, &
     0.52966422253852215131_dp, &
     0.08176876558246862335_dp, &
    -0.52966422253852215131_dp, &
     0.52966422253852215131_dp, &
     0.08176876558246862335_dp, &
    -0.34819753825720418039_dp, &
    -0.34819753825720418039_dp, &
     0.400374091560388519511_dp, &
     0.34819753825720418039_dp, &
    -0.34819753825720418039_dp, &
     0.400374091560388519511_dp, &
     0.34819753825720418039_dp, &
     0.34819753825720418039_dp, &
     0.400374091560388519511_dp, &
    -0.34819753825720418039_dp, &
     0.34819753825720418039_dp, &
     0.400374091560388519511_dp, &
     0.00000000000000000000_dp, &
     0.00000000000000000000_dp, &
     0.83333333333333333333_dp /), &
    (/ 3, 9 /) )

    w(1:order) = w_save(1:order)
    xyz(1:3,1:order) = xyz_save(1:3,1:order)
  end subroutine pyra_unit_o09

  pure subroutine pyra_unit_o13 ( w, xyz ) &
        bind(C, name="pyra_unit_o13")

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
  !    Output, real(dp) W(13), the weights.
  !
  !    Output, real(dp) XYZ(3,13), the abscissas.
  !

    integer(ip), parameter :: order = 13

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(13) = (/ &
     0.063061594202898550725_dp, &
     0.063061594202898550725_dp, &
     0.063061594202898550725_dp, &
     0.063061594202898550725_dp, &
     0.042101946815575556199_dp, &
     0.042101946815575556199_dp, &
     0.042101946815575556199_dp, &
     0.042101946815575556199_dp, &
     0.13172030707666776585_dp, &
     0.13172030707666776585_dp, &
     0.13172030707666776585_dp, &
     0.13172030707666776585_dp, &
     0.05246460761943250889_dp /)
    real(dp), intent(out) :: xyz(3,order)
    real(dp) :: xyz_save(3,13) = reshape ( (/ &
    -0.38510399211870384331_dp, &
    -0.38510399211870384331_dp, &
    0.428571428571428571429_dp, &
     0.38510399211870384331_dp, &
    -0.38510399211870384331_dp, &
    0.428571428571428571429_dp, &
     0.38510399211870384331_dp, &
     0.38510399211870384331_dp, &
    0.428571428571428571429_dp, &
    -0.38510399211870384331_dp, &
     0.38510399211870384331_dp, &
    0.428571428571428571429_dp, &
    -0.40345831960728204766_dp, &
     0.00000000000000000000_dp, &
    0.33928571428571428571_dp,  &
     0.40345831960728204766_dp, &
     0.00000000000000000000_dp, &
    0.33928571428571428571_dp,  &
     0.00000000000000000000_dp, &
    -0.40345831960728204766_dp, &
    0.33928571428571428571_dp,  &
     0.00000000000000000000_dp, &
     0.40345831960728204766_dp, &
    0.33928571428571428571_dp,  &
    -0.53157877436961973359_dp, &
    -0.53157877436961973359_dp, &
    0.08496732026143790850_dp,  &
     0.53157877436961973359_dp, &
    -0.53157877436961973359_dp, &
    0.08496732026143790850_dp,  &
     0.53157877436961973359_dp, &
     0.53157877436961973359_dp, &
    0.08496732026143790850_dp,  &
    -0.53157877436961973359_dp, &
     0.53157877436961973359_dp, &
    0.08496732026143790850_dp,  &
     0.00000000000000000000_dp, &
     0.00000000000000000000_dp, &
    0.76219701803768503595_dp /), &
    (/ 3, 13 /) )

    w(1:order) = w_save(1:order)
    xyz(1:3,1:order) = xyz_save(1:3,1:order)
  end subroutine pyra_unit_o13

  pure subroutine pyra_unit_o18 ( w, xyz ) &
        bind(C, name="pyra_unit_o18")

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
  !    Output, real(dp) W(18), the weights.
  !
  !    Output, real(dp) XYZ(3,18), the abscissas.
  !

    integer(ip), parameter :: order = 18

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(18) = (/ &
     0.023330065296255886709_dp, &
     0.037328104474009418735_dp, &
     0.023330065296255886709_dp, &
     0.037328104474009418735_dp, &
     0.059724967158415069975_dp, &
     0.037328104474009418735_dp, &
     0.023330065296255886709_dp, &
     0.037328104474009418735_dp, &
     0.023330065296255886709_dp, &
     0.05383042853090460712_dp, &
     0.08612868564944737139_dp, &
     0.05383042853090460712_dp, &
     0.08612868564944737139_dp, &
     0.13780589703911579422_dp, &
     0.08612868564944737139_dp, &
     0.05383042853090460712_dp, &
     0.08612868564944737139_dp, &
     0.05383042853090460712_dp /)
    real(dp), intent(out) :: xyz(3,order)
    real(dp) :: xyz_save(3,18) = reshape ( (/ &
    -0.35309846330877704481_dp, &
    -0.35309846330877704481_dp, &
    0.544151844011225288800_dp, &
     0.00000000000000000000_dp, &
    -0.35309846330877704481_dp, &
    0.544151844011225288800_dp, &
     0.35309846330877704481_dp, &
    -0.35309846330877704481_dp, &
    0.544151844011225288800_dp, &
    -0.35309846330877704481_dp, &
     0.00000000000000000000_dp, &
    0.544151844011225288800_dp, &
     0.00000000000000000000_dp, &
     0.00000000000000000000_dp, &
    0.544151844011225288800_dp, &
     0.35309846330877704481_dp, &
     0.00000000000000000000_dp, &
    0.544151844011225288800_dp, &
    -0.35309846330877704481_dp, &
     0.35309846330877704481_dp, &
    0.544151844011225288800_dp, &
     0.00000000000000000000_dp, &
     0.35309846330877704481_dp, &
    0.544151844011225288800_dp, &
     0.35309846330877704481_dp, &
     0.35309846330877704481_dp, &
    0.544151844011225288800_dp, &
    -0.67969709567986745790_dp, &
    -0.67969709567986745790_dp, &
    0.12251482265544137787_dp, &
     0.00000000000000000000_dp, &
    -0.67969709567986745790_dp, &
    0.12251482265544137787_dp, &
     0.67969709567986745790_dp, &
    -0.67969709567986745790_dp, &
    0.12251482265544137787_dp, &
    -0.67969709567986745790_dp, &
     0.00000000000000000000_dp, &
    0.12251482265544137787_dp, &
     0.00000000000000000000_dp, &
     0.00000000000000000000_dp, &
    0.12251482265544137787_dp, &
     0.67969709567986745790_dp, &
     0.00000000000000000000_dp, &
    0.12251482265544137787_dp, &
    -0.67969709567986745790_dp, &
     0.67969709567986745790_dp, &
    0.12251482265544137787_dp, &
     0.00000000000000000000_dp, &
     0.67969709567986745790_dp, &
    0.12251482265544137787_dp, &
     0.67969709567986745790_dp, &
     0.67969709567986745790_dp, &
    0.12251482265544137787_dp /), &
    (/ 3, 18 /) )

    w(1:order) = w_save(1:order)
    xyz(1:3,1:order) = xyz_save(1:3,1:order)
  end subroutine pyra_unit_o18

  pure subroutine pyra_unit_o27 ( w, xyz ) &
        bind(C, name="pyra_unit_o27")

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
  !    Output, real(dp) W(27), the weights.
  !
  !    Output, real(dp) XYZ(3,27), the abscissas.
  !

    integer(ip), parameter :: order = 27

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(27) = (/ &
     0.036374157653908938268_dp, &
     0.05819865224625430123_dp, &
     0.036374157653908938268_dp, &
     0.05819865224625430123_dp, &
     0.09311784359400688197_dp, &
     0.05819865224625430123_dp, &
     0.036374157653908938268_dp, &
     0.05819865224625430123_dp, &
     0.036374157653908938268_dp, &
     0.033853303069413431019_dp, &
     0.054165284911061489631_dp, &
     0.033853303069413431019_dp, &
     0.054165284911061489631_dp, &
     0.08666445585769838341_dp, &
     0.054165284911061489631_dp, &
     0.033853303069413431019_dp, &
     0.054165284911061489631_dp, &
     0.033853303069413431019_dp, &
     0.006933033103838124540_dp, &
     0.011092852966140999264_dp, &
     0.006933033103838124540_dp, &
     0.011092852966140999264_dp, &
     0.017748564745825598822_dp, &
     0.011092852966140999264_dp, &
     0.006933033103838124540_dp, &
     0.011092852966140999264_dp, &
     0.006933033103838124540_dp /)
    real(dp), intent(out) :: xyz(3,order)
    real(dp) :: xyz_save(3,27) = reshape ( (/ &
    -0.7180557413198889387_dp, &
     -0.7180557413198889387_dp, &
     0.07299402407314973216_dp, &
     0.00000000000000000000_dp, &
    -0.7180557413198889387_dp, &
     0.07299402407314973216_dp, &
     0.7180557413198889387_dp, &
     -0.7180557413198889387_dp, &
     0.07299402407314973216_dp, &
    -0.7180557413198889387_dp, &
      0.00000000000000000000_dp, &
    0.07299402407314973216_dp, &
     0.00000000000000000000_dp, &
     0.00000000000000000000_dp, &
    0.07299402407314973216_dp, &
     0.7180557413198889387_dp, &
      0.00000000000000000000_dp, &
    0.07299402407314973216_dp, &
    -0.7180557413198889387_dp, &
      0.7180557413198889387_dp, &
     0.07299402407314973216_dp, &
     0.00000000000000000000_dp, &
     0.7180557413198889387_dp, &
     0.07299402407314973216_dp, &
     0.7180557413198889387_dp, &
      0.7180557413198889387_dp, &
     0.07299402407314973216_dp, &
    -0.50580870785392503961_dp, &
    -0.50580870785392503961_dp, &
    0.34700376603835188472_dp, &
     0.00000000000000000000_dp, &
    -0.50580870785392503961_dp, &
    0.34700376603835188472_dp, &
     0.50580870785392503961_dp, &
    -0.50580870785392503961_dp, &
    0.34700376603835188472_dp, &
    -0.50580870785392503961_dp, &
     0.00000000000000000000_dp, &
    0.34700376603835188472_dp, &
     0.00000000000000000000_dp, &
     0.00000000000000000000_dp, &
    0.34700376603835188472_dp, &
     0.50580870785392503961_dp, &
     0.00000000000000000000_dp, &
    0.34700376603835188472_dp, &
    -0.50580870785392503961_dp, &
     0.50580870785392503961_dp, &
    0.34700376603835188472_dp, &
     0.00000000000000000000_dp, &
     0.50580870785392503961_dp, &
    0.34700376603835188472_dp, &
     0.50580870785392503961_dp, &
     0.50580870785392503961_dp, &
    0.34700376603835188472_dp, &
    -0.22850430565396735360_dp, &
    -0.22850430565396735360_dp, &
    0.70500220988849838312_dp, &
     0.00000000000000000000_dp, &
    -0.22850430565396735360_dp, &
    0.70500220988849838312_dp, &
     0.22850430565396735360_dp, &
    -0.22850430565396735360_dp, &
    0.70500220988849838312_dp, &
    -0.22850430565396735360_dp, &
     0.00000000000000000000_dp, &
    0.70500220988849838312_dp, &
     0.00000000000000000000_dp, &
     0.00000000000000000000_dp, &
    0.70500220988849838312_dp, &
     0.22850430565396735360_dp, &
     0.00000000000000000000_dp, &
    0.70500220988849838312_dp, &
    -0.22850430565396735360_dp, &
     0.22850430565396735360_dp, &
    0.70500220988849838312_dp, &
     0.00000000000000000000_dp, &
     0.22850430565396735360_dp, &
    0.70500220988849838312_dp, &
     0.22850430565396735360_dp, &
     0.22850430565396735360_dp, &
    0.70500220988849838312_dp /), &
    (/ 3, 27 /) )

    w(1:order) = w_save(1:order)
    xyz(1:3,1:order) = xyz_save(1:3,1:order)
  end subroutine pyra_unit_o27

  pure subroutine pyra_unit_o48 ( w, xyz ) &
        bind(C, name="pyra_unit_o48")

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
  !    Output, real(dp) W(48), the weights.
  !
  !    Output, real(dp) XYZ(3,48), the abscissas.
  !

    integer(ip), parameter :: order = 48

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(48) = (/ &
    2.01241939442682455e-002_dp, &
    2.01241939442682455e-002_dp, &
    2.01241939442682455e-002_dp, &
    2.01241939442682455e-002_dp, &
    2.60351137043010779e-002_dp, &
    2.60351137043010779e-002_dp, &
    2.60351137043010779e-002_dp, &
    2.60351137043010779e-002_dp, &
    1.24557795239745531e-002_dp, &
    1.24557795239745531e-002_dp, &
    1.24557795239745531e-002_dp, &
    1.24557795239745531e-002_dp, &
    1.87873998794808156e-003_dp, &
    1.87873998794808156e-003_dp, &
    1.87873998794808156e-003_dp, &
    1.87873998794808156e-003_dp, &
    4.32957927807745280e-002_dp, &
    4.32957927807745280e-002_dp, &
    4.32957927807745280e-002_dp, &
    4.32957927807745280e-002_dp, &
    1.97463249834127288e-002_dp, &
    1.97463249834127288e-002_dp, &
    1.97463249834127288e-002_dp, &
    1.97463249834127288e-002_dp, &
    5.60127223523590526e-002_dp, &
    5.60127223523590526e-002_dp, &
    5.60127223523590526e-002_dp, &
    5.60127223523590526e-002_dp, &
    2.55462562927473852e-002_dp, &
    2.55462562927473852e-002_dp, &
    2.55462562927473852e-002_dp, &
    2.55462562927473852e-002_dp, &
    2.67977366291788643e-002_dp, &
    2.67977366291788643e-002_dp, &
    2.67977366291788643e-002_dp, &
    2.67977366291788643e-002_dp, &
    1.22218992265373354e-002_dp, &
    1.22218992265373354e-002_dp, &
    1.22218992265373354e-002_dp, &
    1.22218992265373354e-002_dp, &
    4.04197740453215038e-003_dp, &
    4.04197740453215038e-003_dp, &
    4.04197740453215038e-003_dp, &
    4.04197740453215038e-003_dp, &
    1.84346316995826843e-003_dp, &
    1.84346316995826843e-003_dp, &
    1.84346316995826843e-003_dp, &
    1.84346316995826843e-003_dp /)
    real(dp), intent(out) :: xyz(3,order)
    real(dp) :: xyz_save(3,48) = reshape ( (/ &
    0.88091731624450909_dp, &
     0.0000000000000000_dp, &
     4.85005494469969989e-02_dp, &
   -0.88091731624450909_dp, &
     0.0000000000000000_dp, &
     4.85005494469969989e-02_dp, &
     0.0000000000000000_dp, &
     0.88091731624450909_dp, &
    4.85005494469969989e-02_dp, &
     0.0000000000000000_dp, &
    -0.88091731624450909_dp, &
    4.85005494469969989e-02_dp, &
    0.70491874112648223_dp, &
     0.0000000000000000_dp, &
     0.23860073755186201_dp, &
   -0.70491874112648223_dp, &
     0.0000000000000000_dp, &
     0.23860073755186201_dp, &
     0.0000000000000000_dp, &
     0.70491874112648223_dp, &
    0.23860073755186201_dp, &
     0.0000000000000000_dp, &
    -0.70491874112648223_dp, &
    0.23860073755186201_dp, &
    0.44712732143189760_dp, &
     0.0000000000000000_dp, &
     0.51704729510436798_dp, &
   -0.44712732143189760_dp, &
     0.0000000000000000_dp, &
     0.51704729510436798_dp, &
     0.0000000000000000_dp, &
     0.44712732143189760_dp, &
    0.51704729510436798_dp, &
     0.0000000000000000_dp, &
    -0.44712732143189760_dp, &
    0.51704729510436798_dp, &
    0.18900486065123448_dp, &
     0.0000000000000000_dp, &
     0.79585141789677305_dp, &
   -0.18900486065123448_dp, &
     0.0000000000000000_dp, &
     0.79585141789677305_dp, &
     0.0000000000000000_dp, &
     0.18900486065123448_dp, &
    0.79585141789677305_dp, &
     0.0000000000000000_dp, &
    -0.18900486065123448_dp, &
    0.79585141789677305_dp, &
    0.36209733410322176_dp, &
     0.36209733410322176_dp, &
    4.85005494469969989e-02_dp, &
   -0.36209733410322176_dp, &
     0.36209733410322176_dp, &
    4.85005494469969989e-02_dp, &
   -0.36209733410322176_dp, &
    -0.36209733410322176_dp, &
    4.85005494469969989e-02_dp, &
    0.36209733410322176_dp, &
    -0.36209733410322176_dp, &
    4.85005494469969989e-02_dp, &
    0.76688932060387538_dp, &
     0.76688932060387538_dp, &
    4.85005494469969989e-02_dp, &
   -0.76688932060387538_dp, &
     0.76688932060387538_dp, &
    4.85005494469969989e-02_dp, &
   -0.76688932060387538_dp, &
    -0.76688932060387538_dp, &
    4.85005494469969989e-02_dp, &
    0.76688932060387538_dp, &
    -0.76688932060387538_dp, &
    4.85005494469969989e-02_dp, &
    0.28975386476618070_dp, &
     0.28975386476618070_dp, &
    0.23860073755186201_dp, &
   -0.28975386476618070_dp, &
     0.28975386476618070_dp, &
    0.23860073755186201_dp, &
   -0.28975386476618070_dp, &
    -0.28975386476618070_dp, &
    0.23860073755186201_dp, &
    0.28975386476618070_dp, &
    -0.28975386476618070_dp, &
    0.23860073755186201_dp, &
    0.61367241226233160_dp, &
     0.61367241226233160_dp, &
    0.23860073755186201_dp, &
   -0.61367241226233160_dp, &
     0.61367241226233160_dp, &
    0.23860073755186201_dp, &
   -0.61367241226233160_dp, &
    -0.61367241226233160_dp, &
    0.23860073755186201_dp, &
    0.61367241226233160_dp, &
    -0.61367241226233160_dp, &
    0.23860073755186201_dp, &
    0.18378979287798017_dp, &
     0.18378979287798017_dp, &
    0.51704729510436798_dp, &
   -0.18378979287798017_dp, &
     0.18378979287798017_dp, &
    0.51704729510436798_dp, &
   -0.18378979287798017_dp, &
    -0.18378979287798017_dp, &
    0.51704729510436798_dp, &
    0.18378979287798017_dp, &
    -0.18378979287798017_dp, &
    0.51704729510436798_dp, &
    0.38925011625173161_dp, &
     0.38925011625173161_dp, &
    0.51704729510436798_dp, &
   -0.38925011625173161_dp, &
     0.38925011625173161_dp, &
    0.51704729510436798_dp, &
   -0.38925011625173161_dp, &
    -0.38925011625173161_dp, &
    0.51704729510436798_dp, &
    0.38925011625173161_dp, &
    -0.38925011625173161_dp, &
    0.51704729510436798_dp, &
    7.76896479525748113e-02_dp, &
     7.76896479525748113e-02_dp, &
    0.79585141789677305_dp, &
   -7.76896479525748113e-02_dp, &
     7.76896479525748113e-02_dp, &
    0.79585141789677305_dp, &
   -7.76896479525748113e-02_dp, &
    -7.76896479525748113e-02_dp, &
    0.79585141789677305_dp, &
    7.76896479525748113e-02_dp, &
    -7.76896479525748113e-02_dp, &
    0.79585141789677305_dp, &
    0.16453962988669860_dp, &
     0.16453962988669860_dp, &
    0.79585141789677305_dp, &
   -0.16453962988669860_dp, &
     0.16453962988669860_dp, &
    0.79585141789677305_dp, &
   -0.16453962988669860_dp, &
    -0.16453962988669860_dp, &
    0.79585141789677305_dp, &
    0.16453962988669860_dp, &
    -0.16453962988669860_dp, &
    0.79585141789677305_dp /), &
    (/ 3, 48 /) )

    w(1:order) = w_save(1:order)
    xyz(1:3,1:order) = xyz_save(1:3,1:order)
  end subroutine pyra_unit_o48

  subroutine pyra_unit_quad_test ( degree_max ) &
        bind(C, name="pyra_unit_quad_test")

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
  !    Input, integer(ip) DEGREE_MAX, the maximum total degree of the
  !    monomials to check.
  !

    integer(ip), parameter :: dim_num = 3

    integer(ip), intent(in), value :: degree_max
    integer(ip) :: expon(dim_num)
    integer(ip) :: h
    logical :: more
    integer(ip) :: order
    real(dp) :: quad
    integer(ip) :: t
    real(dp) :: pyra_unit_volume
    real(dp), allocatable :: v(:)
    real(dp), allocatable :: w(:)
    real(dp), allocatable :: xyz(:,:)

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
  end subroutine pyra_unit_quad_test

  pure function pyra_unit_volume ( ) &
        bind(C, name="pyra_unit_volume")

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
  !    Output, real(dp) PYRA_UNIT_VOLUME, the volume.
  !

    real(dp) :: pyra_unit_volume

    pyra_unit_volume = 4.0_dp / 3.0_dp
  end function pyra_unit_volume

  subroutine quad_unit_monomial ( expon, value ) &
        bind(C, name="quad_unit_monomial")

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
  !    Input, integer(ip) EXPON(2), the exponents.
  !
  !    Output, real(dp) VALUE, the integral of the monomial.
  !

    integer(ip), intent(in) :: expon(2)
    integer(ip) :: i
    real(dp), intent(out) :: value

    value = 1.0_dp

    do i = 1, 2

      if ( mod ( expon(i), 2 ) == 1 ) then
        value = 0.0_dp
      else if ( expon(i) == -1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'QUAD_UNIT_MONOMIAL - Fatal error!'
        write ( *, '(a)' ) '  Exponent of -1 encountered.'
        stop 1
      else
        value = value * 2.0_dp / real ( expon(i) + 1, dp)
      end if

    end do
  end subroutine quad_unit_monomial

  subroutine quad_unit_monomial_test ( degree_max ) &
        bind(C, name="quad_unit_monomial_test")

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
  !    Input, integer(ip) DEGREE_MAX, the maximum total degree of the
  !    monomials to check.
  !

    integer(ip) :: alpha
    integer(ip) :: beta
    integer(ip), intent(in), value :: degree_max
    integer(ip) :: expon(2)
    real(dp) :: quad_unit_volume
    real(dp) :: value

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
  end subroutine quad_unit_monomial_test

  subroutine quad_unit_quad_test ( degree_max ) &
        bind(C, name="quad_unit_quad_test")

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
  !    Input, integer(ip) DEGREE_MAX, the maximum total degree of the
  !    monomials to check.
  !

    integer(ip), parameter :: dim_num = 2

    integer(ip), intent(in), value :: degree_max
    integer(ip) :: expon(dim_num)
    integer(ip) :: h
    integer(ip) :: k
    logical :: more
    integer(ip) :: order
    integer(ip) :: order_1d(dim_num)
    real(dp) :: quad
    real(dp) :: quad_unit_volume
    integer(ip) :: t
    real(dp), allocatable :: v(:)
    real(dp), allocatable :: w(:)
    real(dp), allocatable :: xy(:,:)

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
  end subroutine quad_unit_quad_test

  subroutine quad_unit_rule ( order_1d, w, xy ) &
        bind(C, name="quad_unit_rule")

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
  !    Input, integer(ip) ORDER_1D(2), the order of the rule in
  !    each dimension.  1 <= ORDER_1D(I) <= 5.
  !
  !    Output, real(dp) W(ORDER_1D(1)*ORDER_1D(2)), the weights.
  !
  !    Output, real(dp) XY(2,ORDER_1D(1)*ORDER_1D(2)), the abscissas.
  !

    integer(ip), parameter :: dim_num = 2

    integer(ip) :: dim
    integer(ip) :: order
    integer(ip), intent(in) :: order_1d(2)
    real(dp), intent(out) :: w(order_1D(1)*order_1D(2))
    real(dp), allocatable :: w_1d(:)
    real(dp), allocatable :: x_1d(:)
    real(dp), intent(out) :: xy(2,order_1D(1)*order_1D(2))

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
  end subroutine quad_unit_rule

  pure function quad_unit_volume ( ) &
        bind(C, name="quad_unit_volume")

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
  !    Output, real(dp) QUAD_UNIT_VOLUME, the volume.
  !

    real(dp) :: quad_unit_volume

    quad_unit_volume = 4.0_dp
  end function quad_unit_volume

  pure function r8_choose ( n, k ) &
        bind(C, name="r8_choose")

  !*****************************************************************************80
  !
  !! R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
  !
  !  Discussion:
  !
  !    The value is calculated in such a way as to avoid overflow and
  !    roundoff.  The calculation is done in R8 arithmetic.
  !
  !    The formula used is:
  !
  !      C(N,K) = N! / ( K! * (N-K)! )
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
  !    ML Wolfson, HV Wright,
  !    Algorithm 160:
  !    Combinatorial of M Things Taken N at a Time,
  !    Communications of the ACM,
  !    Volume 6, Number 4, April 1963, page 161.
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, K, are the values of N and K.
  !
  !    Output, real(dp) R8_CHOOSE, the number of combinations of N
  !    things taken K at a time.
  !

    integer(ip) :: i
    integer(ip), intent(in), value :: k
    integer(ip) :: mn
    integer(ip) :: mx
    integer(ip), intent(in), value :: n
    real(dp) :: r8_choose
    real(dp) :: value

    mn = min ( k, n - k )

    if ( mn < 0 ) then

      value = 0.0_dp

    else if ( mn == 0 ) then

      value = 1.0_dp

    else

      mx = max ( k, n - k )
      value = real ( mx + 1, dp)

      do i = 2, mn
        value = ( value * real ( mx + i, dp) ) / real ( i, dp)
      end do

    end if

    r8_choose = value
  end function r8_choose

  pure function r8_mop ( i ) &
        bind(C, name="r8_mop")

  !*****************************************************************************80
  !
  !! R8_MOP returns the I-th power of -1 as an R8.
  !
  !  Discussion:
  !
  !    An R8 is a real(dp) value.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 November 2007
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) I, the power of -1.
  !
  !    Output, real(dp) R8_MOP, the I-th power of -1.
  !

    integer(ip), intent(in), value :: i
    real(dp) :: r8_mop

    if ( mod ( i, 2 ) == 0 ) then
      r8_mop = + 1.0_dp
    else
      r8_mop = - 1.0_dp
    end if
  end function r8_mop

  pure function r8mat_det_4d ( a ) &
        bind(C, name="r8mat_det_4d")

  !*****************************************************************************80
  !
  !! R8MAT_DET_4D computes the determinant of a 4 by 4 matrix.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    16 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real(dp) A(4,4), the matrix whose determinant is desired.
  !
  !    Output, real(dp) R8MAT_DET_4D, the determinant of the matrix.
  !

    real(dp), intent(in) :: a(4,4)
    real(dp) :: r8mat_det_4d

    r8mat_det_4d = &
        a(1,1) * ( &
          a(2,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
        - a(2,3) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
        + a(2,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) ) &
      - a(1,2) * ( &
          a(2,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
        - a(2,3) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
        + a(2,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) ) &
      + a(1,3) * ( &
          a(2,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
        - a(2,2) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
        + a(2,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) ) &
      - a(1,4) * ( &
          a(2,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
        - a(2,2) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
        + a(2,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) )
  end function r8mat_det_4d

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
  !    17 April 2009
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
  !    Local, integer START, the first location of a block of values to set.
  !
  !    Local, integer CONTIG, the number of consecutive values to set.
  !
  !    Local, integer SKIP, the distance from the current value of START
  !    to the next location of a block of values to set.
  !
  !    Local, integer REP, the number of blocks of values to set.
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

  subroutine r8vec_direct_product2 ( factor_index, factor_order, factor_value, &
    factor_num, point_num, w ) &
        bind(C, name="r8vec_direct_product2")

  !*****************************************************************************80
  !
  !! R8VEC_DIRECT_PRODUCT2 creates a direct product of R8VEC's.
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
  !    This routine carries out the task involving the weights W.
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
  !      W(1:4) = ( 2, 3, 5, 7 )
  !
  !    Rule 2:
  !      Order = 3
  !      W(1:3) = ( 11, 13, 17 )
  !
  !    Rule 3:
  !      Order = 2
  !      W(1:2) = ( 19, 23 )
  !
  !    Product Rule:
  !      Order = 24
  !      W(1:24) =
  !        ( 2 * 11 * 19 )
  !        ( 3 * 11 * 19 )
  !        ( 4 * 11 * 19 )
  !        ( 7 * 11 * 19 )
  !        ( 2 * 13 * 19 )
  !        ( 3 * 13 * 19 )
  !        ( 5 * 13 * 19 )
  !        ( 7 * 13 * 19 )
  !        ( 2 * 17 * 19 )
  !        ( 3 * 17 * 19 )
  !        ( 5 * 17 * 19 )
  !        ( 7 * 17 * 19 )
  !        ( 2 * 11 * 23 )
  !        ( 3 * 11 * 23 )
  !        ( 5 * 11 * 23 )
  !        ( 7 * 11 * 23 )
  !        ( 2 * 13 * 23 )
  !        ( 3 * 13 * 23 )
  !        ( 5 * 13 * 23 )
  !        ( 7 * 13 * 23 )
  !        ( 2 * 17 * 23 )
  !        ( 3 * 17 * 23 )
  !        ( 5 * 17 * 23 )
  !        ( 7 * 17 * 23 )
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
  !    Input/output, real(dp) W(POINT_NUM), the elements of the
  !    direct product, which are built up gradually.
  !
  !  Local Parameters:
  !
  !    Local, integer(ip) START, the first location of a block of values
  !    to set.
  !
  !    Local, integer(ip) CONTIG, the number of consecutive values 
  !    to set.
  !
  !    Local, integer SKIP, the distance from the current value of START
  !    to the next location of a block of values to set.
  !
  !    Local, integer REP, the number of blocks of values to set.
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
    real(dp), intent(inout) :: w(point_num)

    if ( factor_index == 1 ) then
      contig = 1
      skip = 1
      rep = point_num
      w(1:point_num) = 1.0_dp
    end if

    rep = rep / factor_order
    skip = skip * factor_order

    do j = 1, factor_order

      start = 1 + ( j - 1 ) * contig

      do k = 1, rep
        w(start:start+contig-1) = w(start:start+contig-1) * factor_value(j)
        start = start + skip
      end do

    end do

    contig = contig * factor_order
  end subroutine r8vec_direct_product2

  subroutine subcomp_next ( n, k, a, more, h, t ) &
        bind(C, name="subcomp_next")

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
  !    Input, integer(ip) N, the integer whose subcompositions
  !    are desired.
  !
  !    Input, integer(ip) K, the number of parts in the subcomposition.
  !
  !    Input/output, integer(ip) A(K), the parts of the subcomposition.
  !
  !    Input/output, logical MORE, set by the user to start the computation,
  !    and by the routine to terminate it.
  !
  !    Input/output, integer H, T, two internal parameters needed for the
  !    computation.  The user should allocate space for these in the calling
  !    program, include them in the calling sequence, but never alter them!
  !

    integer(ip), intent(in), value :: k

    integer(ip), intent(inout) :: a(k)
    integer(ip), intent(inout) :: h
    logical, intent(inout) :: more
    logical, save :: more2 = .false.
    integer(ip), intent(in), value :: n
    integer(ip), save :: n2 = 0
    integer(ip), intent(inout) :: t
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
  end subroutine subcomp_next

  pure subroutine tetr_unit_monomial ( expon, value ) &
        bind(C, name="tetr_unit_monomial")

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
  !    Input, integer(ip) EXPON(3), the exponents.
  !
  !    Output, real(dp) VALUE, the integral of the monomial.
  !

    integer(ip), intent(in) :: expon(3)
    integer(ip) :: i
    integer(ip) :: k
    real(dp), intent(out) :: value
  !
  !  The first computation ends with VALUE = 1.0;
  !
    value = 1.0_dp
  !
  !  The first loop simply calculates 1, so we short circuit it.
  !
  ! k = 0
  !
  ! do i = 1, expon(1)
  !   k = k + 1
  !   value = value * real ( i, dp) / real ( k, dp)
  ! end do

    k = expon(1)
    do i = 1, expon(2)
      k = k + 1
      value = value * real ( i, dp) / real ( k, dp)
    end do

    do i = 1, expon(3)
      k = k + 1
      value = value * real ( i, dp) / real ( k, dp)
    end do

    k = k + 1
    value = value / real ( k, dp)

    k = k + 1
    value = value / real ( k, dp)

    k = k + 1
    value = value / real ( k, dp)
  end subroutine tetr_unit_monomial

  subroutine tetr_unit_monomial_test ( degree_max ) &
        bind(C, name="tetr_unit_monomial_test")

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
  !    Input, integer(ip) DEGREE_MAX, the maximum total degree of the
  !    monomials to check.
  !

    integer(ip) :: alpha
    integer(ip) :: beta
    integer(ip), intent(in), value :: degree_max
    integer(ip) :: expon(3)
    integer(ip) :: gamma
    real(dp) :: tetr_unit_volume
    real(dp) :: value

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
  end subroutine tetr_unit_monomial_test

  pure subroutine tetr_unit_o01 ( w, xyz ) &
        bind(C, name="tetr_unit_o01")

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
  !    Output, real(dp) W(1), the weights.
  !
  !    Output, real(dp) XYZ(3,1), the abscissas.
  !

    integer(ip), parameter :: order = 1

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(1) = (/ &
      1.0000000000000000000_dp /)
    real(dp), intent(out) :: xyz(3,order)
    real(dp) :: xyz_save(3,1) = reshape ( (/ &
      0.25000000000000000000_dp,  0.25000000000000000000_dp,  &
      0.25000000000000000000_dp /), &
    (/ 3, 1 /) )

    w(1:order) = w_save(1:order)
    xyz(1:3,1:order) = xyz_save(1:3,1:order)
  end subroutine tetr_unit_o01

  pure subroutine tetr_unit_o04 ( w, xyz ) &
        bind(C, name="tetr_unit_o04")

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
  !    Output, real(dp) W(4), the weights.
  !
  !    Output, real(dp) XYZ(3,4), the abscissas.
  !

    integer(ip), parameter :: order = 4

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(4) = (/ &
      0.25000000000000000000_dp, &
      0.25000000000000000000_dp, &
      0.25000000000000000000_dp, &
      0.25000000000000000000_dp /)
    real(dp), intent(out) :: xyz(3,order)
    real(dp) :: xyz_save(3,4) = reshape ( (/ &
      0.58541019662496845446_dp, &
      0.13819660112501051518_dp, &
      0.13819660112501051518_dp, &
      0.13819660112501051518_dp, &
      0.58541019662496845446_dp, &
      0.13819660112501051518_dp, &
      0.13819660112501051518_dp, &
      0.13819660112501051518_dp, &
      0.58541019662496845446_dp, &
      0.13819660112501051518_dp, &
      0.13819660112501051518_dp, &
      0.13819660112501051518_dp /), (/ 3, 4 /) )

    w(1:order) = w_save(1:order)
    xyz(1:3,1:order) = xyz_save(1:3,1:order)
  end subroutine tetr_unit_o04

  pure subroutine tetr_unit_o08 ( w, xyz ) &
        bind(C, name="tetr_unit_o08")

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
  !    Output, real(dp) W(8), the weights.
  !
  !    Output, real(dp) XYZ(3,8), the abscissas.
  !

    integer(ip), parameter :: order = 8

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(8) = (/ &
      0.13852796651186214232_dp, &
      0.13852796651186214232_dp, &
      0.13852796651186214232_dp, &
      0.13852796651186214232_dp, &
      0.11147203348813785768_dp, &
      0.11147203348813785768_dp, &
      0.11147203348813785768_dp, &
      0.11147203348813785768_dp /)
    real(dp), intent(out) :: xyz(3,order)
    real(dp) :: xyz_save(3,8) = reshape ( (/ &
      0.015835909865720057993_dp, &
      0.32805469671142664734_dp,  &
      0.32805469671142664734_dp,  &
      0.32805469671142664734_dp,  &
      0.015835909865720057993_dp, &
      0.32805469671142664734_dp,  &
      0.32805469671142664734_dp,  &
      0.32805469671142664734_dp,  &
      0.015835909865720057993_dp, &
      0.32805469671142664734_dp,  &
      0.32805469671142664734_dp,  &
      0.32805469671142664734_dp,  &
      0.67914317820120795168_dp,  &
      0.10695227393293068277_dp,  &
      0.10695227393293068277_dp,  &
      0.10695227393293068277_dp,  &
      0.67914317820120795168_dp,  &
      0.10695227393293068277_dp,  &
      0.10695227393293068277_dp,  &
      0.10695227393293068277_dp,  &
      0.67914317820120795168_dp,  &
      0.10695227393293068277_dp,  &
      0.10695227393293068277_dp,  &
      0.10695227393293068277_dp /), &
    (/ 3, 8 /) )

    w(1:order) = w_save(1:order)
    xyz(1:3,1:order) = xyz_save(1:3,1:order)
  end subroutine tetr_unit_o08

  pure subroutine tetr_unit_o08b ( w, xyz ) &
        bind(C, name="tetr_unit_o08b")

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
  !    Output, real(dp) W(8), the weights.
  !
  !    Output, real(dp) XYZ(3,8), the abscissas.
  !

    integer(ip), parameter :: order = 8

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(8) = (/ &
      0.025000000000000000000_dp, &
      0.025000000000000000000_dp, &
      0.025000000000000000000_dp, &
      0.025000000000000000000_dp, &
      0.22500000000000000000_dp, &
      0.22500000000000000000_dp, &
      0.22500000000000000000_dp, &
      0.22500000000000000000_dp /)
    real(dp), intent(out) :: xyz(3,order)
    real(dp) :: xyz_save(3,8) = reshape ( (/ &
      1.00000000000000000000_dp, &
      0.00000000000000000000_dp, &
      0.00000000000000000000_dp, &
      0.00000000000000000000_dp, &
      1.00000000000000000000_dp, &
      0.00000000000000000000_dp, &
      0.00000000000000000000_dp, &
      0.00000000000000000000_dp, &
      1.00000000000000000000_dp, &
      0.00000000000000000000_dp, &
      0.00000000000000000000_dp, &
      0.00000000000000000000_dp, &
      0.00000000000000000000_dp, &
      0.33333333333333333333_dp, &
      0.33333333333333333333_dp, &
      0.33333333333333333333_dp, &
      0.00000000000000000000_dp, &
      0.33333333333333333333_dp, &
      0.33333333333333333333_dp, &
      0.33333333333333333333_dp, &
      0.00000000000000000000_dp, &
      0.33333333333333333333_dp, &
      0.33333333333333333333_dp, &
      0.33333333333333333333_dp /), &
    (/ 3, 8 /) )

    w(1:order) = w_save(1:order)
    xyz(1:3,1:order) = xyz_save(1:3,1:order)
  end subroutine tetr_unit_o08b

  pure subroutine tetr_unit_o14 ( w, xyz ) &
        bind(C, name="tetr_unit_o14")

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
  !    Output, real(dp) W(ORDER), the weights.
  !
  !    Output, real(dp) XYZ(3,ORDER), the abscissas.
  !

    integer(ip), parameter :: order = 14

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(14) = (/ &
      0.073493043116361949544_dp, &
      0.073493043116361949544_dp, &
      0.073493043116361949544_dp, &
      0.073493043116361949544_dp, &
      0.11268792571801585080_dp, &
      0.11268792571801585080_dp, &
      0.11268792571801585080_dp, &
      0.11268792571801585080_dp, &
      0.042546020777081466438_dp, &
      0.042546020777081466438_dp, &
      0.042546020777081466438_dp, &
      0.042546020777081466438_dp, &
      0.042546020777081466438_dp, &
      0.042546020777081466438_dp /)
    real(dp), intent(out) :: xyz(3,order)
    real(dp) :: xyz_save(3,14) = reshape ( (/ &
      0.72179424906732632079_dp, &
      0.092735250310891226402_dp, &
      0.092735250310891226402_dp, &
      0.092735250310891226402_dp, &
      0.72179424906732632079_dp, &
      0.092735250310891226402_dp, &
      0.092735250310891226402_dp, &
      0.092735250310891226402_dp, &
      0.72179424906732632079_dp, &
      0.092735250310891226402_dp, &
      0.092735250310891226402_dp, &
      0.092735250310891226402_dp, &
      0.067342242210098170608_dp, &
      0.31088591926330060980_dp, &
      0.31088591926330060980_dp, &
      0.31088591926330060980_dp, &
      0.067342242210098170608_dp, &
      0.31088591926330060980_dp, &
      0.31088591926330060980_dp, &
      0.31088591926330060980_dp, &
      0.067342242210098170608_dp, &
      0.31088591926330060980_dp, &
      0.31088591926330060980_dp, &
      0.31088591926330060980_dp, &
      0.045503704125649649492_dp, &
      0.045503704125649649492_dp, &
      0.45449629587435035051_dp, &
      0.045503704125649649492_dp, &
      0.45449629587435035051_dp, &
      0.045503704125649649492_dp, &
      0.045503704125649649492_dp, &
      0.45449629587435035051_dp, &
      0.45449629587435035051_dp, &
      0.45449629587435035051_dp, &
      0.045503704125649649492_dp, &
      0.045503704125649649492_dp, &
      0.45449629587435035051_dp, &
      0.045503704125649649492_dp, &
      0.45449629587435035051_dp, &
      0.45449629587435035051_dp, &
      0.45449629587435035051_dp, &
      0.045503704125649649492_dp /), &
    (/ 3, 14 /) )

    w(1:order) = w_save(1:order)
    xyz(1:3,1:order) = xyz_save(1:3,1:order)
  end subroutine tetr_unit_o14

  pure subroutine tetr_unit_o14b ( w, xyz ) &
        bind(C, name="tetr_unit_o14b")

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
  !    Output, real(dp) W(14), the weights.
  !
  !    Output, real(dp) XYZ(3,14), the abscissas.
  !

    integer(ip), parameter :: order = 14

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(14) = (/ &
      0.13283874668559071814_dp, &
      0.13283874668559071814_dp, &
      0.13283874668559071814_dp, &
      0.13283874668559071814_dp, &
      0.088589824742980710434_dp, &
      0.088589824742980710434_dp, &
      0.088589824742980710434_dp, &
      0.088589824742980710434_dp, &
      0.019047619047619047619_dp, &
      0.019047619047619047619_dp, &
      0.019047619047619047619_dp, &
      0.019047619047619047619_dp, &
      0.019047619047619047619_dp, &
      0.019047619047619047619_dp  /)
    real(dp), intent(out) :: xyz(3,order)
    real(dp) :: xyz_save(3,14) = reshape ( (/ &
      0.056881379520423421748_dp, &
      0.31437287349319219275_dp, &
      0.31437287349319219275_dp, &
      0.31437287349319219275_dp, &
      0.056881379520423421748_dp, &
      0.31437287349319219275_dp, &
      0.31437287349319219275_dp, &
      0.31437287349319219275_dp, &
      0.056881379520423421748_dp, &
      0.31437287349319219275_dp, &
      0.31437287349319219275_dp, &
      0.31437287349319219275_dp, &
      0.69841970432438656092_dp, &
      0.10052676522520447969_dp, &
      0.10052676522520447969_dp, &
      0.10052676522520447969_dp, &
      0.69841970432438656092_dp, &
      0.10052676522520447969_dp, &
      0.10052676522520447969_dp, &
      0.10052676522520447969_dp, &
      0.69841970432438656092_dp, &
      0.10052676522520447969_dp, &
      0.10052676522520447969_dp, &
      0.10052676522520447969_dp, &
      0.50000000000000000000_dp, &
      0.50000000000000000000_dp, &
      0.00000000000000000000_dp, &
      0.50000000000000000000_dp, &
      0.00000000000000000000_dp, &
      0.50000000000000000000_dp, &
      0.50000000000000000000_dp, &
      0.00000000000000000000_dp, &
      0.00000000000000000000_dp, &
      0.00000000000000000000_dp, &
      0.50000000000000000000_dp, &
      0.50000000000000000000_dp, &
      0.00000000000000000000_dp, &
      0.50000000000000000000_dp, &
      0.00000000000000000000_dp, &
      0.00000000000000000000_dp, &
      0.00000000000000000000_dp, &
      0.50000000000000000000_dp /), &
    (/ 3, 14 /) )

    w(1:order) = w_save(1:order)
    xyz(1:3,1:order) = xyz_save(1:3,1:order)
  end subroutine tetr_unit_o14b

  pure subroutine tetr_unit_o15 ( w, xyz ) &
        bind(C, name="tetr_unit_o15")

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
  !    Output, real(dp) W(15), the weights.
  !
  !    Output, real(dp) XYZ(3,15), the abscissas.
  !

    integer(ip), parameter :: order = 15

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(15) = (/ &
      0.071937083779018620010_dp, &
      0.071937083779018620010_dp, &
      0.071937083779018620010_dp, &
      0.071937083779018620010_dp, &
      0.069068207226272385281_dp, &
      0.069068207226272385281_dp, &
      0.069068207226272385281_dp, &
      0.069068207226272385281_dp, &
      0.052910052910052910053_dp, &
      0.052910052910052910053_dp, &
      0.052910052910052910053_dp, &
      0.052910052910052910053_dp, &
      0.052910052910052910053_dp, &
      0.052910052910052910053_dp, &
      0.11851851851851851852_dp /)
    real(dp), intent(out) :: xyz(3,order)
    real(dp) :: xyz_save(3,15) = reshape ( (/ &
      0.72408676584183090163_dp, &
    0.091971078052723032789_dp, &
    0.091971078052723032789_dp, &
      0.091971078052723032789_dp, &
    0.72408676584183090163_dp, &
    0.091971078052723032789_dp, &
      0.091971078052723032789_dp, &
    0.091971078052723032789_dp, &
    0.72408676584183090163_dp, &
      0.091971078052723032789_dp, &
    0.091971078052723032789_dp, &
    0.091971078052723032789_dp, &
      0.040619116511110274837_dp, &
    0.31979362782962990839_dp, &
    0.31979362782962990839_dp, &
      0.31979362782962990839_dp, &
    0.040619116511110274837_dp, &
    0.31979362782962990839_dp, &
      0.31979362782962990839_dp, &
    0.31979362782962990839_dp, &
    0.040619116511110274837_dp, &
      0.31979362782962990839_dp, &
    0.31979362782962990839_dp, &
    0.31979362782962990839_dp, &
      0.44364916731037084426_dp, &
    0.44364916731037084426_dp, &
    0.056350832689629155741_dp, &
      0.44364916731037084426_dp, &
    0.056350832689629155741_dp, &
    0.44364916731037084426_dp, &
      0.44364916731037084426_dp, &
    0.056350832689629155741_dp, &
    0.056350832689629155741_dp, &
      0.056350832689629155741_dp, &
    0.44364916731037084426_dp, &
    0.44364916731037084426_dp, &
      0.056350832689629155741_dp, &
    0.44364916731037084426_dp, &
    0.056350832689629155741_dp, &
      0.056350832689629155741_dp, &
    0.056350832689629155741_dp, &
    0.44364916731037084426_dp, &
      0.25000000000000000000_dp, &
    0.25000000000000000000_dp, &
    0.25000000000000000000_dp /), &
    (/ 3, 15 /) )

    w(1:order) = w_save(1:order)
    xyz(1:3,1:order) = xyz_save(1:3,1:order)
  end subroutine tetr_unit_o15

  pure subroutine tetr_unit_o15b ( w, xyz ) &
        bind(C, name="tetr_unit_o15b")

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
  !    Output, real(dp) W(15), the weights.
  !
  !    Output, real(dp) XYZ(3,15), the abscissas.
  !

    integer(ip), parameter :: order = 15

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(15) = (/ &
      0.036160714285714285714_dp, &
      0.036160714285714285714_dp, &
      0.036160714285714285714_dp, &
      0.036160714285714285714_dp, &
      0.069871494516173816465_dp, &
      0.069871494516173816465_dp, &
      0.069871494516173816465_dp, &
      0.069871494516173816465_dp, &
      0.065694849368318756074_dp, &
      0.065694849368318756074_dp, &
      0.065694849368318756074_dp, &
      0.065694849368318756074_dp, &
      0.065694849368318756074_dp, &
      0.065694849368318756074_dp, &
      0.18170206858253505484_dp /)
    real(dp), intent(out) :: xyz(3,order)
    real(dp) :: xyz_save(3,15) = reshape ( (/ &
      0.00000000000000000000_dp, &
    0.33333333333333333333_dp, &
    0.33333333333333333333_dp, &
      0.33333333333333333333_dp, &
    0.00000000000000000000_dp, &
    0.33333333333333333333_dp, &
      0.33333333333333333333_dp, &
    0.33333333333333333333_dp, &
    0.00000000000000000000_dp, &
      0.33333333333333333333_dp, &
    0.33333333333333333333_dp, &
    0.33333333333333333333_dp, &
      0.72727272727272727273_dp, &
    0.090909090909090909091_dp, &
    0.090909090909090909091_dp, &
      0.090909090909090909091_dp, &
    0.72727272727272727273_dp, &
    0.090909090909090909091_dp, &
      0.090909090909090909091_dp, &
    0.090909090909090909091_dp, &
    0.72727272727272727273_dp, &
      0.090909090909090909091_dp, &
    0.090909090909090909091_dp, &
    0.090909090909090909091_dp, &
      0.43344984642633570176_dp, &
    0.43344984642633570176_dp, &
    0.066550153573664298240_dp, &
      0.43344984642633570176_dp, &
    0.066550153573664298240_dp, &
    0.43344984642633570176_dp, &
      0.43344984642633570176_dp, &
    0.066550153573664298240_dp, &
    0.066550153573664298240_dp, &
      0.066550153573664298240_dp, &
    0.43344984642633570176_dp, &
    0.43344984642633570176_dp, &
      0.066550153573664298240_dp, &
    0.43344984642633570176_dp, &
    0.066550153573664298240_dp, &
      0.066550153573664298240_dp, &
    0.066550153573664298240_dp, &
    0.43344984642633570176_dp, &
      0.25000000000000000000_dp, &
    0.25000000000000000000_dp, &
    0.250000000000000000_dp /), &
    (/ 3, 15 /) )

    w(1:order) = w_save(1:order)
    xyz(1:3,1:order) = xyz_save(1:3,1:order)
  end subroutine tetr_unit_o15b

  pure subroutine tetr_unit_o24 ( w, xyz ) &
        bind(C, name="tetr_unit_o24")

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
  !    Output, real(dp) W(24), the weights.
  !
  !    Output, real(dp) XYZ(3,24), the abscissas.
  !

    integer(ip), parameter :: order = 24

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(24) = (/ &
      0.039922750257869636194_dp, &
      0.039922750257869636194_dp, &
      0.039922750257869636194_dp, &
      0.039922750257869636194_dp, &
      0.010077211055345822612_dp, &
      0.010077211055345822612_dp, &
      0.010077211055345822612_dp, &
      0.010077211055345822612_dp, &
      0.055357181543927398338_dp, &
      0.055357181543927398338_dp, &
      0.055357181543927398338_dp, &
      0.055357181543927398338_dp, &
      0.048214285714285714286_dp, &
      0.048214285714285714286_dp, &
      0.048214285714285714286_dp, &
      0.048214285714285714286_dp, &
      0.048214285714285714286_dp, &
      0.048214285714285714286_dp, &
      0.048214285714285714286_dp, &
      0.048214285714285714286_dp, &
      0.048214285714285714286_dp, &
      0.048214285714285714286_dp, &
      0.048214285714285714286_dp, &
      0.048214285714285714286_dp /)
    real(dp), intent(out) :: xyz(3,order)
    real(dp) :: xyz_save(3,24) = reshape ( (/ &
      0.35619138622025439121_dp, &
      0.21460287125991520293_dp, &
      0.21460287125991520293_dp, &
      0.21460287125991520293_dp, &
      0.35619138622025439121_dp, &
      0.21460287125991520293_dp, &
      0.21460287125991520293_dp, &
      0.21460287125991520293_dp, &
      0.35619138622025439121_dp, &
      0.21460287125991520293_dp, &
      0.21460287125991520293_dp, &
      0.21460287125991520293_dp, &
      0.87797812439616594065_dp, &
      0.040673958534611353116_dp, &
      0.040673958534611353116_dp, &
      0.040673958534611353116_dp, &
      0.87797812439616594065_dp, &
      0.040673958534611353116_dp, &
      0.040673958534611353116_dp, &
      0.040673958534611353116_dp, &
      0.87797812439616594065_dp, &
      0.040673958534611353116_dp, &
      0.040673958534611353116_dp, &
      0.040673958534611353116_dp, &
      0.032986329573173468968_dp, &
      0.32233789014227551034_dp, &
      0.32233789014227551034_dp, &
      0.32233789014227551034_dp, &
      0.032986329573173468968_dp, &
      0.32233789014227551034_dp, &
      0.32233789014227551034_dp, &
      0.32233789014227551034_dp, &
      0.032986329573173468968_dp, &
      0.32233789014227551034_dp, &
      0.32233789014227551034_dp, &
      0.32233789014227551034_dp, &
      0.60300566479164914137_dp, &
      0.26967233145831580803_dp, &
      0.063661001875017525299_dp, &
      0.60300566479164914137_dp, &
      0.063661001875017525299_dp, &
      0.26967233145831580803_dp, &
      0.60300566479164914137_dp, &
      0.063661001875017525299_dp, &
      0.063661001875017525299_dp, &
      0.063661001875017525299_dp, &
      0.60300566479164914137_dp, &
      0.26967233145831580803_dp, &
      0.063661001875017525299_dp, &
      0.60300566479164914137_dp, &
      0.063661001875017525299_dp, &
      0.063661001875017525299_dp, &
      0.063661001875017525299_dp, &
      0.60300566479164914137_dp, &
      0.26967233145831580803_dp, &
      0.60300566479164914137_dp, &
      0.063661001875017525299_dp, &
      0.26967233145831580803_dp, &
      0.063661001875017525299_dp, &
      0.60300566479164914137_dp, &
      0.26967233145831580803_dp, &
      0.063661001875017525299_dp, &
      0.063661001875017525299_dp, &
      0.063661001875017525299_dp, &
      0.26967233145831580803_dp, &
      0.60300566479164914137_dp, &
      0.063661001875017525299_dp, &
      0.26967233145831580803_dp, &
      0.063661001875017525299_dp, &
      0.063661001875017525299_dp, &
      0.063661001875017525299_dp, &
      0.26967233145831580803_dp /), &
    (/ 3, 24 /) )

    w(1:order) = w_save(1:order)
    xyz(1:3,1:order) = xyz_save(1:3,1:order)
  end subroutine tetr_unit_o24

  subroutine tetr_unit_quad_test ( degree_max ) &
        bind(C, name="tetr_unit_quad_test")

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
  !    Input, integer(ip) DEGREE_MAX, the maximum total degree of the
  !    monomials to check.
  !

    integer(ip), parameter :: dim_num = 3

    integer(ip), intent(in), value :: degree_max
    integer(ip) :: expon(dim_num)
    integer(ip) :: h
    logical :: more
    integer(ip) :: order
    real(dp) :: quad
    integer(ip) :: t
    real(dp) :: tetr_unit_volume
    real(dp), allocatable :: v(:)
    real(dp), allocatable :: w(:)
    real(dp), allocatable :: xyz(:,:)

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
  end subroutine tetr_unit_quad_test

  pure function tetr_unit_volume ( ) &
        bind(C, name="tetr_unit_volume")

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
  !    Output, real(dp) TETR_UNIT_VOLUME, the volume.
  !

    real(dp) :: tetr_unit_volume

    tetr_unit_volume = 1.0_dp / 6.0_dp
  end function tetr_unit_volume

  pure subroutine trig_unit_monomial ( expon, value ) &
        bind(C, name="trig_unit_monomial")

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
  !    Input, integer(ip) EXPON(2), the exponents.
  !
  !    Output, real(dp) VALUE, the integral of the monomial.
  !

    integer(ip), intent(in) :: expon(2)
    integer(ip) :: i
    integer(ip) :: k
    real(dp), intent(out) :: value
  !
  !  The first computation ends with VALUE = 1.0;
  !
    value = 1.0_dp

  ! k = 0
  !
  !  The first loop simply computes 1 so we short circuit it!
  !
  ! do i = 1, expon(1)
  !   k = k + 1
  !   value = value * real ( i, dp) / real ( k, dp)
  ! end do

    k = expon(1)

    do i = 1, expon(2)
      k = k + 1
      value = value * real ( i, dp) / real ( k, dp)
    end do

    k = k + 1
    value = value / real ( k, dp)

    k = k + 1
    value = value / real ( k, dp)
  end subroutine trig_unit_monomial

  subroutine trig_unit_monomial_test ( degree_max ) &
        bind(C, name="trig_unit_monomial_test")

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
  !    Input, integer(ip) DEGREE_MAX, the maximum total degree of the
  !    monomials to check.
  !

    integer(ip) :: alpha
    integer(ip) :: beta
    integer(ip), intent(in), value :: degree_max
    integer(ip) :: expon(2)
    real(dp) :: trig_unit_volume
    real(dp) :: value

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
  end subroutine trig_unit_monomial_test

  pure subroutine trig_unit_o01 ( w, xy ) &
        bind(C, name="trig_unit_o01")

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
  !    Output, real(dp) W(1), the weights.
  !
  !    Output, real(dp) XY(2,1), the abscissas.
  !

    integer(ip), parameter :: order = 1

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(1) = (/ &
      1.0_dp /)
    real(dp), intent(out) :: xy(2,order)
    real(dp) :: xy_save(2,1) = reshape ( (/ &
      0.33333333333333333333_dp, &
      0.33333333333333333333_dp /), &
    (/ 2, 1 /) )

    w(1:order) = w_save(1:order)
    xy(1:2,1:order) = xy_save(1:2,1:order)
  end subroutine trig_unit_o01

  pure subroutine trig_unit_o03 ( w, xy ) &
        bind(C, name="trig_unit_o03")

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
  !    Output, real(dp) W(3), the weights.
  !
  !    Output, real(dp) XY(2,3), the abscissas.
  !

    integer(ip), parameter :: order = 3

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(3) = (/ &
      0.33333333333333333333_dp, &
      0.33333333333333333333_dp, &
      0.33333333333333333333_dp /)
    real(dp), intent(out) :: xy(2,order)
    real(dp) :: xy_save(2,3) = reshape ( (/ &
      0.66666666666666666667_dp, &
      0.16666666666666666667_dp, &
      0.16666666666666666667_dp, &
      0.66666666666666666667_dp, &
      0.16666666666666666667_dp, &
      0.16666666666666666667_dp /), &
    (/ 2, 3 /) )

    w(1:order) = w_save(1:order)
    xy(1:2,1:order) = xy_save(1:2,1:order)
  end subroutine trig_unit_o03

  pure subroutine trig_unit_o03b ( w, xy ) &
        bind(C, name="trig_unit_o03b")

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
  !    Output, real(dp) W(3), the weights.
  !
  !    Output, real(dp) XY(2,3), the abscissas.
  !

    integer(ip), parameter :: order = 3

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(3) = (/ &
      0.33333333333333333333_dp, &
      0.33333333333333333333_dp, &
      0.33333333333333333333_dp /)
    real(dp), intent(out) :: xy(2,order)
    real(dp) :: xy_save(2,3) = reshape ( (/ &
      0.0_dp, &
      0.5_dp, &
      0.5_dp, &
      0.0_dp, &
      0.5_dp, &
      0.5_dp /), &
    (/ 2, 3 /) )

    w(1:order) = w_save(1:order)
    xy(1:2,1:order) = xy_save(1:2,1:order)
  end subroutine trig_unit_o03b

  pure subroutine trig_unit_o06 ( w, xy ) &
        bind(C, name="trig_unit_o06")

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
  !    Output, real(dp) W(6), the weights.
  !
  !    Output, real(dp) XY(2,6), the abscissas.
  !

    integer(ip), parameter :: order = 6

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(6) = (/ &
      0.22338158967801146570_dp, &
      0.22338158967801146570_dp, &
      0.22338158967801146570_dp, &
      0.10995174365532186764_dp, &
      0.10995174365532186764_dp, &
      0.10995174365532186764_dp /)
    real(dp), intent(out) :: xy(2,order)
    real(dp) :: xy_save(2,6) = reshape ( (/ &
      0.10810301816807022736_dp, &
      0.44594849091596488632_dp, &
      0.44594849091596488632_dp, &
      0.10810301816807022736_dp, &
      0.44594849091596488632_dp, &
      0.44594849091596488632_dp, &
      0.81684757298045851308_dp, &
      0.091576213509770743460_dp, &
      0.091576213509770743460_dp, &
      0.81684757298045851308_dp, &
      0.091576213509770743460_dp, &
      0.091576213509770743460_dp /), &
    (/ 2, 6 /) )

    w(1:order) = w_save(1:order)
    xy(1:2,1:order) = xy_save(1:2,1:order)
  end subroutine trig_unit_o06

  pure subroutine trig_unit_o06b ( w, xy ) &
        bind(C, name="trig_unit_o06b")

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
  !    Output, real(dp) W(6), the weights.
  !
  !    Output, real(dp) XY(2,6), the abscissas.
  !

    integer(ip), parameter :: order = 6

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(6) = (/ &
      0.30000000000000000000_dp, &
      0.30000000000000000000_dp, &
      0.30000000000000000000_dp, &
      0.033333333333333333333_dp, &
      0.033333333333333333333_dp, &
      0.033333333333333333333_dp /)
    real(dp), intent(out) :: xy(2,order)
    real(dp) :: xy_save(2,6) = reshape ( (/ &
      0.66666666666666666667_dp, &
      0.16666666666666666667_dp, &
      0.16666666666666666667_dp, &
      0.66666666666666666667_dp, &
      0.16666666666666666667_dp, &
      0.16666666666666666667_dp, &
      0.0_dp, &
      0.5_dp, &
      0.5_dp, &
      0.0_dp, &
      0.5_dp, &
      0.5_dp /), &
    (/ 2, 6 /) )

    w(1:order) = w_save(1:order)
    xy(1:2,1:order) = xy_save(1:2,1:order)
  end subroutine trig_unit_o06b

  pure subroutine trig_unit_o07 ( w, xy ) &
        bind(C, name="trig_unit_o07")

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
  !    Output, real(dp) W(7), the weights.
  !
  !    Output, real(dp) XY(2,7), the abscissas.
  !

    integer(ip), parameter :: order = 7

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(7) = (/ &
      0.12593918054482715260_dp, &
      0.12593918054482715260_dp, &
      0.12593918054482715260_dp, &
      0.13239415278850618074_dp, &
      0.13239415278850618074_dp, &
      0.13239415278850618074_dp, &
      0.22500000000000000000_dp /)
    real(dp), intent(out) :: xy(2,order)
    real(dp) :: xy_save(2,7) = reshape ( (/ &
      0.79742698535308732240_dp, &
      0.10128650732345633880_dp, &
      0.10128650732345633880_dp, &
      0.79742698535308732240_dp, &
      0.10128650732345633880_dp, &
      0.10128650732345633880_dp, &
      0.059715871789769820459_dp, &
      0.47014206410511508977_dp, &
      0.47014206410511508977_dp, &
      0.059715871789769820459_dp, &
      0.47014206410511508977_dp, &
      0.47014206410511508977_dp, &
      0.33333333333333333333_dp, &
      0.33333333333333333333_dp /), &
    (/ 2, 7 /) )

    w(1:order) = w_save(1:order)
    xy(1:2,1:order) = xy_save(1:2,1:order)
  end subroutine trig_unit_o07

  pure subroutine trig_unit_o12 ( w, xy ) &
        bind(C, name="trig_unit_o12")

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
  !    Output, real(dp) W(12), the weights.
  !
  !    Output, real(dp) XY(2,12), the abscissas.
  !

    integer(ip), parameter :: order = 12

    real(dp), intent(out) :: w(order)
    real(dp) :: w_save(12) = (/ &
       0.050844906370206816921_dp, &
       0.050844906370206816921_dp, &
       0.050844906370206816921_dp, &
       0.11678627572637936603_dp, &
       0.11678627572637936603_dp, &
       0.11678627572637936603_dp, &
       0.082851075618373575194_dp, &
       0.082851075618373575194_dp, &
       0.082851075618373575194_dp, &
       0.082851075618373575194_dp, &
       0.082851075618373575194_dp, &
       0.082851075618373575194_dp /)
    real(dp), intent(out) :: xy(2,order)
    real(dp) :: xy_save(2,12) = reshape ( (/ &
      0.87382197101699554332_dp, &
      0.063089014491502228340_dp, &
      0.063089014491502228340_dp, &
      0.87382197101699554332_dp, &
      0.063089014491502228340_dp, &
      0.063089014491502228340_dp, &
      0.50142650965817915742_dp, &
      0.24928674517091042129_dp, &
      0.24928674517091042129_dp, &
      0.50142650965817915742_dp, &
      0.24928674517091042129_dp, &
      0.24928674517091042129_dp, &
      0.053145049844816947353_dp, &
      0.31035245103378440542_dp, &
      0.31035245103378440542_dp, &
      0.053145049844816947353_dp, &
      0.053145049844816947353_dp, &
      0.63650249912139864723_dp, &
      0.31035245103378440542_dp, &
      0.63650249912139864723_dp, &
      0.63650249912139864723_dp, &
      0.053145049844816947353_dp, &
      0.63650249912139864723_dp, &
      0.31035245103378440542_dp /), &
    (/ 2, 12 /) )

    w(1:order) = w_save(1:order)
    xy(1:2,1:order) = xy_save(1:2,1:order)
  end subroutine trig_unit_o12

  subroutine trig_unit_quad_test ( degree_max ) &
        bind(C, name="trig_unit_quad_test")

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
  !    Input, integer(ip) DEGREE_MAX, the maximum total degree of the
  !    monomials to check.
  !

    integer(ip), parameter :: dim_num = 2

    integer(ip), intent(in), value :: degree_max
    integer(ip) :: expon(dim_num)
    integer(ip) :: h
    logical :: more
    integer(ip) :: order
    real(dp) :: quad
    integer(ip) :: t
    real(dp) :: trig_unit_volume
    real(dp), allocatable :: v(:)
    real(dp), allocatable :: w(:)
    real(dp), allocatable :: xy(:,:)

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
  end subroutine trig_unit_quad_test

  pure function trig_unit_volume ( ) &
        bind(C, name="trig_unit_volume")

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
  !    Output, real(dp) TRIG_UNIT_VOLUME, the volume.
  !

    real(dp) :: trig_unit_volume

    trig_unit_volume = 0.5_dp
  end function trig_unit_volume

  subroutine wedg_unit_monomial ( expon, value ) &
        bind(C, name="wedg_unit_monomial")

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
  !    Input, integer(ip) EXPON(3), the exponents.
  !
  !    Output, real(dp) VALUE, the integral of the monomial.
  !

    integer(ip), intent(in) :: expon(3)
    integer(ip) :: i
    integer(ip) :: k
    real(dp), intent(out) :: value
  !
  !  The first computation ends with VALUE = 1.0;
  !
    value = 1.0_dp

  ! k = 0
  !
  !  The first loop simply computes 1 so we short circuit it!
  !
  ! do i = 1, expon(1)
  !   k = k + 1
  !   value = value * real ( i, dp) / real ( k, dp)
  ! end do

    k = expon(1)

    do i = 1, expon(2)
      k = k + 1
      value = value * real ( i, dp) / real ( k, dp)
    end do

    k = k + 1
    value = value / real ( k, dp)

    k = k + 1
    value = value / real ( k, dp)
  !
  !  Now account for integration in Z.
  !
    if ( expon(3) == - 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'WEDG_UNIT_MONOMIAL - Fatal error!'
      write ( *, '(a)' ) '  EXPON(3) = -1 is not a legal input.'
      stop 1
    else if ( mod ( expon(3), 2 ) == 1 ) then
      value = 0.0_dp
    else
      value = value * 2.0_dp / real ( expon(3) + 1, dp)
    end if
  end subroutine wedg_unit_monomial

  subroutine wedg_unit_monomial_test ( degree_max ) &
        bind(C, name="wedg_unit_monomial_test")

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
  !    Input, integer(ip) DEGREE_MAX, the maximum total degree of the
  !    monomials to check.
  !

    integer(ip) :: alpha
    integer(ip) :: beta
    integer(ip), intent(in), value :: degree_max
    integer(ip) :: expon(3)
    integer(ip) :: gamma
    real(dp) :: value
    real(dp) :: wedg_unit_volume

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
  end subroutine wedg_unit_monomial_test

  subroutine wedg_unit_quad_test ( degree_max ) &
        bind(C, name="wedg_unit_quad_test")

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
  !    Input, integer(ip) DEGREE_MAX, the maximum total degree of the
  !    monomials to check.
  !

    integer(ip), parameter :: dim_num = 3
    integer(ip), parameter :: test_num = 7

    integer(ip), intent(in), value :: degree_max
    integer(ip) :: expon(dim_num)
    integer(ip) :: h
    integer(ip) :: line_order
    integer(ip) :: line_order_array(test_num) = (/ &
      1, 2, 2, 3, 2, 3, 4 /)
    logical :: more
    integer(ip) :: order
    real(dp) :: quad
    integer(ip) :: t
    integer(ip) :: test
    integer(ip) :: trig_order
    integer(ip) :: trig_order_index
    integer(ip) :: trig_order_array(test_num) = (/ &
      1, 3, -3, 6, -6, 7, 12 /)
    real(dp) :: wedg_unit_volume
    real(dp), allocatable :: v(:)
    real(dp), allocatable :: w(:)
    real(dp), allocatable :: xyz(:,:)

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
  end subroutine wedg_unit_quad_test

  subroutine wedg_unit_rule ( line_order, trig_order, w, xyz ) &
        bind(C, name="wedg_unit_rule")

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
  !    Input, integer(ip) LINE_ORDER, the index of the line rule.
  !    The index of the rule is equal to the order of the rule.
  !    1 <= LINE_ORDER <= 5.
  !
  !    Input, integer(ip) TRIG_ORDER, the indes of the triangle rule.
  !    The index of the rule is 1, 3, -3, 6, -6, 7 or 12.
  !
  !    Output, real(dp) W(LINE_ORDER*abs(TRIG_ORDER)), the weights.
  !
  !    Output, real(dp) XYZ(3,LINE_ORDER*abs(TRIG_ORDER)), the abscissas.
  !

    integer(ip), intent(in), value :: line_order
    integer(ip), intent(in), value :: trig_order

    integer(ip) :: i
    integer(ip) :: j
    integer(ip) :: k
    real(dp) :: line_w(line_order)
    real(dp) :: line_x(line_order)
    real(dp) :: trig_w(abs(trig_order))
    real(dp) :: trig_xy(2,abs(trig_order))
    real(dp), intent(out) :: w(line_order*abs(trig_order))
    real(dp), intent(out) :: xyz(3,line_order*abs(trig_order))

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
  end subroutine wedg_unit_rule

  pure function wedg_unit_volume ( ) &
        bind(C, name="wedg_unit_volume")

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
  !    Output, real(dp) WEDG_UNIT_VOLUME, the volume.
  !

    real(dp) :: wedg_unit_volume

    wedg_unit_volume = 1.0_dp
  end function wedg_unit_volume

  subroutine wedg_unit_write_test ( ) &
        bind(C, name="wedg_unit_write_test")

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

    integer(ip), parameter :: dim_num = 3
    integer(ip), parameter :: rule_num = 7

    integer(ip) :: line_order
    integer(ip) :: line_order_array(rule_num) = (/ &
      1, 2, 2, 3, 2, 3, 4 /)
    integer(ip) :: order
    integer(ip) :: rule
    integer(ip) :: trig_order
    integer(ip) :: trig_order_array(rule_num) = (/ &
      1, 3, -3, 6, -6, 7, 12 /)
    real(dp), allocatable :: w(:)
    character ( len = 255 ) :: w_filename
    real(dp), allocatable :: x(:,:)
    character ( len = 255 ) :: x_filename

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
  end subroutine wedg_unit_write_test

end module felippa_mod
