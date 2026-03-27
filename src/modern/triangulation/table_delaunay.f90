!> table_delaunay — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module table_delaunay_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: diaedg, dtris2, i4_modp, i4_sign, i4_wrap, i4vec_indicator
  public :: lrline, perm_check, perm_inverse, r82vec_permute, r82vec_sort_heap_index_a, swapec
  public :: vbedg

contains

  function diaedg ( x0, y0, x1, y1, x2, y2, x3, y3 ) &
        bind(C, name="diaedg")

  !*****************************************************************************80
  !
  !! DIAEDG chooses a diagonal edge.
  !
  !  Discussion:
  !
  !    The routine determines whether 0--2 or 1--3 is the diagonal edge
  !    that should be chosen, based on the circumcircle criterion, where
  !    (X0,Y0), (X1,Y1), (X2,Y2), (X3,Y3) are the vertices of a simple
  !    quadrilateral in counterclockwise order.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 February 2001
  !
  !  Author:
  !
  !    Original FORTRAN77 version by Barry Joe.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Barry Joe,
  !    GEOMPACK - a software package for the generation of meshes
  !    using geometric algorithms,
  !    Advances in Engineering Software,
  !    Volume 13, pages 325-331, 1991.
  !
  !  Parameters:
  !
  !    Input, real(dp) X0, Y0, X1, Y1, X2, Y2, X3, Y3, the
  !    coordinates of the vertices of a quadrilateral, given in
  !    counter clockwise order.
  !
  !    Output, integer(ip) DIAEDG, chooses a diagonal:
  !    +1, if diagonal edge 02 is chosen;
  !    -1, if diagonal edge 13 is chosen;
  !     0, if the four vertices are cocircular.
  !

    real(dp) :: ca
    real(dp) :: cb
    integer(ip) :: diaedg
    real(dp) :: dx10
    real(dp) :: dx12
    real(dp) :: dx30
    real(dp) :: dx32
    real(dp) :: dy10
    real(dp) :: dy12
    real(dp) :: dy30
    real(dp) :: dy32
    real(dp) :: s
    real(dp) :: tol
    real(dp) :: tola
    real(dp) :: tolb
    real(dp), intent(in), value :: x0                                    !! x-coordinate of vertex 0
    real(dp), intent(in), value :: x1                                    !! x-coordinate of vertex 1
    real(dp), intent(in), value :: x2                                    !! x-coordinate of vertex 2
    real(dp), intent(in), value :: x3                                    !! x-coordinate of vertex 3
    real(dp), intent(in), value :: y0                                    !! y-coordinate of vertex 0
    real(dp), intent(in), value :: y1                                    !! y-coordinate of vertex 1
    real(dp), intent(in), value :: y2                                    !! y-coordinate of vertex 2
    real(dp), intent(in), value :: y3                                    !! y-coordinate of vertex 3

    tol = 100.0_dp * epsilon ( tol )

    dx10 = x1 - x0
    dy10 = y1 - y0
    dx12 = x1 - x2
    dy12 = y1 - y2
    dx30 = x3 - x0
    dy30 = y3 - y0
    dx32 = x3 - x2
    dy32 = y3 - y2

    tola = tol * max ( abs ( dx10 ), abs ( dy10 ), abs ( dx30 ), abs ( dy30 ) )
    tolb = tol * max ( abs ( dx12 ), abs ( dy12 ), abs ( dx32 ), abs ( dy32 ) )

    ca = dx10 * dx30 + dy10 * dy30
    cb = dx12 * dx32 + dy12 * dy32

    if ( tola < ca .and. tolb < cb ) then

      diaedg = -1

    else if ( ca < -tola .and. cb < -tolb ) then

      diaedg = 1

    else

      tola = max ( tola, tolb )
      s = ( dx10 * dy30 - dx30 * dy10 ) * cb + ( dx32 * dy12 - dx12 * dy32 ) * ca

      if ( tola < s ) then
        diaedg = -1
      else if ( s < -tola ) then
        diaedg = 1
      else
        diaedg = 0
      end if

    end if
  end function diaedg

  subroutine dtris2 ( point_num, point_xy, tri_num, tri_vert, tri_nabe ) &
        bind(C, name="dtris2")

  !*****************************************************************************80
  !
  !! DTRIS2 constructs a Delaunay triangulation of 2D vertices.
  !
  !  Discussion:
  !
  !    The routine constructs the Delaunay triangulation of a set of 2D vertices
  !    using an incremental approach and diagonal edge swaps.  Vertices are
  !    first sorted in lexicographically increasing (X,Y) order, and
  !    then are inserted one at a time from outside the convex hull.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    25 August 2001
  !
  !  Author:
  !
  !    Original FORTRAN77 version by Barry Joe.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Barry Joe,
  !    GEOMPACK - a software package for the generation of meshes
  !    using geometric algorithms,
  !    Advances in Engineering Software,
  !    Volume 13, pages 325-331, 1991.
  !
  !  Parameters:
  !
  !    Input, integer(ip) POINT_NUM, the number of vertices.
  !
  !    Input/output, real(dp) POINT_XY(2,POINT_NUM), the coordinates
  !    of the vertices.  On output, the vertices have been sorted into
  !    dictionary order.
  !
  !    Output, integer(ip) TRI_NUM, the number of triangles in the
  !    triangulation; TRI_NUM is equal to 2*POINT_NUM - NB - 2, where NB is the
  !    number of boundary vertices.
  !
  !    Output, integer(ip) TRI_VERT(3,TRI_NUM), the nodes that make up
  !    each triangle.  The elements are indices of POINT_XY.  The vertices of the
  !    triangles are in counter clockwise order.
  !
  !    Output, integer(ip) TRI_NABE(3,TRI_NUM), the triangle neighbor
  !    list.  Positive elements are indices of TIL; negative elements are used
  !    for links of a counter clockwise linked list of boundary edges;
  !    LINK = -(3*I + J-1) where I, J = triangle, edge index; TRI_NABE(J,I) refers
  !    to the neighbor along edge from vertex J to J+1 (mod 3).
  !

    integer(ip), intent(in), value :: point_num                          !! number of vertices

    real(dp) :: cmax
    integer(ip) :: e
    integer(ip) :: i
    integer(ip) :: ierr
    integer(ip) :: indx(point_num)
    integer(ip) :: j
    integer(ip) :: k
    integer(ip) :: l
    integer(ip) :: ledg
    integer(ip) :: lr
    integer(ip) :: lrline
    integer(ip) :: ltri
    integer(ip) :: m
    integer(ip) :: m1
    integer(ip) :: m2
    integer(ip) :: n
    real(dp), intent(inout) :: point_xy(2,point_num)
    integer(ip) :: redg
    integer(ip) :: rtri
    integer(ip) :: stack(point_num)
    integer(ip) :: t
    real(dp) :: tol
    integer(ip) :: top
    integer(ip), intent(out) :: tri_nabe(3,point_num*2)
    integer(ip), intent(out) :: tri_num
    integer(ip), intent(out) :: tri_vert(3,point_num*2)

    tol = 100.0_dp * epsilon ( tol )

    ierr = 0
  !
  !  Sort the vertices by increasing (x,y).
  !
    call r82vec_sort_heap_index_a ( point_num, point_xy, indx )

    call r82vec_permute ( point_num, indx, point_xy )
  !
  !  Make sure that the data points are "reasonably" distinct.
  !
    m1 = 1

    do i = 2, point_num

      m = m1
      m1 = i

      k = 0

      do j = 1, 2

        cmax = max ( abs ( point_xy(j,m) ), abs ( point_xy(j,m1) ) )

        if ( tol * ( cmax + 1.0_dp ) &
             < abs ( point_xy(j,m) - point_xy(j,m1) ) ) then
          k = j
          exit
        end if

      end do

      if ( k == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
        write ( *, '(a,i8)' ) '  Fails for point number I = ', i
        write ( *, '(a,i8)' ) '  M = ', m
        write ( *, '(a,i8)' ) '  M1 = ', m1
        write ( *, '(a,2g14.6)' ) '  X,Y(M)  = ', point_xy(1,m), point_xy(2,m)
        write ( *, '(a,2g14.6)' ) '  X,Y(M1) = ', point_xy(1,m1), point_xy(2,m1)
        ierr = 224
      end if

    end do
  !
  !  Starting from points M1 and M2, search for a third point M that
  !  makes a "healthy" triangle (M1,M2,M)
  !
    m1 = 1
    m2 = 2
    j = 3

    do

      if ( point_num < j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
        ierr = 225
      end if

      m = j

      lr = lrline ( point_xy(1,m), point_xy(2,m), point_xy(1,m1), &
        point_xy(2,m1), point_xy(1,m2), point_xy(2,m2), 0.0_dp )

      if ( lr /= 0 ) then
        exit
      end if

      j = j + 1

    end do
  !
  !  Set up the triangle information for (M1,M2,M), and for any other
  !  triangles you created because points were collinear with M1, M2.
  !
    tri_num = j - 2

    if ( lr == -1 ) then

      tri_vert(1,1) = m1
      tri_vert(2,1) = m2
      tri_vert(3,1) = m
      tri_nabe(3,1) = -3

      do i = 2, tri_num

        m1 = m2
        m2 = i+1
        tri_vert(1,i) = m1
        tri_vert(2,i) = m2
        tri_vert(3,i) = m
        tri_nabe(1,i-1) = -3 * i
        tri_nabe(2,i-1) = i
        tri_nabe(3,i) = i - 1

      end do

      tri_nabe(1,tri_num) = -3 * tri_num - 1
      tri_nabe(2,tri_num) = -5
      ledg = 2
      ltri = tri_num

    else

      tri_vert(1,1) = m2
      tri_vert(2,1) = m1
      tri_vert(3,1) = m
      tri_nabe(1,1) = -4

      do i = 2, tri_num
        m1 = m2
        m2 = i+1
        tri_vert(1,i) = m2
        tri_vert(2,i) = m1
        tri_vert(3,i) = m
        tri_nabe(3,i-1) = i
        tri_nabe(1,i) = -3 * i - 3
        tri_nabe(2,i) = i - 1
      end do

      tri_nabe(3,tri_num) = -3 * tri_num
      tri_nabe(2,1) = -3 * tri_num - 2
      ledg = 2
      ltri = 1

    end if
  !
  !  Insert the vertices one at a time from outside the convex hull,
  !  determine visible boundary edges, and apply diagonal edge swaps until
  !  Delaunay triangulation of vertices (so far) is obtained.
  !
    top = 0

    do i = j+1, point_num

      m = i
      m1 = tri_vert(ledg,ltri)

      if ( ledg <= 2 ) then
        m2 = tri_vert(ledg+1,ltri)
      else
        m2 = tri_vert(1,ltri)
      end if

      lr = lrline ( point_xy(1,m), point_xy(2,m), point_xy(1,m1), &
        point_xy(2,m1), point_xy(1,m2), point_xy(2,m2), 0.0_dp )

      if ( 0 < lr ) then
        rtri = ltri
        redg = ledg
        ltri = 0
      else
        l = -tri_nabe(ledg,ltri)
        rtri = l / 3
        redg = mod(l,3) + 1
      end if

      call vbedg ( point_xy(1,m), point_xy(2,m), point_num, point_xy, tri_num, &
        tri_vert, tri_nabe, ltri, ledg, rtri, redg )

      n = tri_num + 1
      l = -tri_nabe(ledg,ltri)

      do

        t = l / 3
        e = mod ( l, 3 ) + 1
        l = -tri_nabe(e,t)
        m2 = tri_vert(e,t)

        if ( e <= 2 ) then
          m1 = tri_vert(e+1,t)
        else
          m1 = tri_vert(1,t)
        end if

        tri_num = tri_num + 1
        tri_nabe(e,t) = tri_num
        tri_vert(1,tri_num) = m1
        tri_vert(2,tri_num) = m2
        tri_vert(3,tri_num) = m
        tri_nabe(1,tri_num) = t
        tri_nabe(2,tri_num) = tri_num - 1
        tri_nabe(3,tri_num) = tri_num + 1
        top = top + 1

        if ( point_num < top ) then
          ierr = 8
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
          write ( *, '(a)' ) '  Stack overflow.'
        end if

        stack(top) = tri_num

        if ( t == rtri .and. e == redg ) then
          exit
        end if

      end do

      tri_nabe(ledg,ltri) = -3 * n - 1
      tri_nabe(2,n) = -3 * tri_num - 2
      tri_nabe(3,tri_num) = -l
      ltri = n
      ledg = 2

      call swapec ( m, top, ltri, ledg, point_num, point_xy, tri_num, &
        tri_vert, tri_nabe, stack, ierr )

      if ( ierr /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
        write ( *, '(a)' ) '  Error return from SWAPEC.'
      end if

    end do
  !
  !  Now account for the sorting that we did.
  !
    do i = 1, 3
      do j = 1, tri_num
        tri_vert(i,j) = indx ( tri_vert(i,j) )
      end do
    end do

    call perm_inverse ( point_num, indx )

    call r82vec_permute ( point_num, indx, point_xy )
  end subroutine dtris2

  function i4_modp ( i, j ) &
        bind(C, name="i4_modp")

  !*****************************************************************************80
  !
  !! I4_MODP returns the nonnegative remainder of I4 division.
  !
  !  Discussion:
  !
  !    If
  !      NREM = I4_MODP ( I, J )
  !      NMULT = ( I - NREM ) / J
  !    then
  !      I = J * NMULT + NREM
  !    where NREM is always nonnegative.
  !
  !    The MOD function computes a result with the same sign as the
  !    quantity being divided.  Thus, suppose you had an angle A,
  !    and you wanted to ensure that it was between 0 and 360.
  !    Then mod(A,360) would do, if A was positive, but if A
  !    was negative, your result would be between -360 and 0.
  !
  !    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
  !
  !    An I4 is an integer(ip) value.
  !
  !  Example:
  !
  !        I     J     MOD I4_MODP    Factorization
  !
  !      107    50       7       7    107 =  2 *  50 + 7
  !      107   -50       7       7    107 = -2 * -50 + 7
  !     -107    50      -7      43   -107 = -3 *  50 + 43
  !     -107   -50      -7      43   -107 =  3 * -50 + 43
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    02 March 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) I, the number to be divided.
  !
  !    Input, integer(ip) J, the number that divides I.
  !
  !    Output, integer(ip) I4_MODP, the nonnegative remainder when I is
  !    divided by J.
  !

    integer(ip), intent(in), value :: i                                  !! number to be divided
    integer(ip) :: i4_modp
    integer(ip), intent(in), value :: j                                  !! divisor
    integer(ip) :: value

    if ( j == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'I4_MODP - Fatal error!'
      write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
      stop
    end if

    value = mod ( i, j )

    if ( value < 0 ) then
      value = value + abs ( j )
    end if

    i4_modp = value
  end function i4_modp

  function i4_sign ( x ) &
        bind(C, name="i4_sign")

  !*****************************************************************************80
  !
  !! I4_SIGN evaluates the sign of an I4.
  !
  !  Discussion:
  !
  !    An I4 is an integer(ip) value.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 March 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) X, the number whose sign is desired.
  !
  !    Output, integer(ip) I4_SIGN, the sign of X:
  !

    integer(ip) :: i4_sign
    integer(ip), intent(in), value :: x                                  !! number whose sign is desired

    if ( x < 0 ) then
      i4_sign = -1
    else
      i4_sign = +1
    end if
  end function i4_sign

  function i4_wrap ( ival, ilo, ihi ) &
        bind(C, name="i4_wrap")

  !*****************************************************************************80
  !
  !! I4_WRAP forces an I4 to lie between given limits by wrapping.
  !
  !  Discussion:
  !
  !    An I4 is an integer(ip) value.
  !
  !  Example:
  !
  !    ILO = 4, IHI = 8
  !
  !    I  Value
  !
  !    -2     8
  !    -1     4
  !     0     5
  !     1     6
  !     2     7
  !     3     8
  !     4     4
  !     5     5
  !     6     6
  !     7     7
  !     8     8
  !     9     4
  !    10     5
  !    11     6
  !    12     7
  !    13     8
  !    14     4
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 August 2003
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) IVAL, a value.
  !
  !    Input, integer(ip) ILO, IHI, the desired bounds.
  !
  !    Output, integer(ip) I4_WRAP, a "wrapped" version of the value.
  !

    integer(ip) :: i4_modp
    integer(ip) :: i4_wrap
    integer(ip), intent(in), value :: ihi                                !! upper bound
    integer(ip), intent(in), value :: ilo                                !! lower bound
    integer(ip), intent(in), value :: ival                               !! value to wrap
    integer(ip) :: jhi
    integer(ip) :: jlo
    integer(ip) :: value
    integer(ip) :: wide

    jlo = min ( ilo, ihi )
    jhi = max ( ilo, ihi )

    wide = jhi - jlo + 1

    if ( wide == 1 ) then
      value = jlo
    else
      value = jlo + i4_modp ( ival - jlo, wide )
    end if

    i4_wrap = value
  end function i4_wrap

  subroutine i4vec_indicator ( n, a ) &
        bind(C, name="i4vec_indicator")

  !*****************************************************************************80
  !
  !! I4VEC_INDICATOR sets an I4VEC to the indicator vector.
  !
  !  Discussion:
  !
  !    An I4VEC is a vector of I4's.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    01 May 2007
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of elements of A.
  !
  !    Output, integer(ip) A(N), the array to be initialized.
  !

    integer(ip), intent(in), value :: n                                  !! number of elements

    integer(ip), intent(out) :: a(n)
    integer(ip) :: i

    do i = 1, n
      a(i) = i
    end do
  end subroutine i4vec_indicator

  function lrline ( xu, yu, xv1, yv1, xv2, yv2, dv ) &
        bind(C, name="lrline")

  !*****************************************************************************80
  !
  !! LRLINE determines if a point is left of, right or, or on a directed line.
  !
  !  Discussion:
  !
  !    The directed line is parallel to, and at a signed distance DV from
  !    a directed base line from (XV1,YV1) to (XV2,YV2).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 July 2001
  !
  !  Author:
  !
  !    Original FORTRAN77 version by Barry Joe.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Barry Joe,
  !    GEOMPACK - a software package for the generation of meshes
  !    using geometric algorithms,
  !    Advances in Engineering Software,
  !    Volume 13, pages 325-331, 1991.
  !
  !  Parameters:
  !
  !    Input, real(dp) XU, YU, the coordinates of the point whose
  !    position relative to the directed line is to be determined.
  !
  !    Input, real(dp) XV1, YV1, XV2, YV2, the coordinates of two points
  !    that determine the directed base line.
  !
  !    Input, real(dp) DV, the signed distance of the directed line
  !    from the directed base line through the points (XV1,YV1) and (XV2,YV2).
  !    DV is positive for a line to the left of the base line.
  !
  !    Output, integer(ip) LRLINE, the result:
  !    +1, the point is to the right of the directed line;
  !     0, the point is on the directed line;
  !    -1, the point is to the left of the directed line.
  !

    real(dp), intent(in), value :: dv                                    !! signed distance
    real(dp) :: dx
    real(dp) :: dxu
    real(dp) :: dy
    real(dp) :: dyu
    integer(ip) :: lrline
    real(dp) :: t
    real(dp) :: tol
    real(dp) :: tolabs
    real(dp), intent(in), value :: xu                                    !! x-coordinate of test point
    real(dp), intent(in), value :: xv1                                   !! x-coordinate of base point 1
    real(dp), intent(in), value :: xv2                                   !! x-coordinate of base point 2
    real(dp), intent(in), value :: yu                                    !! y-coordinate of test point
    real(dp), intent(in), value :: yv1                                   !! y-coordinate of base point 1
    real(dp), intent(in), value :: yv2                                   !! y-coordinate of base point 2

    tol = 100.0_dp * epsilon ( tol )

    dx = xv2 - xv1
    dy = yv2 - yv1
    dxu = xu - xv1
    dyu = yu - yv1

    tolabs = tol * max ( abs ( dx ), abs ( dy ), abs ( dxu ), &
      abs ( dyu ), abs ( dv ) )

    t = dy * dxu - dx * dyu + dv * sqrt ( dx * dx + dy * dy )

    if ( tolabs < t ) then
      lrline = 1
    else if ( -tolabs <= t ) then
      lrline = 0
    else
      lrline = -1
    end if
  end function lrline

  subroutine perm_check ( n, p, base, ierror ) &
        bind(C, name="perm_check")

  !*****************************************************************************80
  !
  !! PERM_CHECK checks that a vector represents a permutation.
  !
  !  Discussion:
  !
  !    The routine verifies that each of the integers from BASE to
  !    to BASE+N-1 occurs among the N entries of the permutation.
  !
  !    Set the input quantity BASE to 0, if P is a 0-based permutation,
  !    or to 1 if P is a 1-based permutation.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    31 October 2008
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of entries.
  !
  !    Input, integer(ip) P(N), the array to check.
  !
  !    Input, integer(ip) BASE, the index base.
  !
  !    Output, integer(ip) IERROR, error flag.
  !    0, the array represents a permutation.
  !    nonzero, the array does not represent a permutation.  The smallest
  !    missing value is equal to IERROR.
  !

    integer(ip), intent(in), value :: n                                  !! number of entries

    integer(ip), intent(in), value :: base                               !! index base
    integer(ip) :: find
    integer(ip), intent(out) :: ierror
    integer(ip), intent(in) :: p(n)
    integer(ip) :: seek

    ierror = 0

    do seek = base, base + n - 1

      ierror = 1

      do find = 1, n
        if ( p(find) == seek ) then
          ierror = 0
          exit
        end if
      end do

      if ( ierror /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PERM_CHECK - Fatal error!'
        write ( *, '(a)' ) '  The input array does not represent'
        write ( *, '(a)' ) '  a proper permutation.'
        stop
      end if

    end do
  end subroutine perm_check

  subroutine perm_inverse ( n, p ) &
        bind(C, name="perm_inverse")

  !*****************************************************************************80
  !
  !! PERM_INVERSE inverts a permutation "in place".
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    02 January 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of objects being permuted.
  !
  !    Input/output, integer(ip) P(N), the permutation, in standard
  !    index form.  On output, P describes the inverse permutation
  !

    integer(ip), intent(in), value :: n                                  !! number of objects

    integer(ip), parameter :: base = 1
    integer(ip) :: i
    integer(ip) :: i0
    integer(ip) :: i1
    integer(ip) :: i2
    integer(ip) :: i4_sign
    integer(ip) :: ierror
    integer(ip) :: is
    integer(ip), intent(inout) :: p(n)

    if ( n <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PERM_INVERSE - Fatal error!'
      write ( *, '(a,i8)' ) '  Input value of N = ', n
      stop
    end if

    call perm_check ( n, p, base, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PERM_INVERSE - Fatal error!'
      write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
      stop
    end if


    is = 1

    do i = 1, n

      i1 = p(i)

      do while ( i < i1 )
        i2 = p(i1)
        p(i1) = -i2
        i1 = i2
      end do

      is = - i4_sign ( p(i) )
      p(i) = is * abs ( p(i) )

    end do

    do i = 1, n

      i1 = - p(i)

      if ( 0 <= i1 ) then

        i0 = i

        do

          i2 = p(i1)
          p(i1) = i0

          if ( i2 < 0 ) then
            exit
          end if

          i0 = i1
          i1 = i2

        end do

      end if

    end do
  end subroutine perm_inverse

  subroutine r82vec_permute ( n, p, a ) &
        bind(C, name="r82vec_permute")

  !*****************************************************************************80
  !
  !! R82VEC_PERMUTE permutes an R82VEC in place.
  !
  !  Discussion:
  !
  !    An R82VEC is an array of pairs of R8 values.
  !
  !    The same logic can be used to permute an array of objects of any
  !    arithmetic type, or an array of objects of any complexity.  The only
  !    temporary storage required is enough to store a single object.  The number
  !    of data movements made is N + the number of cycles of order 2 or more,
  !    which is never more than N + N/2.
  !
  !  Example:
  !
  !    Input:
  !
  !      N = 5
  !      P = (   2,    4,    5,    1,    3 )
  !      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
  !          (11.0, 22.0, 33.0, 44.0, 55.0 )
  !
  !    Output:
  !
  !      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
  !             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    13 March 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of objects.
  !
  !    Input, integer(ip) P(N), the permutation.  P(I) = J means
  !    that the I-th element of the output array should be the J-th
  !    element of the input array.
  !
  !    Input/output, real(dp) A(2,N), the array to be permuted.
  !

    integer(ip), intent(in), value :: n                                  !! number of objects
    integer(ip), parameter :: dim_num = 2

    real(dp), intent(inout) :: a(dim_num,n)
    real(dp) :: a_temp(dim_num)
    integer(ip), parameter :: base = 1
    integer(ip) :: ierror
    integer(ip) :: iget
    integer(ip) :: iput
    integer(ip) :: istart
    integer(ip), intent(inout) :: p(n)

    call perm_check ( n, p, base, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R82VEC_PERMUTE - Fatal error!'
      write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
      stop
    end if
  !
  !  Search for the next element of the permutation that has not been used.
  !
    do istart = 1, n

      if ( p(istart) < 0 ) then

        cycle

      else if ( p(istart) == istart ) then

        p(istart) = - p(istart)
        cycle

      else

        a_temp(1:dim_num) = a(1:dim_num,istart)
        iget = istart
  !
  !  Copy the new value into the vacated entry.
  !
        do

          iput = iget
          iget = p(iget)

          p(iput) = - p(iput)

          if ( iget < 1 .or. n < iget ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R82VEC_PERMUTE - Fatal error!'
            write ( *, '(a)' ) '  A permutation index is out of range.'
            write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
            stop
          end if

          if ( iget == istart ) then
            a(1:dim_num,iput) = a_temp(1:dim_num)
            exit
          end if

          a(1:dim_num,iput) = a(1:dim_num,iget)

        end do

      end if

    end do
  !
  !  Restore the signs of the entries.
  !
    p(1:n) = - p(1:n)
  end subroutine r82vec_permute

  subroutine r82vec_sort_heap_index_a ( n, a, indx ) &
        bind(C, name="r82vec_sort_heap_index_a")

  !*****************************************************************************80
  !
  !! R82VEC_SORT_HEAP_INDEX_A ascending index heaps an R82VEC.
  !
  !  Discussion:
  !
  !    An R82VEC is an array of R82's.
  !
  !    The sorting is not actually carried out.  Rather an index array is
  !    created which defines the sorting.  This array may be used to sort
  !    or index the array, or to sort or index related arrays keyed on the
  !    original array.
  !
  !    Once the index array is computed, the sorting can be carried out
  !    "implicitly:
  !
  !      A(1:2,INDX(1:N)) is sorted,
  !
  !    or explicitly, by the call
  !
  !      call r82vec_permute ( n, indx, a )
  !
  !    after which A(1:2,I), I = 1 to N is sorted.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    08 December 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of entries in the array.
  !
  !    Input, real(dp) A(2,N), an array to be index-sorted.
  !
  !    Output, integer(ip) INDX(N), the sort index.  The
  !    I-th element of the sorted array is A(1:2,INDX(I)).
  !

    integer(ip), intent(in), value :: n                                  !! number of entries
    integer(ip), parameter :: dim_num = 2

    real(dp), intent(in) :: a(dim_num,n)
    real(dp) :: aval(dim_num)
    integer(ip) :: i
    integer(ip), intent(out) :: indx(n)
    integer(ip) :: indxt
    integer(ip) :: ir
    integer(ip) :: j
    integer(ip) :: l

    if ( n < 1 ) then
    end if

    do i = 1, n
      indx(i) = i
    end do

    if ( n == 1 ) then
    end if

    l = n / 2 + 1
    ir = n

    do

      if ( 1 < l ) then

        l = l - 1
        indxt = indx(l)
        aval(1:dim_num) = a(1:dim_num,indxt)

      else

        indxt = indx(ir)
        aval(1:dim_num) = a(1:dim_num,indxt)
        indx(ir) = indx(1)
        ir = ir - 1

        if ( ir == 1 ) then
          indx(1) = indxt
          exit
        end if

      end if

      i = l
      j = l + l

      do while ( j <= ir )

        if ( j < ir ) then
          if (   a(1,indx(j)) <  a(1,indx(j+1)) .or. &
               ( a(1,indx(j)) == a(1,indx(j+1)) .and. &
                 a(2,indx(j)) <  a(2,indx(j+1)) ) ) then
            j = j + 1
          end if
        end if

        if (   aval(1) <  a(1,indx(j)) .or. &
             ( aval(1) == a(1,indx(j)) .and. &
               aval(2) <  a(2,indx(j)) ) ) then
          indx(i) = indx(j)
          i = j
          j = j + j
        else
          j = ir + 1
        end if

      end do

      indx(i) = indxt

    end do
  end subroutine r82vec_sort_heap_index_a

  subroutine swapec ( i, top, btri, bedg, point_num, point_xy, tri_num, &
    tri_vert, tri_nabe, stack, ierr ) &
        bind(C, name="swapec")

  !*****************************************************************************80
  !
  !! SWAPEC swaps diagonal edges until all triangles are Delaunay.
  !
  !  Discussion:
  !
  !    The routine swaps diagonal edges in a 2D triangulation, based on
  !    the empty circumcircle criterion, until all triangles are Delaunay,
  !    given that I is the index of the new vertex added to the triangulation.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 July 2001
  !
  !  Author:
  !
  !    Original FORTRAN77 version by Barry Joe.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Barry Joe,
  !    GEOMPACK - a software package for the generation of meshes
  !    using geometric algorithms,
  !    Advances in Engineering Software,
  !    Volume 13, pages 325-331, 1991.
  !
  !  Parameters:
  !
  !    Input, integer(ip) I, the index of the new vertex.
  !
  !    Input/output, integer(ip) TOP, the index of the top of the stack.
  !    On output, TOP is zero.
  !
  !    Input/output, integer(ip) BTRI, BEDG; on input, if positive, are
  !    the triangle and edge indices of a boundary edge whose updated indices
  !    must be recorded.  On output, these may be updated because of swaps.
  !
  !    Input, integer(ip) POINT_NUM, the number of points.
  !
  !    Input, real(dp) POINT_XY(2,POINT_NUM), the coordinates
  !    of the points.
  !
  !    Input, integer(ip) TRI_NUM, the number of triangles.
  !
  !    Input/output, integer(ip) TRI_VERT(3,TRI_NUM), the triangle
  !    incidence list.  May be updated on output because of swaps.
  !
  !    Input/output, integer(ip) TRI_NABE(3,TRI_NUM), the triangle
  !    neighbor list; negative values are used for links of the counter-clockwise
  !    linked list of boundary edges;  May be updated on output because of swaps.
  !      LINK = -(3*I + J-1) where I, J = triangle, edge index.
  !
  !    Workspace, integer(ip) STACK(MAXST); on input, entries 1 through
  !    TOP contain the indices of initial triangles (involving vertex I)
  !    put in stack; the edges opposite I should be in interior;  entries
  !    TOP+1 through MAXST are used as a stack.
  !
  !    Output, integer(ip) IERR is set to 8 for abnormal return.
  !

    integer(ip), intent(in), value :: point_num                          !! number of points
    integer(ip), intent(in), value :: tri_num                            !! number of triangles

    integer(ip) :: a
    integer(ip) :: b
    integer(ip), intent(inout) :: bedg
    integer(ip), intent(inout) :: btri
    integer(ip) :: c
    integer(ip) :: diaedg
    integer(ip) :: e
    integer(ip) :: ee
    integer(ip) :: em1
    integer(ip) :: ep1
    integer(ip) :: f
    integer(ip) :: fm1
    integer(ip) :: fp1
    integer(ip), intent(in), value :: i                                  !! index of new vertex
    integer(ip), intent(out) :: ierr
    integer(ip) :: i4_wrap
    integer(ip) :: l
    integer(ip) :: r
    integer(ip) :: s
    integer(ip), intent(inout) :: stack(point_num)
    integer(ip) :: swap
    integer(ip) :: t
    integer(ip), intent(inout) :: top
    integer(ip), intent(inout) :: tri_nabe(3,tri_num)
    integer(ip), intent(inout) :: tri_vert(3,tri_num)
    integer(ip) :: tt
    integer(ip) :: u
    real(dp), intent(in) :: point_xy(2,point_num)
    real(dp) :: x
    real(dp) :: y
  !
  !  Determine whether triangles in stack are Delaunay, and swap
  !  diagonal edge of convex quadrilateral if not.
  !
    x = point_xy(1,i)
    y = point_xy(2,i)

    do

      if ( top <= 0 ) then
        exit
      end if

      t = stack(top)
      top = top - 1

      if ( tri_vert(1,t) == i ) then
        e = 2
        b = tri_vert(3,t)
      else if ( tri_vert(2,t) == i ) then
        e = 3
        b = tri_vert(1,t)
      else
        e = 1
        b = tri_vert(2,t)
      end if

      a = tri_vert(e,t)
      u = tri_nabe(e,t)

      if ( tri_nabe(1,u) == t ) then
        f = 1
        c = tri_vert(3,u)
      else if ( tri_nabe(2,u) == t ) then
        f = 2
        c = tri_vert(1,u)
      else
        f = 3
        c = tri_vert(2,u)
      end if

      swap = diaedg ( x, y, point_xy(1,a), point_xy(2,a), point_xy(1,c), &
        point_xy(2,c), point_xy(1,b), point_xy(2,b) )

      if ( swap == 1 ) then

        em1 = i4_wrap ( e - 1, 1, 3 )
        ep1 = i4_wrap ( e + 1, 1, 3 )
        fm1 = i4_wrap ( f - 1, 1, 3 )
        fp1 = i4_wrap ( f + 1, 1, 3 )

        tri_vert(ep1,t) = c
        tri_vert(fp1,u) = i
        r = tri_nabe(ep1,t)
        s = tri_nabe(fp1,u)
        tri_nabe(ep1,t) = u
        tri_nabe(fp1,u) = t
        tri_nabe(e,t) = s
        tri_nabe(f,u) = r

        if ( 0 < tri_nabe(fm1,u) ) then
          top = top + 1
          stack(top) = u
        end if

        if ( 0 < s ) then

          if ( tri_nabe(1,s) == u ) then
            tri_nabe(1,s) = t
          else if ( tri_nabe(2,s) == u ) then
            tri_nabe(2,s) = t
          else
            tri_nabe(3,s) = t
          end if

          top = top + 1

          if ( point_num < top ) then
            ierr = 8
          end if

          stack(top) = t

        else

          if ( u == btri .and. fp1 == bedg ) then
            btri = t
            bedg = e
          end if

          l = - ( 3 * t + e - 1 )
          tt = t
          ee = em1

          do while ( 0 < tri_nabe(ee,tt) )

            tt = tri_nabe(ee,tt)

            if ( tri_vert(1,tt) == a ) then
              ee = 3
            else if ( tri_vert(2,tt) == a ) then
              ee = 1
            else
              ee = 2
            end if

          end do

          tri_nabe(ee,tt) = l

        end if

        if ( 0 < r ) then

          if ( tri_nabe(1,r) == t ) then
            tri_nabe(1,r) = u
          else if ( tri_nabe(2,r) == t ) then
            tri_nabe(2,r) = u
          else
            tri_nabe(3,r) = u
          end if

        else

          if ( t == btri .and. ep1 == bedg ) then
            btri = u
            bedg = f
          end if

          l = - ( 3 * u + f - 1 )
          tt = u
          ee = fm1

          do while ( 0 < tri_nabe(ee,tt) )

            tt = tri_nabe(ee,tt)

            if ( tri_vert(1,tt) == b ) then
              ee = 3
            else if ( tri_vert(2,tt) == b ) then
              ee = 1
            else
              ee = 2
            end if

          end do

          tri_nabe(ee,tt) = l

        end if

      end if

    end do
  end subroutine swapec

  subroutine vbedg ( x, y, point_num, point_xy, tri_num, tri_vert, tri_nabe, &
    ltri, ledg, rtri, redg ) &
        bind(C, name="vbedg")

  !*****************************************************************************80
  !
  !! VBEDG determines which boundary edges are visible to a point.
  !
  !  Discussion:
  !
  !    The point (X,Y) is assumed to be outside the convex hull of the
  !    region covered by the 2D triangulation.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    25 August 2001
  !
  !  Author:
  !
  !    Original FORTRAN77 version by Barry Joe.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Barry Joe,
  !    GEOMPACK - a software package for the generation of meshes
  !    using geometric algorithms,
  !    Advances in Engineering Software,
  !    Volume 13, pages 325-331, 1991.
  !
  !  Parameters:
  !
  !    Input, real(dp) X, Y, the coordinates of a point outside
  !    the convex hull of the current triangulation.
  !
  !    Input, integer(ip) POINT_NUM, the number of points.
  !
  !    Input, real(dp) POINT_XY(2,POINT_NUM), the coordinates
  !    of the vertices.
  !
  !    Input, integer(ip) TRI_NUM, the number of triangles.
  !
  !    Input, integer(ip) TRI_VERT(3,TRI_NUM), the triangle incidence
  !    list.
  !
  !    Input, integer(ip) TRI_NABE(3,TRI_NUM), the triangle neighbor
  !    list; negative values are used for links of a counter clockwise linked
  !    list of boundary edges;
  !      LINK = -(3*I + J-1) where I, J = triangle, edge index.
  !
  !    Input/output, integer(ip) LTRI, LEDG.  If LTRI /= 0 then these
  !    values are assumed to be already computed and are not changed, else they
  !    are updated.  On output, LTRI is the index of boundary triangle to the
  !    left of the leftmost boundary triangle visible from (X,Y), and LEDG is
  !    the boundary edge of triangle LTRI to the left of the leftmost boundary
  !    edge visible from (X,Y).  1 <= LEDG <= 3.
  !
  !    Input/output, integer(ip) RTRI.  On input, the index of the
  !    boundary triangle to begin the search at.  On output, the index of the
  !    rightmost boundary triangle visible from (X,Y).
  !
  !    Input/output, integer(ip) REDG, the edge of triangle RTRI that
  !    is visible from (X,Y).  1 <= REDG <= 3.
  !

    integer(ip), intent(in), value :: point_num                          !! number of points
    integer(ip), intent(in), value :: tri_num                            !! number of triangles

    integer(ip) :: a
    integer(ip) :: b
    integer(ip) :: e
    integer(ip) :: i4_wrap
    integer(ip) :: l
    logical :: ldone
    integer(ip), intent(inout) :: ledg
    integer(ip) :: lr
    integer(ip) :: lrline
    integer(ip), intent(inout) :: ltri
    real(dp), intent(in) :: point_xy(2,point_num)
    integer(ip), intent(inout) :: redg
    integer(ip), intent(inout) :: rtri
    integer(ip) :: t
    integer(ip), intent(in) :: tri_nabe(3,tri_num)
    integer(ip), intent(in) :: tri_vert(3,tri_num)
    real(dp), intent(in), value :: x                                     !! x-coordinate of test point
    real(dp), intent(in), value :: y                                     !! y-coordinate of test point
  !
  !  Find the rightmost visible boundary edge using links, then possibly
  !  leftmost visible boundary edge using triangle neighbor information.
  !
    if ( ltri == 0 ) then
      ldone = .false.
      ltri = rtri
      ledg = redg
    else
      ldone = .true.
    end if

    do

      l = -tri_nabe(redg,rtri)
      t = l / 3
      e = mod ( l, 3 ) + 1
      a = tri_vert(e,t)

      if ( e <= 2 ) then
        b = tri_vert(e+1,t)
      else
        b = tri_vert(1,t)
      end if

      lr = lrline ( x, y, point_xy(1,a), point_xy(2,a), point_xy(1,b), &
        point_xy(2,b), 0.0_dp )

      if ( lr <= 0 ) then
        exit
      end if

      rtri = t
      redg = e

    end do

    if ( ldone ) then
    end if

    t = ltri
    e = ledg

    do

      b = tri_vert(e,t)
      e = i4_wrap ( e-1, 1, 3 )

      do while ( 0 < tri_nabe(e,t) )

        t = tri_nabe(e,t)

        if ( tri_vert(1,t) == b ) then
          e = 3
        else if ( tri_vert(2,t) == b ) then
          e = 1
        else
          e = 2
        end if

      end do

      a = tri_vert(e,t)

      lr = lrline ( x, y, point_xy(1,a), point_xy(2,a), point_xy(1,b), &
         point_xy(2,b), 0.0_dp )

      if ( lr <= 0 ) then
        exit
      end if

    end do

    ltri = t
    ledg = e
  end subroutine vbedg

end module table_delaunay_mod
