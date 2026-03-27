!> delaunay_lmap_2d — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module delaunay_lmap_2d_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: diaedg_lmap, dtris2_lmap, i4_modp, i4_wrap, i4vec_indicator, i4vec_sort_heap_index_a
  public :: lrline, perm_inv, r82vec_permute, r82vec_sort_heap_index_a, s_index_last, swapec_lmap
  public :: transform_lmap, triangulation_plot_eps, vbedg

contains

  function diaedg_lmap ( x0, y0, x1, y1, x2, y2, x3, y3, matrix ) &
        bind(C, name="diaedg_lmap")

  !*****************************************************************************80
  !
  !! DIAEDG_LMAP chooses a diagonal edge under a linear map.
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
  !    16 April 2004
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
  !    Input, real(dp) MATRIX(2,2), a linear map to be implicitly applied
  !    to the points.
  !
  !    Output, integer(ip) DIAEDG_LMAP, chooses a diagonal:
  !    +1, if diagonal edge 02 is chosen;
  !    -1, if diagonal edge 13 is chosen;
  !     0, if the four vertices are cocircular.
  !

    real(dp) :: ca
    real(dp) :: cb
    integer(ip) :: diaedg_lmap
    real(dp) :: dx10
    real(dp) :: dx12
    real(dp) :: dx30
    real(dp) :: dx32
    real(dp) :: dy10
    real(dp) :: dy12
    real(dp) :: dy30
    real(dp) :: dy32
    real(dp), intent(in) :: matrix(2,2)
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

    call transform_lmap ( matrix, dx10, dy10 )
    call transform_lmap ( matrix, dx12, dy12 )
    call transform_lmap ( matrix, dx30, dy30 )
    call transform_lmap ( matrix, dx32, dy32 )

    tola = tol * max ( abs ( dx10 ), abs ( dy10 ), abs ( dx30 ), abs ( dy30 ) )
    tolb = tol * max ( abs ( dx12 ), abs ( dy12 ), abs ( dx32 ), abs ( dy32 ) )

    ca = dx10 * dx30 + dy10 * dy30
    cb = dx12 * dx32 + dy12 * dy32

    if ( tola < ca .and. tolb < cb ) then

      diaedg_lmap = -1

    else if ( ca < -tola .and. cb < -tolb ) then

      diaedg_lmap = 1

    else

      tola = max ( tola, tolb )
      s = ( dx10 * dy30 - dx30 * dy10 ) * cb + ( dx32 * dy12 - dx12 * dy32 ) * ca

      if ( tola < s ) then
        diaedg_lmap = -1
      else if ( s < -tola ) then
        diaedg_lmap = 1
      else
        diaedg_lmap = 0
      end if

    end if
  end function diaedg_lmap

  subroutine dtris2_lmap ( point_num, point_xy, matrix, tri_num, tri_vert, &
    tri_nabe ) &
        bind(C, name="dtris2_lmap")

  !*****************************************************************************80
  !
  !! DTRIS2_LMAP constructs a Delaunay triangulation of 2D linear mapped vertices.
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
  !    16 April 2004
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
  !    Input, real(dp) POINT_XY(2,POINT_NUM), the coordinates
  !    of the vertices.
  !
  !    Input, real(dp) MATRIX(2,2), the linear map which should
  !    be implicitly applied to the points before the triangulation is done.
  !
  !    Output, integer(ip) TRI_NUM, the number of triangles in the
  !    triangulation; TRI_NUM is equal to 2*POINT_NUM - NB - 2, where NB is the
  !    number of boundary vertices.
  !
  !    Output, integer(ip) TRI_VERT(3,TRI_NUM), the nodes that make up
  !    each triangle.  The elements are indices of POINT_XY.  The vertices of
  !    the triangles are in counter clockwise order.
  !
  !    Output, integer(ip) TRI_NABE(3,TRI_NUM), the triangle neighbor
  !    list.  Positive elements are indices of TIL; negative elements are used
  !    for links of a counter clockwise linked list of boundary edges;
  !    LINK = -(3*I + J-1) where I, J = triangle, edge index; TRI_NABE(J,I)
  !    refers to the neighbor along edge from vertex J to J+1 (mod 3).
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
    real(dp), intent(in) :: matrix(2,2)
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

    call r82vec_permute ( point_num, point_xy, indx )
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
        write ( *, '(a)' ) 'DTRIS2_LMAP - Fatal error!'
        write ( *, '(a,i8)' ) '  Fails for point number I = ', i
        write ( *, '(a,i8)' ) '  M = ', m
        write ( *, '(a,i8)' ) '  M1 = ', m1
        write ( *, '(a,2g14.6)' ) '  X,Y(M)  = ', point_xy(1,m), point_xy(2,m)
        write ( *, '(a,2g14.6)' ) '  X,Y(M1) = ', point_xy(1,m1), point_xy(2,m1)
        ierr = 224
        stop
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
        write ( *, '(a)' ) 'DTRIS2_LMAP - Fatal error!'
        ierr = 225
        stop
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
          write ( *, '(a)' ) 'DTRIS2_LMAP - Fatal error!'
          write ( *, '(a)' ) '  Stack overflow.'
          stop
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

      call swapec_lmap ( m, matrix, top, ltri, ledg, point_num, point_xy, &
        tri_num, tri_vert, tri_nabe, stack, ierr )

      if ( ierr /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DTRIS2_LMAP - Fatal error!'
        write ( *, '(a)' ) '  Error return from SWAPEC_LMAP.'
        stop
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

    call perm_inv ( point_num, indx )

    call r82vec_permute ( point_num, point_xy, indx )
  end subroutine dtris2_lmap

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
  !  Example:
  !
  !        I     J     MOD  I4_MODP    Factorization
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

    if ( j == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'I4_MODP - Fatal error!'
      write ( *, '(a,i8)' ) '  I4_MODP ( I, J ) called with J = ', j
      stop
    end if

    i4_modp = mod ( i, j )

    if ( i4_modp < 0 ) then
      i4_modp = i4_modp + abs ( j )
    end if
  end function i4_modp

  function i4_wrap ( ival, ilo, ihi ) &
        bind(C, name="i4_wrap")

  !*****************************************************************************80
  !
  !! I4_WRAP forces an I4 to lie between given limits by wrapping.
  !
  !  Example:
  !
  !    ILO = 4, IHI = 8
  !
  !    I  I4_WRAP
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
  !    Input, integer(ip) ILO, IHI, the desired bounds for the value.
  !
  !    Output, integer(ip) I4_WRAP, a "wrapped" version of IVAL.
  !

    integer(ip) :: i4_modp
    integer(ip) :: i4_wrap
    integer(ip), intent(in), value :: ihi                                !! upper bound
    integer(ip), intent(in), value :: ilo                                !! lower bound
    integer(ip), intent(in), value :: ival                               !! value to wrap
    integer(ip) :: jhi
    integer(ip) :: jlo
    integer(ip) :: wide

    jlo = min ( ilo, ihi )
    jhi = max ( ilo, ihi )

    wide = jhi - jlo + 1

    if ( wide == 1 ) then
      i4_wrap = jlo
    else
      i4_wrap = jlo + i4_modp ( ival - jlo, wide )
    end if
  end function i4_wrap

  subroutine i4vec_indicator ( n, a ) &
        bind(C, name="i4vec_indicator")

  !*****************************************************************************80
  !
  !! I4VEC_INDICATOR sets an I4VEC to the indicator vector.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    09 November 2000
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

  subroutine i4vec_sort_heap_index_a ( n, a, indx ) &
        bind(C, name="i4vec_sort_heap_index_a")

  !*****************************************************************************80
  !
  !! I4VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an I4VEC.
  !
  !  Discussion:
  !
  !    The sorting is not actually carried out.  Rather an index array is
  !    created which defines the sorting.  This array may be used to sort
  !    or index the array, or to sort or index related arrays keyed on the
  !    original array.
  !
  !    Once the index array is computed, the sorting can be carried out
  !    "implicitly:
  !
  !      A(INDX(I)), I = 1 to N is sorted,
  !
  !    or explicitly, by the call
  !
  !      call I4VEC_PERMUTE ( N, A, INDX )
  !
  !    after which A(I), I = 1 to N is sorted.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    25 September 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of entries in the array.
  !
  !    Input, integer(ip) A(N), an array to be index-sorted.
  !
  !    Output, integer(ip) INDX(N), the sort index.  The
  !    I-th element of the sorted array is A(INDX(I)).
  !

    integer(ip), intent(in), value :: n                                  !! number of entries

    integer(ip), intent(in) :: a(n)
    integer(ip) :: aval
    integer(ip) :: i
    integer(ip), intent(out) :: indx(n)
    integer(ip) :: indxt
    integer(ip) :: ir
    integer(ip) :: j
    integer(ip) :: l

    if ( n <= 1 ) then
    end if

    do i = 1, n
      indx(i) = i
    end do

    l = n / 2 + 1
    ir = n

    do

      if ( 1 < l ) then

        l = l - 1
        indxt = indx(l)
        aval = a(indxt)

      else

        indxt = indx(ir)
        aval = a(indxt)
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
          if ( a(indx(j)) < a(indx(j+1)) ) then
            j = j + 1
          end if
        end if

        if ( aval < a(indx(j)) ) then
          indx(i) = indx(j)
          i = j
          j = j + j
        else
          j = ir + 1
        end if

      end do

      indx(i) = indxt

    end do
  end subroutine i4vec_sort_heap_index_a

  function lrline ( xu, yu, xv1, yv1, xv2, yv2, dv ) &
        bind(C, name="lrline")

  !*****************************************************************************80
  !
  !! LRLINE determines if a point is left of, right or, or on a directed line.
  !
  !  Discussion:
  !
  !    The directed line is paralled to, and at a signed distance DV from
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

  subroutine perm_inv ( n, p ) &
        bind(C, name="perm_inv")

  !*****************************************************************************80
  !
  !! PERM_INV inverts a permutation "in place".
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    25 July 2000
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

    integer(ip) :: i
    integer(ip) :: i0
    integer(ip) :: i1
    integer(ip) :: i2
    integer(ip) :: is
    integer(ip), intent(inout) :: p(n)

    if ( n <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PERM_INV - Fatal error!'
      write ( *, '(a,i8)' ) '  Input value of N = ', n
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

      is = -sign ( 1, p(i) )
      p(i) = sign ( p(i), is )

    end do

    do i = 1, n

      i1 = -p(i)

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
  end subroutine perm_inv

  subroutine r82vec_permute ( n, a, p ) &
        bind(C, name="r82vec_permute")

  !*****************************************************************************80
  !
  !! R82VEC_PERMUTE permutes an R82VEC in place.
  !
  !  Discussion:
  !
  !    This routine permutes an array of real "objects", but the same
  !    logic can be used to permute an array of objects of any arithmetic
  !    type, or an array of objects of any complexity.  The only temporary
  !    storage required is enough to store a single object.  The number
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
  !    11 January 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of objects.
  !
  !    Input/output, real(dp) A(2,N), the array to be permuted.
  !
  !    Input, integer(ip) P(N), the permutation.  P(I) = J means
  !    that the I-th element of the output array should be the J-th
  !    element of the input array.  P must be a legal permutation
  !    of the integers from 1 to N, otherwise the algorithm will
  !    fail catastrophically.
  !

    integer(ip), intent(in), value :: n                                  !! number of objects

    real(dp), intent(inout) :: a(2,n)
    real(dp) :: a_temp(2)
    integer(ip) :: iget
    integer(ip) :: iput
    integer(ip) :: istart
    integer(ip), intent(inout) :: p(n)
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

        a_temp(1:2) = a(1:2,istart)
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
            stop
          end if

          if ( iget == istart ) then
            a(1:2,iput) = a_temp(1:2)
            exit
          end if

          a(1:2,iput) = a(1:2,iget)

        end do

      end if

    end do
  !
  !  Restore the signs of the entries.
  !
    p(1:n) = -p(1:n)
  end subroutine r82vec_permute

  subroutine r82vec_sort_heap_index_a ( n, a, indx ) &
        bind(C, name="r82vec_sort_heap_index_a")

  !*****************************************************************************80
  !
  !! R82VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R82VEC.
  !
  !  Discussion:
  !
  !    The sorting is not actually carried out.  Rather an index array is
  !    created which defines the sorting.  This array may be used to sort
  !    or index the array, or to sort or index related arrays keyed on the
  !    original array.
  !
  !    Once the index array is computed, the sorting can be carried out
  !    "implicitly:
  !
  !      A(1:2,INDX(I)), I = 1 to N is sorted,
  !
  !    or explicitly, by the call
  !
  !      call R82VEC_PERMUTE ( N, A, INDX )
  !
  !    after which A(1:2,I), I = 1 to N is sorted.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    11 January 2004
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

    real(dp), intent(in) :: a(2,n)
    real(dp) :: aval(2)
    integer(ip) :: i
    integer(ip), intent(out) :: indx(n)
    integer(ip) :: indxt
    integer(ip) :: ir
    integer(ip) :: j
    integer(ip) :: l

    if ( n < 1 ) then
    end if

    if ( n == 1 ) then
      indx(1) = 1
    end if

    call i4vec_indicator ( n, indx )

    l = n / 2 + 1
    ir = n

    do

      if ( 1 < l ) then

        l = l - 1
        indxt = indx(l)
        aval(1:2) = a(1:2,indxt)

      else

        indxt = indx(ir)
        aval(1:2) = a(1:2,indxt)
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

  function s_index_last ( s, sub ) &
        bind(C, name="s_index_last")

  !*****************************************************************************80
  !
  !! S_INDEX_LAST finds the LAST occurrence of a given substring.
  !
  !  Discussion:
  !
  !    It returns the location in the string at which the substring SUB is
  !    first found, or 0 if the substring does not occur at all.
  !
  !    The routine is also trailing blank insensitive.  This is very
  !    important for those cases where you have stored information in
  !    larger variables.  If S is of length 80, and SUB is of
  !    length 80, then if S = 'FRED' and SUB = 'RED', a match would
  !    not be reported by the standard FORTRAN INDEX, because it treats
  !    both variables as being 80 characters long!  This routine assumes that
  !    trailing blanks represent garbage!
  !
  !    This means that this routine cannot be used to find, say, the last
  !    occurrence of a substring 'A ', since it assumes the blank space
  !    was not specified by the user, but is, rather, padding by the
  !    system.  However, as a special case, this routine can properly handle
  !    the case where either S or SUB is all blanks.
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
  !    Input, character ( len = * ) S, the string to be searched.
  !
  !    Input, character ( len = * ) SUB, the substring to search for.
  !
  !    Output, integer(ip) S_INDEX_LAST.  0 if SUB does not occur in
  !    the string.  Otherwise S_INDEX_LAST = I, where S(I:I+LENS-1) = SUB,
  !    where LENS is the length of SUB, and is the last place
  !    this happens.
  !

    integer(ip) :: i
    integer(ip) :: j
    integer(ip) :: llen1
    integer(ip) :: llen2
    character ( len = * ), intent(in) :: s                               !! string to search
    integer(ip) :: s_index_last
    character ( len = * ), intent(in) :: sub                             !! substring to find

    s_index_last = 0

    llen1 = len_trim ( s )
    llen2 = len_trim ( sub )
  !
  !  In case S or SUB is blanks, use LEN
  !
    if ( llen1 == 0 ) then
      llen1 = len ( s )
    end if

    if ( llen2 == 0 ) then
      llen2 = len ( sub )
    end if

    if ( llen1 < llen2 ) then
    end if

    do j = 1, llen1+1-llen2

      i = llen1 + 2 - llen2 - j

      if ( s(i:i+llen2-1) == sub ) then
        s_index_last = i
      end if

    end do
  end function s_index_last

  subroutine swapec_lmap ( i, matrix, top, btri, bedg, point_num, point_xy, &
    tri_num, tri_vert, tri_nabe, stack, ierr ) &
        bind(C, name="swapec_lmap")

  !*****************************************************************************80
  !
  !! SWAPEC_LMAP swaps diagonal edges until all triangles are Delaunay.
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
  !    Input, real(dp) MATRIX(2,2), the transformation matrix.
  !
  !    Input/output, integer(ip) TOP, the index of the top of the stack.
  !    On output, TOP is zero.
  !
  !    Input/output, integer(ip) BTRI, BEDG; on input, if positive, are
  !    the triangle and edge indices of a boundary edge whose updated indices
  !    must be recorded.  On output, these may be updated because of swaps.
  !
  !    Input, intger POINT_NUM, the number of points.
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
  !    Workspace, integer STACK(MAXST); on input, entries 1 through TOP
  !    contain the indices of initial triangles (involving vertex I)
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
    integer(ip) :: diaedg_lmap
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
    real(dp), intent(in) :: matrix(2,2)
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

      swap = diaedg_lmap ( x, y, point_xy(1,a), point_xy(2,a), point_xy(1,c), &
        point_xy(2,c), point_xy(1,b), point_xy(2,b), matrix )

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
  end subroutine swapec_lmap

  subroutine transform_lmap ( matrix, dx, dy ) &
        bind(C, name="transform_lmap")

  !*****************************************************************************80
  !
  !! TRANSFORM_LMAP applies a transformation matrix to a vector.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    16 April 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real(dp) MATRIX(2,2), the linear map to be applied
  !    to the data.
  !
  !    Input/output, real(dp) DX, DY, the components of the vector.
  !

    real(dp), intent(inout) :: dx                                        !! x-component of vector
    real(dp) :: dx2
    real(dp), intent(inout) :: dy                                        !! y-component of vector
    real(dp) :: dy2
    real(dp), intent(in) :: matrix(2,2)

    dx2 = matrix(1,1) * dx + matrix(1,2) * dy
    dy2 = matrix(2,1) * dx + matrix(2,2) * dy

    dx = dx2
    dy = dy2
  end subroutine transform_lmap

  subroutine triangulation_plot_eps ( file_name, node_num, node_x, node_y, &
    element_num, element_mask, element_node, title ) &
        bind(C, name="triangulation_plot_eps")

  !*****************************************************************************80
  !
  !! TRIANGULATION_PLOT_EPS creates an EPS file image of a triangulation.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    15 April 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) FILE_NAME, the name of the file to create.
  !
  !    Input, integer(ip) NODE_NUM, the number of nodes.
  !
  !    Input, real(dp) NODE_X(NODE_NUM), NODE_Y(NODE_NUM),
  !    the coordinates of the nodes.
  !
  !    Input, integer(ip) ELEMENT_NUM, the number of elements.
  !
  !    Input, logical ELEMENT_MASK(ELEMENT_NUM), a mask for the elements.
  !
  !    Input, integer(ip) ELEMENT_NODE(3,ELEMENT_NUM), the
  !    element->node data.
  !
  !    Input, character ( len = * ) TITLE, a title for the plot.
  !

    integer(ip), intent(in), value :: element_num                        !! number of elements
    integer(ip), intent(in), value :: node_num                           !! number of nodes

    real(dp) :: ave_x
    real(dp) :: ave_y
    integer(ip), parameter :: circle_size = 3
    real(dp) :: dif
    integer(ip) :: element
    logical, intent(in) :: element_mask(element_num)
    integer(ip), intent(in) :: element_node(3,element_num)
    integer(ip) :: eps_unit
    integer(ip) :: eps_x
    integer(ip) :: eps_y
    character ( len = * ), intent(in) :: file_name                       !! output file name
    integer(ip) :: i
    integer(ip) :: ios
    integer(ip) :: j
    integer(ip) :: local
    integer(ip) :: node
    logical :: node_mask(node_num)
    real(dp), intent(in) :: node_x(node_num)
    real(dp) :: node_x_max
    real(dp) :: node_x_min
    real(dp), intent(in) :: node_y(node_num)
    real(dp) :: node_y_max
    real(dp) :: node_y_min
    real(dp) :: scale
    character ( len = 40 ) :: string
    character ( len = * ), intent(in) :: title                           !! plot title
  !
  !  Determine the range of the unmasked elements.
  !
    node_x_min =  huge ( node_x_min )
    node_x_max = -huge ( node_x_max )
    node_y_min =  huge ( node_y_min )
    node_y_max = -huge ( node_y_max )

    node_mask(1:node_num) = .false.

    do element = 1, element_num
      if ( element_mask(element) ) then
        do j = 1, 3
          node = element_node(j,element)
          node_mask(node) = .true.
          node_x_min = min ( node_x_min, node_x(node) )
          node_x_max = max ( node_x_max, node_x(node) )
          node_y_min = min ( node_y_min, node_y(node) )
          node_y_max = max ( node_y_max, node_y(node) )
        end do
      end if
    end do

    if ( node_y_max - node_y_min < node_x_max - node_x_min ) then
      scale = node_x_max - node_x_min
      dif = ( node_x_max - node_x_min ) - ( node_y_max - node_y_min )
      node_y_max = node_y_max + 0.5_dp * dif
      node_y_min = node_y_min - 0.5_dp * dif
    else
      scale = node_y_max - node_y_min
      dif = ( node_y_max - node_y_min ) - ( node_x_max - node_x_min )
      node_x_max = node_x_max + 0.5_dp * dif
      node_x_min = node_x_min - 0.5_dp * dif
    end if

    call get_unit ( eps_unit )

    open ( unit = eps_unit, file = file_name, status = 'replace', iostat = ios )

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TRIANGULATION_PLOT_EPS - Fatal error!'
      write ( *, '(a)' ) '  Could not open the output EPS file.'
      stop
    end if

    write ( eps_unit, '(a)' ) '%!PS-Adobe-3.0 EPSF-3.0'
    write ( eps_unit, '(a)' ) &
      '%%Creator: triangulation_plot_eps(delaunay_lmap_2d.f90)'
    write ( eps_unit, '(a)' ) '%%Title: ' // trim ( file_name )
    write ( eps_unit, '(a)' ) '%%Pages: 1'
    write ( eps_unit, '(a)' ) '%%BoundingBox:    36    36   576   756'
    write ( eps_unit, '(a)' ) '%%Document-Fonts: Times-Roman'
    write ( eps_unit, '(a)' ) '%%LanguageLevel: 1'
    write ( eps_unit, '(a)' ) '%%EndComments'
    write ( eps_unit, '(a)' ) '%%BeginProlog'
    write ( eps_unit, '(a)' ) '/inch {72 mul} def'
    write ( eps_unit, '(a)' ) '%%EndProlog'
    write ( eps_unit, '(a)' ) '%%Page:      1     1'
    write ( eps_unit, '(a)' ) 'save'
    write ( eps_unit, '(a)' ) '%'
    write ( eps_unit, '(a)' ) '% Set RGB line color.'
    write ( eps_unit, '(a)' ) '%'
    write ( eps_unit, '(a)' ) ' 0.9000 0.9000 0.9000 setrgbcolor'
    write ( eps_unit, '(a)' ) '%'
    write ( eps_unit, '(a)' ) '% Draw a gray border around the page.'
    write ( eps_unit, '(a)' ) '%'
    write ( eps_unit, '(a)' ) 'newpath'
    write ( eps_unit, '(a)' ) '    36   126 moveto'
    write ( eps_unit, '(a)' ) '   576   126 lineto'
    write ( eps_unit, '(a)' ) '   576   666 lineto'
    write ( eps_unit, '(a)' ) '    36   666 lineto'
    write ( eps_unit, '(a)' ) '    36   126 lineto'
    write ( eps_unit, '(a)' ) 'stroke'
    write ( eps_unit, '(a)' ) '%'
    write ( eps_unit, '(a)' ) '% Set RGB line color.'
    write ( eps_unit, '(a)' ) '%'
    write ( eps_unit, '(a)' ) ' 0.0000 0.0000 0.0000 setrgbcolor'

    write ( eps_unit, '(a)' ) '%'
    write ( eps_unit, '(a)' ) '%  Label the plot:'
    write ( eps_unit, '(a)' ) '%'
    write ( eps_unit, '(a)' ) ' 0.0000 0.0000 0.0000 setrgbcolor'
    write ( eps_unit, '(a)' ) '/Times-Roman findfont 0.50 inch scalefont setfont'
    write ( eps_unit, '(a)' ) '    36   666 moveto'
    write ( eps_unit, '(a)' ) '(' // trim ( title ) // ') show'

    write ( eps_unit, '(a)' ) '%'
    write ( eps_unit, '(a)' ) '% Define a clipping polygon'
    write ( eps_unit, '(a)' ) '%'
    write ( eps_unit, '(a)' ) '    36   126 moveto'
    write ( eps_unit, '(a)' ) '   576   126 lineto'
    write ( eps_unit, '(a)' ) '   576   666 lineto'
    write ( eps_unit, '(a)' ) '    36   666 lineto'
    write ( eps_unit, '(a)' ) '    36   126 lineto'
    write ( eps_unit, '(a)' ) 'clip newpath'

    write ( eps_unit, '(a)' ) '%'
    write ( eps_unit, '(a)' ) '%  Draw filled dots at each node:'
    write ( eps_unit, '(a)' ) '%'
    write ( eps_unit, '(a)' ) ' 0.0000 0.0000 0.9000 setrgbcolor'

    do node = 1, node_num

      if ( node_mask(node) ) then

        eps_x = int &
          ( ( node_x_max - node_x(node)              ) *  61.0_dp   &
          + (            + node_x(node) - node_x_min ) * 551.0_dp ) &
          / scale

        eps_y = int &
          ( ( node_y_max - node_y(node)              ) * 151.0_dp   &
          + (              node_y(node) - node_y_min ) * 641.0_dp ) &
          / scale

        write ( eps_unit, '(a,i4,2x,i4,2x,i4,a)' ) &
          'newpath  ', eps_x, eps_y, circle_size, ' 0 360 arc closepath fill'

      end if

    end do
  !
  !  Label the nodes.
  !
    if ( node_num <= 200 ) then

      write ( eps_unit, '(a)' ) '%'
      write ( eps_unit, '(a)' ) '%  Label the nodes:'
      write ( eps_unit, '(a)' ) '%'
      write ( eps_unit, '(a)' ) ' 0.0000 0.0000 1.0000 setrgbcolor'
      write ( eps_unit, '(a)' ) &
        '/Times-Roman findfont 0.20 inch scalefont setfont'

      do node = 1, node_num

        if ( node_mask(node) ) then

          eps_x = int &
            ( ( node_x_max - node_x(node)              ) *  61.0_dp   &
            + (            + node_x(node) - node_x_min ) * 551.0_dp ) &
            / scale

          eps_y = int &
            ( ( node_y_max - node_y(node)              ) * 151.0_dp   &
            + (              node_y(node) - node_y_min ) * 641.0_dp ) &
            / scale

          write ( string, '(i4)' ) node
          string = adjustl ( string )

          write ( eps_unit, '(i4,2x,i4,a)' ) eps_x, eps_y+5, &
            ' moveto (' // trim ( string ) // ') show'

        end if

      end do

    end if

    write ( eps_unit, '(a)' ) '%'
    write ( eps_unit, '(a)' ) '%  Draw the element sides:'
    write ( eps_unit, '(a)' ) '%'
    write ( eps_unit, '(a)' ) ' 0.9000 0.0000 0.0000 setrgbcolor'

    do element = 1, element_num

      if ( .not. element_mask(element) ) then
        cycle
      end if

      local = 1
      node = element_node(local,element)

      eps_x = int &
        ( ( node_x_max - node_x(node)              ) *  61.0_dp   &
        + (            + node_x(node) - node_x_min ) * 551.0_dp ) &
        / scale

      eps_y = int &
        ( ( node_y_max - node_y(node)              ) * 151.0_dp   &
        + (              node_y(node) - node_y_min ) * 641.0_dp ) &
        / scale

      write ( eps_unit, '(a,i4,2x,i4,a)' ) 'newpath ', eps_x, eps_y, ' moveto'

      do

        local = mod ( local, 3 ) + 1
        node = element_node(local,element)

        eps_x = int &
          ( ( node_x_max - node_x(node)              ) *  61.0_dp   &
          + (            + node_x(node) - node_x_min ) * 551.0_dp ) &
          / scale

        eps_y = int &
          ( ( node_y_max - node_y(node)              ) * 151.0_dp   &
          + (              node_y(node) - node_y_min ) * 641.0_dp ) &
          / scale

        write ( eps_unit, '(i4,2x,i4,a)' ) eps_x, eps_y, ' lineto'

        if ( local == 1 ) then
          exit
        end if

      end do

      write ( eps_unit, '(a)' ) 'stroke'

    end do
  !
  !  Label the elements.
  !  (This is usually turned off here, because it clutters up the output.)
  !
    if ( element_num <= 20 ) then

      write ( eps_unit, '(a)' ) '%'
      write ( eps_unit, '(a)' ) '%  Label the elements:'
      write ( eps_unit, '(a)' ) '%'
      write ( eps_unit, '(a)' ) ' 1.0000 0.0000 0.0000 setrgbcolor'
      write ( eps_unit, '(a)' ) &
        '/Times-Roman findfont 0.30 inch scalefont setfont'

      do element = 1, element_num

        if ( .not. element_mask(element) ) then
          cycle
        end if

        ave_x = 0.0_dp
        ave_y = 0.0_dp

        do i = 1, 3

          node = element_node(i,element)

          ave_x = ave_x + node_x(node)
          ave_y = ave_y + node_y(node)

        end do

        ave_x = ave_x / 3.0_dp
        ave_y = ave_y / 3.0_dp

        eps_x = int &
          ( ( node_x_max - ave_x              ) *  61.0_dp   &
          + (            + ave_x - node_x_min ) * 551.0_dp ) &
          / scale

        eps_y = int &
          ( ( node_y_max - ave_y              ) * 151.0_dp   &
          + (              ave_y - node_y_min ) * 641.0_dp ) &
          / scale

        write ( string, '(i4)' ) element
        string = adjustl ( string )

        write ( eps_unit, '(i4,2x,i4,a)' ) eps_x, eps_y, ' moveto (' &
          // trim ( string ) // ') show'

      end do

    end if

    write ( eps_unit, '(a)' ) '%'
    write ( eps_unit, '(a)' ) 'restore showpage'
    write ( eps_unit, '(a)' ) '%'
    write ( eps_unit, '(a)' ) '% End of page'
    write ( eps_unit, '(a)' ) '%'
    write ( eps_unit, '(a)' ) '%%Trailer'
    write ( eps_unit, '(a)' ) '%%EOF'

    close ( unit = eps_unit )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULATION_PLOT_EPS:'
    write ( *, '(a)' ) '  An encapsulated PostScript file was created'
    write ( *, '(a)' ) '  containing an image of the triangulation.'
    write ( *, '(a)' ) '  The file is named "' // trim ( file_name ) // '".'
  end subroutine triangulation_plot_eps

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
  !    19 January 2009
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
  !  Modified:
  !
  !    25 August 2001
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
  !    Input, integer(ip) TRI_VERT(3,TRI_NUM), the triangle
  !    incidence list.
  !
  !    Input, integer(ip) TRI_NABE(3,TRI_NUM), the triangle neighbor
  !    list; negative values are used for links of a counter clockwise linked
  !    list of boundary edges;
  !      LINK = -(3*I + J-1) where I, J = triangle, edge index.
  !
  !    Input/output, integer(ip) LTRI, LEDG.  If LTRI /= 0 then these
  !    values are assumed to be already computed and are not changed, else they
  !    are updated.  On output, LTRI is the index of boundary triangle to the left
  !    of the leftmost boundary triangle visible from (X,Y), and LEDG is the
  !    boundary edge of triangle LTRI to the left of the leftmost boundary edge
  !    visible from (X,Y).  1 <= LEDG <= 3.
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
      e = i4_wrap ( e - 1, 1, 3 )

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

end module delaunay_lmap_2d_mod
