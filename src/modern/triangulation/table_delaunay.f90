!> table_delaunay — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

function diaedg ( x0, y0, x1, y1, x2, y2, x3, y3 )

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
!    Input, double precision X0, Y0, X1, Y1, X2, Y2, X3, Y3, the
!    coordinates of the vertices of a quadrilateral, given in
!    counter clockwise order.
!
!    Output, integer DIAEDG, chooses a diagonal:
!    +1, if diagonal edge 02 is chosen;
!    -1, if diagonal edge 13 is chosen;
!     0, if the four vertices are cocircular.
!
  implicit none

  double precision ca
  double precision cb
  integer diaedg
  double precision dx10
  double precision dx12
  double precision dx30
  double precision dx32
  double precision dy10
  double precision dy12
  double precision dy30
  double precision dy32
  double precision s
  double precision tol
  double precision tola
  double precision tolb
  double precision x0
  double precision x1
  double precision x2
  double precision x3
  double precision y0
  double precision y1
  double precision y2
  double precision y3

  tol = 100.0D+00 * epsilon ( tol )

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
end

subroutine dtris2 ( point_num, point_xy, tri_num, tri_vert, tri_nabe )

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
!    Input, integer POINT_NUM, the number of vertices.
!
!    Input/output, double precision POINT_XY(2,POINT_NUM), the coordinates 
!    of the vertices.  On output, the vertices have been sorted into 
!    dictionary order.
!
!    Output, integer TRI_NUM, the number of triangles in the 
!    triangulation; TRI_NUM is equal to 2*POINT_NUM - NB - 2, where NB is the 
!    number of boundary vertices.
!
!    Output, integer TRI_VERT(3,TRI_NUM), the nodes that make up
!    each triangle.  The elements are indices of POINT_XY.  The vertices of the 
!    triangles are in counter clockwise order.
!
!    Output, integer TRI_NABE(3,TRI_NUM), the triangle neighbor 
!    list.  Positive elements are indices of TIL; negative elements are used 
!    for links of a counter clockwise linked list of boundary edges; 
!    LINK = -(3*I + J-1) where I, J = triangle, edge index; TRI_NABE(J,I) refers
!    to the neighbor along edge from vertex J to J+1 (mod 3).
!
  implicit none

  integer point_num

  double precision cmax
  integer e
  integer i
  integer ierr
  integer indx(point_num)
  integer j
  integer k
  integer l
  integer ledg
  integer lr
  integer lrline
  integer ltri
  integer m
  integer m1
  integer m2
  integer n
  double precision point_xy(2,point_num)
  integer redg
  integer rtri
  integer stack(point_num)
  integer t
  double precision tol
  integer top
  integer tri_nabe(3,point_num*2)
  integer tri_num
  integer tri_vert(3,point_num*2)

  tol = 100.0D+00 * epsilon ( tol )

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

      if ( tol * ( cmax + 1.0D+00 ) &
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
      point_xy(2,m1), point_xy(1,m2), point_xy(2,m2), 0.0D+00 )

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
      point_xy(2,m1), point_xy(1,m2), point_xy(2,m2), 0.0D+00 )

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
end

function lrline ( xu, yu, xv1, yv1, xv2, yv2, dv )

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
!    Input, double precision XU, YU, the coordinates of the point whose
!    position relative to the directed line is to be determined.
!
!    Input, double precision XV1, YV1, XV2, YV2, the coordinates of two points
!    that determine the directed base line.
!
!    Input, double precision DV, the signed distance of the directed line
!    from the directed base line through the points (XV1,YV1) and (XV2,YV2).
!    DV is positive for a line to the left of the base line.
!
!    Output, integer LRLINE, the result:
!    +1, the point is to the right of the directed line;
!     0, the point is on the directed line;
!    -1, the point is to the left of the directed line.
!
  implicit none

  double precision dv
  double precision dx
  double precision dxu
  double precision dy
  double precision dyu
  integer lrline
  double precision t
  double precision tol
  double precision tolabs
  double precision xu
  double precision xv1
  double precision xv2
  double precision yu
  double precision yv1
  double precision yv2

  tol = 100.0D+00 * epsilon ( tol )

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
end

subroutine swapec ( i, top, btri, bedg, point_num, point_xy, tri_num, &
  tri_vert, tri_nabe, stack, ierr )

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
!    Input, integer I, the index of the new vertex.
!
!    Input/output, integer TOP, the index of the top of the stack.
!    On output, TOP is zero.
!
!    Input/output, integer BTRI, BEDG; on input, if positive, are 
!    the triangle and edge indices of a boundary edge whose updated indices
!    must be recorded.  On output, these may be updated because of swaps.
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, double precision POINT_XY(2,POINT_NUM), the coordinates
!    of the points.
!
!    Input, integer TRI_NUM, the number of triangles.
!
!    Input/output, integer TRI_VERT(3,TRI_NUM), the triangle 
!    incidence list.  May be updated on output because of swaps.
!
!    Input/output, integer TRI_NABE(3,TRI_NUM), the triangle 
!    neighbor list; negative values are used for links of the counter-clockwise 
!    linked list of boundary edges;  May be updated on output because of swaps.
!      LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Workspace, integer STACK(MAXST); on input, entries 1 through
!    TOP contain the indices of initial triangles (involving vertex I)
!    put in stack; the edges opposite I should be in interior;  entries
!    TOP+1 through MAXST are used as a stack.
!
!    Output, integer IERR is set to 8 for abnormal return.
!
  implicit none

  integer point_num
  integer tri_num

  integer a
  integer b
  integer bedg
  integer btri
  integer c
  integer diaedg
  integer e
  integer ee
  integer em1
  integer ep1
  integer f
  integer fm1
  integer fp1
  integer i
  integer ierr
  integer i4_wrap
  integer l
  integer r
  integer s
  integer stack(point_num)
  integer swap
  integer t
  integer top
  integer tri_nabe(3,tri_num)
  integer tri_vert(3,tri_num)
  integer tt
  integer u
  double precision point_xy(2,point_num)
  double precision x
  double precision y
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
end

subroutine vbedg ( x, y, point_num, point_xy, tri_num, tri_vert, tri_nabe, &
  ltri, ledg, rtri, redg )

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
!    Input, double precision X, Y, the coordinates of a point outside
!    the convex hull of the current triangulation.
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, double precision POINT_XY(2,POINT_NUM), the coordinates 
!    of the vertices.
!
!    Input, integer TRI_NUM, the number of triangles.
!
!    Input, integer TRI_VERT(3,TRI_NUM), the triangle incidence 
!    list.
!
!    Input, integer TRI_NABE(3,TRI_NUM), the triangle neighbor 
!    list; negative values are used for links of a counter clockwise linked 
!    list of boundary edges;
!      LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Input/output, integer LTRI, LEDG.  If LTRI /= 0 then these 
!    values are assumed to be already computed and are not changed, else they 
!    are updated.  On output, LTRI is the index of boundary triangle to the 
!    left of the leftmost boundary triangle visible from (X,Y), and LEDG is 
!    the boundary edge of triangle LTRI to the left of the leftmost boundary
!    edge visible from (X,Y).  1 <= LEDG <= 3.
!
!    Input/output, integer RTRI.  On input, the index of the 
!    boundary triangle to begin the search at.  On output, the index of the 
!    rightmost boundary triangle visible from (X,Y).
!
!    Input/output, integer REDG, the edge of triangle RTRI that 
!    is visible from (X,Y).  1 <= REDG <= 3.
!
  implicit none

  integer point_num
  integer tri_num

  integer a
  integer b
  integer e
  integer i4_wrap
  integer l
  logical ldone
  integer ledg
  integer lr
  integer lrline
  integer ltri
  double precision point_xy(2,point_num)
  integer redg
  integer rtri
  integer t
  integer tri_nabe(3,tri_num)
  integer tri_vert(3,tri_num)
  double precision x
  double precision y
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
      point_xy(2,b), 0.0D+00 )

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
       point_xy(2,b), 0.0D+00 )

    if ( lr <= 0 ) then
      exit
    end if

  end do

  ltri = t
  ledg = e
end
