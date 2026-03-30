!> sphere_voronoi — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

subroutine stripack_interface ( point_file_name )

!*****************************************************************************80
!
!! STRIPACK_INTERFACE calls STRIPACK routines.
!
!  Discussion:
!
!    10 June 2002: Changed this routine so that it can handle any size
!    problem, by making all arrays allocatable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) POINT_FILE_NAME, the name of the input file.
!
  implicit none

  double precision a
  double precision , allocatable, dimension ( : ) :: ds
  double precision elat
  double precision elon
  integer i
  integer ierror
  integer iunit
  integer , allocatable, dimension ( : ) :: iwk
  integer k
  integer kt
  integer , allocatable, dimension ( :, : ) :: lbtri
  integer , allocatable, dimension ( : ) :: lend
  integer , allocatable, dimension ( : ) :: list
  integer , allocatable, dimension ( : ) :: listc
  integer lnew
  integer lp
  integer lpl
  integer , allocatable, dimension ( : ) :: lptr
  integer , allocatable, dimension ( :, : ) :: ltri
  integer n
  integer na
  integer nb
  integer nn
  double precision norm
  integer nt
  integer ntemp
  logical numbr
  integer nv
  double precision , parameter :: pltsiz = 7.5D+00
  character ( len = * ) point_file_name
  double precision , allocatable, dimension ( : ) :: rc
  integer side_max
  double precision vlat
  double precision vlon
  character ( len = 255 ) voronoi_plot_file_name
  character ( len = 255 ) voronoi_plot_title
  double precision , allocatable, dimension ( : ) :: x
  double precision , allocatable, dimension ( : ) :: xc
  double precision , allocatable, dimension ( : ) :: y
  double precision , allocatable, dimension ( : ) :: yc
  double precision , allocatable, dimension ( : ) :: z
  double precision , allocatable, dimension ( : ) :: zc
!
!  Count the number of lines of (X,Y,Z) data.
!
  call file_row_count ( point_file_name, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACK_INTERFACE - Fatal error!'
    write ( *, '(a)' ) '  The input file has no data.'
  end if
!
!  Allocate everything.
!
  allocate ( ds(1:n) )
  allocate ( iwk(1:2*n) )
  allocate ( lbtri(1:6,1:n) )
  allocate ( lend(1:n) )
  allocate ( list(1:6*(n-2)) )
  allocate ( listc(1:6*(n-2)) )
  allocate ( lptr(1:6*(n-2)) )
  allocate ( ltri(1:9,1:2*(n-2)) )
  allocate ( rc(1:2*(n-2)) )
  allocate ( x(1:n) )
  allocate ( xc(1:2*(n-2)) )
  allocate ( y(1:n) )
  allocate ( yc(1:2*(n-2)) )
  allocate ( z(1:n) )
  allocate ( zc(1:2*(n-2)) )
!
!  Read the (X,Y,Z) data from a file.
!
  call xyz_read ( point_file_name, n, x, y, z, ierror )
!
!  Make sure the data is on the unit sphere.
!
  do i = 1, n
    norm = sqrt ( x(i)**2 + y(i)**2 + z(i)**2 )
    x(i) = x(i) / norm
    y(i) = y(i) / norm
    z(i) = z(i) / norm
  end do
!
!  Create the triangulation.
!
  call trmesh ( n, x, y, z, list, lptr, lend, lnew, iwk, iwk(n+1), ds, ierror )

  if ( ierror == -2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACK_INTERFACE - Fatal error!'
    write ( *, '(a)' ) '  Error in TRMESH.'
    write ( *, '(a)' ) '  The first 3 nodes are collinear.'
    stop
  end if

  if ( 0 < ierror ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACK_INTERFACE - Fatal error!'
    write ( *, '(a)' ) '  Error in TRMESH.'
    write ( *, '(a)' ) '  Duplicate nodes encountered.'
    stop
  end if
!
!  Create a triangle list.
!
  call trlist ( n, list, lptr, lend, 9, nt, ltri, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACK_INTERFACE - Fatal error!'
    write ( *, '(a)' ) '  Error in TRLIST.'
    stop
  end if
!
!  Construct the Voronoi diagram.
!
!  Note that the triangulation data structure is altered if NB > 0.
!
  call crlist ( n, n, x, y, z, list, lend, lptr, lnew, &
    lbtri, listc, nb, xc, yc, zc, rc, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACK_INTERFACE - Fatal error!'
    write ( *, '(a)' ) '  Error in CRLIST.'
    write ( *, '(a,i8)' ) '  IERROR = ', ierror
    stop
  end if
!
!  Count the number of polygons of each size.
!
  call poly_count_max ( n, lend, lptr, listc, side_max )

  call poly_count ( n, lend, lptr, listc, side_max )
!
!  Plot the portion of the Voronoi diagram contained 
!  in the hemisphere centered at E = (ELAT,ELON), where ELAT and ELON
!  are taken to be the center of the range of
!  the nodal latitudes and longitudes.
!
  elat = 0.0D+00
  elon = 0.0D+00
  a = 90.0D+00
  numbr = ( n <= 200 )
  nt = 2 * n - 4

  voronoi_plot_file_name = point_file_name
  call file_name_ext_swap ( voronoi_plot_file_name, 'eps' )

  voronoi_plot_title = '(' // trim ( point_file_name ) // ')'

  call get_unit ( iunit )

  open ( unit = iunit, file = voronoi_plot_file_name )

  call vrplot ( iunit, pltsiz, elat, elon, a, n, x, y, z, nt, listc, &
    lptr, lend, xc, yc, zc, voronoi_plot_title, numbr, ierror )

  close ( unit = iunit )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACK_INTERFACE - Warning!'
    write ( *, '(a,i8)' ) '  VRPLOT returned error code ', ierror
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VRPLOT created the Voronoi plot file: "' // &
    trim ( voronoi_plot_file_name ) // '".'
!
!  Write the XYZ and XYZF files that define the Voronoi information.
!
  call vr_to_xyzf ( n, xc, yc, zc, listc, lptr, lend, point_file_name )
!
!  Free memory.
!
  deallocate ( ds )
  deallocate ( iwk )
  deallocate ( lbtri )
  deallocate ( lend )
  deallocate ( list )
  deallocate ( listc )
  deallocate ( lptr )
  deallocate ( ltri )
  deallocate ( rc )
  deallocate ( x )
  deallocate ( xc )
  deallocate ( y )
  deallocate ( yc )
  deallocate ( z )
  deallocate ( zc )
end

subroutine poly_count ( n, lend, lptr, listc, side_max )

!*****************************************************************************80
!
!! POLY_COUNT counts the number of polygons of each size in the diagram.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of Voronoi polygons.
!
!    Input, integer LEND(N), some kind of pointer.
!
!    Input, integer LPTR(6*(N-2)), some other kind of pointer.
!
!    Input, integer LISTC(6*(N-2)), some other kind of pointer.
!
!    Input, integer SIDE_MAX, the maximum polygonal order.
!
  implicit none

  integer n
  integer side_max

  integer count(0:side_max)
  integer edges
  integer i
  integer kv
  integer lend(n)
  integer listc(6*(n-2))
  integer lp
  integer lpl
  integer lptr(6*(n-2))
  integer n0
  integer sides
  integer vertices

  count(0:side_max) = 0

  edges = 0
  vertices = 0

  do n0 = 1, n

    lpl = lend(n0)

    lp = lpl

    sides = 0

    do

      lp = lptr(lp)
      kv = listc(lp)

      vertices = max ( vertices, kv )
      sides = sides + 1
      edges = edges + 1

      if ( lp == lpl ) then
        exit
      end if

    end do

    if ( sides < 3 ) then
      write ( *, * ) '  Polygon ', n0, ' has ', sides, ' sides.'
    end if

    if ( side_max < sides ) then
      write ( *, * ) '  Polygon ', n0, ' has too many sides!'
    end if

    if ( sides <= 0 ) then
      count(0) = count(0) + 1
    else if ( 0 < sides .and. sides <= side_max ) then
      count(sides) = count(sides) + 1
    end if

  end do

  edges = edges / 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POLY_COUNT'
  write ( *, '(a)' ) '  Number of polygons of each shape.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Faces =    ', n
  write ( *, '(a,i8)' ) '  Vertices = ', vertices
  write ( *, '(a,i8)' ) '  Edges =    ', edges
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Check Euler''s formula:'
  write ( *, '(a,i8)' ) '  F+V-E-2 =  ', n + vertices - edges - 2
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Sides    Number'
  write ( *, '(a)' ) ' '

  do i = 1, side_max
    if ( count(i) /= 0 ) then
      write ( *, '(2x,i8,2x,i8)' ) i, count(i)
    end if
  end do
end

subroutine poly_count_max ( n, lend, lptr, listc, side_max )

!*****************************************************************************80
!
!! POLY_COUNT_MAX determines the highest order polygon.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of Voronoi polygons.
!
!    Input, integer LEND(N), some kind of pointer.
!
!    Input, integer LPTR(6*(N-2)), some other kind of pointer.
!
!    Input, integer LISTC(6*(N-2)), some other kind of pointer.
!
!    Output, integer SIDE_MAX, the highest order polgon.
!
  implicit none

  integer n

  integer i
  integer lend(n)
  integer listc(6*(n-2))
  integer lp
  integer lpl
  integer lptr(6*(n-2))
  integer n0
  integer side_max
  integer sides

  side_max = - 1

  do n0 = 1, n

    lpl = lend(n0)

    lp = lpl

    sides = 0

    do

      lp = lptr(lp)

      sides = sides + 1

      if ( lp == lpl ) then
        exit
      end if

    end do

    side_max = max ( side_max, sides )

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POLY_COUNT_MAX'
  write ( *, '(a,i8)' ) '  Highest order polygon is ', side_max
end

subroutine vr_to_xyzf ( n, xc, yc, zc, listc, lptr, lend, point_file_name )

!*****************************************************************************80
!
!! VR_TO_XYZF makes an XYZF file of Voronoi diagram data.
!
!  Discussion:
!
!    Actually, this routine makes two files:
!
!    an XYZ file of the Voronoi vertices,
!    an XYZF file of the faces formed by the Voronoi vertices.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of nodes.
!
!    Input, double precision XC(2*(N-2)), YC(2*(N-2)), ZC(2*(N-2)),
!    the coordinates  of the Voronoi vertices.
!
!    Input, integer LISTC(6*(N-2)), LPTR(6*(N-2)), LEND(N), 
!    information defining the triangulation, created by TRMESH.
!
!    Input, character ( len = * ) POINT_FILE_NAME, the name of the
!    file containing the point coordinates.
!
  implicit none

  integer n

  logical, parameter :: header = .false.
  integer i
  integer iunit
  integer kv1
  integer kv2
  integer lend(n)
  integer listc(6*(n-2))
  integer lp
  integer lpl
  integer lptr(6*(n-2))
  integer node1
  integer node2
  character ( len = *  ) point_file_name
  character ( len = 40 ) string
  character ( len = 255 ) voronoi_face_file_name
  character ( len = 255 ) voronoi_vertex_file_name
  double precision xc(2*(n-2))
  double precision yc(2*(n-2))
  double precision zc(2*(n-2))
!
!  Write the Voronoi diagram vertex file.
!
  voronoi_vertex_file_name = point_file_name

  call file_name_ext_swap ( voronoi_vertex_file_name, 'voronoi.xyz' )

  call get_unit ( iunit )

  open ( unit = iunit, file = voronoi_vertex_file_name )

  if ( header ) then
    write ( iunit, '(a)' ) '# "' // trim ( voronoi_vertex_file_name ) // '".'
    write ( iunit, '(a)' ) '# created by SPHERE_VORONOI.F90.'
    write ( iunit, '(a)' ) '#'
    write ( iunit, '(a)' ) '# This file contains the coordinates of points that'
    write ( iunit, '(a)' ) '# are the vertices of faces in the Voronoi diagram.'
    write ( iunit, '(a)' ) '#'
  end if

  do i = 1, 2 * ( n - 2 )
    write ( iunit, '(3f8.4)' ) xc(i), yc(i), zc(i)
  end do
  close ( unit = iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VR_TO_XYZF:'
  write ( *, '(a)' ) '  Wrote the Voronoi XYZ file "' &
    // trim ( voronoi_vertex_file_name ) // '".'
!
!  Write the Voronoi diagram face file.
!
  voronoi_face_file_name = point_file_name

  call file_name_ext_swap ( voronoi_face_file_name, 'voronoi.xyzf' )

  call get_unit ( iunit )

  open ( unit = iunit, file = voronoi_face_file_name )

  if ( header ) then
    write ( iunit, '(a)' ) '# "' // trim ( voronoi_face_file_name ) // '".'
    write ( iunit, '(a)' ) '# created by SPHERE_VORONOI.F90.'
    write ( iunit, '(a)' ) '#'
    write ( iunit, '(a)' ) &
      '# This file contains the indices of Voronoi vertices that'
    write ( iunit, '(a)' ) '# form the faces of the Voronoi diagram.'
    write ( iunit, '(a)' ) '#'
    write ( iunit, '(a)' ) '# The coordinates of these vertices are stored in'
    write ( iunit, '(a)' ) '# "' // trim ( voronoi_vertex_file_name ) // '".'
    write ( iunit, '(a)' ) '#'
  end if

  do i = 1, n

    lpl = lend(i)

    kv2 = listc(lpl) 
    write ( iunit, '(i8)', advance = 'no' ) kv2

    lp = lpl

    do

      lp = lptr(lp)
      kv1 = kv2
      kv2 = listc(lp)

      write ( iunit, '(i8)', advance = 'no' ) kv2

      if ( lp == lpl ) then
        write ( iunit, '(i8)', advance = 'yes' ) -1
        exit
      end if

    end do

  end do

  close ( unit = iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VR_TO_XYZF:'
  write ( *, '(a)' ) '  Wrote the Voronoi XYZF file "'  &
    // trim ( voronoi_face_file_name ) // '".'
end
