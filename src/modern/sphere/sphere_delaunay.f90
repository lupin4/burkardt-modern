!> sphere_delaunay — Modern Fortran 2018
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
  character ( len = 255 ) delaunay_plot_file_name
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
  character ( len = 255 ) trplot_title
  double precision vlat
  double precision vlon
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
!  Plot the portion of the triangulation contained 
!  in the hemisphere centered at E = (ELAT,ELON), where ELAT and ELON
!  are taken to be the center of the range of
!  the nodal latitudes and longitudes.
!
  elat = 0.0D+00
  elon = 0.0D+00
  a = 90.0D+00
  numbr = ( n <= 200 )

  trplot_title = '(' // trim ( point_file_name ) // ')'

  delaunay_plot_file_name = point_file_name
  call file_name_ext_swap ( delaunay_plot_file_name, 'eps' )

  call get_unit ( iunit )

  open ( unit = iunit, file = delaunay_plot_file_name )

  call trplot ( iunit, pltsiz, elat, elon, a, n, x, y, z, list, &
    lptr, lend, trplot_title, numbr, ierror )

  close ( unit = iunit )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACK_INTERFACE - Warning!'
    write ( *, '(a,i8)' ) '  TRPLOT returned error code ', ierror
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRPLOT created the triangulation plot file: "' // &
    trim ( delaunay_plot_file_name ) // '".'
!
!  Write the XYZL file that indexes the points that form Delaunay 
!  triangulation lines.
!
  call tr_to_xyzl ( n, x, y, z, list, lptr, lend, point_file_name )
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

subroutine tr_to_xyzl ( n, x, y, z, list, lptr, lend, point_file_name )

!*****************************************************************************80
!
!! TR_TO_XYZL makes an XYZL file of Delaunay triangulation data.
!
!  Discussion:
!
!    The XYZL file, in combination with the original XYZ file, can be used
!    to draw the lines that form the Delaunay triangulation.
!
!    The XYZL file is simply a list of point indices from the XYZ file
!    which are to be connected.  The termination of a line is
!    indicated by an index value of -1.
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
!    Input, double precision X(N), Y(N), Z(N), the coordinates of the nodes.
!
!    Input, integer LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N),
!    information defining the triangulation, created by TRMESH.
!
!    Input, character ( len = * ) POINT_FILE_NAME, the name of the
!    file containing the point coordinates.
!
  implicit none

  integer n

  character ( len = 80 ) delaunay_file_name
  logical, parameter :: header = .false.
  integer iunit
  integer lend(n)
  integer list(6*(n-2))
  integer lp
  integer lpl
  integer lptr(6*(n-2))
  integer node1
  integer node2
  character ( len = *  ) point_file_name
  character ( len = 40 ) string
  double precision x(n)
  double precision y(n)
  double precision z(n)

  delaunay_file_name = point_file_name

  call file_name_ext_swap ( delaunay_file_name, 'xyzl' )

  call get_unit ( iunit )

  open ( unit = iunit, file = delaunay_file_name )
!
!  Write header.
!
  if ( header ) then
    write ( iunit, '(a)' ) '# "' // trim ( delaunay_file_name ) // '".'
    write ( iunit, '(a)' ) '# created by SPHERE_DELAUNAY.F90.'
    write ( iunit, '(a)' ) '#'
    write ( iunit, '(a)' ) &
      '# This file contains the indices of points to be connected'
    write ( iunit, '(a)' ) '# by lines, to form a Delaunay triangulation.'
    write ( iunit, '(a)' ) &
      '# The points are in the file "' // trim ( point_file_name ) // '".'
    write ( iunit, '(a)' ) '#'
    write ( iunit, '(a)' ) &
      '# This file lists the indices of points to be connected'
    write ( iunit, '(a)' ) '# by lines.  A line may include several points.'
    write ( iunit, '(a)' ) '# Each line is terminated by an index of -1.'
    write ( iunit, '(a)' ) '#'
  end if
!
!  List all the line segments that emanate from a point.
!  This will involve listing some segments twice.
!
  do node1 = 1, n

    lpl = lend(node1)
    lp = lpl

    do

      lp = lptr(lp)
      node2 = abs ( list(lp) )

      write ( iunit, '(i8)' ) node1
      write ( iunit, '(i8)' ) node2
      write ( iunit, '(i8)' ) -1

      if ( lp == lpl ) then
        exit
      end if

    end do

  end do

  close ( unit = iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TR_TO_XYZL:'
  write ( *, '(a)' ) '  Wrote the Delaunay XYZL file "' // &
    trim ( delaunay_file_name ) // '".'
end
