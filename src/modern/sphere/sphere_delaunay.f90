!> sphere_delaunay -- Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module sphere_delaunay_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: stripack_interface, tr_to_xyzl, xyz_read

contains

  subroutine stripack_interface ( point_file_name ) &
        bind(C, name="stripack_interface")

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

    character ( len = * ), intent(in) :: point_file_name

    real(dp) :: a
    character ( len = 255 ) :: delaunay_plot_file_name
    real(dp), allocatable, dimension ( : ) :: ds
    real(dp) :: elat
    real(dp) :: elon
    integer(ip) :: i
    integer(ip) :: ierror
    integer(ip) :: iunit
    integer(ip), allocatable, dimension ( : ) :: iwk
    integer(ip) :: k
    integer(ip) :: kt
    integer(ip), allocatable, dimension ( :, : ) :: lbtri
    integer(ip), allocatable, dimension ( : ) :: lend
    integer(ip), allocatable, dimension ( : ) :: list
    integer(ip), allocatable, dimension ( : ) :: listc
    integer(ip) :: lnew
    integer(ip) :: lp
    integer(ip) :: lpl
    integer(ip), allocatable, dimension ( : ) :: lptr
    integer(ip), allocatable, dimension ( :, : ) :: ltri
    integer(ip) :: n
    integer(ip) :: na
    integer(ip) :: nb
    integer(ip) :: nn
    real(dp) :: norm
    integer(ip) :: nt
    integer(ip) :: ntemp
    logical :: numbr
    integer(ip) :: nv
    real(dp), parameter :: pltsiz = 7.5_dp
    real(dp), allocatable, dimension ( : ) :: rc
    integer(ip) :: side_max
    character ( len = 255 ) :: trplot_title
    real(dp) :: vlat
    real(dp) :: vlon
    real(dp), allocatable, dimension ( : ) :: x
    real(dp), allocatable, dimension ( : ) :: xc
    real(dp), allocatable, dimension ( : ) :: y
    real(dp), allocatable, dimension ( : ) :: yc
    real(dp), allocatable, dimension ( : ) :: z
    real(dp), allocatable, dimension ( : ) :: zc
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
    elat = 0.0_dp
    elon = 0.0_dp
    a = 90.0_dp
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
  end subroutine stripack_interface

  subroutine tr_to_xyzl ( n, x, y, z, list, lptr, lend, point_file_name ) &
        bind(C, name="tr_to_xyzl")

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
  !    Input, integer(ip) N, the number of nodes.
  !
  !    Input, real(dp) X(N), Y(N), Z(N), the coordinates of the nodes.
  !
  !    Input, integer(ip) LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N),
  !    information defining the triangulation, created by TRMESH.
  !
  !    Input, character ( len = * ) POINT_FILE_NAME, the name of the
  !    file containing the point coordinates.
  !

    integer(ip), intent(in), value :: n
    real(dp), intent(in) :: x(n)
    real(dp), intent(in) :: y(n)
    real(dp), intent(in) :: z(n)
    integer(ip), intent(in) :: list(6*(n-2))
    integer(ip), intent(in) :: lptr(6*(n-2))
    integer(ip), intent(in) :: lend(n)
    character ( len = * ), intent(in) :: point_file_name

    character ( len = 80 ) :: delaunay_file_name
    logical, parameter :: header = .false.
    integer(ip) :: iunit
    integer(ip) :: lp
    integer(ip) :: lpl
    integer(ip) :: node1
    integer(ip) :: node2
    character ( len = 40 ) :: string

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
  end subroutine tr_to_xyzl

  subroutine xyz_read ( point_file_name, n, x, y, z, ierror ) &
        bind(C, name="xyz_read")

  !*****************************************************************************80
  !
  !! XYZ_READ reads graphics information from an XYZ file.
  !
  !  Discussion:
  !
  !    Comment lines begin with '#";
  !    The XYZ coordinates of a point are written on a single line.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    05 June 2002
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) POINT_FILE_NAME, the name of the input file.
  !
  !    Input, integer(ip) N, the number of points.
  !
  !    Output, real(dp) X(N), Y(N), Z(N), the point coordinates.
  !
  !    Output, integer(ip) IERROR, error flag.
  !    0, no error occurred.
  !    nonzero, an error occurred.
  !

    character ( len = * ), intent(in) :: point_file_name
    integer(ip), intent(in), value :: n
    real(dp), intent(out) :: x(n)
    real(dp), intent(out) :: y(n)
    real(dp), intent(out) :: z(n)
    integer(ip), intent(out) :: ierror

    logical :: done
    integer(ip) :: i
    integer(ip) :: ios
    integer(ip) :: iunit
    integer(ip) :: lchar
    character ( len = 256 ) :: line
    integer(ip) :: n2
    real(dp) :: temp(3)
    integer(ip) :: text_num
    character ( len = 100 ) :: word

    n2 = 0
    ierror = 0
    word = ' '
    text_num = 0

    call get_unit ( iunit )

    open ( unit = iunit, file = point_file_name, status = 'old', iostat = ios )

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'XYZ_READ - Fatal error!'
      write ( *, '(a)' ) '  Could not open the input file.'
      ierror = 1
      stop
    end if
  !
  !  Read a line of text from the file.
  !
    do

      read ( iunit, '(a)', iostat = ios ) line

      if ( ios /= 0 ) then
        exit
      end if

      text_num = text_num + 1
  !
  !  If this line begins with '#' , then it's a comment.  Read a new line.
  !
      if ( line(1:1) == '#' ) then
        cycle
      end if
  !
  !  If this line is blank, then record that information.
  !
      if ( len_trim ( line ) == 0 ) then
        cycle
      end if
  !
  !  This line records a node's coordinates.
  !
      n2 = n2 + 1

      if ( n2 <= n ) then

        done = .true.

        do i = 1, 3

          call word_next_read ( line, word, done )

          call s_to_r8 ( word, temp(i), ierror, lchar )

          if ( ierror /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'XYZ_READ - Fatal error!'
            write ( *, '(a,i8)' ) '  S_TO_R8 returned IERROR = ', ierror
            write ( *, '(a,i8)' ) '  Reading (X,Y,Z) component ', i
            exit
          end if

        end do

        if ( ierror /= 0 ) then
          exit
        end if

        x(n2) = temp(1)
        y(n2) = temp(2)
        z(n2) = temp(3)

      end if

    end do

    close ( unit = iunit )
  !
  !  Report.
  !
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8,a)' ) 'XYZ_READ:'
    write ( *, '(a,i8,a)' ) '  Read ', text_num, ' text lines from ' &
      // trim ( point_file_name )
    write ( *, '(a,i8,a)' ) '  Read ', n2, ' sets of (X,Y,Z) coordinates.'
  end subroutine xyz_read

end module sphere_delaunay_mod
