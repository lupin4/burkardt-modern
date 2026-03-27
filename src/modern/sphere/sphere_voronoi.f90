!> sphere_voronoi -- Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module sphere_voronoi_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: stripack_interface, poly_count, poly_count_max, vr_to_xyzf, xyz_read

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
    real(dp) :: vlat
    real(dp) :: vlon
    character ( len = 255 ) :: voronoi_plot_file_name
    character ( len = 255 ) :: voronoi_plot_title
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
    elat = 0.0_dp
    elon = 0.0_dp
    a = 90.0_dp
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
  end subroutine stripack_interface

  subroutine poly_count ( n, lend, lptr, listc, side_max ) &
        bind(C, name="poly_count")

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
  !    Input, integer(ip) N, the number of Voronoi polygons.
  !
  !    Input, integer(ip) LEND(N), some kind of pointer.
  !
  !    Input, integer(ip) LPTR(6*(N-2)), some other kind of pointer.
  !
  !    Input, integer(ip) LISTC(6*(N-2)), some other kind of pointer.
  !
  !    Input, integer(ip) SIDE_MAX, the maximum polygonal order.
  !

    integer(ip), intent(in), value :: n
    integer(ip), intent(in), value :: side_max
    integer(ip), intent(in) :: lend(n)
    integer(ip), intent(in) :: listc(6*(n-2))
    integer(ip), intent(in) :: lptr(6*(n-2))

    integer(ip) :: count(0:side_max)
    integer(ip) :: edges
    integer(ip) :: i
    integer(ip) :: kv
    integer(ip) :: lp
    integer(ip) :: lpl
    integer(ip) :: n0
    integer(ip) :: sides
    integer(ip) :: vertices

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
  end subroutine poly_count

  subroutine poly_count_max ( n, lend, lptr, listc, side_max ) &
        bind(C, name="poly_count_max")

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
  !    Input, integer(ip) N, the number of Voronoi polygons.
  !
  !    Input, integer(ip) LEND(N), some kind of pointer.
  !
  !    Input, integer(ip) LPTR(6*(N-2)), some other kind of pointer.
  !
  !    Input, integer(ip) LISTC(6*(N-2)), some other kind of pointer.
  !
  !    Output, integer(ip) SIDE_MAX, the highest order polygon.
  !

    integer(ip), intent(in), value :: n
    integer(ip), intent(in) :: lend(n)
    integer(ip), intent(in) :: listc(6*(n-2))
    integer(ip), intent(in) :: lptr(6*(n-2))
    integer(ip), intent(out) :: side_max

    integer(ip) :: i
    integer(ip) :: lp
    integer(ip) :: lpl
    integer(ip) :: n0
    integer(ip) :: sides

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
  end subroutine poly_count_max

  subroutine vr_to_xyzf ( n, xc, yc, zc, listc, lptr, lend, point_file_name ) &
        bind(C, name="vr_to_xyzf")

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
  !    Input, integer(ip) N, the number of nodes.
  !
  !    Input, real(dp) XC(2*(N-2)), YC(2*(N-2)), ZC(2*(N-2)),
  !    the coordinates  of the Voronoi vertices.
  !
  !    Input, integer(ip) LISTC(6*(N-2)), LPTR(6*(N-2)), LEND(N),
  !    information defining the triangulation, created by TRMESH.
  !
  !    Input, character ( len = * ) POINT_FILE_NAME, the name of the
  !    file containing the point coordinates.
  !

    integer(ip), intent(in), value :: n
    real(dp), intent(in) :: xc(2*(n-2))
    real(dp), intent(in) :: yc(2*(n-2))
    real(dp), intent(in) :: zc(2*(n-2))
    integer(ip), intent(in) :: listc(6*(n-2))
    integer(ip), intent(in) :: lptr(6*(n-2))
    integer(ip), intent(in) :: lend(n)
    character ( len = * ), intent(in) :: point_file_name

    logical, parameter :: header = .false.
    integer(ip) :: i
    integer(ip) :: iunit
    integer(ip) :: kv1
    integer(ip) :: kv2
    integer(ip) :: lp
    integer(ip) :: lpl
    integer(ip) :: node1
    integer(ip) :: node2
    character ( len = 40 ) :: string
    character ( len = 255 ) :: voronoi_face_file_name
    character ( len = 255 ) :: voronoi_vertex_file_name
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
  end subroutine vr_to_xyzf

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
    character ( len = 255 ) :: line
    integer(ip) :: n2
    real(dp) :: temp(3)
    integer(ip) :: text_num
    character ( len = 255 ) :: word

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

end module sphere_voronoi_mod
