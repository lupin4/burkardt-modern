!> sphere_voronoi — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module sphere_voronoi_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  implicit none
  private

  public :: stripack_interface, poly_count, poly_count_max, vr_to_xyzf, xyz_read

contains

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

    real(real64) a
    real(real64), allocatable, dimension ( : ) :: ds
    real(real64) elat
    real(real64) elon
    integer(int32) i
    integer(int32) ierror
    integer(int32) iunit
    integer(int32), allocatable, dimension ( : ) :: iwk
    integer(int32) k
    integer(int32) kt
    integer(int32), allocatable, dimension ( :, : ) :: lbtri
    integer(int32), allocatable, dimension ( : ) :: lend
    integer(int32), allocatable, dimension ( : ) :: list
    integer(int32), allocatable, dimension ( : ) :: listc
    integer(int32) lnew
    integer(int32) lp
    integer(int32) lpl
    integer(int32), allocatable, dimension ( : ) :: lptr
    integer(int32), allocatable, dimension ( :, : ) :: ltri
    integer(int32) n
    integer(int32) na
    integer(int32) nb
    integer(int32) nn
    real(real64) norm
    integer(int32) nt
    integer(int32) ntemp
    logical numbr
    integer(int32) nv
    real(real64), parameter :: pltsiz = 7.5e+00_real64
    character ( len = * ) point_file_name
    real(real64), allocatable, dimension ( : ) :: rc
    integer(int32) side_max
    real(real64) vlat
    real(real64) vlon
    character ( len = 255 ) voronoi_plot_file_name
    character ( len = 255 ) voronoi_plot_title
    real(real64), allocatable, dimension ( : ) :: x
    real(real64), allocatable, dimension ( : ) :: xc
    real(real64), allocatable, dimension ( : ) :: y
    real(real64), allocatable, dimension ( : ) :: yc
    real(real64), allocatable, dimension ( : ) :: z
    real(real64), allocatable, dimension ( : ) :: zc
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
    elat = 0.0e+00_real64
    elon = 0.0e+00_real64
    a = 90.0e+00_real64
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
  !    Input, integer(int32) N, the number of Voronoi polygons.
  !
  !    Input, integer(int32) LEND(N), some kind of pointer.
  !
  !    Input, integer(int32) LPTR(6*(N-2)), some other kind of pointer.
  !
  !    Input, integer(int32) LISTC(6*(N-2)), some other kind of pointer.
  !
  !    Input, integer(int32) SIDE_MAX, the maximum polygonal order.
  !

    integer(int32) n
    integer(int32) side_max

    integer(int32) count(0:side_max)
    integer(int32) edges
    integer(int32) i
    integer(int32) kv
    integer(int32) lend(n)
    integer(int32) listc(6*(n-2))
    integer(int32) lp
    integer(int32) lpl
    integer(int32) lptr(6*(n-2))
    integer(int32) n0
    integer(int32) sides
    integer(int32) vertices

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
  !    Input, integer(int32) N, the number of Voronoi polygons.
  !
  !    Input, integer(int32) LEND(N), some kind of pointer.
  !
  !    Input, integer(int32) LPTR(6*(N-2)), some other kind of pointer.
  !
  !    Input, integer(int32) LISTC(6*(N-2)), some other kind of pointer.
  !
  !    Output, integer(int32) SIDE_MAX, the highest order polgon.
  !

    integer(int32) n

    integer(int32) i
    integer(int32) lend(n)
    integer(int32) listc(6*(n-2))
    integer(int32) lp
    integer(int32) lpl
    integer(int32) lptr(6*(n-2))
    integer(int32) n0
    integer(int32) side_max
    integer(int32) sides

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
  !    Input, integer(int32) N, the number of nodes.
  !
  !    Input, real(real64) XC(2*(N-2)), YC(2*(N-2)), ZC(2*(N-2)),
  !    the coordinates  of the Voronoi vertices.
  !
  !    Input, integer(int32) LISTC(6*(N-2)), LPTR(6*(N-2)), LEND(N), 
  !    information defining the triangulation, created by TRMESH.
  !
  !    Input, character ( len = * ) POINT_FILE_NAME, the name of the
  !    file containing the point coordinates.
  !

    integer(int32) n

    logical, parameter :: header = .false.
    integer(int32) i
    integer(int32) iunit
    integer(int32) kv1
    integer(int32) kv2
    integer(int32) lend(n)
    integer(int32) listc(6*(n-2))
    integer(int32) lp
    integer(int32) lpl
    integer(int32) lptr(6*(n-2))
    integer(int32) node1
    integer(int32) node2
    character ( len = *  ) point_file_name
    character ( len = 40 ) string
    character ( len = 255 ) voronoi_face_file_name
    character ( len = 255 ) voronoi_vertex_file_name
    real(real64) xc(2*(n-2))
    real(real64) yc(2*(n-2))
    real(real64) zc(2*(n-2))
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

  subroutine xyz_read ( point_file_name, n, x, y, z, ierror )

  !*****************************************************************************80
  !
  !! XYZ_READ reads graphics information from an XYZ file.
  !
  !  Discussion:
  !
  !    Comment lines begin with '#";
  !    The XYZ coordinates of a point are written on a single line.
  !
  !  Example:
  !
  !     # cube.xyz
  !     #
  !     0 0 0
  !     0 0 1
  !     0 1 0
  !     0 1 1
  !     1 0 0
  !     1 0 1
  !     1 1 0
  !     1 1 1
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
  !    Input, character ( len = * ) XYZ_FILE_NAME, the name of the input file.
  !
  !    Input, integer(int32) N, the number of points.
  !
  !    Output, real(real64) X(N), Y(N), Z(N), the point coordinates.
  !
  !    Output, integer(int32) IERROR, error flag.
  !    0, no error occurred.
  !    nonzero, an error occurred.
  !

    integer(int32) n

    logical done
    integer(int32) i
    integer(int32) ierror
    integer(int32) ios
    integer(int32) iunit
    integer(int32) lchar
    character ( len = 255 ) line
    integer(int32) n2
    real(real64) temp(3)
    integer(int32) text_num
    character ( len = 255 ) word
    real(real64) x(n)
    character ( len = * ) point_file_name
    real(real64) y(n)
    real(real64) z(n)

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
  end

end module sphere_voronoi_mod
