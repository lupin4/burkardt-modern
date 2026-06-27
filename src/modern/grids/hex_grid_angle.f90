!> hex_grid_angle — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

function box_01_contains_point_2d ( p )

!*****************************************************************************80
!
!! BOX_01_CONTAINS_POINT_2D determines if a point is inside the unit box in 2D.
!
!  Discussion:
!
!    A unit box in 2D is the set of points (X,Y) with the property that
!
!      0.0 <= X <= 1.0
!    and
!      0.0 <= Y <= 1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision P(2), the point to be checked.
!
!    Output, logical BOX_01_CONTAINS_POINT_2D, is TRUE if the point is
!    inside the box.
!
  implicit none

  integer , parameter :: dim_num = 2

  logical box_01_contains_point_2d
  double precision p(dim_num)

  box_01_contains_point_2d = &
  ( &
    all ( 0.0D+00 <= p(1:dim_num) ) &
  .and. &
    all (            p(1:dim_num) <= 1.0D+00 ) &
  )
end

function box_contains_point_2d ( box, p )

!*****************************************************************************80
!
!! BOX_CONTAINS_POINT_2D determines if a point is inside a box in 2D.
!
!  Discussion:
!
!    A unit box in 2D is the set of points (X,Y) with the property that
!
!      0.0 <= X <= 1.0
!    and
!      0.0 <= Y <= 1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision BOX(2,2), the lower left and upper right
!    corners of the box.
!
!    Input, double precision P(2), the point to be checked.
!
!    Output, logical BOX_CONTAINS_POINT_2D, is TRUE if the point is
!    inside the box.
!
  implicit none

  integer , parameter :: dim_num = 2

  double precision box(dim_num,2)
  logical box_contains_point_2d
  double precision p(dim_num)

  box_contains_point_2d = &
    all ( box(1:dim_num,1) <= p(1:dim_num) ) .and. &
    all ( p(1:dim_num) <= box(1:dim_num,2) )
end

function cos_deg ( angle )

!*****************************************************************************80
!
!! COS_DEG returns the cosine of an angle given in degrees.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision ANGLE, the angle, in degrees.
!
!    Output, double precision COS_DEG, the cosine of the angle.
!
  implicit none

  double precision angle
  double precision cos_deg
  double precision , parameter :: degrees_to_radians &
    = 3.141592653589793D+00 / 180.0D+00

  cos_deg = cos ( degrees_to_radians * angle )
end

function degrees_to_radians ( degrees )

!*****************************************************************************80
!
!! DEGREES_TO_RADIANS converts an angle measure from degrees to radians.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision DEGREES, the angle measure in degrees.
!
!    Output, double precision DEGREES_TO_RADIANS, the angle measure in radians.
!
  implicit none

  double precision degrees
  double precision degrees_to_radians

  degrees_to_radians = ( degrees / 180.0D+00 ) * 3.141592653589793D+00
end

subroutine hex_grid_angle_01 ( center, angle, h, n, r )

!*****************************************************************************80
!
!! HEX_GRID_ANGLE_01 sets the points in an angled hex grid in the unit box.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision CENTER(2), the center of the grid.
!    This point must be inside the unit square.
!
!    Input, double precision ANGLE, the angle, in degrees, of the grid.
!    Normally, 0 <= ANGLE <= 180, but any value is allowed.
!
!    Input, double precision H, the spacing between neighboring
!    points on a grid line.
!
!    Input, integer N, the number of points of the angled hex grid
!    that are within the unit square.  This value may have been computed
!    by calling HEX_GRID_ANGLE_01_SIZE.
!
!    Output, double precision R(2,N), the grid points.
!
  implicit none

  integer n
  integer , parameter :: dim_num = 2

  double precision angle
  double precision angle2
  logical box_01_contains_point_2d
  double precision center(dim_num)
  double precision cos_deg
  double precision h
  integer i
  integer j
  integer k
  integer layer
  integer layer_size
  double precision point(dim_num)
  double precision r(dim_num,n)
  double precision r8_modp
  double precision sin_deg
  integer size
!
!  Ninny checks.
!
  if ( .not. box_01_contains_point_2d ( center ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEX_GRID_ANGLE_01 - Fatal error!'
    write ( *, '(a)' ) '  The center point of the grid is not'
    write ( *, '(a)' ) '  inside the unit square.'
    write ( *, '(a,2g14.6)' ) '  CENTER = ', center(1:dim_num)
    stop
  end if

  if ( h == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEX_GRID_ANGLE_01 - Fatal error!'
    write ( *, '(a)' ) '  The grid spacing must be nonzero.'
    write ( *, '(a,g14.6)' ) '  H = ', h
    stop
  end if

  layer = 0
  point(1:dim_num) = center(1:dim_num)

  k = 1
  if ( k <= n ) then
    r(1:dim_num,k) = center(1:dim_num)
  end if

  do

    layer = layer + 1

    layer_size = 0

    angle2 = angle
!
!  Compute the first point on the new layer.
!
    point(1:dim_num) = point(1:dim_num) &
      + h * (/ cos_deg ( angle2 ), sin_deg ( angle2 ) /)

    angle2 = r8_modp ( angle2 + 60.0D+00, 360.0D+00 )

    do i = 1, 6

      angle2 = r8_modp ( angle2 + 60.0D+00, 360.0D+00 )

      do j = 1, layer

        point(1:dim_num) = point(1:dim_num) &
          + h * (/ cos_deg ( angle2 ), sin_deg ( angle2 ) /)

        if ( box_01_contains_point_2d ( point ) ) then

          layer_size = layer_size + 1
          k = k + 1

          if ( k <= n ) then
            r(1:dim_num,k) = point(1:dim_num)
          end if

        end if

      end do

    end do

    if ( layer_size == 0 ) then
      exit
    end if

  end do
end

subroutine hex_grid_angle_01_size ( center, angle, h, n )

!*****************************************************************************80
!
!! HEX_GRID_ANGLE_01_SIZE counts points in an angled hex grid in the unit box.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision CENTER(2), the center of the grid.
!    This point must be inside the unit square.
!
!    Input, double precision ANGLE, the angle, in degrees, of the grid.
!    Normally, 0 <= ANGLE <= 180, but any value is allowed.
!
!    Input, double precision H, the spacing between neighboring
!    points on a grid line.
!
!    Output, integer N, the number of points of the angled hex grid
!    that are within the unit square.
!
  implicit none

  integer , parameter :: dim_num = 2

  double precision angle
  double precision angle2
  logical box_01_contains_point_2d
  double precision center(dim_num)
  double precision cos_deg
  double precision h
  integer i
  integer j
  integer layer
  integer layer_size
  integer n
  double precision point(dim_num)
  double precision r8_modp
  double precision sin_deg
!
!  Ninny checks.
!
  if ( .not. box_01_contains_point_2d ( center ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEX_GRID_ANGLE_01_SIZE - Fatal error!'
    write ( *, '(a)' ) '  The center point of the grid is not'
    write ( *, '(a)' ) '  inside the unit square.'
    write ( *, '(a,2g14.6)' ) '  CENTER = ', center(1:dim_num)
    stop
  end if

  if ( h == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEX_GRID_ANGLE_01_SIZE - Fatal error!'
    write ( *, '(a)' ) '  The grid spacing must be nonzero.'
    write ( *, '(a,g14.6)' ) '  H = ', h
    stop
  end if

  n = 0

  layer = 0
  point(1:dim_num) = center(1:dim_num)

  n = 1

  do

    layer = layer + 1

    layer_size = 0

    angle2 = angle
!
!  Compute the first point on the new layer.
!
    point(1:dim_num) = point(1:dim_num) &
      + h * (/ cos_deg ( angle2 ), sin_deg ( angle2 ) /)

    angle2 = r8_modp ( angle2 + 60.0D+00, 360.0D+00 )

    do i = 1, 6

      angle2 = r8_modp ( angle2 + 60.0D+00, 360.0D+00 )

      do j = 1, layer

        point(1:dim_num) = point(1:dim_num) &
          + h * (/ cos_deg ( angle2 ), sin_deg ( angle2 ) /)

        if ( box_01_contains_point_2d ( point ) ) then
          layer_size = layer_size + 1
          n = n + 1
        end if

      end do

    end do

    if ( layer_size == 0 ) then
      exit
    end if

  end do
end

subroutine hex_grid_angle_01_write ( center, angle, h, n, r, file_out_name )

!*****************************************************************************80
!
!! HEX_GRID_ANGLE_01_WRITE writes an angled hex grid dataset to a file.
!
!  Discussion:
!
!    The initial lines of the file are comments, which begin with a
!    '#' character.
!
!    Thereafter, each line of the file contains the 2-dimensional
!    components of the next entry of the dataset.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision CENTER(2), the "center" of the grid.
!
!    Input, double precision ANGLE, the angle of the grid.
!
!    Input, double precision H, the spacing between points on a grid line.
!
!    Input, integer N, the number of elements in the subsequence.
!
!    Input, double precision R(2,N), the points.
!
!    Input, character ( len = * ) FILE_OUT_NAME, the output file name.
!
  implicit none

  integer , parameter :: dim_num = 2
  integer n

  double precision angle
  double precision center(dim_num)
  character ( len = * ) file_out_name
  double precision h
  integer file_out_unit
  integer ios
  integer j
  double precision r(dim_num,n)
  character ( len = 40 ) string

  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEX_GRID_ANGLE_01_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file:'
    write ( *, '(a)' ) '  "' // trim ( file_out_name ) // '".'
    stop
  end if

  write ( file_out_unit, '(a)'        ) '#  ' // trim ( file_out_name )
  write ( file_out_unit, '(a)'        ) &
    '#  created by HEX_GRID_ANGLE_01_WRITE.F90'
  write ( file_out_unit, '(a)'        ) '#'
  write ( file_out_unit, '(a)'        ) '#'
  write ( file_out_unit, '(a,i12)'    ) '#  DIM_NUM = ', dim_num
  write ( file_out_unit, '(a,i12)'    ) '#  N =       ', n
  write ( file_out_unit, '(a,2g14.6)' ) '#  CENTER =  ', center(1:dim_num)
  write ( file_out_unit, '(a,g14.6)'  ) '#  ANGLE =   ', angle
  write ( file_out_unit, '(a,g14.6)'  ) '#  H =       ', h
  write ( file_out_unit, '(a,g14.6)'  ) '#  EPSILON = ', epsilon ( r(1,1) )
  write ( file_out_unit, '(a)'        ) '#'

  do j = 1, n
    write ( file_out_unit, '(2x,f10.6,2x,f10.6)' ) r(1:dim_num,j)
  end do

  close ( unit = file_out_unit )
end

subroutine hex_grid_angle ( box, center, angle, h, n, r )

!*****************************************************************************80
!
!! HEX_GRID_ANGLE sets the points in an angled hex grid in a box.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision BOX(2,2), the lower left and upper right
!    corners of the box.
!
!    Input, double precision CENTER(2), the center of the grid.
!    This point must be inside the unit square.
!
!    Input, double precision ANGLE, the angle, in degrees, of the grid.
!    Normally, 0 <= ANGLE <= 180, but any value is allowed.
!
!    Input, double precision H, the spacing between neighboring
!    points on a grid line.
!
!    Input, integer N, the number of points of the angled hex grid
!    that are within the unit square.  This value may have been computed
!    by calling HEX_GRID_ANGLE_01_SIZE.
!
!    Output, double precision R(2,N), the grid points.
!
  implicit none

  integer n
  integer , parameter :: dim_num = 2

  double precision angle
  double precision angle2
  double precision box(dim_num,2)
  logical box_contains_point_2d
  double precision center(dim_num)
  double precision cos_deg
  double precision h
  integer i
  integer j
  integer k
  integer layer
  integer layer_size
  double precision point(dim_num)
  double precision r(dim_num,n)
  double precision r8_modp
  double precision sin_deg
  integer size
!
!  Ninny checks.
!
  if ( .not. box_contains_point_2d ( box, center ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEX_GRID_ANGLE - Fatal error!'
    write ( *, '(a)' ) '  The center point of the grid is not'
    write ( *, '(a)' ) '  inside the box.'
    write ( *, '(a,2g14.6)' ) '  CENTER = ', center(1:dim_num)
    stop
  end if

  if ( h == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEX_GRID_ANGLE - Fatal error!'
    write ( *, '(a)' ) '  The grid spacing must be nonzero.'
    write ( *, '(a,g14.6)' ) '  H = ', h
    stop
  end if

  layer = 0
  point(1:dim_num) = center(1:dim_num)

  k = 1
  if ( k <= n ) then
    r(1:dim_num,k) = center(1:dim_num)
  end if

  do

    layer = layer + 1

    layer_size = 0

    angle2 = angle
!
!  Compute the first point on the new layer.
!
    point(1:dim_num) = point(1:dim_num) &
      + h * (/ cos_deg ( angle2 ), sin_deg ( angle2 ) /)

    angle2 = r8_modp ( angle2 + 60.0D+00, 360.0D+00 )

    do i = 1, 6

      angle2 = r8_modp ( angle2 + 60.0D+00, 360.0D+00 )

      do j = 1, layer

        point(1:dim_num) = point(1:dim_num) &
          + h * (/ cos_deg ( angle2 ), sin_deg ( angle2 ) /)

        if ( box_contains_point_2d ( box, point ) ) then

          layer_size = layer_size + 1
          k = k + 1

          if ( k <= n ) then
            r(1:dim_num,k) = point(1:dim_num)
          end if

        end if

      end do

    end do

    if ( layer_size == 0 ) then
      exit
    end if

  end do
end

subroutine hex_grid_angle_size ( box, center, angle, h, n )

!*****************************************************************************80
!
!! HEX_GRID_ANGLE_SIZE counts the points in an angled hex grid in a box.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision BOX(2,2), the lower left and upper right
!    corners of the box.
!
!    Input, double precision CENTER(2), the center of the grid.
!    This point must be inside the box
!
!    Input, double precision ANGLE, the angle, in degrees, of the grid.
!    Normally, 0 <= ANGLE <= 180, but any value is allowed.
!
!    Input, double precision H, the spacing between neighboring
!    points on a grid line.
!
!    Output, integer N, the number of points of the angled hex grid
!    that are within the unit square.
!
  implicit none

  integer , parameter :: dim_num = 2

  double precision angle
  double precision angle2
  double precision box(dim_num,2)
  logical box_contains_point_2d
  double precision center(dim_num)
  double precision cos_deg
  double precision h
  integer i
  integer j
  integer layer
  integer layer_size
  integer n
  double precision point(dim_num)
  double precision r8_modp
  double precision sin_deg
!
!  Ninny checks.
!
  if ( .not. box_contains_point_2d ( box, center ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEX_GRID_ANGLE_SIZE - Fatal error!'
    write ( *, '(a)' ) '  The center point of the grid is not'
    write ( *, '(a)' ) '  inside the box.'
    write ( *, '(a,2g14.6)' ) '  CENTER = ', center(1:dim_num)
    stop
  end if

  if ( h == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEX_GRID_ANGLE_SIZE - Fatal error!'
    write ( *, '(a)' ) '  The grid spacing must be nonzero.'
    write ( *, '(a,g14.6)' ) '  H = ', h
    stop
  end if

  n = 0

  layer = 0
  point(1:dim_num) = center(1:dim_num)

  n = 1

  do

    layer = layer + 1

    layer_size = 0

    angle2 = angle
!
!  Compute the first point on the new layer.
!
    point(1:dim_num) = point(1:dim_num) &
      + h * (/ cos_deg ( angle2 ), sin_deg ( angle2 ) /)

    angle2 = r8_modp ( angle2 + 60.0D+00, 360.0D+00 )

    do i = 1, 6

      angle2 = r8_modp ( angle2 + 60.0D+00, 360.0D+00 )

      do j = 1, layer

        point(1:dim_num) = point(1:dim_num) &
          + h * (/ cos_deg ( angle2 ), sin_deg ( angle2 ) /)

        if ( box_contains_point_2d ( box, point ) ) then
          layer_size = layer_size + 1
          n = n + 1
        end if

      end do

    end do

    if ( layer_size == 0 ) then
      exit
    end if

  end do
end

subroutine hex_grid_angle_write ( box, center, angle, h, n, r, file_out_name )

!*****************************************************************************80
!
!! HEX_GRID_ANGLE_WRITE writes an angled hex grid dataset in a box to a file.
!
!  Discussion:
!
!    The initial lines of the file are comments, which begin with a
!    '#' character.
!
!    Thereafter, each line of the file contains the 2-dimensional
!    components of the next entry of the dataset.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision BOX(2,2), the lower left and upper right
!    corners of the box.
!
!    Input, double precision CENTER(2), the "center" of the grid.
!
!    Input, double precision ANGLE, the angle of the grid.
!
!    Input, double precision H, the spacing between points on a grid line.
!
!    Input, integer N, the number of elements in the subsequence.
!
!    Input, double precision R(2,N), the points.
!
!    Input, character ( len = * ) FILE_OUT_NAME, the output file name.
!
  implicit none

  integer , parameter :: dim_num = 2
  integer n

  double precision angle
  double precision box(dim_num,2)
  double precision center(dim_num)
  character ( len = * ) file_out_name
  double precision h
  integer file_out_unit
  integer ios
  integer j
  double precision r(dim_num,n)
  character ( len = 40 ) string

  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEX_GRID_ANGLE_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file:'
    write ( *, '(a)' ) '  "' // trim ( file_out_name ) // '".'
    stop
  end if

  write ( file_out_unit, '(a)'        ) '#  ' // trim ( file_out_name )
  write ( file_out_unit, '(a)'        ) &
    '#  created by HEX_GRID_ANGLE_WRITE.F90'
  write ( file_out_unit, '(a)'        ) '#'
  write ( file_out_unit, '(a)'        ) '#'
  write ( file_out_unit, '(a,i12)'    ) '#  DIM_NUM = ', dim_num
  write ( file_out_unit, '(a,i12)'    ) '#  N =       ', n
  write ( file_out_unit, '(a,2g14.6)' ) '#  BOX_LO =  ', box(1:dim_num,1)
  write ( file_out_unit, '(a,2g14.6)' ) '#  BOX_HI =  ', box(1:dim_num,2)
  write ( file_out_unit, '(a,2g14.6)' ) '#  CENTER =  ', center(1:dim_num)
  write ( file_out_unit, '(a,g14.6)'  ) '#  ANGLE =   ', angle
  write ( file_out_unit, '(a,g14.6)'  ) '#  H =       ', h
  write ( file_out_unit, '(a,g14.6)'  ) '#  EPSILON = ', epsilon ( r(1,1) )
  write ( file_out_unit, '(a)'        ) '#'

  do j = 1, n
    write ( file_out_unit, '(2x,f10.6,2x,f10.6)' ) r(1:dim_num,j)
  end do

  close ( unit = file_out_unit )
end

function r8_modp ( x, y )

!*****************************************************************************80
!
!! R8_MODP returns the nonnegative remainder of R8 division.
!
!  Discussion:
!
!    If
!      REM = R8_MODP ( X, Y )
!      RMULT = ( X - REM ) / Y
!    then
!      X = Y * RMULT + REM
!    where REM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360.0) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, R8_MODP(A,360.0) is between 0 and 360, always.
!
!  Example:
!
!        I         J     MOD R8_MODP  R8_MODP Factorization
!
!      107        50       7       7    107 =  2 *  50 + 7
!      107       -50       7       7    107 = -2 * -50 + 7
!     -107        50      -7      43   -107 = -3 *  50 + 43
!     -107       -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision X, the number to be divided.
!
!    Input, double precision Y, the number that divides X.
!
!    Output, double precision R8_MODP, the nonnegative remainder
!    when X is divided by Y.
!
  implicit none

  double precision r8_modp
  double precision x
  double precision y

  if ( y == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_MODP - Fatal error!'
    write ( *, '(a,g14.6)' ) '  R8_MODP ( X, Y ) called with Y = ', y
    stop
  end if

  r8_modp = mod ( x, y )

  if ( r8_modp < 0.0D+00 ) then
    r8_modp = r8_modp + abs ( y )
  end if
end

function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a double precision value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer k
  double precision r8_uniform_01
  integer seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed) * 4.656612875D-10
end

subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of double precision values.
!
!    For now, the input quantity SEED is an integer variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, double precision R(N), the vector of pseudorandom values.
!
  implicit none

  integer n

  integer i
  integer k
  integer seed
  double precision r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed) * 4.656612875D-10

  end do
end

function sin_deg ( angle )

!*****************************************************************************80
!
!! SIN_DEG returns the sine of an angle given in degrees.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision ANGLE, the angle, in degrees.
!
!    Output, double precision SIN_DEG, the sine of the angle.
!
  implicit none

  double precision angle
  double precision , parameter :: degrees_to_radians &
    = 3.141592653589793D+00 / 180.0D+00
  double precision sin_deg

  sin_deg = sin ( degrees_to_radians * angle )
end
