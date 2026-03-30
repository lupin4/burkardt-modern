!> tetrahedron_properties — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

subroutine tetrahedron_centroid ( tetra, centroid )

!*****************************************************************************80
!
!! TETRAHEDRON_CENTROID computes the centroid of a tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision TETRA(3,4) the tetrahedron vertices.
!
!    Output, double precision CENTROID(3), the coordinates of the centroid.
!
  implicit none

  double precision centroid(3)
  integer i
  double precision tetra(3,4)

  do i = 1, 3
    centroid(i) = sum ( tetra(i,1:4) ) / 4.0D+00
  end do
end

subroutine tetrahedron_circumsphere ( tetra, r, pc )

!*****************************************************************************80
!
!! TETRAHEDRON_CIRCUMSPHERE computes the circumsphere of a tetrahedron.
!
!  Discussion:
!
!    The circumsphere, or circumscribed sphere, of a tetrahedron is the
!    sphere that passes through the four vertices.  The circumsphere is
!    not necessarily the smallest sphere that contains the tetrahedron.
!
!    Surprisingly, the diameter of the sphere can be found by solving
!    a 3 by 3 linear system.  This is because the vectors P2 - P1,
!    P3 - P1 and P4 - P1 are secants of the sphere, and each forms a
!    right triangle with the diameter through P1.  Hence, the dot product of
!    P2 - P1 with that diameter is equal to the square of the length
!    of P2 - P1, and similarly for P3 - P1 and P4 - P1.  This determines
!    the diameter vector originating at P1, and hence the radius and
!    center.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Adrian Bowyer, John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Parameters:
!
!    Input, double precision TETRA(3,4) the tetrahedron vertices.
!
!    Output, double precision R, PC(3), the center of the
!    circumscribed sphere, and its radius.  If the linear system is
!    singular, then R = -1, PC(1:3) = 0.
!
  implicit none

  double precision a(3,4)
  integer i
  integer info
  integer j
  double precision pc(3)
  double precision r
  double precision tetra(3,4)
!
!  Set up the linear system.
!
  a(1:3,1:3) = transpose ( tetra(1:3,2:4) )

  do j = 1, 3
    a(1:3,j) = a(1:3,j) - tetra(j,1)
  end do

  do i = 1, 3
    a(i,4) = sum ( a(i,1:3)**2 )
  end do
!
!  Solve the linear system.
!
  call r8mat_solve ( 3, 1, a, info )
!
!  If the system was singular, return a consolation prize.
!
  if ( info /= 0 ) then
    r = -1.0D+00
    pc(1:3) = 0.0D+00
  end if
!
!  Compute the radius and center.
!
  r = 0.5D+00 * sqrt ( sum ( a(1:3,4)**2 ) )

  pc(1:3) = tetra(1:3,1) + 0.5D+00 * a(1:3,4)
end

subroutine tetrahedron_dihedral_angles ( tetra, angle )

!*****************************************************************************80
!
!! TETRAHEDRON_DIHEDRAL_ANGLES computes dihedral angles of a tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision TETRA(3,4), the vertices of the tetrahedron.
!
!    Output, double precision ANGLE(6), the dihedral angles along the
!    axes AB, AC, AD, BC, BD and CD, respectively.
!
  implicit none

  double precision ab(3)
  double precision abc_normal(3)
  double precision abd_normal(3)
  double precision ac(3)
  double precision acd_normal(3)
  double precision ad(3)
  double precision angle(6)
  double precision bc(3)
  double precision bcd_normal(3)
  double precision bd(3)
  double precision cd(3)
  double precision , parameter :: r8_pi = 3.141592653589793D+00
  double precision tetra(3,4)

  call tetrahedron_edges ( tetra, ab, ac, ad, bc, bd, cd )

  call r8vec_cross_3d ( ac, ab, abc_normal )
  call r8vec_cross_3d ( ab, ad, abd_normal )
  call r8vec_cross_3d ( ad, ac, acd_normal )
  call r8vec_cross_3d ( bc, bd, bcd_normal )

  call r8vec_angle_3d ( abc_normal, abd_normal, angle(1) )
  call r8vec_angle_3d ( abc_normal, acd_normal, angle(2) )
  call r8vec_angle_3d ( abd_normal, acd_normal, angle(3) )
  call r8vec_angle_3d ( abc_normal, bcd_normal, angle(4) )
  call r8vec_angle_3d ( abd_normal, bcd_normal, angle(5) )
  call r8vec_angle_3d ( acd_normal, bcd_normal, angle(6) )

  angle(1:6) = r8_pi - angle(1:6)
end

subroutine tetrahedron_edges ( tetra, ab, ac, ad, bc, bd, cd )

!*****************************************************************************80
!
!! TETRAHEDRON_EDGES computes the edges of a tetrahedron.
!
!  Discussion:
!
!    The vertices are A, B, C, D.  The edge from A to B is denoted by AB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 May 2014
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, double precision TETRA(3,4), the vertices of the tetrahedron.
!
!    Output, double precision AB(3), AC(3), AD(3), BC(3), BD(3), CD(3), 
!    vectors that represent the edges of the tetrahedron.
!
  implicit none

  double precision ab(3)
  double precision ac(3)
  double precision ad(3)
  double precision bc(3)
  double precision bd(3)
  double precision cd(3)
  double precision tetra(3,4)

  ab(1:3) = tetra(1:3,2) - tetra(1:3,1)
  ac(1:3) = tetra(1:3,3) - tetra(1:3,1)
  ad(1:3) = tetra(1:3,4) - tetra(1:3,1)
  bc(1:3) = tetra(1:3,3) - tetra(1:3,2)
  bd(1:3) = tetra(1:3,4) - tetra(1:3,2)
  cd(1:3) = tetra(1:3,4) - tetra(1:3,3)
end

subroutine tetrahedron_edge_length ( tetra, edge_length )

!*****************************************************************************80
!
!! TETRAHEDRON_EDGE_LENGTH returns edge lengths of a tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision TETRA(3,4), the tetrahedron vertices.
!
!    Output, double precision EDGE_LENGTH(6), the length of the edges.
!
  implicit none

  double precision r8vec_length
  double precision edge_length(6)
  integer j1
  integer j2
  integer k
  double precision tetra(3,4)

  k = 0
  do j1 = 1, 3
    do j2 = j1 + 1, 4
      k = k + 1
      edge_length(k) = r8vec_length ( 3, tetra(1:3,j2) - tetra(1:3,j1) )
     end do
  end do
end

subroutine tetrahedron_face_angles ( tetra, angles )

!*****************************************************************************80
!
!! TETRAHEDRON_FACE_ANGLES returns the 12 face angles of a tetrahedron.
!
!  Discussion:
!
!    The tetrahedron has 4 triangular faces.  This routine computes the
!    3 planar angles associated with each face.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision TETRA(3,4) the tetrahedron vertices.
!
!    Output, double precision ANGLES(3,4), the face angles.
!
  implicit none

  double precision angles(3,4)
  double precision tri(3,3)
  double precision tetra(3,4)
!
!  Face 123
!
  tri(1:3,1:3) = tetra(1:3,1:3)
  call triangle_angles_3d ( tri, angles(1:3,1) )
!
!  Face 124
!
  tri(1:3,1:2) = tetra(1:3,1:2)
  tri(1:3,3) = tetra(1:3,4)
  call triangle_angles_3d ( tri, angles(1:3,2) )
!
!  Face 134
!
  tri(1:3,1) = tetra(1:3,1)
  tri(1:3,2:3) = tetra(1:3,3:4)
  call triangle_angles_3d ( tri, angles(1:3,3) )
!
!  Face 234
!
  tri(1:3,1:3) = tetra(1:3,2:4)
  call triangle_angles_3d ( tri, angles(1:3,4) )
end

subroutine tetrahedron_face_areas ( tetra, areas )

!*****************************************************************************80
!
!! TETRAHEDRON_FACE_AREAS returns the 4 face areas of a tetrahedron.
!
!  Discussion:
!
!    The tetrahedron has 4 triangular faces.  This routine computes the
!    area of each face.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision TETRA(3,4) the tetrahedron vertices.
!
!    Output, double precision AREAS(4), the face areas.
!
  implicit none

  double precision areas(4)
  double precision tri(3,3)
  double precision tetra(3,4)
!
!  Face 123
!
  tri(1:3,1:3) = tetra(1:3,1:3)
  call triangle_area_3d ( tri, areas(1) )
!
!  Face 124
!
  tri(1:3,1:2) = tetra(1:3,1:2)
  tri(1:3,3) = tetra(1:3,4)
  call triangle_area_3d ( tri, areas(2) )
!
!  Face 134
!
  tri(1:3,1) = tetra(1:3,1)
  tri(1:3,2:3) = tetra(1:3,3:4)
  call triangle_area_3d ( tri, areas(3) )
!
!  Face 234
!
  tri(1:3,1:3) = tetra(1:3,2:4)
  call triangle_area_3d ( tri, areas(4) )
end

subroutine tetrahedron_insphere ( tetra, r, pc )

!*****************************************************************************80
!
!! TETRAHEDRON_INSPHERE finds the insphere of a tetrahedron.
!
!  Discussion:
!
!    The insphere of a tetrahedron is the inscribed sphere, which touches
!    each face of the tetrahedron at a single point.
!
!    The points of contact are the centroids of the triangular faces
!    of the tetrahedron.  Therefore, the point of contact for a face
!    can be computed as the average of the vertices of that face.
!
!    The sphere can then be determined as the unique sphere through
!    the four given centroids.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Schneider, David Eberly,
!    Geometric Tools for Computer Graphics,
!    Elsevier, 2002,
!    ISBN: 1558605940,
!    LC: T385.G6974.
!
!  Parameters:
!
!    Input, double precision TETRA(3,4), the vertices of the tetrahedron.
!
!    Output, double precision R, PC(3), the radius and the center
!    of the sphere.
!
  implicit none

  double precision b(4,4)
  double precision r8mat_det_4d
  double precision r8vec_length
  double precision gamma
  double precision l123
  double precision l124
  double precision l134
  double precision l234
  double precision n123(1:3)
  double precision n124(1:3)
  double precision n134(1:3)
  double precision n234(1:3)
  double precision pc(1:3)
  double precision r
  double precision tetra(1:3,4)
  double precision v21(1:3)
  double precision v31(1:3)
  double precision v41(1:3)
  double precision v32(1:3)
  double precision v42(1:3)
  double precision v43(1:3)

  call tetrahedron_edges ( tetra, v21, v31, v41, v32, v42, v43 )

  call r8vec_cross_3d ( v21, v31, n123 )
  call r8vec_cross_3d ( v41, v21, n124 )
  call r8vec_cross_3d ( v31, v41, n134 )
  call r8vec_cross_3d ( v42, v32, n234 )

  l123 = r8vec_length ( 3, n123 )
  l124 = r8vec_length ( 3, n124 )
  l134 = r8vec_length ( 3, n134 )
  l234 = r8vec_length ( 3, n234 )

  pc(1:3) = ( l234 * tetra(1:3,1)   &
            + l134 * tetra(1:3,2)   &
            + l124 * tetra(1:3,3)   &
            + l123 * tetra(1:3,4) ) &
            / ( l234 + l134 + l124 + l123 )

  b(1:3,1:4) = tetra(1:3,1:4)
  b(4,1:4) = 1.0D+00

  gamma = abs ( r8mat_det_4d ( b ) )

  r = gamma / ( l234 + l134 + l124 + l123 )
end

subroutine tetrahedron_quality1 ( tetra, quality )

!*****************************************************************************80
!
!! TETRAHEDRON_QUALITY1: "quality" of a tetrahedron.
!
!  Discussion:
!
!    The quality of a tetrahedron is 3 times the ratio of the radius of
!    the inscribed sphere divided by that of the circumscribed sphere.
!
!    An equilateral tetrahredron achieves the maximum possible quality of 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision TETRA(3,4), the tetrahedron vertices.
!
!    Output, double precision QUALITY, the quality of the tetrahedron.
!
  implicit none

  double precision pc(3)
  double precision quality
  double precision r_in
  double precision r_out
  double precision tetra(3,4)

  call tetrahedron_circumsphere ( tetra, r_out, pc )

  call tetrahedron_insphere ( tetra, r_in, pc )

  quality = 3.0D+00 * r_in / r_out
end

subroutine tetrahedron_quality2 ( tetra, quality2 )

!*****************************************************************************80
!
!! TETRAHEDRON_QUALITY2: "quality" of a tetrahedron.
!
!  Discussion:
!
!    The quality measure #2 of a tetrahedron is:
!
!      QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
!
!    where
!
!      RIN = radius of the inscribed sphere;
!      LMAX = length of longest side of the tetrahedron.
!
!    An equilateral tetrahredron achieves the maximum possible quality of 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Qiang Du, Desheng Wang,
!    The Optimal Centroidal Voronoi Tesselations and the Gersho's
!    Conjecture in the Three-Dimensional Space,
!    Computers and Mathematics with Applications,
!    Volume 49, 2005, pages 1355-1373.
!
!  Parameters:
!
!    Input, double precision TETRA(3,4), the tetrahedron vertices.
!
!    Output, double precision QUALITY2, the quality of the tetrahedron.
!
  implicit none

  double precision edge_length(6)
  double precision l_max
  double precision pc(3)
  double precision quality2
  double precision r_in
  double precision tetra(3,4)

  call tetrahedron_edge_length ( tetra, edge_length )

  l_max = maxval ( edge_length(1:6) )

  call tetrahedron_insphere ( tetra, r_in, pc )

  quality2 = 2.0D+00 * sqrt ( 6.0D+00 ) * r_in / l_max
end

subroutine tetrahedron_quality3 ( tetra, quality3 )

!*****************************************************************************80
!
!! TETRAHEDRON_QUALITY3 computes the mean ratio of a tetrahedron.
!
!  Discussion:
!
!    This routine computes QUALITY3, the eigenvalue or mean ratio of
!    a tetrahedron.
!
!      QUALITY3 = 12 * ( 3 * volume )^(2/3) / (sum of squares of edge lengths).
!
!    This value may be used as a shape quality measure for the tetrahedron.
!
!    For an equilateral tetrahedron, the value of this quality measure
!    will be 1.  For any other tetrahedron, the value will be between
!    0 and 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2005
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
!    Input, double precision TETRA(3,4), the vertices of the tetrahedron.
!
!    Output, double precision QUALITY3, the mean ratio of the tetrahedron.
!
  implicit none

  double precision ab(3)
  double precision ac(3)
  double precision ad(3)
  double precision bc(3)
  double precision bd(3)
  double precision cd(3)
  double precision denom
  double precision lab
  double precision lac
  double precision lad
  double precision lbc
  double precision lbd
  double precision lcd
  double precision quality3
  double precision tetra(3,4)
  double precision volume
!
!  Compute the vectors representing the sides of the tetrahedron.
!
  call tetrahedron_edges ( tetra, ab, ac, ad, bc, bd, cd )
!
!  Compute the squares of the lengths of the sides.
!
  lab = sum ( ab(1:3)**2 )
  lac = sum ( ac(1:3)**2 )
  lad = sum ( ad(1:3)**2 )
  lbc = sum ( bc(1:3)**2 )
  lbd = sum ( bd(1:3)**2 )
  lcd = sum ( cd(1:3)**2 )
!
!  Compute the volume.
!
  volume = abs ( &
      ab(1) * ( ac(2) * ad(3) - ac(3) * ad(2) ) &
    + ab(2) * ( ac(3) * ad(1) - ac(1) * ad(3) ) &
    + ab(3) * ( ac(1) * ad(2) - ac(2) * ad(1) ) ) / 6.0D+00

  denom = lab + lac + lad + lbc + lbd + lcd

  if ( denom == 0.0D+00 ) then
    quality3 = 0.0D+00
  else
    quality3 = 12.0D+00 * ( 3.0D+00 * volume )**( 2.0D+00 / 3.0D+00 ) / denom
  end if
end

subroutine tetrahedron_quality4 ( tetra, quality4 )

!*****************************************************************************80
!
!! TETRAHEDRON_QUALITY4 computes the minimum solid angle of a tetrahedron.
!
!  Discussion:
!
!    This routine computes a quality measure for a tetrahedron, based
!    on the sine of half the minimum of the four solid angles.
!
!    The quality measure for an equilateral tetrahedron should be 1,
!    since the solid angles of such a tetrahedron are each equal to pi.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2005
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
!    Input, double precision TETRA(3,4), the vertices of the tetrahedron.
!
!    Output, double precision QUALITY4, the value of the quality measure.
!
  implicit none

  double precision ab(3)
  double precision ac(3)
  double precision ad(3)
  double precision bc(3)
  double precision bd(3)
  double precision cd(3)
  double precision denom
  double precision l1
  double precision l2
  double precision l3
  double precision lab
  double precision lac
  double precision lad
  double precision lbc
  double precision lbd
  double precision lcd
  double precision quality4
  double precision tetra(3,4)
  double precision volume
!
!  Compute the vectors that represent the sides.
!
  call tetrahedron_edges ( tetra, ab, ac, ad, bc, bd, cd )
!
!  Compute the lengths of the sides.
!
  lab = sqrt ( sum ( ab(1:3)**2 ) )
  lac = sqrt ( sum ( ac(1:3)**2 ) )
  lad = sqrt ( sum ( ad(1:3)**2 ) )
  lbc = sqrt ( sum ( bc(1:3)**2 ) )
  lbd = sqrt ( sum ( bd(1:3)**2 ) )
  lcd = sqrt ( sum ( cd(1:3)**2 ) )
!
!  Compute the volume
!
  volume = abs ( &
      ab(1) * ( ac(2) * ad(3) - ac(3) * ad(2) ) &
    + ab(2) * ( ac(3) * ad(1) - ac(1) * ad(3) ) &
    + ab(3) * ( ac(1) * ad(2) - ac(2) * ad(1) ) ) / 6.0D+00

  quality4 = 1.0D+00

  l1 = lab + lac
  l2 = lab + lad
  l3 = lac + lad

  denom = ( l1 + lbc ) * ( l1 - lbc ) &
        * ( l2 + lbd ) * ( l2 - lbd ) &
        * ( l3 + lcd ) * ( l3 - lcd )

  if ( denom <= 0.0D+00 ) then
    quality4 = 0.0D+00
  else
    quality4 = min ( quality4, 12.0D+00 * volume / sqrt ( denom ) )
  end if

  l1 = lab + lbc
  l2 = lab + lbd
  l3 = lbc + lbd

  denom = ( l1 + lac ) * ( l1 - lac ) &
        * ( l2 + lad ) * ( l2 - lad ) &
        * ( l3 + lcd ) * ( l3 - lcd )

  if ( denom <= 0.0D+00 ) then
    quality4 = 0.0D+00
  else
    quality4 = min ( quality4, 12.0D+00 * volume / sqrt ( denom ) )
  end if

  l1 = lac + lbc
  l2 = lac + lcd
  l3 = lbc + lcd

  denom = ( l1 + lab ) * ( l1 - lab ) &
        * ( l2 + lad ) * ( l2 - lad ) &
        * ( l3 + lbd ) * ( l3 - lbd )

  if ( denom <= 0.0D+00 ) then
    quality4 = 0.0D+00
  else
    quality4 = min ( quality4, 12.0D+00 * volume / sqrt ( denom ) )
  end if

  l1 = lad + lbd
  l2 = lad + lcd
  l3 = lbd + lcd

  denom = ( l1 + lab ) * ( l1 - lab ) &
        * ( l2 + lac ) * ( l2 - lac ) &
        * ( l3 + lbc ) * ( l3 - lbc )

  if ( denom <= 0.0D+00 ) then
    quality4 = 0.0D+00
  else
    quality4 = min ( quality4, 12.0D+00 * volume / sqrt ( denom ) )
  end if

  quality4 = quality4 * 1.5D+00 * sqrt ( 6.0D+00 )
end

subroutine tetrahedron_solid_angles ( tetra, angle )

!*****************************************************************************80
!
!! TETRAHEDRON_SOLID_ANGLES computes solid angles of a tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision TETRA(3,4), the vertices of the tetrahedron.
!
!    Output, double precision ANGLE(4), the solid angles.
!
  implicit none

  double precision angle(4)
  double precision dihedral_angles(6)
  double precision , parameter :: r8_pi = 3.141592653589793D+00
  double precision tetra(3,4)

  call tetrahedron_dihedral_angles ( tetra, dihedral_angles )

  angle(1) = dihedral_angles(1) &
           + dihedral_angles(2) &
           + dihedral_angles(3) - r8_pi

  angle(2) = dihedral_angles(1) &
           + dihedral_angles(4) &
           + dihedral_angles(5) - r8_pi

  angle(3) = dihedral_angles(2) &
           + dihedral_angles(4) &
           + dihedral_angles(6) - r8_pi

  angle(4) = dihedral_angles(3) &
           + dihedral_angles(5) &
           + dihedral_angles(6) - r8_pi
end

subroutine tetrahedron_volume ( tetra, volume )

!*****************************************************************************80
!
!! TETRAHEDRON_VOLUME computes the volume of a tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision TETRA(3,4), the vertices of the tetrahedron.
!
!    Output, double precision VOLUME, the volume of the tetrahedron.
!
  implicit none

  double precision a(4,4)
  double precision r8mat_det_4d
  double precision tetra(3,4)
  double precision volume

  a(1:3,1:4) = tetra(1:3,1:4)
  a(4,1:4) = 1.0D+00

  volume = abs ( r8mat_det_4d ( a ) ) / 6.0D+00
end

subroutine triangle_angles_3d ( t, angle )

!*****************************************************************************80
!
!! TRIANGLE_ANGLES_3D computes the angles of a triangle in 3D.
!
!  Discussion:
!
!    The law of cosines is used:
!
!      C * C = A * A + B * B - 2 * A * B * COS ( GAMMA )
!
!    where GAMMA is the angle opposite side C.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision T(3,3), the triangle vertices.
!
!    Output, double precision ANGLE(3), the angles opposite
!    sides P1-P2, P2-P3 and P3-P1, in radians.
!
  implicit none

  double precision a
  double precision angle(3)
  double precision b
  double precision c
  double precision r8_acos
  double precision , parameter :: r8_pi = 3.141592653589793D+00
  double precision t(3,3)
!
!  Compute the length of each side.
!
  a = sqrt ( sum ( ( t(1:3,1) - t(1:3,2) )**2 ) )
  b = sqrt ( sum ( ( t(1:3,2) - t(1:3,3) )**2 ) )
  c = sqrt ( sum ( ( t(1:3,3) - t(1:3,1) )**2 ) )
!
!  Take care of a ridiculous special case.
!
  if ( a == 0.0D+00 .and. b == 0.0D+00 .and. c == 0.0D+00 ) then
    angle(1:3) = 2.0D+00 * r8_pi / 3.0D+00
  end if

  if ( c == 0.0D+00 .or. a == 0.0D+00 ) then
    angle(1) = r8_pi
  else
    angle(1) = r8_acos ( ( c * c + a * a - b * b ) / ( 2.0D+00 * c * a ) )
  end if

  if ( a == 0.0D+00 .or. b == 0.0D+00 ) then
    angle(2) = r8_pi
  else
    angle(2) = r8_acos ( ( a * a + b * b - c * c ) / ( 2.0D+00 * a * b ) )
  end if

  if ( b == 0.0D+00 .or. c == 0.0D+00 ) then
    angle(3) = r8_pi
  else
    angle(3) = r8_acos ( ( b * b + c * c - a * a ) / ( 2.0D+00 * b * c ) )
  end if
end

subroutine triangle_area_3d ( t, area )

!*****************************************************************************80
!
!! TRIANGLE_AREA_3D computes the area of a triangle in 3D.
!
!  Discussion:
!
!    This routine uses the fact that the norm of the cross product
!    of two vectors is the area of the parallelogram they form.
!
!    Therefore, the area of the triangle is half of that value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Adrian Bowyer, John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Parameters:
!
!    Input, double precision T(3,3), the triangle vertices.
!
!    Output, double precision AREA, the area of the triangle.
!
  implicit none

  double precision area
  double precision cross(3)
  double precision t(3,3)
!
!  Compute the cross product vector.
!
  cross(1) = ( t(2,2) - t(2,1) ) * ( t(3,3) - t(3,1) ) &
           - ( t(3,2) - t(3,1) ) * ( t(2,3) - t(2,1) )

  cross(2) = ( t(3,2) - t(3,1) ) * ( t(1,3) - t(1,1) ) &
           - ( t(1,2) - t(1,1) ) * ( t(3,3) - t(3,1) )

  cross(3) = ( t(1,2) - t(1,1) ) * ( t(2,3) - t(2,1) ) &
           - ( t(2,2) - t(2,1) ) * ( t(1,3) - t(1,1) )

  area = 0.5D+00 * sqrt ( sum ( cross(1:3)**2 ) )
end
