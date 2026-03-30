!> tet_mesh_quality — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

subroutine tet_mesh_node_order ( tetra_order, tetra_num, tetra_node, &
  node_num, node_order )

!*****************************************************************************80
!
!! TET_MESH_NODE_ORDER: determine the order of nodes in a tetra mesh.
!
!  Discussion:
!
!    The order of a node is the number of tetrahedrons that use that node
!    as a vertex.
!
!    Tetrahedrons of order 4 or 10 are allowed as input.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TETRA_ORDER, the order of the mesh, either 4 or 10.
!
!    Input, integer TETRA_NUM, the number of tetrahedrons.
!
!    Input, integer TETRA_NODE(TETRA_ORDER,TETRA_NUM), the nodes
!    that make up the tetrahedrons.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Output, integer NODE_ORDER(NODE_NUM), the order of each node.
!
  implicit none

  integer node_num
  integer tetra_num
  integer tetra_order

  integer i
  integer node
  integer node_order(node_num)
  integer tetra
  integer tetra_node(tetra_order,tetra_num)

  node_order(1:node_num) = 0

  do tetra = 1, tetra_num
    do i = 1, tetra_order
      node = tetra_node(i,tetra)
      if ( node < 1 .or. node_num < node ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TET_MESH_NODE_ORDER - Fatal error!'
        write ( *, '(a)' ) '  Illegal entry in TETRA_NODE.'
        stop
      else
        node_order(node) = node_order(node) + 1
      end if
    end do
  end do
end

subroutine tet_mesh_quality1 ( node_num, node_xyz, tetra_order, tetra_num, &
  tetra_node, tetra_quality )

!*****************************************************************************80
!
!! TET_MESH_QUALITY1 returns the quality of each tet in a mesh.
!
!  Discussion:
!
!    The overall tet mesh quality measure is the minimum of the corresponding
!    tetrahedron quality measure, over all tetrahedrons in the tet mesh.
!
!    Although tetrahedrons of order 10 are allowed as input,
!    only the first 4 nodes (presumed to be the vertices) are used
!    in the calculation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, double precision NODE_XYZ(3,NODE_NUM), the nodes.
!
!    Input, integer TETRA_ORDER, the order of the mesh, either 4 or 10.
!
!    Input, integer TETRA_NUM, the number of tetrahedrons.
!
!    Input, integer TETRA_NODE(TETRA_ORDER,TETRA_NUM), the nodes
!    that make up the tetrahedrons.
!
!    Output, double precision TETRA_QUALITY(TETRA_NUM), the quality
!    measure for each tetrahedron.
!
  implicit none

  integer , parameter :: dim_num = 3
  integer node_num
  integer tetra_num
  integer tetra_order

  double precision node_xyz(dim_num,node_num)
  integer tetra
  integer tetra_node(tetra_order,tetra_num)
  double precision tetra_quality(tetra_num)
  double precision tetra_xyz(dim_num,4)

  do tetra = 1, tetra_num

    tetra_xyz(1:dim_num,1:4) = node_xyz(1:dim_num,tetra_node(1:4,tetra))

    call tetrahedron_quality1_3d ( tetra_xyz, tetra_quality(tetra) )

  end do
end

subroutine tet_mesh_quality2 ( node_num, node_xyz, tetra_order, tetra_num, &
  tetra_node, tetra_quality )

!*****************************************************************************80
!
!! TET_MESH_QUALITY2 returns the quality of each tet in a mesh.
!
!  Discussion:
!
!    The overall tet mesh quality measure is the minimum of the
!    corresponding tetrahedron quality measure, over all tetrahedrons in the
!    tet mesh.
!
!    Although tetrahedrons of order 10 are allowed as input,
!    only the first 4 nodes (presumed to be the vertices) are used
!    in the calculation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, double precision NODE_XYZ(3,NODE_NUM), the nodes.
!
!    Input, integer TETRA_ORDER, the order of the mesh, either 4 or 10.
!
!    Input, integer TETRA_NUM, the number of tetrahedrons.
!
!    Input, integer TETRA_NODE(TETRA_ORDER,TETRA_NUM), the nodes
!    that make up the tetrahedrons.
!
!    Output, double precision TETRA_QUALITY(TETRA_NUM), the quality
!    measure for each tetrahedron.
!
  implicit none

  integer , parameter :: dim_num = 3
  integer node_num
  integer tetra_num
  integer tetra_order

  double precision node_xyz(dim_num,node_num)
  integer tetra
  integer tetra_node(tetra_order,tetra_num)
  double precision tetra_quality(tetra_num)
  double precision tetra_xyz(dim_num,4)

  do tetra = 1, tetra_num

    tetra_xyz(1:dim_num,1:4) = node_xyz(1:dim_num,tetra_node(1:4,tetra))

    call tetrahedron_quality2_3d ( tetra_xyz, tetra_quality(tetra) )

  end do
end

subroutine tet_mesh_quality3 ( node_num, node_xyz, tetra_order, tetra_num, &
  tetra_node, tetra_quality )

!*****************************************************************************80
!
!! TET_MESH_QUALITY3 returns the quality of each tet in a mesh.
!
!  Discussion:
!
!    The overall tet mesh quality measure is the minimum of the
!    corresponding tetrahedron quality measure, over all tetrahedrons in the
!    tet mesh.
!
!    Although tetrahedrons of order 10 are allowed as input,
!    only the first 4 nodes (presumed to be the vertices) are used
!    in the calculation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, double precision NODE_XYZ(3,NODE_NUM), the nodes.
!
!    Input, integer TETRA_ORDER, the order of the mesh, either 4 or 10.
!
!    Input, integer TETRA_NUM, the number of tetrahedrons.
!
!    Input, integer TETRA_NODE(TETRA_ORDER,TETRA_NUM), the nodes
!    that make up the tetrahedrons.
!
!    Output, double precision TETRA_QUALITY(TETRA_NUM), the quality
!    measure for each tetrahedron.
!
  implicit none

  integer , parameter :: dim_num = 3
  integer node_num
  integer tetra_num
  integer tetra_order

  double precision node_xyz(dim_num,node_num)
  integer tetra
  integer tetra_node(tetra_order,tetra_num)
  double precision tetra_quality(tetra_num)
  double precision tetra_xyz(dim_num,4)

  do tetra = 1, tetra_num

    tetra_xyz(1:dim_num,1:4) = node_xyz(1:dim_num,tetra_node(1:4,tetra))

    call tetrahedron_quality3_3d ( tetra_xyz, tetra_quality(tetra) )

  end do
end

subroutine tet_mesh_quality4 ( node_num, node_xyz, tetra_order, tetra_num, &
  tetra_node, tetra_quality )

!*****************************************************************************80
!
!! TET_MESH_QUALITY4 returns the quality of each tet in a mesh.
!
!  Discussion:
!
!    The overall tet mesh quality measure is the minimum of the
!    corresponding tetrahedron quality measure, over all tetrahedrons in the
!    tet mesh.
!
!    Although tetrahedrons of order 10 are allowed as input,
!    only the first 4 nodes (presumed to be the vertices) are used
!    in the calculation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, double precision NODE_XYZ(3,NODE_NUM), the nodes.
!
!    Input, integer TETRA_ORDER, the order of the mesh, either 4 or 10.
!
!    Input, integer TETRA_NUM, the number of tetrahedrons.
!
!    Input, integer TETRA_NODE(TETRA_ORDER,TETRA_NUM), the nodes
!    that make up the tetrahedrons.
!
!    Output, double precision TETRA_QUALITY(TETRA_NUM), the quality
!    measure for each tetrahedron.
!
  implicit none

  integer , parameter :: dim_num = 3
  integer node_num
  integer tetra_num
  integer tetra_order

  double precision node_xyz(dim_num,node_num)
  integer tetra
  integer tetra_node(tetra_order,tetra_num)
  double precision tetra_quality(tetra_num)
  double precision tetra_xyz(dim_num,4)

  do tetra = 1, tetra_num

    tetra_xyz(1:dim_num,1:4) = node_xyz(1:dim_num,tetra_node(1:4,tetra))

    call tetrahedron_quality4_3d ( tetra_xyz, tetra_quality(tetra) )

  end do
end

subroutine tet_mesh_quality5 ( node_num, node_xyz, tetra_order, tetra_num, &
  tetra_node, tetra_quality )

!*****************************************************************************80
!
!! TET_MESH_QUALITY5 returns the quality of each tet in a mesh.
!
!  Discussion:
!
!    The overall tet mesh quality measure is the ratio of the minimum
!    tetrahedron volume to the maximum tetrahedron volume.
!
!    Although tetrahedrons of order 10 are allowed as input,
!    only the first 4 nodes (presumed to be the vertices) are used
!    in the calculation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, double precision NODE_XYZ(3,NODE_NUM), the nodes.
!
!    Input, integer TETRA_ORDER, the order of the mesh, either 4 or 10.
!
!    Input, integer TETRA_NUM, the number of tetrahedrons.
!
!    Input, integer TETRA_NODE(TETRA_ORDER,TETRA_NUM), the nodes
!    that make up the tetrahedrons.
!
!    Output, double precision TETRA_QUALITY(TETRA_NUM), the quality
!    measure for each tetrahedron.
!
  implicit none

  integer , parameter :: dim_num = 3
  integer node_num
  integer tetra_num
  integer tetra_order

  double precision node_xyz(dim_num,node_num)
  integer tetra
  integer tetra_node(tetra_order,tetra_num)
  double precision tetra_quality(tetra_num)
  double precision tetra_xyz(dim_num,4)
  double precision volume_max

  do tetra = 1, tetra_num

    tetra_xyz(1:dim_num,1:4) = node_xyz(1:dim_num,tetra_node(1:4,tetra))

    call tetrahedron_volume_3d ( tetra_xyz, tetra_quality(tetra) )

  end do

  volume_max = maxval ( tetra_quality(1:tetra_num) )

  if ( 0.0D+00 < volume_max ) then
    tetra_quality(1:tetra_num) = tetra_quality(1:tetra_num) / volume_max
  end if
end

subroutine tetrahedron_circumsphere_3d ( tetra, r, pc )

!*****************************************************************************80
!
!! TETRAHEDRON_CIRCUMSPHERE_3D computes the circumsphere of a tetrahedron in 3D.
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
!    Adrian Bowyer and John Woodwark,
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

  integer , parameter :: dim_num = 3
  integer , parameter :: rhs_num = 1

  double precision a(dim_num,dim_num+rhs_num)
  integer i
  integer info
  integer j
  double precision pc(dim_num)
  double precision r
  double precision tetra(dim_num,4)
!
!  Set up the linear system.
!
  a(1:dim_num,1:3) = transpose ( tetra(1:dim_num,2:4) )

  do j = 1, dim_num
    a(1:dim_num,j) = a(1:dim_num,j) - tetra(j,1)
  end do

  do i = 1, 3
    a(i,4) = sum ( a(i,1:3)**2 )
  end do
!
!  Solve the linear system.
!
  call r8mat_solve ( dim_num, rhs_num, a, info )
!
!  If the system was singular, return a consolation prize.
!
  if ( info /= 0 ) then
    r = -1.0D+00
    pc(1:dim_num) = 0.0D+00
  end if
!
!  Compute the radius and center.
!
  r = 0.5D+00 * sqrt ( sum ( a(1:dim_num,4)**2 ) )

  pc(1:dim_num) = tetra(1:dim_num,1) + 0.5D+00 * a(1:dim_num,4)
end

subroutine tetrahedron_edge_length_3d ( tetra, edge_length )

!*****************************************************************************80
!
!! TETRAHEDRON_EDGE_LENGTH_3D returns edge lengths of a tetrahedron in 3D.
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

  integer , parameter :: dim_num = 3

  double precision r8vec_length
  double precision edge_length(6)
  integer j1
  integer j2
  integer k
  double precision tetra(dim_num,4)

  k = 0
  do j1 = 1, 3
    do j2 = j1+1, 4
      k = k + 1
      edge_length(k) = r8vec_length ( dim_num, &
        tetra(1:dim_num,j2) - tetra(1:dim_num,j1) )
    end do
  end do
end

subroutine tetrahedron_insphere_3d ( tetra, r, pc )

!*****************************************************************************80
!
!! TETRAHEDRON_INSPHERE_3D finds the insphere of a tetrahedron in 3D.
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
!    David Eberly,
!    Centers of a Simplex,
!    http://www.geometrictools.com
!
!  Parameters:
!
!    Input, double precision TETRA(3,4), the vertices of the tetrahedron.
!
!    Output, double precision R, PC(3), the radius and the center
!    of the sphere.
!
  implicit none

  integer , parameter :: dim_num = 3

  double precision b(4,4)
  double precision r8mat_det_4d
  double precision r8vec_length
  double precision gamma
  double precision l123
  double precision l124
  double precision l134
  double precision l234
  double precision n123(1:dim_num)
  double precision n124(1:dim_num)
  double precision n134(1:dim_num)
  double precision n234(1:dim_num)
  double precision pc(1:dim_num)
  double precision r
  double precision tetra(1:dim_num,4)
  double precision v21(1:dim_num)
  double precision v31(1:dim_num)
  double precision v41(1:dim_num)
  double precision v32(1:dim_num)
  double precision v42(1:dim_num)
  double precision v43(1:dim_num)

  v21(1:dim_num) = tetra(1:dim_num,2) - tetra(1:dim_num,1)
  v31(1:dim_num) = tetra(1:dim_num,3) - tetra(1:dim_num,1)
  v41(1:dim_num) = tetra(1:dim_num,4) - tetra(1:dim_num,1)
  v32(1:dim_num) = tetra(1:dim_num,3) - tetra(1:dim_num,2)
  v42(1:dim_num) = tetra(1:dim_num,4) - tetra(1:dim_num,2)
  v43(1:dim_num) = tetra(1:dim_num,4) - tetra(1:dim_num,3)

  call r8vec_cross_3d ( v21, v31, n123 )
  call r8vec_cross_3d ( v41, v21, n124 )
  call r8vec_cross_3d ( v31, v41, n134 )
  call r8vec_cross_3d ( v42, v32, n234 )

  l123 = r8vec_length ( dim_num, n123 )
  l124 = r8vec_length ( dim_num, n124 )
  l134 = r8vec_length ( dim_num, n134 )
  l234 = r8vec_length ( dim_num, n234 )

  pc(1:dim_num) = ( l234 * tetra(1:dim_num,1)   &
                  + l134 * tetra(1:dim_num,2)   &
                  + l124 * tetra(1:dim_num,3)   &
                  + l123 * tetra(1:dim_num,4) ) &
                / ( l234 + l134 + l124 + l123 )

  b(1:dim_num,1:4) = tetra(1:dim_num,1:4)
  b(4,1:4) = 1.0D+00

  gamma = abs ( r8mat_det_4d ( b ) )

! gamma = abs ( &
!     ( tetra(1,2) * tetra(2,3) * tetra(3,4) &
!     - tetra(1,3) * tetra(2,4) * tetra(3,2) &
!     + tetra(1,4) * tetra(2,2) * tetra(3,3) ) &
!   - ( tetra(1,1) * tetra(2,3) * tetra(3,4) &
!     - tetra(1,3) * tetra(2,4) * tetra(3,1) &
!     + tetra(1,4) * tetra(2,1) * tetra(3,3) ) &
!   + ( tetra(1,1) * tetra(2,2) * tetra(3,4) &
!     - tetra(1,2) * tetra(2,4) * tetra(3,1) &
!     + tetra(1,4) * tetra(2,1) * tetra(3,2) ) &
!   - ( tetra(1,1) * tetra(2,2) * tetra(3,3) &
!     - tetra(1,2) * tetra(2,3) * tetra(3,1) &
!     + tetra(1,3) * tetra(2,1) * tetra(3,2) ) )

  r = gamma / ( l234 + l134 + l124 + l123 )
end

subroutine tetrahedron_quality1_3d ( tetra, quality )

!*****************************************************************************80
!
!! TETRAHEDRON_QUALITY1_3D: "quality" of a tetrahedron in 3D.
!
!  Discussion:
!
!    The quality of a tetrahedron is 3.0 times the ratio of the radius of
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

  integer , parameter :: dim_num = 3

  double precision pc(dim_num)
  double precision quality
  double precision r_in
  double precision r_out
  double precision tetra(dim_num,4)

  call tetrahedron_circumsphere_3d ( tetra, r_out, pc )

  call tetrahedron_insphere_3d ( tetra, r_in, pc )

  quality = 3.0D+00 * r_in / r_out
end

subroutine tetrahedron_quality2_3d ( tetra, quality2 )

!*****************************************************************************80
!
!! TETRAHEDRON_QUALITY2_3D: "quality" of a tetrahedron in 3D.
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
!    Qiang Du and Desheng Wang,
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

  integer , parameter :: dim_num = 3

  double precision edge_length(6)
  double precision l_max
  double precision pc(dim_num)
  double precision quality2
  double precision r_in
  double precision tetra(dim_num,4)

  call tetrahedron_edge_length_3d ( tetra, edge_length )

  l_max = maxval ( edge_length(1:6) )

  call tetrahedron_insphere_3d ( tetra, r_in, pc )

  quality2 = 2.0D+00 * sqrt ( 6.0D+00 ) * r_in / l_max
end

subroutine tetrahedron_quality3_3d ( tetra, quality3 )

!******************************************************************************
!
!! TETRAHEDRON_QUALITY3_3D computes the mean ratio of a tetrahedron.
!
!  Discussion:
!
!    This routine computes QUALITY3, the eigenvalue or mean ratio of
!    a tetrahedron.
!
!      QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of squares of edge lengths).
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

  integer , parameter :: dim_num = 3

  double precision ab(dim_num)
  double precision ac(dim_num)
  double precision ad(dim_num)
  double precision bc(dim_num)
  double precision bd(dim_num)
  double precision cd(dim_num)
  double precision denom
  double precision lab
  double precision lac
  double precision lad
  double precision lbc
  double precision lbd
  double precision lcd
  double precision quality3
  double precision tetra(dim_num,4)
  double precision volume
!
!  Compute the vectors representing the sides of the tetrahedron.
!
  ab(1:3) = tetra(1:dim_num,2) - tetra(1:dim_num,1)
  ac(1:3) = tetra(1:dim_num,3) - tetra(1:dim_num,1)
  ad(1:3) = tetra(1:dim_num,4) - tetra(1:dim_num,1)
  bc(1:3) = tetra(1:dim_num,3) - tetra(1:dim_num,2)
  bd(1:3) = tetra(1:dim_num,4) - tetra(1:dim_num,2)
  cd(1:3) = tetra(1:dim_num,4) - tetra(1:dim_num,3)
!
!  Compute the squares of the lengths of the sides.
!
  lab = sum ( ab(1:dim_num)**2 )
  lac = sum ( ac(1:dim_num)**2 )
  lad = sum ( ad(1:dim_num)**2 )
  lbc = sum ( bc(1:dim_num)**2 )
  lbd = sum ( bd(1:dim_num)**2 )
  lcd = sum ( cd(1:dim_num)**2 )
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

subroutine tetrahedron_quality4_3d ( tetra, quality4 )

!******************************************************************************
!
!! TETRAHEDRON_QUALITY4_3D computes the minimum solid angle of a tetrahedron.
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

  integer , parameter :: dim_num = 3

  double precision a(dim_num)
  double precision ab(dim_num)
  double precision ac(dim_num)
  double precision ad(dim_num)
  double precision b(dim_num)
  double precision bc(dim_num)
  double precision bd(dim_num)
  double precision c(dim_num)
  double precision cd(dim_num)
  double precision d(dim_num)
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
  double precision tetra(dim_num,4)
  double precision volume
!
!  Compute the vectors that represent the sides.
!
  ab(1:dim_num) = tetra(1:dim_num,2) - tetra(1:dim_num,1)
  ac(1:dim_num) = tetra(1:dim_num,3) - tetra(1:dim_num,1)
  ad(1:dim_num) = tetra(1:dim_num,4) - tetra(1:dim_num,1)
  bc(1:dim_num) = tetra(1:dim_num,3) - tetra(1:dim_num,2)
  bd(1:dim_num) = tetra(1:dim_num,4) - tetra(1:dim_num,2)
  cd(1:dim_num) = tetra(1:dim_num,4) - tetra(1:dim_num,3)
!
!  Compute the lengths of the sides.
!
  lab = sqrt ( sum ( ab(1:dim_num)**2 ) )
  lac = sqrt ( sum ( ac(1:dim_num)**2 ) )
  lad = sqrt ( sum ( ad(1:dim_num)**2 ) )
  lbc = sqrt ( sum ( bc(1:dim_num)**2 ) )
  lbd = sqrt ( sum ( bd(1:dim_num)**2 ) )
  lcd = sqrt ( sum ( cd(1:dim_num)**2 ) )
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

subroutine tetrahedron_volume_3d ( tetra, volume )

!*****************************************************************************80
!
!! TETRAHEDRON_VOLUME_3D computes the volume of a tetrahedron in 3D.
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

  integer , parameter :: dim_num = 3

  double precision a(4,4)
  double precision r8mat_det_4d
  double precision tetra(dim_num,4)
  double precision volume

  a(1:dim_num,1:4) = tetra(1:dim_num,1:4)
  a(4,1:4) = 1.0D+00

  volume = abs ( r8mat_det_4d ( a ) ) / 6.0D+00
end
