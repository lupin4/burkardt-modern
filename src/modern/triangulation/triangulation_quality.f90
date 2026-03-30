!> triangulation_quality — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

subroutine bandwidth_mesh ( element_order, element_num, element_node, &
  ml, mu, m )

!*****************************************************************************80
!
!! BANDWIDTH_MESH: bandwidth of finite element mesh.
!
!  Discussion:
!
!    The quantity computed here is the "geometric" bandwidth determined
!    by the finite element mesh alone.
!
!    If a single finite element variable is associated with each node
!    of the mesh, and if the nodes and variables are numbered in the
!    same way, then the geometric bandwidth is the same as the bandwidth
!    of a typical finite element matrix.
!
!    The bandwidth M is defined in terms of the lower and upper bandwidths:
!
!      M = ML + 1 + MU
!
!    where
!
!      ML = maximum distance from any diagonal entry to a nonzero
!      entry in the same row, but earlier column,
!
!      MU = maximum distance from any diagonal entry to a nonzero
!      entry in the same row, but later column.
!
!    Because the finite element node adjacency relationship is symmetric,
!    we are guaranteed that ML = MU.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ELEMENT_ORDER, the order of the elements.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
!    ELEMENT_NODE(I,J) is the global index of local node I in element J.
!
!    Output, integer ML, MU, the lower and upper bandwidths 
!    of the matrix.
!
!    Output, integer M, the bandwidth of the matrix.
!
  implicit none

  integer element_num
  integer element_order

  integer element
  integer element_node(element_order,element_num)
  integer global_i
  integer global_j
  integer local_i
  integer local_j
  integer m
  integer ml
  integer mu

  ml = 0
  mu = 0

  do element = 1, element_num

    do local_i = 1, element_order
      global_i = element_node(local_i,element)

      do local_j = 1, element_order
        global_j = element_node(local_j,element)

        mu = max ( mu, global_j - global_i )
        ml = max ( ml, global_i - global_j )

      end do
    end do
  end do

  m = ml + 1 + mu
end
