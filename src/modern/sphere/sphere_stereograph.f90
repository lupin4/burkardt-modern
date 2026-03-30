!> sphere_stereograph — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

subroutine sphere_stereograph ( m, n, p, q )

!*****************************************************************************80
!
!! SPHERE_STEREOGRAPH computes the stereographic image of points on a sphere.
!
!  Discussion:
!
!    We start with a sphere of radius 1 and center (0,0,0).
!
!    The north pole N = (0,0,1) is the point of tangency to the sphere
!    of a plane, and the south pole S = (0,0,-1) is the focus for the
!    stereographic projection.
!
!    For any point P on the sphere, the stereographic projection Q of the
!    point is defined by drawing the line from S through P, and computing
!    Q as the intersection of this line with the plane.
!
!    Actually, we allow the spatial dimension M to be arbitrary.  Values
!    of M make sense starting with 2.  The north and south poles are
!    selected as the points (0,0,...,+1) and (0,0,...,-1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    C F Marcus,
!    The stereographic projection in vector notation,
!    Mathematics Magazine,
!    Volume 39, Number 2, March 1966, pages 100-102.
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, double precision P(M,N), a set of points on the unit sphere.
!
!    Output, double precision Q(M,N), the coordinates of the
!    image points.
!
  implicit none

  integer m
  integer n

  integer i
  integer j
  double precision p(m,n)
  double precision q(m,n)

  do j = 1, n
    do i = 1, m - 1
      q(i,j) = 2.0D+00 * p(i,j) / ( 1.0D+00 + p(m,j) )
    end do
    q(m,j) = 1.0D+00
  end do
end

subroutine sphere_stereograph_inverse ( m, n, q, p )

!*****************************************************************************80
!
!! SPHERE_STEREOGRAPH_INVERSE computes stereographic preimages of points.
!
!  Discussion:
!
!    We start with a sphere of radius 1 and center (0,0,0).
!
!    The north pole N = (0,0,1) is the point of tangency to the sphere
!    of a plane, and the south pole S = (0,0,-1) is the focus for the
!    stereographic projection.
!
!    For any point Q on the plane, the stereographic inverse projection
!    P of the point is defined by drawing the line from S through Q, and
!    computing P as the intersection of this line with the sphere.
!
!    Actually, we allow the spatial dimension M to be arbitrary.  Values
!    of M make sense starting with 2.  The north and south poles are
!    selected as the points (0,0,...,+1) and (0,0,...,-1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    C F Marcus,
!    The stereographic projection in vector notation,
!    Mathematics Magazine,
!    Volume 39, Number 2, March 1966, pages 100-102.
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, double precision Q(M,N), the points, which are presumed to lie
!    on the plane Z = 1.
!
!    Output, double precision P(M,N), the stereographic
!    inverse projections of the points.
!
  implicit none

  integer m
  integer n

  integer j
  double precision p(m,n)
  double precision q(m,n)
  double precision qn

  do j = 1, n

    qn = sum ( q(1:m-1,j)**2 )

    p(1:m-1,j) = 4.0D+00 * q(1:m-1,j) / ( 4.0D+00 + qn )

    p(m,j) = ( 4.0D+00 - qn ) / ( 4.0D+00 + qn )

  end do
end

subroutine sphere_stereograph2 ( m, n, p, focus, center, q )

!*****************************************************************************80
!
!! SPHERE_STEREOGRAPH2 computes the stereographic image of points on a sphere.
!
!  Discussion:
!
!    We start with a sphere of center C.
!
!    F is a point on the sphere which is the focus of the mapping,
!    and the antipodal point 2*C-F is the point of tangency
!    to the sphere of a plane.
!
!    For any point P on the sphere, the stereographic projection Q of the
!    point is defined by drawing the line from F through P, and computing
!    Q as the intersection of this line with the plane.
!
!    The spatial dimension M is arbitrary, but should be at least 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    C F Marcus,
!    The stereographic projection in vector notation,
!    Mathematics Magazine,
!    Volume 39, Number 2, March 1966, pages 100-102.
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer the number of points.
!
!    Input, double precision P[M*N], a set of points on the unit sphere.
!
!    Input, double precision FOCUS[M], the coordinates of the focus point.
!
!    Input, double precision CENTER[M], the coordinates of the center of 
!    the sphere.
!
!    Output, double precision Q[M*N], the coordinates of the
!    image points,
!
  implicit none

  integer m
  integer n

  double precision center(m)
  double precision cf_dot_pf
  double precision cf_normsq
  double precision focus(m)
  integer i
  integer j
  double precision p(m,n)
  double precision q(m,n)
  double precision s

  do j = 1, n 
    cf_normsq = 0.0D+00
    cf_dot_pf = 0.0D+00
    do i = 1, m
      cf_normsq = cf_normsq + ( center(i) - focus(i) ) ** 2
      cf_dot_pf = cf_dot_pf + ( center(i) - focus(i) ) * ( p(i,j) - focus(i) )
    end do
    s = 2.0D+00 * cf_normsq / cf_dot_pf
    do i = 1, m
      q(i,j) = s * p(i,j) + ( 1.0D+00 - s ) * focus(i)
    end do
  end do
end

subroutine sphere_stereograph2_inverse ( m, n, q, focus, center, p )

!*****************************************************************************80
!
!! SPHERE_STEREOGRAPH2_INVERSE computes stereographic preimages of points.
!
!  Discussion:
!
!    We start with a sphere of center C.
!
!    F is a point on the sphere which is the focus of the mapping,
!    and the antipodal point 2*C-F is the point of tangency
!    to the sphere of a plane.
!
!    For any point Q on the plane, the stereographic inverse projection
!    P of the point is defined by drawing the line from F through Q, and
!    computing P as the intersection of this line with the sphere.
!
!    The spatial dimension M is arbitrary, but should be at least 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    C F Marcus,
!    The stereographic projection in vector notation,
!    Mathematics Magazine,
!    Volume 39, Number 2, March 1966, pages 100-102.
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, double precision Q(M,N), the points, which are presumed to lie
!    on the plane.
!
!    Input, double precision FOCUS(M), the coordinates of the focus point.
!
!    Input, double precision CENTER(M), the coordinates of the center 
!    of the sphere.
!
!    Output, double precision P(M,N), the stereographic
!    inverse projections of the points.
!
  implicit none

  integer m
  integer n

  double precision center(m)
  double precision cf_dot_qf
  double precision focus(m)
  integer i
  integer j
  double precision p(m,n)
  double precision q(m,n)
  double precision qf_normsq
  double precision s

  do j = 1, n

    cf_dot_qf = 0.0D+00
    qf_normsq = 0.0D+00
    do i = 1, m
      cf_dot_qf = cf_dot_qf + ( center(i) - focus(i) ) * ( q(i,j) - focus(i) )
      qf_normsq = qf_normsq + ( q(i,j) - focus(i) ) ** 2
    end do

    s = 2.0D+00 * cf_dot_qf / qf_normsq
    do i = 1, m
      p(i,j) = s * q(i,j) + ( 1.0D+00 - s ) * focus(i)
    end do

  end do
end

subroutine uniform_on_sphere01_map ( dim_num, n, seed, x )

!*****************************************************************************80
!
!! UNIFORM_ON_SPHERE01_MAP maps uniform points onto the unit sphere.
!
!  Discussion:
!
!    The sphere has center 0 and radius 1.
!
!    This procedure is valid for any spatial dimension DIM_NUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Russell Cheng,
!    Random Variate Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998, pages 168.
!
!    Reuven Rubinstein,
!    Monte Carlo Optimization, Simulation, and Sensitivity
!    of Queueing Networks,
!    Krieger, 1992,
!    ISBN: 0894647644,
!    LC: QA298.R79.
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Input, integer N, the number of points.
!
!    Input/output, integer SEED, a seed for the random
!    number generator.
!
!    Output, double precision X(DIM_NUM,N), the points.
!
  implicit none

  integer dim_num
  integer n

  integer j
  double precision norm
  integer seed
  double precision x(dim_num,n)
!
!  Fill a matrix with normally distributed values.
!
  call r8mat_normal_01 ( dim_num, n, seed, x )
!
!  Normalize each column.
!
  do j = 1, n
!
!  Compute the length of the vector.
!
    norm = sqrt ( sum ( x(1:dim_num,j)**2 ) )
!
!  Normalize the vector.
!
    x(1:dim_num,j) = x(1:dim_num,j) / norm

  end do
end
