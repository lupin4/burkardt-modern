!> simplex_coordinates — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

subroutine simplex_coordinates1 ( n, x )

!*****************************************************************************80
!
!! SIMPLEX_COORDINATES1 computes the Cartesian coordinates of simplex vertices.
!
!  Discussion:
!
!    The simplex will have its centroid at 0;
!
!    The sum of the vertices will be zero.
!
!    The distance of each vertex from the origin will be 1.
!
!    The length of each edge will be constant.
!
!    The dot product of the vectors defining any two vertices will be - 1 / N.
!    This also means the angle subtended by the vectors from the origin
!    to any two distinct vertices will be arccos ( - 1 / N ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the spatial dimension.
!
!    Output, double precision X(N,N+1), the coordinates of the vertices
!    of a simplex in N dimensions.  
!
  implicit none

  integer n

  integer i
  integer j
  double precision x(n,n+1)

  x(1:n,1:n+1) = 0.0D+00

  do i = 1, n
!
!  Set X(I,I) so that sum ( X(1:I,I)**2 ) = 1.
!
    x(i,i) = sqrt ( 1.0D+00 - sum ( x(1:i-1,i)**2 ) )
!
!  Set X(I,J) for J = I+1 to N+1 by using the fact that XI dot XJ = - 1 / N 
!
    do j = i + 1, n + 1
      x(i,j) = ( - 1.0D+00 / real ( n) &
        - dot_product ( x(1:i-1,i), x(1:i-1,j) ) ) / x(i,i)
    end do

  end do
end

subroutine simplex_coordinates2 ( n, x )

!*****************************************************************************80
!
!! SIMPLEX_COORDINATES2 computes the Cartesian coordinates of simplex vertices.
!
!  Discussion:
!
!    This routine uses a simple approach to determining the coordinates of
!    the vertices of a regular simplex in n dimensions.
!
!    We want the vertices of the simplex to satisfy the following conditions:
!
!    1) The centroid, or average of the vertices, is 0.
!    2) The distance of each vertex from the centroid is 1.
!       By 1), this is equivalent to requiring that the sum of the squares
!       of the coordinates of any vertex be 1.
!    3) The distance between any pair of vertices is equal (and is not zero.)
!    4) The dot product of any two coordinate vectors for distinct vertices
!       is -1/N; equivalently, the angle subtended by two distinct vertices
!       from the centroid is arccos ( -1/N).
!
!    Note that if we choose the first N vertices to be the columns of the
!    NxN identity matrix, we are almost there.  By symmetry, the last column
!    must have all entries equal to some value A.  Because the square of the
!    distance between the last column and any other column must be 2 (because
!    that's the distance between any pair of columns), we deduce that
!    (A-1)^2 + (N-1)*A^2 = 2, hence A = (1-sqrt(1+N))/N.  Now compute the 
!    centroid C of the vertices, and subtract that, to center the simplex 
!    around the origin.  Finally, compute the norm of one column, and rescale 
!    the matrix of coordinates so each vertex has unit distance from the origin.
!
!    This approach devised by John Burkardt, 19 September 2010.  What,
!    I'm not the first?
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the spatial dimension.
!
!    Output, double precision X(N,N+1), the coordinates of the vertices
!    of a simplex in N dimensions.  
!
  implicit none

  integer n

  double precision a
  double precision c(n)
  integer j
  double precision s
  double precision x(n,n+1)

  x(1:n,1:n+1) = 0.0D+00

  do j = 1, n
    x(j,j) = 1.0D+00
  end do

  a = ( 1.0D+00 - sqrt ( 1.0D+00 + real ( n) ) ) &
    / real ( n)

  x(1:n,n+1) = a
!
!  Now adjust coordinates so the centroid is at zero.
!
  c(1:n) = sum ( x(1:n,1:n+1), dim = 2 ) / real ( n + 1)

  do j = 1, n + 1
    x(1:n,j) = x(1:n,j) - c(1:n)
  end do
!
!  Now scale so each column has norm 1.
!
  s = sqrt ( sum ( x(1:n,1)**2 ) )

  x(1:n,1:n+1) = x(1:n,1:n+1) / s
end

subroutine simplex_volume ( n, x, volume )

!*****************************************************************************80
!
!! SIMPLEX_VOLUME computes the volume of a simplex.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the spatial dimension.
!
!    Input, double precision X(N,N+1), the coordinates of the vertices
!    of a simplex in N dimensions.  
!
!    Output, double precision VOLUME, the volume of the simplex.
!
  implicit none

  integer n

  double precision a(n,n)
  double precision det
  integer i
  integer j
  double precision volume
  double precision x(n,n+1)

  a(1:n,1:n) = x(1:n,1:n)
  do j = 1, n
    a(1:n,j) = a(1:n,j) - x(1:n,n+1)
  end do

  call r8mat_det ( n, a, det )

  volume = abs ( det )
  do i = 1, n
    volume = volume / real ( i)
  end do
end
