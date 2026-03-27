!> simplex_coordinates -- Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module simplex_coordinates_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: r8_factorial, r8mat_det, simplex_coordinates1, simplex_coordinates2, simplex_volume

contains

  pure function r8_factorial ( n ) &
        bind(C, name="r8_factorial")

  !*****************************************************************************80
  !
  !! R8_FACTORIAL computes the factorial of N.
  !
  !  Discussion:
  !
  !    factorial ( N ) = product ( 1 <= I <= N ) I
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    16 January 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the argument of the factorial function.
  !    If N is less than 1, the function value is returned as 1.
  !
  !    Output, real(dp) R8_FACTORIAL, the factorial of N.
  !

    real(dp) :: r8_factorial
    integer(ip), intent(in), value :: n
    integer(ip) :: i

    r8_factorial = 1.0_dp

    do i = 1, n
      r8_factorial = r8_factorial * real ( i, dp)
    end do
  end function r8_factorial

  pure subroutine r8mat_det ( n, a, det ) &
        bind(C, name="r8mat_det")

  !*****************************************************************************80
  !
  !! R8MAT_DET computes the determinant of an R8MAT.
  !
  !  Discussion:
  !
  !    An R8MAT is an array of R8 values.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 December 2004
  !
  !  Author:
  !
  !    Original FORTRAN77 version by Helmut Spaeth.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Helmut Spaeth,
  !    Cluster Analysis Algorithms
  !    for Data Reduction and Classification of Objects,
  !    Ellis Horwood, 1980, page 125-127.
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the order of the matrix.
  !
  !    Input, real(dp) A(N,N), the matrix whose determinant is desired.
  !
  !    Output, real(dp) DET, the determinant of the matrix.
  !

    integer(ip), intent(in), value :: n
    real(dp), intent(in) :: a(n,n)
    real(dp), intent(out) :: det
    real(dp) :: b(n,n)
    integer(ip) :: j
    integer(ip) :: k
    integer(ip) :: m
    integer(ip) :: piv(1)
    real(dp) :: t

    b(1:n,1:n) = a(1:n,1:n)

    det = 1.0_dp

    do k = 1, n

      piv = maxloc ( abs ( b(k:n,k) ) )

      m = piv(1) + k - 1

      if ( m /= k ) then
        det = - det
        t      = b(m,k)
        b(m,k) = b(k,k)
        b(k,k) = t
      end if

      det = det * b(k,k)

      if ( b(k,k) /= 0.0_dp ) then

        b(k+1:n,k) = -b(k+1:n,k) / b(k,k)

        do j = k + 1, n
          if ( m /= k ) then
            t      = b(m,j)
            b(m,j) = b(k,j)
            b(k,j) = t
          end if
          b(k+1:n,j) = b(k+1:n,j) + b(k+1:n,k) * b(k,j)
        end do

      end if

    end do
  end subroutine r8mat_det

  pure subroutine simplex_coordinates1 ( n, x ) &
        bind(C, name="simplex_coordinates1")

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
  !    Input, integer(ip) N, the spatial dimension.
  !
  !    Output, real(dp) X(N,N+1), the coordinates of the vertices
  !    of a simplex in N dimensions.
  !

    integer(ip), intent(in), value :: n
    real(dp), intent(out) :: x(n,n+1)
    integer(ip) :: i
    integer(ip) :: j

    x(1:n,1:n+1) = 0.0_dp

    do i = 1, n
  !
  !  Set X(I,I) so that sum ( X(1:I,I)**2 ) = 1.
  !
      x(i,i) = sqrt ( 1.0_dp - sum ( x(1:i-1,i)**2 ) )
  !
  !  Set X(I,J) for J = I+1 to N+1 by using the fact that XI dot XJ = - 1 / N
  !
      do j = i + 1, n + 1
        x(i,j) = ( - 1.0_dp / real ( n, dp) &
          - dot_product ( x(1:i-1,i), x(1:i-1,j) ) ) / x(i,i)
      end do

    end do
  end subroutine simplex_coordinates1

  pure subroutine simplex_coordinates2 ( n, x ) &
        bind(C, name="simplex_coordinates2")

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
  !    Input, integer(ip) N, the spatial dimension.
  !
  !    Output, real(dp) X(N,N+1), the coordinates of the vertices
  !    of a simplex in N dimensions.
  !

    integer(ip), intent(in), value :: n
    real(dp), intent(out) :: x(n,n+1)
    real(dp) :: a
    real(dp) :: c(n)
    integer(ip) :: j
    real(dp) :: s

    x(1:n,1:n+1) = 0.0_dp

    do j = 1, n
      x(j,j) = 1.0_dp
    end do

    a = ( 1.0_dp - sqrt ( 1.0_dp + real ( n, dp) ) ) &
      / real ( n, dp)

    x(1:n,n+1) = a
  !
  !  Now adjust coordinates so the centroid is at zero.
  !
    c(1:n) = sum ( x(1:n,1:n+1), dim = 2 ) / real ( n + 1, dp)

    do j = 1, n + 1
      x(1:n,j) = x(1:n,j) - c(1:n)
    end do
  !
  !  Now scale so each column has norm 1.
  !
    s = sqrt ( sum ( x(1:n,1)**2 ) )

    x(1:n,1:n+1) = x(1:n,1:n+1) / s
  end subroutine simplex_coordinates2

  pure subroutine simplex_volume ( n, x, volume ) &
        bind(C, name="simplex_volume")

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
  !    Input, integer(ip) N, the spatial dimension.
  !
  !    Input, real(dp) X(N,N+1), the coordinates of the vertices
  !    of a simplex in N dimensions.
  !
  !    Output, real(dp) VOLUME, the volume of the simplex.
  !

    integer(ip), intent(in), value :: n
    real(dp), intent(in) :: x(n,n+1)
    real(dp), intent(out) :: volume
    real(dp) :: a(n,n)
    real(dp) :: det
    integer(ip) :: i
    integer(ip) :: j

    a(1:n,1:n) = x(1:n,1:n)
    do j = 1, n
      a(1:n,j) = a(1:n,j) - x(1:n,n+1)
    end do

    call r8mat_det ( n, a, det )

    volume = abs ( det )
    do i = 1, n
      volume = volume / real ( i, dp)
    end do
  end subroutine simplex_volume

end module simplex_coordinates_mod
