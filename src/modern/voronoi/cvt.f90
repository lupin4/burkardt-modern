!> cvt — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).

module cvt_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use, intrinsic :: iso_c_binding,   only: c_int, c_double, c_float, c_bool
  implicit none
  private

  integer, parameter :: dp = real64
  integer, parameter :: sp = real32
  integer, parameter :: ip = int32

  public :: cvt, cvt_energy, cvt_iterate, cvt_sample, data_read, find_closest
  public :: halham_leap_check, halham_n_check, halham_dim_num_check, halham_seed_check, halham_step_check, halton_base_check
  public :: i4_to_halton_sequence, i4vec_transpose_print, prime, r8mat_uniform_01, s_cap, s_eqi
  public :: tuple_next_fast, user

contains

  subroutine cvt ( dim_num, n, batch, init, sample, sample_num, it_max, &
    it_fixed, seed, r, it_num, it_diff, energy ) &
        bind(C, name="cvt")

  !*****************************************************************************80
  !
  !! CVT computes a Centroidal Voronoi Tessellation.
  !
  !  Discussion:
  !
  !    This routine carries out the CVT iteration.
  !
  !    It initializes the CVT generators, unless the user indicates that
  !    this has already been done.
  !
  !    It sets a flag INITIALIZE that indicates whether the random number
  !    generator needs to be initialized.
  !
  !    It decides whether to use a new seed on each iteration, or to
  !    reuse a seed value.
  !
  !    It then repeatedly calls a routine that takes another step 
  !    of the iteration.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    23 June 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Qiang Du, Vance Faber, Max Gunzburger,
  !    Centroidal Voronoi Tessellations: Applications and Algorithms,
  !    SIAM Review, Volume 41, 1999, pages 637-676.
  !
  !  Parameters:
  !
  !    Input, integer(ip) DIM_NUM, the spatial dimension.
  !
  !    Input, integer(ip) N, the number of Voronoi cells.
  !
  !    Input, integer(ip) BATCH, sets the maximum number of sample points
  !    generated at one time.  It is inefficient to generate the sample
  !    points 1 at a time, but memory intensive to generate them all
  !    at once.  You might set BATCH to min ( SAMPLE_NUM, 10000 ), for instance.
  !    BATCH must be at least 1.
  !
  !    Input, integer(ip) INIT, specifies how the points are to be 
  !    initialized.
  !    -1, 'RANDOM', using FORTRAN RANDOM function;
  !     0, 'UNIFORM', using a simple uniform RNG;
  !     1, 'HALTON', from a Halton sequence;
  !     2, 'GRID', points from a grid;
  !     3, 'USER', call "user" routine;
  !     4, points are already initialized on input.
  !
  !    Input, integer(ip) SAMPLE, specifies how the sampling is done.
  !    -1, 'RANDOM', using FORTRAN RANDOM function;
  !     0, 'UNIFORM', using a simple uniform RNG;
  !     1, 'HALTON', from a Halton sequence;
  !     2, 'GRID', points from a grid;
  !     3, 'USER', call "user" routine.
  !
  !    Input, integer(ip) SAMPLE_NUM, the number of sample points.
  !
  !    Input, integer(ip) IT_MAX, the maximum number of iterations.
  !
  !    Input, integer(ip) IT_FIXED, the maximum number of iterations to 
  !    take with a fixed set of sample points.
  !
  !    Input/output, integer(ip) SEED, the current random number seed.
  !
  !    Input/output, real(dp) R(DIM_NUM,N), the approximate CVT points.
  !    If INIT = 4 on input, then it is assumed that these values have been
  !    initialized.  On output, the CVT iteration has been applied to improve
  !    the value of the points.
  !
  !    Output, integer(ip) IT_NUM, the number of iterations taken.  
  !    Generally, this will be equal to IT_MAX, unless the iteration tolerance was
  !    satisfied early.
  !
  !    Output, real(dp) IT_DIFF, the L2 norm of the difference
  !    between the iterates.
  !
  !    Output, real(dp) ENERGY, the discrete "energy", divided
  !    by the number of sample points.
  !

    integer(ip), intent(in), value :: dim_num
    integer(ip), intent(in), value :: n

    integer(ip), intent(in), value :: batch
    logical, parameter :: DEBUG = .true.
    real(dp), intent(out) :: energy
    integer(ip), intent(in), value :: init
    logical :: initialize
    real(dp), intent(out) :: it_diff
    integer(ip), intent(in), value :: it_fixed
    integer(ip), intent(in), value :: it_max
    integer(ip), intent(out) :: it_num
    real(dp), intent(inout) :: r(dim_num,n)
    integer(ip), intent(in), value :: sample
    integer(ip), intent(in), value :: sample_num
    integer(ip), intent(inout) :: seed
    integer(ip) :: seed_base
    integer(ip) :: seed_init

    if ( batch < 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT - Fatal error!'
      write ( *, '(a)' ) '  The input value BATCH < 1.'
      stop
    end if

    if ( seed <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT - Fatal error!'
      write ( *, '(a)' ) '  The input value SEED <= 0.'
      stop
    end if

    if ( DEBUG ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Step    SEED        L2-Change         Energy'
      write ( *, '(a)' ) ' '
    end if

    it_num = 0
    it_diff = 0.0_dp
    energy = 0.0_dp
    seed_init = seed
  !
  !  Initialize the data unless the user has already done that.
  !
    if ( init /= 4 ) then

      initialize = .true.

      call cvt_sample ( dim_num, n, n, init, initialize, seed, r )

    end if
  !
  !  NOTE THIS IS NEW:...
  !
    call cvt_energy ( dim_num, n, batch, sample, initialize, sample_num, &
      seed, r, energy )

    if ( DEBUG ) then
      write ( *, '(2x,i4,2x,i12,2x,14x,2x,g14.6)' ) &
        it_num, seed_init, energy
    end if
  !
  !  If the initialization and sampling steps use the same random number
  !  scheme, then the sampling scheme does not have to be initialized.
  !
    if ( init == sample ) then
      initialize = .false.
    else
      initialize = .true.
    end if
  !
  !  Carry out the iteration.
  !
    do while ( it_num < it_max )
  !
  !  If it's time to update the seed, save its current value
  !  as the starting value for all iterations in this cycle.
  !  If it's not time to update the seed, restore it to its initial
  !  value for this cycle.
  !
      if ( mod ( it_num, it_fixed ) == 0 ) then
        seed_base = seed
      else
        seed = seed_base
      end if

      it_num = it_num + 1

      seed_init = seed

      call cvt_iterate ( dim_num, n, batch, sample, initialize, sample_num, &
        seed, r, it_diff, energy )

      initialize = .false.

      if ( DEBUG ) then
        write ( *, '(2x,i4,2x,i12,2x,g14.6,2x,g14.6)' ) &
          it_num, seed_init, it_diff, energy
      end if

    end do
  end subroutine cvt

  subroutine cvt_energy ( dim_num, n, batch, sample, initialize, sample_num, &
    seed, r, energy ) &
        bind(C, name="cvt_energy")

  !*****************************************************************************80
  !
  !! CVT_ENERGY computes the CVT energy of a dataset.
  !
  !  Discussion:
  !
  !    For a given number of generators, a CVT is a minimizer (or at least
  !    a local minimizer) of the CVT energy.  During a CVT iteration,
  !    it should generally be the case that the CVT energy decreases from
  !    step to step, and that perturbations or adjustments of an
  !    approximate CVT will almost always have higher CVT energy.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    02 December 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) DIM_NUM, the spatial dimension.
  !
  !    Input, integer(ip) N, the number of generators.
  !
  !    Input, integer(ip) BATCH, the maximum number of sample points 
  !    to generate at one time.
  !
  !    Input, integer(ip) SAMPLE, specifies how the sampling is done.
  !    -1, 'RANDOM', using FORTRAN RANDOM function;
  !     0, 'UNIFORM', using a simple uniform RNG;
  !     1, 'HALTON', from a Halton sequence;
  !     2, 'GRID', points from a grid;
  !     3, 'USER', call "user" routine.
  !
  !    Input, logical INITIALIZE, is TRUE if the pseudorandom process should be
  !    reinitialized.
  !
  !    Input, integer(ip) SAMPLE_NUM, the number of sample points to use.
  !
  !    Input/output, integer(ip) SEED, a seed for the random 
  !    number generator.
  !
  !    Input, real(dp) R(DIM_NUM,N), the coordinates of the points.
  !
  !    Output, real(dp) ENERGY, the estimated CVT energy.
  !

    integer(ip), intent(in), value :: batch
    integer(ip), intent(in), value :: dim_num
    integer(ip), intent(in), value :: n

    real(dp), intent(out) :: energy
    integer(ip) :: get
    integer(ip) :: have
    logical, intent(in), value :: initialize
    integer(ip) :: j
    integer(ip) :: nearest(batch)
    real(dp), intent(in) :: r(dim_num,n)
    real(dp) :: s(dim_num,batch)
    integer(ip), intent(in), value :: sample
    integer(ip), intent(in), value :: sample_num
    integer(ip), intent(inout) :: seed

    have = 0
    energy = 0.0_dp

    do while ( have < sample_num )

      get = min ( sample_num - have, batch )

      call cvt_sample ( dim_num, sample_num, get, sample, initialize, seed, s )

      have = have + get

      call find_closest ( dim_num, n, get, s, r, nearest )

      do j = 1, get

        energy = energy &
          + sum ( ( s(1:dim_num,j) - r(1:dim_num,nearest(j)) )**2 )
      end do

    end do

    energy = energy / real ( sample_num, dp)
  end subroutine cvt_energy

  subroutine cvt_iterate ( dim_num, n, batch, sample, initialize, sample_num, &
    seed, r, it_diff, energy ) &
        bind(C, name="cvt_iterate")

  !*****************************************************************************80
  !
  !! CVT_ITERATE takes one step of the CVT iteration.
  !
  !  Discussion:
  !
  !    The routine is given a set of points, called "generators", which
  !    define a tessellation of the region into Voronoi cells.  Each point
  !    defines a cell.  Each cell, in turn, has a centroid, but it is
  !    unlikely that the centroid and the generator coincide.
  !
  !    Each time this CVT iteration is carried out, an attempt is made
  !    to modify the generators in such a way that they are closer and
  !    closer to being the centroids of the Voronoi cells they generate.
  !
  !    A large number of sample points are generated, and the nearest generator
  !    is determined.  A count is kept of how many points were nearest to each
  !    generator.  Once the sampling is completed, the location of all the
  !    generators is adjusted.  This step should decrease the discrepancy
  !    between the generators and the centroids.
  !
  !    The centroidal Voronoi tessellation minimizes the "energy",
  !    defined to be the integral, over the region, of the square of
  !    the distance between each point in the region and its nearest generator.
  !    The sampling technique supplies a discrete estimate of this
  !    energy.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    14 September 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Qiang Du, Vance Faber, Max Gunzburger,
  !    Centroidal Voronoi Tessellations: Applications and Algorithms,
  !    SIAM Review, Volume 41, 1999, pages 637-676.
  !
  !  Parameters:
  !
  !    Input, integer(ip) DIM_NUM, the spatial dimension.
  !
  !    Input, integer(ip) N, the number of points to generate.
  !
  !    Input, integer(ip) BATCH, sets the maximum number of sample points
  !    generated at one time.  It is inefficient to generate the sample
  !    points 1 at a time, but memory intensive to generate them all
  !    at once.  You might set BATCH to min ( SAMPLE_NUM, 10000 ), for instance.
  !    BATCH must be at least 1.
  !
  !    Input, integer(ip) SAMPLE, specifies how the sampling is done.
  !    -1, 'RANDOM', using FORTRAN RANDOM function;
  !     0, 'UNIFORM', using a simple uniform RNG;
  !     1, 'HALTON', from a Halton sequence;
  !     2, 'GRID', points from a grid;
  !     3, 'USER', call "user" routine.
  !
  !    Input, logical INITIALIZE, is TRUE if the random number generator
  !    should be initialized, because this is the first call to it.
  !
  !    Input, integer(ip) SAMPLE_NUM, the number of sample points.
  !
  !    Input/output, integer(ip) SEED, the random number seed.
  !
  !    Input/output, real(dp) R(DIM_NUM,N), the Voronoi
  !    cell generators.  On output, these have been modified
  !
  !    Output, real(dp) IT_DIFF, the L2 norm of the difference
  !    between the iterates.
  !
  !    Output, real(dp) ENERGY, the discrete "energy", divided
  !    by the number of sample points.
  !

    integer(ip), intent(in), value :: batch
    integer(ip), intent(in), value :: dim_num
    integer(ip), intent(in), value :: n
    integer(ip), intent(in), value :: sample_num

    integer(ip) :: count(n)
    real(dp), intent(out) :: energy
    integer(ip) :: get
    integer(ip) :: have
    logical, intent(out) :: initialize
    real(dp), intent(out) :: it_diff
    integer(ip) :: j
    integer(ip) :: nearest(batch)
    real(dp), intent(inout) :: r(dim_num,n)
    real(dp) :: r2(dim_num,n)
    real(dp) :: s(dim_num,batch)
    integer(ip), intent(in), value :: sample
    integer(ip), intent(inout) :: seed
  !
  !  Take each generator as the first sample point for its region.
  !  This can slightly slow the convergence, but it simplifies the
  !  algorithm by guaranteeing that no region is completely missed
  !  by the sampling.
  !
    energy = 0.0_dp
    r2(1:dim_num,1:n) = r(1:dim_num,1:n)
    count(1:n) = 1
  !
  !  Generate the sampling points S in batches.
  !
    have = 0

    do while ( have < sample_num )

      get = min ( sample_num - have, batch )

      call cvt_sample ( dim_num, sample_num, get, sample, initialize, seed, s )

      initialize = .false.
      have = have + get
  !
  !  Find the index N of the nearest cell generator to each sample point S.
  !
      call find_closest ( dim_num, n, get, s, r, nearest )
  !
  !  Add S to the centroid associated with generator N.
  !
      do j = 1, get
        r2(1:dim_num,nearest(j)) = r2(1:dim_num,nearest(j)) + s(1:dim_num,j)
        energy = energy + sum ( ( r(1:dim_num,nearest(j)) - s(1:dim_num,j) )**2 )
        count(nearest(j)) = count(nearest(j)) + 1
      end do

    end do
  !
  !  Estimate the centroids.
  !
    do j = 1, n
      r2(1:dim_num,j) = r2(1:dim_num,j) / real ( count(j), dp)
    end do
  !
  !  Determine the sum of the distances between the old generators 
  !  and the estimated centroids.
  !
    it_diff = 0.0_dp
    do j = 1, n
      it_diff = it_diff + sqrt ( sum ( ( r2(1:dim_num,j) - r(1:dim_num,j) )**2 ) )
    end do
  !
  !  Replace the generators by the centroids.
  !
    r(1:dim_num,1:n) = r2(1:dim_num,1:n)
  !
  !  Normalize the discrete energy estimate.
  !
    energy = energy / real ( sample_num, dp) 
  end subroutine cvt_iterate

  subroutine cvt_sample ( dim_num, n, n_now, sample, initialize, seed, r ) &
        bind(C, name="cvt_sample")

  !*****************************************************************************80
  !
  !! CVT_SAMPLE returns sample points.
  !
  !  Discussion:
  !
  !    N sample points are to be taken from the region, of spatial
  !    dimension DIM_NUM.  In most cases, this region is a unit box.
  !
  !    These sample points are usually created by a pseudorandom process
  !    for which the points are essentially indexed by a quantity called
  !    SEED.  To get N sample points, we generate values with indices
  !    SEED through SEED+N-1.
  !
  !    It may not be practical to generate all the sample points in a 
  !    single call.  For that reason, the routine allows the user to
  !    request a total of N points, but to require that only N_NOW be
  !    generated now (on this call).  
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    23 June 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) DIM_NUM, the spatial dimension.
  !
  !    Input, integer(ip) N, the number of sample points to be generated.
  !
  !    Input, integer(ip) N_NOW, the number of sample points to be 
  !    generated on this call.  N_NOW must be at least 1.
  !
  !    Input, integer(ip) SAMPLE, specifies how the sampling is done.
  !    -1, 'RANDOM', using FORTRAN RANDOM function;
  !     0, 'UNIFORM', using a simple uniform RNG;
  !     1, 'HALTON', from a Halton sequence;
  !     2, 'GRID', points from a grid;
  !     3, 'USER', from the "user" routine.
  !
  !    Input, logical INITIALIZE, is TRUE if the pseudorandom process should be
  !    reinitialized.
  !
  !    Input/output, integer(ip) SEED, the random number seed.
  !
  !    Output, real(dp) R(DIM_NUM,N_NOW), the sample points.
  !

    integer(ip), intent(in), value :: dim_num
    integer(ip), intent(in), value :: n_now

    real(dp) :: exponent
    integer(ip), allocatable, dimension ( : ) :: halton_base
    integer(ip), allocatable, dimension ( : ) :: halton_leap
    integer(ip), allocatable, dimension ( : ) :: halton_seed
    integer(ip) :: halton_step
    integer(ip) :: i
    logical, intent(in), value :: initialize
    integer(ip) :: j
    integer(ip), intent(in), value :: n
    integer(ip) :: ngrid
    integer(ip) :: prime
    real(dp), intent(out) :: r(dim_num,n_now)
    integer(ip) :: rank
    integer(ip) :: rank_max
    integer(ip), intent(in), value :: sample
    integer(ip), intent(inout) :: seed
    integer(ip), allocatable, dimension ( : ) :: tuple

    if ( n_now < 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_SAMPLE - Fatal error!'
      write ( *, '(a)' ) '  N_NOW < 1.'
      stop
    end if

    if ( sample == -1 ) then

      if ( initialize ) then
        call random_initialize ( seed )
      end if

      call random_number ( harvest = r(1:dim_num,1:n_now) )

      seed = seed + n_now * dim_num

    else if ( sample == 0 ) then

      call r8mat_uniform_01 ( dim_num, n_now, seed, r )

    else if ( sample == 1 ) then

      allocate ( halton_seed(1:dim_num) )
      allocate ( halton_leap(1:dim_num) )
      allocate ( halton_base(1:dim_num) )

      halton_step = seed
      halton_seed(1:dim_num) = 0
      halton_leap(1:dim_num) = 1

      do i = 1, dim_num
        halton_base(i) = prime ( i )
      end do

      call i4_to_halton_sequence ( dim_num, n_now, halton_step, halton_seed, &
        halton_leap, halton_base, r(1:dim_num,1:n_now) )

      deallocate ( halton_seed )
      deallocate ( halton_leap )
      deallocate ( halton_base )

      seed = seed + n_now

    else if ( sample == 2 ) then

      allocate ( tuple(1:dim_num) )

      exponent = 1.0_dp / real ( dim_num, dp)
      ngrid = int ( ( real ( n, dp) )**exponent )
      rank_max = ngrid**dim_num

      if ( rank_max < n ) then
        ngrid = ngrid + 1
        rank_max = ngrid**dim_num
      end if

      if ( initialize ) then
        rank = -1
        call tuple_next_fast ( ngrid, dim_num, rank, tuple )
      end if

      rank = mod ( seed, rank_max )

      do j = 1, n_now
        call tuple_next_fast ( ngrid, dim_num, rank, tuple )
        rank = rank + 1
        rank = mod ( rank, rank_max )
        r(1:dim_num,j) = real ( 2 * tuple(1:dim_num) - 1, dp) &
          / real ( 2 * ngrid, dp)
      end do

      seed = seed + n_now

      deallocate ( tuple )

    else if ( sample == 3 ) then

      call user ( dim_num, n_now, seed, r )

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_SAMPLE - Fatal error!'
      write ( *, '(a,i6,a)' ) '  The value of SAMPLE = ', sample, ' is illegal.'
      stop

    end if
  end subroutine cvt_sample

  subroutine data_read ( file_in_name, dim_num, n, r, success ) &
        bind(C, name="data_read")

  !*****************************************************************************80
  !
  !! DATA_READ reads generator coordinate data from a file.
  !
  !  Discussion:
  !
  !    The file is assumed to contain one record per line.
  !
  !    Records beginning with the '#' character are comments, and are ignored.
  !    Blank lines are also ignored.
  !
  !    Each line that is not ignored is assumed to contain exactly (or at least)
  !    M real numbers, representing the coordinates of a point.
  !
  !    There are assumed to be exactly (or at least) N such records.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    25 July 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) FILE_IN_NAME, the name of the input file.
  !
  !    Input, integer(ip) DIM_NUM, the number of spatial dimensions.
  !
  !    Input, integer(ip) N, the number of points.  The program
  !    will stop reading data once N values have been read.
  !
  !    Output, real(dp) R(DIM_NUM,N), the point coordinates.
  !
  !    Output, logical SUCCESS, is TRUE if the data was read properly.
  !

    integer(ip), intent(in), value :: dim_num
    integer(ip), intent(in), value :: n

    character ( len = * ), intent(in) :: file_in_name
    integer(ip) :: file_in_unit
    integer(ip) :: i
    integer(ip) :: ierror
    integer(ip) :: ios
    character ( len = 255 ) :: line
    real(dp), intent(out) :: r(dim_num,n)
    logical, intent(out) :: success
    real(dp) :: x(dim_num)

    success = .true.

    call get_unit ( file_in_unit )

    open ( unit = file_in_unit, file = file_in_name, status = 'old', &
      iostat = ios )

    if ( ios /= 0 ) then
      success = .false.
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Could not open the input file: ' // &
        trim ( file_in_name )
    end if

    i = 0

    do while ( i < n )

      read ( file_in_unit, '(a)', iostat = ios ) line

      if ( ios /= 0 ) then
        success = .false.
      end if

      if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
        cycle
      end if

      call s_to_r8vec ( line, dim_num, x, ierror )

      if ( ierror /= 0 ) then
        cycle
      end if

      i = i + 1

      r(1:dim_num,i) = x(1:dim_num)

    end do

    close ( unit = file_in_unit )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_READ:'
    write ( *, '(a,i6)' ) '  Read coordinate data from file.'
  end subroutine data_read

  pure subroutine find_closest ( dim_num, n, sample_num, s, r, nearest ) &
        bind(C, name="find_closest")

  !*****************************************************************************80
  !
  !! FIND_CLOSEST finds the nearest R point to each S point.
  !
  !  Discussion:
  !
  !    This routine finds the closest Voronoi cell generator by checking every
  !    one.  For problems with many cells, this process can take the bulk
  !    of the CPU time.  Other approaches, which group the cell generators into
  !    bins, can run faster by a large factor.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    02 August 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) DIM_NUM, the spatial dimension.
  !
  !    Input, integer(ip) N, the number of cell generators.
  !
  !    Input, integer(ip) SAMPLE_NUM, the number of sample points.
  !
  !    Input, real(dp) S(DIM_NUM,SAMPLE_NUM), the points to be checked.
  !
  !    Input, real(dp) R(DIM_NUM,N), the cell generators.
  !
  !    Output, integer(ip) NEAREST(SAMPLE_NUM), the index of the nearest 
  !    cell generators.
  !

    integer(ip), intent(in), value :: dim_num
    integer(ip), intent(in), value :: n
    integer(ip), intent(in), value :: sample_num

    real(dp) :: dist_sq_min
    real(dp) :: dist_sq
    integer(ip) :: jr
    integer(ip) :: js
    integer(ip), intent(out) :: nearest(sample_num)
    real(dp), intent(in) :: r(dim_num,n)
    real(dp), intent(in) :: s(dim_num,sample_num)

    do js = 1, sample_num

      dist_sq_min = huge ( dist_sq_min )
      nearest(js) = -1

      do jr = 1, n

        dist_sq = sum ( ( r(1:dim_num,jr) - s(1:dim_num,js) )**2 )

        if ( dist_sq < dist_sq_min ) then
          dist_sq_min = dist_sq
          nearest(js) = jr
        end if

      end do

    end do
  end subroutine find_closest

  function halham_leap_check ( dim_num, leap ) &
        bind(C, name="halham_leap_check")

  !*****************************************************************************80
  !
  !! HALHAM_LEAP_CHECK checks LEAP for a Halton or Hammersley sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    21 September 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) DIM_NUM, the spatial dimension.
  !
  !    Input, integer(ip) LEAP(DIM_NUM), the leap vector.
  !
  !    Output, logical, HALHAM_LEAP_CHECK, true if LEAP is legal.
  !

    integer(ip), intent(in), value :: dim_num

    logical :: halham_leap_check
    integer(ip), intent(in) :: leap(dim_num)

    if ( any ( leap(1:dim_num) < 1 ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HALHAM_LEAP_CHECK - Fatal error!'
      write ( *, '(a)' ) '  Some entry of LEAP < 1!'
      write ( *, '(a)' ) ' '
      call i4vec_transpose_print ( dim_num, leap, 'LEAP:  ' )
      halham_leap_check = .false.
    else
      halham_leap_check = .true.
    end if
  end function halham_leap_check

  function halham_n_check ( n ) &
        bind(C, name="halham_n_check")

  !*****************************************************************************80
  !
  !! HALHAM_N_CHECK checks N for a Halton or Hammersley sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    21 September 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the spatial dimension.
  !
  !    Output, logical HALHAM_N_CHECK, true if N is legal.
  !

    logical :: halham_n_check
    integer(ip), intent(in), value :: n

    if ( n < 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HALHAM_N_CHECK - Fatal error!'
      write ( *, '(a)' ) '  N < 1.'
      write ( *, '(a,i12)' ) '  N = ', n
      halham_n_check = .false.
    else
      halham_n_check = .true.
    end if
  end function halham_n_check

  function halham_dim_num_check ( dim_num ) &
        bind(C, name="halham_dim_num_check")

  !*****************************************************************************80
  !
  !! HALHAM_DIM_NUM_CHECK checks DIM_NUM for a Halton or Hammersley sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    21 September 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) DIM_NUM, the spatial dimension.
  !
  !    Output, logical HALHAM_DIM_NUM_CHECK, true if DIM_NUM is legal.
  !

    logical :: halham_dim_num_check
    integer(ip), intent(in), value :: dim_num

    if ( dim_num < 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HALHAM_DIM_NUM_CHECK - Fatal error!'
      write ( *, '(a)' ) '  DIM_NUM < 1.'
      write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
      halham_dim_num_check = .false.
    else
      halham_dim_num_check = .true.
    end if
  end function halham_dim_num_check

  function halham_seed_check ( dim_num, seed ) &
        bind(C, name="halham_seed_check")

  !*****************************************************************************80
  !
  !! HALHAM_SEED_CHECK checks SEED for a Halton or Hammersley sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    21 September 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) DIM_NUM, the spatial dimension.
  !
  !    Input, integer(ip) SEED(DIM_NUM), the seed vector.
  !
  !    Output, logical, HALHAM_SEED_CHECK, true if SEED is legal.
  !

    integer(ip), intent(in), value :: dim_num

    logical :: halham_seed_check
    integer(ip), intent(in) :: seed(dim_num)

    if ( any ( seed(1:dim_num) < 0 ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HALHAM_SEED_CHECK - Fatal error!'
      write ( *, '(a)' ) '  Some entry of SEED < 0!'
      write ( *, '(a)' ) ' '
      call i4vec_transpose_print ( dim_num, seed, 'SEED:  ' )
      halham_seed_check = .false.
    else
      halham_seed_check = .true.
    end if
  end function halham_seed_check

  function halham_step_check ( step ) &
        bind(C, name="halham_step_check")

  !*****************************************************************************80
  !
  !! HALHAM_STEP_CHECK checks STEP for a Halton or Hammersley sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    21 September 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) STEP, the index of the subsequence element.
  !
  !    Output, logical HALHAM_STEP_CHECK, true if STEP is legal.
  !

    logical :: halham_step_check
    integer(ip), intent(in), value :: step

    if ( step < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HALHAM_STEP_CHECK - Fatal error!'
      write ( *, '(a)' ) '  STEP < 0.'
      write ( *, '(a,i12)' ) '  STEP = ', step
      halham_step_check = .false.
    else
      halham_step_check = .true.
    end if
  end function halham_step_check

  function halton_base_check ( dim_num, base ) &
        bind(C, name="halton_base_check")

  !*****************************************************************************80
  !
  !! HALTON_BASE_CHECK checks BASE for a Halton sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    21 September 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) DIM_NUM, the spatial dimension.
  !
  !    Input, integer(ip) BASE(DIM_NUM), the bases.
  !
  !    Output, logical, HALTON_BASE_CHECK, true if BASE is legal.
  !

    integer(ip), intent(in), value :: dim_num

    integer(ip), intent(in) :: base(dim_num)
    logical :: halton_base_check

    if ( any ( base(1:dim_num) <= 1 ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HALTON_BASE_CHECK - Fatal error!'
      write ( *, '(a)' ) '  Some entry of BASE is <= 1!'
      write ( *, '(a)' ) ' '
      call i4vec_transpose_print ( dim_num, base, 'BASE:  ' )
      halton_base_check = .false.
    else
      halton_base_check = .true.
    end if
  end function halton_base_check

  subroutine i4_to_halton_sequence ( dim_num, n, step, seed, leap, base, r ) &
        bind(C, name="i4_to_halton_sequence")

  !*****************************************************************************80
  !
  !! I4_TO_HALTON_SEQUENCE computes N elements of a leaped Halton subsequence.
  !
  !  Discussion:
  !
  !    The DIM_NUM-dimensional Halton sequence is really DIM_NUM separate
  !    sequences, each generated by a particular base.
  !
  !    This routine selects elements of a "leaped" subsequence of the
  !    Halton sequence.  The subsequence elements are indexed by a
  !    quantity called STEP, which starts at 0.  The STEP-th subsequence
  !    element is simply element
  !
  !      SEED(1:DIM_NUM) + STEP * LEAP(1:DIM_NUM)
  !
  !    of the original Halton sequence.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    21 September 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    J H Halton,
  !    On the efficiency of certain quasi-random sequences of points
  !    in evaluating multi-dimensional integrals,
  !    Numerische Mathematik,
  !    Volume 2, 1960, pages 84-90.
  !
  !    J H Halton and G B Smith,
  !    Algorithm 247: Radical-Inverse Quasi-Random Point Sequence,
  !    Communications of the ACM,
  !    Volume 7, 1964, pages 701-702.
  !
  !    Ladislav Kocis and William Whiten,
  !    Computational Investigations of Low-Discrepancy Sequences,
  !    ACM Transactions on Mathematical Software,
  !    Volume 23, Number 2, 1997, pages 266-294.
  !
  !  Parameters:
  !
  !    Input, integer(ip) DIM_NUM, the spatial dimension.
  !    1 <= DIM_NUM is required.
  !
  !    Input, integer(ip) N, the number of elements of the sequence.
  !
  !    Input, integer(ip) STEP, the index of the subsequence element.
  !    0 <= STEP is required.
  !
  !    Input, integer(ip) SEED(DIM_NUM), the Halton sequence index 
  !    corresponding to STEP = 0.
  !
  !    Input, integer(ip) LEAP(DIM_NUM), the succesive jumps in the 
  !    Halton sequence.
  !
  !    Input, integer(ip) BASE(DIM_NUM), the Halton bases.
  !
  !    Output, real(dp) R(DIM_NUM,N), the next N elements of the
  !    leaped Halton subsequence, beginning with element STEP.
  !

    integer(ip), intent(in), value :: dim_num
    integer(ip), intent(in), value :: n

    integer(ip), intent(in) :: base(dim_num)
    real(dp) :: base_inv
    integer(ip) :: digit(n)
    logical :: halham_leap_check
    logical :: halham_n_check
    logical :: halham_dim_num_check
    logical :: halham_seed_check
    logical :: halham_step_check
    logical :: halton_base_check
    integer(ip) :: i
    integer(ip) :: j
    integer(ip), intent(in) :: leap(dim_num)
    real(dp), intent(out) :: r(dim_num,n)
    integer(ip), intent(in) :: seed(dim_num)
    integer(ip) :: seed2(n)
    integer(ip), intent(in), value :: step
  !
  !  Check the input.
  !
    if ( .not. halham_dim_num_check ( dim_num ) ) then
      stop
    end if

    if ( .not. halham_n_check ( n ) ) then
      stop
    end if

    if ( .not. halham_step_check ( step ) ) then
      stop
    end if

    if ( .not. halham_seed_check ( dim_num, seed ) ) then
      stop
    end if

    if ( .not. halham_leap_check ( dim_num, leap ) ) then
      stop
    end if

    if ( .not. halton_base_check ( dim_num, base ) ) then
      stop
    end if
  !
  !  Calculate the data.
  !
    r(1:dim_num,1:n) = 0.0_dp

    do i = 1, dim_num

      do j = 1, n
        seed2(j) = seed(i) + ( step + j - 1 ) * leap(i)
      end do

      base_inv = real ( 1.0_dp, dp) / real ( base(i), dp)

      do while ( any ( seed2(1:n) /= 0 ) )
        digit(1:n) = mod ( seed2(1:n), base(i) )
        r(i,1:n) = r(i,1:n) + real ( digit(1:n), dp) * base_inv
        base_inv = base_inv / real ( base(i), dp)
        seed2(1:n) = seed2(1:n) / base(i)
      end do

    end do
  end subroutine i4_to_halton_sequence

  subroutine i4vec_transpose_print ( n, a, title ) &
        bind(C, name="i4vec_transpose_print")

  !*****************************************************************************80
  !
  !! I4VEC_TRANSPOSE_PRINT prints an integer vector "transposed".
  !
  !  Example:
  !
  !    A = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 /)
  !    TITLE = 'My vector:  '
  !
  !    My vector:      1    2    3    4    5
  !                    6    7    8    9   10
  !                   11
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    03 July 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the number of components of the vector.
  !
  !    Input, integer(ip) A(N), the vector to be printed.
  !
  !    Input, character ( len = * ) TITLE, a title to be printed first.
  !    TITLE may be blank.
  !

    integer(ip), intent(in), value :: n

    integer(ip), intent(in) :: a(n)
    integer(ip) :: ihi
    integer(ip) :: ilo
    character ( len = 11 ) :: string
    character ( len = * ), intent(in) :: title
    integer(ip) :: title_len

    if ( 0 < len_trim ( title ) ) then

      title_len = len ( title )

      write ( string, '(a,i3,a)' ) '(', title_len, 'x,5i12)'

      do ilo = 1, n, 5
        ihi = min ( ilo + 5 - 1, n )
        if ( ilo == 1 ) then
          write ( *, '(a, 5i12)' ) title, a(ilo:ihi)
        else
          write ( *, string      )        a(ilo:ihi)
        end if
      end do

    else

      do ilo = 1, n, 5
        ihi = min ( ilo + 5 - 1, n )
        write ( *, '(5i12)' ) a(ilo:ihi)
      end do

    end if
  end subroutine i4vec_transpose_print

  function prime ( n ) &
        bind(C, name="prime")

  !*****************************************************************************80
  !
  !! PRIME returns any of the first PRIME_MAX prime numbers.
  !
  !  Discussion:
  !
  !    PRIME_MAX is 1600, and the largest prime stored is 13499.
  !
  !    Thanks to Bart Vandewoestyne for pointing out a typo, 18 February 2005.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    18 February 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Milton Abramowitz and Irene Stegun,
  !    Handbook of Mathematical Functions,
  !    US Department of Commerce, 1964, pages 870-873.
  !
  !    Daniel Zwillinger,
  !    CRC Standard Mathematical Tables and Formulae,
  !    30th Edition,
  !    CRC Press, 1996, pages 95-98.
  !
  !  Parameters:
  !
  !    Input, integer(ip) N, the index of the desired prime number.
  !    In general, is should be true that 0 <= N <= PRIME_MAX.
  !    N = -1 returns PRIME_MAX, the index of the largest prime available.
  !    N = 0 is legal, returning PRIME = 1.
  !
  !    Output, integer(ip) PRIME, the N-th prime.  If N is out of range,
  !    PRIME is returned as -1.
  !

    integer(ip), parameter :: prime_max = 1600

    integer(ip), save :: icall = 0
    integer(ip), intent(in), value :: n
    integer(ip), save, dimension ( prime_max ) :: npvec
    integer(ip) :: prime

    if ( icall == 0 ) then

      icall = 1

      npvec(1:100) = (/ &
          2,    3,    5,    7,   11,   13,   17,   19,   23,   29, &
         31,   37,   41,   43,   47,   53,   59,   61,   67,   71, &
         73,   79,   83,   89,   97,  101,  103,  107,  109,  113, &
        127,  131,  137,  139,  149,  151,  157,  163,  167,  173, &
        179,  181,  191,  193,  197,  199,  211,  223,  227,  229, &
        233,  239,  241,  251,  257,  263,  269,  271,  277,  281, &
        283,  293,  307,  311,  313,  317,  331,  337,  347,  349, &
        353,  359,  367,  373,  379,  383,  389,  397,  401,  409, &
        419,  421,  431,  433,  439,  443,  449,  457,  461,  463, &
        467,  479,  487,  491,  499,  503,  509,  521,  523,  541 /)

      npvec(101:200) = (/ &
        547,  557,  563,  569,  571,  577,  587,  593,  599,  601, &
        607,  613,  617,  619,  631,  641,  643,  647,  653,  659, &
        661,  673,  677,  683,  691,  701,  709,  719,  727,  733, &
        739,  743,  751,  757,  761,  769,  773,  787,  797,  809, &
        811,  821,  823,  827,  829,  839,  853,  857,  859,  863, &
        877,  881,  883,  887,  907,  911,  919,  929,  937,  941, &
        947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013, &
       1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, &
       1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, &
       1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223 /)

      npvec(201:300) = (/ &
       1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, &
       1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, &
       1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, &
       1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, &
       1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, &
       1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, &
       1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, &
       1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, &
       1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, &
       1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987 /)

      npvec(301:400) = (/ &
       1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, &
       2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, &
       2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, &
       2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, &
       2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, &
       2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, &
       2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, &
       2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, &
       2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, &
       2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741 /)

      npvec(401:500) = (/ &
       2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, &
       2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, &
       2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, &
       3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, &
       3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, &
       3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, &
       3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, &
       3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, &
       3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, &
       3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571 /)

      npvec(501:600) = (/ &
       3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, &
       3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, &
       3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, &
       3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, &
       3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, &
       4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, &
       4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, &
       4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, &
       4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, &
       4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409 /)

      npvec(601:700) = (/ &
       4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, &
       4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, &
       4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, &
       4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, &
       4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, &
       4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, &
       4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, &
       5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, &
       5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, &
       5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279 /)

      npvec(701:800) = (/ &
       5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, &
       5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, &
       5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, &
       5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639, &
       5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, &
       5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, &
       5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, &
       5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, &
       5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, &
       6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133 /)

      npvec(801:900) = (/ &
       6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, &
       6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, &
       6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, &
       6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, &
       6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, &
       6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, &
       6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, &
       6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, &
       6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, &
       6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997 /)

      npvec(901:1000) = (/ &
       7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, &
       7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, &
       7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, &
       7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, &
       7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, &
       7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, &
       7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, &
       7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, &
       7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, &
       7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919 /)

      npvec(1001:1100) = (/ &
       7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017, &
       8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111, &
       8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219, &
       8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291, &
       8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387, &
       8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501, &
       8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597, &
       8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677, &
       8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, &
       8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831 /)

      npvec(1101:1200) = (/ &
       8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929, &
       8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011, &
       9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109, &
       9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199, &
       9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283, &
       9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377, &
       9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439, &
       9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533, &
       9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631, &
       9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733 /)

      npvec(1201:1300) = (/ &
       9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811, &
       9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887, &
       9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973,10007, &
      10009,10037,10039,10061,10067,10069,10079,10091,10093,10099, &
      10103,10111,10133,10139,10141,10151,10159,10163,10169,10177, &
      10181,10193,10211,10223,10243,10247,10253,10259,10267,10271, &
      10273,10289,10301,10303,10313,10321,10331,10333,10337,10343, &
      10357,10369,10391,10399,10427,10429,10433,10453,10457,10459, &
      10463,10477,10487,10499,10501,10513,10529,10531,10559,10567, &
      10589,10597,10601,10607,10613,10627,10631,10639,10651,10657 /)

      npvec(1301:1400) = (/ &
      10663,10667,10687,10691,10709,10711,10723,10729,10733,10739, &
      10753,10771,10781,10789,10799,10831,10837,10847,10853,10859, &
      10861,10867,10883,10889,10891,10903,10909,10937,10939,10949, &
      10957,10973,10979,10987,10993,11003,11027,11047,11057,11059, &
      11069,11071,11083,11087,11093,11113,11117,11119,11131,11149, &
      11159,11161,11171,11173,11177,11197,11213,11239,11243,11251, &
      11257,11261,11273,11279,11287,11299,11311,11317,11321,11329, &
      11351,11353,11369,11383,11393,11399,11411,11423,11437,11443, &
      11447,11467,11471,11483,11489,11491,11497,11503,11519,11527, &
      11549,11551,11579,11587,11593,11597,11617,11621,11633,11657 /)

      npvec(1401:1500) = (/ &
      11677,11681,11689,11699,11701,11717,11719,11731,11743,11777, &
      11779,11783,11789,11801,11807,11813,11821,11827,11831,11833, &
      11839,11863,11867,11887,11897,11903,11909,11923,11927,11933, &
      11939,11941,11953,11959,11969,11971,11981,11987,12007,12011, &
      12037,12041,12043,12049,12071,12073,12097,12101,12107,12109, &
      12113,12119,12143,12149,12157,12161,12163,12197,12203,12211, &
      12227,12239,12241,12251,12253,12263,12269,12277,12281,12289, &
      12301,12323,12329,12343,12347,12373,12377,12379,12391,12401, &
      12409,12413,12421,12433,12437,12451,12457,12473,12479,12487, &
      12491,12497,12503,12511,12517,12527,12539,12541,12547,12553 /)

     npvec(1501:1600) = (/ &
      12569,12577,12583,12589,12601,12611,12613,12619,12637,12641, &
      12647,12653,12659,12671,12689,12697,12703,12713,12721,12739, &
      12743,12757,12763,12781,12791,12799,12809,12821,12823,12829, &
      12841,12853,12889,12893,12899,12907,12911,12917,12919,12923, &
      12941,12953,12959,12967,12973,12979,12983,13001,13003,13007, &
      13009,13033,13037,13043,13049,13063,13093,13099,13103,13109, &
      13121,13127,13147,13151,13159,13163,13171,13177,13183,13187, &
      13217,13219,13229,13241,13249,13259,13267,13291,13297,13309, &
      13313,13327,13331,13337,13339,13367,13381,13397,13399,13411, &
      13417,13421,13441,13451,13457,13463,13469,13477,13487,13499 /)

    end if

    if ( n == -1 ) then
      prime = prime_max
    else if ( n == 0 ) then
      prime = 1
    else if ( n <= prime_max ) then
      prime = npvec(n)
    else
      prime = -1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PRIME - Fatal error!'
      write ( *, '(a,i6)' ) '  Illegal prime index N = ', n
      write ( *, '(a,i6)' ) '  N should be between 1 and PRIME_MAX =', prime_max
      stop
    end if
  end function prime

  pure subroutine r8mat_uniform_01 ( m, n, seed, r ) &
        bind(C, name="r8mat_uniform_01")

  !*****************************************************************************80
  !
  !! R8MAT_UNIFORM_01 fills an array with pseudorandom numbers.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    11 August 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Paul Bratley, Bennett Fox, L E Schrage,
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
  !  Parameters:
  !
  !    Input, integer(ip)M, N, the number of rows and columns in 
  !    the array.
  !
  !    Input/output, integer(ip)SEED, the "seed" value, which should 
  !    NOT be 0.  On output, SEED has been updated.
  !
  !    Output, real(dp) R(M,N), the array of pseudorandom values.
  !

    integer(ip), intent(in), value :: m
    integer(ip), intent(in), value :: n

    integer(ip) :: i
    integer(ip) :: j
    integer(ip) :: k
    integer(ip), intent(out) :: seed
    real(dp), intent(out) :: r(m,n)

    do j = 1, n

      do i = 1, m

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed < 0 ) then
          seed = seed + 2147483647
        end if

        r(i,j) = real ( seed, dp) * 4.656612875E-10

      end do
    end do
  end subroutine r8mat_uniform_01

  subroutine s_cap ( s ) &
        bind(C, name="s_cap")

  !*****************************************************************************80
  !
  !! S_CAP replaces any lowercase letters by uppercase ones in a string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    28 June 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character ( len = * ) S, the string to be transformed.
  !

    character :: c
    integer(ip) :: i
    integer :: nchar
    character ( len = * ), intent(inout) :: s

    nchar = len_trim ( s )

    do i = 1, nchar

      c = s(i:i)
      call ch_cap ( c )
      s(i:i) = c

    end do
  end subroutine s_cap

  function s_eqi ( s1, s2 ) &
        bind(C, name="s_eqi")

  !*****************************************************************************80
  !
  !! S_EQI is a case insensitive comparison of two strings for equality.
  !
  !  Example:
  !
  !    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, S2, the strings to compare.
  !
  !    Output, logical S_EQI, the result of the comparison.
  !

    character :: c1
    character :: c2
    integer(ip) :: i
    integer(ip) :: len1
    integer(ip) :: len2
    integer(ip) :: lenc
    logical :: s_eqi
    character ( len = * ), intent(in) :: s1
    character ( len = * ), intent(in) :: s2

    len1 = len ( s1 )
    len2 = len ( s2 )
    lenc = min ( len1, len2 )

    s_eqi = .false.

    do i = 1, lenc

      c1 = s1(i:i)
      c2 = s2(i:i)
      call ch_cap ( c1 )
      call ch_cap ( c2 )

      if ( c1 /= c2 ) then
      end if

    end do

    do i = lenc + 1, len1
      if ( s1(i:i) /= ' ' ) then
      end if
    end do

    do i = lenc + 1, len2
      if ( s2(i:i) /= ' ' ) then
      end if
    end do

    s_eqi = .true.
  end function s_eqi

  subroutine tuple_next_fast ( m, n, rank, x ) &
        bind(C, name="tuple_next_fast")

  !*****************************************************************************80
  !
  !! TUPLE_NEXT_FAST computes the next element of a tuple space, "fast".
  !
  !  Discussion:
  !
  !    The elements are N vectors.  Each entry is constrained to lie
  !    between 1 and M.  The elements are produced one at a time.
  !    The first element is
  !      (1,1,...,1)
  !    and the last element is
  !      (M,M,...,M)
  !    Intermediate elements are produced in lexicographic order.
  !
  !    This code was written as a possibly faster version of TUPLE_NEXT.
  !
  !  Example:
  !
  !    N = 2,
  !    M = 3
  !
  !    INPUT        OUTPUT
  !    -------      -------
  !    Rank          X
  !    ----          ----
  !   -1            -1 -1
  !
  !    0             1  1
  !    1             1  2
  !    2             1  3
  !    3             2  1
  !    4             2  2
  !    5             2  3
  !    6             3  1
  !    7             3  2
  !    8             3  3
  !    9             1  1
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    11 August 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) M, the maximum entry in any component.
  !    M must be greater than 0.
  !
  !    Input, integer(ip) N, the number of components.
  !    N must be greater than 0.
  !
  !    Input, integer(ip) RANK, indicates the rank of the tuple.
  !    Typically, 0 <= RANK < N**M.  Values of RANK greater than
  !    N**M are legal and meaningful; they are equivalent to the
  !    corresponding value mod (N**M).  If RANK < 0, this indicates 
  !    that this is the first call for the given values of (M,N).  
  !    Initialization is done, and X is set to a dummy value.
  !
  !    Output, integer(ip) X(N), the next tuple, or a dummy value if
  !    initialization has just been done.
  !

    integer(ip), intent(in), value :: n

    integer(ip), save, allocatable, dimension ( : ) :: base
    integer(ip) :: i
    integer(ip), intent(in), value :: m
    integer(ip), intent(in), value :: rank
    integer(ip), intent(out) :: x(n)

    if ( rank < 0 ) then

      if ( m <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TUPLE_NEXT_FAST - Fatal error!'
        write ( *, '(a)' ) '  The value M <= 0 is not allowed.'
        write ( *, '(a,i6)' ) '  M = ', m
        stop
      end if

      if ( n <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TUPLE_NEXT_FAST - Fatal error!'
        write ( *, '(a)' ) '  The value N <= 0 is not allowed.'
        write ( *, '(a,i6)' ) '  N = ', n
        stop
      end if

      if ( allocated ( base ) ) then
        deallocate ( base )
      end if
      allocate ( base(1:n) )

      base(n) = 1
      do i = n-1, 1, -1
        base(i) = base(i+1) * m
      end do

      x(1:n) = -1

    else

      x(1:n) = mod ( rank / base(1:n), m ) + 1

    end if
  end subroutine tuple_next_fast

  subroutine user ( dim_num, n, seed, r ) &
        bind(C, name="user")

  !*****************************************************************************80
  !
  !! USER samples points in a user-specified region with given density.
  !
  !  Discussion:
  !
  !    This routine can be used to 
  !
  !    * specify an interesting initial configuration for the data,
  !      by specifing that USER be used for initialization (INIT = 3);
  !
  !    * specify the shape of the computational region, by specifying
  !      that sample points are to be generated by this routine, 
  !      (SAMPLE = 3) and then returning sample points uniformly at random.
  !
  !    * specify the distribution or density function, by specifying
  !      that sample points are to be generated by this routine, 
  !      (SAMPLE = 3 ) and then returning sample points according to a 
  !      given probability density function.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    23 June 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer(ip) DIM_NUM, the spatial dimension.
  !
  !    Input, integer(ip) N, the number of sample points desired.
  !
  !    Input/output, integer(ip) SEED, the "seed" value.  The meaning and
  !    use of this variable is up to the user.
  !
  !    Output, real(dp) R(DIM_NUM,N), an array of sample points from
  !    the region.
  !

    integer(ip), intent(in), value :: dim_num
    integer(ip), intent(in), value :: n

    real(dp) :: angle(n)
    integer(ip) :: i
    integer(ip) :: j
    integer(ip) :: k
    real(dp), parameter :: pi = 3.141592653589793_dp
    integer(ip), intent(inout) :: seed
    real(dp), intent(out) :: r(dim_num,n)
    real(dp) :: radius(n)
  !
  !  We sample points in the unit circle uniformly.
  !
    call random_number ( harvest = angle(1:n) )
    angle(1:n) = 2.0_dp * pi * angle(1:n)

    call random_number ( harvest = radius(1:n) )
    radius(1:n) = sqrt ( radius(1:n) )

    r(1,1:n) = radius(1:n) * cos ( angle(1:n) )
    r(2,1:n) = radius(1:n) * sin ( angle(1:n) )
  end subroutine user

end module cvt_mod
