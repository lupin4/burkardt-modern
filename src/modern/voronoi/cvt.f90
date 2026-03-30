!> cvt — Modern Fortran 2018
!>
!> Modernized from John Burkardt's original (GNU LGPL).
!> Standalone routines (no module wrapping) for clean C symbol names.

subroutine cvt ( dim_num, n, batch, init, sample, sample_num, it_max, &
  it_fixed, seed, r, it_num, it_diff, energy )

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
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer N, the number of Voronoi cells.
!
!    Input, integer BATCH, sets the maximum number of sample points
!    generated at one time.  It is inefficient to generate the sample
!    points 1 at a time, but memory intensive to generate them all
!    at once.  You might set BATCH to min ( SAMPLE_NUM, 10000 ), for instance.
!    BATCH must be at least 1.
!
!    Input, integer INIT, specifies how the points are to be 
!    initialized.
!    -1, 'RANDOM', using FORTRAN RANDOM function;
!     0, 'UNIFORM', using a simple uniform RNG;
!     1, 'HALTON', from a Halton sequence;
!     2, 'GRID', points from a grid;
!     3, 'USER', call "user" routine;
!     4, points are already initialized on input.
!
!    Input, integer SAMPLE, specifies how the sampling is done.
!    -1, 'RANDOM', using FORTRAN RANDOM function;
!     0, 'UNIFORM', using a simple uniform RNG;
!     1, 'HALTON', from a Halton sequence;
!     2, 'GRID', points from a grid;
!     3, 'USER', call "user" routine.
!
!    Input, integer SAMPLE_NUM, the number of sample points.
!
!    Input, integer IT_MAX, the maximum number of iterations.
!
!    Input, integer IT_FIXED, the maximum number of iterations to 
!    take with a fixed set of sample points.
!
!    Input/output, integer SEED, the current random number seed.
!
!    Input/output, double precision R(DIM_NUM,N), the approximate CVT points.
!    If INIT = 4 on input, then it is assumed that these values have been
!    initialized.  On output, the CVT iteration has been applied to improve
!    the value of the points.
!
!    Output, integer IT_NUM, the number of iterations taken.  
!    Generally, this will be equal to IT_MAX, unless the iteration tolerance was
!    satisfied early.
!
!    Output, double precision IT_DIFF, the L2 norm of the difference
!    between the iterates.
!
!    Output, double precision ENERGY, the discrete "energy", divided
!    by the number of sample points.
!
  implicit none

  integer dim_num
  integer n

  integer batch
  logical, parameter :: DEBUG = .true.
  double precision energy
  integer init
  logical initialize
  double precision it_diff
  integer it_fixed
  integer it_max
  integer it_num
  double precision r(dim_num,n)
  integer sample
  integer sample_num
  integer seed
  integer seed_base
  integer seed_init

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
  it_diff = 0.0D+00
  energy = 0.0D+00
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
end

subroutine cvt_energy ( dim_num, n, batch, sample, initialize, sample_num, &
  seed, r, energy )

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
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer N, the number of generators.
!
!    Input, integer BATCH, the maximum number of sample points 
!    to generate at one time.
!
!    Input, integer SAMPLE, specifies how the sampling is done.
!    -1, 'RANDOM', using FORTRAN RANDOM function;
!     0, 'UNIFORM', using a simple uniform RNG;
!     1, 'HALTON', from a Halton sequence;
!     2, 'GRID', points from a grid;
!     3, 'USER', call "user" routine.
!
!    Input, logical INITIALIZE, is TRUE if the pseudorandom process should be
!    reinitialized.
!
!    Input, integer SAMPLE_NUM, the number of sample points to use.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Input, double precision R(DIM_NUM,N), the coordinates of the points.
!
!    Output, double precision ENERGY, the estimated CVT energy.
!
  implicit none

  integer batch
  integer dim_num
  integer n

  double precision energy
  integer get
  integer have
  logical initialize
  integer j
  integer nearest(batch)
  double precision r(dim_num,n)
  double precision s(dim_num,batch)
  integer sample
  integer sample_num
  integer seed

  have = 0
  energy = 0.0D+00

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

  energy = energy / real ( sample_num)
end

subroutine cvt_iterate ( dim_num, n, batch, sample, initialize, sample_num, &
  seed, r, it_diff, energy )

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
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer N, the number of points to generate.
!
!    Input, integer BATCH, sets the maximum number of sample points
!    generated at one time.  It is inefficient to generate the sample
!    points 1 at a time, but memory intensive to generate them all
!    at once.  You might set BATCH to min ( SAMPLE_NUM, 10000 ), for instance.
!    BATCH must be at least 1.
!
!    Input, integer SAMPLE, specifies how the sampling is done.
!    -1, 'RANDOM', using FORTRAN RANDOM function;
!     0, 'UNIFORM', using a simple uniform RNG;
!     1, 'HALTON', from a Halton sequence;
!     2, 'GRID', points from a grid;
!     3, 'USER', call "user" routine.
!
!    Input, logical INITIALIZE, is TRUE if the random number generator
!    should be initialized, because this is the first call to it.
!
!    Input, integer SAMPLE_NUM, the number of sample points.
!
!    Input/output, integer SEED, the random number seed.
!
!    Input/output, double precision R(DIM_NUM,N), the Voronoi
!    cell generators.  On output, these have been modified
!
!    Output, double precision IT_DIFF, the L2 norm of the difference
!    between the iterates.
!
!    Output, double precision ENERGY, the discrete "energy", divided
!    by the number of sample points.
!
  implicit none

  integer batch
  integer dim_num
  integer n
  integer sample_num

  integer count(n)
  double precision energy
  integer get
  integer have
  logical initialize
  double precision it_diff
  integer j
  integer nearest(batch)
  double precision r(dim_num,n)
  double precision r2(dim_num,n)
  double precision s(dim_num,batch)
  integer sample
  integer seed
!
!  Take each generator as the first sample point for its region.
!  This can slightly slow the convergence, but it simplifies the
!  algorithm by guaranteeing that no region is completely missed
!  by the sampling.
!
  energy = 0.0D+00
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
    r2(1:dim_num,j) = r2(1:dim_num,j) / real ( count(j))
  end do
!
!  Determine the sum of the distances between the old generators 
!  and the estimated centroids.
!
  it_diff = 0.0D+00
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
  energy = energy / real ( sample_num) 
end

subroutine cvt_sample ( dim_num, n, n_now, sample, initialize, seed, r )

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
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer N, the number of sample points to be generated.
!
!    Input, integer N_NOW, the number of sample points to be 
!    generated on this call.  N_NOW must be at least 1.
!
!    Input, integer SAMPLE, specifies how the sampling is done.
!    -1, 'RANDOM', using FORTRAN RANDOM function;
!     0, 'UNIFORM', using a simple uniform RNG;
!     1, 'HALTON', from a Halton sequence;
!     2, 'GRID', points from a grid;
!     3, 'USER', from the "user" routine.
!
!    Input, logical INITIALIZE, is TRUE if the pseudorandom process should be
!    reinitialized.
!
!    Input/output, integer SEED, the random number seed.
!
!    Output, double precision R(DIM_NUM,N_NOW), the sample points.
!
  implicit none

  integer dim_num
  integer n_now

  double precision exponent
  integer , allocatable, dimension ( : ) :: halton_base
  integer , allocatable, dimension ( : ) :: halton_leap
  integer , allocatable, dimension ( : ) :: halton_seed
  integer halton_step
  integer i
  logical initialize
  integer j
  integer n
  integer ngrid
  integer prime
  double precision r(dim_num,n_now)
  integer rank
  integer rank_max
  integer sample
  integer seed
  integer , allocatable, dimension ( : ) :: tuple

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

    exponent = 1.0D+00 / real ( dim_num)
    ngrid = int ( ( real ( n) )**exponent )
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
      r(1:dim_num,j) = real ( 2 * tuple(1:dim_num) - 1) &
        / real ( 2 * ngrid)
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
end

subroutine data_read ( file_in_name, dim_num, n, r, success )

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
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer N, the number of points.  The program
!    will stop reading data once N values have been read.
!
!    Output, double precision R(DIM_NUM,N), the point coordinates.
!
!    Output, logical SUCCESS, is TRUE if the data was read properly.
!
  implicit none

  integer dim_num
  integer n

  character ( len = * ) file_in_name
  integer file_in_unit
  integer i
  integer ierror
  integer ios
  character ( len = 255 ) line
  double precision r(dim_num,n)
  logical success
  double precision x(dim_num)

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
end

subroutine find_closest ( dim_num, n, sample_num, s, r, nearest )

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
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer N, the number of cell generators.
!
!    Input, integer SAMPLE_NUM, the number of sample points.
!
!    Input, double precision S(DIM_NUM,SAMPLE_NUM), the points to be checked.
!
!    Input, double precision R(DIM_NUM,N), the cell generators.
!
!    Output, integer NEAREST(SAMPLE_NUM), the index of the nearest 
!    cell generators.
!
  implicit none

  integer dim_num
  integer n
  integer sample_num

  double precision dist_sq_min
  double precision dist_sq
  integer jr
  integer js
  integer nearest(sample_num)
  double precision r(dim_num,n)
  double precision s(dim_num,sample_num)

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
end

function halham_leap_check ( dim_num, leap )

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
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer LEAP(DIM_NUM), the leap vector.
!
!    Output, logical, HALHAM_LEAP_CHECK, true if LEAP is legal.
!
  implicit none

  integer dim_num

  logical halham_leap_check
  integer leap(dim_num)

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
end

function halham_n_check ( n )

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
!    Input, integer N, the spatial dimension.
!
!    Output, logical HALHAM_N_CHECK, true if N is legal.
!
  implicit none

  logical halham_n_check
  integer n

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HALHAM_N_CHECK - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    write ( *, '(a,i12)' ) '  N = ', n
    halham_n_check = .false.
  else
    halham_n_check = .true.
  end if
end

function halham_dim_num_check ( dim_num )

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
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Output, logical HALHAM_DIM_NUM_CHECK, true if DIM_NUM is legal.
!
  implicit none

  logical halham_dim_num_check
  integer dim_num

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HALHAM_DIM_NUM_CHECK - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    write ( *, '(a,i12)' ) '  DIM_NUM = ', dim_num
    halham_dim_num_check = .false.
  else
    halham_dim_num_check = .true.
  end if
end

function halham_seed_check ( dim_num, seed )

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
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer SEED(DIM_NUM), the seed vector.
!
!    Output, logical, HALHAM_SEED_CHECK, true if SEED is legal.
!
  implicit none

  integer dim_num

  logical halham_seed_check
  integer seed(dim_num)

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
end

function halham_step_check ( step )

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
!    Input, integer STEP, the index of the subsequence element.
!
!    Output, logical HALHAM_STEP_CHECK, true if STEP is legal.
!
  implicit none

  logical halham_step_check
  integer step

  if ( step < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HALHAM_STEP_CHECK - Fatal error!'
    write ( *, '(a)' ) '  STEP < 0.'
    write ( *, '(a,i12)' ) '  STEP = ', step
    halham_step_check = .false.
  else
    halham_step_check = .true.
  end if
end

function halton_base_check ( dim_num, base )

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
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer BASE(DIM_NUM), the bases.
!
!    Output, logical, HALTON_BASE_CHECK, true if BASE is legal.
!
  implicit none

  integer dim_num

  integer base(dim_num)
  logical halton_base_check

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
end

subroutine s_cap ( s )

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
  implicit none

  character c
  integer i
  integer nchar
  character ( len = * ) s

  nchar = len_trim ( s )

  do i = 1, nchar

    c = s(i:i)
    call ch_cap ( c )
    s(i:i) = c

  end do
end

function s_eqi ( s1, s2 )

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
  implicit none

  character c1
  character c2
  integer i
  integer len1
  integer len2
  integer lenc
  logical s_eqi
  character ( len = * ) s1
  character ( len = * ) s2

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
end

subroutine tuple_next_fast ( m, n, rank, x )

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
!    Input, integer M, the maximum entry in any component.
!    M must be greater than 0.
!
!    Input, integer N, the number of components.
!    N must be greater than 0.
!
!    Input, integer RANK, indicates the rank of the tuple.
!    Typically, 0 <= RANK < N**M.  Values of RANK greater than
!    N**M are legal and meaningful; they are equivalent to the
!    corresponding value mod (N**M).  If RANK < 0, this indicates 
!    that this is the first call for the given values of (M,N).  
!    Initialization is done, and X is set to a dummy value.
!
!    Output, integer X(N), the next tuple, or a dummy value if
!    initialization has just been done.
!
  implicit none

  integer n

  integer , save, allocatable, dimension ( : ) :: base
  integer i
  integer m
  integer rank
  integer x(n)

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
end

subroutine user ( dim_num, n, seed, r )

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
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer N, the number of sample points desired.
!
!    Input/output, integer SEED, the "seed" value.  The meaning and
!    use of this variable is up to the user.
!
!    Output, double precision R(DIM_NUM,N), an array of sample points from
!    the region.
!
  implicit none

  integer dim_num
  integer n

  double precision angle(n)
  integer i
  integer j
  integer k
  double precision , parameter :: pi = 3.141592653589793D+00
  integer seed
  double precision r(dim_num,n)
  double precision radius(n)
!
!  We sample points in the unit circle uniformly.
!
  call random_number ( harvest = angle(1:n) )
  angle(1:n) = 2.0D+00 * pi * angle(1:n)

  call random_number ( harvest = radius(1:n) )
  radius(1:n) = sqrt ( radius(1:n) )

  r(1,1:n) = radius(1:n) * cos ( angle(1:n) )
  r(2,1:n) = radius(1:n) * sin ( angle(1:n) )
end
