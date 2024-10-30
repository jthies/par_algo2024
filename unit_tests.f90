!! This file is part of Fortuno.
! Licensed under the BSD-2-Clause Plus Patent license.
! SPDX-License-Identifier: BSD-2-Clause-Patent

module m_unit_tests
  use m_dotprod
  use fortuno_coarray, only : as_char, test => coa_pure_case_item, context => coa_context,&
      & is_equal, test_item
  implicit none

  private
  public :: get_unit_tests

real, parameter :: invalid_entry = -256.0

!! helper function to test equality (up to floating-point accuracy) of floats of different kind
interface is_float_equal
  module procedure is_float4_equal, is_float8_equal
end interface is_float_equal

! set to .True. to enable debugging output
logical, parameter :: verbose = .True.

contains

  ! Returns the tests from this module
  function get_unit_tests() result(testitems)
    type(test_item), allocatable :: testitems(:)

    testitems = [&
        test("dotprod: all variants", test_dotprod_variants), &
        test("sorting: quicksort", test_quicksort), &
        test("sorting: merge", test_merge), &
        test("sorting: mergeparts", test_mergeparts), &
        test("sorting: mergeparts with odd nparts", test_mergeparts_odd_nparts), &
        test("sorting: sort_coarray", test_sort_coa), &
        test('primes: simple_sieve', test_simple_sieve), &
        test('primes: simple_sieve', test_parallel_sieve) &
    ]

  end function get_unit_tests

!!!!!!!!!!!!!!!!!!!!!!
!! Helper functions !!
!!!!!!!!!!!!!!!!!!!!!!

  logical function is_float4_equal(a, b)
  implicit none
  real(kind=4), intent(in) :: a, b
  is_float4_equal = (abs(a-b)/abs(a+b) < epsilon(a))
  end function is_float4_equal

  logical function is_float8_equal(a, b)
  implicit none
  real(kind=8), intent(in) :: a, b
  is_float8_equal = (abs(a-b)/abs(a+b) < epsilon(a))
  end function is_float8_equal


  logical function is_sorted(x)

  implicit none

  real, dimension(:), intent(in) :: x
  integer :: is_decreasing, i, N

    N = size(X)
    is_decreasing = 0
    do i=1,N-1
      if (x(i+1)<x(i)) then
        is_decreasing = is_decreasing + 1
      end if
    end do


  is_sorted = (is_decreasing==0)

  end function is_sorted

  ! checks both local and global sorting, and if there are
  ! any invalid entries (-256.0), which we use to initialize
  ! the extra space before calling sort_coarray below
  logical function is_sorted_coa(x, nloc)

  implicit none

  real, dimension(:), intent(in) :: x[*]
  integer, intent(in) :: nloc
  
  integer :: me, nproc
  logical :: locally_sorted, globally_sorted, invalid_entries
  integer, save :: co_nloc[*]

  me = this_image()
  nproc = num_images()

  co_nloc = nloc
  sync all

  invalid_entries = any(X(1:nloc) == invalid_entry)

  ! check local sorting
  locally_sorted = is_sorted(X(1:nloc)[me])

  globally_sorted = .true.
  if (me>1) then
      globally_sorted = globally_sorted .and. &
          (X(1)[me] >= X(co_nloc[me-1])[me-1])
  end if

  if (me<nproc) then
      globally_sorted = globally_sorted .and. &
          (X(nloc)[me] <= X(1)[me+1])
  end if

  if (invalid_entries) then
      write(*,*) 'invalid elements in local array (nloc_after may be incorrect)'
  end if
  if (.not. locally_sorted .and. verbose) then
      write(*,*) 'elements not locally sorted on image ',this_image()
      write(*,*) X(:)[me]
  end if
  if (.not. globally_sorted .and. verbose) then
      write(*,*) 'elements not globally sorted'
  end if
  is_sorted_coa = locally_sorted .and. globally_sorted .and. (.not. invalid_entries)
end function is_sorted_coa
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Actual unit tests       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine test_dotprod_variants(ctx)
  
    class(context), intent(inout) :: ctx

    integer, parameter :: N=100
    real(kind=8), save :: x(N)[*], y(N)[*]
    integer :: i

    real(kind=8) :: dot_xy, s_cosum, s_allgather, s_gatherbcast, s_butterfly

    do i=1,N
        x(i) = dble(i+N*(this_image()-1))
        y(i) = dble(this_image())/x(i)
    end do

    ! expected value: On each process we get N*this_image(), so in total it's sum_p p*N = N * (p+1)(p/2)
    dot_xy = N*(num_images()+1)*(dble(num_images())/2.d0)

    s_cosum = dot_cosum(x,y)
    s_allgather = dot_allgather(x,y)
    s_gatherbcast = dot_gatherbcast(x,y)
    s_butterfly = dot_butterfly(x,y)

    ! THEN each rank must contain source rank's value
    call ctx%check(is_float_equal(s_cosum, dot_xy),msg='dot_cosum')
    call ctx%check(is_float_equal(s_allgather, dot_xy), msg='dot_allgather')
    call ctx%check(is_float_equal(s_gatherbcast, dot_xy), msg='dot_gatherbcast')
    call ctx%check(is_float_equal(s_butterfly, dot_xy), msg='dot_butterfly')

  end subroutine test_dotprod_variants


  subroutine test_quicksort(ctx)
    use :: m_sorting, only: quicksort
    class(context), intent(inout) :: ctx

    integer, parameter :: N=121
    real :: x(N)
    integer :: i

    call random_init(repeatable=.True., image_distinct=.True.)
    call random_number(x)

    call quicksort(x, 1, N)

    ! THEN each rank must contain source rank's value
    call ctx%check(is_sorted(X))

  end subroutine test_quicksort

  subroutine test_merge(ctx)

    use :: m_sorting, only: merge
    class(context), intent(inout) :: ctx

    integer, parameter :: a=1, b=6, c=10
    integer :: i
    real, dimension(10) :: x = (/1.0, 3.0, 5.0, 7.0, 9.0, &
                                 2.0, 4.0, 6.0, 8.0, 10.0/)
    real, dimension(10) :: tmp

    call merge(x, a, b, c, tmp)
    do i=1,10
        call ctx%check(is_float_equal(X(i), float(i)), 'Expected x(i)=i after "merge"')
    end do
    call ctx%check(is_sorted(X))

  end subroutine test_merge

  subroutine test_mergeparts(ctx)

    use :: m_sorting, only: mergeparts
    class(context), intent(inout) :: ctx

    real, dimension(10) :: x = (/3.0,6.0,9.0, 1.0,7.0, 4.0,8.0,10.0, 2.0,5.0/)
    ! chosen s.t. X(start(i):start(i+1)-1) is sorted for i=1:4
    integer, dimension(5) :: start = (/1, 4, 6, 9, 11/)
    if (verbose) then
        write(*,'(A,10I3)') 'original x:',int(x)
        write(*,'(A,5I3)') 'start:',start
    end if
    call mergeparts(X, start)
    if (verbose) then
        write(*,'(A,10I3)') 'x after mergeparts:',int(x)
    end if
    call ctx%check(is_sorted(X))

  end subroutine test_mergeparts

  subroutine test_mergeparts_odd_nparts(ctx)

    use :: m_sorting, only: mergeparts
    class(context), intent(inout) :: ctx

    real, dimension(13) :: x = (/3.0,6.0,9.0, 1.0,11.0, 4.0,8.0,12.0, 2.0,7.0, 5.0,10.0,13.0/)
    ! chosen s.t. X(start(i):start(i+1)-1) is sorted for i=1:5
    integer, dimension(6) :: start = (/1, 4, 6, 9, 11,14/)
    if (verbose) then
        write(*,'(A,13I3)') 'original x:',int(x)
        write(*,'(A,6I3)') 'start:',start
    end if
    call mergeparts(X, start)
    if (verbose) then
        write(*,'(A,13I3)') 'x after mergeparts:',int(x)
    end if
    call ctx%check(is_sorted(X))

  end subroutine test_mergeparts_odd_nparts

  ! this is a helper that initializes X linearly and
  ! set nloc_after=size(X) so taht we can start a new test.
  subroutine init_coarray(X, nloc, nloc_after)

  implicit none

  real, dimension(:), codimension[*], intent(out) :: X
  integer, intent(in) :: nloc
  integer, intent(out) :: nloc_after
  integer :: i

  nloc_after = 2*nloc+num_images()

  if (size(X)/=nloc_after) then
      stop 'bad array passed to init_coarray for unit testing.'
  end if

    do i=1,nloc
        X(i) = real((this_image()-1)*nloc+i)
    end do

    X(nloc+1:size(X)) = invalid_entry

  end subroutine init_coarray

  subroutine test_sort_coa(ctx)
    use :: m_sorting, only: sort_coarray
    class(context), intent(inout) :: ctx

    real, dimension(:), allocatable :: x[:]
    real :: x2
    integer :: nloc, nloc_after, bs, me, nproc, i, j, k, idx

    me=this_image()
    nproc=num_images()
    bs=3
    nloc=nproc*nproc*bs
    nloc_after = 2*nloc+nproc

    allocate(X(nloc_after)[*])

    call random_init(repeatable=.True., image_distinct=.True.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Test 1: input is sorted !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call init_coarray(X, nloc, nloc_after)
    call sort_coarray(X, nloc, nloc_after)
    call ctx%check(is_sorted_coa(X, nloc_after), msg='case 1: sort_coarray with globally sorted input')
    call ctx%check(nloc==nloc_after, msg='case 1: nloc should stay the same.')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Test 2: input is only locally permuted !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call init_coarray(X, nloc, nloc_after)
    call random_number(X(1:nloc))
    do i=1,nloc
        X(i) = X(i) + real(me-1)
    end do
    call sort_coarray(X, nloc, nloc_after)
    call ctx%check(is_sorted_coa(X, nloc_after), msg='case 2: sort_coarray with only local sorting required')
    call ctx%check(nloc==nloc_after, msg='case 2: nloc should stay the same.')


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Test 3: input is such tthat complete blocks are exchanged !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call init_coarray(X, nloc, nloc_after)
    j = modulo(me, nproc)+1
    do i=1,nloc
        X(i) = real((j-1)*nloc + i)
    end do
    call sort_coarray(X, nloc, nloc_after)
    call ctx%check(is_sorted_coa(X, nloc_after), msg='case 3: sort_coarray with full blocks swapped')
    call ctx%check(nloc==nloc_after, msg='case 3: nloc should stay the same.')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Test 4: a random array !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call init_coarray(X, nloc, nloc_after)
    call random_number(X(1:nloc))
    X(:) = X(:)*10000.0
    call sort_coarray(X, nloc, nloc_after)
    call ctx%check(is_sorted_coa(X, nloc_after), msg='case 3: sort_coarray with random array')

  end subroutine test_sort_coa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PRIMES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_simple_sieve(ctx)
    use m_primes
    implicit none
    class(context), intent(inout) :: ctx

    integer(idx), parameter :: ntests = 3
    integer(idx), parameter, dimension(ntests) :: Ns= (/100, 1000, 10000/)
    integer(idx), parameter, dimension(ntests) :: expected_nprimes = (/25, 168, 1229/)
    
    integer(idx) :: N, nprimes, nprimes_max
    integer(idx), dimension(Ns(ntests)) :: primes
    logical, dimension(Ns(ntests)) :: sieve
    integer(idx) :: k

    do k=1,ntests
      N = Ns(k)
      nprimes_max = 5*int(dble(N)/floor(log(dble(N))))
      sieve(:) = .true.
      primes(:) = -1
      call simple_sieve(sieve(1:N), primes(1:nprimes_max), nprimes)
      call ctx%check(is_equal(int(nprimes,kind=4), int(expected_nprimes(k),kind=4)))
    end do
  end subroutine test_simple_sieve

  subroutine test_parallel_sieve(ctx)
    use m_primes
    implicit none
    class(context), intent(inout) :: ctx

    integer(idx), parameter :: ntests = 3
    integer(idx), parameter, dimension(ntests) :: Ns= (/9973, 20000, 1000000/)
    integer(idx), parameter, dimension(ntests) :: expected_nprimes = (/1229, 2269, 78498/)
    
    integer(idx) :: N, nprimes
    integer(idx), dimension(:), allocatable :: primes[:]
    integer(idx) :: k

    do k=1,ntests
      N = Ns(k)
      call parallel_sieve(N, primes, nprimes)
      call ctx%check(is_equal(int(nprimes,kind=4), int(expected_nprimes(k),kind=4)))
      deallocate(primes)
    end do
  end subroutine test_parallel_sieve
  
end module m_unit_tests
