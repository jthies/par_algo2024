module m_primes

implicit none

! if this variable is .true., debugging output is produced
#ifdef DEBUG
logical, parameter :: verbose = .true.
#else
logical, parameter :: verbose = .false.
#endif

! this determines the size of an integer used in this module.
! default integer(4) will only take us up to N=2 billion, which
! is not enough to justify using a supercomputer :)
integer, parameter :: idx=8

contains

subroutine simple_sieve(sieve, primes, nprimes)
logical, dimension(:), intent(inout) :: sieve
integer(idx), dimension(:), intent(out) :: primes
integer(idx), intent(out) :: nprimes

integer(idx) :: N, sqN, i, j

N = ubound(sieve,1)

nprimes = 0
sqN = int(sqrt(dble(N)))

do i=2,sqN
    if (sieve(i)) then
        nprimes = nprimes+1
        primes(nprimes) = i
        do j=i*i, N, i
            sieve(j) = .false.
        end do
    end if
end do
call collect_primes(sieve(sqN+1:N), sqN+1, N, primes, nprimes)

end subroutine simple_sieve

! Assuming that sieve(1) refers to number 'a',
! for every element p in pprimes, 
! set entries offset:p:end to false.,
! where 'offset' is the smallest multiple of p larger than min(a,p^2)
subroutine simple_filter(a, sieve, primes)

implicit none

integer(idx), intent(in) :: a
logical, dimension(:), intent(inout) :: sieve
integer(idx), dimension(:), intent(in) :: primes

integer(idx) :: i, offset, p, b

b = size(sieve)

do i=1,size(primes)
  p = primes(i)
  offset = a
  do while (modulo(offset,p)>0)
      offset = offset+1
  end do
  offset = max(p*p,offset)-a+1
  sieve(offset:b:p) = .false.
end do

end subroutine simple_filter

! This subroutine does the same as 'simple_filter' above,
! but it uses cache-blocking to achieve much higher performance
subroutine fast_filter(a, sieve, primes)
implicit none
integer(idx), intent(in) :: a
logical, dimension(:), intent(inout) :: sieve
integer(idx), dimension(:), intent(in) :: primes

integer(idx) :: n, nprimes
integer(idx) :: i
integer(idx) :: chunk, imax

chunk = 1000 ! cover the range [a,b] with this step size

n = size(sieve)
nprimes = size(primes)

do i=1,n,chunk
    imax=min(i+chunk,n)
    call simple_filter(a+i-1, sieve(i:imax), primes)
end do

end subroutine fast_filter

! given a filtered sieve(a:b), insert any positions still .true. into primes(nprimes+1:nprimes'),
! where nprimes is the value on input, and nprimes' on output.
subroutine collect_primes(sieve, a, b, primes, nprimes)

logical, dimension(:), intent(in) :: sieve
integer(idx), intent(in) :: a, b
integer(idx), dimension(:), intent(inout) :: primes
integer(idx), intent(inout) :: nprimes

integer(idx) :: i

do i=1, b-a+1
    if (sieve(i)) then
        nprimes=nprimes+1
        primes(nprimes) = a+i-1
    end if
end do

end subroutine collect_primes

! given an integer N, computes all primes <= N
! in parallel. On input, primes(:)[:] should be
! unallocated. On output, it will contain the full list
! of prime numbers on every process. nprimes will contain the
! number of primes found, s.t., primes(nprimes) is the largest
! prime <= N.
subroutine parallel_sieve(N, primes, nprimes)

integer(idx), intent(in) :: N
integer(idx), intent(out) :: nprimes
integer(idx), allocatable :: primes(:)[:]
integer(idx) :: nloc, nprimes_old, nprimes_max
integer(idx), dimension(:), allocatable :: nprimes_new[:]
integer(idx) :: first_batch
integer :: proc_i
integer(idx) :: P, Q, i, imin, imax
logical, dimension(:), allocatable :: sieve
integer(idx), dimension(:), allocatable :: new_primes[:]
character(len=12) :: arg
real(kind=8) :: t0, t1

nprimes_max = 5*int(dble(N)/floor(log(dble(N))))
nloc = N/num_images()

first_batch = min(1000_8, nloc)

imin = (this_image()-1)*nloc+1
imax = this_image()*nloc
if (this_image()==num_images()) imax = this_image()*nloc + modulo(N, int(num_images(),8))

allocate(sieve(imin:imax))
allocate(primes(nprimes_max)[*])
allocate(nprimes_new(num_images())[*])
allocate(new_primes(nprimes_max)[*])

sieve(:) = .true.

nprimes_old = 0

if (this_image() == 1) then
    if (verbose) write(*,'(A,I0,A,I0,A)') 'SUPERSTEP 0: image 1 computes all primes in range [',imin,',',imin+first_batch-1,']'
    call simple_sieve(sieve(imin:imin+first_batch-1), primes, nprimes)
    if (verbose) write(*,'(A,I0,A)') 'Found ',nprimes, ' initial primes.'
end if

if (verbose .and. this_image()==1) write(*,*) 'SUPERSTEP 1: broadcast initial primes.'
call co_broadcast(nprimes, 1)
call co_broadcast(primes(1:nprimes), 1)

do while (.true.)
    ! Filter with all newly found primes in previous step
    P=primes(nprimes_old+1)
    Q=primes(nprimes)
    ! we can skip the filtering completely if none of
    ! the primes would lead to any cross-outs
    !if (Q*Q<=imax) then
    if (.true.) then
        if (verbose) then
            write(*,'(A,I0,A,I0,A,I0,A)') 'SUPERSTEP 2: P',this_image(), ' filter range [',imin,',', imax,']'
            write(*,'(A,I0,A,I0,A)')      '             using primes in range [',P,',',Q,']'
        end if
        call simple_filter(imin, sieve, primes(nprimes_old+1:nprimes))
    end if
    ! Collect new primes if you are the active process.
    ! The active process is the owner of the last prime number identified,
    ! because he is most likely the owner of the next (few).
    ! If we have filtered the whole range with all prime numbers up to P, then
    ! the prime numbers up to P^2 are "laid bare" (only prime numbers
    ! remain in sieve(P+1:P^2).
    P = max(imin, primes(nprimes)+1)
    Q = min(imax,primes(nprimes)*primes(nprimes))
    nprimes_old = nprimes
    nprimes_new(:) = 0
    sync all
    if (verbose .and. P<=Q) write(*,'(A,I0,A,I0,A,I0,A)') 'SUPERSTEP 3: P',this_image(),' collects primes in range [',P+1,',',Q,']'
    call collect_primes(sieve(P:Q), P, Q, new_primes, nprimes_new(this_image()))
    do proc_i=1,num_images()
        nprimes_new(this_image())[proc_i] = nprimes_new(this_image())
    end do
    sync all
    if (verbose .and. this_image()==1) write(*,*) 'nprimes_new: ', nprimes_new(:)
    do proc_i=1,num_images()
        if (nprimes_new(proc_i)>0) then
            primes(nprimes+1:nprimes+nprimes_new(proc_i)) = new_primes(1:nprimes_new(proc_i))[proc_i]
            nprimes = nprimes + nprimes_new(proc_i)
        end if
    end do
    sync all
    
    P = primes(nprimes)
    Q = primes(nprimes_old)
    if (verbose .and. this_image()==1) write(*,'(A,I0,A,I0)') 'nprimes=',nprimes, ' max prime:',P
    if (Q*Q>=N) exit
end do

do while (primes(nprimes)>n)
    nprimes=nprimes-1
end do

deallocate(sieve)
! it is up to the caller to deallocate primes
end subroutine parallel_sieve

end module m_primes

