! Find all prime numbers smaller than an input integer n.
! This is Exercise 1.7 (a-c) on page 71 of the book by Rob Bisseling.
!
! a) When can we stop? Once we have crossed out the multiples up to sqrt(n).
!    Note also that when crossing out multiples of j, we can start from j*j.
! b) see below.
! c) The Prime Number Theorem states that the number of primes p <= N is approximately N/log(N).
!    Consequently, if M=sqrt(N), there are M/log(M) numbers x for which we have to cross out their multiples,
!    starting from x^2 every time. We write:
!
!       C = \sum_{p<N: p is prime} (N-p^2)/p = \psum N/p - \psum p 
!             ^=:\psum
!
program par_sieve

use m_benchmarks, only: wtime
use m_primes

implicit none
#ifdef DEBUG
logical, parameter :: verbose = .true.
#else
logical, parameter :: verbose = .false.
#endif
integer(idx) :: N, nloc, nprimes, nprimes_old, nprimes_max
integer(idx), dimension(:), allocatable :: nprimes_new[:]
integer(idx) :: first_batch
integer :: proc_i
integer(idx) :: P, Q, i, imin, imax
logical, dimension(:), allocatable :: sieve
integer(idx), dimension(:), allocatable :: new_primes[:], primes[:]
character(len=12) :: arg
real(kind=8) :: t0, t1

if (command_argument_count() > 0) then
    call get_command_argument(1,arg)
    read(arg,*) N
else
    if (this_image()==1) write(*,*) 'Using N=100 for testing. To change this, pass N on the command-line.'
    N = 100
end if

! make sure N is divisible by num_images()
N = N + modulo(N, int(num_images(),8))

nprimes_max = 5*int(dble(N)/floor(log(dble(N))))
nloc = N/num_images()
!first_batch = min(100_8, nloc)
first_batch = 2_8

imin = (this_image()-1)*nloc+1
imax = this_image()*nloc

allocate(sieve(imin:imax))
allocate(primes(nprimes_max)[*])
allocate(nprimes_new(num_images())[*])
allocate(new_primes(nprimes_max)[*])

t0 = wtime()

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
        call filter_range(imin, sieve, primes(nprimes_old+1:nprimes))
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

sync all
t1 = wtime()

if (this_image()==1) then
  write(*,'(A,I0, A, I0)') 'Number of primes smaller than ',N, ': ',nprimes
  if (nprimes<=10000) then
    do i=1,nprimes
      write(*,'(I0)') primes(i)
    end do
  end if
  write(*,'(A,G12.4, A)') 'Elapsed time: ',t1-t0,' seconds.'
end if

!do i=imin, imax
!    if (sieve(i)) then
!        write(*,*) 'TROET ',i
!    end if
!end do

deallocate(sieve)
deallocate(primes)

end program par_sieve



