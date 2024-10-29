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
program seq_sieve

use m_benchmarks, only: wtime
use m_primes

implicit none

logical, parameter :: verbose = .true.
integer(idx) :: N, nloc, nprimes, nprimes_old, nprimes_max
integer(idx) :: batch_size, active_proc
integer(idx) :: P, Q, i, imin, imax
logical, dimension(:), allocatable :: sieve
integer(idx), dimension(:), allocatable :: primes[:]
character(len=12) :: arg
real(kind=8) :: t0, t1

if (command_argument_count() > 0) then
    call get_command_argument(1,arg)
    read(arg,*) N
else
    write(*,*) 'Using N=100 for testing. To change this, pass N on the command-line.'
    N = 100
end if

! make sure N is divisible by num_images()
N = N + modulo(N, int(num_images(),8))

nprimes_max = 5*int(dble(N)/floor(log(dble(N))))
nloc = N/num_images()
batch_size = min(100_8, nloc)

imin = (this_image()-1)*nloc+1
imax = this_image()*nloc

allocate(sieve(imin:imax))
allocate(primes(nprimes_max)[*])

t0 = wtime()

sieve(:) = .true.

active_proc = 1
nprimes_old = 0

if (this_image() == active_proc) then
    if (verbose) write(*,'(A,I0,A,I0,A,I0,A)') 'SUPERSTEP 0: image ',active_proc,' computes all primes in range [',imin,',',imin+batch_size-1,']'
    call simple_sieve(sieve(imin:imin+batch_size-1), primes, nprimes)
    if (verbose) write(*,'(A,I0,A)') 'Found ',nprimes, ' initial primes.'
end if

if (verbose .and. this_image()==1) write(*,*) 'SUPERSTEP 1: broadcast initial primes.'
call co_broadcast(nprimes, active_proc)
call co_broadcast(primes(1:nprimes), active_proc)

do while (.true.)
    ! Filter with all newly found primes in previous step
    if (verbose) then
        write(*,'(A,I0,A,I0,A,I0,A)') 'SUPERSTEP 2: P',this_image(), ' filter range [',imin,',', imax,']'
        write(*,'(A,I0,A,I0,A)') '        using primes in range [',primes(nprimes_old+1),',',primes(nprimes),']'
    end if
    if (verbose) then
        write(*,*) 'current nprimes: ',nprimes
        write(*,*) 'current nprimes_old: ',nprimes_old
        write(*,*) 'current primes: ',primes(1:nprimes)
    end if
    call filter_range(imin, sieve, primes(nprimes_old+1:nprimes))

    ! Collect new primes if you are the active process.
    ! The active process is the owner of the last prime number identified,
    ! because he is most likely the owner of the next (few).
    ! If we have filtered the whole range with all prime numbers up to P, then
    ! the prime numbers up to P^2 are "laid bare" (only prime numbers
    ! remain in sieve(P+1:P^2).
    P = primes(nprimes)
    active_proc = (P+1)/nloc+1
    Q = min(P*P, active_proc*nloc)
    nprimes_old = nprimes
    if (this_image()==active_proc) then
        if (verbose) write(*,'(A,I0,A,I0,A,I0,A)') 'SUPERSTEP 3: P',this_image(),' collects primes in range [',P+1,',',Q,']'
        call collect_primes(sieve(P+1:Q), P+1, Q, primes, nprimes)
    end if
    call co_broadcast(nprimes, active_proc)
    call co_broadcast(primes(nprimes_old+1:nprimes), active_proc)
    P = primes(nprimes)
    if (verbose .and. this_image()==0) write(*,'(A,I0,A,I0)') 'nprimes=',nprimes, ' max prime:',P
    if (Q>=N) exit
end do

do while (primes(nprimes)>n)
    nprimes=nprimes-1
end do

sync all
t1 = wtime()

if (this_image()==1) then
  write(*,'(A,I0, A, I0)') 'Number of primes smaller than ',N, ': ',nprimes
  if (nprimes<=100) then
    do i=1,nprimes
      write(*,'(I0)') primes(i)
    end do
  end if
  write(*,'(A,G12.4, A)') 'Elapsed time: ',t1-t0,' seconds.'
end if

deallocate(sieve)
deallocate(primes)

end program seq_sieve



