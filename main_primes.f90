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
integer(idx) :: N, nprimes
integer(idx), dimension(:), allocatable :: primes[:]
character(len=12) :: arg
real(kind=8) :: t0, t1
integer(idx) :: i

if (command_argument_count() > 0) then
    call get_command_argument(1,arg)
    read(arg,*) N
else
    if (this_image()==1) write(*,*) 'Using N=100 for testing. To change this, pass N on the command-line.'
    N = 100
end if

t0 = wtime()
call parallel_sieve(N, primes, nprimes)
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

deallocate(primes)

end program par_sieve



