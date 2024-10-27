module m_primes

implicit none

contains

subroutine simple_sieve(sieve, primes, nprimes)
logical, dimension(:), intent(inout) :: sieve
integer, dimension(:), intent(out) :: primes
integer, intent(out) :: nprimes

integer :: N, sqN, i, j

N = size(sieve)

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
do i=sqN+1,N
    if (sieve(i)) then
        nprimes = nprimes+1
        primes(nprimes) = i
    end if
end do

end subroutine simple_sieve

end module m_primes

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

integer :: N, nprimes, nprimes_max, i
logical, dimension(:), allocatable :: sieve
integer, dimension(:), allocatable :: primes
character(len=12) :: arg
real(kind=8) :: t0, t1

if (command_argument_count() > 0) then
    call get_command_argument(1,arg)
    read(arg,*) N
else
    write(*,*) 'Using N=100 for testing. To change this, pass N on the command-line.'
    N = 100
end if

nprimes_max = 5*int(dble(N)/floor(log(dble(N))))

allocate(sieve(N))
allocate(primes(nprimes_max))

t0 = wtime()

sieve(:) = .true.

call simple_sieve(sieve, primes, nprimes)

write(*,'(A,I0, A, I0)') 'Number of primes smaller than ',N, ': ',nprimes
if (nprimes<=100) then
  do i=1,nprimes
    write(*,'(I0)') primes(i)
  end do
end if

t1 = wtime()
write(*,'(A,G12.4, A)') 'Elapsed time: ',t1-t0,' seconds.'

deallocate(sieve)
deallocate(primes)

end program seq_sieve



