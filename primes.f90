module m_primes

implicit none

integer, parameter :: idx=8

contains

subroutine simple_sieve(sieve, primes, nprimes)
logical, dimension(:), intent(inout) :: sieve
integer(idx), dimension(:), intent(out) :: primes
integer(idx), intent(out) :: nprimes

integer(idx) :: N, sqN, i, j

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
call collect_primes(sieve, sqN+1, N, primes, nprimes)

end subroutine simple_sieve

! for sieve(a:b), element a:p:b to false.
subroutine filter_range(a, b, p, sieve)

implicit none

integer(idx), intent(in) :: a, b
integer(idx), dimension(:), intent(in) :: p
logical, dimension(:), intent(inout) :: sieve

integer(idx) :: i, offset

do i=1,size(p)
  offset = a-modulo(a, p(i))
  sieve(max(a,p(i)*p(i)):b:p(i)) = .false.
end do

end subroutine filter_range

! given a filtered sieve(a:b), insert any positions still .true. into primes(nprimes+1:nprimes'),
! where nprimes is the value on input, and nprimes' on output.
subroutine collect_primes(sieve, a, b, primes, nprimes)

logical, dimension(:), intent(in) :: sieve
integer(idx), intent(in) :: a, b
integer(idx), dimension(:), intent(inout) :: primes
integer(idx), intent(inout) :: nprimes

integer(idx) :: i

do i=a, b
    if (sieve(i)) then
        nprimes=nprimes+1
        primes(nprimes) = i
    end if
end do

end subroutine collect_primes

end module m_primes

