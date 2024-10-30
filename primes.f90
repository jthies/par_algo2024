module m_primes

implicit none

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
subroutine filter_range(a, sieve, primes)

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

end subroutine filter_range

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

end module m_primes

