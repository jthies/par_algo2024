module m_dotprod

implicit none

contains

! compute the dot product of two distributed vectors x and y
! using the Fortran intrinsics dot_product (serial) and
! co_sum.
real(8) function dot_cosum(x, y)

implicit none

real(kind=8), dimension(:), intent(in) :: x, y

dot_cosum = dot_product(x, y)
call co_sum(dot_cosum)

end function dot_cosum

! compute the dot product in the way described in the book
! by Rob Bisseling: everyone does their local work (using again
! the dot_product intrinsic), and then shares their local contribution
! with everyone else in a single communication superstep ('allgather' pattern)
real(8) function dot_allgather(x, y)

implicit none
real(kind=8), dimension(:), intent(in) :: x, y
real(kind=8), save, allocatable :: local_dots(:)[:]
integer :: me, np, i

np = num_images()
me = this_image()

if ((.not. allocated(local_dots)) .or. (size(local_dots) /= np)) then
    allocate(local_dots(np)[*])
end if

local_dots(me) = dot_product(x, y)

do i=1,np
    local_dots(me)[i] = local_dots(me)
end do

sync all

dot_allgather = sum(local_dots)

end function dot_allgather

! This is your TODO: Write a variant that collects the local contributions
! only on image 1 ("gather" pattern), does the final summation there and 
! sends the final result to everyone else ("broadcast" pattern).
real(8) function dot_gatherbcast(x, y)

implicit none
real(kind=8), dimension(:), intent(in) :: x, y

! Your job: You can either define a coarray of dimension np and use put operations
! (as in the allgather variant), or you can define only a scalar coarray variable
! and an array on image 1, and then use 'get' operations.

dot_gatherbcast = 0.d0

end function dot_gatherbcast

! Simple implementation of the butterfly reduction algorithm that
! will only work for num_images() a power of 2.
real(8) function dot_butterfly(x, y)

implicit none
real(kind=8), dimension(:), intent(in) :: x, y
! note: the standard requires coarrays inside functions to be
! persistent in one way or the other so that one doesn't accidently
! recreate the coarray many times (which causes an implicit sync)
real(kind=8), save :: s[*], t[*]
integer :: me, np
integer :: nlev, l, other

np = num_images()
me = this_image()

if (np /= 2**(log(dble(np))/log(2.d0))) then
    stop 'The butterfly_dot variant only works for num_images() a power of 2'
end if

s = dot_product(x, y)

!
! Butterfly reduction: E.g., for np=4 we get
!
!  s1  s2  s3  s4
!    \/      \/
!    /\      /\
! s12 s12 s34  s34
!    \___/
!    /   \ (and similar for p2 <-> P4
nlev = log(dble(np))/log(2.d0)

do l=0, nlev-1
    other = ieor(me-1, 2**l)+1
    t[other] = s
    sync all
    !write(*,'(A,I1,A,I1,A,I1,A,G8.2,A,G8.2)') 'level ',l,' ',me,'->',other, ', s=',s, ', t=',t
    s=s+t
    sync all
end do

dot_butterfly = s

end function dot_butterfly

end module m_dotprod
