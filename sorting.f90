module m_sorting

implicit none

contains

recursive subroutine quicksort(x, lo, hi)

implicit none

real, dimension(:), intent(inout) :: x
integer :: lo, hi, i, N

N = size(X)
i = split(X, lo, hi)
if (i-1 > lo) then
    call quicksort(x, lo, i-1);
end if
if (i+1 < hi) then
    call quicksort(x, i+1, hi);
end if
end subroutine quicksort

integer function split(x, lo, hi)
implicit none
real, dimension(:), intent(inout) :: x
integer, intent(in) :: lo, hi
integer :: i, j, piv
real :: val, rnd

! pick piv , with lo ≤ piv ≤ hi;
call random_number(rnd)
piv = lo+anint(rnd*(hi-lo))
val = X(piv)
call swap(x(piv), x(hi))
i = lo
do j = lo, hi-1
    if (x(j) < val) then
        call swap(x(i), x(j))
        i = i+1
    end if
end do
call swap(x(i) , x(hi) )
split = i
end function split

subroutine swap(x, y)
implicit none
real, intent(inout) :: x, y
real :: tmp
tmp=x
x=y
y=tmp
end subroutine swap

! given an array X and integers a <= b <= c s.t.
! X(a:b-1) and X(b:c) are sorted already, sorts
! X(a:c) using merge-sort.
! tmp should be a real array allocated with at least c-a+1 elements,
! it's contents are arbitrary on input and output.
subroutine merge(x, a, b, c, tmp)

real, dimension(:), intent(inout) :: X
integer, intent(in) :: a, b, c

integer :: i, j, k
real, dimension(:), intent(out) :: tmp

if (a >= b .or. b > c) then ! one of the parts is empty
    return
end if

i = a !index for range [a,b-1]
j = b !index for range [b,c]
k = 1 !index for tmp
     
! Compare values as long as no part is empty
do while (i < b .and. j <= c)
    if (x(i) < x(j)) then
        tmp(k) = x(i)
        i = i+1
    else
        tmp(k) = x(j)
        j = j+1
    end if
    k = k+1
end do

! Copy the remaining values
if (i >= b) then
    tmp(k:k+(c-j)) = x(j:c)
end if
if (j > c) then
    tmp(k:k+(b-i)) = x(i:b)
end if
! Copy the values back into x
X(a:c) = tmp(1:c-a+1)

end subroutine merge


! This subroutine sorts the array x for
! the index range start(1) <= i < start(p), p:=size(start).
!
! On input, subranges start(k) <= i < start(k+1) are assumed
! to have been sorted already, for k = 1,2,...,p-1.
!
subroutine mergeparts(x, start)

implicit none

real, dimension(:), intent(inout) :: x
integer, dimension(:), intent(in) :: start

integer :: p, nparts, i, k
integer, allocatable, dimension(:) :: cstart
real, allocatable, dimension(:) :: tmp

p = size(start)

allocate(tmp(start(p)-start(1)+1))

! Initialize current data
nparts= p-1  ! current number of parts
allocate(cstart(p)) ! current start
cstart(:) = start(:)

do while (nparts > 1)
    !write(*,*) 'nparts=',nparts
    !write(*,*) 'cstart=',cstart(1:nparts+1)
    ! Merge pairs of parts
    do i=1,nparts/2
        !write(*,*) 'merge ', cstart(2*i-1), cstart(2*i),cstart(2*i+1)
        call merge(x, cstart(2*i-1), cstart(2*i),cstart(2*i+1)-1, tmp)
    end do
    do i=0,nparts/2-1
        cstart(i+1)=cstart(2*i+1)
    end do
    if (modulo(nparts,2)==1) then
        cstart((nparts+1)/2)=cstart(nparts)
    end if
    nparts = (nparts+1)/2
    cstart(nparts+1) = start(p)
end do

end subroutine mergeparts

!! 
subroutine sort_coarray(X, nloc, nloc_out)

implicit none

integer, intent(in) :: nloc
integer, intent(inout) :: nloc_out
real(kind=4), dimension(nloc_out), codimension[*], intent(inout) :: x
integer :: me, nproc, i, t, count, offset

real, dimension(:), allocatable, save :: sample[:]
real, dimension(:), allocatable :: all_samples, splitval, x_buf
integer, dimension(:), allocatable :: start

! counts(s)[t] is the number of elements process t receives from process s
! offset(s)[t] indicates that process t will obtain those values starting from
! position offset(s)[t] on process s.
integer, allocatable, save :: counts(:)[:], offsets(:)[:]

me = this_image()
nproc = num_images()

! We assume that the local block size is a multiple of nproc^2 for this algorithm,
if (modulo(nloc, nproc*nproc) /= 0) then
    stop 'sort_coarray only works if the local block size is a multiple of nproc^2'
end if

! note: Fortran auto-deallocates arrays at the end of a subroutine,
! so whenever we call this we have to allocate the arrays:
allocate(start(nproc+1))
allocate(splitval(nproc+1))
allocate(all_samples(nproc*nproc))

! co-arrays are different because allocation/deallocation
! overhead is substantial: Therefore, the standard requires
! them to be declared 'save', which means that they persist
! between calls.
if (.not. allocated(sample)) then
    allocate(sample(nproc)[*])
    allocate(counts(nproc)[*])
    allocate(offsets(nproc)[*])
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUPERSTEP 0: sort local block and create samples !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call quicksort(x, 1, nloc)

do i=1,nproc
    sample(i) = x((i-1)*(nloc/nproc)+1)
end do

sync all

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUPERSTEP 1-2: allgather the samples and sort them !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do t=1,nproc
    all_samples((t-1)*nproc+1:t*nproc) = sample(1:nproc)[t]
    start(t) = (t-1)*nproc+1
end do
start(nproc+1) = nproc*nproc+1

! note that we don't need to sync because we used a 'get'

!! Sort the samples: We can exploit that the intervals start[t]:start[t+1]-1
! are already sorted by using mergeparts instead of quicksort.
call mergeparts(all_samples, start)

!! Create splitters: elements s.t. splitval(t) <= X(i) < splitval(t+1) are sent to process t
do t= 1,nproc
    splitval(t) = all_samples((t-1)*nproc+1)
end do

! 'huge' returns the largest possible value for the data type of x
splitval(nproc+1) = huge(x)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUPERSTEP 3: Split the local block and send the resulting parts      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

i=1

do t=1,nproc
    offset= i ! index of first value to be sent
    count = 0 ! number of elements to be sent to process t
    do while (i<=nloc .and. x(i)<splitval(t+1))
        count=count+1
        i=i+1
    end do
    ! this way the other process knows how much memory to allocate
    ! and get from me:
    counts(me)[t] = count
    ! and the offset where it should get from:
    offsets(me)[t] = offset
end do

sync all

! get the values from all other processes and store them in the vararray
count = sum(counts(:))
allocate(x_buf(count))
i=1
do t=1,nproc
    x_buf(i:i+counts(t)-1) = X(offsets(t):offsets(t)+counts(t)-1)[t]
    i=i+counts(t)
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUPERSTEP 4: Sort the received parts !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

sync all

start(1) = 1
do t=1,nproc
    start(t+1) = start(t) + counts(t)
end do

nloc_out = count
x(1:count) = x_buf
call mergeparts(x, start)

end subroutine sort_coarray

end module m_sorting
