program sorting

use m_sorting, only: sort_coarray
use m_benchmarks, only: wtime
use m_sorting, only: sort_coarray
use iso_c_binding, only: tab => c_horizontal_tab

implicit none

integer :: nloc, nloc_after, offset, co_nloc[*]
real(kind=4), allocatable :: x(:)[:], x_orig(:)[:]
real(kind=4), allocatable :: x_gathered(:), x_orig_gathered(:)

character(len=12) :: arg
integer :: me, np, i, j, p
real(kind=8) :: t0, t1

me = this_image()
np = num_images()

if (command_argument_count() > 0) then
    call get_command_argument(1,arg)
    read(arg,*) nloc
else if (me==1) then
    write(*,*) 'Using approximately 100 elements per process by default. To change this, pass N on the command-line.'
    nloc = 100
end if

! make sure nloc is divisible by np^2 (our algorithm needs this)
nloc = ceiling(dble(nloc)/dble(np*np))*np*np

if (me==1) then
    write(*,'(A,I4,A,I8,A)') 'Sort random real array on ',np,' processes with ',nloc, ' elements per process.'
end if

! The sort_coarray function requires us to allocate additional space
! so that the output can have a variable number of elements per process:
nloc_after = 2*nloc+np
allocate(x(nloc_after)[*])
allocate(x_orig(nloc)[*])

! create different random numbers every time the program is run, and
! different ones per processor (image)
call random_init(repeatable=.True., image_distinct=.True.)
call random_number(X_orig)
X_orig = X_orig*100
x(1:nloc) = x_orig(:)

sync all
t0 = wtime()
call sort_coarray(x, nloc, nloc_after)
sync all
t1 = wtime()


! ordered output of X by the root process
if (me==1) then
    write(*,'(A,E16.8,A)') 'Elapsed time: ', t1-t0, ' seconds.'
end if

co_nloc = nloc_after
sync all
if (nloc<=144 .and. me==1) then
    offset=0
    do p=1,me-1
        offset = offset + co_nloc[p]
    end do
    allocate(x_gathered(nloc*np))
    allocate(x_orig_gathered(nloc*np))
    j=1
    do p=1,np
        do i=1,nloc
            x_orig_gathered((p-1)*nloc+i) = x_orig(i)[p]
        end do
        do i=1,co_nloc[p]
            x_gathered(j) = x(i)[p]
            j=j+1
        end do
    end do
end if
sync all
if (nloc<=144 .and. me==1) then
    do i=1,np*nloc
        write(*,'(I4, A, E16.8, A, F16.8)') i, tab, X_orig_gathered(i), tab, X_gathered(i)
    end do
end if

end program sorting
