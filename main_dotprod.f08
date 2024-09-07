!!!!!!!!!!!!!!!!!!!!!!
!! HELPER FUNCTIONS !!
!!!!!!!!!!!!!!!!!!!!!!

! wtime(): Get elapsed (wallclock) time in seconds since an arbitrary time in the past
! Typical usage would be:
!
! real(kind=8) :: t0, t1
! 
! t0=wtime()
! ... do some work ...
! t1=wtime()
! write(*,*) 'Elapsed time: ',t1-t0, ' seconds.'
!
real(kind=8) function wtime()

implicit none

integer(kind=8) :: count, rate

call system_clock(count, rate)

wtime = dble(count)/dble(rate)
end function wtime

! given a function like dot_cosum and vectors x, y, run ntimes
! times the operation dot_cosum(x,y) and return the average run time.
real(8) function timeit(dot_foo, x, y, ntimes)

implicit none

interface
  real(8) function wtime()
  end function wtime
  real(8) function dot_foo(x, y)
    real(8), dimension(:), intent(in) :: x, y
  end function dot_foo
end interface

real(kind=8), dimension(:), intent(in) :: x, y
integer(kind=8), intent(in) :: ntimes

integer(kind=8) :: i
real(kind=8) :: t0, t1, s

t0 = wtime()
do i=1,ntimes
  s = dot_foo(x, y)
end do
t1 = wtime()

! for checking the corretness of the implementation,
! print the result of the dot product
if (this_image()==1) then
    write(*,*) 'Result of dotprod variant: ',s
end if

timeit = (t1-t0)/dble(ntimes)

end function timeit


!! Program for testing and benchmarking the dot product variants in dotprod.f08.
!!
!! Driver that uses a fixed N obtained from the command-line,
!! allocates a **total** of N elements across all processes,
!! and prints one line of output for the different variants.
program main_dotprod

use iso_c_binding, only: tab => c_horizontal_tab
use m_dotprod

implicit none

! In order to use the functions above, they have to be declared here.
! To circumvent this, we could put them in a module and add a "use"
! statement, like for the dotprod variants.
interface
real(kind=8) function wtime()
end function wtime
real(kind=8) function timeit(dot_foo,x, y, ntimes)

! one of our input arguments is a function, hence again - an interface declaration
! so the compiler knows it's argument and return types.
interface
  real(8) function dot_foo(x, y)
    real(8), dimension(:), intent(in) :: x, y
  end function dot_foo
end interface
real(kind=8), dimension(:), intent(in) :: x, y
integer(kind=8), intent(in) :: ntimes
end function timeit
end interface

integer(kind=8) :: N, ntimes, k
integer :: np, me, i, nloc

real(kind=8) :: s, t0, t1, rel_err, my_rel_err
real(kind=8), dimension(:), allocatable :: x, y

integer :: num_args
character(len=12) :: arg

! timing results for the three variants
real(kind=8) :: t(4)

np = num_images()
me = this_image()

num_args = command_argument_count()
N = 1000000
if (num_args>0) then
    call get_command_argument(1,arg)
    read(arg,*) N
end if

ntimes = max(50_8, 10000000_8/N)
nloc = N/np
if (me <= modulo(N,int(np,kind=8))) then
    nloc = nloc + 1
end if

allocate(x(nloc), y(nloc))

do k=1,nloc
    x(k) = 1.d0*dble(k)
    y(k) = 1.d0/dble(k)
end do

! Run benchmarks and print results for different dot-prodcuct variants
t(1) = timeit(dot_cosum, x, y, ntimes)
t(2) = timeit(dot_allgather, x, y, ntimes)
t(3) = timeit(dot_gatherbcast, x, y, ntimes)
t(4) = timeit(dot_butterfly, x, y, ntimes)

if (this_image()==1) then
    write(*,'(A6, A1, A8,A1,A16,A1,A16,A1,A16,A1,A16)') '%nproc', tab, 'N', tab, 't_cosum', tab, 't_allgather', tab, 't_gattherbcast', tab, 't_butterfly'
    write(*, '(I3, A1, I8, A1, E16.8, A1, E16.8, A1, E16.8, A1, E16.8)') np, tab, N, tab, t(1), tab, t(2), tab, t(3),tab, t(4)
end if

deallocate(x, y)

end program main_dotprod
