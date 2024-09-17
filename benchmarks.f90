module m_benchmarks

implicit none

CONTAINS

! get elapsed (wallclock) time in seconds since an arbitrary time in the past
real(8) function wtime()

implicit none

integer(kind=8) :: count, rate

call system_clock(count, rate)

wtime = dble(count)/dble(rate)
end function wtime

! If the file 'filename' exists and contains a NAMELIST
! with variables bsp_L and bsp_G (real(8), dimension(48)),
! try to read bsp_l L(num_images()) and bsp_G(num_images()) and return them
! as via the function arguments L and G
subroutine read_bsp_params(filename, L, G)

implicit none

character(*), intent(in) :: filename
real(kind=8), intent(out) :: G, L

integer :: npmax
real(kind=8), allocatable, dimension(:) :: bsp_G, bsp_L

integer :: fu

namelist  /bsp_npmax/ npmax
namelist /bsp_params/ bsp_L, bsp_G

open(newunit=fu, file=filename,action='read', status='old')
read(unit=fu, nml=bsp_npmax)
allocate(bsp_G(npmax), bsp_L(npmax))
read(unit=fu, nml=bsp_params)
close(fu)

G=0.d0
L=0.d0

if (num_images() <= npmax) then
  G = bsp_G(num_images())
  L = bsp_L(num_images())
end if

if (G==0.d0 .or. L==0.d0) then
    if (this_image()==1) then
        write(*,'(A,I3,A,A)') 'Values G and/or L for ',num_images(),' processes not available in file ',filename
    end if
end if

end subroutine read_bsp_params

! measures bandwidth bw [GB/s\] and flop rate R [GFlop/s] of an axpy operation
! (y=a*x+y, two loads, one store) of given size N, averaged over ntimes instances.
! This is a sequential benchmark executed by the calling process.
subroutine bench_axpy(N, ntimes, bw, R)

implicit none

integer(kind=8), intent(in) :: N
integer(kind=8), intent(in)         :: ntimes
real(kind=8), intent(out) :: bw, R

real(kind=8), dimension(N)  :: x, y
real(kind=8) :: a, t0, t1
integer(kind=8) :: i, k

t0 = wtime()

do k=1,ntimes
!$omp simd
    do i=1,N
        y(i) = a*x(i) + y(i)
    end do
end do

t1 = wtime()

bw = (24.D0*dble(N)*dble(ntimes)) / (t1-t0) * 1.0e-9
R  = (2.D0 *dble(N)*dble(ntimes)) / (t1-t0) * 1.0e-9

end subroutine bench_axpy

! Measure the average wallclock-time for an "h-relation" (as defined in the BSP model by Rob Bisseling).
! For a given number of processes P and integer h>=0, an h-relation is a communication pattern where
! each process j sends h distict one-word (8 bytes here) messages to up to min(h,P) processes
real(8) function bench_hrel(h,ntimes)

implicit none

integer(kind=8), intent(in) :: h, ntimes

integer(kind=8), dimension(h) :: sendbuf
! note: we let remote processors p each write to a contiguous
! memory space recvbuf(:,p) with leading dimension a multiple of the
! cache-line length. That way, no false sharing can occur (that is,
! no synchronization because cache-line elements are updated by multiple
! processors).
integer(kind=8), allocatable :: recvbuf(:,:)[:]

integer :: i, j, np, me
integer, dimension(h) :: dest, destidx
real(kind=8) :: t0, t1

np = num_images()
me = this_image()

! let every process write to a memory location aligned to 10 cache-line (80 element) blocks
! to avoid unwanted synchronization due to cache-line sharing. Note that this
! is different from the 'interleaved' writing in the Bisseling implementation of bspbench,
! where adjacent processes write to adjancent cache-line elements and the performnace is much
! worse on modern CPUs.
allocate(recvbuf(h+modulo(h,80_8), np)[*])

! Initialize communication pattern
do i=1,h
    ! arbitrary data to send
    sendbuf(i)= i;
    ! note: mod(x,0) is not defined in Fortran and leads to an FPE
    if (np==1) then
        dest(i) = 1
    else
        ! destination processor is one of the p-1 others
        dest(i) = modulo(me+1+modulo(i,(np-1)),np)+1
    end if
end do

! Measure time of ntimes h-relations
sync all
t0 = wtime()
do j=1,ntimes
    do i=1,h
        recvbuf(i,me)[dest(i)] = sendbuf(i)
    end do
    sync all
end do
t1 = wtime()

bench_hrel = (t1-t0)/dble(ntimes)

end function bench_hrel

end module m_benchmarks
