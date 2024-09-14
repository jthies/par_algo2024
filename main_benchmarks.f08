
! helper function to solve an over-determined linear system with n rows and 2 columns in the least-squares sense.
! The algorithm is simple and not numerically robust, it is only used by our benchmark driver
! to put a straight line through some data points and should not be used for other numerical
! calculations
subroutine least_squares(n, A, x, b)

implicit none

integer(kind=8), intent(in) :: n
real(kind=8), dimension(n, 2), intent(in) :: A
real(kind=8), dimension(2), intent(out) :: x
real(kind=8), dimension(n), intent(in) :: b

real(kind=8), dimension(2, 2) :: G
real(kind=8) :: piv

G = matmul(transpose(A), A)
x = matmul(transpose(A), b)

! Solve (A'A)x = A'b in-place
piv = G(2,1)/G(1,1)
G(2,:) = G(2,:) - piv * G(1,:)
X(2) = x(2) - piv*x(1)
x(2) = x(2) / G(2,2)
x(1) = (x(1) - G(1,2)*x(2))/G(1,1)

end subroutine least_squares

!! This program runs a series of classical benchmarks:
!! - axpy for an increasing vector length,
!!   reporting memory bandwidth and flop rate. These parameters
!!   define the maximum attainable performance according to the roofline model.
!! - bsp-style communication benchmark to determine bandwidth and latency parameters.
!!
!! The parallel implementation is based on coarrays.
!!
program main_benchmarks

use iso_c_binding, only: tab => c_horizontal_tab
use m_benchmarks

implicit none

! start with N=Nmin (i.e., 1kB here)
integer(kind=8), parameter :: Nmin = 128
! double nsteps times (up to 1GB here)
integer, parameter :: nsteps = 20
integer(kind=8) :: N, ntimes
integer :: i, j
real(8) :: bw, R
integer :: nproc

! measure h-relations for 1 <= h <= maxh
integer(kind=8), parameter :: maxh = 2048
integer(kind=8) :: h
real(kind=8), dimension(0:maxh) :: t_hrel
real(kind=8), dimension(0:maxh, 2) :: A_hrel
! will contain BSP model parameters L (latency) and G ("gap" or inverse bandwidth)
real(kind=8), dimension(2) :: p_hrel

nproc = num_images()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! AXPY benchmarks for floating-point performance and memory bandwidth  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

N = Nmin

if (this_image()==1) then
    write(*,'(A6, A1, A8,A1,A16,A1,A16,A1,A16)') '%nproc', tab, 'N', tab, 'mem [GB]', tab, 'bw [GB/s]', tab, 'R [GFlop/s]'
    write(*,'(A6,I4,A8)') 'datBR{',nproc,'} = [...'
end if

do i=1,nsteps
    ntimes = max(50_8, (Nmin*10000_8)/N)
    call bench_axpy(N, ntimes, bw, R)
    sync all
    call co_sum(bw, result_image=1)
    call co_sum(R, result_image=1)

    if (this_image()==1) then
        write(*,'(I5, A1, I8, A1, G16.8,A1,G16.8,A1,G16.8)') nproc, tab, N, tab, dble(N)*16.0e-9, tab, bw, tab, R
        !write(*,*) N, dble(N)*16.0e-9, bw, R
    end if
    N = N * 2
end do

if (this_image()==1) then
    write(*,*) '];'
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! BSP-benchmarks for communication model (latency and bandwidth)       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ntimes=50

do h=0,maxh
    A_hrel(h, 1) = 1.0_8
    A_hrel(h, 2) = h
    !if (this_image()==1) then
    !   write(*,'(A,I4,A)') 'benchmark ',h,'-relations...'
    !end if
    t_hrel(h) = bench_hrel(h, ntimes)
end do

call least_squares(maxh+1, A_hrel, p_hrel, t_hrel)

if (this_image()==1) then

    write(*,'(A6, A1, A4, A1, A5)') '%nproc', tab, 'h', tab, 't [s]'
    write(*,'(A6,I4,A8)') 'dat_h{',nproc,'} = [...'

    do h=0, maxh
        write(*,'(I6, A1, I4, A1, G16.8)') nproc, tab, h, tab, t_hrel(h)
    end do
    write(*,*) '];'
    write(*,*) '%--------------------------------------------------------------------'
    write(*,'(A2,I4,A4,G8.2,A9,G8.2,A6)') 'L(',num_images(),') = ',p_hrel(1),'; %[s]'
    write(*,'(A2,I4,A4,G8.2,A9,G8.2,A6)') 'G(',num_images(),') = ',p_hrel(2),'; %[s]'
end if
end program main_benchmarks

