program hello

implicit none

! This was an example from one of the online tutorials.
! It does not work -- the output order is left to the process management of the
! cluster (SLURM in this case), and on DelftBlue it is always serialized at the end

write(*,'(A,I2,A,I2)') 'Hello from image',this_image(),'out of',num_images()
sync all
write(*,'(A,I2,A,I2)') 'Goodbye from image',this_image(),'out of',num_images()

end program hello
