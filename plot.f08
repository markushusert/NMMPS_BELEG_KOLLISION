module plot
  use type_particle,only:Particle
  use setup,only:NP
  use statistics, only:current_timestep, E_kin, E_pot, E_disp
  implicit none

  contains
      subroutine write_plot_file(array_of_particles)
          implicit none
          type(Particle),dimension(:), intent(inout)::array_of_particles
          integer:: counter
          character*8 timestep_string
          character*64 filename, filename_energy
          logical, save :: first_call = .TRUE.
          character*32,parameter::resultdir="results"


          print *,"writing plot output at timestep",current_timestep
          
          write(timestep_string,"(I6.6)") current_timestep
          filename="plot_output_"//trim(timestep_string)//".txt"
          filename_energy="plot_output_energy"//".txt"
          ! output data into a file

          if(first_call) then
            first_call = .FALSE.
            call system('mkdir '//trim(resultdir))
            call chdir(trim(resultdir))
            open(unit=21, file=filename_energy, status = "replace")
          else
            open(unit=21, file=filename_energy,status = "old",Position = "append")
          end if
          write(21,*) current_timestep, E_kin, E_pot, E_disp
          close (unit=21)

          open(unit=20,file=filename, status = 'replace')
          do counter=1,NP
              write(20,*) array_of_particles(counter)%Position
          end do
          close(20)

          call chdir("../")
      end subroutine write_plot_file
end module plot
