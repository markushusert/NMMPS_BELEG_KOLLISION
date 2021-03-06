module plot
  use type_particle,only:Particle,array_of_particles
  use setup,only:NP,output_inkr,mass_mean
  use solid_dynamics,only:global_energies
  use statistics, only:current_timestep, E_kin,E_num, E_pot, E_disp,counted_collisions
  implicit none

  contains
      subroutine write_plot_file()
          implicit none
          integer:: counter
          integer massflag
          character*8 timestep_string
          character*64 filename, filename_energy
          logical, save :: first_call = .TRUE.
          character*32,parameter::resultdir="results"

          call global_energies()

          print *,"writing plot output at timestep",current_timestep
          
          write(timestep_string,"(I6.6)") current_timestep
          filename="plot_output_"//trim(timestep_string)//".txt"
          filename_energy="plot_output_energy"//".txt"
          ! output data into a file

          if(first_call) then
            first_call = .FALSE.
            call system('mkdir '//trim(resultdir))
            call system('cp starting_positions.txt '//trim(resultdir))
            call system('cp config.txt '//trim(resultdir))
            call chdir(trim(resultdir))
            open(unit=21, file=filename_energy, status = "replace")
          else
            call chdir(trim(resultdir))
            open(unit=21, file=filename_energy,status = "old",Position = "append")
          end if
          write(21,*) current_timestep, E_kin, E_pot, E_disp ,E_num,E_disp+E_num,&
          E_kin+ E_pot+ E_disp+E_num,counted_collisions
          close (unit=21)

          if (mod(current_timestep,output_inkr).eq.0) then
            open(unit=20,file=filename, status = 'replace')
            do counter=1,NP
                if (array_of_particles(counter)%masse.gt.mass_mean) then
                  massflag=1
                else
                  massflag=0
                end if
                write(20,*) array_of_particles(counter)%Position,array_of_particles(counter)%radius,&
                array_of_particles(counter)%Velocity,massflag
            end do
            close(20)
          end if

          call chdir("../")
      end subroutine write_plot_file
      
      
end module plot
