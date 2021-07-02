module solid_dynamics
    use type_particle,only: Particle
    use setup,only:g,dt,NP,areawidth
    use statistics,only:E_kin,E_pot,E_disp
    implicit none
    integer i

    contains
    subroutine move(t_inkrement,particle_array)
        real,intent(in)::t_inkrement  !duration to advance by
        Type(Particle),intent(inout)::particle_array(:)

        !TODO update position of all particles in particle-array based on t_inkrement
        !Yvi
        do while (i<=size(particle_array))
          particle_array(i)%Position = particle_array(i)%Position + particle_array(i)%Velocity * t_inkrement !update position
          call teleport_particles(particle_array)
        end do
        !end Yvi

    end subroutine move

    subroutine time_integration(array_of_particles)
        implicit none
        type(Particle),dimension(:), intent(inout)::array_of_particles
        integer:: counter
        real:: absolute_velocity
        real:: deltaE_kin
        do counter=1,NP
          array_of_particles(counter)%velocity(2) = real(array_of_particles(counter)%velocity(2))-g*DT
          absolute_velocity = norm2(array_of_particles(counter)%velocity)
          deltaE_kin = 0.5*array_of_particles(counter)%masse*absolute_velocity**2
          E_kin = E_kin + deltaE_kin
        end do
    end subroutine time_integration

    subroutine teleport_particles(array_of_particles)
        implicit none
        type(Particle),dimension(:), intent(inout)::array_of_particles
        integer:: counter
        do counter=1,NP
          if(real(array_of_particles(counter)%position(1)) < 0) then
          array_of_particles(counter)%position(1) = real(array_of_particles(counter)%position(1)) + areawidth
          endif
          if(real(array_of_particles(counter)%position(1)) > areawidth) then
          array_of_particles(counter)%position(1) = real(array_of_particles(counter)%position(1)) - areawidth
          endif
        end do
    end subroutine teleport_particles

end module solid_dynamics
