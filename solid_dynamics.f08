module solid_dynamics
    use type_particle,only: Particle,array_of_particles
    use setup,only:g,dt,NP,areawidth,NB,use_activation
    use statistics,only:E_kin,E_pot,E_disp
    use compile_constants
    implicit none
    integer i
    
    contains
    subroutine move(t_inkrement)
        DOUBLE PRECISION,intent(in)::t_inkrement  !duration to advance by

        !TODO update position of all particles in particle-array based on t_inkrement
        !Yvi
        do i=nb+1,Np
          if (array_of_particles(i)%active) then
            array_of_particles(i)%Position = array_of_particles(i)%Position + array_of_particles(i)%Velocity * t_inkrement !update position
          end if
        end do
        call teleport_particles()
        !end Yvi

    end subroutine move

    subroutine time_integration(ratio)
        implicit none
        real,intent(in)::ratio
        DOUBLE PRECISIOn delta_v(dim)
        integer:: counter
        do counter=NB+1,Np
          delta_v=0.0d0
          if (array_of_particles(counter)%active) then
            delta_v(dim)=-g*DT
            delta_v=delta_v+array_of_particles(counter)%Force/array_of_particles(counter)%masse*dt
            delta_v=delta_v*ratio
            array_of_particles(counter)%velocity=array_of_particles(counter)%velocity+delta_v
          end if
        end do
    end subroutine time_integration

    subroutine update_active_status()
      integer iter_part
      DOUBLE PRECISION E_kin_single
      DOUBLE PRECISION E_pot_single
      DOUBLE PRECISION E_ges
      DOUBLE PRECISION E_lim

      if (.not.use_activation) then
        return
      end if
      
      do iter_part=nb+1,np
        E_lim=g*array_of_particles(iter_part)%masse*&
        (array_of_particles(iter_part)%radius+array_of_particles(1)%radius)
        E_kin_single = (array_of_particles(iter_part)%velocity(1)**2 + array_of_particles(iter_part)%velocity(2)**2)&
         * 0.5 * array_of_particles(iter_part)%masse
        E_pot_single = array_of_particles(iter_part)%masse * g * array_of_particles(iter_part)%position(2)
        E_ges = E_kin_single + E_pot_single
        if (E_ges < E_lim) then
          array_of_particles(iter_part)%active=.false.
          array_of_particles(iter_part)%Velocity=0.0d0
        end if
      end do
    end subroutine update_active_status
    subroutine teleport_particles()
        implicit none
        integer:: counter,iter_dim
        do counter=1,NP
          do iter_dim=1,dim-1
            if(array_of_particles(counter)%position(iter_dim) < 0) then
              array_of_particles(counter)%position(iter_dim) =&
               array_of_particles(counter)%position(iter_dim) + areawidth(iter_dim)
              endif
              if(array_of_particles(counter)%position(iter_dim) > areawidth(iter_dim)) then
              array_of_particles(counter)%position(iter_dim) =&
              array_of_particles(counter)%position(iter_dim) - areawidth(iter_dim)
              endif
          end do
         
        end do
    end subroutine teleport_particles

    subroutine global_energies()
      implicit none
      DOUBLE PRECISION absolute_velocity
      integer counter

      E_kin=0.0d0
      E_pot=0.0d0
      do counter=Nb+1,Np!do not loop over bootom spheres
        absolute_velocity = norm2(array_of_particles(counter)%velocity)
        E_kin = E_kin + 0.5*array_of_particles(counter)%masse*absolute_velocity**2
        E_pot=E_pot+g*array_of_particles(counter)%masse*array_of_particles(counter)%position(dim)
      end do
    end subroutine global_energies

end module solid_dynamics
