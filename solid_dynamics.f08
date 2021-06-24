module solid_dynamics
    use type_particle,only: Particle
    use setup,only:g,dt
    implicit none

    contains
    subroutine move(t_inkrement,particle_array)
        real,intent(in)::t_inkrement  !duration to advance by
        Type(Particle),intent(inout)::particle_array(:)
        
        !TODO update position of all particles in particle-array based on t_inkrement

    end subroutine move

    subroutine time_integration(particle_array)
        Type(Particle),intent(inout)::particle_array(:)
        
        !TODO update velocity of all particles based on timestep dt and g

    end subroutine time_integration

    subroutine teleport_particles(particle_array)
        Type(Particle),intent(inout)::particle_array(:)
        
        !TODO check for all particles if the have left the area of our simulation
        !if so then "teleport" to the other side

    end subroutine teleport_particles

end module solid_dynamics