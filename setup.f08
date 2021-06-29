module setup
    use type_particle,only: Particle,array_of_particles
    use collision_list_module,only: init_list
    implicit none

    !---------Variables
    integer NP          !number of particles
    integer NB          !number of particles at the bottom

    real RP             !radius of particles
    real Rhop           !density of particles
    real k_stoss        !Stoßzahl
    real g              !Erdbeschleunigung
    real dt             !timestep
    real t_ges          !duration

    contains
        subroutine user_input()
            
            !TODO set variables defined above
            t_ges=1.0
            dt=0.1
            Np=100!Set number of particles (hardcoded for now)
            call alloclist()
            !TODO fill attributes of particles in array_of_particles

        end subroutine user_input
        subroutine alloclist()
            !allocates dynamically sized arrays
            allocate(array_of_particles(NP))
            call init_list()
        end subroutine alloclist
end module setup