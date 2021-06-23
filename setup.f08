module setup
    use type_particle
    use collisions
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
            !set Variables
            !currently just hardcoded some user-inputs
            Np=100!Set number of particles
            call alloclist()

        end subroutine user_input
        subroutine alloclist()
            !allocates dynamically sized arrays
            allocate(array_of_particles(NP))
            allocate(collisionpartners(2,NP**2))
            allocate(tab_list(NP**2))
        end subroutine alloclist
end module setup