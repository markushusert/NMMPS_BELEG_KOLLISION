module settings
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
    
end module settings