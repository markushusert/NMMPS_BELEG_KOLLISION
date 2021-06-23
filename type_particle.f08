module type_particle
    !purpose: contains definition of our defined type: particle
    !and stores array of particles
    use compile_constants
    implicit none

    type Particle
        real radius
        real,DIMENSION(dim)::Velocity,Position
    end type Particle
    
    type(Particle),dimension(:),ALLOCATABLE::array_of_particles
end module type_particle