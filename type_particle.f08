module type_particle
    !purpose: contains definition of our defined type: particle
    !and stores array of particles
    use compile_constants
    implicit none

    type Particle
        DOUBLE PRECISION radius,masse
        DOUBLE PRECISION,DIMENSION(dim)::Velocity,Position,Force
        logical active,use_soft_sphere
    end type Particle
    
    type(Particle),dimension(:),pointer::array_of_particles
end module type_particle