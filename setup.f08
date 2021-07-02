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
    real areawidth      !width of simulation area
    real radius
    real,parameter::pi=3.14159265
    character use_config
    integer i
    integer n
    integer, allocatable :: seed(:)
    type(Particle) :: bottom_particle, top_particle
    namelist /general/  Np, & !number of particles
                        Nb, & !number of particles at bottom
                        Rp, & !radius of particles
                        k_stoss, & !coefficient of restitution
                        radius, & ! radius of the spheres
                        t_ges, &
                        dt,&
                        Rhop,&
                        g
    real                 :: random_value_x, random_value_y

    contains
        subroutine user_input()
            
            call random_seed(size = n)
           allocate(seed(n))
           call random_seed(get = seed)
            !print *, 'use config.txt for initialization? y/n'
            !read *, use_config
            if (.true. .or. use_config == "y") then
              open (1, file='config.txt', action = 'read')
                read(1,general)
              close (1)
            else
              print *, 'Enter number of particles'
              read(*,*) Np
              print *, 'Enter number of particles at the bottom'
              read(*,*) Nb
              print *, 'Enter radius of particles'
              read(*,*) Rp
              print *, 'Enter coefficient of restitution'
              read(*,*) k_stoss
              print *, 'Enter radius'
              read(*,*) radius
              print *, 'Enter tges'
              read(*,*) t_ges
              print *, 'Enter dt'
              read(*,*) dt
              print *, 'Enter Rhop'
              read(*,*) Rhop
              print *, 'Enter g'
              read(*,*) g
            end if
            ! end Yvi
            areawidth=radius*Nb
            !Yvi
            call alloclist()
            !TODO fill attributes of particles in array_of_particles
            do i=1,Np
              if (i<=Nb) then
                bottom_particle%velocity=[0.0,0.0]
                bottom_particle%position=[i-1+Rp, Rp]
                bottom_particle%radius=radius
                bottom_particle%masse=top_particle%radius**3*pi*4/3*Rhop
                array_of_particles(i) = bottom_particle
                
              else
                call random_number(random_value_x)
                call random_number(random_value_y)
                top_particle%velocity=[0.0,0.0]
                top_particle%position=[Rp*Nb*random_value_x, Rp*Nb*random_value_y]
                top_particle%radius=radius
                top_particle%masse=top_particle%radius**3*pi*4/3*Rhop
                
                array_of_particles(i) = top_particle
              end if
              !insert Überlappung detection
              !end Yvi
            end do

        end subroutine user_input
        subroutine alloclist()
            !allocates dynamically sized arrays
            allocate(array_of_particles(NP))
            call init_list()
        end subroutine alloclist
end module setup
