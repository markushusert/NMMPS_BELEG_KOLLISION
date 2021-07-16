module setup
    use type_particle,only: Particle,array_of_particles
    use collision_list_module,only: init_list
    use compile_constants
    Use, intrinsic :: iso_fortran_env, Only : iostat_end
    USE ieee_arithmetic
    implicit none

    !---------Variables
    integer NP          !number of particles
    integer NB          !number of particles at the bottom
    integer NB1d

    real radius             !radius of particles
    real Rhop           !density of particles
    DOUBLE PRECISION E_Lim
    DOUBLE PRECISION k_stoss        !Stoßzahl
    DOUBLE PRECISION g              !Erdbeschleunigung
    DOUBLE PRECISION dt             !timestep
    DOUBLE PRECISION t_ges          !duration
    DOUBLE PRECISION rho_layer(10)   !rho for eeach layer
    DOUBLE PRECISION rho_to_use
    DOUBLE precision mass_mean
    logical use_activation
    integer mod_collision_calc
    real areawidth(dim-1)      !width of simulation area
    real v_init
    real,parameter::pi=3.14159265
    real,parameter::wiggle_room=1.0d0!how many radii do the spheres have to be placed in
    real,parameter::rel_intervall=wiggle_room+2.0d0
    character use_config
    integer i
    integer iter_1,iter_2
    integer output_inkr
    integer iter_layer
    integer n,n_sphere_in_layer,n_layer
    real height
    integer, allocatable :: seed(:)
    logical read_start_positions,turn_of_hard_sphere
    integer use_soft_sphere
    real length_intervall,start_intervall
    type(Particle) :: bottom_particle, top_particle
    namelist /general/  n_sphere_in_layer, & !number of particles
                        Nb, & !number of particles at bottom
                        radius, & !radius of particles
                        k_stoss, & !coefficient of restitution
                        t_ges, &
                        dt,&
                        Rhop,&
                        g,&
                        height,&
                        n_layer, &
                        v_init,&
                        read_start_positions,&
                        output_inkr,&
                        rho_layer,&
                        use_activation,&
                        mod_collision_calc,&
                        use_soft_sphere,&
                        turn_of_hard_sphere
    real                 :: random_value_x, random_value_y
    real infi
    
    
    

    contains
        subroutine user_input()
          logical file_exists
          integer iter_dim
          integer iter_dim_temp
          integer compare_offset
            
          mass_mean=0.0d0
            call random_seed(size = n)
           allocate(seed(n))
           call random_seed(get = seed)
            !print *, 'use config.txt for initialization? y/n'
            !read *, use_config
            if (.true. .or. use_config == "y") then
              output_inkr=1!default value
              use_activation=.false. !default
              turn_of_hard_sphere=.false.
              use_soft_sphere=0
              mod_collision_calc=0 !default
              open (1, file='config.txt', action = 'read')
                read(1,general)
              close (1)
              NB1d=NB
              NB=NB1d**(dim-1)
            else
              print *, 'Enter number of particles in a layer'
              read(*,*) n_sphere_in_layer
              print *, 'Enter number of particles at the bottom'
              read(*,*) Nb
              print *, 'Enter radius of particles'
              read(*,*) radius
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
            NP=n_layer*n_sphere_in_layer**(dim-1)+NB
            areawidth=radius*Nb1d*2
            if (n_layer*2.gt.NB1d) then
              print *,"there cannot be more than",NB1d/2,"spheres in a layer"
              call exit()
            end if
            
            !TODO fill attributes of particles in array_of_particles
            IF (ieee_support_inf(infi)) THEN
              infi = ieee_value(infi,  ieee_positive_inf)
            END IF
            !check if start_position file exists
            INQUIRE(FILE=init_filename, EXIST=file_exists)
            if (.not.read_start_positions.or. .not.file_exists) then 
              !Yvi
              call alloclist()
              print *,"randomly generating particle positions"
              i=0
              do iter_2=0,NB1d-1
                if (dim.eq.3 .or.iter_2.eq.0) then
                  do iter_1=0,NB1d-1
                    i=i+1
                    bottom_particle%velocity=0.0
                    bottom_particle%position=0.0d0
                    bottom_particle%position(1)=bottom_particle%position(1)+iter_1*2*radius
                    bottom_particle%position(2)=bottom_particle%position(2)+iter_2*2*radius
                    bottom_particle%radius=radius
                    bottom_particle%masse=infi
                    bottom_particle%active=.false.
                    bottom_particle%use_soft_sphere=.true. !set to true since both partners need to be true for soft sphere
                    array_of_particles(i) = bottom_particle
                  end do
                end if
              end do   

              do iter_layer=0,n_layer-1
                if (rho_layer(iter_layer+1).gt.0.0d0) then
                  rho_to_use=rho_layer(iter_layer+1)
                else
                  rho_to_use=Rhop
                end if
                do iter_2=0,n_sphere_in_layer-1
                  if (dim.eq.3 .or.iter_2.eq.0) then
                    do iter_1=0,n_sphere_in_layer-1
                      i=i+1
                      top_particle%velocity=0.0d0
                      top_particle%velocity(dim)=v_init
                      top_particle%position=0.0d0
                      top_particle%position(dim)=height+iter_layer*radius*3
                      top_particle%active=.true.
                      do iter_dim=1,dim-1
                        call random_number(random_value_x)
                        if (iter_dim.eq.1) then
                          iter_dim_temp=iter_1
                          compare_offset=1
                        else
                          iter_dim_temp=iter_2
                          compare_offset=n_sphere_in_layer
                        end if
                        length_intervall=areawidth(iter_dim)/n_sphere_in_layer-2*radius
                        start_intervall=(length_intervall+2*radius)*iter_dim_temp
                        
                        !each particle gets an intervall of NB/n_sphere_in_layer radii in x-direction where it can be placed wherever
                        top_particle%position(iter_dim)=start_intervall+length_intervall*random_value_x

                        bottom_particle=array_of_particles(i-compare_offset)
                        if (abs(top_particle%Position(iter_dim)-bottom_particle%Position(iter_dim)).lt.2*radius) then
                          print *,"error in random creation"
                          print *,"x-position",top_particle%Position(iter_dim),bottom_particle%Position(iter_dim)
                          call exit()
                        end if
                      end do
                      top_particle%radius=radius
                      top_particle%masse=top_particle%radius**3*pi*4/3*rho_to_use
                      mass_mean=mass_mean+top_particle%masse/(NP-NB)
                      array_of_particles(i) = top_particle
                    end do
                  end if
                end do
              end do
                  
              call write_init_file()
            else
              call read_init_file()
            end if
        end subroutine user_input

        subroutine write_init_file()
          integer iter_particle
          type(particle),pointer::p
          open(25,file=init_filename,status="replace")
          do iter_particle=1,Np
            p=>array_of_particles(iter_particle)
            write(25,*) p%masse,p%radius,p%Velocity,p%Position
          end do

          close(25,status="keep")
        end subroutine

        subroutine read_init_file()
          integer iter_particle
          type(particle),pointer::p
          character*18 massestring
          open(25,file=init_filename,status="old")
          call determine_number_particles(25)
          rewind(25) !go back to start of file
          call alloclist()
          do iter_particle=1,Np
            p=>array_of_particles(iter_particle)
            read(25,*) massestring,p%radius,p%Velocity,p%Position
            if (trim(massestring)=="Infinity") then
              p%masse=infi
              p%active=.false.
            else
              read(massestring,*)  p%masse
              p%active=.true.
            end if
            p%Force=0.0d0
            if (iter_particle.gt.nb) mass_mean=mass_mean+p%masse/(NP-NB)            
          end do
          close(25,status="keep")
        end subroutine read_init_file
        subroutine determine_number_particles(filenumber)
          integer,intent(in)::filenumber
          integer error
          type(particle) p
          NP=0
          NB=0
          do
            read(filenumber,*,iostat=error) p%masse,p%radius,p%Velocity,p%Position
            select case(error)
            case(0)
              NP=NP+1
              if (.not. ieee_is_finite(p%masse)) then
                NB=NB+1
              end if
            case(iostat_end)
              exit
            case default
              print *,"error in start-position file"
              call exit()
            end select
          end do

        end subroutine determine_number_particles
        subroutine alloclist()
          integer iter
            !allocates dynamically sized arrays
            allocate(array_of_particles(NP))
            call init_list()
            do iter=1,np
              array_of_particles(iter)%use_soft_sphere=.false.
              array_of_particles(iter)%Force=0.0d0
              array_of_particles(iter)%active=.true.
            end do
        end subroutine alloclist
end module setup
