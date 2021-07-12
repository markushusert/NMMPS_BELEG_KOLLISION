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

    real radius             !radius of particles
    real Rhop           !density of particles
    real k_stoss        !Stoßzahl
    real g              !Erdbeschleunigung
    real dt             !timestep
    real t_ges          !duration
    real areawidth      !width of simulation area
    real v_init
    real,parameter::pi=3.14159265
    character use_config
    integer i
    integer iter_layer,iter_place
    integer n,n_sphere_in_layer,n_layer
    real height
    integer, allocatable :: seed(:)
    logical read_start_positions
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
                        read_start_positions
    real                 :: random_value_x, random_value_y
    real infi
    
    
    

    contains
        subroutine user_input()
          logical file_exists
          
            
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
            NP=n_layer*n_sphere_in_layer+NB
            areawidth=radius*Nb*2
            if (n_layer.gt.NB/2) then
              print *,"there cannot be more than NB/2 spheres in a layer"
              call exit()
            end if
            !Yvi
            call alloclist()
            !TODO fill attributes of particles in array_of_particles
            IF (ieee_support_inf(infi)) THEN
              infi = ieee_value(infi,  ieee_positive_inf)
            END IF
            !check if start_position file exists
            INQUIRE(FILE=init_filename, EXIST=file_exists)
            if (.not.read_start_positions.or. .not.file_exists) then 
              print *,"randomly generating particle positions"      
              do i=1,NB
                bottom_particle%velocity=[0.0,0.0]
                bottom_particle%position=[(i-1)*2*radius, 0.0]
                bottom_particle%radius=radius
                bottom_particle%masse=infi
                array_of_particles(i) = bottom_particle
              end do
              i=nb
              do iter_layer=0,n_layer-1
                do iter_place=0,n_sphere_in_layer-1
                  i=i+1
                  call random_number(random_value_x)

                  top_particle%velocity=[0.0,v_init]
                  start_intervall=radius*2*real(NB)/n_sphere_in_layer*iter_place
                  length_intervall=radius*2*real(NB)/n_sphere_in_layer
                  !each particle gets an intervall of NB/n_sphere_in_layer radii in x-direction where it can be placed wherever
                  top_particle%position=[start_intervall+length_intervall*random_value_x, height+iter_layer*radius*3]
                  top_particle%radius=radius
                  top_particle%masse=top_particle%radius**3*pi*4/3*Rhop
                  array_of_particles(i) = top_particle
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

          do iter_particle=1,Np
            p=>array_of_particles(iter_particle)
            read(25,*) massestring,p%radius,p%Velocity,p%Position
            if (trim(massestring)=="Infinity") then
              p%masse=infi
            else
              read(massestring,*)  p%masse
            end if
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
            !allocates dynamically sized arrays
            allocate(array_of_particles(NP))
            call init_list()
        end subroutine alloclist
end module setup
