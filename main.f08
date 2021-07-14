program main
    
    use setup
    use type_particle
    use collisions
    use solid_dynamics
    use statistics
    use plot
    use linked_list_demo
    USE ieee_arithmetic

    implicit none
    integer::n_steps,iter_step,counter
    type(Particle)::colliding_particles(2)
    type(collision)::next_crash
    logical any_collisions_left
    DOUBLE PRECISION,target:: acctim
    real,pointer::temp
    real dummy
    
    call user_input()
    n_steps=floor(t_ges/dt)
    print *,"starting program, calculating",n_steps,"timesteps"
   
    !-----------MAIN-LOOP
    do iter_step=1,n_steps
        acctim=0.0d0
        counted_collisions=0
        call time_integration(0.5)

        call collision_detection()
        !call print_list(get_collision_list())

        !-----------EVENT-DRIVEN-HARD-SPHERE
        do 
            !print *,"next 3 colls"
            !call print_list(get_collision_list(),3)
            next_crash=next_collision(any_collisions_left)
            

            if (.not.any_collisions_left) then
                exit !break out of loop
            end if
            counted_collisions=counted_collisions+1

            call print_collision(next_crash)
            
            call move(next_crash%time-acctim)

            acctim=next_crash%time
            
            colliding_particles(:)=array_of_particles(next_crash%partners)

            call collision_calculation(colliding_particles)

            !colliding_particles was a seperate copy of our array, so we need to update array as well
            array_of_particles(next_crash%partners)=colliding_particles

            call collision_update(next_crash%partners,acctim)

            dummy=0.0
        end do

        call move(DT-acctim) !MOVE REMAINING TIME
        
        call update_active_status()
        time=time+dt
        current_timestep=current_timestep+1
        
        call time_integration(0.5)

        do counter=1,np
            if (array_of_particles(counter)%Position(dim).lt.0.0) then
                print *,"particle",counter,"fell through"
            end if
        end do

        call write_plot_file()
    end do

    if (.false.) then
        print *,"running linked list demo"
        !call  execute_demo()
        call collision_list_demo
    end if
    print *,"----------------programm finished-------------------"
end program main

