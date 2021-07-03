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
    integer::n_steps,iter_step
    type(Particle)::colliding_particles(2)
    type(collision)::next_crash
    logical any_collisions_left
    real,target:: acctim
    real,pointer::temp
    
    call user_input()
    n_steps=floor(t_ges/dt)
    print *,"starting program, calculating",n_steps,"timesteps"
   
    !-----------MAIN-LOOP
    do iter_step=1,n_steps
        acctim=0.0
        call time_integration()

        call collision_detection()

        !-----------EVENT-DRIVEN-HARD-SPHERE
        do 
            next_crash=next_collision(any_collisions_left)

            if (.not.any_collisions_left) then
                exit !break out of loop
            end if
            
            call move(next_crash%time-acctim)

            acctim=next_crash%time
            
            colliding_particles(:)=array_of_particles(next_crash%partners)

            call collision_calculation(colliding_particles)

            call collision_update(next_crash%partners,acctim)

        end do
        
        time=time+dt
        current_timestep=current_timestep+1
        
        call write_plot_file()
    end do

    if (.false.) then
        print *,"running linked list demo"
        !call  execute_demo()
        call collision_list_demo
    end if
    print *,"----------------programm finished-------------------"
end program main

