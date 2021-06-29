program main
    
    use setup
    use type_particle
    use collisions
    use solid_dynamics
    use statistics
    use plot
    use linked_list_demo

    implicit none
    integer::n_steps,iter_step
    type(Particle)::colliding_particles(2)
    type(collision)::next_crash
    logical any_collisions_left
    real acctim

    !-----------SETUP
    print *,"starting program"
    call user_input()
    n_steps=floor(t_ges/dt)

    !-----------MAIN-LOOP
    do iter_step=1,n_steps
        acctim=0.0
        call time_integration(array_of_particles)

        call collision_detection(array_of_particles)

        !-----------EVENT-DRIVEN-HARD-SPHERE
        do 
            next_crash=next_collision(any_collisions_left)

            if (.not.any_collisions_left) then
                exit !break out of loop
            end if
            
            call move(next_crash%time-acctim,array_of_particles)

            acctim=next_crash%time
            
            colliding_particles(:)=array_of_particles(next_crash%partners)

            call collision_calculation(colliding_particles)

            call collision_update(array_of_particles,next_crash%partners)

        end do
        
        time=time+dt
        
        call write_plot_file(array_of_particles)
    end do

    if (.true.) then
        print *,"running linked list demo"
        !call  execute_demo()
        call collision_list_demo
    end if
    
end program main

