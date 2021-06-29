program main
    
    use setup
    use type_particle
    use collisions
    use solid_dynamics
    use plot
    use linked_list_demo

    implicit none
    integer::n_steps,iter_step
    type(Particle)::colliding_particles(2)
    type(collision)::next_crash
    integer id_of_colliding_particles(2)
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

        call sort_list()

        next_crash=next_collision(any_collisions_left)

        !-----------EVENT-DRIVEN-HARD-SPHERE
        do while (any_collisions_left)
            
            call move(next_crash%time-acctim,array_of_particles)

            acctim=next_crash%time
            
            colliding_particles(:)=array_of_particles(id_of_colliding_particles)

            call collision_calculation(colliding_particles)

            call collision_update(array_of_particles,id_of_colliding_particles)

            next_crash=next_collision(any_collisions_left)
        end do
        
        call write_plot_file(array_of_particles)
    end do

    if (.true.) then
        print *,"running linked list demo"
        !call  execute_demo()
        call collision_list_demo
    end if
    
end program main

