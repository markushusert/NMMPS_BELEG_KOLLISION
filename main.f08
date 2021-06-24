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
    integer id_of_colliding_particles(2)
    real acctim
    real time_next_collision !measured from the start of the timestep
    !so that we do not need to update collision timeings after each move

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

        time_next_collision=10000.0!insert first node of linked list here

        !-----------EVENT-DRIVEN-HARD-SPHERE
        do while (time_next_collision.le.dt)
            call move(time_next_collision-acctim,array_of_particles)

            acctim=time_next_collision
            
            id_of_colliding_particles=collisionpartners(:,1)

            colliding_particles(:)=array_of_particles(id_of_colliding_particles)

            call collision_calculation(colliding_particles)

            call collision_update(array_of_particles,id_of_colliding_particles)

            time_next_collision=tab_list(1)
        end do
        
        call write_plot_file(array_of_particles)
    end do

    if (.true.) then
        print *,"running linked list demo"
        call  execute_demo()
    end if
    
end program main

