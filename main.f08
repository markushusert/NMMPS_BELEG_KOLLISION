program main
    
    use setup
    use type_particle
    use collisions

    implicit none
    print *,"NP beofre user input",NP
    call user_input()
    print *,"NP after user input",NP
    call collision_detection(array_of_particles)

    
end program main

