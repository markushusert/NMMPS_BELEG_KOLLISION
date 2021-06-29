module collisions
    use type_particle,only: Particle
    use setup,only:k_stoss,Rhop
    use collision_list_module, only:collision,collision_list_demo,next_collision
    implicit none

    

    contains
        subroutine collision_detection(part_array)
            type(Particle),intent(in)::part_array(:)
            
            !TODO fill collisionpartners, and tab_list based on particles in part_array
            !think about teleporting!

        end subroutine collision_detection

        subroutine collision_update(part_array,id_of_colliding_parts)
            type(Particle),intent(in)::part_array(:)
            integer,intent(in)::id_of_colliding_parts(2)
            !TODO:
            !1 remove colliding_particles from collision list
            !2 compute new collisions of particle 1 and 2
            !3 insert new collisions in collision list

        end subroutine collision_update

        subroutine sort_list()

            !TODO sort collisionpartners and tab_list chronologically
            !info collisionpartners and tab_list are accesible inside of this subroutine
            !as global variables

        end subroutine

        subroutine collision_calculation(colliding_particles)
            Type(Particle),intent(inout)::colliding_particles(2)

            !TODO calculate new velocities of Pcolliding_particles based on the Stoßgesetz
        end subroutine collision_calculation

end module collisions