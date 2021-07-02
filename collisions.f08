module collisions
    use type_particle,only: Particle
    use setup,only:k_stoss,Rhop,NP,dt
    use compile_constants
    use collision_list_module, only:collision,collision_list_demo&
    ,next_collision,get_collision_list,clear_list,get_next,access_ptr,delete_next,insert_data
    use linked_list, only:list_t

    implicit none

    

    contains
        subroutine collision_detection(part_array)
            !TODO fill collisionpartners, and tab_list based on particles in part_array
            !think about teleporting!
            type(Particle),intent(in)::part_array(:)
            
            
            integer iter_1,iter_2

            

            call clear_list()

            do iter_1=1,NP
                do iter_2=1,NP
                    if (.not.iter_2.eq.iter_1)then
                        call process_particles(part_array,iter_1,iter_2,0.0)
                    end if
                end do
            end do

        end subroutine collision_detection

        subroutine process_particles(part_array,iter_1,iter_2,acctim)
            type(Particle),intent(in)::part_array(:)
            integer,intent(IN)::iter_1,iter_2
            real,intent(in)::acctim

            type(collision)::current_collision
            type(list_t),pointer:: startpoint_list
            real collision_time

            startpoint_list=>get_collision_list()
            collision_time=calc_collision_time(part_array(iter_1),part_array(iter_1))+acctim
            if (collision_time.lt.dt) then
                current_collision%time=collision_time
                current_collision%partners=[iter_1,iter_2]
                call insert_collision_to_list(startpoint_list,current_collision)
            end if
        end subroutine process_particles

        function calc_collision_time(particle1,particle2) result(time)
            real time
            type(Particle),intent(in)::particle1,particle2
            real,dimension(dim)::rel_vel,rel_pos
            real rel_vel_norm
            real rel_pos_norm
            real sum_radius

            sum_radius=particle1%radius+particle2%radius
            rel_vel=particle1%Velocity-particle2%Velocity
            rel_pos=particle1%Position-particle2%Position

            rel_vel_norm=sqrt(sum(rel_vel**2))
            rel_pos_norm=sqrt(sum(rel_pos**2))

            time=(-dot_product(rel_vel,rel_pos)&
            -sqrt(dot_product(rel_vel,rel_pos)**2)-rel_vel_norm**2*(rel_pos_norm**2-sum_radius**2))&
            /rel_vel_norm**2
        end function calc_collision_time

        subroutine insert_collision_to_list(startpoint_list,current_collision)
            type(list_t),pointer::startpoint_list
            type(list_t),pointer::current_element,next_element
            type(collision),intent(in)::current_collision
            type(collision),pointer::collision_in_list
            logical next_ele_exists

            current_element=>startpoint_list
            
            do
                next_element=>get_next(current_element,next_ele_exists)
                if (next_ele_exists) then
                    collision_in_list=>access_ptr(next_element)
                    if (collision_in_list%time.lt.current_collision%time) then
                        current_element=>next_element
                        cycle!do not insert element yet
                    end if
                end if
                call insert_data(current_element,current_collision)
                exit
            end do

        end subroutine insert_collision_to_list

        subroutine collision_update(part_array,id_of_colliding_parts,acctim)
            type(Particle),intent(in)::part_array(:)
            integer,intent(in)::id_of_colliding_parts(2)
            real,intent(in):: acctim
            !TODO:
            !1 remove colliding_particles from collision list
            !2 compute new collisions of particle 1 and 2
            !3 insert new collisions in collision list

            call remove_ids_from_list(id_of_colliding_parts)
            call process_particles(part_array,id_of_colliding_parts(1),id_of_colliding_parts(2),acctim)

        end subroutine collision_update

        subroutine remove_ids_from_list(id_of_colliding_parts)
            integer,intent(in)::id_of_colliding_parts(:)
            type(list_t),pointer::current_node,next_node
            type(collision),pointer::current_collision
            logical next_element_exists
            integer iter

            current_node=>get_collision_list()

            do
                next_node=>get_next(current_node,next_element_exists)
                if (next_element_exists) then
                    current_collision=>access_ptr(next_node)
                    do iter=1,2
                        if (any(id_of_colliding_parts.eq.current_collision%partners(iter))) then
                            call delete_next(current_node)
                        end if
                    end do
                else
                    exit
                end if
                current_node=>next_node
            end do
        end subroutine remove_ids_from_list

        subroutine sort_list()

            !TODO sort collisionpartners and tab_list chronologically
            !info collisionpartners and tab_list are accesible inside of this subroutine
            !as global variables

        end subroutine

        subroutine collision_calculation(colliding_particles)
            Type(Particle),intent(inout)::colliding_particles(2)

            !TODO calculate new velocities of Pcolliding_particles based on the Stoßgesetz
            !Yvi
            real, dimension(2) :: vector12, vector21,velocity_vector12
            vector12 = colliding_particles(2)%Position - colliding_particles(1)%Position
            velocity_vector12 = colliding_particles(2)%Velocity - colliding_particles(1)%Velocity
            !source: Ouyang: Particle-motion-resolved discrete model for simulating gas-solid fluidization, 1999
            colliding_particles(1)%Velocity =  colliding_particles(1)%Velocity + &
            1.0/2.0 * (1+k_stoss) *  (velocity_vector12 * vector12)/(norm2(vector12))**2 * vector12
            colliding_particles(2)%Velocity =  colliding_particles(2)%Velocity -&
             1.0/2.0 * (1+k_stoss) *  (velocity_vector12 * vector12)/(norm2(vector12))**2 * vector12

             !calculate energy_loss
            !end Yvi
        end subroutine collision_calculation

end module collisions
