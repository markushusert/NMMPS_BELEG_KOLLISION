module collisions
    use type_particle,only: Particle,array_of_particles
    use setup,only:k_stoss,Rhop,NP,dt,areawidth,NB
    use compile_constants
    use statistics,only:E_disp
    use ieee_arithmetic
    use collision_list_module, only:collision,collision_list_demo&
    ,next_collision,get_collision_list,clear_list,get_next,access_ptr,delete_next,insert_data&
    ,number_collisions,print_list,print_collision
    use linked_list, only:list_t
    

    implicit none

    

    contains
        subroutine collision_detection()
            !TODO fill collisionpartners, and tab_list based on particles in part_array
            !think about teleporting!            
            integer iter_1,iter_2
            call clear_list()

            do iter_1=NB+1,NP
                call check_particle_against_all_others(iter_1,0.0)
            end do

            if (.true.) then
                !call print_list(get_collision_list())
            end if
        end subroutine collision_detection
        subroutine check_particle_against_all_others(iter,acctim)
            integer,intent(in)::iter
            real,intent(in)::acctim
            integer iter2
            do iter2=NB+1,np
                if(.not.iter2.eq.iter) then
                    call process_particles(iter,iter2,acctim)
                end if
            end do
        end subroutine
        subroutine process_particles(iter_1,iter_2,acctim)
            integer,intent(IN)::iter_1,iter_2
            real,intent(in)::acctim

            type(collision)::current_collision
            type(list_t),pointer:: startpoint_list
            real collision_time

            startpoint_list=>get_collision_list()
            collision_time=calc_collision_time([array_of_particles(iter_1),array_of_particles(iter_2)])+acctim
            if (collision_time.lt.dt) then
                
                current_collision%time=collision_time
                current_collision%partners=[iter_1,iter_2]
                !print *,"collision between particles",iter_1,iter_2,"at time",collision_time
                call insert_collision_to_list(startpoint_list,current_collision)
            end if
        end subroutine process_particles

        function calc_collision_time(colliding_particles) result(time)
            real time
            type(Particle),intent(in)::colliding_particles(2)
            real,dimension(dim)::rel_vel,current_rel_pos
            real,dimension(dim,9)::rel_pos
            real rel_vel_norm
            real rel_pos_norm
            real sum_radius
            real temp_coll_time
            integer iter_config

            sum_radius=colliding_particles(1)%radius+colliding_particles(2)%radius
            rel_vel=colliding_particles(2)%Velocity-colliding_particles(1)%Velocity
            rel_pos=calc_rel_pos(colliding_particles)

            !set time to nan by default
            time=ieee_value(time, ieee_quiet_nan)

            do iter_config=1,9
                current_rel_pos=rel_pos(:,iter_config)

                if (norm2(current_rel_pos).eq.0.0) then
                    exit
                end if

                rel_vel_norm=sqrt(sum(rel_vel**2))
                rel_pos_norm=sqrt(sum(current_rel_pos**2))
                temp_coll_time=(-dot_product(rel_vel,current_rel_pos)&
                -sqrt(dot_product(rel_vel,current_rel_pos)**2-rel_vel_norm**2*(rel_pos_norm**2-sum_radius**2)))&
                /rel_vel_norm**2
                if (.not. (ieee_is_nan(temp_coll_time).or.temp_coll_time.lt.0.0)) then
                    !not also includes the case that time is nan
                    if (.not.temp_coll_time.gt.time) then
                        time=temp_coll_time
                    end if
                end if
                
            end do
        end function calc_collision_time

        function calc_rel_pos(colliding_particles)
            type(Particle),intent(IN)::colliding_particles(2)
            real,dimension(dim,9)::calc_rel_pos
            real,dimension(dim)::rel_pos_temp
            integer,DIMENSION(dim)::direction_vector
            integer iter_dim1,iter_dim2,counter


            !returns distance between 2 particles in 3 different konfigurations
            !1. no teleportation
            !2. part1 is teleported forward
            !3. part2 is teleported
            !this is needed to consider teleportation when determining collisions
            !the first entry calc_rel_pos(:,1) is guaranteed to be the shortest connection
            counter=0
            calc_rel_pos=0.0
            do iter_dim1=-1,1
                do iter_dim2=-1,1
                    if (dim.eq.3 .or. iter_dim2.eq.0) then
                        !only vary the second dimension if we have a 3d calculation
                        counter=counter+1
                        direction_vector=[0,0]
                        direction_vector(1)=iter_dim1
                        direction_vector(2)=iter_dim2
                        rel_pos_temp=colliding_particles(2)%Position-colliding_particles(1)%Position+areawidth*direction_vector
                        if (counter.eq.1)then
                            calc_rel_pos(:,1)=rel_pos_temp
                        else
                            if (norm2(rel_pos_temp).lt.norm2(calc_rel_pos(:,1))) then
                                !if distance using the current teleportations is smallest distance yet, put in first position
                                calc_rel_pos(:,counter)=calc_rel_pos(:,1)
                                calc_rel_pos(:,1)=rel_pos_temp
                            else
                                calc_rel_pos(:,counter)=rel_pos_temp
                            end if
                        end if
                    end if

                end do
            end do
        end function calc_rel_pos

        function check_list(startpoint_list,wrong_member)
            type(list_t),pointer::startpoint_list
            integer,optional,intent(out)::wrong_member
            type(list_t),pointer::current_element,next_element
            logical check_list
            logical next_ele_exists
            integer*8 used_adresses(1000)
            integer counter

            counter=0
            current_element=>startpoint_list
            used_adresses=0
            check_list=.true.
            do
                next_element=>get_next(current_element,next_ele_exists)
                if (next_ele_exists) then
                    counter=counter+1
                    if (any(used_adresses(1:counter).eq.loc(next_element))) then
                        print *,"found error at collision",counter
                        !call print_list(startpoint_list)
                        check_list=.false.

                        if (present(wrong_member)) then
                            wrong_member=counter
                        end if
                        return
                    end if
                    used_adresses(counter)=loc(next_element)
                    current_element=>next_element
                else
                    exit
                end if
                
            end do
        end function check_list
        
        subroutine insert_collision_to_list(startpoint_list,current_collision)
            type(list_t),pointer::startpoint_list
            type(list_t),pointer::current_element,next_element
            type(collision),intent(in)::current_collision
            type(collision),pointer::collision_in_list
            logical next_ele_exists,list_correct
            integer wrong_member
            integer counter

            current_element=>startpoint_list
            number_collisions=number_collisions+1
            counter=0
            do
                counter=counter+1
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
            list_correct=check_list(startpoint_list,wrong_member)
            if (.not.list_correct) then
                call print_list(startpoint_list,wrong_member)
            end if
            !call print_list(startpoint_list)
        end subroutine insert_collision_to_list

        subroutine collision_update(id_of_colliding_parts,acctim)
            integer,intent(in)::id_of_colliding_parts(2)
            real,intent(in):: acctim
            !TODO:
            !1 remove colliding_particles from collision list
            !2 compute new collisions of particle 1 and 2
            !3 insert new collisions in collision list
            
            call remove_ids_from_list(id_of_colliding_parts)
            !print *,"after removing"
            !call print_list(get_collision_list(),3)
            call check_particle_against_all_others(id_of_colliding_parts(1),acctim)
            !print *,"after update1"
            !call print_list(get_collision_list(),3)
            call check_particle_against_all_others(id_of_colliding_parts(2),acctim)
            !print *,"after update2"
            !call print_list(get_collision_list(),3)
        end subroutine collision_update

        subroutine remove_ids_from_list(id_of_colliding_parts)
            integer,intent(in)::id_of_colliding_parts(:)
            type(list_t),pointer::current_node,next_node,temp_node
            type(collision),pointer::current_collision
            logical next_element_exists
            integer iter,counter
            logical,parameter:: debug=.false.
            logical delete_node

            !debug=.false.
            if (id_of_colliding_parts(1).eq.61 .and. id_of_colliding_parts(2).eq.81) then
                !debug=.true.
            end if
            current_node=>get_collision_list()
            
            if (debug) then
                print *,"number_collisions",number_collisions
                call print_list(get_collision_list())
            end if
            counter=0
            do
                counter=counter+1
                next_node=>get_next(current_node,next_element_exists)
                if (next_element_exists) then
                    
                    current_collision=>access_ptr(next_node)
                    if (debug) then
                        !call print_list(current_node,3)
                        print *,counter,"comparing collision-id",current_collision%partners
                    end if
                    delete_node=.false.
                    do iter=1,2
                        if (any(id_of_colliding_parts.eq.current_collision%partners(iter))) then
                            delete_node=.true.
                            
                            exit
                        end if
                    end do

                    if (delete_node) then
                        call delete_next(current_node,temp_node)
                        number_collisions=number_collisions-1
                    else
                        current_node=>get_next(current_node,next_element_exists)
                    end if      
                    
                else
                    exit
                end if
                
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
            real, dimension(2) :: vector12,velocity_vector12
            real,dimension(dim,9)::vector12_temp

            vector12_temp = calc_rel_pos(colliding_particles)
            !shortest possible collision time
            vector12=vector12_temp(:,1)
            vector12=vector12/norm2(vector12)!normalise distance vector
            velocity_vector12 = colliding_particles(2)%Velocity - colliding_particles(1)%Velocity
            !source: Ouyang: Particle-motion-resolved discrete model for simulating gas-solid fluidization, 1999
            colliding_particles(1)%Velocity =  colliding_particles(1)%Velocity + &
            1.0/(1.0+colliding_particles(1)%masse/colliding_particles(2)%masse)&
            * (1+k_stoss) *  dot_product(velocity_vector12 ,vector12)* vector12
            colliding_particles(2)%Velocity =  colliding_particles(2)%Velocity - &
            1.0/(1.0+colliding_particles(2)%masse/colliding_particles(1)%masse)&
            * (1+k_stoss) *  dot_product(velocity_vector12 ,vector12)* vector12

             !calculate energy_loss
            E_disp=E_disp+&
            (1-k_stoss**2)*0.5+dot_product(velocity_vector12 ,vector12)**2&
            /(1/colliding_particles(1)%masse+1/colliding_particles(2)%masse)
            !end Yvi
        end subroutine collision_calculation

end module collisions
