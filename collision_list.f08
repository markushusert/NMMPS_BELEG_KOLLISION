module collision_list_module
    use linked_list
    implicit none
    !integer,dimension(:,:),ALLOCATABLE::collisionpartners
    !real,dimension(:),ALLOCATABLE::tab_list

    type :: collision
        integer,DIMENSION(2)::partners
        real time!measured from start of timestep
    end type collision

    ! A trick to allow us to store pointers in the list
    type :: collision_ptr
        type(collision), pointer :: p
    end type collision_ptr

    type(list_t), pointer,PRIVATE :: collision_list => null()
    type(list_t), pointer :: last_collision => null()

    integer number_collisions
    private ::access_ptr!function noch nicht implementier, brauchen wir wohl auch nicht

    contains
        subroutine clear_list()
            !clears linked lists for tab and collision partners
            call list_free(list_next(collision_list))
            number_collisions=0
        end subroutine clear_list

        function next_collision(found_coll)
            !returns the collision which is the next one to occur
            type(collision)::next_collision
            type(list_t),pointer::first_ele
            logical found_coll

            first_ele=>get_next(collision_list,found_coll)
            if (found_coll) then
                next_collision=access_copy(first_ele)
            end if

        end function next_collision

        function access_copy(list_element) result(data)
            !creates a copy of the collision stored in list_element
            type(list_t), pointer,intent(in)::list_element
            type(collision)::data
            type(collision),pointer::temp
            temp=>access_ptr(list_element)
            data=temp
            
        end function access_copy

        function access_ptr(list_element) result(data)
            !creates a link to the collision stored in list_element
            !if you manipulate data returned by this subroutine, the original list will change too
            type(list_t), pointer,intent(in)::list_element
            type(collision),pointer::data
            integer, dimension(:), pointer :: raw_data
            type(collision_ptr)::data_ptr

            raw_data=>list_get(list_element)
            if (associated(raw_data)) then
                data_ptr=Transfer(raw_data,data_ptr)
            else
                print *,"list-element contains no value"
                if (loc(list_element).eq.loc(collision_list)) then
                    print *,"do not try to access the starting point of the list"
                end if
                call exit()
            end if
            
            data=>data_ptr%p
        end function access_ptr
        subroutine set_data(data,list_element)
            !inserts the collision data after the element list_element
            type(collision),intent(in)::data
            type(collision_ptr),pointer::data_ptr
            type(list_t),pointer,intent(inout)::list_element

            allocate(data_ptr%p)
            data_ptr%p=data
            call list_put(list_element,DATA=transfer(data_ptr, list_data))

        end subroutine set_data

        subroutine insert_data(list_element,data)
            !inserts the collision data after the element list_element
            type(collision),intent(in),value::data
            type(list_t),pointer,intent(inout)::list_element
            type(collision_ptr)::data_ptr

            allocate(data_ptr%p)
            data_ptr%p=data
            call list_insert(list_element,DATA=transfer(data_ptr, list_data))

        end subroutine insert_data

        subroutine insert_element(element_to_insert,list_element)
            !inserts the element_to_insert after list_element
            type(list_t),pointer,intent(inout)::element_to_insert
            type(list_t),pointer,intent(inout)::list_element
            
            call list_insert_node(list_element,element_to_insert)

        end subroutine insert_element

        function cut_next(list_element,found_next) result(next_element)
            !returns a pointer to the next element of the list
            !and cuts it out of the linked list
            !this function can help with moving an element
            !if found next is false the result must not be used since no next element exists
            !(end of list reached)
            type(list_t), pointer,intent(in)::list_element
            type(list_t), pointer::next_element
            logical found_next
 
            next_element=>list_cut_next(list_element)

            found_next=ASSOCIATED(next_element)
        end function cut_next

        subroutine print_list(list_element)
            !prints all elements in linked list starting after list_element
            type(list_t), pointer,intent(in)::list_element
            type(list_t), pointer::next_element,current_element,temp
            type(collision),pointer:: data
            logical nextflag,tempflag

            current_element=>list_element
            do 
                next_element=>get_next(current_element,nextflag)
                if (.not.nextflag) then
                    !last element reached
                    exit
                end if
                current_element=>next_element

                data=>access_ptr(current_element)
                if (.not.ASSOCIATED(data)) then
                    print *,"data has not been associated"
                    call exit()
                end if
                temp=>get_next(current_element,nextflag)
                print *,"element at position",loc(current_element),"next_element",loc(tempflag)
                call print_collision(data)
            end do
        end subroutine print_list
        subroutine print_collision(data)
            type(collision),intent(in):: data
            print *,"kollision between partikles",data%partners(1),"and",data%partners(2)&
            ,"at time",data%time
        end subroutine print_collision
        function get_next(list_element,found_next) result(next_element)
            !returns a pointer to the next element of the list
            !if found next is false the result must not be used since no next element exists
            !(end of list reached)
            type(list_t), pointer,intent(in)::list_element
            type(list_t), pointer::next_element
            logical found_next
     
            next_element=>list_next(list_element)
            found_next=ASSOCIATED(next_element)
        end function get_next

        subroutine init_list()
            !inits the start node for both linked lists
            !start list on start_node, do not put data to start node
            if (associated(collision_list)) then
                print *,"ERROR, only initialise the collisionlist once"
                call exit()
            else
                call list_init(collision_list)
                number_collisions=0
            end if
            
        end subroutine init_list
        function get_collision_list() result (startpoint)
            type(list_t),pointer::startpoint
            startpoint=>collision_list
        end function
        subroutine collision_list_demo()
            !this demo does the following
            !   -initialise the collision list
            !   -add 3 collisions to it
            !   -remove the 2d collision
            !   -modify the 3rd collision
            !   -acces the 1s collision
            type(collision),dimension(3)::collision_array
            type(collision),pointer::my_coll
            type(collision) next_collision_obj
            type(list_t),pointer::start_of_list,current_node
            logical flag
            integer i

            !fill collision array
            collision_array(1)%partners=(/1,2/)
            collision_array(2)%partners=(/3,4/)
            collision_array(3)%partners=(/5,6/)
            collision_array(1)%time=0.5
            collision_array(2)%time=0.3
            collision_array(3)%time=0.2

            
            start_of_list=>get_collision_list()
            current_node=>start_of_list

            do i=1,3
                call insert_data(current_node,collision_array(i))
                current_node=>list_next(current_node)
            end do
            print *,"filled linked list"
            
            call print_list(start_of_list)

            print *,"removing second element"

            call list_delnext(list_next(start_of_list))

            call print_list(start_of_list)

            print *,"modifying now second element"

            my_coll=>access_ptr(list_next(list_next(start_of_list)))

            my_coll%time=2.0

            call print_list(start_of_list)

            print *,"fetching first collision"

            next_collision_obj=next_collision(flag)

            if (flag) then
                print *,"first colliison object received"
                call print_collision(next_collision_obj)
            end if

            next_collision_obj%partners=(/10,20/)

            print *,"inserting new collision at first place"

            call insert_data(start_of_list,next_collision_obj)

            call print_list(start_of_list)
        end subroutine collision_list_demo
end module collision_list_module