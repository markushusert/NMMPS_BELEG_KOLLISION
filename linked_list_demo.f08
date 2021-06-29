module linked_list_demo

    use linked_list

    implicit none
    type :: pair_demo
        integer,DIMENSION(2)::partners
    end type pair_demo

    ! A trick to allow us to store pointers in the list
    type :: pair_demo_ptr
        type(pair_demo), pointer :: p
    end type pair_demo_ptr

    
    

    contains
    subroutine execute_demo()
        type(list_t),pointer::start_node
        type(list_t),pointer::next_node
        integer iter
        type(pair_demo),dimension(3),target::collision_partners !dummy collision partners we want to store
        type(pair_demo),target::new_pair_demo
        type(pair_demo_ptr)::pointer_to_collision_partners
        collision_partners(1)%partners=[1,2]
        collision_partners(2)%partners=[3,4]
        collision_partners(3)%partners=[5,6]
        !start list on start_node, do not put data to start node
        call list_init(start_node)
        
        next_node=>start_node

        print *,"filling linked list with dummy data"
        do iter=1,3
            !make pointer point to data we want to store
            pointer_to_collision_partners%p=>collision_partners(iter)
            print *,"storing value",pointer_to_collision_partners%p%partners,"in position",iter

            !store data in linked list
            call list_insert(next_node,DATA=transfer(pointer_to_collision_partners, list_data))
            !move next_node to node we just created (optional)
            next_node=>list_next(next_node)
        end do

        print *,"printing values of linked list"
        next_node=>list_next(start_node)
        iter=0
        do while (associated(next_node))!as long as a new node is found

            iter=iter+1
            !make pointer point to data we want to store
            pointer_to_collision_partners=Transfer(list_get(next_node),pointer_to_collision_partners)
            print *,"reading value",pointer_to_collision_partners%p%partners,"in position",iter

            next_node=>list_next(next_node)
        end do

        print *,"manipulating list"

        print *,"changing first partner of first pair_demo to 7"
        new_pair_demo%partners=[7,2]
        pointer_to_collision_partners%p=>new_pair_demo
        call list_put(list_next(start_node),transfer(pointer_to_collision_partners,list_data))


        print *,"removing second node from list"
        next_node=>start_node
        !to get to second node, we have to iterate from the beginning
        do iter=1,2-1
            next_node=>list_next(next_node)
        end do

        !next node is the first data node here, because subroutine list_pop deletes the next node after its argument
        call list_delnext(next_node)
        !call test()

        print *,"printing values of list after manipulation"
        iter=0
        do while (associated(next_node))!as long as a new node is found
            iter=iter+1
            !make pointer point to data we want to store
            pointer_to_collision_partners=Transfer(list_get(next_node),pointer_to_collision_partners)
            print *,"reading value",pointer_to_collision_partners%p%partners,"in position",iter
            next_node=>list_next(next_node)
        end do
        
         print *,"die Bedienung der linked list ist aktuell noch etwas unelegant"
         print *,"ich schreibe als erstes functionen, damit wir die elemente der linked list einfacher erreichen können"
         print *,"den datentypen würde ich trotzdem sehr gerne verwenden, da er gut für unser problem geeignet ist"
    end subroutine execute_demo

end module linked_list_demo
