module collision_list
    use linked_list
    implicit none
    !integer,dimension(:,:),ALLOCATABLE::collisionpartners
    !real,dimension(:),ALLOCATABLE::tab_list

    type :: pair
        integer,DIMENSION(2)::partners
    end type pair

    ! A trick to allow us to store pointers in the list
    type :: pair_ptr
        type(pair), pointer :: p
    end type pair_ptr

    ! A trick to allow us to store pointers in the list
    type :: real_ptr
        real, pointer :: p
    end type real_ptr

    type(list_t), pointer :: tab_list => null()
    type(list_t), pointer :: collisionpartners => null()
     

    contains
        subroutine clear_lists()
            !clears linked lists for tab and collision partners
        end subroutine clear_lists

        subroutine init_lists()
            !inits the start node for both linked lists
        end subroutine init_lists
end module collision_list