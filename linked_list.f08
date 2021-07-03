! A generic linked list object
!from http://fortranwiki.org/fortran/show/gen_list
module linked_list
    implicit none
  
    private
  
    public :: list_t
    public :: list_data
    public :: list_init
    public :: list_free
    public :: list_insert
    public :: list_put
    public :: list_get
    public :: list_next
    public :: list_delnext
    public :: list_insert_node
    public :: list_cut_next
  
    ! A public variable to use as a MOLD for transfer()
    integer, dimension(:), allocatable :: list_data
  
    ! Linked list node data type
    type :: list_t
       private
       integer, dimension(:), pointer :: data => null()
       type(list_t), pointer :: next => null()
    end type list_t
  
  contains
  
    ! Initialize a head node SELF and optionally store the provided DATA.
    subroutine list_init(self, data)
      type(list_t), pointer :: self
      integer, dimension(:), intent(in), optional :: data
  
      allocate(self)
      nullify(self%next)
  
      if (present(data)) then
         allocate(self%data(size(data)))
         self%data = data
      else
         nullify(self%data)
      end if
    end subroutine list_init
  
    ! Free the entire list and all data, beginning at SELF
    subroutine list_free(self)
      type(list_t), pointer :: self
      type(list_t), pointer :: current
      type(list_t), pointer :: next
  
      current => self
      do while (associated(current))
         next => current%next
         if (associated(current%data)) then
            deallocate(current%data)
            nullify(current%data)
         end if
         deallocate(current)
         nullify(current)
         current => next
      end do
    end subroutine list_free

    ! Free the node following after self
    subroutine list_delnext(self,next_ele)
        type(list_t), pointer :: self
        type(list_t), pointer :: current
        type(list_t), pointer :: next
        type(list_t), pointer :: next_next
        type(list_t), pointer,optional ::next_ele
    
        !get nodes 1 place and 2 places after self
        current => self
        next => current%next
        next_next =>next%next
        
        !delete node 1 place after self
        if (associated(next%data)) then
            deallocate(next%data)
            nullify(next%data)
        end if
        deallocate(next)
        nullify(next)

        !fill the gap with node 2 places after
        current%next=>next_next
        if (present(next_ele)) then
            next_ele=>current%next
        end if

    end subroutine list_delnext
 
    ! Return the next node after SELF
    function list_next(self) result(next_node)
      type(list_t), pointer :: self
      type(list_t), pointer :: next_node
      next_node => self%next
    end function list_next
    
    function list_cut_next(self) result(next_node)
      type(list_t), pointer :: self
      type(list_t), pointer :: next_node
      next_node => self%next

      !remove next_node by having self point to its next 
      self%next=>next_node%next
    end function list_cut_next

    subroutine list_insert_node(self,node_to_insert)
      type(list_t), pointer,intent(inout) :: self
      type(list_t), pointer,intent(inout) :: node_to_insert
      type(list_t), pointer :: next_node

      !if a node was following after self, make node_to_insert point to it
      if (associated(self%next)) then
         next_node=>self%next
         node_to_insert%next=>next_node
      end if
      self%next=>node_to_insert
   end subroutine list_insert_node
    ! Insert a list node after SELF containing DATA (optional)
    subroutine list_insert(self, data)
      type(list_t), pointer :: self
      integer, dimension(:), intent(in), optional :: data
      type(list_t), pointer :: next
  
      allocate(next)
  
      if (present(data)) then
         allocate(next%data(size(data)))
         next%data = data
      else
         nullify(next%data)
      end if
  
      next%next => self%next
      self%next => next
    end subroutine list_insert
  
    ! Store the encoded DATA in list node SELF
    subroutine list_put(self, data)
      type(list_t), pointer :: self
      integer, dimension(:), intent(in) :: data
  
      if (associated(self%data)) then
         deallocate(self%data)
         nullify(self%data)
      end if
      
      allocate(self%data(size(data)))
      self%data = data
    end subroutine list_put
  
    ! Return the DATA stored in the node SELF
    function list_get(self) result(data)
      type(list_t), pointer :: self
      integer, dimension(:), pointer :: data
      data => self%data
    end function list_get
  
end module linked_list