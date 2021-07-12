module compile_constants
    !purpose: define constants known at compile time
    integer,parameter:: dim=3  !choose 2D or 3D simulation
    character*64,parameter::init_filename="starting_positions.txt"
end module compile_constants