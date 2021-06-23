module collisions
    use type_particle

    implicit none

    integer,dimension(:,:),ALLOCATABLE::collisionpartners
    real,dimension(:),ALLOCATABLE::tab_list
    contains
        subroutine collision_detection(part_array)
            type(Particle),intent(in)::part_array(:)
            print *,"location of argument",loc(part_array)
            print *,"location of global",loc(array_of_particles)
        end subroutine collision_detection

end module collisions