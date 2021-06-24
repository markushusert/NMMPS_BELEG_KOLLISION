module plot
    use type_particle,only:Particle
    implicit none

    contains
        subroutine write_plot_file(particle_array)
            type(Particle),intent(in)::particle_array(:)

        end subroutine write_plot_file
end module plot