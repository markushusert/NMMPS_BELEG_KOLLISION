program main
    
    use settings
    use setup

    implicit none
    print *,"NP beofre user input",NP
    call user_input()
    print *,"NP after user input",NP
end program main

