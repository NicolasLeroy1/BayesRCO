program rng_test
    use mod_defs, only: dp
    use mod_random
    implicit none

    integer :: i, n
    integer, dimension(:), allocatable :: seed
    character(len=20) :: arg
    integer :: user_seed

    ! Initialize the seed from the command line
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read(arg, *) user_seed
    else
        user_seed = 12345
    end if

    call random_seed(size=n)
    allocate(seed(n))
    ! Seeding logic matching Fortran 13 / libgfortran default
    ! Let's just use user_seed for all elements or something predictable
    seed = [(user_seed + i - 1, i = 1, n)]
    call random_seed(put=seed)

    print *, "Uniforms:"
    do i = 1, 10
        print "(F20.16)", rand_uniform(0.0_dp, 1.0_dp)
    end do

    print *, "Normals:"
    do i = 1, 10
        print "(F20.16)", rand_normal(0.0_dp, 1.0_dp)
    end do

    print *, "Gammas (shape=2.0, scale=1.0):"
    do i = 1, 10
        print "(F20.16)", rand_gamma(2.0_dp, 1.0_dp)
    end do

end program rng_test
