program test_mcmc_utils
    use mod_defs, only: dp
    use mod_data, only: ModelConfig, GenomicData, MCMCState, MCMCStorage
    use mod_io
    use mod_stats
    use mod_mcmc_utils
    use mod_random, only: rand_chi_square, rand_normal, rand_scaled_inverse_chi_square, rdirichlet
    implicit none

    ! Local declarations
    type(ModelConfig) :: config
    type(GenomicData), target :: gdata
    type(MCMCState), target :: mstate
    type(MCMCStorage) :: mstore

    integer :: i, j, n
    integer, dimension(:), allocatable :: seed

    ! Seed RNG for reproducibility
    call random_seed(size=n)
    allocate(seed(n))
    seed = [(12345 + i - 1, i = 1, n)]
    call random_seed(put=seed)

    ! Setup minimal config
    config%ncat = 2
    config%ndist = 4
    config%VCE = .true.
    config%dfvara = 4.0d0
    config%dfvare = 4.0d0
    config%vara_ap = 0.01d0
    config%vare_ap = 0.01d0

    ! Setup minimal gdata
    gdata%nloci = 10
    gdata%nind = 5
    gdata%nt = 5

    ! Allocate arrays
    allocate(gdata%X(gdata%nt, gdata%nloci))
    allocate(gdata%xpx(gdata%nloci))
    allocate(gdata%nannot(gdata%nloci))
    allocate(gdata%C(gdata%nloci, config%ncat))
    allocate(gdata%permvec(gdata%nloci))
    allocate(gdata%why(gdata%nind))
    allocate(gdata%trains(gdata%nind))
    allocate(gdata%snptracker(gdata%nloci, config%ncat))

    allocate(mstate%g(gdata%nloci))
    allocate(mstate%yadj(gdata%nt))
    allocate(mstate%gp(config%ndist))
    allocate(mstate%gpin(config%ndist))
    allocate(mstate%p(config%ndist, config%ncat))
    allocate(mstate%log_p(config%ndist, config%ncat))
    allocate(mstate%log_gp(config%ndist))
    allocate(mstate%vare_gp(config%ndist))
    allocate(mstate%delta(config%ndist))
    allocate(mstate%dirx(config%ndist))
    allocate(mstate%snpindist(config%ndist, config%ncat))
    allocate(mstate%varindist(config%ndist, config%ncat))
    allocate(mstate%ytemp(gdata%nloci))

    allocate(mstore%gstore(gdata%nloci))
    allocate(mstore%pstore(config%ndist, config%ncat))
    allocate(mstore%snpstore(config%ndist, config%ncat))
    allocate(mstore%varstore(config%ndist, config%ncat))
    allocate(mstore%varistore(gdata%nloci))
    allocate(mstore%varustore(gdata%nloci))
    allocate(mstore%mu_vare_store(4))
    allocate(mstore%indiststore(gdata%nloci, config%ndist))

    ! Initialize test values
    gdata%trains = 0
    gdata%why = [1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0]
    mstate%nnind = 5.0d0
    mstate%vara = 1.0d0
    mstate%vare = 1.0d0
    mstate%mu = 0.0d0
    mstate%gpin = [0.0d0, 0.0001d0, 0.001d0, 0.01d0]
    mstate%delta = [1.0d0, 1.0d0, 1.0d0, 1.0d0]

    ! Fill X with simple test pattern
    do i = 1, gdata%nt
        do j = 1, gdata%nloci
            gdata%X(i, j) = dble(i + j) / 10.0d0
        end do
    end do

    ! Fill C with test pattern (alternating annotations)
    gdata%C = 0
    do j = 1, gdata%nloci
        gdata%C(j, mod(j-1, config%ncat) + 1) = 1
    end do

    ! Test mcmc_init_common
    call mcmc_init_common(config, gdata, mstore)

    print *, "xpx:"
    do j = 1, gdata%nloci
        print "(E25.16)", gdata%xpx(j)
    end do

    print *, "nannot:"
    do j = 1, gdata%nloci
        print "(I4)", gdata%nannot(j)
    end do

    ! Test mcmc_start_values_common
    call mcmc_start_values_common(config, gdata, mstate)

    print *, "yhat:"
    print "(E25.16)", mstate%yhat

    print *, "vary:"
    print "(E25.16)", mstate%vary

    print *, "g_init:"
    do j = 1, gdata%nloci
        print "(E25.16)", mstate%g(j)
    end do

    print *, "p_init:"
    do j = 1, config%ncat
        do i = 1, config%ndist
            print "(E25.16)", mstate%p(i, j)
        end do
    end do

    ! Test mcmc_iteration_pre_common
    mstate%rep = 1
    call mcmc_iteration_pre_common(config, mstate)

    print *, "mu_after_pre:"
    print "(E25.16)", mstate%mu

    print *, "vare_after_pre:"
    print "(E25.16)", mstate%vare

    print *, "log_gp:"
    do i = 2, config%ndist
        print "(E25.16)", mstate%log_gp(i)
    end do

    print *, "vare_gp:"
    do i = 2, config%ndist
        print "(E25.16)", mstate%vare_gp(i)
    end do

end program test_mcmc_utils
