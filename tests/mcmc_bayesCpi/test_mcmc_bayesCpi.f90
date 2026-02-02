program test_mcmc_bayesCpi
    use mod_defs, only: dp
    use mod_data, only: ModelConfig, GenomicData, MCMCState, MCMCStorage
    use mod_random
    use mod_mcmc_utils
    use mod_mcmc_bayesCpi
    implicit none

    ! Local declarations
    type(ModelConfig) :: config
    type(GenomicData), target :: gdata
    type(MCMCState), target :: mstate
    type(MCMCStorage) :: mstore

    integer :: i, j, n, k
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
    config%msize = 0

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
    allocate(gdata%vsnptrack(gdata%nloci))
    allocate(gdata%gannot(gdata%nloci, config%ncat))
    allocate(gdata%includedloci(gdata%nloci))

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
    allocate(mstate%s(config%ndist))
    allocate(mstate%stemp(config%ndist))

    allocate(mstore%gstore(gdata%nloci))
    allocate(mstore%pstore(config%ndist, config%ncat))
    allocate(mstore%annotstore(gdata%nloci, config%ncat))

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

    ! Fill C with test pattern
    gdata%C = 0
    do j = 1, gdata%nloci
        gdata%C(j, mod(j-1, config%ncat) + 1) = 1
    end do

    ! Initialize permvec
    do k = 1, gdata%nloci
        gdata%permvec(k) = k
    end do

    ! Compute nannot
    do j = 1, gdata%nloci
        gdata%nannot(j) = sum(gdata%C(j, :))
    end do

    ! Compute xpx
    do j = 1, gdata%nloci
        gdata%xpx(j) = dot_product(gdata%X(:, j), gdata%X(:, j))
    end do

    ! Initialize p
    do j = 1, config%ncat
        mstate%p(1, j) = 0.5d0
        mstate%p(2, j) = 0.25d0
        mstate%p(3, j) = 0.15d0
        mstate%p(4, j) = 0.10d0
    end do

    ! Initialize log_p
    do j = 1, config%ncat
        do i = 1, config%ndist
            mstate%log_p(i, j) = dlog(mstate%p(i, j))
        end do
    end do

    ! Initialize gp and vare_gp
    mstate%gp = mstate%gpin * mstate%vara
    do i = 2, config%ndist
        mstate%log_gp(i) = dlog(mstate%gp(i))
        mstate%vare_gp(i) = mstate%vare / mstate%gp(i)
    end do

    ! Initialize g and gannot
    mstate%g = 0.1d0
    gdata%gannot = 0.0d0

    ! Initialize yadj
    mstate%yadj = gdata%why(1:gdata%nt) - 3.0d0

    ! Run init
    call mcmc_bayesCpi_init(config, gdata, mstate, mstore)

    print *, "gannot_after_init:"
    do j = 1, gdata%nloci
        do i = 1, config%ncat
            print "(F20.16)", gdata%gannot(j, i)
        end do
    end do

    print *, "snptracker_after_init:"
    do j = 1, gdata%nloci
        do i = 1, config%ncat
            print "(I4)", gdata%snptracker(j, i)
        end do
    end do

    ! Run one iteration of bayesCpi kernel
    mstate%rep = 1
    call mcmc_bayesCpi_kernel(config, gdata, mstate, mstore)

    print *, "g_after_kernel:"
    do j = 1, gdata%nloci
        print "(F20.16)", mstate%g(j)
    end do

    print *, "included_after_kernel:"
    print "(I6)", mstate%included

end program test_mcmc_bayesCpi
