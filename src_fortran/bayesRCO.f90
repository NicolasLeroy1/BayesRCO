program bayesR
    use mod_types
    use mod_io
    use mod_mcmc
    implicit none

    type(BayesConfig) :: config
    type(BayesData) :: data
    type(BayesModel) :: model
    integer :: narg, i

    ! Parse minimal arguments (simplified driver)
    narg = command_argument_count()
    i = 1
    do while (i <= narg)
        call process_arg(i, config)
    end do

    ! Initialize config gpin/delta if not set
    if (.not. allocated(config%gpin)) then
        allocate(config%gpin(4))
        config%gpin = [0.0_wp, 0.0001_wp, 0.001_wp, 0.01_wp]
    end if
    if (.not. allocated(config%delta)) then
        allocate(config%delta(size(config%gpin)))
        config%delta = 1.0_wp
    end if

    print *, "BayesRCO Refactored Fortran Version"
    
    call get_data_size(config, data)
    call load_phenotypes(config, data)
    call load_genotypes(config, data)
    call load_categories(config, data)
    call xcenter(data)

    ! Write frequencies to file to match expected output for verification
    block
        integer :: u, j
        open(newunit=u, file=trim(config%outprefix)//'.frq', status='unknown')
        do j = 1, data%nloci
            write(u, '(F10.6)') data%freq(j)
        end do
        close(u)
    end block

    call init_model(config, data, model)
    call run_mcmc_driver(config, data, model)

    print *, "Done."

contains

    subroutine process_arg(idx, cfg)
        integer, intent(inout) :: idx
        type(BayesConfig), intent(inout) :: cfg
        character(len=256) :: arg, val
        
        call get_command_argument(idx, arg)
        if (arg == '-bfile') then
            call get_command_argument(idx+1, cfg%inprefix)
            idx = idx + 2
        else if (arg == '-out') then
            call get_command_argument(idx+1, cfg%outprefix)
            idx = idx + 2
        else if (arg == '-numit') then
            call get_command_argument(idx+1, val); read(val, *) cfg%numit
            idx = idx + 2
        else if (arg == '-burnin') then
            call get_command_argument(idx+1, val); read(val, *) cfg%burnin
            idx = idx + 2
        else if (arg == '-thin') then
            call get_command_argument(idx+1, val); read(val, *) cfg%thin
            idx = idx + 2
        else if (arg == '-seed') then
            call get_command_argument(idx+1, val); read(val, *) cfg%seed
            idx = idx + 2
        else if (arg == '-catfile') then
            call get_command_argument(idx+1, cfg%catRC)
            idx = idx + 2
        else
            idx = idx + 1
        end if
    end subroutine process_arg

end program bayesR
