program bayesR
    use mod_defs, only: dp
    use mod_data
    use mod_cmd
    use mod_io
    use mod_stats
    use mod_mcmc
    use mod_standardize
    implicit none

    integer :: i, j
    character(len=8)  :: cdate
    character(len=10) :: ctime, ci, ca, cj

    call date_and_time(date=cdate, time=ctime)
    call parse()

    call get_size(config, gdata)
    call load_phenos_plink(config, gdata)
    call allocate_data(config, gdata, mstate, mstore)
    call parse_priors()
    call load_categories(config, gdata)

    if (config%mcmc) then
        open(newunit=config%unit_log, file=config%logfil, status='unknown', form='formatted')
        write(config%unit_log, 901) 'Program BayesRCO'
        write(config%unit_log, 908) 'Run started at', cdate(1:4), cdate(5:6), cdate(7:8), ctime(1:2), ctime(3:4), ctime(5:6)
        write(config%unit_log, 902) 'Prefix for input files', trim(config%inprefix)
        write(config%unit_log, 902) 'Prefix for output files', trim(config%outprefix)
        write(config%unit_log, 903) 'Phenotype column', config%trait_pos
        write(config%unit_log, 903) 'No. of loci', gdata%nloci
        write(config%unit_log, 903) 'No. of individuals', gdata%nind
        write(config%unit_log, 903) 'No. of training individuals', gdata%nt
        write(config%unit_log, 906) 'Prior Vara', mstate%vara, config%dfvara
        write(config%unit_log, 906) 'Prior Vare', mstate%vare, config%dfvare
        write(config%unit_log, 903) 'Model size', config%msize
        write(config%unit_log, 903) 'No. of cycles', config%numit
        write(config%unit_log, 903) 'Burnin ', config%burnin
        write(config%unit_log, 903) 'Thinning rate', config%thin
        write(config%unit_log, 903) 'No. of mixtures', config%ndist
        write(config%unit_log, 905) 'Variance of dist ', mstate%gpin
        write(config%unit_log, 905) 'Dirichlet prior', mstate%delta
        write(config%unit_log, 903) 'Seed ', config%seed1
        write(config%unit_log, 909) 'SNP output ', config%snpout
        write(config%unit_log, 909) 'Cat output', config%cat
        write(config%unit_log, 909) 'Beta output', config%beta
        write(config%unit_log, 903) 'No. of SNP categories ', config%ncat
        call flush(config%unit_log)
    end if

    call load_snp_binary(config, gdata)
    call xcenter(config, gdata)
    call init_random_seed(config)

    if (config%mcmc) then
        mstate%nnind = dble(gdata%nt)
        if (config%snpout) then
            open(newunit=config%unit_loc, file=config%locfil, status='unknown', action='write')
        end if
        if (config%cat) then
            open(newunit=config%unit_cat, file=config%catfil, status='unknown', action='write')
        end if
        if (config%beta) then
            open(newunit=config%unit_beta, file=config%betafil, status='unknown', action='write')
        end if
        open(newunit=config%unit_hyp, file=config%hypfil, status='unknown', form='formatted')
        write(config%unit_hyp, '(2(A10,1x),2(A12,1x),A7)', advance='no') 'Replicate', 'Nsnp', 'Va', 'Ve', ' '
        do j = 1, config%ncat
            write(cj, '(I8)') j
            cj = adjustl(cj)
            do i = 1, config%ndist
                write(ci, '(I8)') i
                ci = adjustl(ci)
                ca = "Nk"//trim(ci)//"_"//trim(cj)
                write(config%unit_hyp, '(A10,1x)', advance="no") ca
            end do
        end do
        do j = 1, config%ncat
            write(cj, '(I8)') j
            cj = adjustl(cj)
            do i = 1, config%ndist
                write(ci, '(I8)') i
                ci = adjustl(ci)
                ca = "Vk"//trim(ci)//"_"//trim(cj)
                write(config%unit_hyp, '(A12)', advance="no") ca
            end do
        end do
        write(config%unit_hyp, *)
        
        if (config%dfvara < -2.0d0) then
            config%VCE = .false.
            mstate%yhat = sum(gdata%why, mask=gdata%trains == 0) / mstate%nnind
            mstate%vary = sum((gdata%why - mstate%yhat) * (gdata%why - mstate%yhat), mask=gdata%trains == 0) / &
                (mstate%nnind - 1.0d0)
            mstate%vara = mstate%vara * mstate%vary
        else
            config%VCE = .true.
            config%vara_ap = mstate%vara
            config%vare_ap = mstate%vare
            if (config%dfvara == -2.0d0) config%vara_ap = 0.0d0
            if (config%dfvare == -2.0d0) config%vare_ap = 0.0d0
        end if

        call run_mcmc(config, gdata, mstate, mstore)

    else
        open(newunit=config%unit_log, file=config%logfil, status='unknown', form='formatted')
        write(config%unit_log, 901) 'Program BayesR'
        write(config%unit_log, 908) 'Run started at', cdate(1:4), cdate(5:6), cdate(7:8), ctime(1:2), ctime(3:4), ctime(5:6)
        write(config%unit_log, 902) 'Prefix for input files', trim(config%inprefix)
        write(config%unit_log, 902) 'Prefix for output files', trim(config%outprefix)
        write(config%unit_log, 903) 'Phenotype column', config%trait_pos
        write(config%unit_log, 903) 'No. of loci', gdata%nloci
        write(config%unit_log, 903) 'No. of individuals', gdata%nind
        write(config%unit_log, 903) 'No. of individuals to predict', gdata%nt
        call load_param(config, gdata, mstore, mstate)
        call compute_dgv(gdata, mstate, mstore)
        call write_dgv(config, gdata)
    end if

    call date_and_time(date=cdate, time=ctime)
    write(config%unit_log, 908) 'Run ended at', cdate(1:4), cdate(5:6), cdate(7:8), ctime(1:2), ctime(3:4), ctime(5:6)
    close(config%unit_log)

901 format(a)
902 format(a, t30, ': ', a)
903 format(a, t30, '= ', i8)
905 format(a, t30, '= ', 10f10.5)
906 format(a, t30, '= ', 2f10.6)
908 format(a20, 1x, a4, '-', a2, '-', a2, ' ', a2, ':', a2, ':', a2)
909 format(a, t30, '= ', l)

end program bayesR
