module mod_mcmc_mixture
    use mod_defs, only: dp, LOG_UPPER_LIMIT
    use mod_data
    use mod_random
    use mod_io
    use mod_stats
    use mod_mcmc_utils
    implicit none

contains

    subroutine mcmc_mixture()
        integer :: i, j, k, kk, jj, snploc, l, aa, rep
        real(dp) :: r, skk, sk, ssculm, maxs, maxtemp, detV, clike
        logical :: overflow
        character(len=10) :: ci, ca, cj

        print *, 'mixture'
        print *, 'mixture'
        call mcmc_init_common()
        gdata%vsnptrack = 2
        
        call mcmc_start_values_common()
        
        ! Fill the annotation vector for 1 annotation SNP
        do k = 1, gdata%nloci
            if (gdata%nannot(k) == 1) then
                do j = 1, config%ncat
                    if (gdata%C(k, j) == 1) then
                        gdata%a(k) = j
                    end if
                end do
            end if
        end do

        each_cycle: do rep = 1, config%numit
            mstate%rep = rep
            call mcmc_iteration_pre_common()
            if (config%permute) then
                call permutate(gdata%permvec, gdata%nloci)
            end if
            gdata%snptracker = 0
            ! Choose the annotation of each SNP for the iteration
            do k = 1, gdata%nloci
                snploc = gdata%permvec(k)
                if (gdata%nannot(snploc) > 1) then
                    z => gdata%X(:, snploc)
                    mstate%zz = gdata%xpx(snploc)
                    mstate%zz_vare = mstate%zz / mstate%vare
                    mstate%gk = mstate%g(snploc)
                    gdata%atemp = 0
                    if (mstate%rep /= 1) then
                        gdata%atemp(gdata%a(snploc)) = 1
                    end if
                    mstate%dira = dble(gdata%C(snploc, 1:config%ncat)) + dble(gdata%atemp)
                    mstate%pia = rdirichlet2(config%ncat, mstate%dira)
                    if (gdata%vsnptrack(snploc) > 1) then
                        mstate%ytemp = mstate%yadj + z * mstate%gk
                    end if
                    mstate%rhs = dot_product(mstate%ytemp, z)
                    mstate%lhs = mstate%zz / mstate%vare
                    mstate%ss = 0.0d0
                    mstate%maxs = 0.0d0
                    do i = 2, config%ndist
                        mstate%uhat = mstate%rhs / (mstate%zz + mstate%vare_gp(i))
                        mstate%maxtemp = 0.5d0 * mstate%uhat * mstate%rhs / mstate%vare
                        if (mstate%maxtemp > mstate%maxs) then
                            mstate%maxs = mstate%maxtemp
                        end if
                    end do
                    do j = 1, config%ncat
                        if (gdata%C(snploc, j) == 1) then
                            mstate%ss(j) = mstate%p(1, j) * dexp(-mstate%maxs)
                            do kk = 2, config%ndist
                                mstate%detV = mstate%gp(kk) * mstate%zz_vare + 1.0d0
                                mstate%uhat = mstate%rhs / (mstate%zz + mstate%vare_gp(kk))
                                mstate%ss(j) = mstate%ss(j) + mstate%p(kk, j) * mstate%detV**(-0.5d0) * dexp(0.5d0 * mstate%uhat * mstate%rhs / mstate%vare - mstate%maxs)
                            end do
                            mstate%ss(j) = dlog(mstate%pia(j)) + dlog(mstate%ss(j))
                        end if
                    end do
                    mstate%sstemp = 0.0d0
                    do kk = 1, config%ncat
                        skk = 0.0d0
                        if (gdata%C(snploc, kk) == 1) then
                            skk = mstate%ss(kk)
                            sk = 0.0d0
                            overflow = .false.
                            do l = 1, config%ncat
                                if (gdata%C(snploc, l) == 1) then
                                    if (l == kk) cycle
                                    clike = mstate%ss(l) - skk
                                    if (clike < -LOG_UPPER_LIMIT) then ! underflow
                                        cycle
                                    else if (clike > LOG_UPPER_LIMIT) then
                                        overflow = .true.
                                        exit
                                    end if
                                    sk = sk + dexp(clike)
                                end if
                            end do
                            if (overflow) then
                                mstate%sstemp(kk) = 0.0d0
                            else
                                mstate%sstemp(kk) = 1.0d0 / (1.0d0 + sk)
                            end if
                        end if
                    end do
                    ssculm = 0.0d0
                    call random_number(r)
                    config%annotflag = 1
                    do kk = 1, config%ncat
                        if (gdata%C(snploc, kk) == 1) then
                            ssculm = ssculm + mstate%sstemp(kk)
                            if (r < ssculm) then
                                config%annotflag = kk
                                exit
                            end if
                        end if
                    end do
                    gdata%a(snploc) = config%annotflag
                end if
            end do

            do k = 1, gdata%nloci
                snploc = gdata%permvec(k)
                j = gdata%a(snploc)
                z => gdata%X(:, snploc)
                mstate%zz = gdata%xpx(snploc)
                mstate%zz_vare = mstate%zz / mstate%vare
                mstate%gk = mstate%g(snploc)
                if (gdata%vsnptrack(snploc) > 1) then
                    mstate%yadj = mstate%yadj + z * mstate%gk
                end if
                mstate%rhs = dot_product(mstate%yadj, z)
                mstate%lhs = mstate%zz / mstate%vare
                mstate%s(1) = mstate%log_p(1, j)
                do kk = 2, config%ndist
                    mstate%logdetV = dlog(mstate%gp(kk) * mstate%zz_vare + 1.0d0)
                    mstate%uhat = mstate%rhs / (mstate%zz + mstate%vare_gp(kk))
                    mstate%s(kk) = -0.5d0 * (mstate%logdetV - (mstate%rhs * mstate%uhat / mstate%vare)) + mstate%log_p(kk, j)
                end do
                mstate%stemp = 0.0d0
                do kk = 1, config%ndist
                    skk = mstate%s(kk)
                    sk = 0.0d0
                    overflow = .false.
                    do l = 1, config%ndist
                        if (l == kk) cycle
                        clike = mstate%s(l) - skk
                        if (clike < -LOG_UPPER_LIMIT) then ! underflow
                            cycle
                        else if (clike > LOG_UPPER_LIMIT) then
                            overflow = .true.
                            exit
                        end if
                        sk = sk + dexp(clike)
                    end do
                    if (overflow) then
                        mstate%stemp(kk) = 0.0d0
                    else
                        mstate%stemp(kk) = 1.0d0 / (1.0d0 + sk)
                    end if
                end do
                ssculm = 0.0d0
                call random_number(r)
                config%indistflag = 1
                do kk = 1, config%ndist
                    ssculm = ssculm + mstate%stemp(kk)
                    if (r < ssculm) then
                        config%indistflag = kk
                        exit
                    end if
                end do
                gdata%snptracker(snploc, j) = config%indistflag
                gdata%vsnptrack(snploc) = config%indistflag
                if (config%indistflag == 1) then
                    mstate%gk = 0.0d0
                else
                    mstate%v1 = mstate%zz + mstate%vare / mstate%gp(config%indistflag)
                    mstate%gk = rand_normal(mstate%rhs / mstate%v1, dsqrt(mstate%vare / mstate%v1))
                    mstate%yadj = mstate%yadj - z * mstate%gk
                    mstate%included = mstate%included + 1
                end if
                mstate%g(snploc) = mstate%gk
                if (config%msize > 0 .and. mstate%rep > config%mrep) then
                    if (mstate%included >= config%msize) exit
                end if
            end do ! each loci

            do j = 1, config%ncat
                do i = 1, config%ndist
                    mstate%snpindist(i, j) = count(gdata%snptracker(:, j) == i)
                    mstate%varindist(i, j) = sum(mstate%g * mstate%g, mask=gdata%snptracker(:, j) == i)
                end do
            end do
            mstate%included = gdata%nloci - sum(mstate%snpindist(1, :))

            call mcmc_update_hypers_common(config%ndist)

            if (mod(mstate%rep, config%thin) == 0) then
                if (mstate%rep > config%burnin) then
                    do i = 1, gdata%nloci
                        jj = gdata%vsnptrack(i)
                        mstore%indiststore(i, jj) = mstore%indiststore(i, jj) + 1.0d0
                        aa = gdata%a(i)
                        mstore%annotstore(i, aa) = mstore%annotstore(i, aa) + 1.0d0
                    end do
                    call mcmc_save_samples_common('(i10,1x,i10,1x,2(E15.7,1x),100(i10,1x))', '(100E15.7,1x)')
                end if
            end if
            ! re-calibrate residuals
            if (mod(mstate%rep, 1000) == 0) then
                ! call compute_residuals
            end if
        end do each_cycle

        call mcmc_calculate_posterior_means()
        call output_model()
        mstate%mu = mstore%mu_vare_store(1)
        call compute_dgv()
        call write_dgv()
    end subroutine mcmc_mixture

end module mod_mcmc_mixture
