module mod_mcmc_additive
    use mod_defs, only: dp, LOG_UPPER_LIMIT
    use mod_data
    use mod_random
    use mod_io
    use mod_stats
    use mod_mcmc_utils
    implicit none

contains

    subroutine mcmc_additive()
        integer :: i, j, k, kk, jj, snploc, l, kcat, rep
        real(dp) :: r, skk, sk, ssculm, clike
        logical :: overflow
        character(len=10) :: ci, ca, cj

        print *, 'additive'
        print *, 'additive'
        call mcmc_init_common()
        mstore%annotstore = dble(gdata%C(:, 1:config%ncat))
        gdata%gannot = 0.0d0
        do j = 1, config%ncat
            where (gdata%C(:, j) == 1) gdata%snptracker(:, j) = 2
        end do
        
        call mcmc_start_values_common()
        do kcat = 1, config%ncat
            gdata%permannot(kcat) = kcat
        end do
        call compute_residuals()
        each_cycle2: do rep = 1, config%numit
            mstate%rep = rep
            call mcmc_iteration_pre_common()
            if (config%permute) then
                call permutate(gdata%permvec, gdata%nloci)
            end if
            call permutate(gdata%permannot, config%ncat)
            do k = 1, gdata%nloci
                snploc = gdata%permvec(k)
                mstate%gk = mstate%g(snploc)
                z => gdata%X(:, snploc)
                mstate%zz = gdata%xpx(snploc)
                mstate%zz_vare = mstate%zz / mstate%vare
                if (gdata%vsnptrack(snploc) > 1) then
                    mstate%yadj = mstate%yadj + z * mstate%gk
                end if
                mstate%rhs = dot_product(mstate%yadj, z)
                do kcat = 1, config%ncat
                    j = gdata%permannot(kcat)
                    mstate%log_p(1, j) = dlog(mstate%p(1, j))
                    do i = 2, config%ndist
                        mstate%log_p(i, j) = dlog(mstate%p(i, j))
                    end do
                    if (gdata%C(snploc, j) == 1) then
                        z => gdata%X(:, snploc)
                        mstate%zz = gdata%xpx(snploc)
                        mstate%zz_vare = mstate%zz / mstate%vare
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
                        if (config%indistflag == 1) then
                            mstate%gk = 0.0d0
                        else
                            mstate%v1 = mstate%zz + mstate%vare / mstate%gp(config%indistflag)
                            mstate%gk = rand_normal(mstate%rhs / mstate%v1, dsqrt(mstate%vare / mstate%v1))
                            mstate%yadj = mstate%yadj - z * mstate%gk
                        end if
                        gdata%gannot(snploc, j) = mstate%gk
                        if (config%msize > 0 .and. mstate%rep > config%mrep) then
                            if (mstate%included >= config%msize) exit
                        end if
                    end if ! is the loci in the categorie
                end do ! each categorie
            end do ! each loci

                ! Sum loci effects for each iteration
                mstate%g = sum(gdata%gannot, dim=2)
                gdata%vsnptrack = maxval(gdata%snptracker, dim=2)

                do j = 1, config%ncat
                    do i = 1, config%ndist
                        mstate%snpindist(i, j) = count(gdata%snptracker(:, j) == i)
                        mstate%varindist(i, j) = sum(gdata%gannot(:, j) * gdata%gannot(:, j), mask=gdata%snptracker(:, j) == i)
                    end do
                end do
                ! How many loci included?
                gdata%includedloci = 0.0d0
                do i = 1, gdata%nloci
                    if (gdata%vsnptrack(i) > 1) then
                        gdata%includedloci(i) = 1.0d0
                    end if
                end do
                mstate%included = int(sum(gdata%includedloci))

            call mcmc_update_hypers_common(config%ndist)

            if (mod(mstate%rep, config%thin) == 0) then
                if (mstate%rep > config%burnin) then
                    do i = 1, gdata%nloci
                        do j = 1, config%ncat
                            jj = gdata%snptracker(i, j)
                            if (jj > 0) then
                                mstore%indiststore(i, jj) = mstore%indiststore(i, jj) + 1.0d0
                            end if
                        end do
                    end do
                    call mcmc_save_samples_common('(i10,1x,i10,1x,2(E15.7,1x),100(i10,1x))', '(100E15.7,1x)')
                end if
            end if
            ! re-calibrate residuals
            if (mod(mstate%rep, 1000) == 0) then
                ! call compute_residuals
            end if
        end do each_cycle2

        call mcmc_calculate_posterior_means()
        do i = 1, gdata%nloci
             if (gdata%nannot(i) > 1) then
                 mstore%indiststore(i, :) = mstore%indiststore(i, :) / gdata%nannot(i)
             end if
        end do
        call output_model()
        mstate%mu = mstore%mu_vare_store(1)
        call compute_dgv()
        call write_dgv()
    end subroutine mcmc_additive

end module mod_mcmc_additive
