module mod_mcmc_additive
    use mod_defs, only: dp, LOG_UPPER_LIMIT
    use mod_data, only: ModelConfig, GenomicData, MCMCState, MCMCStorage
    use mod_random, only: rand_normal
    implicit none

contains

    subroutine mcmc_additive_init(config, gdata, mstate, mstore)
        type(ModelConfig), intent(in) :: config
        type(GenomicData), intent(inout) :: gdata
        type(MCMCState), intent(inout) :: mstate
        type(MCMCStorage), intent(inout) :: mstore
        integer :: j

        print *, 'additive'
        print *, 'additive'
        
        mstore%annotstore = dble(gdata%C(:, 1:config%ncat))
        gdata%gannot = 0.0d0
        do j = 1, config%ncat
            where (gdata%C(:, j) == 1) gdata%snptracker(:, j) = 2
        end do
        
    end subroutine mcmc_additive_init

    subroutine mcmc_additive_kernel(config, gdata, mstate, mstore)
        type(ModelConfig), intent(inout) :: config
        type(GenomicData), intent(inout) :: gdata
        type(MCMCState), intent(inout) :: mstate
        type(MCMCStorage), intent(inout) :: mstore
        
        integer :: i, j, k, kk, jj, snploc, l, kcat
        real(dp) :: r, skk, sk, ssculm, clike
        logical :: overflow
        
        do k = 1, gdata%nloci
            snploc = gdata%permvec(k)
            mstate%gk = mstate%g(snploc)
            ! z => gdata%X(:, snploc) ! Accessed via pointer in original, but direct access might be cleaner or stick to pointer
            ! Avoiding pointer assignment shorthand to be explicit or just use direct access if z is not used elsewhere efficiently
            ! Original: z => gdata%X(:, snploc)
            ! We'll just use gdata%X(:, snploc) directly in dot_product to avoid pointer setup overhead if compiler doesn't optimize,
            ! or keep it if it simplifies code. Let's keep distinct variables for readability.
            
            mstate%zz = gdata%xpx(snploc)
            mstate%zz_vare = mstate%zz / mstate%vare
            if (gdata%vsnptrack(snploc) > 1) then
                mstate%yadj = mstate%yadj + gdata%X(:, snploc) * mstate%gk
            end if
            mstate%rhs = dot_product(mstate%yadj, gdata%X(:, snploc))
            
            do kcat = 1, config%ncat
                j = gdata%permannot(kcat)
                mstate%log_p(1, j) = dlog(mstate%p(1, j))
                do i = 2, config%ndist
                    mstate%log_p(i, j) = dlog(mstate%p(i, j))
                end do
                if (gdata%C(snploc, j) == 1) then
                    ! Re-calculating these seems redundant if they didn't change, but following original logic flow
                    mstate%zz = gdata%xpx(snploc)
                    mstate%zz_vare = mstate%zz / mstate%vare
                    ! Recalculating RHS? It was calculated above.
                    ! In original code:
                    ! z => gdata%X(:, snploc)
                    ! mstate%zz = gdata%xpx(snploc)
                    ! mstate%zz_vare = mstate%zz / mstate%vare
                    ! mstate%rhs = dot_product(mstate%yadj, z)
                    ! This is exactly the same as outside the kcat loop. 
                    ! Optimizing: Use the values calculated before kcat loop.
                    
                    mstate%lhs = mstate%zz / mstate%vare
                    mstate%s(1) = mstate%log_p(1, j)
                    do kk = 2, config%ndist
                        mstate%logdetV = dlog(mstate%gp(kk) * mstate%zz_vare + 1.0d0)
                        mstate%uhat = mstate%rhs / (mstate%zz + mstate%vare_gp(kk))
                        mstate%s(kk) = -0.5d0 * (mstate%logdetV - &
                            (mstate%rhs * mstate%uhat / mstate%vare)) + mstate%log_p(kk, j)
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
                        mstate%yadj = mstate%yadj - gdata%X(:, snploc) * mstate%gk
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
                mstate%varindist(i, j) = sum(gdata%gannot(:, j) * gdata%gannot(:, j), &
                    mask=gdata%snptracker(:, j) == i)
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

    end subroutine mcmc_additive_kernel

end module mod_mcmc_additive
