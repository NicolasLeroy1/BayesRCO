module mod_mcmc_utils
    use mod_defs, only: dp
    ! use mod_data ! Removed global import
    use mod_data, only: ModelConfig, GenomicData, MCMCState, MCMCStorage
    use mod_io, only: output_beta
    use mod_stats, only: compute_residuals
    implicit none

contains

    subroutine mcmc_save_samples_common(fmt1, fmt2, config, gdata, mstate, mstore)
        character(len=*), intent(in) :: fmt1, fmt2
        type(ModelConfig), intent(in) :: config
        type(GenomicData), intent(in) :: gdata
        type(MCMCState), intent(inout) :: mstate
        type(MCMCStorage), intent(inout) :: mstore
        
        mstate%counter = mstate%counter + 1
        mstore%gstore = mstore%gstore + mstate%g
        mstore%pstore = mstore%pstore + mstate%p
        mstore%varistore = mstore%varistore + mstate%g**2
        mstore%mu_vare_store(1) = mstore%mu_vare_store(1) + mstate%mu
        mstore%mu_vare_store(2) = mstore%mu_vare_store(2) + mstate%included
        mstore%mu_vare_store(3) = mstore%mu_vare_store(3) + mstate%vara
        mstore%mu_vare_store(4) = mstore%mu_vare_store(4) + mstate%vare
        mstore%varstore = mstore%varstore + mstate%varindist
        mstore%snpstore = mstore%snpstore + mstate%snpindist
        if (mstate%counter > 1) then
            mstore%varustore = mstore%varustore + (mstate%counter * mstate%g - mstore%gstore)**2 / &
                (mstate%counter * (mstate%counter - 1))
        end if
        
        write(config%unit_hyp, fmt1, advance='no') mstate%rep, mstate%included, mstate%vara, mstate%vare, mstate%snpindist
        write(config%unit_hyp, fmt2) mstate%varindist
        call flush(config%unit_hyp)
        if (config%beta) call output_beta(config, mstate, gdata)
        ! if(config%cat) call output_cat 
    end subroutine mcmc_save_samples_common

    subroutine mcmc_init_common(config, gdata, mstore)
        type(ModelConfig), intent(in) :: config
        type(GenomicData), intent(inout) :: gdata
        type(MCMCStorage), intent(inout) :: mstore
        integer :: i
        mstore%pstore = 0.0d0
        mstore%gstore = 0.0d0
        mstore%varistore = 0.0d0
        mstore%mu_vare_store = 0.0d0
        mstore%snpstore = 0.0d0
        mstore%indiststore = 0.0d0
        mstore%varstore = 0.0d0
        mstore%varustore = 0.0d0
        gdata%snptracker = 0
        gdata%xpx = 0.0d0
        do i = 1, gdata%nloci
            gdata%xpx(i) = dot_product(gdata%X(:, i), gdata%X(:, i))
        end do
        gdata%nannot = sum(gdata%C(:, 1:config%ncat), dim=2)
    end subroutine mcmc_init_common

    subroutine mcmc_start_values_common(config, gdata, mstate)
        type(ModelConfig), intent(in) :: config
        type(GenomicData), intent(inout) :: gdata
        type(MCMCState), intent(inout) :: mstate
        integer :: j, k
        mstate%mu = 1.0d0
        mstate%yadj = 0.0d0
        mstate%yhat = sum(gdata%why, mask=gdata%trains == 0) / mstate%nnind
        mstate%vary = sum((gdata%why - mstate%yhat) * (gdata%why - mstate%yhat), mask=gdata%trains == 0) / (mstate%nnind - 1.0d0)
        mstate%gp = mstate%gpin * mstate%vara
        mstate%scale = 0.0d0
        do j = 1, config%ncat
            mstate%p(1, j) = 0.5d0
            mstate%p(2:config%ndist, j) = 1.0d0 / mstate%gpin(2:config%ndist)
            mstate%p(2:config%ndist, j) = 0.5 * mstate%p(2:config%ndist, j) / sum(mstate%p(2:config%ndist, j))
        end do
        mstate%g = dsqrt(mstate%vara / (0.5 * dble(gdata%nloci)))
        do k = 1, gdata%nloci
            gdata%permvec(k) = k
        end do
        call compute_residuals(gdata, mstate)
    end subroutine mcmc_start_values_common

    subroutine mcmc_iteration_pre_common(config, mstate)
        use mod_random, only: rand_chi_square, rand_normal
        type(ModelConfig), intent(in) :: config
        type(MCMCState), intent(inout) :: mstate
        integer :: i
        mstate%included = 0
        if (.not. config%VCE) then
            mstate%vare = dot_product(mstate%yadj, mstate%yadj) / rand_chi_square(mstate%nnind + 3.0d0)
        end if
        mstate%yadj = mstate%yadj + mstate%mu
        mstate%mu = rand_normal(sum(mstate%yadj) / mstate%nnind, dsqrt(mstate%vare / mstate%nnind))
        mstate%yadj = mstate%yadj - mstate%mu
        do i = 2, config%ndist
            mstate%log_gp(i) = dlog(mstate%gp(i))
            mstate%vare_gp(i) = mstate%vare / mstate%gp(i)
        end do
    end subroutine mcmc_iteration_pre_common

    subroutine mcmc_update_hypers_common(nc, config, mstate)
        use mod_random, only: rand_scaled_inverse_chi_square, rdirichlet
        integer, intent(in) :: nc
        type(ModelConfig), intent(in) :: config
        type(MCMCState), intent(inout) :: mstate
        integer :: i, j
        if (config%VCE) then
            mstate%scale = (dble(mstate%included) * sum(mstate%g**2) + config%vara_ap * config%dfvara) / &
                (config%dfvara + dble(mstate%included))
            mstate%vara = rand_scaled_inverse_chi_square(dble(mstate%included) + config%dfvara, mstate%scale)
            if (nc == 2) then ! BayesCpi variant
                 mstate%gp(2) = mstate%vara / mstate%included
            else
                 mstate%gp = mstate%gpin * mstate%vara
            end if
            mstate%vare = (dot_product(mstate%yadj, mstate%yadj) + config%vare_ap * config%dfvare) / (mstate%nnind + config%dfvare)
            mstate%vare = rand_scaled_inverse_chi_square(mstate%nnind + config%dfvare, mstate%vare)
        end if

        do j = 1, config%ncat
            mstate%dirx = dble(mstate%snpindist(:, j)) + mstate%delta
            mstate%p(:, j) = rdirichlet(config%ndist, mstate%dirx)
            mstate%log_p(1, j) = dlog(mstate%p(1, j))
            do i = 2, config%ndist
                mstate%log_p(i, j) = dlog(mstate%p(i, j))
            end do
        end do
    end subroutine mcmc_update_hypers_common

    subroutine mcmc_calculate_posterior_means(gdata, mstate, mstore)
        type(GenomicData), intent(in) :: gdata
        type(MCMCState), intent(in) :: mstate
        type(MCMCStorage), intent(inout) :: mstore
        integer :: i
        mstore%gstore = mstore%gstore / mstate%counter
        mstore%pstore = mstore%pstore / mstate%counter
        mstore%mu_vare_store = mstore%mu_vare_store / mstate%counter
        mstore%varstore = mstore%varstore / mstate%counter
        mstore%snpstore = mstore%snpstore / mstate%counter
        mstore%varustore = mstore%varustore / mstate%counter
        mstore%varistore = mstore%varistore / mstate%counter
        do i = 1, gdata%nloci
            mstore%indiststore(i, :) = mstore%indiststore(i, :) / mstate%counter
            mstore%annotstore(i, :) = mstore%annotstore(i, :) / mstate%counter
        end do
    end subroutine mcmc_calculate_posterior_means

end module mod_mcmc_utils
