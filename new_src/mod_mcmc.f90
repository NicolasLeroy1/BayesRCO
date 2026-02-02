module mod_mcmc
    use mod_defs, only: dp
    use mod_data, only: ModelConfig, GenomicData, MCMCState, MCMCStorage
    ! We import the specific modules to access their Init and Kernel routines
    use mod_mcmc_mixture
    use mod_mcmc_additive
    use mod_mcmc_bayesCpi
    use mod_mcmc_utils
    use mod_stats, only: compute_dgv, permutate, compute_residuals
    use mod_io, only: output_model, write_dgv
    implicit none

contains

    subroutine run_mcmc(config, gdata, mstate, mstore)
        type(ModelConfig), intent(inout) :: config ! inout because annotflag etc might be modified in kernels
        type(GenomicData), intent(inout) :: gdata
        type(MCMCState), intent(inout) :: mstate
        type(MCMCStorage), intent(inout) :: mstore
        integer :: i, j, jj

        ! 1. Common Initialization
        call mcmc_init_common(config, gdata, mstore)

        ! 2. Specific Initialization
        if (config%nobayesCpi) then
            if (config%mixture) then
                call mcmc_mixture_init(config, gdata, mstate)
            else
                call mcmc_additive_init(config, gdata, mstate, mstore)
            end if
        else
            call mcmc_bayesCpi_init(config, gdata, mstate, mstore)
        end if

        ! 3. Start Values & Residuals
        call mcmc_start_values_common(config, gdata, mstate)
        if (config%nobayesCpi .and. .not. config%mixture) then
             ! Additive specific extra init: permannot
             do j = 1, config%ncat
                 gdata%permannot(j) = j
             end do
        end if
        call compute_residuals(gdata, mstate)

        ! 4. MCMC Loop
        do i = 1, config%numit
             mstate%rep = i
             
             ! Pre-iteration updates (vara, vare, mu, etc.)
             call mcmc_iteration_pre_common(config, mstate)
             
             if (config%permute) then
                 call permutate(gdata%permvec, gdata%nloci)
             end if
             
             ! Additive specific: permute annotation order
             if (config%nobayesCpi .and. .not. config%mixture) then
                 call permutate(gdata%permannot, config%ncat)
             end if

             ! Kernel Execution
             if (config%nobayesCpi) then
                if (config%mixture) then
                    call mcmc_mixture_kernel(config, gdata, mstate, mstore)
                else
                    call mcmc_additive_kernel(config, gdata, mstate, mstore)
                end if
            else
                call mcmc_bayesCpi_kernel(config, gdata, mstate, mstore)
            end if

            ! Update Hypers (common)
            call mcmc_update_hypers_common(config%ndist, config, mstate)

            ! Saving Samples
            if (mod(mstate%rep, config%thin) == 0) then
                if (mstate%rep > config%burnin) then
                    
                    ! Specific Post-Save Logic: Update indiststore / annotstore
                    ! This logic was slightly different for each, but can be unified or checked
                    ! Additive/BayesCpi/Mixture all seem to update indiststore and annotstore similarly
                    ! but Mixture uses vsnptrack and a(i), others use snptracker loop.
                    
                    if (config%nobayesCpi .and. config%mixture) then
                         ! Mixture specific saving logic
                         do j = 1, gdata%nloci
                            jj = gdata%vsnptrack(j)
                            mstore%indiststore(j, jj) = mstore%indiststore(j, jj) + 1.0d0
                            jj = gdata%a(j) ! Re-using jj for annotation index
                            mstore%annotstore(j, jj) = mstore%annotstore(j, jj) + 1.0d0
                        end do
                        call mcmc_save_samples_common('(i10,1x,i10,1x,2(E15.7,1x),100(i10,1x))', '(100E15.7,1x)', &
                            config, gdata, mstate, mstore)
                        
                    elseif (config%nobayesCpi .and. .not. config%mixture) then
                        ! Additive specific saving logic
                        do j = 1, gdata%nloci
                            do jj = 1, config%ncat
                                if (gdata%snptracker(j, jj) > 0) then
                                    mstore%indiststore(j, gdata%snptracker(j, jj)) = &
                                        mstore%indiststore(j, gdata%snptracker(j, jj)) + 1.0d0
                                end if
                            end do
                        end do
                        call mcmc_save_samples_common('(i10,1x,i10,1x,2(E15.7,1x),100(i10,1x))', '(100E15.7,1x)', &
                            config, gdata, mstate, mstore)
                        
                    else 
                        ! BayesCpi specific saving logic
                        do j = 1, gdata%nloci
                            do jj = 1, config%ncat
                                if (gdata%snptracker(j, jj) > 0) then
                                    mstore%indiststore(j, gdata%snptracker(j, jj)) = &
                                        mstore%indiststore(j, gdata%snptracker(j, jj)) + 1.0d0
                                end if
                            end do
                        end do
                        call mcmc_save_samples_common('(i10,1x,i10,1x,2(E15.7,1x),20(i10,1x))', '(20E15.7,1x)', &
                            config, gdata, mstate, mstore)
                    end if
                end if
            end if
            
            ! recalibrate residuals (if enabled)
            if (mod(mstate%rep, 1000) == 0) then
                ! call compute_residuals
            end if

        end do

        ! 5. Post-Processing
        call mcmc_calculate_posterior_means(gdata, mstate, mstore)
        
        ! Specific Post-Processing
        if (config%nobayesCpi .and. .not. config%mixture) then
             do i = 1, gdata%nloci
                 if (gdata%nannot(i) > 1) then
                     mstore%indiststore(i, :) = mstore%indiststore(i, :) / gdata%nannot(i)
                 end if
            end do
        elseif (.not. config%nobayesCpi) then
             do i = 1, gdata%nloci
                 if (gdata%nannot(i) > 1) then
                     mstore%indiststore(i, :) = mstore%indiststore(i, :) / gdata%nannot(i)
                 end if
            end do
        end if
        
        call output_model(config, gdata, mstore)
        mstate%mu = mstore%mu_vare_store(1)
        call compute_dgv(gdata, mstate, mstore)
        call write_dgv(config, gdata)
        
    end subroutine run_mcmc

end module mod_mcmc
