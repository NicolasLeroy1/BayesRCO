module mod_mcmc
    use mod_defs, only: dp
    use mod_mcmc_mixture
    use mod_mcmc_additive
    use mod_mcmc_bayesCpi
    implicit none

contains

    subroutine run_mcmc()
        if (config%nobayesCpi) then
            if (config%mixture) then
                call mcmc_mixture()
            else
                call mcmc_additive()
            end if
        else
            call mcmc_bayesCpi()
        end if
    end subroutine run_mcmc

end module mod_mcmc
