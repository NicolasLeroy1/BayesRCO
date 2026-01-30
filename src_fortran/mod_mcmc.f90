module mod_mcmc
    use mod_types
    use RDistributions
    implicit none

contains

    subroutine init_model(config, data, model)
        type(BayesConfig), intent(in) :: config
        type(BayesData), intent(in) :: data
        type(BayesModel), intent(inout) :: model
        integer :: i, j

        model%mu = 1.0_wp
        model%vara = config%vara
        model%vare = config%vare
        
        ! Estimate vary
        model%vary = sum((data%why - sum(data%why, mask=data%trains==0)/data%nt)**2, mask=data%trains==0) / (data%nt - 1.0_wp)
        
        if (config%dfvara < -2.0_wp) then
            model%vara = config%vara * model%vary
        end if

        allocate(model%g(data%nloci))
        allocate(model%p(config%ndist, config%ncat))
        allocate(model%gp(config%ndist))
        allocate(model%yadj(data%nt))
        allocate(model%xpx(data%nloci))
        allocate(model%snptracker(data%nloci, config%ncat))
        allocate(model%a(data%nloci))

        model%g = sqrt(model%vara / (0.5_wp * data%nloci))
        model%gp = config%gpin * model%vara
        
        do j = 1, config%ncat
            model%p(1, j) = 0.5_wp
            model%p(2:, j) = 1.0_wp / config%gpin(2:)
            model%p(2:, j) = 0.5_wp * model%p(2:, j) / sum(model%p(2:, j))
        end do

        do i = 1, data%nloci
            model%xpx(i) = dot_product(data%X(:, i), data%X(:, i))
        end do

        ! Initial residuals
        model%yadj = data%why(pack([(i, i=1,size(data%trains))], data%trains==0)) - matmul(data%X, model%g) - model%mu
    end subroutine init_model

    subroutine run_mcmc_driver(config, data, model)
        type(BayesConfig), intent(in) :: config
        type(BayesData), intent(in) :: data
        type(BayesModel), intent(inout) :: model
        integer :: rep, k, i, kk, choice, cat_idx
        real(wp) :: rhs, zz, gk, v1, r, cum
        real(wp) :: s(config%ndist), stemp(config%ndist), vare_gp(config%ndist)
        real(wp) :: logdetV, uhat, max_s
        
        ! Tracking variables for hyperparameter updates
        integer :: mc(config%ndist, config%ncat)
        real(wp) :: gh(config%ndist, config%ncat)
        integer :: unit_hyp
        integer, allocatable :: permvec(:)
        integer :: snploc
        
        ! Scratch variables for BayesRCpi
        integer :: active_cats(100), n_active, i_ac, n_opts, selected_idx
        real(wp) :: log_probs(400), prob_opts(400), sum_s, log_lik
        integer :: map_cat(400), map_dist(400)

        open(newunit=unit_hyp, file=trim(config%outprefix)//'.hyp', status='replace')
        write(unit_hyp, '(A)') "Iteration mu vara vare"

        allocate(permvec(data%nloci))
        do k = 1, data%nloci
            permvec(k) = k
        end do

        do rep = 1, config%numit
            call permutate(permvec)
            
            ! Reset trackers
            mc = 0
            gh = 0.0_wp
            ! 1. Update vare
            model%vare = dot_product(model%yadj, model%yadj) / rand_chi_square(real(data%nt + 3, wp))

            ! 2. Update mu
            model%yadj = model%yadj + model%mu
            model%mu = rand_normal(sum(model%yadj)/data%nt, sqrt(model%vare/data%nt))
            model%yadj = model%yadj - model%mu

            vare_gp(2:) = model%vare / model%gp(2:)

            ! 3. Update SNP effects
            do kk = 1, data%nloci
                k = permvec(kk)
                zz = model%xpx(k)
                gk = model%g(k)
                
                if (gk /= 0.0_wp) then
                    model%yadj = model%yadj + data%X(:, k) * gk
                end if

                rhs = dot_product(model%yadj, data%X(:, k))
                
                ! Identify active categories`
                n_active = 0
                do i = 1, config%ncat
                    if (data%C(k, i) == 1) then
                        n_active = n_active + 1
                        active_cats(n_active) = i
                    end if
                end do
                
                if (n_active == 0) then
                    n_active = 1
                    active_cats(1) = 1 ! Default to first cat if none
                end if

                ! Calculate joint probabilities
                n_opts = 0
                do i_ac = 1, n_active
                    cat_idx = active_cats(i_ac)
                    
                    ! d=1 (Null)
                    n_opts = n_opts + 1
                    log_probs(n_opts) = log(model%p(1, cat_idx))
                    map_cat(n_opts) = cat_idx
                    map_dist(n_opts) = 1

                    ! d > 1
                    do i = 2, config%ndist
                        logdetV = log(model%gp(i) * zz / model%vare + 1.0_wp)
                        uhat = rhs / (zz + vare_gp(i))
                        log_lik = -0.5_wp * (logdetV - (rhs * uhat / model%vare))
                        
                        n_opts = n_opts + 1
                        log_probs(n_opts) = log_lik + log(model%p(i, cat_idx))
                        map_cat(n_opts) = cat_idx
                        map_dist(n_opts) = i
                    end do
                end do

                ! Log-sum-exp
                max_s = log_probs(1)
                do i = 2, n_opts
                    if (log_probs(i) > max_s) max_s = log_probs(i)
                end do

                sum_s = 0.0_wp
                do i = 1, n_opts
                    prob_opts(i) = exp(log_probs(i) - max_s)
                    sum_s = sum_s + prob_opts(i)
                end do

                ! Sample
                call random_number(r)
                r = r * sum_s
                cum = 0.0_wp
                selected_idx = n_opts
                do i = 1, n_opts
                    cum = cum + prob_opts(i)
                    if (r < cum) then
                        selected_idx = i
                        exit
                    end if
                end do

                choice = map_dist(selected_idx)
                cat_idx = map_cat(selected_idx)

                if (choice == 1) then
                    gk = 0.0_wp
                else
                    v1 = zz + vare_gp(choice)
                    gk = rand_normal(rhs/v1, sqrt(model%vare/v1))
                    model%yadj = model%yadj - data%X(:, k) * gk
                    
                    gh(choice, cat_idx) = gh(choice, cat_idx) + gk**2 / config%gpin(choice)
                end if
                
                model%g(k) = gk
                model%snptracker(k, cat_idx) = choice
                mc(choice, cat_idx) = mc(choice, cat_idx) + 1
            end do
            
            ! 4. Update hyperparameters
            call update_vara(config, model, mc, gh)
            call update_pi(config, model, mc)

            write(unit_hyp, '(I0, 1X, F12.6, 1X, F12.6, 1X, F12.6)') rep, model%mu, model%vara, model%vare

            if (mod(rep, 10) == 0 .or. rep == config%numit) then
                print '(A,I0,A,I0,A,F10.6)', "Iteration ", rep, ": included ", count(model%g /= 0.0_wp), ", vare ", model%vare
            end if
        end do
        close(unit_hyp)
    end subroutine run_mcmc_driver

    subroutine update_vara(config, model, mc, gh)
        type(BayesConfig), intent(in) :: config
        type(BayesModel), intent(inout) :: model
        integer, intent(in) :: mc(config%ndist, config%ncat)
        real(wp), intent(in) :: gh(config%ndist, config%ncat)
        
        real(wp) :: total_gh, shape, scale, ans
        integer :: total_mc, j
        
        ! If dfvara is negative (but not -2 which might mean something else in some contexts, 
        ! assuming simple flag for now), we might skip. 
        ! But standard BayesR usually updates if dfvara is provided or uses flat prior.
        ! Here we assume we always update unless specifically fixed.
        
        total_mc = sum(mc(2:, :)) ! Count SNPs in non-zero variance classes
        total_gh = sum(gh(2:, :)) ! Sum of g^2 / sigma_k^2
        
        shape = 0.5_wp * (total_mc + config%dfvara)
        ! If dfvara is "flat" (e.g. 4.0 scale 0 ?) 
        ! Using the logic: prior scale * prior df + sample sum of squares
        scale = 0.5_wp * (config%vara * config%dfvara + total_gh)
        
        ! Handle flat prior case (dfvara < 0 usually means flat in some implementations, 
        ! or dfvara -2 means infinite variance? 
        ! Let's follow a standard update: 
        ! If prior is uninformative (dfvara -> 0), shape ~ mc/2, scale ~ gh/2
        
        if (config%dfvara > 0.0_wp) then
             ans = rand_inverse_gamma(shape, scale) 
             ! Note: rand_inverse_gamma(shape, scale) in RDistributions takes (shape, scale)
             ! but in math 1/gamma(shape, rate). 
             ! Check RDistributions definition: returns 1 / rand_gamma(shape, 1/scale)
             ! So scale passed is indeed the scale parameter (beta) of Inverse Gamma.
             ! The definition in RDistributions is: ans = 1.0_wp / rand_gamma(shape, 1.0_wp / scale)
             ! This implies scale argument IS the scale parameter.
        else
             ! Flat prior
             shape = 0.5_wp * total_mc
             scale = 0.5_wp * total_gh
             if (shape > 0.0_wp) then
                ans = rand_inverse_gamma(shape, scale)
             else
                ans = model%vara ! Keep if no info
             end if
        end if
        
        model%vara = ans
        
        ! Update dependent variables
        model%gp = config%gpin * model%vara
        
    end subroutine update_vara

    subroutine update_pi(config, model, mc)
        type(BayesConfig), intent(in) :: config
        type(BayesModel), intent(inout) :: model
        integer, intent(in) :: mc(config%ndist, config%ncat)
        
        integer :: j, i
        real(wp) :: alpha(config%ndist), res(config%ndist)
        
        ! Update pi separately for each category (RCO)
        ! Or jointly if model doesn't support RCO fully yet?
        ! The struct has p(ndist, ncat). So we can update per cat.
        
        do j = 1, config%ncat
            ! Alpha = prior (delta) + counts
            ! Using config%delta (assumed scalar or array). 
            ! mod_types says real(wp), allocatable :: delta(:)
            
            do i = 1, config%ndist
                alpha(i) = config%delta(i) + mc(i, j)
            end do
            
            res = rdirichlet(config%ndist, alpha)
            model%p(:, j) = res
        end do
        
    end subroutine update_pi

    subroutine permutate(v)
        integer, dimension(:), intent(inout) :: v
        integer :: n, i, k, temp
        real :: r
        n = size(v)
        do i = n, 2, -1
            call random_number(r)
            k = int(r * i) + 1
            temp = v(i)
            v(i) = v(k)
            v(k) = temp
        end do
    end subroutine permutate

end module mod_mcmc
