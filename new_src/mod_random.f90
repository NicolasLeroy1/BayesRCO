! Bayesian hierarchical models for complex trait analysis using a mixture of 
! normal distributions of SNP effects 
! Copyright (C) 2021 Fanny Mollandin
! Copyright (C) 2014 Gerhard Moser

module mod_random
    use mod_defs, only: dp, PI
    implicit none

contains

    function rand_uniform(a, b) result(c)
        real(dp), intent(in) :: a, b
        real(dp) :: c, temp
        call random_number(temp)
        c = a + temp * (b - a)
    end function rand_uniform

    function rand_normal(mean, stdev) result(c)
        real(dp), intent(in) :: mean, stdev
        real(dp) :: c, r, theta, temp(2)
        if (stdev <= 0.0d0) then
            c = mean
        else
            call random_number(temp)
            r = (-2.0d0 * log(temp(1)))**0.5d0
            theta = 2.0d0 * PI * temp(2)
            c = mean + stdev * r * sin(theta)
        end if
    end function rand_normal

    function rand_exponential(mean) result(c)
        real(dp), intent(in) :: mean
        real(dp) :: c, temp
        if (mean <= 0.0d0) then
            write(*, *) "mean must be positive"
        else
            call random_number(temp)
            c = -mean * log(temp)
        end if
    end function rand_exponential

    recursive function rand_gamma(shape, scale) result(ans)
        real(dp), intent(in) :: shape, scale
        real(dp) :: ans, u, x, xsq, d, c, v, g, w
        if (shape <= 0.0d0) then
            write(*, *) "Shape parameter must be positive"
        end if
        if (scale <= 0.0d0) then
            write(*, *) "Scale parameter must be positive"
        end if

        if (shape >= 1.0d0) then
            d = shape - 1.0d0 / 3.0d0
            c = 1.0d0 / sqrt(9.0d0 * d)
            do
                x = rand_normal(0.0d0, 1.0d0)
                v = 1.0d0 + c * x
                do while (v <= 0.0d0)
                    x = rand_normal(0.0d0, 1.0d0)
                    v = 1.0d0 + c * x
                end do
                v = v * v * v
                call random_number(u)
                xsq = x * x
                if ((u < 1.0d0 - .0331d0 * xsq * xsq) .or. (log(u) < 0.5d0 * xsq + d * (1.0d0 - v + log(v)))) then
                    ans = scale * d * v
                    return
                end if
            end do
        else
            g = rand_gamma(shape + 1.0d0, 1.0d0)
            call random_number(w)
            ans = scale * g * (w**(1.0d0 / shape))
            return
        end if
    end function rand_gamma

    function rand_chi_square(dof) result(ans)
        real(dp), intent(in) :: dof
        real(dp) :: ans
        ans = rand_gamma(0.5d0 * dof, 2.0d0)
    end function rand_chi_square

    function rand_scaled_inverse_chi_square(dof, scale) result(ans)
        real(dp), intent(in) :: dof, scale
        real(dp) :: ans
        ans = rand_inverse_gamma(0.5d0 * dof, 0.5d0 * dof * scale)
    end function rand_scaled_inverse_chi_square

    function rand_inverse_gamma(shape, scale) result(ans)
        real(dp), intent(in) :: shape, scale
        real(dp) :: ans
        ans = 1.0d0 / rand_gamma(shape, 1.0d0 / scale)
    end function rand_inverse_gamma

    function rand_weibull(shape, scale) result(ans)
        real(dp), intent(in) :: shape, scale
        real(dp) :: ans, temp
        if (shape <= 0.0d0) then
            write(*, *) "Shape parameter must be positive"
        end if
        if (scale <= 0.0d0) then
            write(*, *) "Scale parameter must be positive"
        end if
        call random_number(temp)
        ans = scale * (-log(temp))**(1.0d0 / shape)
    end function rand_weibull

    function rand_cauchy(median, scale) result(ans)
        real(dp), intent(in) :: median, scale
        real(dp) :: ans, p
        if (scale <= 0.0d0) then
            write(*, *) "Scale parameter must be positive"
        end if
        call random_number(p)
        ans = median + scale * tan(PI * (p - 0.5d0))
    end function rand_cauchy

    function rand_student_t(dof) result(ans)
        real(dp), intent(in) :: dof
        real(dp) :: ans, y1, y2
        if (dof <= 0.d0) then
            write(*, *) "Degrees of freedom must be positive"
        end if
        y1 = rand_normal(0.0d0, 1.0d0)
        y2 = rand_chi_square(dof)
        ans = y1 / sqrt(y2 / dof)
    end function rand_student_t

    function rand_laplace(mean, scale) result(ans)
        real(dp), intent(in) :: mean, scale
        real(dp) :: ans, u
        if (scale <= 0.0d0) then
            write(*, *) "Scale parameter must be positive"
        end if
        call random_number(u)
        if (u < 0.5d0) then
            ans = mean + scale * log(2.0d0 * u)
        else
            ans = mean - scale * log(2.0d0 * (1.0d0 - u))
        end if
    end function rand_laplace

    function rand_log_normal(mu, sigma) result(ans)
        real(dp), intent(in) :: mu, sigma
        real(dp) :: ans
        ans = exp(rand_normal(mu, sigma))
    end function rand_log_normal

    function rand_beta(a, b) result(ans)
        real(dp), intent(in) :: a, b
        real(dp) :: ans, u, v
        if ((a <= 0.0d0) .or. (b <= 0.0d0)) then
            write(*, *) "Beta parameters must be positive"
        end if
        u = rand_gamma(a, 1.0d0)
        v = rand_gamma(b, 1.0d0)
        ans = u / (u + v)
    end function rand_beta

    function rdirichlet(n, irx) result(x)
        integer, intent(in) :: n
        real(dp), intent(in) :: irx(n)
        real(dp) :: x(n), sx
        integer :: i
        do i = 1, n
            x(i) = rand_gamma(irx(i), 1.0d0)
        end do
        sx = sum(x)
        x = x / sx
    end function rdirichlet

    function rdirichlet2(n, irx) result(x)
        integer, intent(in) :: n
        real(dp), intent(in) :: irx(n)
        real(dp) :: x(n), sx
        integer :: i
        do i = 1, n
            if (irx(i) /= 0) then
                x(i) = rand_gamma(irx(i), 1.0d0)
            else
                x(i) = 0.0d0
            end if
        end do
        sx = sum(x)
        x = x / sx
    end function rdirichlet2

end module mod_random
