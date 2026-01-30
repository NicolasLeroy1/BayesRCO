module RDistributions
    use iso_fortran_env, only: wp => real64
    implicit none

    real(wp), parameter :: PI = 3.141592653589793238462_wp

contains

    function rand_uniform(a, b) result(c)
        real(wp), intent(in) :: a, b
        real(wp) :: c, temp
        call random_number(temp)
        c = a + temp * (b - a)
    end function rand_uniform

    function rand_normal(mean, stdev) result(c)
        real(wp), intent(in) :: mean, stdev
        real(wp) :: c, r, theta, temp(2)
        if (stdev <= 0.0_wp) then
            c = mean
        else
            call random_number(temp)
            r = sqrt(-2.0_wp * log(temp(1)))
            theta = 2.0_wp * PI * temp(2)
            c = mean + stdev * r * sin(theta)
        end if
    end function rand_normal

    function rand_exponential(mean) result(c)
        real(wp), intent(in) :: mean
        real(wp) :: c, temp
        if (mean <= 0.0_wp) then
            print *, "mean must be positive"
            c = 0.0_wp
        else
            call random_number(temp)
            c = -mean * log(temp)
        end if
    end function rand_exponential

    recursive function rand_gamma(shape, scale) result(ans)
        real(wp), intent(in) :: shape, scale
        real(wp) :: ans, u, w, d, c, x, xsq, g, v
        
        if (shape >= 1.0_wp) then
            d = shape - 1.0_wp/3.0_wp
            c = 1.0_wp/sqrt(9.0_wp*d)
            do
                x = rand_normal(0.0_wp, 1.0_wp)
                v = 1.0_wp + c*x
                do while (v <= 0.0_wp)
                    x = rand_normal(0.0_wp, 1.0_wp)
                    v = 1.0_wp + c*x
                end do
                v = v*v*v
                call random_number(u)
                xsq = x*x
                if ((u < 1.0_wp - 0.0331_wp*xsq*xsq) .or. (log(u) < 0.5_wp*xsq + d*(1.0_wp - v + log(v)))) then
                    ans = scale * d * v
                    return
                end if
            end do
        else
            g = rand_gamma(shape + 1.0_wp, 1.0_wp)
            call random_number(w)
            ans = scale * g * (w**(1.0_wp/shape))
        end if
    end function rand_gamma

    function rand_chi_square(dof) result(ans)
        real(wp), intent(in) :: dof
        real(wp) :: ans
        ans = rand_gamma(0.5_wp * dof, 2.0_wp)
    end function rand_chi_square

    function rand_inverse_gamma(shape, scale) result(ans)
        real(wp), intent(in) :: shape, scale
        real(wp) :: ans
        ans = 1.0_wp / rand_gamma(shape, 1.0_wp / scale)
    end function rand_inverse_gamma

    function rand_scaled_inverse_chi_square(dof, scale) result(ans)
        real(wp), intent(in) :: dof, scale
        real(wp) :: ans
        ans = rand_inverse_gamma(0.5_wp * dof, 0.5_wp * dof * scale)
    end function rand_scaled_inverse_chi_square

    function rand_beta(a, b) result(ans)
        real(wp), intent(in) :: a, b
        real(wp) :: ans, u, v
        u = rand_gamma(a, 1.0_wp)
        v = rand_gamma(b, 1.0_wp)
        ans = u / (u + v)
    end function rand_beta

    function rdirichlet(n, alpha) result(x)
        integer, intent(in) :: n
        real(wp), intent(in) :: alpha(n)
        real(wp) :: x(n)
        integer :: i
        real(wp) :: sx
        do i = 1, n
            x(i) = rand_gamma(alpha(i), 1.0_wp)
        end do
        sx = sum(x)
        if (sx > 0.0_wp) then
            x = x / sx
        else
            x = 1.0_wp / n
        end if
    end function rdirichlet

end module RDistributions
