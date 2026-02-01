module mod_stats
    use mod_defs, only: dp, MISSING_VALUE
    use mod_data
    implicit none

contains

    subroutine permutate(p, n)
        integer, intent(in) :: n
        integer, intent(inout) :: p(n)
        integer :: k, j, i, ipj, itemp, m
        real(dp) :: u(100)
        do i = 1, n
            p(i) = i
        end do
        ! generate up to 100 u(0,1) numbers at a time.
        do i = 1, n, 100
            m = min(n - i + 1, 100)
            call random_number(u)
            do j = 1, m
                ipj = i + j - 1
                k = int(u(j) * (n - ipj + 1)) + ipj
                itemp = p(ipj)
                p(ipj) = p(k)
                p(k) = itemp
            end do
        end do
    end subroutine permutate

    subroutine compute_dgv()
        integer :: i, tr
        tr = 0
        gdata%pred = MISSING_VALUE
        do i = 1, gdata%nind
            if (gdata%trains(i) == 0) then
                tr = tr + 1
                gdata%pred(i) = mstate%mu + dot_product(gdata%X(tr, 1:gdata%nloci), mstore%gstore(1:gdata%nloci))
            end if
        end do
    end subroutine compute_dgv

    subroutine compute_residuals()
        integer :: tr, i
        tr = 0
        do i = 1, gdata%nind
            if (gdata%trains(i) == 0) then
                tr = tr + 1
                mstate%yadj(tr) = gdata%why(i) - dot_product(gdata%X(tr, 1:gdata%nloci), mstate%g) - mstate%mu
            end if
        end do
    end subroutine compute_residuals

end module mod_stats
