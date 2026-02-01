module mod_io
    use mod_defs, only: dp, MISSING_VALUE
    use mod_data
    implicit none

contains

    subroutine get_size()
        integer :: unit_phen, unit_bim
        gdata%nind = 0
        open(newunit=unit_phen, file=trim(config%phenfil), status='old', form='formatted')
        do
            read(unit_phen, *, iostat=ios)
            if (ios /= 0) exit
            gdata%nind = gdata%nind + 1
        end do
        close(unit_phen)

        gdata%nloci = 0
        open(newunit=unit_bim, file=trim(config%bimfil), status='old', form='formatted')
        do
            read(unit_bim, *, iostat=ios)
            if (ios /= 0) exit
            gdata%nloci = gdata%nloci + 1
        end do
        close(unit_bim)
    end subroutine get_size

    subroutine load_phenos_plink()
        character(len=1024) :: str
        character(len=20) :: why_str
        integer :: pos1, pos2, n, i, unit_in
 
        allocate(gdata%trains(gdata%nind), gdata%why(gdata%nind))
        open(newunit=unit_in, file=trim(config%phenfil), status='old', form='formatted')
        do i = 1, gdata%nind
            read(unit_in, '(a)') str
            pos1 = 1
            n = 0
            do
                pos2 = index(str(pos1:), " ")
                if (pos2 == 0) then
                    n = n + 1
                    why_str = str(pos1:)
                    exit
                end if
                n = n + 1
                why_str = str(pos1:pos1 + pos2 - 2)
                if (n == config%trait_pos + 5) then
                    exit
                end if
                pos1 = pos2 + pos1
            end do
            why_str = trim(why_str)
            if (why_str /= 'NA') then
                gdata%trains(i) = 0
                read(why_str, *) gdata%why(i)
            else
                gdata%trains(i) = 1
                gdata%why(i) = MISSING_VALUE
            end if
        end do
        close(unit_in)
    end subroutine load_phenos_plink

    subroutine load_snp_binary()
        integer :: i, j, k, tr, unit_gen
        integer(1) :: b1
        real(dp), dimension(0:3) :: igen
        real(dp) :: val
        igen(0) = 0.0d0
        igen(1) = 3.0d0
        igen(2) = 1.0d0
        igen(3) = 2.0d0
 
        open(newunit=unit_gen, file=trim(config%genfil), status='old', access='stream', form='unformatted')
        read(unit_gen) b1
        ! Magic number checks omitted for brevity but should be here if needed
        read(unit_gen) b1
        read(unit_gen) b1
 
        do j = 1, gdata%nloci
            k = 0
            tr = 1
            do i = 1, gdata%nind
                if (k == 0) read(unit_gen) b1
                val = igen(ibits(b1, k, 2))
                if (gdata%trains(i) == 0) then
                    gdata%X(tr, j) = val
                    tr = tr + 1
                end if
                k = mod(k + 2, 8)
            end do
        end do
        close(unit_gen)
    end subroutine load_snp_binary

    subroutine init_random_seed()
        integer :: i, n, clock_local
        integer, dimension(:), allocatable :: seed

        call random_seed(size=n)
        allocate(seed(n))

        if (config%seed1 /= 0) then
            do i = 1, n
                seed(i) = abs(config%seed1) + (i - 1)
            end do
        else
            call system_clock(count=clock_local)
            seed = clock_local + 37 * (/(i - 1, i = 1, n)/)
        end if
        call random_seed(put=seed)
        deallocate(seed)
    end subroutine init_random_seed

    subroutine allocate_data()
        if (.not. config%mcmc) then
            where (gdata%trains == 0) gdata%trains = 3
            where (gdata%trains == 1) gdata%trains = 0
            where (gdata%trains == 3) gdata%trains = 1
        end if
        gdata%nt = count(gdata%trains == 0)
        
        allocate(gdata%pred(gdata%nind), mstate%gpin(config%ndist), mstate%gp(config%ndist), &
            mstate%p(config%ndist, config%ncat), gdata%permannot(config%ncat), &
            gdata%X(gdata%nt, gdata%nloci), mstate%delta(config%ndist), mstate%dirx(config%ndist), &
            mstate%g(gdata%nloci), mstate%yadj(gdata%nt), &
            mstate%snpindist(config%ndist, config%ncat), mstate%varindist(config%ndist, config%ncat), &
            z(gdata%nt), mstate%s(config%ndist), mstate%stemp(config%ndist), mstate%sstemp(config%ncat), &
            gdata%xpx(gdata%nloci), mstore%gstore(gdata%nloci), mstore%snpstore(config%ndist, config%ncat), &
            mstore%varstore(config%ndist, config%ncat), mstore%pstore(config%ndist, config%ncat), &
            mstore%indiststore(gdata%nloci, config%ndist), mstore%mu_vare_store(4), gdata%freqstore(gdata%nloci), &
            gdata%permvec(gdata%nloci), gdata%snptracker(gdata%nloci, config%ncat), &
            mstate%log_p(config%ndist, config%ncat), mstate%log_gp(config%ndist), &
            mstate%vare_gp(config%ndist), mstore%varustore(gdata%nloci), mstore%varistore(gdata%nloci), &
            gdata%vsnptrack(gdata%nloci), gdata%nannot(gdata%nloci), gdata%a(gdata%nloci), mstate%ss(config%ncat), &
            gdata%includedloci(gdata%nloci), mstate%pia(config%ncat), mstate%ytemp(gdata%nloci), &
            mstate%dira(config%ncat), gdata%atemp(config%ncat), mstore%annotstore(gdata%nloci, config%ncat), &
            gdata%gannot(gdata%nloci, config%ncat), stat=ios)
            
        if (ios /= 0) then
            stop 'Unable to allocate required storage for data'
        end if
    end subroutine allocate_data

    subroutine load_param()
        character(len=100) :: str
        character(len=10) :: dum
        integer :: i, nc, ios_local, unit_param, unit_mod
        real(dp), dimension(:), allocatable :: gtemp
 
        open(newunit=unit_param, file=config%paramfil, status='old', form='formatted')
        read(unit_param, '(a)') str
        nc = config%ndist + 1
        allocate(gtemp(nc))
        do i = 1, gdata%nloci
            read(unit_param, *, iostat=ios_local) gtemp(1:nc)
            mstore%gstore(i) = gtemp(nc)
            if (ios_local /= 0) then
                write(*, *) 'Error reading file ', adjustl(config%paramfil)
            end if
        end do
        close(unit_param)
        open(newunit=unit_mod, file=config%modfil, status='old', form='formatted')
        read(unit_mod, *) dum, mstate%mu
        close(unit_mod)
    end subroutine load_param

    subroutine load_categories()
        integer :: i, j, ios_local, unit_cat
        open(newunit=unit_cat, file=trim(config%catRC), iostat=ios_local, status='old', form="formatted")
        if (ios_local /= 0) stop 'Cannot open file! '
        allocate(gdata%C(1:gdata%nloci, 1:config%ncat + 1))
        do i = 1, gdata%nloci
            read(unit_cat, *) (gdata%C(i, j), j = 1, config%ncat)
        end do
        close(unit_cat)
    end subroutine load_categories

    subroutine write_dgv()
        integer :: i, unit_dgv
        character(len=2) :: missvalue
        missvalue = 'NA'
        open(newunit=unit_dgv, file=config%mbvfil, status='unknown', form='formatted')
        do i = 1, gdata%nind
            if (gdata%trains(i) == 0) then
                write(unit_dgv, '(E15.7)') gdata%pred(i)
            else
                write(unit_dgv, '(a2)') missvalue
            end if
        end do
        close(unit_dgv)
    end subroutine write_dgv

    subroutine output_model()
        integer :: i, j, unit_param, unit_mod
        character(len=20) :: ci, ca, cj
 
        open(newunit=unit_param, file=config%paramfil, status='unknown', iostat=ios)
        if (ios /= 0) then
            write(*, *) 'Error opening ', trim(config%paramfil)
            stop
        end if
        do i = 1, config%ndist
            write(ci, '(I8)') i
            ci = adjustl(ci)
            ca = "PIP"//trim(ci)
            write(unit_param, '(2(A5))', advance="no") ' ', ca
        end do
        write(unit_param, '(2(A5))', advance="no") ' ', 'beta'
        do i = 1, config%ncat
            write(ci, '(I8)') i
            ci = adjustl(ci)
            ca = "PAIP"//trim(ci)
            write(unit_param, '(2(A5))', advance="no") ' ', ca
        end do
        write(unit_param, '(2(A5))', advance="no") ' ', 'Vbeta'
        write(unit_param, *) ' ', 'Vi'
        do i = 1, gdata%nloci
            write(unit_param, '(50(E15.7,1X))') mstore%indiststore(i, :), mstore%gstore(i), &
                mstore%annotstore(i, :), mstore%varustore(i), mstore%varistore(i)
        end do
        close(unit_param)
 
        open(newunit=unit_mod, file=config%modfil, status='unknown', iostat=ios)
        if (ios /= 0) then
            write(*, *) 'Error opening ', trim(config%modfil)
            stop
        end if
        write(unit_mod, 800) 'Mean', mstore%mu_vare_store(1)
        write(unit_mod, 800) 'Nsnp', mstore%mu_vare_store(2)
        write(unit_mod, 800) 'Va', mstore%mu_vare_store(3)
        write(unit_mod, 800) 'Ve', mstore%mu_vare_store(4)
        do j = 1, config%ncat
            write(cj, '(I8)') j
            cj = adjustl(cj)
            do i = 1, config%ndist
                write(ci, '(I8)') i
                ci = adjustl(ci)
                ca = "Nk"//trim(adjustl(ci))//'_'//trim(cj)
                write(unit_mod, 800) trim(ca), mstore%snpstore(i, j)
            end do
        end do
        do j = 1, config%ncat
            write(cj, '(I8)') j
            cj = adjustl(cj)
            do i = 1, config%ndist
                write(ci, '(I8)') i
                ci = adjustl(ci)
                ca = "Pk"//trim(adjustl(ci))//'_'//trim(cj)
                write(unit_mod, 800) trim(ca), mstore%pstore(i, j)
            end do
        end do
        do j = 1, config%ncat
            write(cj, '(I8)') j
            cj = adjustl(cj)
            do i = 1, config%ndist
                write(ci, '(I8)') i
                ci = adjustl(ci)
                ca = "Vk"//trim(adjustl(ci))//'_'//trim(cj)
                write(unit_mod, 800) trim(ca), mstore%varstore(i, j)
            end do
        end do
        close(unit_mod)
 800     format(a, t10, E15.7)
    end subroutine output_model

    subroutine output_beta()
        integer :: i
        character(len=30) :: betai
        do i = 1, gdata%nloci
            write(betai, '(e15.6)') mstate%g(i)**2
            betai = adjustl(betai)
            write(config%unit_beta, '(1X,a)', advance='no') trim(betai)
        end do
        write(config%unit_beta, *)
    end subroutine output_beta

end module mod_io
