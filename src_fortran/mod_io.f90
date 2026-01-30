module mod_io
    use mod_types
    implicit none

contains

    subroutine get_data_size(config, data)
        type(BayesConfig), intent(in) :: config
        type(BayesData), intent(inout) :: data
        integer :: u, ios
        character(len=1024) :: line

        data%nind = 0
        open(newunit=u, file=trim(config%inprefix)//'.fam', status='old', iostat=ios)
        if (ios /= 0) then
            print *, "Error opening .fam file"
            stop
        end if
        do
            read(u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) > 0) data%nind = data%nind + 1
        end do
        close(u)

        data%nloci = 0
        open(newunit=u, file=trim(config%inprefix)//'.bim', status='old', iostat=ios)
        if (ios /= 0) then
            print *, "Error opening .bim file"
            stop
        end if
        do
            read(u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) > 0) data%nloci = data%nloci + 1
        end do
        close(u)
    end subroutine get_data_size

    subroutine load_phenotypes(config, data)
        type(BayesConfig), intent(in) :: config
        type(BayesData), intent(inout) :: data
        integer :: u, i, ios, n, pos1, pos2
        character(len=8192) :: str
        character(len=64) :: val_str

        allocate(data%why(data%nind), data%trains(data%nind))
        data%nt = 0

        open(newunit=u, file=trim(config%inprefix)//'.fam', status='old')
        do i = 1, data%nind
            read(u, '(a)') str
            pos1 = 1
            n = 0
            do
                pos2 = index(str(pos1:), " ")
                if (pos2 == 0) then
                    n = n + 1
                    val_str = str(pos1:)
                    exit
                end if
                n = n + 1
                val_str = str(pos1:pos1+pos2-2)
                if (n == config%trait_pos + 5) exit
                pos1 = pos1 + pos2
                do while (str(pos1:pos1) == " ")
                    pos1 = pos1 + 1
                end do
            end do
            
            val_str = adjustl(val_str)
            if (trim(val_str) == 'NA') then
                data%trains(i) = 1
                data%why(i) = -999.0_wp
            else
                read(val_str, *) data%why(i)
                data%trains(i) = 0
                data%nt = data%nt + 1
            end if
        end do
        close(u)
    end subroutine load_phenotypes

    subroutine load_genotypes(config, data)
        type(BayesConfig), intent(in) :: config
        type(BayesData), intent(inout) :: data
        integer :: u, i, j, k, tr, ios
        integer(kind=1) :: b1
        real(wp) :: igen(0:3)
        
        igen(0) = 0.0_wp
        igen(1) = 3.0_wp ! missing
        igen(2) = 1.0_wp
        igen(3) = 2.0_wp

        allocate(data%X(data%nt, data%nloci))

        open(newunit=u, file=trim(config%inprefix)//'.bed', status='old', access='stream', form='unformatted')
        read(u) b1 ! magic 1
        read(u) b1 ! magic 2
        read(u) b1 ! mode
        
        do j = 1, data%nloci
            k = 0
            tr = 1
            do i = 1, data%nind
                if (k == 0) read(u) b1
                if (data%trains(i) == 0) then
                    data%X(tr, j) = igen(ibits(b1, k, 2))
                    tr = tr + 1
                end if
                k = mod(k + 2, 8)
            end do
        end do
        close(u)
    end subroutine load_genotypes

    subroutine load_categories(config, data)
        type(BayesConfig), intent(in) :: config
        type(BayesData), intent(inout) :: data
        integer :: u, i, j, ios

        allocate(data%C(data%nloci, config%ncat))
        data%C = 0

        if (len_trim(config%catRC) == 0) return

        open(newunit=u, file=trim(config%catRC), status='old', iostat=ios)
        if (ios /= 0) return

        do i = 1, data%nloci
            read(u, *, iostat=ios) (data%C(i, j), j=1, config%ncat)
            if (ios /= 0) exit
        end do
        close(u)
    end subroutine load_categories

    subroutine xcenter(data)
        type(BayesData), intent(inout) :: data
        integer :: j, nomiss
        real(wp) :: q, mean, sd
        real(wp), allocatable :: xtemp(:)

        allocate(data%freq(data%nloci))
        allocate(xtemp(data%nt))

        do j = 1, data%nloci
            xtemp = data%X(:, j)
            nomiss = count(xtemp < 3.0_wp)
            if (nomiss > 0) then
                q = sum(xtemp, mask=xtemp < 3.0_wp) / (2.0_wp * nomiss)
            else
                q = 0.5_wp
            end if
            data%freq(j) = q

            if (q <= 0.0_wp .or. q >= 1.0_wp) then
                data%X(:, j) = 0.0_wp
            else
                mean = 2.0_wp * q
                sd = sqrt(2.0_wp * q * (1.0_wp - q))
                where (xtemp > 2.0_wp) xtemp = mean
                data%X(:, j) = (xtemp - mean) / sd
            end if
        end do
    end subroutine xcenter

end module mod_io
