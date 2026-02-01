module mod_cmd
    use mod_defs, only: dp
    use mod_data
    implicit none

    integer :: narg, nopt
    character(len=1024) :: arg
    character(len=1024), allocatable, dimension(:) :: cmd_line

    integer, parameter :: &
        a_int = 1, &
        a_float = 2, &
        a_char = 3, &
        a_flag = 4

    type param
        integer, allocatable, dimension(:) :: pos
        character(len=20), allocatable, dimension(:) :: key
        character(len=20), allocatable, dimension(:) :: argtype
        character(len=100), allocatable, dimension(:) :: desc
        integer, allocatable, dimension(:) :: kind
        character(len=200), allocatable, dimension(:) :: default
    end type param

    type(param) :: register

contains

    subroutine init_register()
        nopt = 28
        allocate(register%pos(nopt))
        allocate(register%key(nopt))
        allocate(register%argtype(nopt))
        allocate(register%desc(nopt))
        allocate(register%kind(nopt))
        allocate(register%default(nopt))

        call include_option(1, '-bfile', '[prefix]', 'prefix PLINK binary files', a_char, '')
        call include_option(2, '-out', '[prefix]', 'prefix for output', a_char, '')
        call include_option(3, '-n', '[num]', 'phenotype column', a_int, '1')
        call include_option(4, '-vara', '[num]', 'SNP variance prior', a_float, '0.01')
        call include_option(5, '-vare', '[num]', 'error variance prior', a_float, '0.01')
        call include_option(6, '-dfvara', '[num]', 'degrees of freedom Va', a_float, '-2.0')
        call include_option(7, '-dfvare', '[num]', 'degrees of freedom Ve', a_float, '-2.0')
        call include_option(8, '-delta', '[num]', 'prior for Dirichlet', a_float, '1.0')
        call include_option(9, '-msize', '[num]', 'number of SNPs in reduced update', a_int, '0')
        call include_option(10, '-mrep', '[num]', 'number of full cycles in reduced update', a_int, '5000')
        call include_option(11, '-numit', '[num]', 'length of MCMC chain', a_int, '50000')
        call include_option(12, '-burnin', '[num]', 'burnin steps', a_int, '20000')
        call include_option(13, '-thin', '[num]', 'thinning rate', a_int, '10')
        call include_option(14, '-ndist', '[num]', 'number of mixture distributions', a_int, '4')
        call include_option(15, '-gpin', '[num]', 'effect sizes of mixtures (% x Va)', a_float, '0.0,0.0001,0.001,0.01')
        call include_option(16, '-seed', '[num]', 'initial value for random number', a_int, '0')
        call include_option(17, '-predict', '[flag]', 'perform prediction', a_flag, 'f')
        call include_option(18, '-snpout', '[flag]', 'output detailed SNP info', a_flag, 'f')
        call include_option(19, '-permute', '[flag]', 'permute order of SNP', a_flag, 'f')
        call include_option(20, '-model', '[filename]', 'model summary file (for prediction) ', a_char, '')
        call include_option(21, '-freq', '[filename]', 'SNP frequency file (for prediction)', a_char, '')
        call include_option(22, '-param', '[filename]', 'SNP effect file (for prediction)', a_char, '')
        call include_option(23, '-cat', '[flag]', 'output SNP categories per iteration', a_flag, 'f')
        call include_option(24, '-beta', '[flag]', 'output SNP effect per iteration', a_flag, 'f')
        call include_option(25, '-ncat', '[num]', 'number of SNP categories', a_int, '1')
        call include_option(26, '-catfile', '[filename]', 'SNP categories file', a_char, '')
        call include_option(27, '-additive', '[flag]', 'additive annotations', a_flag, 'f')
        call include_option(28, '-bayesCpi', '[flag]', 'run bayesCpi', a_flag, 'f')
    end subroutine init_register

    subroutine include_option(pos, key, argtype, desc, kind, default)
        integer, intent(in) :: pos, kind
        character(len=*), intent(in) :: key, argtype, desc, default
        register%pos(pos) = pos
        register%key(pos) = trim(key)
        register%argtype(pos) = trim(argtype)
        register%desc(pos) = trim(desc)
        register%kind(pos) = kind
        register%default(pos) = trim(default)
    end subroutine include_option

    function str_match(str1, str2) result(match)
        character(*), intent(in) :: str1, str2
        logical :: match
        match = (trim(str1) == trim(str2))
    end function str_match

    function ntokens(line) result(toks)
        character(*), intent(in) :: line
        integer :: toks, i, n
        character(len=1) :: separator
        separator = ','
        toks = 0
        i = 1
        n = len_trim(line)
        do while (i <= n)
            do while (i <= n .and. line(i:i) == separator)
                i = i + 1
            end do
            if (i > n) exit
            toks = toks + 1
            do while (i <= n .and. line(i:i) /= separator)
                i = i + 1
            end do
        end do
    end function ntokens

    subroutine tokenize(str, sep, dim, val)
        character(len=*), intent(in) :: str
        character(len=1), intent(in) :: sep
        integer, intent(in) :: dim
        character(len=128), dimension(dim), intent(out) :: val
        integer :: pos1, pos2, n
        pos1 = 1
        n = 0
        do
            pos2 = index(str(pos1:), sep)
            if (pos2 == 0) then
                n = n + 1
                val(n) = str(pos1:)
                exit
            end if
            n = n + 1
            val(n) = str(pos1:pos1 + pos2 - 2)
            pos1 = pos2 + pos1
        end do
    end subroutine tokenize

    subroutine get_cmdLine()
        integer :: i
        narg = command_argument_count()
        allocate(cmd_line(narg))
        do i = 1, narg
            call get_command_argument(i, cmd_line(i))
        end do
    end subroutine get_cmdLine

    function cast_int(value) result(r)
        character(len=*), intent(in) :: value
        integer :: r
        read(value, *) r
    end function cast_int

    function cast_float(value) result(r)
        character(len=*), intent(in) :: value
        real :: r
        read(value, *) r
    end function cast_float

    function cast_logical(value) result(r)
        character(len=*), intent(in) :: value
        logical :: r
        r = .false.
        if (str_match(trim(value), 't')) r = .true.
    end function cast_logical

    function is_key(key) result(r)
        character(len=*), intent(in) :: key
        logical :: r
        r = .false.
        if (key(1:1) == '-') r = .true.
    end function is_key

    subroutine update_register()
        integer :: i, k, kind
        character(len=1024) :: key
        i = 1
        do while (i <= narg)
            key = trim(cmd_line(i))
            k = 1
            do while (k <= nopt)
                if (str_match(key, trim(register%key(k)))) then
                    kind = register%kind(k)
                    if (kind == a_flag) then
                        register%default(k) = 't'
                        i = i + 1
                    else if (i == narg) then
                        stop 'ERROR: Problem parsing the command line arguments'
                    else if (is_key(trim(cmd_line(i + 1)))) then
                        stop 'ERROR: Problem parsing the command line arguments'
                    else
                        register%default(k) = trim(cmd_line(i + 1))
                        i = i + 2
                    end if
                    exit
                end if
                k = k + 1
            end do
            if (k > nopt) then
                stop 'ERROR: Unknown command line option'
            end if
        end do
    end subroutine update_register

    subroutine parse_help()
        integer :: i, k
        if (narg == 0) then
            print *, 'Usage: bayesR [options]'
            stop
        end if
        do i = 1, narg
            if (str_match(trim(cmd_line(i)), '-h')) then
                print *, 'Usage: bayesR [options]'
                do k = 1, nopt
                    print *, trim(register%key(k)), ' ', trim(register%desc(k))
                end do
                stop
            end if
        end do
    end subroutine parse_help

    subroutine parse()
        call init_register()
        call get_cmdLine()
        call parse_help()
        call update_register()
        call parse_out()
        call parse_method()
        call parse_plink()
        call parse_categorie()
        call parse_ncat()
        call parse_trait_pos()
        call parse_predict()
        call parse_ndist()
        call parse_initialise()
        call parse_bayesCpi()
    end subroutine parse

    subroutine parse_out()
        integer :: i
        do i = 1, nopt
            if (str_match(trim(register%key(i)), '-out')) then
                config%outprefix = trim(register%default(i))
                config%logfil = trim(config%outprefix)//'.log'
                config%freqfil = trim(config%outprefix)//'.frq'
                config%mbvfil = trim(config%outprefix)//'.gv'
                config%hypfil = trim(config%outprefix)//'.hyp'
                config%locfil = trim(config%outprefix)//'.snp'
                config%catfil = trim(config%outprefix)//'.catit'
                config%betafil = trim(config%outprefix)//'.beta'
                config%modfil = trim(config%outprefix)//'.model'
                config%paramfil = trim(config%outprefix)//'.param'
                exit
            end if
        end do
    end subroutine parse_out

    subroutine parse_method()
        integer :: i
        config%mixture = .true.
        do i = 1, nopt
            if (str_match(trim(register%key(i)), '-additive')) then
                if (cast_logical(register%default(i))) config%mixture = .false.
            end if
        end do
    end subroutine parse_method

    subroutine parse_plink()
        integer :: i
        logical :: fileExist
        do i = 1, nopt
            if (str_match(trim(register%key(i)), '-bfile')) then
                config%inprefix = trim(register%default(i))
                config%genfil = trim(config%inprefix)//'.bed'
                config%phenfil = trim(config%inprefix)//'.fam'
                config%bimfil = trim(config%inprefix)//'.bim'
                inquire(file=config%genfil, exist=fileExist)
                if (.not. fileExist) stop 'file not found'
                exit
            end if
        end do
    end subroutine parse_plink

    subroutine parse_categorie()
        integer :: i
        logical :: fileExist
        do i = 1, nopt
            if (str_match(trim(register%key(i)), '-catfile')) then
                config%catRC = trim(register%default(i))
                inquire(file=config%catRC, exist=fileExist)
                if (.not. fileExist) stop 'file not found'
            end if
        end do
    end subroutine parse_categorie

    subroutine parse_trait_pos()
        integer :: i
        do i = 1, nopt
            if (str_match(trim(register%key(i)), '-n')) then
                config%trait_pos = cast_int(register%default(i))
            end if
        end do
    end subroutine parse_trait_pos

    subroutine parse_ncat()
        integer :: i
        do i = 1, nopt
            if (str_match(trim(register%key(i)), '-ncat')) then
                config%ncat = cast_int(register%default(i))
            end if
        end do
    end subroutine parse_ncat

    subroutine parse_predict()
        integer :: i
        logical :: flag, fileExist
        config%mcmc = .true.
        do i = 1, nopt
            if (str_match(trim(register%key(i)), '-predict')) then
                if (cast_logical(register%default(i))) config%mcmc = .false.
            end if
        end do
        if (.not. config%mcmc) then
            do i = 1, nopt
                if (str_match(trim(register%key(i)), '-model')) then
                    config%modfil = trim(register%default(i))
                else if (str_match(trim(register%key(i)), '-freq')) then
                    config%freqfil = trim(register%default(i))
                else if (str_match(trim(register%key(i)), '-param')) then
                    config%paramfil = trim(register%default(i))
                end if
            end do
        end if
    end subroutine parse_predict

    subroutine parse_ndist()
        integer :: i
        do i = 1, nopt
            if (str_match(trim(register%key(i)), '-ndist')) then
                config%ndist = cast_int(register%default(i))
            end if
        end do
    end subroutine parse_ndist

    subroutine parse_priors()
        integer :: i, k, nitem
        character(len=128), dimension(:), allocatable :: c_string
        allocate(c_string(config%ndist))
        do i = 1, nopt
            if (str_match(trim(register%key(i)), '-vara')) then
                mstate%vara = cast_float(register%default(i))
            else if (str_match(trim(register%key(i)), '-vare')) then
                mstate%vare = cast_float(register%default(i))
            else if (str_match(trim(register%key(i)), '-dfvara')) then
                config%dfvara = cast_float(register%default(i))
            else if (str_match(trim(register%key(i)), '-dfvare')) then
                config%dfvare = cast_float(register%default(i))
            else if (str_match(trim(register%key(i)), '-gpin')) then
                nitem = ntokens(register%default(i))
                call tokenize(trim(register%default(i)), ',', config%ndist, c_string)
                do k = 1, config%ndist
                    mstate%gpin(k) = cast_float(trim(c_string(k)))
                end do
            else if (str_match(trim(register%key(i)), '-delta')) then
                nitem = ntokens(register%default(i))
                if (nitem == 1) then
                    mstate%delta = cast_float(trim(register%default(i)))
                else
                    call tokenize(trim(register%default(i)), ',', config%ndist, c_string)
                    do k = 1, config%ndist
                        mstate%delta(k) = cast_float(trim(c_string(k)))
                    end do
                end if
            end if
        end do
    end subroutine parse_priors

    subroutine parse_initialise()
        integer :: i
        do i = 1, nopt
            if (str_match(trim(register%key(i)), '-msize')) then
                config%msize = cast_int(register%default(i))
                config%permute = .true.
            else if (str_match(trim(register%key(i)), '-mrep')) then
                config%mrep = cast_int(register%default(i))
            else if (str_match(trim(register%key(i)), '-numit')) then
                config%numit = cast_int(register%default(i))
            else if (str_match(trim(register%key(i)), '-burnin')) then
                config%burnin = cast_int(register%default(i))
            else if (str_match(trim(register%key(i)), '-thin')) then
                config%thin = cast_int(register%default(i))
            else if (str_match(trim(register%key(i)), '-seed')) then
                config%seed1 = cast_int(register%default(i))
            else if (str_match(trim(register%key(i)), '-snpout')) then
                config%snpout = cast_logical(register%default(i))
            else if (str_match(trim(register%key(i)), '-beta')) then
                config%beta = cast_logical(register%default(i))
            else if (str_match(trim(register%key(i)), '-cat')) then
                config%cat = cast_logical(register%default(i))
            else if (str_match(trim(register%key(i)), '-permute')) then
                config%permute = cast_logical(register%default(i))
            end if
        end do
    end subroutine parse_initialise

    subroutine parse_bayesCpi()
        integer :: i
        config%nobayesCpi = .true.
        do i = 1, nopt
            if (str_match(trim(register%key(i)), '-bayesCpi')) then
                if (cast_logical(register%default(i))) config%nobayesCpi = .false.
            end if
        end do
    end subroutine parse_bayesCpi

end module mod_cmd
