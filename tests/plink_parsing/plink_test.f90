program plink_test
    use mod_defs, only: dp, MISSING_VALUE
    use mod_data
    use mod_io
    implicit none

    integer :: i, j

    config%inprefix = "test_data"
    config%phenfil = "test_data.fam"
    config%bimfil = "test_data.bim"
    config%genfil = "test_data.bed"
    config%trait_pos = 1
    config%mcmc = .true.

    call get_size(config, gdata)
    call load_phenos_plink(config, gdata)
    
    ! Allocate X
    gdata%nt = count(gdata%trains == 0)
    allocate(gdata%X(gdata%nt, gdata%nloci))
    
    call load_snp_binary(config, gdata)

    print *, "Phenotypes (why):"
    do i = 1, gdata%nind
        if (gdata%why(i) == MISSING_VALUE) then
            print *, "NA"
        else
            print "(F20.6)", gdata%why(i)
        end if
    end do

    print *, "Genotypes (X):"
    do i = 1, gdata%nt
        do j = 1, gdata%nloci
            print "(F5.1)", gdata%X(i, j)
        end do
    end do

end program plink_test
