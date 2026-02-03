program test_standardize
    use mod_defs, only: dp, GENOTYPE_MISSING_THRESHOLD
    use mod_data
    use mod_io
    use mod_standardize
    implicit none

    type(ModelConfig) :: config_local
    type(GenomicData) :: gdata_local
    integer :: i, j

    config_local%inprefix = "test_data"
    config_local%phenfil = "test_data.fam"
    config_local%bimfil = "test_data.bim"
    config_local%genfil = "test_data.bed"
    config_local%freqfil = "test_data.frq"
    config_local%trait_pos = 1
    config_local%mcmc = .true.

    call get_size(config_local, gdata_local)
    call load_phenos_plink(config_local, gdata_local)
    
    ! Allocate X
    gdata_local%nt = count(gdata_local%trains == 0)
    allocate(gdata_local%X(gdata_local%nt, gdata_local%nloci))
    allocate(gdata_local%freqstore(gdata_local%nloci))
    
    call load_snp_binary(config_local, gdata_local)
    
    ! Run xcenter
    call xcenter(config_local, gdata_local)

    print *, "Allele Frequencies (freqstore):"
    do j = 1, gdata_local%nloci
        print "(E25.16)", gdata_local%freqstore(j)
    end do

    print *, "Standardized Genotypes (X):"
    do i = 1, gdata_local%nt
        do j = 1, gdata_local%nloci
            print "(E25.16)", gdata_local%X(i, j)
        end do
    end do

end program test_standardize
