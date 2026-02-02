module mod_standardize
    use mod_defs, only: dp, GENOTYPE_MISSING_THRESHOLD
    use mod_data, only: ModelConfig, GenomicData
    implicit none

contains

    subroutine xcenter(config, gdata)
        type(ModelConfig), intent(in) :: config
        type(GenomicData), intent(inout) :: gdata
        integer :: j, nomiss
        real(dp) :: q, qtest
        real(dp), dimension(gdata%nt) :: xtemp
        if (config%mcmc) then
            do j = 1, gdata%nloci
                xtemp = gdata%X(:, j)
                nomiss = count(xtemp < GENOTYPE_MISSING_THRESHOLD)
                q = sum(xtemp, mask=xtemp < GENOTYPE_MISSING_THRESHOLD) / (2.0d0 * nomiss)
                if (q == 1.0d0 .or. q == 0.0d0) then
                    gdata%X(:, j) = 0.0d0
                else
                    where (xtemp > 2.0d0) xtemp = 2.0d0 * q
                    gdata%X(:, j) = (xtemp - 2.0d0 * q) / sqrt(2.0d0 * q * (1.0d0 - q))
                end if
                gdata%freqstore(j) = q
            end do
            open(45, file=trim(config%freqfil), status='unknown')
            do j = 1, gdata%nloci
                write(45, '(F10.6)') gdata%freqstore(j)
            end do
            close(45, status='keep')
        else
            open(45, file=trim(config%freqfil), status='unknown')
            do j = 1, gdata%nloci
                read(45, '(E15.7)') gdata%freqstore(j)
            end do
            close(45, status='keep')
            do j = 1, gdata%nloci
                q = gdata%freqstore(j)
                if (q == 1.0d0 .or. q == 0.0d0) then
                    gdata%X(:, j) = 0.0d0
                else
                    xtemp = gdata%X(:, j)
                    nomiss = count(xtemp < GENOTYPE_MISSING_THRESHOLD)
                    qtest = real(sum(xtemp, mask=xtemp < GENOTYPE_MISSING_THRESHOLD), dp) / (2.0d0 * nomiss)
                    where (xtemp > 2.0d0) xtemp = 2.0d0 * qtest
                    gdata%X(:, j) = (xtemp - 2.0d0 * q) / sqrt(2.0d0 * q * (1.0d0 - q))
                end if
            end do
        end if
    end subroutine xcenter

end module mod_standardize
