module mod_data
    use mod_defs, only: dp
    implicit none

    type :: ModelConfig
        integer :: numit, burnin, thin, ndist, seed1, trait_pos, ncat, msize, mrep
        integer :: indistflag, burn, annotflag
        integer :: unit_log, unit_hyp, unit_loc, unit_cat, unit_beta
        logical :: mcmc, snpout, permute, cat, beta, mixture, nobayesCpi, VCE
        real(dp) :: dfvara, dfvare, vara_ap, vare_ap
        character(len=200) :: genfil, phenfil, bimfil, inprefix, outprefix
        character(len=200) :: logfil, freqfil, mbvfil, hypfil, locfil
        character(len=200) :: modfil, paramfil, betafil, catfil, catRC
    end type ModelConfig

    type :: GenomicData
        integer :: nloci, nind, nt, ntrain, ntest, nref
        real(dp), dimension(:), allocatable :: why, pred, freqstore, xpx, includedloci
        real(dp), dimension(:,:), allocatable :: X
        integer, dimension(:,:), allocatable :: C
        integer, dimension(:), allocatable :: nannot, a, vsnptrack, trains, permvec, permannot, atemp
        integer, dimension(:,:), allocatable :: snptracker
        real(dp), dimension(:,:), allocatable :: gannot
        real(dp) :: msep, bhat, ahat, corr
    end type GenomicData

    type :: MCMCState
        real(dp) :: mu, vara, vare, scale, yhat, vary, nnind
        integer :: included, counter, rep
        real(dp), dimension(:), allocatable :: gp, gpin, delta, g, yadj
        real(dp), dimension(:), allocatable :: dirx, pia, ytemp, dira, s, stemp, sstemp, ss
        real(dp), dimension(:,:), allocatable :: varindist, p, log_p
        integer, dimension(:,:), allocatable :: snpindist
        real(dp), dimension(:), allocatable :: log_gp, vare_gp
        real(dp) :: zz, rhs, lhs, v1, gk, zz_vare
        real(dp) :: logdetV, uhat, total_ssq, detV, maxs, maxtemp
        real(dp) :: xhat, sk, skk, r, ssculm, clike
    end type MCMCState

    type :: MCMCStorage
        real(dp), dimension(:), allocatable :: gstore, mu_vare_store
        real(dp), dimension(:), allocatable :: varustore, varistore
        real(dp), dimension(:,:), allocatable :: snpstore, varstore
        real(dp), dimension(:,:), allocatable :: indiststore, pstore, annotstore
    end type MCMCStorage

    type(ModelConfig) :: config
    type(GenomicData), target :: gdata
    type(MCMCState), target   :: mstate
    type(MCMCStorage) :: mstore

    integer :: ios, clock
    real(dp), pointer :: z(:)

end module mod_data
