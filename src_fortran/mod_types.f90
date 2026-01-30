module mod_types
    use iso_fortran_env, only: wp => real64
    implicit none

    type BayesConfig
        character(len=256) :: inprefix = 'input'
        character(len=256) :: outprefix = 'output'
        character(len=256) :: catRC = ''
        
        integer :: numit = 1000
        integer :: burnin = 100
        integer :: thin = 1
        integer :: seed = 0
        integer :: trait_pos = 1
        integer :: ndist = 4
        integer :: ncat = 1
        
        real(wp) :: vara = 0.01_wp
        real(wp) :: vare = 0.01_wp
        real(wp) :: dfvara = -2.0_wp
        real(wp) :: dfvare = -2.0_wp
        
        real(wp), allocatable :: gpin(:)
        real(wp), allocatable :: delta(:)
        
        logical :: mcmc = .true.
        logical :: snpout = .false.
        logical :: permute = .false.
        logical :: cat = .false.
        logical :: beta = .false.
        logical :: additive = .false.
        logical :: bayesCpi = .false.
    end type BayesConfig

    type BayesData
        integer :: nind = 0
        integer :: nloci = 0
        integer :: nt = 0
        
        real(wp), allocatable :: why(:)    ! [nind]
        integer, allocatable :: trains(:)  ! [nind]
        real(wp), allocatable :: X(:,:)    ! [nt, nloci]
        integer, allocatable :: C(:,:)     ! [nloci, ncat]
        real(wp), allocatable :: freq(:)   ! [nloci]
    end type BayesData

    type BayesModel
        real(wp) :: mu = 1.0_wp
        real(wp) :: vara = 0.0_wp
        real(wp) :: vare = 0.0_wp
        real(wp) :: vary = 0.0_wp
        
        real(wp), allocatable :: g(:)      ! [nloci]
        real(wp), allocatable :: p(:,:)    ! [ndist, ncat]
        real(wp), allocatable :: gp(:)     ! [ndist]
        real(wp), allocatable :: yadj(:)   ! [nt]
        real(wp), allocatable :: xpx(:)    ! [nloci]
        real(wp), allocatable :: gannot(:,:) ! [nloci, ncat]
        
        integer, allocatable :: snptracker(:,:) ! [nloci, ncat]
        integer, allocatable :: a(:)             ! [nloci] (selected category)
    end type BayesModel

end module mod_types
