module mod_defs
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    real(dp), parameter :: PI = 3.141592653589793238462d0
    real(dp), parameter :: LOG_UPPER_LIMIT = 700.0d0
    real(dp), parameter :: MISSING_VALUE = -9999.0d0
    real(dp), parameter :: GENOTYPE_MISSING_THRESHOLD = 3.0d0
end module mod_defs
