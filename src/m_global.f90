module m_global
  use iso_fortran_env, only: r8 => real64
  implicit none
  real(kind=r8), parameter :: pi = 3.141592653589793_r8
  real(kind=r8), parameter :: degenerate = 2._r8 ! degree of degeneration
  real(kind=r8), parameter :: atm = 1.E5_r8 ! standard atmosphere: 100 kPa
  real(kind=r8), parameter :: mass_e = 9.1093837015E-31_r8 ! electron rest mass 
  real(kind=r8), parameter :: planck_const = 6.62607015E-34_r8 ! Planck constant
  real(kind=r8), parameter :: Boltzmann_const = 1.380649E-23_r8 ! Boltzmann constant
  real(kind=r8), parameter :: gas_const = 8.31446261815324_r8 ! gas constant

end module