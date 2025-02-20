module OneDim_QM_NI_defs

  use OneDim_QM_NI_kinds, only: dp

  implicit none
  private

  real(dp), parameter, public :: pi = acos(-1.0_dp)

  complex(dp), parameter, public :: cmplx_0 = cmplx(0.0_dp, 0.0_dp, dp)
  complex(dp), parameter, public :: cmplx_1 = cmplx(1.0_dp, 0.0_dp, dp)
  complex(dp), parameter, public :: cmplx_i = cmplx(0.0_dp, 1.0_dp, dp)

  !-\frac{\hbar^2}{2m_e} in eV*\r{A}^2, with m_e the electron mass.
  real(dp), parameter, public :: electronic_kinetic_prefactor = -3.8099821433_dp
  !-\frac{\hbar^2}{2AMU} in peV*(\mu m)^2, with AMU the atomic mass unit.
  real(dp), parameter, public :: atomic_kinetic_prefactor = -20.9007965444_dp

  !(e/\hbar)*10^{-15}. The conversion factor to convert energies
  !\hbar\omega in eV to angular frequencies \omega in fs^{-1},
  !or \hbar\omega in peV to angular frequencies \omega in ms^{-1}.
  real(dp), parameter, public :: time_conversion_factor = 1.519267449_dp

end module OneDim_QM_NI_defs
