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

  public :: set_metric

contains

  function set_metric(N, st, fs) result(u)
    !Spetialized metric to account for boundary conditions.
    integer, intent(in) :: N
    logical, intent(in) :: st, fs
    complex(dp) :: u(N, N)

    integer :: i

    u = cmplx_0
    do i = 1, N
      u(i, i) = cmplx_1
    enddo

    !Do we impose boundary conditions at start and finish?
    !If so, set the corresponding entry to (numerical) zero.
    if (st) u(1, 1) = epsilon(1.0_dp)
    if (fs) u(N, N) = epsilon(1.0_dp)

  end function set_metric

end module OneDim_QM_NI_defs
