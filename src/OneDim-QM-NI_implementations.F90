module OneDim_QM_NI_implementations

  use OneDim_QM_NI_kinds, only: dp
  use OneDim_QM_NI_defs, only: time_conversion_factor

  implicit none
  private

  public :: particle_in_a_box_pot
  public :: harmonic_oscillator_pot
  public :: free_particle_hamil
  public :: dipole_cosine_hamil

contains

  function particle_in_a_box_pot(position, args) result(u)
    real(dp), intent(in) :: position
    real(dp), intent(in) :: args(:)

    real(dp) :: u

    u = 0.0_dp
  end function particle_in_a_box_pot

  function harmonic_oscillator_pot(position, args) result(u)
    real(dp), intent(in) :: position
    real(dp), intent(in) :: args(:)

    real(dp) :: u
    !args(1) stands for the "spring constant".
    u = 0.5_dp*args(1)*position**2
  end function harmonic_oscillator_pot

  function free_particle_hamil(H, pos, time, targs) result(u)
    complex(dp), intent(in) :: H(:, :), pos(:, :)
    real(dp), intent(in) :: time
    real(dp), intent(in) :: targs(:)

    complex(dp) :: u(size(H(:, 1)), size(H(1, :)))

    u = H
  end function free_particle_hamil

  function dipole_cosine_hamil(H, pos, time, targs) result(u)
    complex(dp), intent(in) :: H(:, :), pos(:, :)
    real(dp), intent(in) :: time
    real(dp), intent(in) :: targs(:)

    complex(dp) :: u(size(H(:, 1)), size(H(1, :)))
    !targs(1) stands for e*amplitude, measured in eV/A or peV/\mu m.
    !targs(2) stands for the driving frequency in either eV-s or peV-s.
    u = H + targs(1)*cos(targs(2)*time_conversion_factor*time)*pos
  end function dipole_cosine_hamil

end module OneDim_QM_NI_implementations
