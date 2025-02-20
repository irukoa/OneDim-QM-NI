program Basic_Functionality

  use Iso_Fortran_ENV, only: stdout => output_unit, stderr => error_unit

  use OneDim_QM_NI_kinds, only: dp
  use OneDim_QM_NI_defs, only: cmplx_i, pi, &
    electronic_kinetic_prefactor, atomic_kinetic_prefactor, &
    time_conversion_factor
  use OneDim_QM_NI

  implicit none

  type(time_dependent_system) :: td_test

  integer :: i, N, Nt
  character(len=1024) :: aux
  real(dp) :: rnd_tm
  integer :: rnd_tm_instant
  complex(dp) :: tevop_ref
  real(dp) :: tol = 1.0E3_dp*epsilon(1.0_dp)

  real(dp) :: reference_atm_eig(3) = [205.45894233699917_dp, 821.83374557932848_dp, 1849.1183384979399_dp], &
              reference_elec_eig(3) = [37.452874091464579_dp, 149.81112746850590_dp, 337.07365341513753_dp]
  real(dp) :: reference_pos = cmplx(-0.19492613539931294_dp, 0.0_dp, dp)

  call random_seed()

  N = 1000
  Nt = 1000
  write (unit=aux, fmt="(I0)") N

  write (unit=stderr, fmt="(A)") "Starting test program..."
  write (unit=stderr, fmt="(A)") "Consider the particle in a box system with N = "//trim(adjustl(aux))//" states."
  write (unit=stderr, fmt="(A)") "Details: atomic system, 1AMU, [0, 1]\mu m range."
  write (unit=stderr, fmt="(A)") ""

  call td_test%init(name="Atomic_PiaB_Test", &
                    type="atomic", &
                    mass=1.0_dp, &
                    start=0.0_dp, &
                    finish=1.0_dp, &
                    steps=N, &
                    args=[1.0_dp], &
                    potential=particle_in_a_box_pot)

  write (unit=stderr, fmt="(A)") "Eigenvalues obtained. Comparing with reference..."
  do i = 1, 3
    if (abs(td_test%eig(i) - reference_atm_eig(i)) > tol) error stop "FAIL: Eigenvalue mismatch."
  enddo
  write (unit=stderr, fmt="(A)") "PASS"
  write (unit=stderr, fmt="(A)") ""
  write (unit=stderr, fmt="(A)") "INFO: Comparison between analytic and numeric first excitation energy:"
  write (unit=aux, fmt="(F15.8)") - atomic_kinetic_prefactor*(pi**2)*(2**2 - 1**2)
  write (unit=stderr, fmt="(A)") "Analytic: "//trim(adjustl(aux))//" neV."
  write (unit=aux, fmt="(F15.8)") td_test%eig(2) - td_test%eig(1)
  write (unit=stderr, fmt="(A)") "Numeric: "//trim(adjustl(aux))//" neV."
  write (unit=stderr, fmt="(A)") ""

  write (unit=stderr, fmt="(A)") "Comparing position operator..."
  if (abs(td_test%pos(2, 3) - reference_pos) > tol) error stop "FAIL: Value mismatch."
  write (unit=stderr, fmt="(A)") "PASS"
  write (unit=stderr, fmt="(A)") ""
  write (unit=stderr, fmt="(A)") "INFO: Comparison between analytic and numeric position operator's (2, 3) element:"
  write (unit=aux, fmt="(F15.8)") - (8.0_dp*1.0_dp/pi**2)*(2*3)/(2**2 - 3**2)**2
  write (unit=stderr, fmt="(A)") "Analytic: "//trim(adjustl(aux))//" \mu m."
  write (unit=aux, fmt="(F15.8)") real(td_test%pos(2, 3), dp)
  write (unit=stderr, fmt="(A)") "Numeric: "//trim(adjustl(aux))//" \mu m."
  write (unit=stderr, fmt="(A)") ""

  write (unit=stderr, fmt="(A)") "Repeating for: electronic system, 1m_e, [0, 1]\r{A} range."
  write (unit=stderr, fmt="(A)") ""

  call td_test%init(name="Electronic_PiaB_Test", &
                    type="electronic", &
                    mass=1.0_dp, &
                    start=0.0_dp, &
                    finish=1.0_dp, &
                    steps=N, &
                    args=[1.0_dp], &
                    potential=particle_in_a_box_pot)

  write (unit=stderr, fmt="(A)") "Eigenvalues obtained. Comparing with reference..."
  do i = 1, 3
    if (abs(td_test%eig(i) - reference_elec_eig(i)) > tol) error stop "FAIL: Eigenvalue mismatch."
  enddo
  write (unit=stderr, fmt="(A)") "PASS"
  write (unit=stderr, fmt="(A)") ""
  write (unit=stderr, fmt="(A)") "INFO: Comparison between analytic and numeric first excitation energy:"
  write (unit=aux, fmt="(F15.8)") - electronic_kinetic_prefactor*(pi**2)*(2**2 - 1**2)
  write (unit=stderr, fmt="(A)") "Analytic: "//trim(adjustl(aux))//" eV."
  write (unit=aux, fmt="(F15.8)") td_test%eig(2) - td_test%eig(1)
  write (unit=stderr, fmt="(A)") "Numeric: "//trim(adjustl(aux))//" eV."
  write (unit=stderr, fmt="(A)") ""

  write (unit=stderr, fmt="(A)") "Comparing position operator..."
  if (abs(td_test%pos(2, 3) - reference_pos) > tol) error stop "FAIL: Value mismatch."
  write (unit=stderr, fmt="(A)") "PASS"
  write (unit=stderr, fmt="(A)") ""
  write (unit=stderr, fmt="(A)") "INFO: Comparison between analytic and numeric position operator's (2, 3) element:"
  write (unit=aux, fmt="(F15.8)") - (8.0_dp*1.0_dp/pi**2)*(2*3)/(2**2 - 3**2)**2
  write (unit=stderr, fmt="(A)") "Analytic: "//trim(adjustl(aux))//" \r{A}."
  write (unit=aux, fmt="(F15.8)") real(td_test%pos(2, 3), dp)
  write (unit=stderr, fmt="(A)") "Numeric: "//trim(adjustl(aux))//" \r{A}."
  write (unit=stderr, fmt="(A)") ""

  write (unit=stderr, fmt="(A)") "Start testing free dynamical evolution."
  write (unit=stderr, fmt="(A)") "Details: electronic system, 1m_e, [0, 1]\r{A} range,"
  write (unit=stderr, fmt="(A)") "selected 2 states, starting with state #3,"
  write (unit=stderr, fmt="(A)") "Time propagation from 0 to 1 eV^{-1}."
  write (unit=stderr, fmt="(A)") ""

  call td_test%init_td(t_start=0.0_dp, &
                       t_end=1.0_dp, &
                       t_steps=Nt, &
                       t_args=[1.0_dp], &
                       Ht=free_particle_hamil, &
                       selected_states=2, &
                       selected_states_start=3)

  write (unit=stderr, fmt="(A)") "Sampling..."
  write (unit=stderr, fmt="(A)") ""

  call td_test%get_tevop(parallel=.false.)

  write (unit=stderr, fmt="(A)") "Time evolution operator obtained in serial. Comparing with random reference..."
  call random_number(rnd_tm)
  rnd_tm = 1.0_dp + real(Nt - 1, dp)*rnd_tm
  rnd_tm_instant = floor(rnd_tm)
  tevop_ref = exp(-cmplx_i*td_test%eig(3)*time_conversion_factor*real(rnd_tm_instant - 1, dp)/real(Nt - 1, dp))
  if (abs(td_test%tevop(1, 1, rnd_tm_instant) - tevop_ref) > tol) error stop "FAIL: Value mismatch."
  write (unit=stderr, fmt="(A)") "PASS"
  write (unit=stderr, fmt="(A)") ""

  call td_test%get_tevop(parallel=.true.)

  write (unit=stderr, fmt="(A)") "Time evolution operator obtained in parallel. Comparing with random reference..."
  call random_number(rnd_tm)
  rnd_tm = 1.0_dp + real(Nt - 1, dp)*rnd_tm
  rnd_tm_instant = floor(rnd_tm)
  tevop_ref = exp(-cmplx_i*td_test%eig(3)*time_conversion_factor*real(rnd_tm_instant - 1, dp)/real(Nt - 1, dp))
  if (abs(td_test%tevop(1, 1, rnd_tm_instant) - tevop_ref) > tol) error stop "FAIL: Value mismatch."
  write (unit=stderr, fmt="(A)") "PASS"
  write (unit=stderr, fmt="(A)") ""

  write (unit=stderr, fmt="(A)") "Testing finished."

end program Basic_Functionality
