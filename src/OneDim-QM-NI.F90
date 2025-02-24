module OneDim_QM_NI

  use Iso_Fortran_ENV, only: stdout => output_unit
  use OMP_LIB

  use OneDim_QM_NI_kinds, only: dp
  use OneDim_QM_NI_defs, only: pi, cmplx_0, cmplx_1, cmplx_i, &
    electronic_kinetic_prefactor, atomic_kinetic_prefactor, &
    time_conversion_factor
  use OneDim_QM_NI_utility, only: diagonalize, inverse, deg_list, schur, SVD, expsh, logu
  use OneDim_QM_NI_implementations, only: particle_in_a_box_pot, harmonic_oscillator_pot, &
    free_particle_hamil, dipole_cosine_hamil

  implicit none
  private

  real(dp), parameter :: underflow_tolerance = 1.0E3*epsilon(1.0_dp)

  !Definitions in this module.
  public :: steady_state_system
  public :: time_dependent_system

  !Implementations of potentials/Hamiltonians.
  public :: particle_in_a_box_pot
  public :: harmonic_oscillator_pot
  public :: free_particle_hamil
  public :: dipole_cosine_hamil

  !Utilities and external constants.
  public :: diagonalize
  public :: inverse
  public :: deg_list
  public :: schur
  public :: SVD
  public :: expsh
  public :: logu

  public :: cmplx_0, cmplx_1, cmplx_i
  public :: pi
  public :: electronic_kinetic_prefactor
  public :: atomic_kinetic_prefactor
  public :: time_conversion_factor

  public :: dp

  abstract interface
    function frml(position, args) result(u)
      !Represents a potential function in the position basis.
      import dp
      real(dp), intent(in) :: position
      real(dp), intent(in) :: args(:)

      real(dp) :: u
    end function frml
  end interface

  abstract interface
    function tdep_hamil(H, pos, time, targs) result(u)
      !Represents a time-dependent Hamiltonian in the Hamiltonian basis.
      import dp
      complex(dp), intent(in) :: H(:, :), pos(:, :)
      real(dp), intent(in) :: time
      real(dp), intent(in) :: targs(:)

      complex(dp) :: u(size(H(:, 1)), size(H(1, :)))
    end function tdep_hamil
  end interface

  type :: steady_state_system
    character(len=120) :: name
    integer :: tp
    real(dp) :: mass
    real(dp) :: start, finish                             !1-Dimensional interval.
    integer :: number_of_states                           !Number of discretization points in the interval (coincides with the number of states).
    complex(dp), allocatable :: Hx(:, :)                  !Hamiltonian in the position basis.
    complex(dp), allocatable :: Hh(:, :)                  !Hamiltonian in the Hamiltonian basis (diagonal).
    real(dp), allocatable :: eig(:)                       !Energy levels.
    complex(dp), allocatable :: rot(:, :)                 !Rotation from position to Hamiltonian basis. Normalized to correspond to wave-functions.
    complex(dp), allocatable :: pos(:, :)                 !Position matrix elements in the Hamiltonian basis.
    procedure(frml), nopass, pointer :: formula => null() !Implementation of the potential function.
    real(dp), allocatable :: args(:)                      !External arguments in potential function.
    logical :: initialized = .false.
  contains
    procedure, pass(self) :: init => construcor_st
  end type

  type, extends(steady_state_system) :: time_dependent_system
    real(dp) :: time_start, time_end                                !Time interval.
    integer :: time_steps                                           !Number of discretization points.
    integer :: selected_states                                      !Selected states (for approximations).
    integer :: selected_states_start = 1                            !First chosen state.
    procedure(tdep_hamil), nopass, pointer :: hamiltonian => null() !Implementation of the time dependent Hamiltonian.
    real(dp), allocatable :: t_args(:)                              !External arguments in time dependent Hamiltonian.
    complex(dp), allocatable :: tevop(:, :, :)                      !Time evolution operator.
    logical :: sampling_ready = .false.
  contains
    procedure, pass(self) :: init_td => construcor_td
    procedure, pass(self) :: get_tevop => td_sampler
  end type

contains

  subroutine construcor_st(self, name, type, mass, start, finish, steps, args, potential, silent)
    class(steady_state_system), intent(out) :: self

    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: type
    real(dp), intent(in) :: mass
    real(dp), intent(in) :: start, finish
    integer, intent(in) :: steps
    real(dp), intent(in) :: args(:)
    procedure(frml) :: potential
    logical, optional, intent(in) :: silent

    real(dp) :: prefactor

    character(len=1024) :: errormsg
    integer :: istat
    logical :: not_silent = .true.

    self%name = name

    select case (type)
    case ("electronic")
      !Expected: potential function in eVs, mass in units
      !of the electron mass, lengts in Angstroms.
      prefactor = electronic_kinetic_prefactor/mass
      self%tp = 0
    case ("atomic")
      !Expected: potential function in pico eVs, mass in units
      !of AMUs, lengths in micro meters.
      prefactor = atomic_kinetic_prefactor/mass
      self%tp = 1
    case default
      error stop "SS System specification error: type is neither of <<electronic>>, <<atomic>>."
    end select

    if (present(silent)) not_silent = .not. silent

    if (not_silent) write (unit=stdout, fmt="(A)") "Initializing SS System <<"//trim(adjustl(self%name))//">>."

    if (start >= finish) error stop "SS System specification error: start >= finish."
    self%start = start
    self%finish = finish
    if (steps < 3) error stop "SS System specification error: steps < 3."
    self%number_of_states = steps

    if (not_silent) write (unit=stdout, fmt="(A)") "start, finish, steps verified..."

    allocate (self%Hx(self%number_of_states, self%number_of_states), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "SS System allocation failure: Hx. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif
    allocate (self%Hh(self%number_of_states, self%number_of_states), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "SS System allocation failure: Hh. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif
    allocate (self%eig(self%number_of_states), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "SS System allocation failure: eig. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif
    allocate (self%rot(self%number_of_states, self%number_of_states), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "SS System allocation failure: rot. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif
    allocate (self%pos(self%number_of_states, self%number_of_states), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "SS System allocation failure: pos. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif

    if (not_silent) write (unit=stdout, fmt="(A)") "Data allocation complete..."

    allocate (self%args(size(args)), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "SS System allocation failure: args. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif
    self%args = args

    if (not_silent) write (unit=stdout, fmt="(A)") "Argument allocation complete..."

    self%formula => potential

    if (not_silent) write (unit=stdout, fmt="(A)") "Potential function linked..."
    self%Hx = Hamil(start=self%start, &
                    finish=self%finish, &
                    steps=self%number_of_states, &
                    formula=self%formula, &
                    prefactor=prefactor, &
                    args=self%args)

    if (not_silent) write (unit=stdout, fmt="(A)") "Hamiltonian matrix calculated (position representation)..."

    self%pos = Pos(start=self%start, &
                   finish=self%finish, &
                   steps=self%number_of_states)

    if (not_silent) write (unit=stdout, fmt="(A)") "Position matrix calculated (position representation)..."
    if (not_silent) write (unit=stdout, fmt="(A)") "Diagonalizing..."

    call diagonalize(matrix=self%Hx, P=self%rot, eig=self%eig, D=self%Hh)

    if (not_silent) write (unit=stdout, fmt="(A)") "Transforming position matrix to Hamiltonian representation..."

    self%pos = matmul(self%rot, matmul(self%pos, conjg(transpose(self%rot))))

    if (not_silent) write (unit=stdout, fmt="(A)") "Normalizing eigenvectors to interpret them as wave-functions..."
    !We normalize rot such that (finish - start)*sum(rot(i, :)*rot(j, :))/(N-1) = \delta_{ij}.
    !This way, rot represents the wavefunction, since the integral of the wavefunction
    !in the [start, finish] range is defined as L*sum(rot(i, :)*rot(j, :))/(N-1).
    self%rot = self%rot*sqrt(real(self%number_of_states - 1, dp))/sqrt(self%finish - self%start)

    if (not_silent) write (unit=stdout, fmt="(A)") "Done."
    if (not_silent) write (unit=stdout, fmt="(A)") ""

    self%initialized = .true.

  end subroutine construcor_st

  function Hamil(start, finish, steps, formula, prefactor, args) result(u)
    !Function to obtain the Hamiltonian in the position representation, Hx(x).
    real(dp), intent(in) :: start, finish
    integer, intent(in) :: steps
    procedure(frml) :: formula
    real(dp), intent(in) :: prefactor
    real(dp), intent(in) :: args(:)

    complex(dp) :: u(steps, steps)

    integer :: i
    real(dp) :: position
    real(dp) :: spacing, &
                k_fac

    spacing = (finish - start)/real(steps - 1, dp)
    k_fac = prefactor/(spacing**2)

    if (abs(k_fac) < underflow_tolerance) error stop "SS System specification error: underflow, too large N"

    u = cmplx_0

    u(1, 1) = -2*k_fac + formula(start, args)
    u(1, 2) = k_fac
    do i = 2, steps - 1
      position = start + spacing*real(i - 1, dp)
      u(i, i + 1) = k_fac
      u(i, i) = -2*k_fac + formula(position, args)
      u(i, i - 1) = k_fac
    enddo
    u(steps, steps - 1) = k_fac
    u(steps, steps) = -2*k_fac + formula(finish, args)
  end function Hamil

  function Pos(start, finish, steps) result(u)
    !Function to obtain the position in the defining representation, x(x).
    real(dp), intent(in) :: start, finish
    integer, intent(in) :: steps

    complex(dp) :: u(steps, steps)

    integer :: i
    real(dp) :: spacing

    u = cmplx_0
    spacing = (finish - start)/real(steps - 1, dp)

    do i = 1, steps
      u(i, i) = cmplx_1*(start + spacing*real(i - 1, dp))
    enddo

  end function Pos

  subroutine construcor_td(self, t_start, t_end, t_steps, t_args, Ht, selected_states, selected_states_start, silent)
    class(time_dependent_system), intent(inout) :: self

    real(dp), intent(in) :: t_start, t_end
    integer, intent(in) :: t_steps
    real(dp), intent(in) :: t_args(:)
    procedure(tdep_hamil) :: Ht
    integer, optional, intent(in) :: selected_states, selected_states_start
    logical, optional, intent(in) :: silent

    character(len=1024) :: errormsg, aux1, aux2
    integer :: istat
    logical :: not_silent = .true.

    if (.not. self%initialized) error stop "TD System specification error: SS System not initialized."
    if (self%sampling_ready) then
      write (unit=stdout, fmt="(A)") "Warning: TD System <<"//trim(adjustl(self%name))//">> was ready for sampling."
      write (unit=stdout, fmt="(A)") "Deallocating TD related allocatables for re-execution..."
      deallocate (self%tevop, self%t_args)
      write (unit=stdout, fmt="(A)") "Done."
    endif

    if (present(silent)) not_silent = .not. silent

    if (not_silent) write (unit=stdout, fmt="(A)") "Initializing TD System <<"//trim(adjustl(self%name))//">>."

    if (t_start >= t_end) error stop "TD System specification error: t_start >= t_end."
    self%time_start = t_start
    self%time_end = t_end

    if (t_steps < 1) error stop "TD System specification error: t_steps < 1."
    self%time_steps = t_steps

    if (not_silent) write (unit=stdout, fmt="(A)") "t_start, t_end, t_steps verified..."

    if (present(selected_states)) then
      if (selected_states < 1) &
        error stop "TD System specification error: selected_states < 1."
      if (selected_states > self%number_of_states) &
        error stop "TD System specification error: selected_states < max_number_of_states."
      self%selected_states = selected_states
      if (present(selected_states_start)) then
        if (selected_states_start < 1) &
          error stop "TD System specification error: selected_states_start < 1."
        if (selected_states_start > self%number_of_states) &
          error stop "TD System specification error: selected_states_start < max_number_of_states."
        self%selected_states_start = selected_states_start
      endif
    endif

    if ((present(selected_states_start)) .and. (.not. present(selected_states))) &
      error stop "TD System specification error: selected_states_start present but selected_states is not."

    write (unit=aux1, fmt="(I0)") self%selected_states
    write (unit=aux2, fmt="(I0)") self%selected_states_start

    if (not_silent) write (unit=stdout, fmt="(A)") &
      "Selecting "//trim(adjustl(aux1))//" states, starting at state #"//trim(adjustl(aux2))//"..."

    allocate (self%tevop(self%selected_states, self%selected_states, self%time_steps), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "SS System allocation failure: tevop. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif

    if (not_silent) write (unit=stdout, fmt="(A)") "Data allocation complete..."

    allocate (self%t_args(size(t_args)), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "TD System allocation failure: t_args. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif
    self%t_args = t_args

    if (not_silent) write (unit=stdout, fmt="(A)") "Argument allocation complete..."

    self%hamiltonian => Ht

    if (not_silent) write (unit=stdout, fmt="(A)") "Time dependent Hamiltonian function linked..."

    self%sampling_ready = .true.

    if (not_silent) write (unit=stdout, fmt="(A)") "TD System ready for sampling."
    if (not_silent) write (unit=stdout, fmt="(A)") ""

  end subroutine construcor_td

  subroutine td_sampler(self, parallel)
    class(time_dependent_system), intent(inout) :: self
    logical, optional, intent(in) :: parallel

    integer :: i
    real(dp) :: dt
    logical :: run_parallel = .true.

    integer :: thr
    character(len=1024) :: aux1, aux2

    if (.not. self%initialized) error stop "TD System sampling error: SS System not initialized."
    if (.not. self%sampling_ready) error stop "TD System sampling error: TD System not ready for sampling."

    write (unit=stdout, fmt="(A)") "Starting the sampling of the time evolution operator..."

    if (present(parallel)) run_parallel = parallel

    dt = (self%time_end - self%time_start)/real(self%time_steps - 1, dp)

    self%tevop(:, :, 1) = cmplx_0

    do i = 1, self%selected_states
      self%tevop(i, i, 1) = cmplx_1
    enddo

    write (unit=aux1, fmt="(I0)") self%time_steps

    if (run_parallel) then
      write (unit=stdout, fmt="(A)") "Parallel calculation of the time dependent hamiltonian:"
!$OMP     PARALLEL
      thr = OMP_GET_NUM_THREADS()
      write (unit=aux2, fmt="(I0)") thr
!$OMP     SINGLE
      write (unit=stdout, fmt="(A)") &
        "Distributing "//trim(adjustl(aux1))//" processes across "//trim(adjustl(aux2))//" OMP threads."
!$OMP     END SINGLE
!$OMP     DO
      do i = 2, self%time_steps
        self%tevop(:, :, i) = self%hamiltonian(H=self%Hh(self%selected_states_start: &
                                                         self%selected_states_start + self%selected_states - 1, &
                                                         self%selected_states_start: &
                                                         self%selected_states_start + self%selected_states - 1), &
                                               pos=self%pos(self%selected_states_start: &
                                                            self%selected_states_start + self%selected_states - 1, &
                                                            self%selected_states_start: &
                                                            self%selected_states_start + self%selected_states - 1), &
                                               time=self%time_start + dt*real(i - 1, dp), &
                                               targs=self%t_args)
        self%tevop(:, :, i) = expsh(-cmplx_i*dt*time_conversion_factor*self%tevop(:, :, i))
      enddo
!$OMP     END DO
!$OMP     END PARALLEL
    else
      write (unit=stdout, fmt="(A)") "Sequential calculation of the time dependent hamiltonian."
      do i = 2, self%time_steps
        self%tevop(:, :, i) = self%hamiltonian(H=self%Hh(self%selected_states_start: &
                                                         self%selected_states_start + self%selected_states - 1, &
                                                         self%selected_states_start: &
                                                         self%selected_states_start + self%selected_states - 1), &
                                               pos=self%pos(self%selected_states_start: &
                                                            self%selected_states_start + self%selected_states - 1, &
                                                            self%selected_states_start: &
                                                            self%selected_states_start + self%selected_states - 1), &
                                               time=self%time_start + dt*real(i - 1, dp), &
                                               targs=self%t_args)
        self%tevop(:, :, i) = expsh(-cmplx_i*dt*time_conversion_factor*self%tevop(:, :, i))
      enddo
    endif

    write (unit=stdout, fmt="(A)") "Sampling done, accumulating..."
    do i = 2, self%time_steps
      self%tevop(:, :, i) = matmul(self%tevop(:, :, i - 1), self%tevop(:, :, i))
    enddo
    write (unit=stdout, fmt="(A)") "Done."
    write (unit=stdout, fmt="(A)") ""

  end subroutine td_sampler

end module OneDim_QM_NI
