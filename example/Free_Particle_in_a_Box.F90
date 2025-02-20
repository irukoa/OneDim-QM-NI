program Free_Particle_in_a_Box

  use Iso_Fortran_ENV, only: stdout => output_unit

  use OneDim_QM_NI_kinds, only: dp
  use OneDim_QM_NI_defs, only: cmplx_i, pi, electronic_kinetic_prefactor
  use OneDim_QM_NI, only: time_dependent_system

  implicit none

  !In this example I will show you how I simulate a simple quantum system's
  !steady state properties and dynamics. I know that this file is full
  !of comments and that can overwhelming. However, I believe that reading
  !carefully these comments will help you understand whats going on,
  !gaining insight on the numerical scheme.
  !Also, I recommend that you read the terminal log to understand
  !what the code is doing.

  !The system we are considering is really simple, a free particle in
  !a box of size L, so recall these results:
  !Levels: E_n = \hbar^2 \pi^2 n^2 / (2m L^2),
  !Position matrix elements: r_{nm} = -8 L n m / (pi^2 (n^2 - m^2)^2) if n+m is odd
  !and 0 otherwise.

  !I am confident that this simple example can serve as a tutorial to get used
  !to the code and that soon you will be dealing with much more complex
  !things with ease.

  !Cheers!,
  !Álvaro

  type(time_dependent_system) :: PiaB !Object containing the system.

  !---CONTROL PANEL---!
  !These parameters apply to the definition of the steady state system.
  real(dp) :: box_length = 1.0_dp

  integer :: N = 100 !Discretization points/number of states.

  !These parameters apply to the dynamical system.
  real(dp) :: time_interval = 1.0_dp
  integer :: Nt = 1000                        !Discretization points for the time evolution operator.
  integer :: number_of_selected_levels = 2, & !This says that we will consider 2 levels for dynamics,
             starting_level = 1               !starting at level #1, meaning that we will consider states #1 and #2.
  !-------------------!

  !Auxiliary.
  integer :: i, out
  real(dp) :: x, aux1, aux2

  !This is the way to initialize a steady state system:
  call PiaB%init(name="Electronic_PiaB", &              !we give it a name,
                 type="electronic", &                   !we specify if it is an electronic system or an atomic one(*),
                 mass=1.0_dp, &                         !we give it a mass,
                 start=-box_length/2, &                 !we set the starting and ending points of the one dimensional line,
                 finish=box_length/2, &
                 steps=N, &                             !we set the number of discretization steps, which coincide with the number of states,
                 args=[1.0_dp], &                       !we provide any "external arguments" appearing in the potential (in this case, a dummy array),
                 potential=particle_in_a_box_potential) !and we point the program to our implementation of the potential (see the contains section).

  !(*) If the system is "electronic", the potential must be expressed in eV-s, the mass in units of the electron mass m_e and all lenghts are in \r{A}-s.
  !If the system is "atomic", the potential must be expressed in peV-s, the mass in units of the atomic mass unit, AMU, and all lenghts are in \mu m-s.

  !When running the above line, the program will write to elements PiaB%eig(:), PiaB%rot(:, :), PiaB%pos(:, :), which contain:
  !the energy levels in ascending order, the wavefunctions \psi_n(x_i) = PiaB%rot(i, n), and the position operator in the basis
  !which makes the Hamiltonian diagonal, x_{nm} = PiaB%pos(n, m), respectively.

  !First, since the particle in a box can be solved analitically, we will compare the accuracy of the levels.
  !We do this by writing to file energy differences between contiguous levels.
  open (newunit=out, action="write", status="unknown", file="PiaB_Levels.dat")
  do i = 2, N
    aux1 = PiaB%eig(i) - PiaB%eig(i - 1)                                             !Numerical result,
    aux2 = -(pi**2)*electronic_kinetic_prefactor*((i)**2 - (i - 1)**2)/box_length**2 !analitical result (the kinetic prefactor is -\frac{\hbar^2}{2m_e} in eV*\r{A}^2),
    write (unit=out, fmt="(I0, 2xI0, 3(2xF15.8))") i - 1, i, aux1, aux2, &
      100*abs(aux1 - aux2)/abs(aux2)                                                 !and relative error.
  enddo
  close (unit=out)
  !We can either examine the file or use gnuplot to see how the numerics relate to analitic results.
  !We observe that only the first ~8 levels provide errors below 5%. By repeating the calculation by setting
  !N = 1000 in the control panel, we see that the error decreases significantly, given that the first
  !~169 levels provide errors below 5%. The takeaway is that a large N has to be set in order to
  !obtain accurate results, and even in that case, we can only consider a handful of well described
  !low energy states as accurate. This will be important when we consider dynamics, since we will
  !be truncating the state space.

  !We can also plot wavefunctions:
  open (newunit=out, action="write", status="unknown", file="PiaB_GroundState.dat")
  do i = 1, N
    x = -box_length/2 + box_length*real(i - 1, dp)/real(N - 1, dp)
    write (unit=out, fmt="(F15.8, 2xF15.8)") x, real(PiaB%rot(i, 1), dp)
  enddo
  close (unit=out)
  !this is nodeless, as expected,
  open (newunit=out, action="write", status="unknown", file="PiaB_ExcitedState.dat")
  do i = 1, N
    x = -box_length/2 + box_length*real(i - 1, dp)/real(N - 1, dp)
    write (unit=out, fmt="(F15.8, 2xF15.8)") x, real(PiaB%rot(i, 2), dp)
  enddo
  close (unit=out)
  !and this contains a node.

  !So far so good, we know that the results are accurate and can be improved by considering larger and larger N-s.
  !Just one thing, the computational time increases as N^3 due to a step involving matrix diagonalization,
  !so take that into account when setting N, also, if N is really really large the program will complain
  !due to underflow errors involving the numerical shceme.

  !Now, on to dynamics. Since we will consider that the particle is free,
  !the time dependence of each level is straightforward, just exp(-i*E_n*t/\hbar).
  !We will verify this.

  !This is the way we initialize the dynamical simulation:
  call PiaB%init_td(t_start=0.0_dp, &                            !We proide a time interval in either fs-s or ms-s(*),
                    t_end=time_interval, &
                    t_steps=Nt, &                                !we set the number of discretization points of the time interval,
                    t_args=[1.0_dp], &                           !consider any external arguments (again, in this case, none),
                    Ht=non_interacting_particle_hamiltonian, &   !point the program to the implementation of H(t),
                    selected_states=number_of_selected_levels, & !and lastly, select a subset of the N states(**),
                    selected_states_start=starting_level)        !from starting_level to starting_level + number_of_selected_levels - 1.

  !(*) If the system is "electronic", time intervals are expressed in fs-s and if the system is "atomic", in ms-s,
  !(**) Given that highly energetic states are not accurately described, in practice we will consider just a handful of accurately described states.
  !It is our responsibility to consider N large enough such that states [starting_level, starting_level + number_of_selected_levels - 1] are accurately
  !described.

  !This will perform the dynamical simulation by employing a exponential integrator (the accuratest thing we can get). The larger Nt, the more accurate
  !time propagation. However, take care, computational time scales as Nt*(number_of_selected_levels)^3. Usually, you will need to find
  !a compromise between Nt and number_of_selected_levels, or calculations become large really fast.
  call PiaB%get_tevop(parallel=.true.)
  !I have added a parallelization option, I advise you to use it, it will speed things up by a factor of 10 in most computers,
  !but you can also disable it by setting "parallel = .false.".

  !This will write into U(t_k, t_start) = PiaB%tevop(:, :, k), which is the time evolution operator in the basis that makes the Hamiltonian diagonal,
  !PiaB%tevop(:, :, Nt) propagates a state from t_start to t_end. Consider the first excited state, with energy PiaB%eig(2). We know that after the
  !the considered interval of 1fs, the phase acquired by a state in that state is:
  write (unit=stdout, fmt="(A)") "Analitically accumulated phase (Re, Im):"
  write (unit=stdout, fmt="(F15.8, 2xF15.8)") &
    real(exp(-cmplx_i*PiaB%eig(2)*1.519267449_dp*1.0_dp), dp), aimag(exp(-cmplx_i*PiaB%eig(2)*1.519267449_dp*1.0_dp))
  !with the factor 1.519267449 = e*10^{-15}/\hbar the conversion factor to pass from eV to fs^{-1}. As for the numerical time evolution operator,
  write (unit=stdout, fmt="(A)") "Numerically accumulated phase (Re, Im):"
  write (unit=stdout, fmt="(F15.8, 2xF15.8)") real(PiaB%tevop(2, 2, Nt), dp), aimag(PiaB%tevop(2, 2, Nt))
  write (unit=stdout, fmt="(A)") ""

  !See you later in the next example, where I consider the experimental potential.

contains

  function particle_in_a_box_potential(position, args) result(u)
    !Express the potential in terms of the position x and external arguments.
    real(dp), intent(in) :: position
    real(dp), intent(in) :: args(:)

    real(dp) :: u

    u = 0.0_dp
  end function particle_in_a_box_potential

  function non_interacting_particle_hamiltonian(H, pos, time, targs) result(u)
    !Express the time dependent Hamiltonian in terms of the Hamiltonian, position operator,
    !time variable and external arguments.
    complex(dp), intent(in) :: H(:, :), pos(:, :)
    real(dp), intent(in) :: time
    real(dp), intent(in) :: targs(:)

    complex(dp) :: u(size(H(:, 1)), size(H(1, :)))

    u = H
  end function non_interacting_particle_hamiltonian

end program
