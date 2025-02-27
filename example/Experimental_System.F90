program Experimental_System

  use Iso_Fortran_ENV, only: stdout => output_unit

  use OneDim_QM_NI, only: dp, cmplx_1, pi, &
    time_dependent_system, SVD

  implicit none

  !This is the implementation of the experimental set-up.
  !It is a simplified version of the program I created when
  !I was in Firenze, but the fundamentals are here.

  !First, we define the constant to pass from nK to peV
  !and give the potassium mass in AMUs.
  real(dp), parameter :: kb_over_e_in_peV_per_nK = 0.08617333262_dp
  real(dp), parameter :: potassium_mass_in_AMU = 39.0983_dp

  type(time_dependent_system) :: Exp_Sys

  !In this case, the control panel contains the experimentally tuneable parameters,
  !and the number of states. This is the input experimental people gave me.

  !---CONTROL PANEL---!
  real(dp), parameter :: box_length = 10.0_dp, &
                         param_nn = 20.0_dp

  real(dp) :: leading_amplitude = 350*kb_over_e_in_peV_per_nK, &
              amp2_over_leading = 300.0_dp/350.0_dp, &
              amp3_over_leading = 270.0_dp/350.0_dp, &
              leading_wavelength = 4*box_length/param_nn, &
              wav2_over_leading = param_nn/(param_nn + 1.0_dp), &
              wav3_over_leading = param_nn/(param_nn - 1.0_dp), &
              spring_constant = 0.0_dp, &
              force_constant = 0.0_dp

  integer, parameter :: N = 1025
  !-------------------!

  !Definitions to obtain the tunnelling coefficient.
  real(dp) :: variance = 0.1_dp*box_length
  complex(dp) :: locL(N), locR(N)
  complex(dp) :: overlap(2, 2), U(2, 2), V(2, 2), &
                 Hw(2, 2)
  real(dp) :: eig(2)

  !Auxiliary.
  integer :: i, out
  character(len=120) :: aux
  real(dp) :: x

  !I define an array containing the arguments. This will be passed by the call to initialize
  !the steady state system to the function implementing the potential in the contains section.
  real(dp) :: pot_args(8)

  pot_args = [leading_amplitude, leading_amplitude*amp2_over_leading, leading_amplitude*amp3_over_leading, &
              leading_wavelength, leading_wavelength*wav2_over_leading, leading_wavelength*wav3_over_leading, &
              spring_constant, force_constant]

  !Plot the potential.
  open (newunit=out, action="write", status="unknown", file="Exp_Potential.dat")
  do i = 1, N
    x = -box_length/2 + box_length*real(i - 1, dp)/real(N - 1, dp)
    write (unit=out, fmt="(F15.8, 2xF15.8)") x, &
      experimental_potential(position=x, args=pot_args)
  enddo
  close (unit=out)

  !Obtain the levels, wavefunctions...
  call Exp_Sys%init(name="Experimental_config.", &
                    type="atomic", &
                    mass=potassium_mass_in_AMU, &
                    start=-box_length/2, &
                    finish=box_length/2, &
                    steps=N, &
                    args=pot_args, &
                    potential=experimental_potential)

  !Plot the energy differences.
  open (newunit=out, action="write", status="unknown", file="Exp_Levels.dat")
  do i = 2, N
    write (unit=out, fmt="(I0, 2xI0, 2xF15.8)") i - 1, i, Exp_Sys%eig(i) - Exp_Sys%eig(i - 1)
  enddo
  close (unit=out)

  !---------------------SOME THOUGHTS.
  !Now comes my doubt No. 1: the potential implemented here is not similar to that appearing in the
  !poster, but it is what they have in the MATLAB programs.

  !And the doubt No. 2: the spectrum is continuum-like after the first few levels
  !(the differences between states are almost constant),
  !this is nice, but I do not know how many states to include in a future dynamics’ simulation.
  !My guess may be just the first 2~3, since there is a nice 0.8850 peV difference between states #1 and #2
  !and a really small 0.00323361 peV difference between states #2 and #3, making them almost degenerate.
  !Then, we have a big jump of 2.530 peV-s from level #3 to #4. If we take only the first 2 levels,
  !we can simulate everything as a 2-level system, and we can define \sigma states and everything we need.

  !Remark: the energy h*f = 2.53050062 peV corresponds to a frequency of 612 Hz-s, which I believe is
  !an experimentally achievable driving frequency.

  !Doubt No. 3: regarding dynamics, I do not remember the formula of the time dependent Hamiltonian, which we need,
  !but as long as it was expressed in terms of the Hamiltonian, the position and external arguments everything
  !should be doable easily.
  !--------------------- END THOUGHTS.

  !---------------------TUNNELLING CALCULATION.

  !Now, I will obtain the tunnelling coefficient after truncating the basis to just the first 2 states.
  !Note that the tunnelling coefficients are intrinsically gauge dependent quantities,
  !but let me present what I believe is the correct procedure to obtain them:

  !First, I define left and right localized states, |L'> and |R'>, using Gaussian functions, each localized in a side of the well.
  do i = 1, N
    x = -box_length/2 + box_length*real(i - 1, dp)/real(N - 1, dp)
    locL(i) = cmplx_1*exp(-0.5_dp*(x + 0.25_dp*box_length)**2/variance)
    locR(i) = cmplx_1*exp(-0.5_dp*(x - 0.25_dp*box_length)**2/variance)
  enddo
  !I normalize them,
  locL = (locL/sqrt(real(sum(locL*locL), dp)))*sqrt(real(N - 1, dp))/sqrt(box_length)
  locR = (locR/sqrt(real(sum(locR*locR), dp)))*sqrt(real(N - 1, dp))/sqrt(box_length)

  !and plot them together with the wavefunctions corresponding to the first two levels, |1>, |2>.
  open (newunit=out, action="write", status="unknown", file="Exp_WFCs_and_loc.dat")
  do i = 1, N
    x = -box_length/2 + box_length*real(i - 1, dp)/real(N - 1, dp)
    write (unit=out, fmt="(F15.8, 4(2xF15.8))") x, real(Exp_Sys%rot(i, 1), dp), real(Exp_Sys%rot(i, 2), dp), &
      real(locL(i), dp), real(locR(i), dp)
  enddo

  !Ideally, the basis used to express the Hamiltonian should be like this, so a state would be composed of a superposition of
  !|L'> and |R'> states. The tunnelling coefficient would then be just <L'|H|R'>.
  !In reality, this is not possible since these states do not span the same subspace as the
  !|1>, |2> states do, |L'><L'| + |R'><R'| =/ |1><1| + |2><2| =/ Id.

  !So, I define a matrix containing the products:
  ! <1|L'> <2|L'>
  ! <1|R'> <2|R'>
  overlap(1, 1) = box_length*sum(locL*conjg(Exp_Sys%rot(:, 1)))/real(N - 1, dp)
  overlap(1, 2) = box_length*sum(locL*conjg(Exp_Sys%rot(:, 2)))/real(N - 1, dp)
  overlap(2, 1) = box_length*sum(locR*conjg(Exp_Sys%rot(:, 1)))/real(N - 1, dp)
  overlap(2, 2) = box_length*sum(locR*conjg(Exp_Sys%rot(:, 2)))/real(N - 1, dp)
  !This matrix measures to what degree the basis {|L'>, |R'>} span the same space as
  !the {|1>, |2>} do. In other words, how "aligned" these are in the infinite dimensional
  !state space. If they were aligned, this matrix would be unitary, else no.

  !At this point I perform a singular value decomposition (SVD) of this matrix. This process
  !decomposes a matrix as a product USV*, with U, V unitary and S diagonal with real positive entries.
  call SVD(matrix=overlap, U=U, V=V, eig=eig)
  !The diagonal entries of S, the eigenvalues, measure the alignment of the two bases.
  !Ideally, both eigenvalues should be 1, I will show them in the terminal:
  write (unit=stdout, fmt="(A)"), "SVD Decomposition:"
  write (unit=aux, fmt="(F15.8)") eig(1)
  write (unit=stdout, fmt="(A)"), "Eigenvalue #1 = "//trim(adjustl(aux))//"."
  write (unit=aux, fmt="(F15.8)") eig(2)
  write (unit=stdout, fmt="(A)"), "Eigenvalue #2 = "//trim(adjustl(aux))//"."
  write (unit=stdout, fmt="(A)"), ""
  !As you see, not even close! The takeaway is that the |L'> and |R'> states cannot
  !span the {|1>, |2>} basis. However, we can define an alternative basis, {|L>, |R>},
  !which spans the same space as the {|1>, |2>} basis and is as similar as possible
  !to the previous {|L'>, |R'>} basis. This is done by considering the following:
  overlap = matmul(U, conjg(transpose(V)))
  !Notice that the eigenvalues of the matrix S do not appear here, but U, and V do.
  !This assures that the newly defined overlap matrix is unitary.

  !Now, we extract the newly defined left and right localized states, |L> and |R>,
  !in terms of the original |1>, |2> eigenstates of the Hamiltonian.
  locL = overlap(1, 1)*Exp_Sys%rot(:, 1) + overlap(1, 2)*Exp_Sys%rot(:, 2)
  locR = overlap(2, 1)*Exp_Sys%rot(:, 1) + overlap(2, 2)*Exp_Sys%rot(:, 2)
  !As I mentioned, these are as similar as possible to the localized states |L'> and
  !|R'> while spanning the same space as the {|1>, |2>} basis.
  !I will plot the probability density:
  open (newunit=out, action="write", status="unknown", file="Exp_Loc_PDs.dat")
  do i = 1, N
    x = -box_length/2 + box_length*real(i - 1, dp)/real(N - 1, dp)
    write (unit=out, fmt="(F15.8, 2(2xF15.8))") x, abs(locL(i))**2, abs(locR(i))**2
  enddo
  !As you see, these are not as localized as the Gaussians, but this is what
  !the most localized basis made up only of the two states |1>, |2> looks like.
  !In a sense, this basis is the "Wannier gauge" of this simple system.
  !All that is left is to rotate the diagonal Hamiltonian to this localized basis.
  Hw = matmul(overlap, matmul(Exp_Sys%Hh(1:2, 1:2), conjg(transpose(overlap))))
  !The off-diagonal elements are the tunnelling coefficients.
  write (unit=stdout, fmt="(A)"), "Tunnelling coefficient:"
  write (unit=aux, fmt="(F15.8)") real(Hw(1, 2), dp)
  write (unit=stdout, fmt="(A)"), "[Re] = "//trim(adjustl(aux))//" peV."
  write (unit=aux, fmt="(F15.8)") aimag(Hw(1, 2))
  write (unit=stdout, fmt="(A)"), "[Im] = "//trim(adjustl(aux))//" peV."
  write (unit=stdout, fmt="(A)"), ""

  !---------------------END TUNNELLING CALCULATION.

contains

  function experimental_potential(position, args) result(u)
    !The "args" appearing here refers to the experimental arguments.
    !I gave the values above. Notice that the order of arguments is
    !important, as well as the number of arguments.
    real(dp), intent(in) :: position
    real(dp), intent(in) :: args(:)

    real(dp) :: u

    u = args(1)*(sin(2*pi*position/args(4)))**2 + &
        args(2)*(sin(2*pi*position/args(5)))**2 + &
        args(3)*(sin(2*pi*position/args(6)))**2 + &
        0.5_dp*args(7)*(position)**2 + &
        args(8)*position
  end function experimental_potential

end program Experimental_System
