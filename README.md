[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)
# OneDim-QM-NI
### One Dimensional Quantum Mechanics - Numerical Implementation.

This is a tiny modern Fortran library for illustration and teaching purposes. Its aim is to study Fortran and one dimensional quantum systems.

<details>
  <summary>Details: working principle</summary>

## Working principle

### Steady state systems.

The library solves the time independent Schroedinger equation in the position basis,

$$
\frac{-\hbar^2}{2m}\frac{\partial^2 \Psi(x)}{\partial x^2} + V(x)\Psi(x) = E_n \Psi(x).
$$

The laplace operator is considered by its finite difference approximation,

$$
\frac{\partial^2 \Psi(x)}{\partial x^2} = \frac{\Psi(x+\epsilon) -2\Psi(x) + \Psi(x-\epsilon)}{\epsilon^2}.
$$

This leads to a straightforward matrix representation of the Hamiltonian operator in the position basis. By employing the discretization

$$
x \to x_i = x_s + x_f - x_s \frac{i-1}{N-1},\;\epsilon = \frac{x_f - x_s}{N-1}
$$

we obtain

$$
\begin{pmatrix}
-2\alpha + V(x_1)& \alpha & 0 & \cdots& \cdots& \cdots\\
\alpha & -2\alpha + V(x_2)& \alpha & 0 & \cdots & \cdots \\
0 & \alpha & -2\alpha + V(x_3)& \alpha & 0 & \cdots \\
0 & 0 & \ddots& \ddots& \ddots
\end{pmatrix}\begin{pmatrix}
\Psi(x_1) \\
\vdots \\
\Psi(x_N)
\end{pmatrix}
$$

$$ = E_n\begin{pmatrix}
\Psi(x_1) \\
\vdots \\
\Psi(x_N)
\end{pmatrix},
$$
with

$$
\alpha = \frac{-\hbar^2}{2m\epsilon^2}.
$$

Diagonalizing, we obtain the eigenvalues and eigenstates to pass from the position basis to the Hamiltonian basis.

### Time dependent systems.

The performs dynamical simulations by considering the Suzuki-Trotter expansion of the time evolution operator,

$$
\hat{U}(t_k, t_{\text{st}}) = \prod_{j=1}^{k}\text{exp}\left[-\frac{i}{\hbar} \delta t \hat{H}(t_j) \right]+ \mathcal{O}(\delta t^2),
$$

for a number of discretization points $N_t$ and an instant $t_k$, $k\in[1, N_t]$,

$$
t_j = t_{\text{st}} + T\frac{j-1}{N_t-1} = t_{\text{st}} + \delta t (j-1).
$$

In this case, $\hat{H}(t)$ is the time dependent Hamiltonian. We consider that it can be expressed in terms of the steady state Hamiltonian operator and the position operator. In practice, $\hat{H}(t)$ is constructed in the Hamiltonian basis.

</details>

# API

The following derived types are defined in the module `OneDim_QM_NI`:
``` fortran
type, public :: steady_state_system
type, extends(steady_state_system), public :: time_dependent_system
```

## `type(steady_state_system) :: a`
Procedure | Description | Parameters
--- | --- | ---
Constructor subroutine. <br /> Usage: <br /> `call a%init(name, type, mass, start, finish, steps, args, potential[, silent])` | Initializes an steady state system and calculates its eigenvalues. | `character(len=*), intent(in) :: name`: name of the system. <br />`character(len=*), intent(in) :: type`: either `electronic` or `atomic`. If `electronic`, assumes energies in eV, lengths in $\text{\r{A}}$-s, masses in units of the electron mass $m_e$ and time scales in fs-s. If `atomic`, assumes energies in peV, lengths in $\mu$m-s, masses in units of the atomic mass unit AMU and time scales in ms-s. <br /> `real(dp), intent(in) :: mass`: mass of the particle in electronic or atomic units. <br /> `real(dp), intent(in) :: start, finish`: lenght of the 1-dimensional space in electronic or atomic units. <br />`integer, intent(in) :: steps`: number of discretization steps of the interval and number of quantum levels.<br />`real(dp), intent(in) :: args(:)`: Arguments passed to the potential function. <br />`procedure(frml) :: potential`: implementation of the potential function in electronic or atomic units. Must comply with interface [`frml`](#intfc_frml).<br />`logical, optional, intent(in) :: silent`: if `.true.`, the program will not write to terminal.

Component | Description
--- | ---
`a%Hx(N, N)` | Hamiltonian in the position basis. `N` is the number of states and discretization points.
`a%Hh(N, N)` | Hamiltonian in the Hamiltonian basis.
`a%eig(N)` | Energy eigenvalues.
`a%rot(N, N)` | Rotation from position to Hamiltonian basis.
`a%pos(N, N)` | Position matrix elements in the Hamiltonian basis.

<a id="intfc_frml"></a>
### Interface `frml`:
``` fortran
abstract interface
  function frml(position, args) result(u)
    !Represents a potential function in the position basis.
    import dp
    real(dp), intent(in) :: position
    real(dp), intent(in) :: args(:)

    real(dp) :: u
  end function frml
end interface
```

## `type(time_dependent_system) :: b`
Aditionally to all routines and components of `type(steady_state_system)`, contains:
Procedure | Description | Parameters
--- | --- | ---
Constructor subroutine. <br /> Usage: <br /> `call b%init_td(t_start, t_end, t_steps, t_args, Ht, [selected_states, selected_states_start, silent])` | Initializes a dynamics simulation. | `real(dp), intent(in) :: t_start, t_end`: lenght of the time interval in electronic or atomic units. <br />`integer, intent(in) :: steps`: number of discretization steps in the time interval. <br />`real(dp), intent(in) :: t_args(:)`: Arguments passed to the time dependent Hamiltonian function. <br />`procedure(tdep_hamil) :: Ht`: implementation of the time dependent Hamiltonian function in electronic or atomic units. Must comply with interface [`tdep_hamil`](#intfc_tdep_hamil).<br />`integer, optional, intent(in) :: selected_states, selected_states_start`: for approximations. If specified, the dynamics simulation it will only account for the states [`selected_states_start`, `selected_states_start + selected_states - 1`].<br />`logical, optional, intent(in) :: silent`: if `.true.`, the program will not write to terminal.
Dynamical simulation starter. <br /> Usage: <br /> `call b%get_tevop()` | Performs the dynamics simulation, granting acess to the `b%tevop(:, :, :)` time evolution operator.

Component | Description
--- | ---
`b%tevop(N, N, Nt)` | Time evolution operator in the Hamiltonian basis. `N` is the number of selected states (`selected_states` if it was specified in the constructor or the number of discretization points otherwise). `Nt` is the number of discretization points in the time interval. `b%tevop(:, :, ik)` is a matrix representing $\hat{U}(t_k, t_s)$, where $t_k = t_s + (t_e - t_s) \times (k-1)/(N_t-1)$.


<a id="intfc_tdep_hamil"></a>
### Interface `tdep_hamil`:
``` fortran
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
```

# Usage

The repository includes two examples. Please read them carefully. These can be run by using
```
fpm run --example Free_Particle_in_a_Box
```
```
fpm run --example Experimental_System
```

# Build

An automated build is available for [Fortran Package Manager](https://fpm.fortran-lang.org/) users. This is the recommended way to build and use OneDim-QM-NI in your projects. You can add OneDim-QM-NI to your project dependencies by including

```
[dependencies]
OneDim-QM-NI = { git="https://github.com/irukoa/OneDim-QM-NI.git" }
```
in the `fpm.toml` file.
