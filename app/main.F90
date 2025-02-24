program main

  use Iso_Fortran_ENV, only: stdout => output_unit

  use OneDim_QM_NI_kinds
  use OneDim_QM_NI_defs
  use OneDim_QM_NI_implementations
  use OneDim_QM_NI

  implicit none

  character(len=12) :: version = "v1.0.2"

  write (unit=stdout, fmt="(A)") &
    "Auxiliary library 'One Dimensional Quantum Mechanics: Numerical Implementations', "//trim(adjustl(version))//" built."

end program main
