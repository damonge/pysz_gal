MODULE global_var
  integer :: nl
  double precision :: mass_bias, alp_p, beta
  double precision :: z1 = 1.d-5, z2 = 4d0
  integer :: nz_dndz
  double precision :: Mmin = 1d10, Mmax = 5d15
  integer :: nz = 1001, nz_2h = 100
  integer :: nm = 51
  integer :: nz_tll = 2001, nm_tll = 201
  integer :: flag_nu = 0 ! 0:'tot', 1:'cb'
  integer :: flag_tll = 0 ! 0: not calc, 1: calc
  integer :: pk_nz, pk_nk
  double precision :: pi = 3.14159265359d0
  double precision :: Mcut, M1, kappa, sigma_Ncen, alp_Nsat ! HOD params
  double precision :: rmax, rgs ! satellite dist. params
  integer :: flag_cst = 0 ! 0: Sa ́nchez-Conde & Prada (2014)
  double precision, allocatable :: dndz_arr(:)
END MODULE global_var
