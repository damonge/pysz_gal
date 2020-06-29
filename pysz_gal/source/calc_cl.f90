#define MAXSIZE 4096
SUBROUTINE calc_cl(h0_in, obh2_in, och2_in, mnu_in, mass_bias_in,&
                   Mcut_in, M1_in, kappa_in, sigma_Ncen_in, alp_Nsat_in,&
                   rmax_in, rgs_in,&
                   pk_nk_in, pk_nz_in, k_arr, pk_arr,&
                   dndz_in, nz_dndz_in,&
                   z1_in,z2_in,z1_ng,z2_ng,&
                   nl_in, ell_arr, cl_gg, cl_gy, tll, ng, flag_nu_in, flag_tll_in,&
                   nm_in, nz_in)
  !$ USE omp_lib
  use cosmo
  use global_var
  use linearpk_z
  use sigma_z
  use angular_distance
  use mod_ptilde
  IMPLICIT none
  integer, intent(IN) :: pk_nk_in, pk_nz_in, nl_in, nz_dndz_in
  double precision, intent(IN) :: h0_in, obh2_in, och2_in, mnu_in
  double precision, intent(IN) :: mass_bias_in
  double precision, intent(IN) :: Mcut_in, M1_in, kappa_in, sigma_Ncen_in, alp_Nsat_in
  double precision, intent(IN) :: rmax_in, rgs_in
  double precision, intent(IN) :: z1_in, z2_in, z1_ng, z2_ng
  double precision, dimension(0:nz_dndz_in-1) :: dndz_in
  double precision, dimension(0:pk_nk_in-1,0:pk_nz_in-1), intent(IN) :: k_arr, pk_arr
  double precision, dimension(0:nl_in-1), intent(INOUT) :: ell_arr
  double precision, dimension(0:nl_in-1,0:1), intent(INOUT) :: cl_gg, cl_gy
  double precision, dimension(0:nl_in*2-1,0:nl_in*2-1), intent(INOUT) :: tll
  double precision, intent(INOUT) :: ng
  integer, intent(IN) :: flag_nu_in, flag_tll_in
  integer, intent(IN) :: nm_in, nz_in
  double precision :: calc_ng
  external calc_ng

  ! precision of mass and z integration
  nz = nz_in
  nm = nm_in

  nl = nl_in
  z1 = z1_in
  z2 = z2_in

  ! flag for neutrino prescription
  flag_nu = flag_nu_in

  ! Calc Tll or not
  flag_tll = flag_tll_in

  ! read dNdz
  nz_dndz = nz_dndz_in
  allocate(dndz_arr(nz_dndz))
  dndz_arr(1:nz_dndz) = dndz_in(0:nz_dndz-1)

  ! read in linear P(k,z)
  pk_nk = pk_nk_in
  pk_nz = pk_nz_in
  call open_linearpk(pk_nk,pk_nz,k_arr,pk_arr)

  ! input parameters
  !! cosmological parameters
  h0 = h0_in
  obh2 = obh2_in
  och2 = och2_in
  om0_cb = (obh2+och2)/h0**2
  mnu = mnu_in
  onu = mnu/93.14/h0**2
  ode0 = 1.d0-om0_cb-onu
  w = -1d0

  !! tSZ parameters
  mass_bias = mass_bias_in
  ! fixed parameters
  alp_p = 0.12d0
  beta = 0.d0

  !! galaxy parameters
  Mcut = Mcut_in
  M1 = M1_in
  kappa = kappa_in
  sigma_Ncen = sigma_Ncen_in
  alp_Nsat = alp_Nsat_in
  rmax = rmax_in
  rgs = rgs_in

  !!!! preface !!!!
  !! ptilde
  call setup_ptilde

  !! fit sigma^2(R) to Chebyshev polynomials
  call compute_sigma2

  !! compute and tabulate da(z) 
  call setup_da
  !!!!!!!!!!!!!!!!!

  ! calculate Cls ang Ng(z1 < z < z2)
  !$OMP barrier
  call calc_cl_gg_gy(ell_arr,cl_gg,cl_gy,tll)
  ng = calc_ng(z1_ng,z2_ng)
  !$OMP barrier

  call close_linearpk
  call close_sigma
  call close_ptilde
  deallocate(dndz_arr)
  !$OMP barrier

!===============================================================
CONTAINS
!===============================================================
  SUBROUTINE calc_cl_gg_gy(ell_arr,cl_gg,cl_gy,tll)
    !$ USE omp_lib
    USE global_var
    IMPLICIT none
    double precision, intent(INOUT) :: cl_gg(0:nl-1,0:1)
    double precision, intent(INOUT) :: cl_gy(0:nl-1,0:1)
    double precision, intent(INOUT) :: tll(0:nl*2-1,0:nl*2-1)
    double precision, intent(IN) :: ell_arr(:)
    double precision, allocatable :: cls_th_1h(:,:), cls_th_2h(:,:), tll_th(:,:,:)
    integer :: i, j, nth, ith
    double precision :: lnx, lnx1, lnx2, dlnx
    double precision :: lnM, lnM1, lnM2, dlnM
    double precision :: fac1
    double precision :: intg_1h(nl*2), intg_2h(nl*2), intg_tll(nl*2,nl*2)

    cl_gg(:,:) = 0.d0
    cl_gy(:,:) = 0.d0
    intg_tll(:,:) = 0.d0

    lnx1=dlog(1d0+z1); lnx2 = dlog(1d0+z2)
    dlnx = (lnx2-lnx1)/nz
  
    lnM1=dlog(Mmin); lnM2=dlog(Mmax)
    dlnM = (lnM2-lnM1)/nm
    
    ! 1-halo term
    !$OMP parallel private(i,j,lnx,fac1,lnM,intg_1h,intg_2h,ith), shared(cls_th_1h,cls_th_2h,nth,dlnx)
    nth = omp_get_num_threads()
    ith = omp_get_thread_num()
    !$OMP single
    allocate(cls_th_1h(nl*2,nth),cls_th_2h(nl*2,nth))
    cls_th_1h(:,:) = 0.d0
    cls_th_2h(:,:) = 0.d0
    !$OMP end single
    !$OMP do
    do j = 1, nz+1
      lnx = lnx1+dlnx*(j-1)
      fac1 = 1.d0
      if (j == 1 .or. j == nz+1) fac1 = 0.5d0
      ! mass integration
      ! trapezoidal
      intg_1h = (integrand_1h(lnM1,lnx,ell_arr)+integrand_1h(lnM2,lnx,ell_arr))*(0.5d0*dlnM)
      do i = 2, nm
        lnM = lnM1+dlnM*(i-1)
        intg_1h = intg_1h+integrand_1h(lnM,lnx,ell_arr)*dlnM
      end do
      ! ! simpson
      ! intg_1h = integrand_1h(lnM1,lnx,ell_arr)+integrand_1h(lnM2,lnx,ell_arr)
      ! do i = 1, nm, 2
      !   intg_1h = intg_1h+4.d0*integrand_1h(lnM1+dlnM*i,lnx,ell_arr)
      ! end do
      ! do i = 2, nm, 2
      !   intg_1h = intg_1h+2.d0*integrand_1h(lnM1+dlnM*i,lnx,ell_arr)
      ! end do
      ! intg_1h = intg_1h*dlnM/3d0
      ! qgaus
      ! call qgaus2_n20_arr_gy(integrand_1h,lnM1,lnM2,intg_1h,lnx,ell_arr)
      cls_th_1h(:,ith+1) = cls_th_1h(:,ith+1)+intg_1h*dlnx*fac1
    end do 
    !$OMP end do
    !$OMP barrier

    !$OMP barrier
    ! 2-halo term
    ! dlnx = (lnx2-lnx1)/nz_2h
    !$OMP do
    ! do j = 1, nz_2h+1
    do j = 1, nz+1
      lnx = lnx1+dlnx*(j-1)
      fac1 = 1.d0
      if (j == 1 .or. j == nz+1) fac1 = 0.5d0
      call calc_integrand_2h(lnx,ell_arr,intg_2h)
      cls_th_2h(:,ith+1) = cls_th_2h(:,ith+1)+intg_2h*dlnx*fac1
    end do
    !$OMP end do
    !$OMP barrier

    !$OMP end parallel

    do i = 1, nl
      ! 1h term
      cl_gg(i-1,0) = sum(cls_th_1h(i,1:nth))
      cl_gy(i-1,0) = sum(cls_th_1h(i+nl,1:nth))
      ! 2h term
      cl_gg(i-1,1) = sum(cls_th_2h(i,1:nth))
      cl_gy(i-1,1) = sum(cls_th_2h(i+nl,1:nth))
    end do

    !$OMP single
    deallocate(cls_th_1h,cls_th_2h)
    !$OMP end single

    ! Tll
    if (flag_tll == 1) then
      !$OMP parallel private(i,j,lnx,fac1,lnM,intg_tll,ith), shared(tll_th,nth)
      nth = omp_get_num_threads()
      ith = omp_get_thread_num()
      !$OMP single
      allocate(tll_th(nl*2,nl*2,nth))
      tll_th(:,:,:) = 0.d0
      !$OMP end single
      !$OMP do
      do j = 1, nz+1
        lnx = lnx1+dlnx*(j-1)
        fac1 = 1.d0
        if (j == 1 .or. j == nz+1) fac1 = 0.5d0
        ! mass integration
        ! trapezoidal
        intg_tll = (integrand_tll(lnM1,lnx,ell_arr)+integrand_tll(lnM2,lnx,ell_arr))*(0.5d0*dlnM)
        do i = 2, nm
          lnM = lnM1+dlnM*(i-1)
          intg_tll = intg_tll+integrand_tll(lnM,lnx,ell_arr)*dlnM
        end do
        tll_th(:,:,ith+1) = tll_th(:,:,ith+1)+intg_tll*dlnx*fac1
      end do
      !$OMP end do
      !$OMP barrier
      !$OMP end parallel

      do i = 1, nl*2
        do j = 1, nl*2
          tll(i-1,j-1) = sum(tll_th(i,j,1:nth))
        end do
      end do

      !$OMP single
      deallocate(tll_th)
      !$OMP end single
 
    end if

  END SUBROUTINE calc_cl_gg_gy
!!===============================================================
  FUNCTION integrand_1h(lnM200c, lnx, ell_arr)
    USE global_var
    USE cosmo
    USE mf_module
    USE angular_distance
    IMPLICIT none
    double precision :: integrand_1h(nl*2)
    double precision, intent(IN) :: lnx, lnM200c, ell_arr(nl)
    double precision :: ngz, uy, ug, ell
    double precision :: z, dvdz, fac, da
    double precision :: calc_uy, calc_ug, calc_ngz
    double precision :: M500c, lnM500c
    integer :: i
    external calc_uy, calc_ug, calc_ngz, da

    integrand_1h(:) = 0d0

    z = dexp(lnx)-1d0
    call M200_to_Mdel(dexp(lnM200c),z,500d0,M500c)
    lnM500c = dlog(M500c)

    ngz = calc_ngz(z)
    do i = 1, nl
      ell = ell_arr(i)
      uy = calc_uy(lnx,lnM500c,ell)
      ug = calc_ug(lnx,lnM200c,ell)
      integrand_1h(i) = ug*ug/ngz**2d0
      integrand_1h(i+nl) = ug*uy/ngz
    end do

    z = dexp(lnx)-1d0
    dvdz = (1d0+z)**2d0*da(z)**2d0*clight/Ez(z) ! h^-3 Mpc^3
    fac = dvdz*(1+z)*dndlnMh_200c_T08(lnM200c,z)
    !fac = dvdz*(1+z)*dndlnMh_500c_T08(lnM500c,z)

    integrand_1h = integrand_1h*fac

    return  
  END FUNCTION integrand_1h
!===============================================================
  SUBROUTINE calc_integrand_2h(lnx, ell_arr, intg_2h)
    USE global_var
    USE cosmo
    USE mf_module
    USE angular_distance
    USE linearpk_z
    IMPLICIT none
    double precision, intent(IN) :: lnx, ell_arr(nl)
    double precision, intent(INOUT) :: intg_2h(nl*2)
    double precision :: ell
    double precision :: z, dvdz, fac, da, chi
    double precision :: by, bg, lnM1, lnM2, lnM, dlnM
    double precision :: integrand_by, integrand_bg
    double precision :: linear_pk
    double precision :: pk1, pk2, pk_z, dz
    double precision :: calc_ngz, ngz
    integer :: il, iz, im
    external linear_pk, integrand_by, integrand_bg, calc_ngz

    lnM1 = dlog(Mmin); lnM2 = dlog(Mmax)

    intg_2h = 0d0
    z = dexp(lnx)-1d0
    dz = (z2-z1)/(pk_nz-1)
    iz = int((z-z1)/dz)

    ngz = calc_ngz(z)

    do il = 1, nl
      ell = ell_arr(il)

      ! SZ term
      call qgaus2(integrand_by,lnM1,lnM2,by,z,ell) ! Gaussian quadrature, SZ bias

      ! galaxy term
      dlnM = (lnM2-lnM1)/nm
      bg = (integrand_bg(lnM1,z,ell)+integrand_bg(lnM2,z,ell))*(0.5d0*dlnM)
      do im = 2, nm-1
        lnM = lnM1+dlnM*(im-1)
        bg = bg+integrand_bg(lnM,z,ell)*dlnM
     end do
  
      dvdz=(1d0+z)**2d0*da(z)**2d0*clight/Ez(z)
      chi = da(z)*(1.d0+z)
      if (iz < pk_nz-1) then
        pk1 = linear_pk((ell+0.5)/chi,iz+1)
        pk2 = linear_pk((ell+0.5)/chi,iz+2)
        ! pk_z = (pk2-pk1)/dz*(z-dz*iz)+pk1
        pk_z = (pk2-pk1)/dz*(z-(z1+dz*iz))+pk1
      else
        pk_z = linear_pk((ell+0.5)/chi,pk_nz) 
      end if

      fac = pk_z*dvdz*(1d0+z)
  
      intg_2h(il) = bg*bg*fac/ngz**2
      intg_2h(il+nl) = bg*by*fac/ngz
    end do

    return  
  END SUBROUTINE calc_integrand_2h
!===============================================================
 FUNCTION integrand_tll(lnM200c, lnx, ell_arr)
    USE global_var
    USE cosmo
    USE mf_module
    USE angular_distance
    IMPLICIT none
    double precision :: integrand_tll(nl*2,nl*2), integrand_1h(nl*2,1)
    double precision, intent(IN) :: lnx, lnM200c, ell_arr(nl)
    double precision :: uy, ug, ngz, ell
    double precision :: z, dvdz, fac, da
    double precision :: calc_uy, calc_ug, calc_ngz
    double precision :: M500c, lnM500c
    integer :: i
    external calc_uy, calc_ug, calc_ngz, da

    integrand_tll(:,:) = 0d0

    z = dexp(lnx)-1d0
    ngz = calc_ngz(z)
    call M200_to_Mdel(dexp(lnM200c),z,500d0,M500c)
    lnM500c = dlog(M500c)

    do i = 1, nl
      ell = ell_arr(i)
      uy = calc_uy(lnx,lnM500c,ell)
      ug = calc_ug(lnx,lnM200c,ell)
      integrand_1h(i,1) = ug*ug/ngz**2
      integrand_1h(i+nl,1) = ug*uy/ngz
    end do

    integrand_tll = matmul(integrand_1h,transpose(integrand_1h))

    z = dexp(lnx)-1d0
    dvdz = (1d0+z)**2d0*da(z)**2d0*clight/Ez(z) ! h^-3 Mpc^3
    fac = dvdz*(1+z)*dndlnMh_500c_T08(lnM500c,z)

    integrand_tll = integrand_tll*fac

    return  
  END FUNCTION integrand_tll
!===============================================================
END SUBROUTINE calc_cl
