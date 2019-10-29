!!$================================================================
  DOUBLE PRECISION FUNCTION calc_ug(lnx, lnM200c, ell)
    USE global_var
    USE cosmo
    USE angular_distance
    USE HOD
    IMPLICIT none
    double precision, intent(IN) :: lnx, lnM200c, ell
    double precision :: M200c, z
    double precision :: dvdz, da, chi, dchidz
    double precision :: k, u3d, u_sat, dNdz
    external u_sat, dNdz, da

    M200c = dexp(lnM200c)
    z = dexp(lnx)-1d0

    dvdz=(1d0+z)**2d0*da(z)**2d0*2998d0/Ez(z) ! h^-3 Mpc^3
    chi = da(z)*(1.d0+z)
    dchidz = dvdz/chi/chi

    k = ell/chi
    u3d = u_sat(k,M200c,z)
  
    calc_ug = (1./chi/chi)*(dNdz(z)/dchidz) &
            *abs(2.d0*Nsat(M200c)*u3d+(Nsat(M200c)*u3d)**2.d0)**0.5d0

    return
  END FUNCTION calc_ug
!!$================================================================
  DOUBLE PRECISION FUNCTION calc_uy(lnx, lnM500c, ell)
    USE global_var
    USE angular_distance
    USE cosmo
    USE mod_ptilde
    IMPLICIT none
    double precision, intent(IN) :: lnx, lnM500c, ell
    double precision :: M500c, r500, l500
    double precision :: rhoc
    double precision :: da, rombint1
    double precision :: z
    double precision :: ptilde, P0=6.41d0
    external pgnfw, rombint1

    M500c = dexp(lnM500c)
    z = dexp(lnx)-1d0
  
    ! compute omega and rhoc
    rhoc=2.775d11*Ez(z)**2d0 ! critical density in units of h^2 M_sun/Mpc^3
  
    M500c = M500c/(mass_bias*(1.d0+z)**beta)
    r500=(3d0*M500c/4d0/pi/500d0/rhoc)**(1d0/3d0) ! h^-1 Mpc
    l500=da(z)/r500
    calc_uy = 1.65d0*(h0/0.7d0)**2d0*Ez(z)**(8d0/3d0) &
            *(M500c/3d14/0.7d0)**(2d0/3d0+alp_p) &
            *P0*(0.7d0/h0)**1.5d0*ptilde(ell/l500) &
            /0.5176d0*(4d0*3.14159d0)/l500**2d0*(r500/h0)
    calc_uy = 1.04d-4*(calc_uy/50d0)  
    return
  END FUNCTION calc_uy
!!$================================================================
  DOUBLE PRECISION FUNCTION integrand_bg(lnM200c,z,ell)
  ! scale-dependent galaxy bias bg(k,z)
    use cosmo
    use sigma_z
    use mf_module
    use HOD
    use global_var
    use angular_distance
    IMPLICIT none
    double precision, intent(IN) :: lnM200c, z, ell
    double precision :: M200c, M500c, lnM500c, omega
    double precision :: dvdz, chi, da, dchidz
    double precision :: u_sat, dNdz
    external u_sat, dNdz

    M200c = dexp(lnM200c)
    call M200_to_Mdel(M200c,z,500d0,M500c)
    lnM500c = dlog(M500c)
    omega = (om0_cb+onu)*(1d0+z)**3d0/Ez(z)**2d0

    dvdz=(1d0+z)**2d0*da(z)**2d0*2998d0/Ez(z)
    chi = da(z)*(1.d0+z)
    dchidz = dvdz/chi/chi

    integrand_bg = (Ncen(M200c)+Nsat(M200c)*u_sat(ell/chi,M200c,z)) &
                   *bl_delta(lnnu_500c(lnM500c,z),500d0/omega) &
                   *dndlnMh_500c_T08(lnM500c,z)
    integrand_bg = integrand_bg*(1./chi/chi)*(dNdz(z)/dchidz)

  END FUNCTION integrand_bg
!!$================================================================
  DOUBLE PRECISION FUNCTION integrand_by(lnM200c,z,ell)
    use mf_module
    use global_var
    IMPLICIT none
    double precision, intent(IN) :: lnM200c,z,ell
    double precision :: M500c, lnM500c
    double precision :: calc_uy, lnx, omega
  
    call M200_to_Mdel(dexp(lnM200c),z,500d0,M500c)
    lnM500c = dlog(M500c)
    omega = (om0_cb+onu)*(1d0+z)**3d0/Ez(z)**2d0

    lnx = dlog(z+1d0)
    integrand_by = calc_uy(lnx,lnM500c,ell) &
                  *bl_delta(lnnu_500c(lnM500c,z),500d0/omega) &
                  *dndlnMh_500c_T08(lnM500c,z)
  
  END FUNCTION integrand_by
!!$================================================================
   DOUBLE PRECISION FUNCTION u_sat(k,M200c,z)
    USE cosmo
    USE global_var
    IMPLICIT none
    double precision, intent(IN) :: k, M200c, z
    double precision :: a, c, kr, rhoc, R200c
    double precision :: rmax_g
    real :: si, ci, si1c, ci1c, x
    double precision :: c200
    external c200
  
    c = c200(M200c,z)
    rhoc=2.775d11*Ez(z)**2d0 ! critical density in units of h^2 M_sun/Mpc^3
    R200c = (3d0*M200c/4d0/pi/200d0/rhoc)**(1d0/3d0) ! R200, h^-1 Mpc
  
    rmax_g = rmax*R200c
    a = rmax/rgs*c
  
    kr = k*rmax_g/a
    x = kr
    call cisi(x,ci,si)
    x = (1.d0+a)*kr
    call cisi(x,ci1c,si1c)
    u_sat = (sin(kr)*(si1c-si) &
           -sin(a*kr)/(1.d0+a)/kr &
           +cos(kr)*(ci1c-ci)) &
           /(log(1.d0+a)-a/(1.d0+a))
    return
  END FUNCTION u_sat
!!$================================================================
  DOUBLE PRECISION FUNCTION pgnfw(x,y) ! x=r/r500 & y=l/l500
  ! Ref: Planck Collaboration, A&A, 550, A131 (2013)
    USE cosmo
    IMPLICIT none
    double precision :: x,y
    double precision :: c=1.81d0,g=0.31d0,a=1.33d0,b=4.13d0,P0=6.41d0 ! Planck 2013
    double precision :: sinc
    external sinc
  
    pgnfw=P0*(0.7d0/h0)**1.5d0/(c*x)**g/(1d0+(c*x)**a)**((b-g)/a) &
         *x**2d0*sinc(y*x) !dsin(y*x)/(y*x)
    return
  END FUNCTION pgnfw
!!$================================================================
  DOUBLE PRECISION FUNCTION sinc(x)
    double precision :: x
    if (abs(x)<=1d-2) then
      sinc = 1.-x*x/6.+x**4/120.
    else
      sinc = dsin(x)/x
    endif
    return
  END FUNCTION sinc
!!$================================================================
  DOUBLE PRECISION FUNCTION integrand_ngz(lnM200c,z)
    use global_var
    use mf_module
    use HOD
    IMPLICIT none
    double precision, intent(IN) :: lnM200c,z
    double precision :: M200c, M500c

    M200c = dexp(lnM200c)
    call M200_to_Mdel(M200c,z,500d0,M500c)
    integrand_ngz = (Ncen(M200c)+Nsat(M200c))*dndlnMh_500c_T08(dlog(M500c),z)
    return
  END FUNCTION integrand_ngz
!!$================================================================
  DOUBLE PRECISION FUNCTION calc_ngz(z)
    IMPLICIT none
    double precision, intent(IN) :: z
    double precision :: ngz
    double precision :: lnM1, lnM2, lnM200c, dlnM
    double precision :: Mmin=1d10, Mmax=5d15
    integer :: i,n=51
    double precision :: integrand_ngz
    external integrand_ngz

    lnM1=dlog(Mmin) ! minimum mass, h^-1 Msun
    lnM2=dlog(Mmax) ! maximum mass, h^-1 Msun

    ! integrate by the trapezoidal rule
    dlnM=(lnM2-lnM1)/dble(n-1)
    ngz=0.d0
    ngz=integrand_ngz(lnM1,z)*(0.5d0*dlnM)
    ngz=ngz+integrand_ngz(lnM2,z)*(0.5d0*dlnM)
    do i=2, n-1
      lnM200c = lnM1+dble(i-1)*dlnM
      ngz = ngz+integrand_ngz(lnM200c,z)*dlnM
    enddo
    calc_ngz = ngz
    return
  END FUNCTION calc_ngz
!!$================================================================
  DOUBLE PRECISION FUNCTION dNdz(z)
    USE global_var
    IMPLICIT none
    double precision, intent(IN) :: z
    integer :: iz

    iz = int((z-z1)/((z2-z1)/nz_dndz))
    if (iz >= 0 .and. iz < nz_dndz) then
      dNdz = dndz_arr(iz)
    else
      dNdz = 0d0
    endif

    return
  END FUNCTION dNdz
!!$================================================================
  DOUBLE PRECISION FUNCTION c200(M200c,z)
    USE global_var
    IMPLICIT none
    double precision, intent(IN) :: M200c, z
    double precision :: c_arr(6) = (/37.5153d0,-1.5093d0,1.636d-2,3.66d-4,-2.89237d-5,5.32d-7/)
    integer :: i
  
    c200 = 0.d0
    if (flag_cst == 0) then 
      ! Sa ́nchez-Conde & Prada (2014) 
      do i=0,5
        c200 = c200+c_arr(i+1)*(log(M200c))**i
      end do
      c200 = c200/(1.d0+z)
    else if (flag_cst == 1) then
      ! Concentration parameter from Duffy et al. (2008)
      c200 = 5.71d0*(M200c/2d12)**(-0.084)/(1d0+z)**0.47
    else
      print *, 'ERROR: invalid value of flag_cst=',flag_cst
      stop
    endif
  
    return
  END FUNCTION c200
!!$================================================================
SUBROUTINE M200_to_Mdel(M200c,z,delc,Mdel)
  use cosmo
  IMPLICIT none
  double precision, intent(IN) :: M200c, z, delc
  double precision, intent(INOUT) :: Mdel
  double precision :: c, c200
  external c200

  c = c200(M200c,z)
  Mdel = delc/200.d0*M200c &
        *(c*xfunc(delc/200.d0*ffunc(1.0/c)))**(-3.d0)
  return
!!$================================================================
CONTAINS
!!$================================================================
  DOUBLE PRECISION FUNCTION ffunc(x) 
    double precision, intent(in) :: x
    ffunc = x**3.d0*(log(1.0+1.0/x)-1.0/(1.0+x))
    return
  END FUNCTION ffunc
!!$================================================================
  DOUBLE PRECISION FUNCTION xfunc(x)
    double precision, intent(in) :: x
    double precision :: a1 = 0.5116d0,a2=-0.4283d0,a3=-3.13d-3,a4=-3.52d-5
    double precision :: p
    p = a2 + a3*log(x)+a4*log(x)**2.d0
    xfunc = (a1*(x**(2.d0*p))+(3.0/4.0)**2)**(-0.5d0)+2.d0*x
    return
  END FUNCTION xfunc
!!$================================================================
END SUBROUTINE M200_to_Mdel
!!$================================================================
