MODULE HOD
  USE cosmo
  IMPLICIT none
!!====================================================================
  CONTAINS
!!====================================================================
  DOUBLE PRECISION FUNCTION Ncen(M200c)
    USE cosmo
    USE global_var
    IMPLICIT none
    double precision, intent(IN) :: M200c
    !Ncen = 5d-1*erfc(log(Mcut/M200c)/2**0.5/sigma_Ncen)
    Ncen = 5d-1*erfc(log(Mcut/M200c)/sigma_Ncen)
    return
  END FUNCTION Ncen
!!====================================================================
  DOUBLE PRECISION FUNCTION Nsat(M200c)
    USE cosmo
    USE global_var
    IMPLICIT none
    double precision, intent(IN) :: M200c
 
    Nsat = 0.d0
    if (M200c > kappa*Mcut) then
      Nsat = ((M200c-kappa*Mcut)/M1)**alp_Nsat
      !Nsat = Nsat*Ncen(M200c)
    end if
    return
  END FUNCTION Nsat
!!====================================================================
END MODULE HOD

