!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   Inequality technology model
!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, IJAC, ICP(*)
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
  DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
  DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM),DFDP(NDIM,*)

  DOUBLE PRECISION x, a, alpha, c, gamma, k, m, mu, rho, s, tau, ext

  x = U(1)
  a = U(2)

  alpha   = PAR(1) != 0.2
  c       = PAR(2) != 30
  gamma   = PAR(3) != 0.3
  k       = PAR(4) != 1
  m       = PAR(5) != 100
  mu      = PAR(6) != 0.1
  rho     = PAR(7) != 0.5
  s       = PAR(8) != 0.1
  tau     = PAR(9) != 10


  ext = k*((rho*a**20)/(a**20 + 50.**20) + 1 - rho)*a**(alpha + 0.2)*x**gamma
  F(1) = tau*(x)*(x - c)*(1 - x/m) - 100*ext
  F(2) = a**(alpha + alpha + 0.2)*s*k**2*((rho*a**20)/(a**20 + 50.**20) + 1 - rho)**2*x**gamma - a*mu


END SUBROUTINE FUNC

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

SUBROUTINE STPNT(NDIM,U,PAR,T)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
  DOUBLE PRECISION, INTENT(IN) :: T

! Initialize the equation parameters
PAR(1) = 0.2
PAR(2) = 30.
PAR(3) = 0.3
PAR(4) = 1.
PAR(5) = 100.
PAR(6) = 0.1
PAR(7) = 0.5
PAR(8) = 0.1
PAR(9) = 10.

! Initialize the solution
  U(1) = 99.7155
  U(2) = 0.9861025

END SUBROUTINE STPNT

!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE BCND
END SUBROUTINE BCND

SUBROUTINE ICND
END SUBROUTINE ICND

SUBROUTINE FOPT
END SUBROUTINE FOPT

SUBROUTINE PVLS
END SUBROUTINE PVLS
