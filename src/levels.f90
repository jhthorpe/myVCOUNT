MODULE levels
  USE recur

CONTAINS

!------------------------------------------------------------
! level_gen()
!       - generates E,v lists from w and X
!------------------------------------------------------------
! N             : int, nubmer of vibratinal modes
! T             : real*8, temperature
! w             : 1D real*8, list of vibrational freq in cm-1
! X             : 2D real*8, array of coupling constants
! E             : 1D real*8, list of vibrational energies
! v             : 2D int, list of vibrational states
! tolE          : real*8, tol for limit to exp(-E*B)
! tolX          : real*8, tol for when to use harmonic
! M             : int, number of states found
! error         : bool, exit if this becomes true

SUBROUTINE level_gen(N,T,w,X,E,v,tolE,tolX,M,Z,E0,B,error)

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: E
  INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: v
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: X
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: w
  REAL(KIND=8), INTENT(INOUT) :: Z,E0,B
  REAL(KIND=8), INTENT(IN) :: T,tolE,tolX
  INTEGER, INTENT(INOUT) :: M
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N

  REAL(KIND=8) :: kb
  INTEGER :: treat

  error = .FALSE.
  kb = 8.6173303E-5 * 8065.5
  B = 1.0D0/(kb*T)

  WRITE(*,*) "---------------------------------------------------"
  WRITE(*,*) "Begining level counting"

  treat = level_treat(X,tolX)
  IF (treat .EQ. 0) THEN
   !CALL level_VPT2()
    IF (error) RETURN
  ELSE 
    CALL level_HO(N,w(0:N-1),v,E,tolE,B,M,Z,E0,error)
    IF (error) RETURN
  END IF

  WRITE(*,*) 

END SUBROUTINE level_gen

!------------------------------------------------------------
! level_HO
!       -generates energy levels using harmonic oscillator
!        approximation
!------------------------------------------------------------
! Variables
! N             : int, nubmer of vibratinal modes
! w             : 1D real*8, list of vibrational freq in cm-1
! v             : 2D int, list of vibrational states
! E             : 1D real*8, list of vibrational energies
! tolE          : real*8, max total energy to consider 
! B             : real*8, beta, 1/kb*T
! error         : bool, exit if this becomes true

SUBROUTINE level_HO(N,w,v,E,tolE,B,M,Z,E0,error)
  USE recur
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: E
  INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: v
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: w 
  REAL(KIND=8), INTENT(INOUT) :: Z,E0
  REAL(KIND=8), INTENT(IN) :: tolE,B
  INTEGER, INTENT(INOUT) :: M
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N

  INTEGER, DIMENSION(0:N-1) :: vc
  REAL(KIND=8) :: Emax
  INTEGER :: key, maxkey
  LOGICAL :: oob

  error = .FALSE.
  
  E0 = 0.5*SUM(w(0:N-1)) !harmonic ZPE  
  Emax = -1.0*LOG(tolE)/B + E0 !max energy relative to E0
  !Emax = -1.0*LOG(tolE)/B !max energy relative to E0
  WRITE(*,*) "Generating harmonic vib states below (cm-1): ", Emax
  WRITE(*,*) "β = ", B

  maxkey = 128  !2^7just initial guess
  !maxkey = 1 !2^0, testing purposes
  ALLOCATE(v(0:maxkey-1,0:N-1))
  ALLOCATE(E(0:maxkey-1))

  E(0) = 0.0D0
  v(0,0:N-1) = 0
  Z = 0.0
  key = 0
  vc = 0
  oob = .FALSE.
  CALL recur_HO(N-1,N,w,0.0D0,vc(0:N-1),Emax,key,maxkey,E,v,Z,E0,B,oob,error) 
  WRITE(*,*) 
  WRITE(*,*) "There are ", key, "vib states below Emax"
  !WRITE(*,*) "v0, v1, energy"
  !DO i=0,key-1
  !  WRITE(*,*) v(i,0:N-1), E(i)
  !END DO
  M = key
  !WRITE(*,*) "TESTING TESTING TESTING"

END SUBROUTINE level_HO

!------------------------------------------------------------
! level_treat
!       -determines if we should use harmonic or anharmonic
!        treatment
!       - returns 0 if anharmonic, 1 if harmonic
!------------------------------------------------------------
! Variables
! N             : int, number of vib modes
! X             : 2D real*8, list of couplings
! tolX          : real*8, tol for anharmonic vs harmonic  

INTEGER FUNCTION level_treat(X,tolX)
  IMPLICIT NONE
  
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: X
  REAL(KIND=8), INTENT(IN) :: tolX

  IF (MAXVAL(X) .GT. tolX .OR. MINVAL(X) .LT. -tolX) THEN
    WRITE(*,*) 
    WRITE(*,*) "Using VPT2 treatment"
    level_treat = 0
  ELSE 
    WRITE(*,*) 
    WRITE(*,*) "No Xij greater than", tolX
    WRITE(*,*) "Using harmonic treatment"
    level_treat = 1
  END IF 

  

END FUNCTION level_treat
!------------------------------------------------------------

END MODULE levels

