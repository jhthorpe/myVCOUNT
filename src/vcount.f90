!------------------------------------------------------------
! vcount
!       - program that counts vibrational
!         states of a molecule as a fucntion
!         of temperature
!       
!       - This program should only be used for 
!         low temps, otherwise the algorithm
!         for iterating through VPT2 states
!         is less than efficient 
!------------------------------------------------------------
! Variables
! v             : 2D int, list of vibrational states
! E             : 1D real*8, list of vibrational energies
! N             : int, nubmer of vibratinal modes
! tol           : real*8, tollerance for when to stop 
! w             : 1D real*8, list of vibrational freq in cm-1
! X             : 2D real*8, array of coupling constants
! T             : real*8, temperature
! M             : int, number of states found
! error         : bool, exit if this becomes true

PROGRAM vcount
  USE input
  USE levels
  USE recur
  USE sort
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: X
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: E,w
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: v
  REAL(KIND=8) :: T,B,tolE,tolX,Z,E0
  INTEGER :: N,M
  LOGICAL :: error

  WRITE(*,*) "==================================================="
  WRITE(*,*) 
  WRITE(*,*) "Starting vibrational counting program"

  error = .FALSE.
  tolE = 1.0D-6
  tolX = 1.0D-6

  CALL read_input(N,T,w,X,error)
  IF (error) THEN
    CALL mem_clean(X,E,w,v)
    STOP 1 
  END IF

  CALL level_gen(N,T,w,X,E,v,tolE,tolX,M,Z,E0,B,error)
  IF (error) THEN
    CALL mem_clean(X,E,w,v)
    STOP 2
  END IF

  CALL sort_print(N,M,Z,B,E0,E,v,error)
  IF (error) THEN
    CALL mem_clean(X,E,w,v)
    STOP 2
  END IF

  CALL mem_clean(X,E,w,v)
  WRITE(*,*) 
  WRITE(*,*) "==================================================="

CONTAINS

!------------------------------------------------------------
! mem_clean
!       -cleans up memory on early exit
!------------------------------------------------------------
SUBROUTINE mem_clean(X,E,w,v)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: X
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: E,w
  INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: v
  
  IF (ALLOCATED(X)) DEALLOCATE(X)
  IF (ALLOCATED(E)) DEALLOCATE(E)
  IF (ALLOCATED(w)) DEALLOCATE(w)
  IF (ALLOCATED(v)) DEALLOCATE(v)

END SUBROUTINE mem_clean
!------------------------------------------------------------


END PROGRAM vcount
!------------------------------------------------------------
