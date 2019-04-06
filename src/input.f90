MODULE input

CONTAINS

!------------------------------------------------------------
! read_input
!       -reads input... dur 
!------------------------------------------------------------
! Variables
! v             : 2D int, list of vibrational states
! N             : int, nubmer of vibratinal modes
! w             : 1D real*8, list of vibrational freq in cm-1
! X             : 2D real*8, array of coupling constants
! error         : bool, exit if this becomes true

SUBROUTINE read_input(N,T,w,X,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: X
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: w
  REAL(KIND=8),INTENT(INOUT) :: T
  INTEGER, INTENT(INOUT) :: N
  LOGICAL, INTENT(INOUT)  :: error

  CHARACTER(LEN=1024) :: word
  INTEGER :: i,j
  LOGICAL :: exists

  error = .FALSE.

  WRITE(*,*) "Reading input from vcount.dat"

  INQUIRE(file='vcount.dat',EXIST=exists)
  IF (.NOT. exists) THEN
    WRITE(*,*) "Could not find the input file: vcount.dat"
    error = .TRUE.
    RETURN
  END IF

  OPEN(file='vcount.dat',unit=100,status='old')
  READ(100,*) word, N
  READ(100,*) word, T

  ALLOCATE(w(0:N-1))
  ALLOCATE(X(0:N-1,0:N-1))

  READ(100,*)
  READ(100,*)
  DO i=0,N-1
    READ(100,*) w(i)
  END DO

  READ(100,*)
  READ(100,*)
  DO i=0,N-1
    READ(100,*) X(i,0:i)
    DO j=0,i-1
      X(j,i) = X(i,j)
    END DO
  END DO
  CLOSE(unit=100)

  WRITE(*,*)
  WRITE(*,*) "Number of vib modes         :", N
  WRITE(*,*) "Temperature (K)             :", T
  WRITE(*,*)
  WRITE(*,*) "Harmonic frequencies (cm-1)"
  DO i=0,N-1
    WRITE(*,*) w(i)
  END DO
  WRITE(*,*)
  WRITE(*,*) "Anharmonic coupling matrix"
  DO i=0,N-1
    WRITE(*,*) X(i,0:N-1)
  END DO

  WRITE(*,*)

END SUBROUTINE read_input
!------------------------------------------------------------

END MODULE input

!------------------------------------------------------------
