MODULE sort
  

CONTAINS

!------------------------------------------------------------
! sort_print
!       -Sorts and prints E and v 
!        in increasing energies relative to E0
!------------------------------------------------------------
! Variables
! N             : int, nubmer of vib modes
! M             : int, number of states to sort
! Z             : real*8, Z of partion function
! B             : real*8, beta of partition function
! E             : 1D real*8, list of energies
! v             : 2D int, list of vibrational states
! error         : bool, true on exit if error

SUBROUTINE sort_print(N,M,Z,B,E0,E,v,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: E
  INTEGER, DIMENSION(0:,0:), INTENT(INOUT) :: v
  REAL(KIND=8), INTENT(IN) :: Z,B,E0
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N,M

  INTEGER :: i

  error  = .FALSE.

  !WRITE(*,*) "M is", M
  !WRITE(*,*) "Priot to sort, E is..."
  !WRITE(*,*) E(0:M-1)
  CALL quicksort(0,M-1,N,E(0:M-1),v(0:M-1,0:N-1))

  WRITE(*,*) " B is", B

  
  !WRITE(*,*) EXP(-(10.-E0))
  !WRITE(*,*) EXP(-B*(10.-E0))
  !WRITE(*,*) EXP(-B*(10.-E0))/Z
  WRITE(*,*) "vib state, energy, population"
  DO i=0,M-1
    WRITE(*,*) v(i,0:N-1),E(i), EXP(-B*E(i))/Z
  END DO 

END SUBROUTINE sort_print

!------------------------------------------------------------
! quicksort
!       - quicksort alrogithm that will affect both E and v
!------------------------------------------------------------
! Variables
! lo            : int, lowest index
! hi            : int, highest index
! N             : int, length of vector part of v
! E             : 1D real*8, values to be sorted
! v             : 2D int, vecotors to be sorted along E

RECURSIVE SUBROUTINE quicksort(lo,hi,N,E,v)
  IMPLICIT NONE
 
  INTEGER, DIMENSION(lo:hi,0:N-1), INTENT(INOUT) :: v
  REAL(KIND=8), DIMENSION(lo:hi), INTENT(INOUT) :: E
  INTEGER,INTENT(IN) :: lo,hi,N

  INTEGER :: p

  IF (lo .LT. hi) THEN
    CALL partition(p,lo,hi,N,E(lo:hi),v(lo:hi,0:N-1))
    CALL quicksort(lo,p-1,N,E(lo:p-1),v(lo:p-1,0:N-1))
    CALL quicksort(p+1,hi,N,E(p+1:hi),v(p+1:hi,0:N-1))
  END IF

END SUBROUTINE quicksort

!------------------------------------------------------------
!  partition
!       - partition funciton for quicksort
!------------------------------------------------------------
! Variables
! lo, hi        : int, lowest/heighest indices
! N             : int, number of vib modes
! E             : 1D real*8, list of values to sort along
! v             : 2D int, list of vectors to get sorted as well

SUBROUTINE partition(p,lo,hi,N,E,v)
  IMPLICIT NONE

  INTEGER, DIMENSION(lo:hi,0:N-1), INTENT(INOUT) :: v
  REAL(KIND=8), DIMENSION(lo:hi), INTENT(INOUT) :: E
  INTEGER, INTENT(INOUT) :: p
  INTEGER,INTENT(IN) :: lo,hi,N

  REAL(KIND=8) :: pivot,temp
  INTEGER, DIMENSION(0:N-1) :: v1
  INTEGER :: i,j

  !pivot = E((lo + hi)/2)
  pivot = E(hi)
  i = lo
  DO j=lo, hi-1
    IF (E(j) .LT. pivot) THEN
      temp = E(j) 
      E(j) = E(i) 
      E(i) = temp
      v1(0:N-1) = v(j,0:N-1) 
      v(j,0:N-1) = v(i,0:N-1)
      v(i,0:N-1) = v1(0:N-1)
      i = i + 1
    END IF 
  END DO
  temp = E(hi) 
  E(hi) = E(i) 
  E(i) = temp
  v1(0:N-1) = v(hi,0:N-1) 
  v(hi,0:N-1) = v(i,0:N-1)
  v(i,0:N-1) = v1(0:N-1)
  p = i
  RETURN

  !WRITE(*,*) "Pivot is", E((lo+hi)/2), (lo+hi)/2
  !i = lo 
  !j = hi 
!
 ! DO WHILE (E(i) .LT. pivot .AND. i .LT. 100) 
 !   i = i + 1
 ! END DO
!
!  DO WHILE (E(j) .GT. pivot .AND. j .GT. 0)
!    j = j - 1
!  END DO
!
!  IF (i .GT. 99 .OR. j .LT. 0) STOP "fuck"
!
!  IF (i .GE. j) THEN
!    p = j
!    RETURN
!
!  ELSE
!    p = (lo+hi)/2
!    temp = E(j)
!    E(j) = E(i) 
!    E(i) = temp
!    v1(0:N-1) = v(j,0:N-1) 
!    v(j,0:N-1) = v(i,0:N-1)
!    v(i,0:N-1) = v1(0:N-1)
!    RETURN
!  END IF

END SUBROUTINE partition

!------------------------------------------------------------



END MODULE sort
