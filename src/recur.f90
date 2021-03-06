MODULE recur
  IMPLICIT NONE


CONTAINS
!-------------------------------------------------------------------
! recur_HO
!       -recursively tablulates energies of HO levels
!-------------------------------------------------------------------
! Variables
! idx           : int, index of which oscillator we're on
! N             : int, total number of oscillators 
! w             : 1D real*8, list of vibraitonal frequencies
! Ec            : real*8, current energy
! vc            : 1D int, current list of quantum numbers
! Emax          : real*8, largest energy above E0 we can get
! key           : int, key of last level found 
! maxkey        : int, max value of key allowed, for now
! E             : 1D real*8, list of keyied energies
! v             : 2D int, array of keyied quantum vectors
! Z             : real*8, running sum
! oob           : bool, true on exit if cannot be below Emax 
! error         : bool, true on exit if problem 

RECURSIVE SUBROUTINE recur_HO(idx,N,w,Ec,vc,Emax,key,maxkey,E,v,Z,E0,B,oob,error)

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: E
  INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: v
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: w
  INTEGER, DIMENSION(0:), INTENT(INOUT) :: vc
  REAL(KIND=8), INTENT(INOUT) :: Z
  REAL(KIND=8), INTENT(IN) :: Ec, Emax, E0,B
  INTEGER, INTENT(INOUT) :: key,maxkey
  LOGICAL, INTENT(INOUT) :: oob,error
  INTEGER, INTENT(IN) :: idx,N

  REAL(KIND=8) :: En
  INTEGER :: vmax,i      

  oob = .FALSE.
  error = .FALSE.

  !WRITE(*,*) 
  !WRITE(*,*) "-----------------------------"
  !WRITE(*,*) "recur_HO called"
  !WRITE(*,*) "idx = ", idx
  !WRITE(*,*) "Z = ", Z
  !WRITE(*,*) "Ec = ", Ec
  !WRITE(*,*) "vc..." 
  !DO i=0,N-1
  !  WRITE(*,*) vc(i)
  !END DO
 
  IF (idx .EQ. 0) THEN !if we're at the last index to check
    vmax = FLOOR((Emax - Ec)/w(idx))

    IF (vmax .LT. 0) THEN !if we cannot fill the level to below Emax
      oob = .TRUE.

    ELSE !we have at least one level we can fill
      DO i=0,vmax
        IF (key .EQ. maxkey) CALL grow(maxkey,N,E,v,error)
        IF (error) RETURN
        En = Ec + w(idx)*i !this assumes vc(idx) is 0
        E(key) = En
        vc(idx) = i
        v(key,0:N-1) = vc(0:N-1)
        Z = Z + EXP(-B*(En))
        key = key + 1
        !WRITE(*,*) "idx, i", idx, i 
      END DO

    END IF

  ELSE !we are not at the last index, add in data and call again
    vmax = FLOOR((Emax - Ec)/w(idx)) !this assumes vc(idx) is 0 
    
    IF (vmax .LT. 0) THEN !we're out of bounds, go back
      oob = .TRUE.

    ELSE
      DO i=0,vmax
         !WRITE(*,*) "idx, i", idx, i 
        En = Ec + w(idx)*i
        vc(idx) = i
        CALL recur_HO(idx-1,N,w(0:N-1),En,vc(0:N-1),Emax,key,maxkey,E,&
                      v,Z,E0,B,oob,error)
        IF (oob) EXIT
      END DO 

    END IF

  END IF

  !we've finished this index, reset and return 
  vc(idx) = 0
  RETURN

END SUBROUTINE recur_HO

!-------------------------------------------------------------------
! recur_VPT2
!       -recursively tablulates energies of VTP2 levels below Emax
!-------------------------------------------------------------------
! Variables
! idx           : int, index of which oscillator we're on
! N             : int, total number of oscillators 
! w             : 1D real*8, list of vibraitonal frequencies
! X             : 2D real*8, list of coupling constants
! Ec            : real*8, current energy
! vc            : 1D int, current list of quantum numbers
! Emax          : real*8, largest energy above E0 we can get
! key           : int, key of last level found 
! maxkey        : int, max value of key allowed, for now
! E             : 1D real*8, list of keyied energies
! v             : 2D int, array of keyied quantum vectors
! Z             : real*8, running sum
! oob           : bool, true on exit if cannot be below Emax 
! error         : bool, true on exit if problem 

RECURSIVE SUBROUTINE recur_VPT2(idx,N,w,X,Ec,vc,Emax,key,maxkey,E,v,Z,E0,B,oob,error)

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: E
  INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: v
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: X
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: w
  INTEGER, DIMENSION(0:), INTENT(INOUT) :: vc
  REAL(KIND=8), INTENT(INOUT) :: Z
  REAL(KIND=8), INTENT(IN) :: Ec, Emax, E0,B
  INTEGER, INTENT(INOUT) :: key,maxkey
  LOGICAL, INTENT(INOUT) :: oob,error
  INTEGER, INTENT(IN) :: idx,N

  REAL(KIND=8) :: En,temp
  INTEGER :: vmax,i,j      

  oob = .FALSE.
  error = .FALSE.

  !calculate vmax
  temp = 0 
  DO j=0,idx-1
    temp = temp + X(idx,j)*(vc(j)+0.5)
  END DO
  DO j=idx,N-1
    temp = temp + X(j,idx)*(vc(j)+0.5)
  END DO
  temp = temp + w(idx)
  vmax = FLOOR((Emax - Ec)/temp)
 
  IF (idx .EQ. 0) THEN !if we're at the last index to check

    IF (vmax .LT. 0) THEN !if we cannot fill the level to below Emax
      oob = .TRUE.

    ELSE !we have at least one level we can fill
      DO i=0,vmax
        IF (key .EQ. maxkey) CALL grow(maxkey,N,E,v,error)
        IF (error) RETURN
 
        !get En, assuming vc = 0
        En = 0.0
        DO j=0,idx-1
          En = En + X(idx,j)*(vc(j)+0.5)*i
        END DO
        En = En + X(idx,idx)*((i+0.5)*(i+0.5)-0.25)
        DO j=idx+1,N-1
          En = En + X(j,idx)*(vc(j)+0.5)*i 
        END DO
        En = En + w(idx)*i + Ec
     
        E(key) = En
        vc(idx) = i
        v(key,0:N-1) = vc(0:N-1)
        Z = Z + EXP(-B*(En))
        key = key + 1
      END DO

    END IF

  ELSE !we are not at the last index, add in data and call again
    
    IF (vmax .LT. 0) THEN !we're out of bounds, go back
      oob = .TRUE.

    ELSE
      DO i=0,vmax
       
        !get En, assuming vc(idx) = 0 
        En = 0.0
        DO j=0,idx-1
          En = En + X(idx,j)*(vc(j)+0.5)*i
        END DO
        En = En + X(idx,idx)*((i+0.5)*(i+0.5)-0.25)
        DO j=idx+1,N-1
          En = En + X(j,idx)*(vc(j)+0.5)*i 
        END DO
        En = En + w(idx)*i + Ec

        vc(idx) = i
        CALL recur_VPT2(idx-1,N,w(0:N-1),X(0:N-1,0:N-1),En,vc(0:N-1),&
                      Emax,key,maxkey,E,v,Z,E0,B,oob,error)
        IF (oob) EXIT
      END DO 

    END IF

  END IF

  !we've finished this index, reset and return 
  vc(idx) = 0
  RETURN

END SUBROUTINE recur_VPT2

!-------------------------------------------------------------------
! grow 
!       -grows E and v arrays to 2*maxkey
!       -amusingly, this very simple subroutine takes up more
!        space than the fancy ones above...
!-------------------------------------------------------------------
! Variables
! maxkey        : int, maxkey and max index to keep
! N             : int, number of vib modes
! E             : 1D real*8, array to double and copy 
! v             : 2D int, array to double and copy
! error         : bool, true on exit if error

SUBROUTINE grow(maxkey,N,E,v,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: E
  INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: v
  INTEGER, INTENT(INOUT) :: maxkey
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: A
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: B
  INTEGER :: astat

  error = .FALSE.
  !WRITE(*,*) "Dynamically growing tables, ",maxkey,"->", 2*maxkey

  ! move over E
  ALLOCATE(A(0:2*maxkey-1),STAT=astat)
  IF (astat .NE. 0) THEN
    WRITE(*,*) "@recur:grow  -- Could not allocate new tables"
    error = .TRUE. 
    RETURN
  END IF

  A(0:maxkey-1) = E(0:maxkey-1)

  DEALLOCATE(E)
  ALLOCATE(E(0:2*maxkey-1),STAT=astat)
  IF (astat .NE. 0) THEN
    WRITE(*,*) "@recur:grow  -- Could not allocate new tables"
    error = .TRUE. 
    RETURN
  END IF
  E(0:2*maxkey-1) = A(0:2*maxkey-1)
  DEALLOCATE(A)

  ! move over v
  ALLOCATE(B(0:2*maxkey-1,0:N-1),STAT=astat)
  IF (astat .NE. 0) THEN
    WRITE(*,*) "@recur:grow  -- Could not allocate new tables"
    error = .TRUE. 
    RETURN
  END IF

  B(0:maxkey-1,0:N-1) = v(0:maxkey-1,0:N-1)

  DEALLOCATE(v)
  ALLOCATE(v(0:2*maxkey-1,0:maxkey-1),STAT=astat)
  IF (astat .NE. 0) THEN
    WRITE(*,*) "@recur:grow  -- Could not allocate new tables"
    error = .TRUE. 
    RETURN
  END IF
  v(0:2*maxkey-1,0:N-1) = B(0:2*maxkey-1,0:N-1)
  DEALLOCATE(B)
  
  maxkey = 2*maxkey 

END SUBROUTINE grow
!-------------------------------------------------------------------

END MODULE recur
