!------------------------------------------------------------------------------
! routine to read in setting file sys.in
!------------------------------------------------------------------------------
MODULE READIN

IMPLICIT NONE

CONTAINS

SUBROUTINE READER(latt, nstrain, delta, code, md, npt,polynmord)
USE NRTYPE

INTEGER :: latt, nstrain, polynmord, i, idum, N,ierr
REAL(DP) :: temperature, rdum
REAL(DP), DIMENSION(:), ALLOCATABLE :: delta
CHARACTER(len=8) :: code, md,npt
CHARACTER :: charac
COMPLEX(DPC) :: cdum
LOGICAL :: lopen, ldum
INTEGER, PARAMETER :: iu=12


lopen = .false.
OPEN(iu,file='elastic.in',status='old')
!------------------------------------------------------------------------------
! lattice type
!------------------------------------------------------------------------------
latt = 1
CALL DOREAD(lopen,'elastic.in',iu,'latt','=', '#', ';', 'I', &
     latt, rdum, cdum, ldum, charac, N, 1, ierr)
IF(((ierr/=0) .and. (ierr/=3)) .or. ((ierr==0) .and. (N<1))) THEN
  WRITE(*,*) "ERROR reading item ''latt'' from file sys.in."
  GOTO 150
END IF
IF(ierr==3) WRITE(*,*) "Warning : lattice type 'latt' not defined, taken as &
  Triclinic."
!------------------------------------------------------------------------------
! number of strains 
!------------------------------------------------------------------------------
nstrain = 4

CALL DOREAD(lopen,'elastic.in',iu,'nstrain','=', '#', ';', 'I', &
     nstrain, rdum, cdum, ldum, charac, N, 1, ierr)
IF((ierr/=0)  .or. ((ierr==0) .and. (N<1))) THEN
  WRITE(*,*) "Warning : ERROR reading item ''nstrain'' from file sys.in."
  GOTO 150
END IF

!------------------------------------------------------------------------------
! delta
!------------------------------------------------------------------------------
ALLOCATE(delta(nstrain))

DO i=1, nstrain
  delta(i) = 0.0
END DO

CALL DOREAD(lopen,'elastic.in',iu,'delta','=', '#', ';', 'F', &
     idum, delta, cdum, ldum, charac, N, i, ierr)
IF((ierr/=0) .or. ((ierr==0) .and. (N<1))) THEN
  WRITE(*,*) "Warning : ERROR reading item ''delta'' from file sys.in."
  GOTO 150
END IF

!------------------------------------------------------------------------------
! code 
!------------------------------------------------------------------------------

code = "VASP"

CALL DOREAD(lopen,'elastic.in',iu,'code','=', '#', ';', 'S', &
     idum, rdum, cdum, ldum, code, N, 8, ierr)
IF(((ierr/=0) .and. (ierr/=3)) .or. ((ierr==0) .and. (N<1))) THEN
  WRITE(*,*) "ERROR reading item ''code'' from file sys.in."
  GOTO 150
END IF
IF(ierr==3) WRITE(*,*) "Warning : Calc. method 'code' not defined, taken as &
  VASP."

!------------------------------------------------------------------------------
! polynmord
!------------------------------------------------------------------------------
polynmord = 1
CALL DOREAD(lopen,'elastic.in',iu,'polynmord','=', '#', ';', 'I', &
     polynmord, rdum, cdum, ldum, charac, N, 1, ierr)
IF(((ierr/=0) .and. (ierr/=3)) .or. ((ierr==0) .and. (N<1))) THEN
  WRITE(*,*) "ERROR reading item ''polynmord'' from file sys.in."
  GOTO 150
END IF
IF(ierr==3) WRITE(*,*) "Warning : Fitting order 'polynmord' not defined, &
  taken as 1."

!------------------------------------------------------------------------------
! type of calculation
!------------------------------------------------------------------------------

md = "no" 
CALL DOREAD(lopen,'elastic.in',iu,'md','=', '#', ';', 'S', &
     idum, rdum, cdum, ldum, md, N, 8, ierr)
IF(((ierr/=0) .and. (ierr/=3)) .or. ((ierr==0) .and. (N<1))) THEN
  WRITE(*,*) "ERROR reading item ''md'' from file sys.in."
  GOTO 150
END IF
IF(ierr==3) WRITE(*,*) "Warning : Calc. type 'md' not defined, taken as no." 

npt = "no"
CALL DOREAD(lopen,'elastic.in',iu,'npt','=', '#', ';', 'S', &
     idum, rdum, cdum, ldum, npt, N, 8, ierr)
IF(((ierr/=0) .and. (ierr/=3)) .or. ((ierr==0) .and. (N<1))) THEN
  WRITE(*,*) "ERROR reading item ''npt'' from file sys.in."
  GOTO 150
END IF
IF(ierr==3) WRITE(*,*) "Warning : Calc. type 'NPT' not defined, taken as no."


!---------------------------------------------

CLOSE(iu)

RETURN

150 CONTINUE
WRITE(*,*) ierr, N

STOP

END SUBROUTINE READER

END MODULE READIN
