!------------------------------------------------------------------------------
! extract stresses from subdirectories
!------------------------------------------------------------------------------

SUBROUTINE  EXTRACT(code,md,nstrain,numdist,stress,devt)

USE NRTYPE

IMPLICIT NONE

INTEGER :: i,j,nstrain, numdist,dimens,lsys

REAL(DP), intent(out) :: stress(0:nstrain,1:numdist,1:6) 
!REAL(DP), DIMENSION(:), ALLOCATABLE :: stresspass,devt
REAL(DP), DIMENSION(:), ALLOCATABLE :: stresspass
REAL(DP), intent(out) :: devt(0:nstrain,1:numdist,1:6)

CHARACTER(len=20) :: fnam, dirnam
CHARACTER(len=80) :: cmd,line
CHARACTER(len=8) :: code, md

LOGICAL :: found
INTEGER :: OU
! initilize 
ou = 62
dimens = 6

ALLOCATE(stresspass(dimens))
!ALLOCATE(devt(dimens))

OPEN(ou,file='stress.out',status='new')
OPEN(unit=66,file='devt.dat',status='new')
devt = 0.0
stress = 0.0 ; stresspass = 0
!------------------------------------------------------------------------------
! vasp
!------------------------------------------------------------------------------

IF(trim(adjustl(code))=='vasp' .or. trim(adjustl(code))=='Vasp' .or. &
   trim(adjustl(code))=='VASP' .or. trim(adjustl(code))=='V') THEN
  IF(trim(adjustl(md))=='no' .or. trim(adjustl(md))=='No' .or. &
     trim(adjustl(md))=='NO' .or. trim(adjustl(md))=='N') THEN
!stresses from static calc.
    DO i = 0, nstrain
      DO j = 1, numdist
        WRITE(dirnam,'(A1,I1,A1,I1)') "s",i,"d",j
        IF(I==0) WRITE(dirnam,'(A1,I1,A1,I1)') "s",i,"d",j-1
        cmd="grep "//"'Total CPU'"//" "//trim(dirnam)//"/OUTCAR >dummy.ch"
        CALL SYSTEM(cmd,lsys);OPEN(unit=65,file='dummy.ch',status='unknown') 
        READ(65,*,IOSTAT=lsys); CLOSE(65) ; CALL SYSTEM("rm dummy.ch")
!check OUTCAR
        IF(lsys/=0) WRITE(*,*) "Warning : Calculation improperly terminated &
          in directory: ", dirnam
        cmd="grep "//"'in kB'"//" "//trim(dirnam)//"/OUTCAR "//"|tail -1 &
             |awk '{print $3,$4,$5,$6,$7,$8}' >dummy.dat"
        CALL SYSTEM(cmd,lsys)
        OPEN(unit=64,file='dummy.dat',status='unknown')        
        READ(64,*,IOSTAT=lsys) stress(i,j,:)
        IF(lsys==0) THEN
          CLOSE(64); CALL SYSTEM("rm dummy.dat")
          WRITE(ou,'(2I4,6F14.4)') i,j, stress(i,j,:)
        ELSE
          WRITE(*,*) "Stress reading failed in directory: ",dirnam 
        END IF
        WRITE(66,*) i,j,devt(i,j,:) ! this redundant line used to create devt.dat
        IF(i==0) EXIT
      END DO
    END DO
!stresses from md calc.
  ELSE
    DO i = 0, nstrain
      DO j = 1, numdist
        WRITE(dirnam,'(A1,I1,A1,I1)') "s",i,"d",j
        IF(I==0) WRITE(dirnam,'(A1,I1,A1,I1)') "s",i,"d",j-1
!check OUTCAR
        cmd="grep "//"'Total CPU'"//" "//trim(dirnam)//"/OUTCAR >dummy.ch"
        CALL SYSTEM(cmd,lsys);OPEN(unit=65,file='dummy.ch',status='unknown')
        READ(65,*,IOSTAT=lsys); CLOSE(65) ; CALL SYSTEM("rm dummy.ch")
        IF(lsys/=0) WRITE(*,*) "Warning : Calculation improperly terminated &
          in directory: ", dirnam
        fnam=trim(dirnam)//".dat"
        cmd="grep "//"'Total+kin'"//" "//trim(dirnam)//"/OUTCAR "// &
             "|awk '{print $2,$3,$4,$5,$6,$7}' >"//fnam
        CALL SYSTEM(cmd,lsys)
!calc. average, devt
        CALL BLDEV(fnam,dimens,stresspass,devt(i,j,:))
        WRITE(ou,*) i,j, stresspass
        WRITE(66,*) i,j, devt(i,j,:)
        IF(i==0) EXIT
      END DO
    END DO
  END IF
ELSE
  WRITE(*,*) "'code = vasp' is not well defined in sys.in."
END IF

!------------------------------------------------------------------------------
! QE
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! GULP
!------------------------------------------------------------------------------

CLOSE(ou)
CLOSE(66)

END SUBROUTINE EXTRACT

