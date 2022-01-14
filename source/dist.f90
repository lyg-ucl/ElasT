!------------------------------------------------------------------------------
! generate distorted structures from POSCAR 
!                                                         
! Input ::  POSCAR  vasp5 format                          
! Output::  POSCAR_s1d1   s-strain, d-dsitortion          
!------------------------------------------------------------------------------
SUBROUTINE  DIST(latt,nstrain,numdist,delta)

USE NRTYPE

IMPLICIT  NONE

INTEGER :: latt,nstrain,nr,maxatom,i,j,k,m,numdist

REAL(dp), DIMENSION(3,3) :: vect,newvect
REAL(dp), intent(out)  :: delta(nstrain)
REAL(dp), DIMENSION(8,6,3,3) :: defvect
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: atmcord  ! data array

CHARACTER(len=20) :: fnam 
CHARACTER(len=255) :: line
CHARACTER(len=120) :: argmt

LOGICAL  :: found

INTEGER  :: iu, ou
iu=11 ; ou=63

!remove old files
INQUIRE(file="POSCAR_s1d1",EXIST=found)
IF(found)  CALL SYSTEM("rm POSCAR_*")

!------------------------------------------------------------------------------
!Loop over strain
!------------------------------------------------------------------------------
DO j = 1, nstrain
  ! distortion matrix
  defvect(1,1,1,1) = 1 + delta(j)  
  defvect(1,1,1,2) = 0
  defvect(1,1,1,3) = 0
  defvect(1,1,2,1) = 0
  defvect(1,1,2,2) = 1
  defvect(1,1,2,3) = 0
  defvect(1,1,3,1) = 0
  defvect(1,1,3,2) = 0
  defvect(1,1,3,3) = 1
  defvect(1,2,1,1) = 1
  defvect(1,2,1,2) = 0
  defvect(1,2,1,3) = 0
  defvect(1,2,2,1) = 0
  defvect(1,2,2,2) = 1 + delta(j)
  defvect(1,2,2,3) = 0
  defvect(1,2,3,1) = 0
  defvect(1,2,3,2) = 0
  defvect(1,2,3,3) = 1
  defvect(1,3,1,1) = 1
  defvect(1,3,1,2) = 0
  defvect(1,3,1,3) = 0
  defvect(1,3,2,1) = 0
  defvect(1,3,2,2) = 1
  defvect(1,3,2,3) = 0
  defvect(1,3,3,1) = 0
  defvect(1,3,3,2) = 0
  defvect(1,3,3,3) = 1 + delta(j)
  defvect(1,4,1,1) = 1
  defvect(1,4,1,2) = 0
  defvect(1,4,1,3) = 0
  defvect(1,4,2,1) = 0
  defvect(1,4,2,2) = 1
  defvect(1,4,2,3) = delta(j) / 2
  defvect(1,4,3,1) = 0
  defvect(1,4,3,2) = delta(j) / 2
  defvect(1,4,3,3) = 1
  defvect(1,5,1,1) = 1
  defvect(1,5,1,2) = 0
  defvect(1,5,1,3) = delta(j) / 2
  defvect(1,5,2,1) = 0
  defvect(1,5,2,2) = 1
  defvect(1,5,2,3) = 0
  defvect(1,5,3,1) = delta(j) / 2
  defvect(1,5,3,2) = 0
  defvect(1,5,3,3) = 1 
  defvect(1,6,1,1) = 1
  defvect(1,6,1,2) = delta(j) / 2
  defvect(1,6,1,3) = 0
  defvect(1,6,2,1) = delta(j) / 2
  defvect(1,6,2,2) = 1
  defvect(1,6,2,3) = 0
  defvect(1,6,3,1) = 0
  defvect(1,6,3,2) = 0
  defvect(1,6,3,3) = 1 

  defvect(2,1,1,1) = 1 + delta(j)
  defvect(2,1,1,2) = 0
  defvect(2,1,1,3) = 0
  defvect(2,1,2,1) = 0
  defvect(2,1,2,2) = 1
  defvect(2,1,2,3) = delta(j) / 2
  defvect(2,1,3,1) = 0
  defvect(2,1,3,2) = delta(j) / 2
  defvect(2,1,3,3) = 1
  defvect(2,2,1,1) = 1
  defvect(2,2,1,2) = delta(j) / 2
  defvect(2,2,1,3) = 0
  defvect(2,2,2,1) = delta(j) / 2
  defvect(2,2,2,2) = 1 + delta(j) 
  defvect(2,2,2,3) = 0
  defvect(2,2,3,1) = 0
  defvect(2,2,3,2) = 0
  defvect(2,2,3,3) = 1
  defvect(2,3,1,1) = 1
  defvect(2,3,1,2) = 0
  defvect(2,3,1,3) = delta(j) / 2
  defvect(2,3,2,1) = 0
  defvect(2,3,2,2) = 1
  defvect(2,3,2,3) = 0
  defvect(2,3,3,1) = delta(j) / 2
  defvect(2,3,3,2) = 0
  defvect(2,3,3,3) = 1
  defvect(2,4,1,1) = 1
  defvect(2,4,1,2) = 0
  defvect(2,4,1,3) = 0
  defvect(2,4,2,1) = 0
  defvect(2,4,2,2) = 1
  defvect(2,4,2,3) = 0
  defvect(2,4,3,1) = 0
  defvect(2,4,3,2) = 0
  defvect(2,4,3,3) = 1 + delta(j) 

  defvect(3,1,1,1) = 1 + delta(j) / 2
  defvect(3,1,1,2) = 0
  defvect(3,1,1,3) = 0
  defvect(3,1,2,1) = 0
  defvect(3,1,2,2) = 1 + delta(j) / 2
  defvect(3,1,2,3) = delta(j) / 2
  defvect(3,1,3,1) = 0
  defvect(3,1,3,2) = delta(j) / 2
  defvect(3,1,3,3) = 1
  defvect(3,2,1,1) = 1 - delta(j) / 2
  defvect(3,2,1,2) = 0
  defvect(3,2,1,3) = delta(j) / 2
  defvect(3,2,2,1) = 0
  defvect(3,2,2,2) = 1 + delta(j) / 2
  defvect(3,2,2,3) = 0
  defvect(3,2,3,1) = delta(j) / 2
  defvect(3,2,3,2) = 0
  defvect(3,2,3,3) = 1
  defvect(3,3,1,1) = 1
  defvect(3,3,1,2) = delta(j) / 2
  defvect(3,3,1,3) = 0
  defvect(3,3,2,1) = delta(j) / 2
  defvect(3,3,2,2) = 1
  defvect(3,3,2,3) = 0
  defvect(3,3,3,1) = 0
  defvect(3,3,3,2) = 0
  defvect(3,3,3,3) = 1 + delta(j)

  defvect(4,1,1,1) = 1 + delta(j)
  defvect(4,1,1,2) = delta(j) / 2
  defvect(4,1,1,3) = 0
  defvect(4,1,2,1) = delta(j) / 2
  defvect(4,1,2,2) = 1
  defvect(4,1,2,3) = 0
  defvect(4,1,3,1) = 0
  defvect(4,1,3,2) = 0
  defvect(4,1,3,3) = 1 + delta(j)
  defvect(4,2,1,1) = 1
  defvect(4,2,1,2) = 0
  defvect(4,2,1,3) = 0
  defvect(4,2,2,1) = 0
  defvect(4,2,2,2) = 1 + delta(j) / 2
  defvect(4,2,2,3) = delta(j) / 2
  defvect(4,2,3,1) = 0
  defvect(4,2,3,2) = delta(j) / 2
  defvect(4,2,3,3) = 1 + delta(j)

  defvect(5,1,1,1) = 1 
  defvect(5,1,1,2) = 0 
  defvect(5,1,1,3) = 0 
  defvect(5,1,2,1) = 0 
  defvect(5,1,2,2) = 1 + delta(j) 
  defvect(5,1,2,3) = delta(j) / 2 
  defvect(5,1,3,1) = 0 
  defvect(5,1,3,2) = delta(j) / 2 
  defvect(5,1,3,3) = 1 
  defvect(5,2,1,1) = 1 
  defvect(5,2,1,2) = delta(j) / 2 
  defvect(5,2,1,3) = 0 
  defvect(5,2,2,1) = delta(j) / 2 
  defvect(5,2,2,2) = 1 
  defvect(5,2,2,3) = 0 
  defvect(5,2,3,1) = 0 
  defvect(5,2,3,2) = 0 
  defvect(5,2,3,3) = 1 + delta(j) 

  defvect(6,1,1,1) = 1 + delta(j)
  defvect(6,1,1,2) = 0
  defvect(6,1,1,3) = 0
  defvect(6,1,2,1) = 0
  defvect(6,1,2,2) = 1
  defvect(6,1,2,3) = delta(j) / 2
  defvect(6,1,3,1) = 0
  defvect(6,1,3,2) = delta(j) / 2
  defvect(6,1,3,3) = 1
  defvect(6,2,1,1) = 1
  defvect(6,2,1,2) = 0
  defvect(6,2,1,3) = 0
  defvect(6,2,2,1) = 0
  defvect(6,2,2,2) = 1
  defvect(6,2,2,3) = 0
  defvect(6,2,3,1) = 0
  defvect(6,2,3,2) = 0
  defvect(6,2,3,3) = 1 + delta(j)

  defvect(7,1,1,1) = 1 + delta(j)
  defvect(7,1,1,2) = delta(j) / 2
  defvect(7,1,1,3) = 0
  defvect(7,1,2,1) = delta(j) / 2
  defvect(7,1,2,2) = 1
  defvect(7,1,2,3) = 0
  defvect(7,1,3,1) = 0
  defvect(7,1,3,2) = 0
  defvect(7,1,3,3) = 1
! loop over distortions
  DO i = 1, numdist
     WRITE(fnam,"(A8,I1,A1,I1)") "POSCAR_s", j,"d",i
     open(iu, file='POSCAR', status='old')
     open(ou, file = fnam, status='new')
     DO k = 1, 2
      READ(iu,'(A255)') line
      WRITE(ou,'(A20)') line
     END DO
     DO k = 1, 3
        READ(iu,*) vect(k,:)
     END DO
     DO k = 1, 3
        DO m = 1, 3
           newvect(k,m) = vect(k,1)*defvect(latt,i,1,m) + &
           vect(k,2)*defvect(latt,i,2,m) + vect(k,3)*defvect(latt,i,3,m)
        END DO
     END DO
     DO k = 1, 3
        WRITE(ou,*) (newvect(k,m), m=1,3)
     END DO
     READ(iu,'(A80)') line
     WRITE(ou,'(A80)') line
     READ(iu,'(A80)') line
     WRITE(ou,'(A80)') line
     CALL getatmnr(line, maxatom)  ! subrooutine to transfer string to number
     READ(iu,'(A80)') line
     WRITE(ou,'(A80)') line
     ALLOCATE(atmcord(maxatom,3))
     DO nr = 1, maxatom
        READ(iu,*) atmcord(nr,:)
        WRITE(ou,*) atmcord(nr,:)
     END DO
     DEALLOCATE(atmcord)
     close(iu)
     close(ou)
  END DO
END DO

CONTAINS
 SUBROUTINE getatmnr(line,atomnr)
  CHARACTER(len=255),intent(IN):: line
  INTEGER,intent(OUT) ::atomnr
  INTEGER :: ios,atnr,atnrpos
  CHARACTER(len=255):: line2,atomnrch

 ios=0
 atomnr=0
 line2=line

 DO WHILE (ios==0)
    read(line2,*,IOSTAT=ios) atnr
    IF (ios==0) then
       atomnr=atomnr+atnr
       write(atomnrch,*) atnr
       atomnrch=ADJUSTL(atomnrch)
       atnrpos=Index(line2,trim(atomnrch))+len_trim(atomnrch)
! Index finds the location of a substring in another string,returns 0 if
! not found
! len_trim counts the number of characters including blanks
! this sentence locates the point in the line
       line2=""
       line2(atnrpos:255)=line(atnrpos:255)
       IF (len_trim(line2)<=0) then
          ios=10
       END IF
    END IF
 END DO

 END SUBROUTINE getatmnr

END SUBROUTINE  DIST

