!------------------------------------------------------------------------------
! this routine calculate average, and error, dev. from blocking method
!------------------------------------------------------------------------------

SUBROUTINE BLDEV(fnam,ntype,average,devt)

USE NRTYPE

IMPLICIT NONE

INTEGER :: n,m,ntype,nprim,i,j,k,ierr

REAL(dp), DIMENSION(:,:), ALLOCATABLE :: array1,array2,array3,array4
REAL(dp), DIMENSION(:), ALLOCATABLE ::  c
REAL(dp), intent(out) ::  average(ntype)
REAL(dp), optional,intent(out)  ::  devt(ntype)
CHARACTER(len=10) :: fnam

INTEGER :: IU
IU = 25
OPEN(iu,file=trim(adjustl(fnam)), status='old')

!calc. number of lines  
ierr=0
n=0
DO WHILE (ierr .eq. 0)
  READ(iu,*,IOSTAT=ierr)
  n = n + 1
END DO
n = n -1
!write(*,'(A16,I10)') "number of type: ", ntype
!write(*,'(A16,I10)') "number of data: ", n 


ALLOCATE(array1(ntype,n))
ALLOCATE(array2(ntype,n/2))

ALLOCATE(c(ntype))
!ALLOCATE(devt(ntype))

! read in data and calc. average
average=0.0
REWIND(iu)
DO i=1,n
  READ(iu,*) array1(1:ntype,i)
  average(1:ntype) = average(1:ntype) + array1(1:ntype,i)
END DO
average = average / n

! WRITE(*,*) "Average :             ", average

! blockment data
OPEN(116, file ='devtmp.dat', status='replace')
DO i=1, n / 2
  array2(1:ntype,i) = 0.5*(array1(1:ntype,2*i-1) + array1(1:ntype,2*i))
  c(1:ntype) = c(1:ntype) + (array2(1:ntype,i) - average(1:ntype))**2
  WRITE(116,*) array2(1:ntype,i)
END DO
c = 2*c / n
CLOSE(116)
devt = sqrt( c/(n/2-1) )

 j=1
 DO WHILE (m.GT.4)
   OPEN(116, file ='devtmp.dat', status='old')
   m=n/2
   ALLOCATE(array3(ntype,m))
   ALLOCATE(array4(ntype,m/2))
   DO i=1, m
     READ(116,*) array3(1:ntype,i)
   ENDDO
   CLOSE(116)
   OPEN(116, file ='devtmp.dat', status='replace')
   DO i=1, m / 2
     array4(1:ntype,i) = 0.5*(array3(1:ntype,2*i-1) + array3(1:ntype,2*i))
     c(1:ntype) = c(1:ntype) + (array4(1:ntype,i) - average(1:ntype))**2
     WRITE(116,*) array4(1:ntype,i) 
   ENDDO
   c = 2*c / m
   CLOSE(116)
   devt(:) = SQRT(c(:)/(m-1))
!   ddelta(:) = (delta(:,2)-delta(:,1))/delta(:,1)
!   mxddev = MAXVAL(ddelta)
!   WRITE(*,*) "Block size :         ", j
!   WRITE(*,*) "Stastical error :    ", delta(:,2),ddelta(:)
!   WRITE(*,*) "Stastical error (%): ", delta(:,2)/average
   DEALLOCATE(array3)
   DEALLOCATE(array4)
   n=m ; j=j+1  !; delta(:,1) = delta(:,2)
 END DO


DEALLOCATE(array1)
DEALLOCATE(array2)
!DEALLOCATE(average)
DEALLOCATE(C)
!DEALLOCATE(devt)

CLOSE(25)

RETURN

END SUBROUTINE BLDEV
