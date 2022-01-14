!------------------------------------------------------------------------------
! routine to calculate elastic constants from npt 
!------------------------------------------------------------------------------
SUBROUTINE SSNPT(latt,cijmatrix)

USE NRTYPE
USE CALC, only : getatmnr

IMPLICIT NONE

integer :: i,j,k,ierr,latt,totatom,totstep
integer :: m,n,NRHS, LWMAX, LDA, LDB, INFO, LWORK,RANK
real(dp), dimension(:,:,:), allocatable :: vect,p,cmatrix
real(dp), dimension(:,:), allocatable :: strain,stress,A,B
real(dp), dimension(6,6), intent(inout) :: cijmatrix
real(dp), dimension(:), allocatable :: WORK,iwork,S
real(dp) :: lattice(3,3), RCOND
INTEGER, dimension(:),allocatable ::  atnr
character(len=80) :: fnam
character(len=255) :: line
LOGICAL  :: exist

!------------------------------------------------------------------------------
! read total steps
!------------------------------------------------------------------------------
ierr=0;i=0
open(14,file='XDATCAR') 
do while(ierr==0)
  read(14,'(A255)',iostat=ierr) line
  i=i+1
  if(i==7) then 
    call getatmnr(line,totatom,atnr)
  endif
end do
!write(*,*) totatom
n=i-1
totstep=n/(8+totatom)

allocate(vect(totstep,3,3))
allocate(p(totstep,3,3))
allocate(strain(totstep,6))
allocate(stress(totstep,6))
allocate(cmatrix(totstep,6,6))

!------------------------------------------------------------------------------
! read lattice and calc strain 
!------------------------------------------------------------------------------
rewind(14)
read(14,*) line
read(14,*) line
do i=1,totstep
   read(14,*) vect(i,1,:)
   read(14,*) vect(i,2,:)
   read(14,*) vect(i,3,:)
   if(i==totstep) exit
   do j=1,totatom+5
     read(14,*)
   end do
end do
close(14)
do i=1,3
  do j=1,3
    lattice(i,j)=sum(vect(:,i,j))/totstep
  enddo
enddo
lattice=ainv(lattice)

! P=R^-1*R'
do i=1,totstep
  p(i,:,:) = matmul(lattice,vect(i,:,:))
  strain(i,1)=p(i,1,1)-1
  strain(i,2)=p(i,2,2)-1
  strain(i,3)=p(i,3,3)-1
  strain(i,4)=p(i,2,3)+p(i,3,2)
  strain(i,5)=p(i,1,3)+p(i,3,1)
  strain(i,6)=p(i,1,2)+p(i,2,1)
end do

!------------------------------------------------------------------------------
! read stresses
!------------------------------------------------------------------------------
INQUIRE(file="p.dat",EXIST=exist) 
IF(exist)  CALL SYSTEM("rm p.dat")
line="grep "//"'Total+kin'"//" "//"OUTCAR "// &
             "|awk '{print $2,$3,$4,$5,$6,$7}' >"//"p.dat"
CALL SYSTEM(line,ierr)

open(15,file='p.dat',status='old')
do i=1,totstep
  read(15,*) stress(i,:)
enddo
close(15)
do i=1,6
  stress(:,i) = stress(:,i) -sum(stress(:,i))/totstep
enddo 

!------------------------------------------------------------------------------
! prepare linear equations 
!------------------------------------------------------------------------------
allocate(B(6*totstep,1))
! allocate number of unknown cij as per to system
SELECT CASE(latt)
  CASE(1) ! Triclinic
    allocate(A(totstep*6,21))
  CASE(2) ! Monoclinic
    allocate(A(totstep*6,13))
  CASE(3) ! Orthorhombic
    allocate(A(totstep*6,9))
  CASE(4) ! Tetragonal
    allocate(A(totstep*6,7))
  CASE(5) !Trigonal(3,-3)
    allocate(A(totstep*6,7))
  CASE(6) ! hexagonal
    allocate(A(totstep*6,6))
  CASE(7) ! cubic
    allocate(A(totstep*6,3))
END SELECT

do i=1,totstep
  B((6*(i-1)+1):(6*(i-1)+3),1) = stress(i,1:3)
  B(6*(i-1)+4,1) = stress(i,5)
  B(6*(i-1)+5,1) = stress(i,6)
  B(6*(i-1)+6,1) = stress(i,4)
  SELECT CASE(latt)
    CASE(1)
      A((6*(i-1)+1):(6*(i-1)+6),:) = 0
      A(6*(i-1)+1,1:6) = strain(i,1:6)
      A(6*(i-1)+2,2) = strain(i,1) ; A(6*(i-1)+2,7:11) = strain(i,2:6)
      A(6*(i-1)+3,3) = strain(i,1) ; A(6*(i-1)+3,8) = strain(i,2) 
      A(6*(i-1)+3,12:15) = strain(i,3:6)
      A(6*(i-1)+4,4) = strain(i,1) ; A(6*(i-1)+4,9) = strain(i,2) 
      A(6*(i-1)+4,13) = strain(i,3); A(6*(i-1)+4,16:18) = strain(i,4:6)
      A(6*(i-1)+5,5) = strain(i,1) ; A(6*(i-1)+5,10) = strain(i,2)
      A(6*(i-1)+5,14) = strain(i,3) ; A(6*(i-1)+5,17) = strain(i,4)
      A(6*(i-1)+5,19:20) = strain(i,5:6)
      A(6*(i-1)+6,6) = strain(i,1) ; A(6*(i-1)+6,11) = strain(i,2)
      A(6*(i-1)+6,15) = strain(i,3) ; A(6*(i-1)+6,18) = strain(i,4)
      A(6*(i-1)+5,20:21) = strain(i,5:6)
    CASE(2)
      A((6*(i-1)+1):(6*(i-1)+6),:) = 0
      A(6*(i-1)+1,1) = strain(i,1) ;  A(6*(i-1)+1,2) = strain(i,2)
      A(6*(i-1)+1,3) = strain(i,3) ;  A(6*(i-1)+1,4) = strain(i,5)
      A(6*(i-1)+2,2) = strain(i,1) ;  A(6*(i-1)+2,5) = strain(i,2)
      A(6*(i-1)+2,6) = strain(i,3) ;  A(6*(i-1)+2,7) = strain(i,5)
      A(6*(i-1)+3,3) = strain(i,1) ;  A(6*(i-1)+3,6) = strain(i,2)
      A(6*(i-1)+3,8) = strain(i,3) ;  A(6*(i-1)+3,9) = strain(i,5)
      A(6*(i-1)+4,10) = strain(i,4) ;  A(6*(i-1)+4,11) = strain(i,6)
      A(6*(i-1)+5,4) = strain(i,1) ;  A(6*(i-1)+5,7) = strain(i,2)
      A(6*(i-1)+5,9) = strain(i,3) ;  A(6*(i-1)+5,12) = strain(i,5)
      A(6*(i-1)+6,11) = strain(i,4) ;  A(6*(i-1)+6,13) = strain(i,6)
    CASE(3)
      A((6*(i-1)+1):(6*(i-1)+6),:) = 0
      A(6*(i-1)+1,1:3) = strain(i,1:3) 
      A(6*(i-1)+2,2) = strain(i,1) ; A(6*(i-1)+2,4:5) = strain(i,2:3)
      A(6*(i-1)+3,3) = strain(i,1) ; A(6*(i-1)+3,5:6) = strain(i,2:3)
      A(6*(i-1)+4,7) = strain(i,4) ; A(6*(i-1)+5,8) = strain(i,5)
      A(6*(i-1)+6,9) = strain(i,6)
    CASE(4)
      A((6*(i-1)+1):(6*(i-1)+6),:) = 0
      A(6*(i-1)+1,1:3) = strain(i,1:3) ; A(6*(i-1)+1,7) = strain(i,6)
      A(6*(i-1)+2,2) = strain(i,1) ; A(6*(i-1)+2,1) = strain(i,2) 
      A(6*(i-1)+2,3) = strain(i,3) ; A(6*(i-1)+2,7) = 0-strain(i,6)
      A(6*(i-1)+3,3) = strain(i,1)+strain(i,2);A(6*(i-1)+3,4) = strain(i,3)
      A(6*(i-1)+4,5) = strain(i,4) ;A(6*(i-1)+5,5) = strain(i,5)
      A(6*(i-1)+6,6) = strain(i,6)
    CASE(5)
      A((6*(i-1)+1):(6*(i-1)+6),:) = 0 
      A(6*(i-1)+1,1:5) = strain(i,1:5)
      A(6*(i-1)+2,1) = strain(i,2) ; A(6*(i-1)+2,2) = strain(i,1)
      A(6*(i-1)+2,3) = strain(i,3) ; A(6*(i-1)+2,4) =0- strain(i,4)
      A(6*(i-1)+2,5) =0- strain(i,5) 
      A(6*(i-1)+3,3) = strain(i,1) + strain(i,2) 
      A(6*(i-1)+3,6) = strain(i,3)
      A(6*(i-1)+4,4) = strain(i,1) - strain(i,2)
      A(6*(i-1)+4,5) =0- strain(i,6) ; A(6*(i-1)+4,7) = strain(i,4)
      A(6*(i-1)+5,4) = strain(i,6) ; A(6*(i-1)+5,7) = strain(i,5)
      A(6*(i-1)+5,5) = strain(i,1) - strain(i,2)
      A(6*(i-1)+6,1) = strain(i,6)/2 ; A(6*(i-1)+6,2) =0- strain(i,6)/2
      A(6*(i-1)+6,4) = strain(i,5) ; A(6*(i-1)+6,5) =0- strain(i,4)
    CASE(6)
      A(6*(i-1)+1,1) = strain(i,1) ;  A(6*(i-1)+1,2) = strain(i,2)
      A(6*(i-1)+1,3) = strain(i,3) ;  A(6*(i-1)+1,4:6) = 0
      A(6*(i-1)+2,1) = strain(i,2) ;  A(6*(i-1)+2,2) = strain(i,1)
      A(6*(i-1)+2,3) = strain(i,3) ;  A(6*(i-1)+2,4:6) = 0 
      A(6*(i-1)+3,1:2) = 0 ; A(6*(i-1)+3,5:6) = 0
      A(6*(i-1)+3,3)=strain(i,1)+strain(i,2) ; A(6*(i-1)+3,4)=strain(i,3)
      A(6*(i-1)+4,:) = 0 ; A(6*(i-1)+4,4)=strain(i,4)
      A(6*(i-1)+5,:) = 0 ; A(6*(i-1)+5,5)=strain(i,5)
      A(6*(i-1)+6,:) = 0 ; A(6*(i-1)+6,6)=strain(i,6)
    CASE(7)
      A(6*(i-1)+1,1) = strain(i,1) ;  A(6*(i-1)+1,2) = strain(i,2)+strain(i,3)
      A(6*(i-1)+1,3) = 0
      A(6*(i-1)+2,1) = strain(i,2) ;  A(6*(i-1)+2,2) = strain(i,1)+strain(i,3)
      A(6*(i-1)+2,3) = 0
      A(6*(i-1)+3,1) = strain(i,3) ; A(6*(i-1)+3,2) = strain(i,1)+strain(i,2)
      A(6*(i-1)+3,3) = 0
      A(6*(i-1)+4,:) = 0 ; A(6*(i-1)+4,3)=strain(i,4)
      A(6*(i-1)+5,:) = 0 ; A(6*(i-1)+5,3)=strain(i,5)
      A(6*(i-1)+6,:) = 0 ; A(6*(i-1)+6,3)=strain(i,6)
  END SELECT

end do

!------------------------------------------------------------------------------
! solve strain-stress overdetermined linear equations
!------------------------------------------------------------------------------
LWMAX=10000 ; ALLOCATE(WORK(LWMAX))
M=SIZE(A,1) ; N=SIZE(A,2) ; NRHS=1 ; LDA=M ; LDB=M ; RCOND = -1.0
allocate(s(M)) ; allocate(iwork(14*M))
! Query the optimal workspace.
LWORK = -1
CALL DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK,LWORK, IWORK, INFO )
LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
! Solve the equations A*X = B
CALL DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK,LWORK, IWORK, INFO )
! Check for convergence
IF( INFO.GT.0 ) THEN
  WRITE(*,*)'The algorithm computing SVD failed to converge;'
  WRITE(*,*)'the least squares solution could not be computed.'
  STOP
END IF

!write(*,*) B(1:4,1)
!------------------------------------------------------------------------------
! map cij
!------------------------------------------------------------------------------
cijmatrix=0
SELECT CASE(latt)
  CASE(1)
    cijmatrix(1,1)=B(1,1);cijmatrix(1,2)=B(2,1);cijmatrix(1,3)=B(3,1)
    cijmatrix(1,4)=B(4,1);cijmatrix(1,5)=B(5,1);cijmatrix(1,6)=B(6,1)
    cijmatrix(2,2)=B(7,1);cijmatrix(2,3)=B(8,1);cijmatrix(2,4)=B(9,1)
    cijmatrix(2,5)=B(10,1);cijmatrix(2,6)=B(11,1);cijmatrix(3,3)=B(12,1)
    cijmatrix(3,4)=B(13,1);cijmatrix(3,5)=B(14,1);cijmatrix(3,6)=B(15,1)
    cijmatrix(4,4)=B(16,1);cijmatrix(4,5)=B(17,1);cijmatrix(4,6)=B(18,1)
    cijmatrix(5,5)=B(19,1);cijmatrix(5,6)=B(20,1);cijmatrix(6,6)=B(21,1)
    cijmatrix(2,1)=cijmatrix(1,2);cijmatrix(3,1)=cijmatrix(1,3)
    cijmatrix(3,2)=cijmatrix(2,3);cijmatrix(4,1)=cijmatrix(1,4)
    cijmatrix(4,2)=cijmatrix(2,4);cijmatrix(4,3)=cijmatrix(3,4)
    cijmatrix(5,1)=cijmatrix(1,5);cijmatrix(5,2)=cijmatrix(2,5)
    cijmatrix(5,3)=cijmatrix(3,5);cijmatrix(5,4)=cijmatrix(4,5)
    cijmatrix(6,1)=cijmatrix(1,6);cijmatrix(6,2)=cijmatrix(2,6)
    cijmatrix(6,3)=cijmatrix(3,6);cijmatrix(6,4)=cijmatrix(4,6)
    cijmatrix(6,5)=cijmatrix(5,6)
  CASE(2)
    cijmatrix(1,1)=B(1,1);cijmatrix(1,2)=B(2,1);cijmatrix(1,3)=B(3,1)
    cijmatrix(1,5)=B(4,1);cijmatrix(2,2)=B(5,1);cijmatrix(2,3)=B(6,1)
    cijmatrix(2,5)=B(7,1);cijmatrix(3,3)=B(8,1);cijmatrix(3,5)=B(9,1)
    cijmatrix(4,4)=B(10,1);cijmatrix(4,6)=B(11,1);cijmatrix(5,5)=B(12,1)
    cijmatrix(6,6)=B(13,1)
    cijmatrix(2,1)=cijmatrix(1,2);cijmatrix(3,1)=cijmatrix(1,3)
    cijmatrix(3,2)=cijmatrix(2,3);cijmatrix(5,1)=cijmatrix(1,5)
    cijmatrix(5,2)=cijmatrix(2,5);cijmatrix(5,3)=cijmatrix(3,5)
    cijmatrix(6,4)=cijmatrix(4,6)
  CASE(3)
    cijmatrix(1,1)=B(1,1);cijmatrix(1,2)=B(2,1);cijmatrix(1,3)=B(3,1)
    cijmatrix(2,2)=B(4,1);cijmatrix(2,3)=B(5,1);cijmatrix(3,3)=B(6,1)
    cijmatrix(4,4)=B(7,1);cijmatrix(5,5)=B(8,1);cijmatrix(6,6)=B(9,1)
    cijmatrix(2,1)=cijmatrix(1,2);cijmatrix(3,1)=cijmatrix(1,3)
    cijmatrix(3,2)=cijmatrix(2,3)
  CASE(4)
    cijmatrix(1,1)=B(1,1);cijmatrix(1,2)=B(2,1);cijmatrix(1,3)=B(3,1)
    cijmatrix(3,3)=B(4,1);cijmatrix(4,4)=B(5,1);cijmatrix(6,6)=B(6,1)
    cijmatrix(1,6)=B(7,1)
    cijmatrix(2,1)=cijmatrix(1,2);cijmatrix(3,1)=cijmatrix(1,3)
    cijmatrix(2,3)=cijmatrix(1,3);cijmatrix(3,2)=cijmatrix(1,3)
    cijmatrix(2,6)=0-cijmatrix(1,6) ; cijmatrix(5,5)=cijmatrix(4,4)
  CASE(5)
    cijmatrix(1,1)=B(1,1);cijmatrix(1,2)=B(2,1);cijmatrix(1,3)=B(3,1)
    cijmatrix(1,4)=B(4,1);cijmatrix(1,5)=B(5,1);cijmatrix(3,3)=B(6,1)
    cijmatrix(4,4)=B(7,1)
    cijmatrix(2,1)=cijmatrix(1,2);cijmatrix(3,1)=cijmatrix(1,3)
    cijmatrix(2,3)=cijmatrix(1,3);cijmatrix(3,2)=cijmatrix(1,3)
    cijmatrix(2,4)=0-cijmatrix(1,4);cijmatrix(2,5)=0-cijmatrix(1,5)
    cijmatrix(4,1)=cijmatrix(1,4)
    cijmatrix(4,2)=0-cijmatrix(1,4);cijmatrix(4,6)=0-cijmatrix(1,5)
    cijmatrix(5,1)=cijmatrix(1,5);cijmatrix(5,2)=0-cijmatrix(1,5)
    cijmatrix(5,5)=cijmatrix(4,4);cijmatrix(5,6)=cijmatrix(1,4)
    cijmatrix(6,4)=0-cijmatrix(1,5);cijmatrix(6,5)=cijmatrix(1,4)
    cijmatrix(6,6)=(cijmatrix(1,1)-cijmatrix(1,2))/2
  CASE(6)
    cijmatrix(1,1)=B(1,1);cijmatrix(2,2)=B(1,1);cijmatrix(1,3)=B(3,1)
    cijmatrix(1,2)=B(2,1);cijmatrix(3,3)=B(4,1);cijmatrix(4,4)=B(5,1)
    cijmatrix(6,6)=B(6,1)
    cijmatrix(2,1)=cijmatrix(1,2);cijmatrix(3,1)=cijmatrix(1,3)
    cijmatrix(2,3)=cijmatrix(1,3);cijmatrix(3,2)=cijmatrix(1,3)
    cijmatrix(5,5)=cijmatrix(4,4)
  CASE(7)
    cijmatrix(1,1)=B(1,1);cijmatrix(2,2)=B(1,1);cijmatrix(3,3)=B(1,1)
    cijmatrix(1,2:3)=B(2,1);cijmatrix(2,1)=B(2,1);cijmatrix(2,3)=B(2,1)
    cijmatrix(3,1:2)=B(2,1);cijmatrix(4,4)=B(3,1);cijmatrix(5,5)=B(3,1)
    cijmatrix(6,6)=B(3,1)
END SELECT
! change unit kbar to GPa
cijmatrix = cijmatrix*(-0.1)


WRITE(*,*) "Elastic Stiffness Constant Matrix : "
DO i=1,6
  WRITE(*,'(6F16.8)') cijmatrix(i,:)
END DO

!------------------------------------------------------------------------------
! done
!------------------------------------------------------------------------------
deallocate(vect)
deallocate(p)
deallocate(strain)
deallocate(stress)
deallocate(cmatrix)

RETURN

contains
! function to inverse array A
FUNCTION AINV(A)
real(dp), dimension(:,:), intent(in) :: A
real(dp), dimension(size(A,1),size(A,2)) :: Ainv

real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
integer, dimension(size(A,1)) :: ipiv   ! pivot indices
integer :: n, info

! External procedures defined in LAPACK
!        external DGETRF
!        external DGETRI

! Store A in Ainv to prevent it from being overwritten by LAPACK
Ainv = A
n = size(A,1)

! DGETRF computes an LU factorization of a general M-by-N matrix A
! using partial pivoting with row interchanges.
call DGETRF(n, n, Ainv, n, ipiv, info)

if (info /= 0) then
  stop 'Matrix is numerically singular!'
end if

! DGETRI computes the inverse of a matrix using the LU factorization
! computed by DGETRF.
call DGETRI(n, Ainv, n, ipiv, work, n, info)

if (info /= 0) then
   stop 'Matrix inversion failed!'
end if
end function AINV

END SUBROUTINE SSNPT
