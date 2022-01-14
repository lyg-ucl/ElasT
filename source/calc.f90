!------------------------------------------------------------------------------
! routine to calc. sigma matrix, elastic properties 
! References:
! 1) check Born mechanical stability. M. Born and K. Huang, Dynamics Theory of
!    Crystal Lattices (Oxford University Press, 1954).
! 2) Félix Mouhat and François-Xavier Coudert, Phys. Rev. B 90, 224104 
!------------------------------------------------------------------------------
MODULE CALC

USE NRTYPE

IMPLICIT NONE

! used by rhocalc, getelemtnm, getatmnr 
CHARACTER(len=2), dimension(:), ALLOCATABLE :: elemnm
INTEGER, DIMENSION(:), ALLOCATABLE          :: atnr

CONTAINS

SUBROUTINE RHOCALC(fnam,rho)
INTEGER                  :: i,j,ntype,nr,maxatom
REAL(DP)                 :: a, b, c, scaler, vol, mass, alpha, beta, gamma
REAL(DP), INTENT(OUT)    :: rho
REAL(DP), DIMENSION(:), ALLOCATABLE :: elemms
REAL(DP), DIMENSION(3,3) :: vect
CHARACTER(len=20)        :: fnam, title
CHARACTER(len=255)       :: line

INTEGER                  :: IU
iu=11

OPEN(iu,file=trim(adjustl(fnam)),status='old')
READ(iu,*) 
READ(iu,*) scaler
DO i=1,3
  READ(iu,*) vect(i,:)
END DO
a = sqrt(vect(1,1)*vect(1,1) + vect(1,2)*vect(1,2) + vect(1,3)*vect(1,3))
b = sqrt(vect(2,1)*vect(2,1) + vect(2,2)*vect(2,2) + vect(2,3)*vect(2,3))
c = sqrt(vect(3,1)*vect(3,1) + vect(3,2)*vect(3,2) + vect(3,3)*vect(3,3))
alpha = ACOS((vect(2,1)*vect(3,1) + vect(2,2)*vect(3,2) + vect(2,3)*vect(3,3))/(b*c))
beta = ACOS((vect(1,1)*vect(3,1) + vect(1,2)*vect(3,2) + vect(1,3)*vect(3,3))/(a*c))
gamma = ACOS((vect(1,1)*vect(2,1) + vect(1,2)*vect(2,2) + vect(1,3)*vect(2,3))/(a*b))
vol = a*b*c*scaler*scaler*scaler* & 
    sqrt(1 - COS(alpha)*COS(alpha) - COS(beta)*COS(beta) - &
    COS(gamma)*COS(gamma) - 2*COS(alpha)*COS(beta)*COS(gamma))

READ(iu,'(A255)') line
CALL getelemtnm(line,elemnm,ntype)      ! subroutine to get element name and number of kinds of element

READ(iu,'(A255)') line
CALL getatmnr(line, maxatom, atnr)  ! subrooutine to transfer string to number

CLOSE(iu)

ALLOCATE(elemms(ntype))
DO i=1,ntype
  DO j=1,109
    IF(trim(adjustl(elemnm(i)))==nmarray(j)) elemms(i)=msarray(j)
  END DO
END DO

mass = 0.0
DO i=1, ntype
  mass = mass +atnr(i)*elemms(i)
END DO
rho = mass/vol*1.6605389

RETURN

END SUBROUTINE RHOCALC

!------------------------------------------------------------------------------
! routine to calculate sigma matrix
! sigma_ij = Cij * strain
!------------------------------------------------------------------------------
SUBROUTINE SIGMACALC(latt,sigmatrix,stress,devt,sigmdevt)

INTEGER :: i,j,latt
REAL(DP), DIMENSION(:,:,:), INTENT(OUT) :: sigmatrix,sigmdevt
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: stress,devt


DO i=0, size(stress,1)-1
  DO j=1, size(stress,2)
    IF(i==0) EXIT
    stress(i,j,:) = stress(0,1,:) - stress(i,j,:)
  END DO 
END DO


!initiate 
sigmatrix = 0 ; sigmdevt = 0
DO i=1, size(stress,1)-1
  SELECT CASE(latt)
!in array stress, yz(sigma4)->position5 zx(sig5)->6 xy(sig6)->4
!hexagonal
    CASE(6)
      sigmatrix(i,1,1) = stress(i,1,1)
      sigmatrix(i,1,2) = stress(i,1,2)
      sigmatrix(i,1,3) = (stress(i,1,3) + stress(i,2,1) + stress(i,2,2))/3
      sigmatrix(i,2,1) = stress(i,1,2) 
      sigmatrix(i,2,2) = stress(i,1,1) 
      sigmatrix(i,2,3) = (stress(i,1,3) + stress(i,2,1) + stress(i,2,2))/3
      sigmatrix(i,3,1) = (stress(i,1,3) + stress(i,2,1) + stress(i,2,2))/3
      sigmatrix(i,3,2) = (stress(i,1,3) + stress(i,2,1) + stress(i,2,2))/3
      sigmatrix(i,3,3) = stress(i,2,3) 
      sigmatrix(i,4,4) = stress(i,1,5)  
      sigmatrix(i,5,5) = stress(i,1,5)
      sigmatrix(i,6,6) = (sigmatrix(i,1,1) - sigmatrix(i,1,2))/2
!devt
      sigmdevt(i,1,1) = devt(i,1,1)
      sigmdevt(i,1,2) = devt(i,1,2)
      sigmdevt(i,1,3) = (devt(i,1,3) + devt(i,2,1) + devt(i,2,2))/3
      sigmdevt(i,2,1) = devt(i,1,2)
      sigmdevt(i,2,2) = devt(i,1,1)
      sigmdevt(i,2,3) = (devt(i,1,3) + devt(i,2,1) + devt(i,2,2))/3
      sigmdevt(i,3,1) = (devt(i,1,3) + devt(i,2,1) + devt(i,2,2))/3
      sigmdevt(i,3,2) = (devt(i,1,3) + devt(i,2,1) + devt(i,2,2))/3
      sigmdevt(i,3,3) = devt(i,2,3)
      sigmdevt(i,4,4) = devt(i,1,5)
      sigmdevt(i,5,5) = devt(i,1,5)
      sigmdevt(i,6,6) = (sigmdevt(i,1,1) + sigmdevt(i,1,2))/2
!cubic
    CASE(7)
      sigmatrix(i,1,1) = stress(i,1,1)
      sigmatrix(i,1,2) = stress(i,1,2)
      sigmatrix(i,4,4) = stress(i,1,4)
      sigmatrix(i,1,3) = sigmatrix(i,1,2)
      sigmatrix(i,2,1) = sigmatrix(i,1,2)
      sigmatrix(i,2,3) = sigmatrix(i,1,2)
      sigmatrix(i,3,1) = sigmatrix(i,1,2)
      sigmatrix(i,3,2) = sigmatrix(i,1,2)
      sigmatrix(i,2,2) = sigmatrix(i,1,1)
      sigmatrix(i,3,3) = sigmatrix(i,1,1)
      sigmatrix(i,5,5) = sigmatrix(i,4,4)
      sigmatrix(i,6,6) = sigmatrix(i,4,4)
!devt
      sigmdevt(i,1,1) = devt(i,1,1)
      sigmdevt(i,1,2) = devt(i,1,2)
      sigmdevt(i,4,4) = devt(i,1,4)
      sigmdevt(i,1,3) = sigmdevt(i,1,2)
      sigmdevt(i,2,1) = sigmdevt(i,1,2)
      sigmdevt(i,2,3) = sigmdevt(i,1,2)
      sigmdevt(i,3,1) = sigmdevt(i,1,2)
      sigmdevt(i,3,2) = sigmdevt(i,1,2)
      sigmdevt(i,2,2) = sigmdevt(i,1,1)
      sigmdevt(i,3,3) = sigmdevt(i,1,1)
      sigmdevt(i,5,5) = sigmdevt(i,4,4)
      sigmdevt(i,6,6) = sigmdevt(i,4,4)
!triclinic
    CASE(1)
      sigmatrix(i,1,1) = stress(i,1,1)
      sigmatrix(i,1,2) = (stress(i,1,2)+stress(i,2,1))/2
      sigmatrix(i,1,3) = (stress(i,1,3)+stress(i,3,1))/2
      sigmatrix(i,1,4) = (stress(i,1,5)+stress(i,4,1))/2
      sigmatrix(i,1,5) = (stress(i,1,6)+stress(i,5,1))/2
      sigmatrix(i,1,6) = (stress(i,1,4)+stress(i,6,1))/2
      sigmatrix(i,2,1) = sigmatrix(i,1,2)
      sigmatrix(i,2,2) = stress(i,2,2)
      sigmatrix(i,2,3) = (stress(i,2,3)+stress(i,3,2))/2
      sigmatrix(i,2,4) = (stress(i,2,5)+stress(i,4,2))/2
      sigmatrix(i,2,5) = (stress(i,2,6)+stress(i,5,2))/2
      sigmatrix(i,2,6) = (stress(i,2,4)+stress(i,6,2))/2
      sigmatrix(i,3,1) = sigmatrix(i,1,3)
      sigmatrix(i,3,2) = sigmatrix(i,2,3)
      sigmatrix(i,3,3) = stress(i,3,3)
      sigmatrix(i,3,4) = (stress(i,3,5)+stress(i,4,3))/2
      sigmatrix(i,3,5) = (stress(i,3,6)+stress(i,5,3))/2
      sigmatrix(i,3,6) = (stress(i,3,4)+stress(i,6,3))/2
      sigmatrix(i,4,1) = sigmatrix(i,1,4)
      sigmatrix(i,4,2) = sigmatrix(i,2,4)
      sigmatrix(i,4,3) = sigmatrix(i,3,4)
      sigmatrix(i,4,4) = stress(i,4,5)
      sigmatrix(i,4,5) = (stress(i,4,6)+stress(i,5,4))/2
      sigmatrix(i,4,6) = (stress(i,4,4)+stress(i,6,4))/2
      sigmatrix(i,5,1) = sigmatrix(i,1,5)
      sigmatrix(i,5,2) = sigmatrix(i,2,5)
      sigmatrix(i,5,3) = sigmatrix(i,3,5)
      sigmatrix(i,5,4) = sigmatrix(i,4,5)
      sigmatrix(i,5,5) = stress(i,5,6)
      sigmatrix(i,5,6) = (stress(i,5,4)+stress(i,6,6))/2
      sigmatrix(i,6,1) = sigmatrix(i,1,6)
      sigmatrix(i,6,2) = sigmatrix(i,2,6)
      sigmatrix(i,6,3) = sigmatrix(i,3,6)
      sigmatrix(i,6,4) = sigmatrix(i,4,6)
      sigmatrix(i,6,5) = sigmatrix(i,5,6)
      sigmatrix(i,6,6) = stress(i,6,4)
!devt
      sigmdevt(i,1,1) = devt(i,1,1)
      sigmdevt(i,1,2) = (devt(i,1,2)+devt(i,2,1))/2
      sigmdevt(i,1,3) = (devt(i,1,3)+devt(i,3,1))/2
      sigmdevt(i,1,4) = (devt(i,1,5)+devt(i,4,1))/2
      sigmdevt(i,1,5) = (devt(i,1,6)+devt(i,5,1))/2
      sigmdevt(i,1,6) = (devt(i,1,4)+devt(i,6,1))/2
      sigmdevt(i,2,1) = sigmdevt(i,1,2)
      sigmdevt(i,2,2) = devt(i,2,2)
      sigmdevt(i,2,3) = (devt(i,2,3)+devt(i,3,2))/2
      sigmdevt(i,2,4) = (devt(i,2,5)+devt(i,4,2))/2
      sigmdevt(i,2,5) = (devt(i,2,6)+devt(i,5,2))/2
      sigmdevt(i,2,6) = (devt(i,2,4)+devt(i,6,2))/2
      sigmdevt(i,3,1) = sigmdevt(i,1,3)
      sigmdevt(i,3,2) = sigmdevt(i,2,3)
      sigmdevt(i,3,3) = devt(i,3,3)
      sigmdevt(i,3,4) = (devt(i,3,5)+devt(i,4,3))/2
      sigmdevt(i,3,5) = (devt(i,3,6)+devt(i,5,3))/2
      sigmdevt(i,3,6) = (devt(i,3,4)+devt(i,6,3))/2
      sigmdevt(i,4,1) = sigmdevt(i,1,4)
      sigmdevt(i,4,2) = sigmdevt(i,2,4)
      sigmdevt(i,4,3) = sigmdevt(i,3,4)
      sigmdevt(i,4,4) = devt(i,4,4)
      sigmdevt(i,4,5) = (devt(i,4,6)+devt(i,5,4))/2
      sigmdevt(i,4,6) = (devt(i,4,4)+devt(i,6,4))/2
      sigmdevt(i,5,1) = sigmdevt(i,1,5)
      sigmdevt(i,5,2) = sigmdevt(i,2,5)
      sigmdevt(i,5,3) = sigmdevt(i,3,5)
      sigmdevt(i,5,4) = sigmdevt(i,4,5)
      sigmdevt(i,5,5) = devt(i,5,6)
      sigmdevt(i,5,6) = (devt(i,5,4)+devt(i,6,6))/2
      sigmdevt(i,6,1) = sigmdevt(i,1,6)
      sigmdevt(i,6,2) = sigmdevt(i,2,6)
      sigmdevt(i,6,3) = sigmdevt(i,3,6)
      sigmdevt(i,6,4) = sigmdevt(i,4,6)
      sigmdevt(i,6,5) = sigmdevt(i,5,6)
      sigmdevt(i,6,6) = devt(i,6,4)
!monoclinic
    CASE(2)
      sigmatrix(i,1,1) = stress(i,1,1)
      sigmatrix(i,1,2) = (stress(i,1,2)+stress(i,2,1))/2
      sigmatrix(i,1,3) = (stress(i,1,3)+stress(i,4,1))/2
      sigmatrix(i,1,5) = (stress(i,3,1)+stress(i,1,6))/2
      sigmatrix(i,2,2) = stress(i,2,2)
      sigmatrix(i,2,3) = (stress(i,2,3)+stress(i,4,2))/2
      sigmatrix(i,2,5) = (stress(i,3,2)+stress(i,2,6))/2 
      sigmatrix(i,3,3) = stress(i,4,3)
      sigmatrix(i,3,5) = (stress(i,3,3)+stress(i,4,6))/2
      sigmatrix(i,4,4) = stress(i,1,5)
      sigmatrix(i,4,6) = (stress(i,1,4)+stress(i,2,5))/2
      sigmatrix(i,5,5) = stress(i,3,6)
      sigmatrix(i,6,6) = stress(i,2,4)
      sigmatrix(i,2,1) = sigmatrix(i,1,2)
      sigmatrix(i,3,1) = sigmatrix(i,1,3)
      sigmatrix(i,3,2) = sigmatrix(i,2,3)
      sigmatrix(i,6,4) = sigmatrix(i,4,6)
      sigmatrix(i,5,1) = sigmatrix(i,1,5)
      sigmatrix(i,5,2) = sigmatrix(i,2,5)
      sigmatrix(i,5,3) = sigmatrix(i,3,5)
!devt
      sigmdevt(i,1,1) = devt(i,1,1)
      sigmdevt(i,1,2) = (devt(i,1,2)+devt(i,2,1))/2
      sigmdevt(i,1,3) = (devt(i,1,3)+devt(i,4,1))/2
      sigmdevt(i,1,5) = (devt(i,3,1)+devt(i,1,6))/2
      sigmdevt(i,2,2) = devt(i,2,2)
      sigmdevt(i,2,3) = (devt(i,2,3)+devt(i,4,2))/2
      sigmdevt(i,2,5) = (devt(i,3,2)+devt(i,2,6))/2
      sigmdevt(i,3,3) = devt(i,4,3)
      sigmdevt(i,3,5) = (devt(i,3,3)+devt(i,4,6))/2
      sigmdevt(i,4,4) = devt(i,1,5)
      sigmdevt(i,4,6) = (devt(i,1,4)+devt(i,2,5))/2
      sigmdevt(i,5,5) = devt(i,3,6)
      sigmdevt(i,6,6) = devt(i,2,4)
      sigmdevt(i,2,1) = sigmdevt(i,1,2)
      sigmdevt(i,3,1) = sigmdevt(i,1,3)
      sigmdevt(i,3,2) = sigmdevt(i,2,3)
      sigmdevt(i,6,4) = sigmdevt(i,4,6)
      sigmdevt(i,5,1) = sigmdevt(i,1,5)
      sigmdevt(i,5,2) = sigmdevt(i,2,5)
      sigmdevt(i,5,3) = sigmdevt(i,3,5)
!orthorhombic
    CASE(3)
      sigmatrix(i,1,1) = stress(i,1,1) - stress(i,2,1)
      sigmatrix(i,1,2) = stress(i,1,1) + stress(i,2,1)
      sigmatrix(i,1,3) = stress(i,1,3) - stress(i,2,3)
      sigmatrix(i,2,2) = stress(i,1,2) + stress(i,2,2)
      sigmatrix(i,2,3) = stress(i,1,3) + stress(i,2,3)
      sigmatrix(i,3,3) = stress(i,3,3)
      sigmatrix(i,4,4) = stress(i,1,5)
      sigmatrix(i,5,5) = stress(i,2,6)
      sigmatrix(i,6,6) = stress(i,3,4)
      sigmatrix(i,2,1) = sigmatrix(i,1,2)
      sigmatrix(i,3,1) = sigmatrix(i,1,3)
      sigmatrix(i,3,2) = sigmatrix(i,2,3)
!devt
      sigmdevt(i,1,1) = devt(i,1,1) + devt(i,2,1)
      sigmdevt(i,1,2) = devt(i,1,1) + devt(i,2,1)
      sigmdevt(i,1,3) = devt(i,1,3) + devt(i,2,3)
      sigmdevt(i,2,2) = devt(i,1,2) + devt(i,2,2)
      sigmdevt(i,2,3) = devt(i,1,3) + devt(i,2,3)
      sigmdevt(i,3,3) = devt(i,3,3)
      sigmdevt(i,4,4) = devt(i,1,5)
      sigmdevt(i,5,5) = devt(i,2,6)
      sigmdevt(i,6,6) = devt(i,3,4)
      sigmdevt(i,2,1) = sigmdevt(i,1,2)
      sigmdevt(i,3,1) = sigmdevt(i,1,3)
      sigmdevt(i,3,2) = sigmdevt(i,2,3)
!tetragonal
    CASE(4)
      sigmatrix(i,1,3) = (stress(i,1,3)-stress(i,2,3)) + (stress(i,2,1) + &
                          stress(i,2,2) - (stress(i,1,1)+stress(i,1,2))/2)/2
      sigmatrix(i,4,4) = stress(i,2,5)
      sigmatrix(i,6,6) = stress(i,1,4)
      sigmatrix(i,3,3) = 2 * stress(i,2,3) - stress(i,1,3)
      sigmatrix(i,1,1) = 2 * stress(i,2,2) - 2 * sigmatrix(i,1,3)
      sigmatrix(i,1,2) = 2 * stress(i,2,1) - 2 * sigmatrix(i,1,3)
      sigmatrix(i,1,6) = (stress(i,1,1) - stress(i,1,2) - sigmatrix(i,1,1) & 
                         + sigmatrix(i,1,2))/2
      sigmatrix(i,2,6) = 0 - sigmatrix(i,1,6) 
      sigmatrix(i,2,1) = sigmatrix(i,1,2) 
      sigmatrix(i,2,2) = sigmatrix(i,1,1)
      sigmatrix(i,2,3) = sigmatrix(i,1,3) 
      sigmatrix(i,3,1) = sigmatrix(i,1,3) 
      sigmatrix(i,3,2) = sigmatrix(i,1,3) 
      sigmatrix(i,5,5) = sigmatrix(i,4,4) 
      sigmatrix(i,6,1) = sigmatrix(i,1,6)
      sigmatrix(i,6,2) = 0 - sigmatrix(i,1,6)
!devt
      sigmdevt(i,1,3) = devt(i,2,2) + devt(i,2,1) + devt(i,1,3) &
                         + devt(i,2,3) + (devt(i,1,1)+devt(i,1,2))/2
      sigmdevt(i,4,4) = devt(i,2,5)
      sigmdevt(i,6,6) = devt(i,1,4)
      sigmdevt(i,3,3) = 2 * devt(i,2,3) + devt(i,1,3)
      sigmdevt(i,1,1) = 2 * devt(i,2,2) + 2 * sigmdevt(i,1,3)
      sigmdevt(i,1,2) = 2 * devt(i,2,1) + 2 * sigmdevt(i,1,3)
      sigmdevt(i,1,6) = (devt(i,1,1) + devt(i,1,2) + sigmdevt(i,1,1) &
                         + sigmdevt(i,1,2))/2
      sigmdevt(i,2,6) = 0 + sigmdevt(i,1,6)
      sigmdevt(i,2,1) = sigmdevt(i,1,2)
      sigmdevt(i,2,3) = sigmdevt(i,1,3)
      sigmdevt(i,3,1) = sigmdevt(i,1,3)
      sigmdevt(i,3,2) = sigmdevt(i,1,3)
      sigmdevt(i,5,5) = sigmdevt(i,4,4)
!triagonal
    CASE(5)
      sigmatrix(i,1,1) = (stress(i,2,4) + stress(i,2,6) + stress(i,1,1) &
                         + stress(i,1,2)*3/2)/2
      sigmatrix(i,1,2) = stress(i,1,1) - stress(i,2,6)
      sigmatrix(i,1,3) = (stress(i,1,3) + stress(i,2,1) + stress(i,2,2))/3
      sigmatrix(i,1,4) = stress(i,2,6)
      sigmatrix(i,1,5) = 0 - (stress(i,1,6) + stress(i,1,4) + &
                            stress(i,2,5))/3 
      sigmatrix(i,3,3) = stress(i,2,3)
      sigmatrix(i,4,4) = stress(i,1,5) + stress(i,2,6) 
      sigmatrix(i,2,1) = sigmatrix(i,1,2) 
      sigmatrix(i,2,2) = sigmatrix(i,1,1) 
      sigmatrix(i,2,3) = sigmatrix(i,1,3)
      sigmatrix(i,3,1) = sigmatrix(i,1,3)
      sigmatrix(i,3,2) = sigmatrix(i,1,3)
      sigmatrix(i,2,4) = 0 - sigmatrix(i,1,4)
      sigmatrix(i,4,1) = sigmatrix(i,1,4)
      sigmatrix(i,4,2) = 0 - sigmatrix(i,1,4)
      sigmatrix(i,5,6) = sigmatrix(i,1,4)
      sigmatrix(i,6,5) = sigmatrix(i,1,4)
      sigmatrix(i,2,5) = 0 - sigmatrix(i,1,5)
      sigmatrix(i,4,6) = 0 - sigmatrix(i,1,5)
      sigmatrix(i,5,2) = 0 - sigmatrix(i,1,5)
      sigmatrix(i,6,4) = 0 - sigmatrix(i,1,5)
      sigmatrix(i,5,1) = sigmatrix(i,1,5)
      sigmatrix(i,5,5) = sigmatrix(i,4,4)
      sigmatrix(i,6,6) = (sigmatrix(i,1,1) - sigmatrix(i,1,2))/2
!devt
      sigmdevt(i,1,1) = (devt(i,2,4) + devt(i,2,6) + devt(i,1,1) &
                         + devt(i,1,2)*3/2)/2
      sigmdevt(i,1,2) = devt(i,1,1) + devt(i,2,6)
      sigmdevt(i,1,3) = (devt(i,1,3) + devt(i,2,1) + devt(i,2,2))/3
      sigmdevt(i,1,4) = devt(i,2,6)
      sigmdevt(i,1,5) = 0 + (devt(i,1,6) + devt(i,1,4) + & 
                            devt(i,2,5))/3
      sigmdevt(i,3,3) = devt(i,2,3)
      sigmdevt(i,4,4) = devt(i,1,5) + devt(i,2,6)
      sigmdevt(i,2,1) = sigmdevt(i,1,2)
      sigmdevt(i,2,2) = sigmdevt(i,1,1)
      sigmdevt(i,2,3) = sigmdevt(i,1,3)
      sigmdevt(i,3,1) = sigmdevt(i,1,3)
      sigmdevt(i,3,2) = sigmdevt(i,1,3)
      sigmdevt(i,2,4) = 0 + sigmdevt(i,1,4)
      sigmdevt(i,4,1) = sigmdevt(i,1,4)
      sigmdevt(i,4,2) = 0 + sigmdevt(i,1,4)
      sigmdevt(i,5,6) = sigmdevt(i,1,4)
      sigmdevt(i,6,5) = sigmdevt(i,1,4)
      sigmdevt(i,2,5) = 0 + sigmdevt(i,1,5)
      sigmdevt(i,4,6) = 0 + sigmdevt(i,1,5)
      sigmdevt(i,5,2) = 0 + sigmdevt(i,1,5)
      sigmdevt(i,6,4) = 0 + sigmdevt(i,1,5)
      sigmdevt(i,5,1) = sigmdevt(i,1,5)
      sigmdevt(i,5,5) = sigmdevt(i,4,4)
      sigmdevt(i,6,6) = (sigmdevt(i,1,1) + sigmdevt(i,1,2))/2
  END SELECT
END DO

!do i=1, 4
!  do j=1, 6
!    write(*,*) i, j, sigmatrix(i,j,:)
!  end do
!end do

RETURN

END SUBROUTINE SIGMACALC

!------------------------------------------------------------------------------
SUBROUTINE ELASCALC(latt,cij,rho,born,vrh)

INTEGER                  :: i,j,latt
REAL(DP), DIMENSION(6,6) :: cij, sij
REAL(DP)                 :: rho
LOGICAL,  INTENT(OUT)    :: born
REAL(DP), INTENT(OUT)    :: vrh(3,7)

! elastic compliance constants
sij = AINV(cij)
WRITE(*,*) "Elastic Compliance Constant Matrix : "
DO i=1,6
  WRITE(*,'(6F16.8)') sij(i,:)
END DO

! bulk modulus
vrh(1,1) = ( cij(1,1)+cij(2,2)+cij(3,3) + 2*(cij(1,2)+cij(2,3)+cij(3,1)) ) / 9
vrh(2,1) = 1 / ( sij(1,1)+sij(2,2)+sij(3,3) + 2*(sij(1,2)+sij(2,3)+sij(3,1)) )
vrh(3,1) = ( vrh(1,1) + vrh(2,1) ) / 2
!shear modulus
vrh(1,2) = ( cij(1,1)+cij(2,2)+cij(3,3) - (cij(1,2)+cij(2,3)+cij(3,1)) + &
             3*(cij(4,4)+cij(5,5)+cij(6,6)) ) / 15
vrh(2,2) = 15 / ( 4*(sij(1,1)+sij(2,2)+sij(3,3)) - 4*(sij(1,2)+sij(2,3)+ &
            sij(3,1)) + 3*(sij(4,4)+sij(5,5)+sij(6,6)) )
vrh(3,2) = ( vrh(1,2) + vrh(2,2) ) / 2
! Young's modulus
vrh(:,3) = 9 * vrh(:,2) * vrh(:,1) / (vrh(:,2)+3*vrh(:,1) )
! Poisson's ratio
vrh(:,4) = ( 3*vrh(:,1) - 2*vrh(:,2) ) /( 3*vrh(:,1) + vrh(:,2) ) / 2
! shear velocity
vrh(:,5) = SQRT(vrh(:,2)/rho)
! compressional velocity
vrh(:,6) = SQRT( (vrh(:,1)+vrh(:,2)*4/3) / rho )
! Debye velocity
vrh(:,7) = 1/(vrh(:,6)**3) + 2/(vrh(:,5)**3) 
vrh(:,7) = 3/vrh(:,7)
vrh(:,7) = vrh(:,7)**0.333333333  ! **(1/3) will not work since fortran will do **(1 and /3 
! check Born mechanical stability. M. Born and K. Huang, Dynamics Theory of Crystal Lattices (Oxford
! University Press, 1954).
SELECT CASE(latt)
  CASE(1,2) !triclinic and monoclinic
    DO i=1,6
      DO j=1,6
        IF(cij(i,j)**2<cij(i,i)*cij(j,j) .and. cij(i,j)>0) THEN
          born= .TRUE.
         ELSE 
           born=.FALSE.
         END IF
      END DO
    END DO
  CASE(3) !orthorhombic 
    IF( cij(1,1)>0 .and. cij(1,2)**2 < cij(1,1)*cij(2,2) .and. cij(4,4)>0 &
      .and. cij(5,5)>0 .and. cij(6,6)>0 .and. (cij(1,1)*cij(2,2)*cij(3,3) + &
      2*cij(1,2)*cij(1,3)*cij(2,3))>(cij(1,1)*cij(2,3)**2 + & 
      cij(2,2)*cij(1,3)**2 + cij(3,3)*cij(1,2)**2) .and. cij(4,4)>0 .and. &
      cij(6,6)>0 ) THEN
      born= .TRUE.
    ELSE 
      born=.FALSE.
    END IF
  CASE(4) ! tetrahedral
    IF((cij(1,1)>ABS(cij(1,2))) .and. ( 2*cij(1,3)**2 < (cij(3,3)*(cij(1,1)+cij(1,2))) ) & 
        .and. cij(4,4)>0 .and. (2*cij(1,6)**2)<(cij(6,6)*(cij(1,1)-cij(1,2))) ) THEN
      born= .TRUE.
    ELSE 
      born=.FALSE.
    END IF
  CASE(5) ! triagonal 3, -3
    IF((cij(1,1)>ABS(cij(1,2))) .and. ( 2*cij(1,3)**2 < (cij(3,3)*(cij(1,1)+cij(1,2))) ) &
        .and. cij(4,4)>0 .and. (cij(1,4)**2+cij(1,5)**2) < cij(4,4)*cij(6,6) ) THEN
      born= .TRUE.
    ELSE 
      born=.FALSE.
    END IF
  CASE(6) !hexagonal
    IF( cij(1,1)>ABS(cij(1,2)) .and. 2*cij(1,3)**2 < cij(3,3)*(cij(1,1)+cij(1,2)) &
        .and. cij(4,4)>0 .and. cij(6,6)>0 ) THEN
      born= .TRUE.
    ELSE 
      born=.FALSE.
    END IF
  CASE(7) !cubic
    IF((cij(1,1)-cij(1,2))>0 .and. (cij(1,1)+2*cij(1,2))>0 .and. cij(4,4)>0) THEN
      born= .TRUE.
    ELSE 
      born=.FALSE.
    END IF
END SELECT

RETURN

CONTAINS
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


END SUBROUTINE ELASCALC



 ! get the element name
SUBROUTINE getelemtnm(line,elemnm,ntype)
CHARACTER(len=255),intent(IN):: line
CHARACTER(len=255):: line2,line3
CHARACTER(len=2),dimension(:),ALLOCATABLE, intent(OUT) :: elemnm
CHARACTER(len=2),dimension(:),ALLOCATABLE :: elemt
INTEGER, intent(OUT) :: ntype 
INTEGER :: ios, i, pos

line2=line
ios=0
i = 0
ALLOCATE(elemt(20))

DO WHILE (ios==0)
  i = i + 1
  READ(line2,*,IOSTAT=ios) elemt(i)
  IF (ios==0) then
    elemt(i) = ADJUSTL(elemt(i))
    pos = Index(line2,trim(elemt(i))) + len_trim(elemt(i))
    line2=""
    line2(pos:255)=line(pos:255)
    IF (len_trim(line2)<=0) then
      ios=10
    END IF
  END IF
END DO
ntype = i 
ALLOCATE(elemnm(ntype))
DO WHILE (i>0)
   elemnm(i) = elemt(i)
   i = i - 1
END DO
END SUBROUTINE getelemtnm

SUBROUTINE getatmnr(line,maxatom,atnr)
CHARACTER(len=255),intent(IN):: line
INTEGER,intent(OUT) :: maxatom
INTEGER :: j, ios,nratoms,atnrpos
INTEGER,DIMENSION(:),ALLOCATABLE,intent(OUT) :: atnr
INTEGER,DIMENSION(:),ALLOCATABLE :: atnrtmp
CHARACTER(len=255):: line2,atl

ios=0
nratoms=0
line2=line
j = 0
ALLOCATE(atnrtmp(20))

DO WHILE (ios==0)
   j = j + 1
   read(line2,*,IOSTAT=ios) atnrtmp(j)
   IF (ios==0) then
      nratoms=nratoms+atnrtmp(j)
      write(atl,*) atnrtmp(j)
      atl=ADJUSTL(atl)
      atnrpos=Index(line2,trim(atl))+len_trim(atl)
! Index finds the location of a substring in another string,returns 0 if not found
! len_trim counts the number of characters including blanks
! this sentence locates the point in the line
      line2=""
      line2(atnrpos:255)=line(atnrpos:255)
      IF (len_trim(line2)<=0) then
        ios=10
      END IF
   END IF
END DO
maxatom=nratoms
ALLOCATE(atnr(j))
DO WHILE (j>0)
  atnr(j) = atnrtmp(j)   
  j = j - 1
END DO

RETURN

END SUBROUTINE getatmnr

END MODULE CALC
