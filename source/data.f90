!------------------------------------------------------------------------------
!Definitions of the main variables used by elasticity
!------------------------------------------------------------------------------

MODULE DATA

USE NRTYPE

IMPLICIT NONE

INTEGER :: i,j,latt,nstrain,numdist,polynmord,lsys

REAL(DP) :: dummyr, rho
REAL(DP), dimension(:), allocatable :: delta, x,y,polycoef
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: stress, sigmatrix,devt,sigmdevt
REAL(DP), DIMENSION(6,6) :: cijmatrix
REAL(DP), DIMENSION(3,7) :: vrh

CHARACTER(len=20) :: fnam, dirnam
CHARACTER(len=80) :: cmd,line
CHARACTER(len=8) :: code,md,npt

LOGICAL  :: lopen, found, born


END MODULE DATA

