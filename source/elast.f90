!--------------------------------------------------------------------------
! Program: ElasT.1.2
! Copyright (C) Yunguo Li @ USTC                         
! Email :: liyunguo@ustc.edu.cn           
! Function :: calculate elastic properties via the strain-strss method
!             by using VASP (QE, GULP in development)
!
! Input ::  elastic.in, POSCAR, VASP input files
! Output::  POSCAR_s*, stress.out, devt.dat
!--------------------------------------------------------------------------

PROGRAM  ELAST

USE DATA
USE NRTYPE
USE READIN
USE CALC
USE FITTING

IMPLICIT  NONE

CHARACTER(80), PARAMETER :: version='            ElasT, VERSION 1.2, &
                                     Copyright (C) 2021, Y. Li, USTC'
CHARACTER(80), PARAMETER :: citation='        Y. Li, et al., Computer &
                                     Physics Communications, 2022, 108280'

CHARACTER(80), PARAMETER :: license='    This program comes with ABSOLUTELY &
                                     NO WARRANTY, see LICENSE.'
!WRITE(*,*)'----------------------------------------------------------------'
WRITE(*,*)''
WRITE(*,*)'                        ____|   |                 __ __|        '
WRITE(*,*)'                        __|     |    _` |    __|     |          '
WRITE(*,*)'                        |       |   (   |  \__ \     |          '
WRITE(*,*)'                       _____|  _|  \__,_|  ____/    _|          '
WRITE(*,*)''
WRITE(*,*) version
!WRITE(*,*) citation
WRITE(*,*) license
WRITE(*,*)'----------------------------------------------------------------'


!------------------------------------------------------------------------------
! read lattice type, numbser of strain, strains, calc. code, 
! type of calc., polynominal order
!------------------------------------------------------------------------------

CALL READER(latt, nstrain, delta, code, md, npt,polynmord)

SELECT CASE(latt)
  CASE(1)
    WRITE(*,*) "Working System            : Triclinic"
    numdist = 6
  CASE(2)
    WRITE(*,*) "Working System            : Monoclinic"
    numdist = 4
  CASE(3)
    WRITE(*,*) "Working System            : Orthorhombic"
    numdist = 3
  CASE(4)
    WRITE(*,*) "Working System            : Tetragonal"
    numdist = 2
  CASE(5)
    WRITE(*,*) "Working System            : Trigonal(3,-3)"
    numdist = 2
!  CASE(6)
!    WRITE(*,*) "Working System            : Trigonal(32,3m,-32/m)"
!    numdist = 2
  CASE(6)
    WRITE(*,*) "Working System            : Hexagonal"
    numdist = 2
  CASE(7)
    WRITE(*,*) "Working System            : Cubic"
    numdist = 1
END SELECT 


WRITE(*,'(A29,A8)') " Calculations Using Code   : ", ADJUSTL(code)
IF(trim(adjustl(md))=='no' .or. trim(adjustl(md))=='No' .or. &
   trim(adjustl(md))=='NO' .or. trim(adjustl(md))=='N') THEN
  WRITE(*,'(A29,A6)') " Type of Calculations      : ", "Static"
ELSE
  WRITE(*,'(A29,A2)') " Type of Calculations      : ", "MD"
END IF

IF(trim(adjustl(npt))=='yes' .or. trim(adjustl(npt))=='Yes' .or. &
    trim(adjustl(npt))=='YES' .or. trim(adjustl(npt))=='Y') THEN
  WRITE(*,*)'----------------------------------------------------------------'
  CALL SSNPT(latt,cijmatrix)
  GOTO 130
END IF

WRITE(*,'(A29,I4)') " Number of Distortions     :  ", numdist
WRITE(*,'(A29,I4)') " Number of Applied Strains : ", nstrain
WRITE(*,'(A29)',advance='no') " Applied Strains are       : "

DO j = 1, nstrain
   WRITE(*,'(F6.3,A3)', advance='no') delta(j),'   '
   IF(j==nstrain) write(*,*)
END DO
WRITE(*,*)'----------------------------------------------------------------'

!------------------------------------------------------------------------------
! check existences of elastic.out,stress.out, calc. directories
!------------------------------------------------------------------------------

INQUIRE(file="elastic.out",EXIST=found)
IF(found) THEN
  WRITE(*,*) "Analysis already done! If else please remove elastic.out." 
  STOP
END IF 

INQUIRE(file="stress.out",EXIST=found)
IF(found) THEN
  WRITE(*,*) "Stresses already collected! If else please remove &
              &stress.out and devt.dat."
  GOTO 120
END IF 

! structure and output files  will be checked
DO i = 0, nstrain
  DO j = 1, numdist
    WRITE(fnam,'(A6)') "POSCAR"
    WRITE(dirnam,'(A1,I1,A1,I1)') "s",i,"d",j
    IF(i==0)  WRITE(dirnam,'(A1,I1,A1,I1)') "s",i,"d",j-1
    cmd=trim(dirnam)//"/"//trim(fnam)
    INQUIRE(file=trim(cmd),EXIST=found)
    cmd=trim(dirnam)//"/"//"OUTCAR"
    INQUIRE(file=trim(cmd),EXIST=born)
    IF(i==0) EXIT
  END DO
  IF(.not.found) EXIT
  IF(.not.born) EXIT
END DO

IF(found) THEN
  WRITE(*,*) "Calculation directories already created."
  IF(born) THEN
    WRITE(*,*) "Calculations found completed." 
    GOTO 110
  ELSE 
    WRITE(*,*) "Calculations found incompleted, Please perform calculations."
    STOP
  END IF
END IF

!-------------------------------------------------------------------------------
! create distorted structures
!-------------------------------------------------------------------------------

CALL DIST(latt,nstrain,numdist,delta)

!-------------------------------------------------------------------------------
! set up calculations
!-------------------------------------------------------------------------------

DO i = 0, nstrain
  DO j = 1, numdist
! 'if' construct to treat s0d0 dir.
    IF(i==0) THEN
      WRITE(dirnam,'(A1,I1,A1,I1)') "s",i,"d",j-1
      WRITE(fnam,'(A6)') "POSCAR"
    ELSE
      WRITE(fnam,'(A8,I1,A1,I1)') "POSCAR_s",i,"d",j
      WRITE(dirnam,'(A1,I1,A1,I1)') "s",i,"d",j
    END IF
    cmd="mkdir "//trim(dirnam)
    CALL SYSTEM(cmd, lsys)
    IF(lsys/=0) WRITE(*,*) "Creating directories failed"
    cmd="cp "//fnam//" "//trim(dirnam)//"/POSCAR" 
    CALL SYSTEM(cmd, lsys)
    IF(lsys/=0) WRITE(*,*) "Copying POSCAR_* into directories failed"
    cmd="cp INCAR POTCAR KPOINTS "//trim(dirnam)//"/"
    CALL SYSTEM(cmd, lsys)
    IF(lsys/=0) WRITE(*,*) "Copying POTCAR INCAR KPOINTS into directories & 
                            &failed" 
    IF(i==0) EXIT
  END DO
END DO
STOP

!-------------------------------------------------------------------------------
! extract stresses
!-------------------------------------------------------------------------------

!extract stress
110 WRITE(*,*) "Extracting ..."

ALLOCATE(stress(0:nstrain,1:numdist,6))
ALLOCATE(devt(0:nstrain,1:numdist,6))
CALL EXTRACT(code,md,nstrain,numdist,stress,devt)


!check consistence of number of stresses and number of strains
!if inconsistent, stop progress; read stress()
120 WRITE(*,'(A41)') " Checking consistence of stresses ... ..." 

IF(.not. ALLOCATED(stress)) ALLOCATE(stress(0:nstrain,1:numdist,6))
IF(.not. ALLOCATED(devt)) ALLOCATE(devt(0:nstrain,1:numdist,6))
OPEN(unit=62,file='stress.out',status='old')
OPEN(unit=26,file='devt.dat',status='old')
DO i=0, nstrain
  DO j=1, numdist
    READ(62,*,IOSTAT=lsys) dummyr,dummyr,stress(i,j,:)
    write(*,'(6F16.8)') stress(i,j,:)
    READ(26,*) dummyr,dummyr,devt(i,j,:)
    IF(i==0) EXIT
  END DO
END DO  
CLOSE(62)
CLOSE(26)

IF(lsys/=0) THEN 
  WRITE(*,'(A55)') " Collected No. of stresses incorrect, chech stress.out."
  STOP
ELSE
  WRITE(*,'(A25)') " No. of stresses correct!"  
END IF

!-------------------------------------------------------------------------------
! calculate and fit cij
!-------------------------------------------------------------------------------
! calc. sigma stress for each matrix element
ALLOCATE(sigmatrix(nstrain,6,6))
ALLOCATE(sigmdevt(nstrain,6,6))
CALL SIGMACALC(latt,sigmatrix,stress,devt,sigmdevt)

! convert sigma unit from kB to GPa
sigmatrix = sigmatrix / 10 ; sigmdevt = sigmdevt / 10

! fitting cij 
ALLOCATE(x(nstrain)) ; ALLOCATE(y(nstrain)) ; ALLOCATE(polycoef(polynmord+1)) 

IF(polynmord>(nstrain-1)) WRITE(*,*) "Warning : fitting order bigger than No. of &
  data, overfitting may be expected."

DO i=1,6
  DO j=1,6
    x=delta
    y=sigmatrix(1:nstrain,i,j) 
    polycoef = POLYFIT(x,y,polynmord) 
    cijmatrix(i,j) = polycoef(2) 
    END DO
END DO

! fitting deviation

DO i=1,6
  DO j=1,6
    x=delta
    y=sigmdevt(1:nstrain,i,j)
    polycoef = POLYFIT(x,y,polynmord)
    sigmdevt(:,i,j) = polycoef(2)
    END DO
END DO

! write cij, dev.
WRITE(*,*) "Elastic Stiffness Constant Matrix : "
DO i=1,6
  WRITE(*,'(6F16.8)') cijmatrix(i,:)
END DO

WRITE(*,*) "Deviation of Elastic Stiffness Constant : "
DO i=1,6
  WRITE(*,'(6F16.8)') sum(abs(sigmdevt(:,i,1)))*10, sum(abs(sigmdevt(:,i,2)))*10, &
             sum(abs(sigmdevt(:,i,3)))*10, sum(abs(sigmdevt(:,i,4)))*10, &
            sum(abs(sigmdevt(:,i,5)))*10, sum(abs(sigmdevt(:,i,6)))*10
END DO

DEALLOCATE(stress)    ; DEALLOCATE(devt)
DEALLOCATE(sigmatrix) ; DEALLOCATE(sigmdevt)
DEALLOCATE(x) ; DEALLOCATE(y) ; DEALLOCATE(polycoef)

!-------------------------------------------------------------------------------
! calc. elastic properties
!-------------------------------------------------------------------------------

130 vrh = 0.0
WRITE(fnam,'(A6)') "POSCAR" ; CALL RHOCALC(fnam,rho)
CALL ELASCALC(latt,cijmatrix,rho,born,vrh)
WRITE(*,*) "Born Mechanical Stability : ", born
WRITE(*,'(A28,F13.8)') " Density (g/cm^3)          :",rho
WRITE(*,*) "Approx.  Bulk Modul(GPa)  Shear Modul(GPa)   Young's &
            &Modul(GPa)   Poisson's Ratio      Vs(km/s)        Vp(km/s)           &      
            &Vd(km/s)"
WRITE(*,'(A10,7(F12.4,6X))') " Voigt :  ",vrh(1,:)
WRITE(*,'(A10,7(F12.4,6X))') " Reuss :  ",vrh(2,:) 
WRITE(*,'(A10,7(F12.4,6X))') " Hill  :  ",vrh(3,:)

!-------------------------------------------------------------------------------
! end of prog
!-------------------------------------------------------------------------------

STOP

END PROGRAM ELAST
