MODULE nrtype
!------------------------------------------------------------------------------
! Data type definitions and constant settings
!------------------------------------------------------------------------------

  INTEGER, PARAMETER :: I4B = selected_int_kind(9) 
  INTEGER, PARAMETER :: I2B = selected_int_kind(4) 
  INTEGER, PARAMETER :: I1B = selected_int_kind(2) 

  INTEGER, PARAMETER :: SP = kind(1.0) 
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(10) 

  INTEGER, PARAMETER :: SPC = kind( (1.0,1.0) ) 
  INTEGER, PARAMETER :: DPC = kind( (1.0_dp,1.0_dp) ) 

  INTEGER, PARAMETER :: LGT = kind( .true. ) 

  REAL(DP), PARAMETER :: PI = 3.141592653589793238462643383279502884197_dp
  REAL(DP), PARAMETER :: TWOPI = 2._dp*PI

  COMPLEX(DPC), PARAMETER :: im = (0._dp,1._dp)


  !  Some important Parameters, to convert to a.u.
  !  - AUTOA  = 1. a.u. in Angstroem
  !  - RYTOEV = 1 Ry in Ev
  !  - EVTOJ  = 1 eV in Joule
  !  - AMTOKG = 1 atomic mass unit ("proton mass") in kg
  !  - BOLKEV = Boltzmanns constant in eV/K
  !  - BOLK   = Boltzmanns constant in Joule/K

  REAL(DP), PARAMETER :: EVTOJ=1.60217733E-19_dp,AMTOKG=1.6605402E-27_dp, &
       BOLKEV=8.6173857E-5_dp,BOLK=BOLKEV*EVTOJ, HPLANK=6.6262E-34_dp, &
       zero = 1.d-6, convert_thz_to_cm = 33.357, &
       convert_thz_to_meV = 4.1357

  CHARACTER(len=2), DIMENSION(109) :: NMARRAY=(/ 'H ','He','Li','Be','B ',&
  'C ','N ','O ','F ','Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar','K ',&
  'Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As',&
  'Se','Br','Kr','Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag',&
  'Cd','In','Sn','Sb','Te','I ','Xe','Cs','Be','La','Ce','Pr','Nd','Pm',&
  'Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re',&
  'Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',&
  'Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',&
  'Rf','Db','Sg','Bh','Hs','Mt' /)

  REAL, DIMENSION(109)    :: MSARRAY=(/ 1.0079,4.0026,6.941,9.0122,10.811,&
  12.0107,14.0067,15.9994,18.9984,20.1797,22.9897,24.305,26.9815,28.0855,&
  30.9738,32.065,35.453,39.948,39.0983,40.078,44.9559,47.867,50.9415,&
  51.9961,54.938,55.845,58.9332,58.6934,63.546,65.39,69.723,72.64,74.9216,&
  78.96,79.904,83.8,85.4678,87.62,88.9059,91.224,92.9064,95.94,98.0,101.07,& 
  102.9055,106.42,107.8682,112.411,114.818,118.71,121.76,127.6,126.9045,&
  131.293,132.9055,137.327,138.9055,140.116,140.9077,144.24,145.0,150.36,&
  151.964,157.25,158.9253,162.5,164.9303,167.259,168.9342,173.04,174.967,&
  178.49,180.9479,183.84,186.207,190.23,192.217,195.078,196.9665,200.59,&
  204.3833,207.2,208.9804,209.0,210.0,222.0,223.0,226.0,227.0,232.0381,&
  231.0359,238.0289,237.0,244.0,243.0,247.0,247.0,251.0,252.0,257.0,258.0,&
  259.0,262.0,261.0,262.0,266.0,264.0,277.0,268.0 /)


END MODULE nrtype

