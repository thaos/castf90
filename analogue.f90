!> © LSCE – Laboratory related to CEA/DRF – CNRS – UVSQ, 
!! Sabine Radanovics (sabine.radanovics@lsce.ipsl.fr) andPascal Yiou (pascal.yiou@lsce.ipsl.fr)
!! This source code is part of the CASTf90 software IDDN.FR.001.030008.000.S.P.2016.000.20700
!!
!! This software is governed by the CeCILL license under French law and abiding by the rules of distribution 
!! of free software. You can use, modify and / or redistribute the software under the terms of the 
!! CeCILL license as circulated by CEA, CNRS and INRIA at the following URL "http://www.cecill.info".
!!
!! As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by 
!! the license, users are provided only with a limited warranty and the software's author, 
!! the holder of the economic rights, and the successive licensors have only limited liability.
!!
!! In this respect, the user's attention is drawn to the risks associated with loading, using, 
!! modifying and/or developing or reproducing the software by the user in light of its specific status 
!! of free software, that may mean that it is complicated to manipulate, and that also therefore means 
!! that it is reserved for developers and experienced professionals having in-depth computer knowledge. 
!! Users are therefore encouraged to load and test the software's suitability as regards their requirements 
!! in conditions enabling the security of their systems and/or data to be ensured and, more generally, 
!! to use and operate it in the same conditions as regards security.
!!
!! The fact that you are presently reading this means that you have had knowledge of the CeCILL license 
!! and that you accept its terms.
!!
PROGRAM analogue
!! Program uses IBM extension SYSTEM()
USE write_file
USE routines

IMPLICIT NONE
TYPE (config_type) :: configs
TYPE (dims_type) :: dim_archi
TYPE (dims_type) :: dim_archi_tmp
TYPE (dims_type) :: dim_sim
TYPE (dims_type) :: dim_sim_tmp
CHARACTER(8) :: clockdate
CHARACTER(10) :: clocktime
INTEGER, ALLOCATABLE :: analogue_dates(:,:)
REAL(8), ALLOCATABLE :: distances(:,:)
REAL, ALLOCATABLE :: spatial_corr(:,:)
INTEGER, ALLOCATABLE :: dates_sim(:)
CHARACTER(4000) :: headerline
CHARACTER(50) :: formatstring
!CHARACTER(50) :: headerformat
CHARACTER(11), ALLOCATABLE :: headerpieces(:)
INTEGER :: ia, it
INTEGER :: arg_count
CHARACTER(50) :: configfilename
INTEGER :: maxit
INTEGER :: mem
INTEGER :: archimem
INTEGER :: simmem
INTEGER :: maxchunk
INTEGER :: numberofachunks
INTEGER :: numberofschunks
INTEGER :: ichunk
INTEGER :: schunk
INTEGER :: tstart
INTEGER, ALLOCATABLE :: index_vec(:)
INTEGER, ALLOCATABLE :: sorted_ind(:)
INTEGER :: anastart
INTEGER :: anaend
INTEGER :: timestart
INTEGER :: timeend
INTEGER :: splitfactor

CALL DATE_AND_TIME(clockdate, clocktime)
PRINT*, clockdate, ' ', clocktime
! how many arguments were given to the file?
arg_count = COMMAND_ARGUMENT_COUNT()
SELECT CASE (arg_count)
 CASE (0)
  configfilename = "config.txt"
 CASE DEFAULT
  call GET_COMMAND_ARGUMENT(1, configfilename)
END SELECT
PRINT*, configfilename
! read configuration file
configs = get_configuration(TRIM(configfilename))
IF (.NOT. configs%param%silent) THEN
 PRINT*, "read config"
 PRINT*, TRIM(configs%files%outputfile)
 PRINT*, TRIM(configs%files%archivefile)
END IF
! get memory information
CALL SYSTEM('grep MemTotal /proc/meminfo |cut -c11-25 > fort.24')
READ(24,*) mem
PRINT*, mem
! get dimensions of the input data files
dim_archi = get_dims(TRIM(configs%files%archivefile))
print*, dim_archi
dim_sim = get_dims(TRIM(configs%files%simulationfile))
print*, dim_sim
IF (.NOT. configs%param%silent) PRINT*, "got dimensions"
! calculate memory needed to read the archive file (in kB, because memory in meminfo is given in kB)
archimem = dim_archi%lon_dim * dim_archi%lat_dim * dim_archi%time_dim * 0.008 
simmem = dim_sim%lon_dim * dim_sim%lat_dim * dim_sim%time_dim * 0.008 
IF (mem < 12000000) THEN
 splitfactor=4
ELSE IF (mem < 36000000) THEN
 splitfactor=8
ELSE 
 splitfactor = 8
END IF
IF (archimem+simmem > mem/REAL(splitfactor) ) THEN
 PRINT*, 'big archive and/or simulation'
 PRINT*, 'archive memory= ', archimem, 'simumlation memory= ', simmem, 'system memory= ', mem
 maxchunk = mem/(2*splitfactor) *125/dim_archi%lon_dim/dim_archi%lat_dim ! 125 8byte values can be stored per kb
 numberofachunks = CEILING(dim_archi%time_dim/REAL(maxchunk))
 numberofschunks = CEILING(dim_sim%time_dim/REAL(maxchunk))
 PRINT*, 'maxchunk= ', maxchunk, 'numberofachunks= ', numberofachunks, 'numberofschunks= ' , numberofschunks
ELSE
 numberofachunks = 1
 numberofschunks = 1
END IF

 IF (configs%param%calccor) THEN
    ALLOCATE(headerpieces(3*configs%param%nanalog))
 ELSE
    ALLOCATE(headerpieces(2*configs%param%nanalog))
 END IF

IF (numberofachunks == 1 .AND. numberofschunks == 1) THEN
 PRINT*, 'one chunk' 
 
 ALLOCATE(analogue_dates(configs%param%nanalog, dim_sim%time_dim), &
   & distances(configs%param%nanalog, dim_sim%time_dim), &
   & spatial_corr(configs%param%nanalog, dim_sim%time_dim), dates_sim(dim_sim%time_dim))
! call the main routine that reads the data and triggers the calculation
 CALL mainsub(dim_archi, dim_sim, configs%param%nanalog, analogue_dates, &
  & distances, spatial_corr, configs%param%silent, configs%param%varname, &
  & configs%files%archivefile, configs%files%simulationfile, &
  & configs%param%seacyc, configs%files%seacycfilebase, &
  & configs%files%seacycfilesim, configs%param%cycsmooth, &
  & configs%param%calccor, configs%param%distfun, configs%param%seasonwin, &
  & configs%param%timewin, dates_sim)
ELSE ! if the calculation has to be done in chunks
 ALLOCATE(analogue_dates(configs%param%nanalog*numberofachunks, dim_sim%time_dim), &
   & distances(configs%param%nanalog*numberofachunks, dim_sim%time_dim), &
   & spatial_corr(configs%param%nanalog*numberofachunks, dim_sim%time_dim), &
   & index_vec(configs%param%nanalog*numberofachunks), dates_sim(dim_sim%time_dim), &
   & sorted_ind(configs%param%nanalog))
 dim_archi_tmp = dim_archi
 dim_sim_tmp = dim_sim
 DO schunk =1, numberofschunks ! simulation chunk loop
!  prepare the simulation start and count numbers
  IF (schunk == numberofschunks) THEN
    dim_sim_tmp%time_dim = dim_sim%time_dim - (numberofschunks-1)*maxchunk
  ELSE
    dim_sim_tmp%time_dim = maxchunk
  END IF
  tstart = (schunk-1)*maxchunk + 1
  dim_sim_tmp%ncstart = (/1,1, tstart/)
  dim_sim_tmp%nccount = (/dim_sim_tmp%lon_dim, dim_sim_tmp%lat_dim, dim_sim_tmp%time_dim/)
  timestart = tstart
  timeend = tstart+dim_sim_tmp%time_dim-1
  PRINT*, dim_sim
  PRINT*, dim_sim_tmp
  PRINT*, timestart, timeend
  !STOP
  DO ichunk = 1, numberofachunks ! archive chunk loop
! prepare the archive start and count numbers
   IF (ichunk == numberofachunks) THEN
    dim_archi_tmp%time_dim = dim_archi%time_dim - (numberofachunks-1)*maxchunk
   ELSE
    dim_archi_tmp%time_dim = maxchunk
   END IF
   tstart = (ichunk-1)*maxchunk + 1
  !PRINT*, 'tstart = ', tstart
   dim_archi_tmp%ncstart = (/1,1, tstart/)
   dim_archi_tmp%nccount = (/dim_archi_tmp%lon_dim, dim_archi_tmp%lat_dim, dim_archi_tmp%time_dim/)
   anastart = (ichunk-1)*configs%param%nanalog +1
   anaend = ichunk*configs%param%nanalog
! calculate analogues in different archive sections and put them to different array sections
   CALL mainsub(dim_archi_tmp, dim_sim_tmp, configs%param%nanalog, &
    & analogue_dates(anastart:anaend,timestart:timeend), &
    & distances(anastart:anaend,timestart:timeend), spatial_corr(anastart:anaend,timestart:timeend), &
    & configs%param%silent, configs%param%varname, &
    & configs%files%archivefile, configs%files%simulationfile, &
    & configs%param%seacyc, configs%files%seacycfilebase, &
    & configs%files%seacycfilesim, configs%param%cycsmooth, &
    & configs%param%calccor, configs%param%distfun, configs%param%seasonwin, &
    & configs%param%timewin, dates_sim(timestart:timeend))
  PRINT*, 'end mainsub', ichunk
  END DO
 END DO
 ! do final sorting
 index_vec = (/ (ia, ia=1, configs%param%nanalog*numberofachunks) /)
 DO it = 1, dim_sim%time_dim
  CALL sort_index(index_vec, distances(:,it), &
    & configs%param%nanalog*numberofachunks, sorted_ind, &
    & distances(1:configs%param%nanalog,it), configs%param%nanalog)
  analogue_dates(1:configs%param%nanalog,it) = analogue_dates(sorted_ind,it)
  spatial_corr(1:configs%param%nanalog,it) = spatial_corr(sorted_ind,it)
 END DO
 PRINT*, 'end final sorting'
END IF

SELECT CASE (TRIM(configs%param%oformat))
 CASE (".nc")
  maxit = dim_sim%time_dim-configs%param%timewin+1
  ! write output with correlations
  IF (configs%param%calccor) THEN
   CALL write_dates_dists(TRIM(configs%files%outputfile), maxit, &
    & configs%param%nanalog, dates_sim(1:maxit), analogue_dates(1:configs%param%nanalog, 1:maxit), &
    & distances(1:configs%param%nanalog, 1:maxit), &
    &  TRIM(configs%param%distfun), configs%atts, spatial_corr(1:configs%param%nanalog, 1:maxit))
  ELSE
   CALL write_dates_dists(TRIM(configs%files%outputfile), maxit, &
    & configs%param%nanalog, dates_sim(1:maxit), analogue_dates(1:configs%param%nanalog, 1:maxit), &
    & distances(1:configs%param%nanalog, 1:maxit), &
    &  TRIM(configs%param%distfun), configs%atts)
  END IF
 CASE DEFAULT
 ! write output with correlations
  IF (configs%param%calccor) THEN
!PRINT*, 'construct header'
 ! construct header
   DO ia=1,configs%param%nanalog
    IF (ia < 10 ) THEN
!PRINT*,'ia= ', ia, configs%param%nanalog
     WRITE(headerpieces(ia), '(A, I1)') 'date.an', ia
     WRITE(headerpieces(ia+configs%param%nanalog), '(A, I1)') 'dis', ia
     WRITE(headerpieces(ia+2*configs%param%nanalog), '(A, I1)') 'cor', ia
    ELSE IF (ia < 100) THEN
!PRINT*, 'ia= ', ia
     WRITE(headerpieces(ia), '(A, I2)') 'date.an', ia
     WRITE(headerpieces(ia+configs%param%nanalog), '(A, I2)') 'dis', ia
     WRITE(headerpieces(ia+2*configs%param%nanalog), '(A, I2)') 'cor', ia
    ELSE 
     WRITE(headerpieces(ia), '(A, I3)') 'date.an', ia
     WRITE(headerpieces(ia+configs%param%nanalog), '(A, I3)') 'dis', ia
     WRITE(headerpieces(ia+2*configs%param%nanalog), '(A, I3)') 'cor', ia
    END IF
   END DO 
! write header
! WRITE(headerformat,'(A,I3.3,A)') '(',configs%param%nanalog*3 + 1,'A)'
! PRINT*, headerformat, size(headerpieces)
   WRITE(headerline,*) 'date ', headerpieces
 PRINT*, 'headerline written' 
   WRITE(formatstring,'(A,3(I3.3,A))') '(I9,',configs%param%nanalog, 'I9,', &
    & configs%param%nanalog,'F19.3,',configs%param%nanalog,'F13.7 )'
 PRINT*, formatstring
   OPEN(22,FILE=TRIM(configs%files%outputfile))
    WRITE(22,'(A)') TRIM(headerline)
  IF (.NOT. configs%param%silent)  PRINT*, 'header written'
! write data
    DO it=1,dim_sim%time_dim-configs%param%timewin+1
     WRITE(22,TRIM(formatstring)) &
      & dates_sim(it), analogue_dates(1:configs%param%nanalog, it), &
      & distances(1:configs%param%nanalog, it), spatial_corr(1:configs%param%nanalog, it)
    END DO
   CLOSE(22)
! write output without correlations
  ELSE
! construct header
   DO ia=1,configs%param%nanalog
    IF (ia < 10) THEN
     WRITE(headerpieces(ia), '(A, I1)') 'date.an', ia
     WRITE(headerpieces(ia+configs%param%nanalog), '(A, I1)') 'dis', ia
    ELSE IF (ia < 100 ) THEN
     WRITE(headerpieces(ia), '(A, I2)') 'date.an', ia
     WRITE(headerpieces(ia+configs%param%nanalog), '(A, I2)') 'dis', ia
    ELSE
     WRITE(headerpieces(ia), '(A, I3)') 'date.an', ia
     WRITE(headerpieces(ia+configs%param%nanalog), '(A, I3)') 'dis', ia
    END IF
   END DO 
! write header
   WRITE(headerline,*) 'date ', headerpieces
   WRITE(formatstring,'(A,3(I3.3,A))') '(I9,',configs%param%nanalog, 'I9,', configs%param%nanalog,'F19.3 )'
   OPEN(22,FILE=TRIM(configs%files%outputfile))
    WRITE(22,'(A)') headerline
! write data
    DO it=1,dim_sim%time_dim-configs%param%timewin+1
     WRITE(22,TRIM(formatstring)) &
      & dates_sim(it), analogue_dates(1:configs%param%nanalog, it), &
      & distances(1:configs%param%nanalog, it)
    END DO
   CLOSE(22)
  END IF
END SELECT

CALL DATE_AND_TIME(clockdate, clocktime)
PRINT*, clockdate, ' ', clocktime

END PROGRAM analogue
