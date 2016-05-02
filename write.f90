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
MODULE write_file
USE config
USE netcdf

CONTAINS
!> write output to netcdf file 
SUBROUTINE write_dates_dists(filename, numberofsimdays, numberofanalogues, &
 & simdates, anadates, distance, distname, atts, rankcor)
IMPLICIT NONE
CHARACTER(*) :: filename
INTEGER :: numberofsimdays
INTEGER :: numberofanalogues
INTEGER :: simdates(numberofsimdays)
INTEGER :: anadates(numberofanalogues, numberofsimdays)
REAL(8) :: distance(numberofanalogues, numberofsimdays)
CHARACTER(*) :: distname
TYPE(atts_type) :: atts
REAL, OPTIONAL :: rankcor(numberofanalogues, numberofsimdays)

INTEGER :: state
INTEGER :: ncid
INTEGER :: nana_dimid
INTEGER :: time_dimid
INTEGER :: nana_varid
INTEGER :: t_varid
INTEGER :: anadate_varid
INTEGER :: dist_varid
INTEGER :: cor_varid
INTEGER :: m
 
state = NF90_CREATE(TRIM(filename), NF90_HDF5, ncid) ! create netCDF dataset: enter define mode
 IF (state /= NF90_NOERR) THEN
  PRINT*, 'error: failed to create netcdf file'
 ELSE
 ! define dimensions
  state = NF90_DEF_DIM(ncid, 'numberofanalogues', numberofanalogues, nana_dimid) ! define dimensions: from name and length
  IF (state /= NF90_NOERR) CALL print_err(state)
  state = NF90_DEF_DIM(ncid, 'time', NF90_UNLIMITED, time_dimid) ! define time dimension
   IF (state /= NF90_NOERR) CALL print_err(state)
 ! define variables 
  state = NF90_DEF_VAR(ncid, 'numberofanalogues', NF90_INT, nana_dimid, nana_varid)! define variables: from name, type, dims
   IF (state /= NF90_NOERR) CALL print_err(state)
  state = NF90_DEF_VAR(ncid, 'time', NF90_INT, time_dimid, t_varid) ! define variables: from name, type, dims
   IF (state /= NF90_NOERR) CALL print_err(state)
  state = NF90_DEF_VAR(ncid, 'anadates', NF90_INT, (/nana_dimid, time_dimid/), anadate_varid)! define variables: from name, type, dims
   IF (state /= NF90_NOERR) CALL print_err(state)
  state = NF90_DEF_VAR(ncid, 'distance', NF90_FLOAT, (/nana_dimid, time_dimid/), dist_varid) ! define variables: from name, type, dims
   IF (state /= NF90_NOERR) CALL print_err(state)
  IF (PRESENT(rankcor)) THEN
   state = NF90_DEF_VAR(ncid, 'correlation', NF90_FLOAT, (/nana_dimid, time_dimid/), cor_varid) ! define variables: from name, type, dims
    IF (state /= NF90_NOERR) CALL print_err(state)
  END IF
 ! attributes
  state = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'title', 'CASTf90 analogue dates') ! assign attribute values
   IF (state /= NF90_NOERR) CALL print_err(state)
  state = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'institution', 'LSCE') ! assign attribute values
   IF (state /= NF90_NOERR) CALL print_err(state)
  state = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'data_source_simulation', TRIM(atts%simsource)) ! assign attribute values
   IF (state /= NF90_NOERR) CALL print_err(state)
  state = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'data_source_archive', TRIM(atts%archisource)) ! assign attribute values
 IF (state /= NF90_NOERR) CALL print_err(state)
  state = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'archive_period', TRIM(atts%archiperiod)) ! assign attribute values
 IF (state /= NF90_NOERR) CALL print_err(state)
  state = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'predictor_variable', TRIM(atts%predictorvar))! assign attribute values
  IF (state /= NF90_NOERR) CALL print_err(state)
  state = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'predictor_domain', TRIM(atts%predictordom))! assign attribute values
  IF (state /= NF90_NOERR) CALL print_err(state)
  state = NF90_PUT_ATT(ncid, t_varid, 'long_name', 'time') ! assign attribute values
   IF (state /= NF90_NOERR) CALL print_err(state)
  state = NF90_PUT_ATT(ncid, anadate_varid, 'long_name', 'analogue dates') ! assign attribute values
  IF (state /= NF90_NOERR) CALL print_err(state)
  state = NF90_PUT_ATT(ncid, dist_varid, 'long_name', TRIM(distname)) ! assign attribute values
  IF (state /= NF90_NOERR) CALL print_err(state)
  IF (PRESENT(rankcor)) THEN
   state = NF90_PUT_ATT(ncid, cor_varid, 'long_name', 'rank correlation') ! assign attribute values
   IF (state /= NF90_NOERR) CALL print_err(state)
  END IF
  state = NF90_ENDDEF(ncid) ! end definitions: leave define mode
  IF (state /= NF90_NOERR) CALL print_err(state)
  state = NF90_PUT_VAR(ncid, nana_varid, (/ (m, m=1,numberofanalogues) /)) ! provide values for variable
   IF (state /= NF90_NOERR) CALL print_err(state)
  state = NF90_PUT_VAR(ncid, t_varid, simdates)
   IF (state /= NF90_NOERR) CALL print_err(state)
  state = NF90_PUT_VAR(ncid, anadate_varid, anadates)
  IF (state /= NF90_NOERR) CALL print_err(state)
  state = NF90_PUT_VAR(ncid, dist_varid, distance)
   IF (state /= NF90_NOERR) CALL print_err(state)
  IF (PRESENT(rankcor)) THEN
   state = NF90_PUT_VAR(ncid, cor_varid, rankcor)
   IF (state /= NF90_NOERR) CALL print_err(state)
  END IF
 state =  NF90_CLOSE(ncid)! close: save new netCDF dataset
 END IF 
END SUBROUTINE write_dates_dists

!********************************************

subroutine print_err(state)
IMPLICIT NONE
integer, intent ( in) :: state

print *, trim(nf90_strerror(state))
stop "Stopped"

end subroutine print_err

!*******************8
END MODULE write_file

