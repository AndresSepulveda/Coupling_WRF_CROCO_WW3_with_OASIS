# Coupling_WRF_CROCO_WW3_with_OASIS
Notes on coupling WRF, CROCO, and WW3 using OASIS for a UBUNTU 18.04 machine

A mix of "Documentation for coupling with OASIS in CROCO, WRF, WW3" by Swen JULLIEN, Gildas CAMBON (March 7, 2018) and
personal experience.

# OASIS

Donwload OASIS-MCT version 3

```console
cd oasis
svn checkout https://oasis3mct.cerfacs.fr/svn/branches/OASIS3-MCT_3.0_branch
cd $HOME/oasis/OASIS3 -MCT_3.0 _branch/oasis3 -mct/util/make_dir
make realclean -f TopMakefileOasis3 > oasis_clean.out
make -f TopMakefileOasis3 > oasis_make.out
```
Check directory

```console
compile_oa3-mct
```
and file oasis_make.out to see if everything went ok


# OASIS - namcouple file

```console
#########################################################################
# This is a typical input file for OASIS3-MCT.
# Keywords used in previous versions of OASIS3 
# but now obsolete are marked "Not used"
# Don't hesitate to ask precisions or make suggestions (oasishelp@cerfacs.fr). 
#
# Any line beginning with # is ignored. Blank lines are not allowed.
#
#########################################################################
#
# NFIELDS: total number of fields being exchanged
 $NFIELDS
 9
#########################################################################
# NBMODEL: number of models and their names (6 characters) 
 $NBMODEL
 2  wrfexe  crocox
###########################################################################
# RUNTIME: total simulated time for the actual run in seconds (<I8)
 $RUNTIME
 86400
###########################################################################
# NLOGPRT: debug and time statistics informations printed in log file 
#          First number: 0 one log file for master, and one for other procs
#                        1 one log file for master, and one for other errors
#                        2 one file per proc with normal diagnostics
#                        5 as 2 + initial debug info
#                        10 as 5 + routine calling tree
#                        12 as 10 + some routine calling notes
#                        15 as 12 + even more debug diagnostics
#                        20 as 15 + some extra runtime analysis
#                        30 full debug information
#          Second number: time statistics
#          		 0 nothing calculated
#          		 1 one file for proc 0 and min/max of other procs
#          		 2 as 1 + one file per proc
#          		 3 as 2 + proc 0 writes all procs results in its file
 $NLOGPRT
 1 1
###########################################################################
# Beginning of fields exchange definition
 $STRINGS
#
# For each exchanged field:
#
# line 1: field in sending model, field in target model, unused, coupling 
#         period, number of transformation, restart file, field status
# line 2: nb of pts for sending model grid (without halo) first dim, and second dim,
#         for target grid first dim, and second dim, sending model grid name, target 
#         model grid name, lag = time step of sending model
# line 3: sending model grid periodical (P) or regional (R), and nb of overlapping 
#         points, target model grid periodical (P) or regional (R), and number of
#         overlapping points
# line 4: list of transformations performed
# line 5: parameters for each transformation
#
# See the correspondances between variables in models and in OASIS:
# Note: for CROCO and WRF nesting capability is useable in coupled 
#       mode. For CROCO the domain in defined by the last number 
#       of coupled field name. For WRF, WRF domain is defined by
#       the number after WRF_d, and the domain of the coupled model
#       (CROCO for example) is given by EXT_d in coupled field name 
#
# |--------------------------------------------------------------|
# | Possibly sent fields by CROCO:                 CROCO | OASIS |
# |--------------------------------------------------------------|
# |     t(:,:,N,nnew,itemp)  |    SRMSSTV0                       |
# |                   zeta   |    SRMSSHV0                       |
# |     u v (at rho points)  |    SRMVOCE0 SRMUOCE0              |
# |--------------------------------------------------------------|
# | Possibly received fields by CROCO:            CROCO | OASIS  |
# |--------------------------------------------------------------|
# |                  srflx   |    RRMSRFL0                       |
# |       stflx(:,:,isalt)   |    RRMEVPR0                       |
# |      stflx(,:,:,itemp)   |    RRMSTFL0                       |
# |                  sustr   |    RRMTAUX0                       |
# |                  svstr   |    RRMTAUY0                       |
# |                  smstr   |    RRMTAUM0                       |
# |                  whrm    |    RRM__HS0                       |
# |                  wfrq    |    RRMT0M10                       |
# |                  wdrx    |    RRMCDIR0                       |
# |                  wdre    |    RRMSDIR0                       |
# |--------------------------------------------------------------|
# | Possibly sent fields by WW3:                    WW3 | OASIS  |
# |--------------------------------------------------------------|
# |            not defined   |    WW3_ODRY                       |
# |                   T0M1   |    WW3_T0M1                       |
# |                     HS   |    WW3__OHS                       |
# |                    DIR   |    WW3_CDIR WW3_SDIR              |
# |                    BHD   |    WW3__BHD                       |
# |                    TWO   |    WW3_TWOX WW3_TWOY              |
# |                    UBR   |    WW3__UBR                       |
# |                    FOC   |    WW3__FOC                       |
# |                    TAW   |    WW3_TAWX WW3_TAWY              |
# |                     LM   |    WW3___LM                       |
# |                    CUR   |    WW3_WSSU WW3_WSSV              |
# |                    CHA   |    WW3__CHA                       |
# |                     HS   |    WW3__AHS                       |
# |                     FP   |    WW3___FP                       |
# |--------------------------------------------------------------|
# | Possibly received fields by WW3:                WW3 | OASIS  |
# |--------------------------------------------------------------|
# |            not defined   |    WW3_OWDH WW3_OWDU WW3_OWDV     |
# |                    SSH   |    WW3__SSH                       |
# |                    CUR   |    WW3_OSSU WW3_OSSV              |
# |                    WND   |    WW3__U10 WW3__V10              |
# |--------------------------------------------------------------|
# | Possibly sent fields by WRF:                    WRF | OASIS  |
# |--------------------------------------------------------------|
# |                GSW   |    WRF_d01_EXT_d01_SURF_NET_SOLAR     |
# |        QFX-(RAINCV                                           |
# |       +RAINNCV)/DT   |    WRF_d01_EXT_d01_EVAP-PRECIP        |
# |   GLW-STBOLT*EMISS                                           |
# |     *SST**4-LH-HFX   |    WRF_d01_EXT_d01_SURF_NET_NON-SOLAR |
# | taut * u_uo / wspd   |    WRF_d01_EXT_d01_TAUX               |
# | taut * u_uo / wspd   |    WRF_d01_EXT_d01_TAUY               |
# |               taut   |    WRF_d01_EXT_d01_TAUMOD             |
# |               u_uo   |    WRF_d01_EXT_d01_U_01               |
# |               v_vo   |    WRF_d01_EXT_d01_V_01               |
# |--------------------------------------------------------------|
# | Possibly received fields by WRF:                WRF | OASIS  |
# |--------------------------------------------------------------|
# |                    SST   |    WRF_d01_EXT_d01_SST            |
# |                   UOCE   |    WRF_d01_EXT_d01_UOCE           |
# |                   VOCE   |    WRF_d01_EXT_d01_VOCE           |
# |               CHA_COEF   |    WRF_d01_EXT_d01_CHA_COEF       |
# |--------------------------------------------------------------|
#
#                     ------------------------------------
#                       WRF (wrfexe) ==> CROCO (crocox)
#                     ------------------------------------
#
#~~~~~~~~~~~
# TAUX : zonal stress (N.m-2)
#~~~~~~~~~~~
WRF_d01_EXT_d01_TAUX RRMTAUX0 1 3600 1 atm.nc EXPORTED
56 50 41 42 atmt ocnt  LAG=150
R  0  R  0
SCRIPR
DISTWGT LR SCALAR LATLON 1 4
#
#~~~~~~~~~~~
# TAUY : meridional stress (N.m-2)
#~~~~~~~~~~~
WRF_d01_EXT_d01_TAUY RRMTAUY0 1 3600 1 atm.nc EXPORTED
56 50 41 42 atmt ocnt  LAG=150
R  0  R  0
SCRIPR
DISTWGT LR SCALAR LATLON 1 4
#
#~~~~~~~~~~~
# TAUMOD : stress module (N.m-2)
#~~~~~~~~~~~
WRF_d01_EXT_d01_TAUMOD RRMTAUM0 1 3600 1 atm.nc EXPORTED
56 50 41 42 atmt ocnt  LAG=150
R  0  R  0
SCRIPR
DISTWGT LR SCALAR LATLON 1 4
#
#~~~~~~~~~~~
# EVAP-PRECIP : E-P flux (kg.m-2.s-1)
#~~~~~~~~~~~
WRF_d01_EXT_d01_EVAP-PRECIP RRMEVPR0 1 3600 1 atm.nc EXPORTED
56 50 41 42 atmt ocnt  LAG=150
R  0  R  0
SCRIPR
DISTWGT LR SCALAR LATLON 1 4
#
#~~~~~~~~~~~
# NET SOLAR FLUX (W.m-2)
#~~~~~~~~~~~
WRF_d01_EXT_d01_SURF_NET_SOLAR RRMSRFL0 1 3600 1 atm.nc EXPORTED
56 50 41 42 atmt ocnt  LAG=150
R  0  R  0
SCRIPR
DISTWGT LR SCALAR LATLON 1 4
#
#~~~~~~~~~~~
# NET NON-SOLAR FLUX (W.m-2)
#~~~~~~~~~~~
WRF_d01_EXT_d01_SURF_NET_NON-SOLAR RRMSTFL0 1 3600 1 atm.nc EXPORTED
56 50 41 42 atmt ocnt  LAG=150
R  0  R  0
SCRIPR
DISTWGT LR SCALAR LATLON 1 4
#
#                     ------------------------------------
#                       CROCO (crocox) ==> WRF (wrfexe)
#                     ------------------------------------
#
#~~~~~~~~~~~
# SST (K)
#~~~~~~~~~~~
SRMSSTV0 WRF_d01_EXT_d01_SST 1 3600 1 oce.nc EXPORTED
41 42 56 50 ocnt atmt LAG=3600
R  0  R  0
SCRIPR
DISTWGT LR SCALAR LATLON 1 4
#
#~~~~~~~~~~~
# UOCE : zonal current (m.s-1)
#~~~~~~~~~~~
SRMUOCE0 WRF_d01_EXT_d01_UOCE 1 3600 1 oce.nc EXPORTED
41 42 56 50 ocnt atmt LAG=3600
R  0  R  0
SCRIPR
DISTWGT LR SCALAR LATLON 1 4
#
#~~~~~~~~~~~
# VOCE : meridonal current (m.s-1)
#~~~~~~~~~~~
SRMVOCE0 WRF_d01_EXT_d01_VOCE 1 3600 1 oce.nc EXPORTED
41 42 56 50 ocnt atmt LAG=3600
R  0  R  0
SCRIPR
DISTWGT LR SCALAR LATLON 1 4
#
###########################################################################
$END
```
# OASIS - Create OASIS grid files form WRF

```console
#!/bin/bash
set -x
## ----------------------------------------------------------------------------- #
## - Create grids.nc, masks.nc, files from WRF for oasis                       - #
## - because call to oasis_grid function not yet implemented in WRF            - #
##                                                                             - #
## Mandatory inputs:                                                           - #
##  - a file from WRF containing lon,lat,mask (with full path)                 - #
##  - the output destination directory                                         - #
##                                                                             - #
## ----------------------------------------------------------------------------- #
#
# Further Information:   
# http://www.croco-ocean.org
#  
# This file is part of CROCOTOOLS
#
# CROCOTOOLS is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation; either version 2 of the License,
# or (at your option) any later version.
#
# CROCOTOOLS is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA
#
# Copyright (c) 2018 S. Jullien
# swen.jullien@ifremer.fr
## ----------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------ #

gridfile=$1
mydir=$2

# ------------------------------------------------------------------------------ #

echo '*******************************************'
echo 'START script create_oasis_grids_for_wrf.sh'
echo '*******************************************'
echo ' '

# First check if inputs are ok
if [[ -z $gridfile ]] || [[ -z $mydir ]] ; then
    echo 'ERROR: inputs are not correctly specified.'
    echo '       this script needs 2 inputs:'
    echo '       - a file from WRF containing the mask, lon and lat (with full path)'
    echo '       - the output destination directory'
    echo ' Exit...'
    echo ' '
    exit 1
fi
# ------------------------------------------------------------------------------ #

mytmpgrd=$mydir/grd_tmp.nc
grdfile=$mydir/grids.wrf.nc
mskfile=$mydir/masks.wrf.nc

wrflon=XLONG
wrflat=XLAT
wrfmask=LANDMASK

# Extract lon,lat,mask
echo '---> Extract '${wrflon}', '${wrflat}', and '${wrfmask}' variables...'
ncks -O -v ${wrflon},${wrflat},${wrfmask} -d Time,0 $gridfile ${mytmpgrd}

# Put them on the stag grid
echo '---> Put them on the stag grid' 
./to_wrf_stag_grid.sh ${mytmpgrd} ${mytmpgrd}

# remove time dimension
echo '---> Remove time dimension...'
ncwa -O -a Time ${mytmpgrd} ${mytmpgrd}

# compute the last lon and lat
Nlon=`ncdump -h $gridfile | grep "west_east = " | cut -d '=' -f2 | cut -d ';' -f1`
Nlon=${Nlon// /}
Nlat=`ncdump -h $gridfile | grep "south_north = " | cut -d '=' -f2 | cut -d ';' -f1`
Nlat=${Nlat// /}
Nlonstag=$(($Nlon + 1))
Nlatstag=$(($Nlat + 1))
Nlonm1=$(($Nlon - 1))
Nlatm1=$(($Nlat - 1))
echo '---> compute the last lon...'
ncap2 -F -O -s "${wrflon}(:,$Nlonstag)=${wrflon}(:,$Nlon)+(${wrflon}(:,$Nlon)-${wrflon}(:,$Nlonm1))" ${mytmpgrd} ${mytmpgrd}
echo '---> compute the last lat...'
ncap2 -F -O -s "${wrflat}($Nlatstag,:)=${wrflat}($Nlat,:)+(${wrflat}($Nlat,:)-${wrflat}($Nlatm1,:))" ${mytmpgrd} ${mytmpgrd}

# change mask from float to integer
echo '---> Change mask from float to integer...'
ncap2 -O -s "${wrfmask}=int(${wrfmask})" ${mytmpgrd} ${mytmpgrd}

# rename dimensions
echo '---> rename dimensions...'
ncrename -d west_east,x_atmt -d south_north,y_atmt ${mytmpgrd}
# rename variables
echo '---> Rename variables...'
ncrename -v ${wrfmask},atmt.msk  -v ${wrflat},atmt.lat ${mytmpgrd} 
ncrename -v ${wrflon},atmt.lon ${mytmpgrd} 

# create grid file
echo '---> Ceate grid file...'
echo '======================='
ncks -O -v atmt.lon,atmt.lat ${mytmpgrd} ${grdfile}
ncatted -h -O -a ,global,d,, ${grdfile} ${grdfile}
ncatted -h -O -a ,atmt.lon,d,, ${grdfile} ${grdfile}
ncatted -h -O -a ,atmt.lat,d,, ${grdfile} ${grdfile}

# create mask file
echo '---> Create mask file...'
echo '========================='
ncks -O -v atmt.msk ${mytmpgrd} ${mskfile}
ncatted -h -O -a ,global,d,, ${mskfile} ${mskfile}
ncatted -h -O -a ,atmt.msk,d,, ${mskfile} ${mskfile}

rm ${mytmpgrd}

echo 'DONE: grids.wrf.nc and masks.wrf.nc have been created in '$mydir
echo ' '
```

# OASIS - Create restart from calm conditions

```console
#!/bin/bash -e

## ----------------------------------------------------------------------------- #
## - Create restart file for oasis                                             - #
## - with all variables set to 0                                               - #
##                                                                             - #
## Mandatory inputs:                                                           - #
##  - a file from this model containing the mask (with full path)              - #
##  - the oasis restart file name (with full path)                             - #
##  - the model: wrf, croco, or ww3 cases are accepted                         - #
##  - the list of variables that have to be generated in this restart file     - #
##                                                                             - #
## ----------------------------------------------------------------------------- #
#
# Further Information:   
# http://www.croco-ocean.org
#  
# This file is part of CROCOTOOLS
#
# CROCOTOOLS is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation; either version 2 of the License,
# or (at your option) any later version.
#
# CROCOTOOLS is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA
#
# Copyright (c) 2018 S. Jullien
# swen.jullien@ifremer.fr
## ----------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------ #

filein=$1
fileout=$2
model=$3
varlist=$4

# ------------------------------------------------------------------------------ #
echo '**********************************************************'
echo 'START script create_oasis_restart_from_calm_conditions.sh'
echo '**********************************************************'
echo ' '

# First check if inputs are ok
if [[ -z $filein ]] || [[ -z $fileout ]] || [[ -z $model ]] || [[ -z $varlist ]] ; then
    echo 'ERROR: inputs are not correctly specified.'
    echo '       this script needs 4 inputs:'
    echo '       - a file from this model containing the mask (with full path)'
    echo '       - the oasis restart file name (with full path)'
    echo '       - the model: wrf croco or ww3 cases are accepted'
    echo '       - the list of variables that have to be generated in this restart file'
    echo ' Exit...'
    echo ' '
    exit 1
fi 
# ------------------------------------------------------------------------------ #

mydir=$(dirname "$fileout")
filetmp=$mydir/rst_tmp.nc

echo 'Initialize restart file '$fileout' to 0 for variables: '$varlist
echo '========================================================================='

if [ $model == wrf ] ; then

    # Extract mask
    echo '---> Extract LANDMASK variable...'
    ncks -O -v LANDMASK -d Time,0 $filein ${filetmp}
    
    # Put it on the stag grid
    echo '---> Put it on the stag grid' 
    ./to_wrf_stag_grid.sh ${filetmp} ${filetmp}
 
    # remove time dimension
    echo '---> Remove time dimension...'
    ncwa -O -a Time ${filetmp} ${filetmp}
    
    # set the variable to 0 and rename it var0
    echo '---> Set the variable to 0 and rename it var0...'
    ncap2 -O -v -s "var0=LANDMASK*0" ${filetmp} ${filetmp}
    
elif  [ $model == croco ] ; then

    # Extract dimensions
    xi_rho=`ncdump -h $filein | grep "xi_rho = " | cut -d '=' -f2 | cut -d ';' -f1`
    xi_rho=${xi_rho// /}
    eta_rho=`ncdump -h $filein | grep "eta_rho = " | cut -d '=' -f2 | cut -d ';' -f1`
    eta_rho=${eta_rho// /}

    # Extract mask
    echo '---> Extract mask_rho variable (only interior grid indices, i.e. in fortran convention : 2:end-1)...'
    ncks -O -F -d xi_rho,2,$((${xi_rho}-1)) -d eta_rho,2,$((${eta_rho}-1)) -v mask_rho $filein ${filetmp}

    # set the variable to 0 and rename it var0
    echo '---> Set the variable to 0 and rename it var0...'
    ncap2 -O -v -s "var0=mask_rho*0" ${filetmp} ${filetmp}

elif  [ $model == ww3 ] ; then
    
    # Extract mask
    echo '---> Extract MAPSTA variable...'
    ncks -O -3 -v MAPSTA $filein ${filetmp}

    # set the variable to 0 and rename it var0
    echo '---> Set the variable to 0 and rename it var0...'
    ncap2 -O -v -s "var0=MAPSTA*0" ${filetmp} ${filetmp}

else
    
    echo 'ERROR: '$model' case is not implemented yet. Exit...'
    echo ' '
    exit 1

fi # model 
    
# START LOOP on varlist #
#-----------------------#
for var in $varlist ; do

    echo ' '
    echo '================================'
    echo 'Create variable: '$var'...'
    echo '================================'
    ncks -A -v var0 ${filetmp} $fileout
    ncrename -v var0,$var $fileout

done
rm ${filetmp}

# remove global attributes
ncatted -h -O -a ,global,d,, $fileout $fileout

echo 'DONE for '$varlist' variables => initialized to 0'
echo ' '
```



# WRF

# WRF - Create CPLMASK

Use something similar to this program named make_CPLMASK.m

```console
nc=netcdf('wrfinput_d01','w');
lm=nc{'LANDMASK'}(:,:,:);
lk=nc{'LAKEMASK'}(:,:,:);
cp=nc{'CPLMASK'}(:,:,:,:);
indx0=find(lm==0);
indxl=find(lk==1);
cp(indx0)=1;
cp(indxl)=0;
nc{'CPLMASK'}(:,1,:,:)=cp;
close(nc)
```

# WRF - namelist.input

To add variable SST update

in &time_control

```console
io_form_auxinput4        = 2,
auxinput4_interval       = 60,
auxinput4_inname         = "wrflowinp_d<domain>",
```

in &domains

```console
num_ext_model_couple_dom = 1,
```

in &physics

```console
sst_update               = 1,
```



# CROCO

Define MPI and OA_COUPLING in cppdefs.h

# CROCO - cppdefs.h

```console
                      /* Parallelization */
# undef  OPENMP
# define  MPI
                      /* I/O server */
# undef  XIOS
                      /* Non-hydrostatic option */
# undef  NBQ
                      /* Nesting */
# undef  AGRIF
# undef  AGRIF_2WAY
                      /* OA and OW Coupling via OASIS (MPI) */
# define  OA_COUPLING
# undef  OW_COUPLING
```

# CROCO - jobcomp

```console
NETCDFLIB="-L/home/mosa/libraries/netcdf/lib -lnetcdf"
NETCDFINC="-I/home/mosa/libraries/netcdf/include"
NETCDFLIB=$(nf-config --flibs)
NETCDFINC=-I$(nf-config --includedir)
#
# set MPI directories if needed
#
MPIF90="/usr/bin/mpif90"
MPILIB=""
MPIINC=""

#
# set OASIS-MCT (or OASIS3) directories if needed
#
PRISM_ROOT_DIR=/home/mosa/compile_oa3-mct

```

# WW3
