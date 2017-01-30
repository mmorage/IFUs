#!/usr/bin/env python
 
import sys,os,string

import astropy
from astropy.io import fits
import numpy as np
import pyfits
import scipy
from scipy import ndimage
import pyraf
from pyraf import iraf
from iraf import gemini, gmos, gemtools, gfreduce, unlearn,astutil,rv,gffindblocks

imagen_in=sys.argv[1]
SBias=sys.argv[2]
ARC=sys.argv[3]       
FLAT=sys.argv[4]
TWIL=sys.argv[5]

iraf.unlearn('gfreduce')


gmos.logfile="example.log"
gemtools.logfile="example.log"

# Load packages:
 

#seteo de gfextract!


gmos.gfextract.line='100'
gmos.gfextract.nsum='10'
gfreduce.line=100

gmos.gfextract.line=100
gmos.gfextract.nsum=10
gfreduce.line=100

#Running BIAS and OVERSCAN in all images aiming at the correction of the bad columns on the CCDs

iraf.gfreduce(inimages=imagen_in,outpref='default',fl_inter='no',fl_vardq='yes',fl_addmdf='yes',fl_over='yes',fl_trim='no',fl_gscrrej='no',fl_wavtran='no',fl_skysub='no',trace='no',recenter='no',bias=SBias,fl_bias='yes',fl_gsappwave='no',fl_fluxcal='no',fl_fulldq='no',fl_fixgap='no',line=100,nsum=10)

iraf.gfreduce(inimages=ARC,outpref='default',fl_inter='no',fl_vardq='yes',fl_addmdf='yes',fl_over='yes',fl_trim='no',fl_gscrrej='no',fl_wavtran='no',fl_skysub='no',trace='no',recenter='no',bias=SBias,fl_bias='yes',fl_gsappwave='no',fl_fluxcal='no',fl_fulldq='no',fl_fixgap='no',line=100,nsum=10)

iraf.gfreduce(inimages=FLAT,outpref='default',fl_inter='no',fl_vardq='yes',fl_addmdf='yes',fl_over='yes',fl_trim='no',fl_gscrrej='no',fl_wavtran='no',fl_skysub='no',trace='no',recenter='no',bias=SBias,fl_bias='yes',fl_gsappwave='no',fl_fluxcal='no',fl_fulldq='no',fl_fixgap='no',line=100,nsum=10)

 
#run cl script from German Gimeno that corrects the bad columns on the Amplifiers. 


###############
# ## OLD
# ##
# ##iraf.gfreduce(inimages=FLAT,outpref='default',fl_inter='no',fl_vardq='yes',fl_addmdf='yes',fl_over='yes',fl_trim='no',fl_gscrrej='no',fl_wavtran='no',fl_skysub='no',trace='no',recenter='no',bias=SBias,fl_bias='yes',fl_gsappwave='no',fl_fluxcal='no',fl_fulldq='no',fl_fixgap='no',line=100,nsum=10)
# ##
# ##
#########################


FLATg='rrg'+FLAT[:-5]+'_tmp.fits'
#FLATg='grg'+FLAT[:-5]+'.fits'

# SE USA!!
iraf.gfreduce(inimages=FLATg,outpref='REF_FLAT_',fl_inter='yes',fl_vardq='no',fl_addmdf='no',fl_over='no',fl_trim='no',fl_gscrrej='no',fl_wavtran='no',fl_skysub='no',trace='yes',recenter='yes',bias=SBias,fl_bias='no',fl_gsappwave='no',fl_fluxcal='no',fl_fulldq='no',fl_fixgap='no',line=100,nsum=10)

 
#################################################################################################################### 
# OLD
#
#  ###Scattered light subtraction
#
#  ##image= imagen bias y trimmed... i.e rgS.....
#  ##extspc= FLAT reference
#  ##mask=output
#  
#  #iraf.gffindblocks(image='rg'+FLAT,extspec='REF_FLAT_'+FLAT,mask=FLAT[:-5]+'_GAPS')#
#  #
#  
#  ##iraf.gfscatsub.outimage=""
#  ##iraf.gfscatsub.prefix="b"#
#  
#  ##iraf.gfscatsub.xorder="5,9,5"   # try 5,5,9,5,5,5 for new GMOS-N CCDs
#  ##iraf.gfscatsub.yorder="5,7,5"   # try 5,5,9,5,5,5 for new GMOS-N CCDs
#  ##iraf.gfscatsub.cross='yes'#
#  
#  ####iraf.gfscatsub(image='REF_FLAT_'+FLAT,mask='GAPS_'+FLAT,)
#  
#  
#  #iraf.gfscatsub(image='rg'+FLAT[:-5],prefix='b_',mask=FLAT[:-5]+'_GAPS',xorder="3",yorder="3",cross='yes')
###################################################################################################################

###ARC: 

ARCg='rrg'+ARC[:-5]+'_tmp.fits'

iraf.gfreduce(inimages=ARCg,outpref='ARC_',fl_inter='no',fl_vardq='no',fl_addmdf='no',fl_over='no',fl_trim='no',fl_gscrrej='no',fl_wavtran='no',fl_skysub='no',trace='no',recenter='no',fl_extract='yes',bias=SBias,fl_bias='no',orde='1',weight="none",fl_gsappwave='no',fl_fluxcal='no',ref='REF_FLAT_'+FLATg,fl_fulldq='yes')

#wavelenght calibration
iraf.gswavelength(inimages='ARC_'+ARCg,coordlist='CuAr_GMOS.dat',fl_inter='yes')
 
##Apply QE to the FLATS
#################
# OLD
#####iraf.gqecorr(inimages='b_rg'+FLAT,refimages='ARC_'+ARC[:-5],qecorr_data="gmosQEfactors.dat") 
#
#######################
FLATg='rrg'+FLAT[:-5]+'_tmp.fits'

iraf.gqecorr(inimages=FLATg,outpref='q_',refimages='ARC_'+ARCg[:-5],qecorr_data="gmosQEfactors.dat") 


##############################################################################
#  ###Repitiendo la extraccion con el flat corregido por scattered qb_rg...
# OLD
#  ####iraf.gfreduce(inimages='qtest_rg'+FLAT,outpref='FLAT_',fl_inter='no',fl_vardq='no',fl_addmdf='no',fl_over='no',fl_trim='no',fl_gscrrej='no',fl_wavtran='no',fl_skysub='no',trace='no',recenter='no',fl_extract='yes',bias=SBias,fl_bias='yes',fl_gsappwave='yes',fl_fluxcal='no',ref='REF_FLAT_'+FLAT,fl_fulldq='yes',fl_fixgap='no')
#################################################################################

iraf.gfreduce(inimages='q_'+FLATg,outpref='FLAT_',fl_inter='no',fl_vardq='no',fl_addmdf='no',fl_over='no',fl_trim='no',fl_gscrrej='no',fl_wavtran='no',fl_skysub='no',trace='no',recenter='no',fl_extract='yes',bias=SBias,fl_bias='yes',fl_gsappwave='yes',fl_fluxcal='no',ref='REF_FLAT_'+FLATg,fl_fulldq='yes',fl_fixgap='no')

###haciendo el responso... normalizando el flat

####################################################
##  RESPONSE. IT works sometimes, It deppends on the data
##
##  iraf.gfresponse(inimage='FLAT_qb_rg'+FLAT,outimage='RESP_'+FLAT,skyimage="",fl_inter='yes',function='spline3',order='45',sample="*",wavtraname='ARC_'+ARC[:-5])
##  iraf.gfresponse(inimage='FLAT_q_rg'+FLAT,outimage='RESP_'+FLAT,skyimage="",fl_inter='yes',function='spline3',order='45',sample="*",wavtraname='ARC_'+ARC[:-5]) 
##
##########################################################

#OK Check the iraf version here: The first one is for the standard gfresponse, the second one for the developer version from Gemini

#iraf.gfresponse(inimage='FLAT_q_'+FLATg,outimage='RESP_'+FLATg,skyimage="",fl_inter='yes',function='spline3',order='45',sample="*") 
iraf.gfresponse(inimage='FLAT_q_'+FLATg,outimage='RESP_'+FLATg,skyimage="",fl_inter='yes',function='spline3',order='45',sample="*",wavtraname='ARC_'+ARCg[:-5]) 



###############################################################
# OLD
# #### Preparando para remover scattered light from the science
# #### Primero corriendo gfreduce en Science solamente con rg ,es decir gprepared bias trimed y 
# 
# ###iraf.gfreduce(inimages=imagen_in,fl_inter='no',fl_vardq='yes',fl_addmdf='yes',fl_over='yes',fl_trim='yes',fl_gscrrej='no',fl_wavtran='no',fl_skysub='no',bias=SBias,fl_bias='yes',fl_gsappwave='no',fl_fluxcal='no')
# 
# ###problemas con el xorder 3#
# 
# ###iraf.gfscatsub(image='rg'+imagen_in[:-5],prefix='b',mask=FLAT[:-5]+'_GAPS',xorder="1",yorder="3",cross='yes')
# 
# 
# ##output b
# 
########################################################################

#RUN CCDCLEANING

#see and run limpia_cosmicos_parall.py (require lacosmic from Van Dokkum and parallel python)

imagen_ing='rrg'+imagen_in[:-5]+'_tmp_CC2.fits'

 
#Quantum efficiency correction  to science (QE correction):

iraf.gqecorr(inimages=imagen_ing[:-5],outpref='q_',refimages='ARC_'+ARCg[:-5],qecorr_data='gmosQEfactors.dat')  
 

#TOTAL EXTRACTION!!!

FLATg='rrg'+FLAT[:-5]+'_tmp.fits'
ARCg='rrg'+ARC[:-5]+'_tmp.fits'

iraf.gfreduce(inimages='q_'+imagen_ing,outpref='SC_',fl_inter='no',fl_vardq='no',fl_addmdf='no',fl_over='no',fl_trim='no',fl_gscrrej='no',fl_extract='yes',fl_gsappwave='yes',fl_wavtran='yes',fl_skysub='no',trace='no',recenter='no',bias=SBias,fl_bias='yes',orde='1',weight="none",fl_fluxcal='no',ref='FLAT_q_'+FLATg,response='RESP_'+FLATg,wavtraname='ARC_'+ARCg[:-5],fl_fulldq='yes') 

#####
# Sky subtraction needs to be done. Not in this pipeline, since needs carefully checking....

