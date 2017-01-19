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

#Running BIAS and OVERSCAN PARA CORREGIR DEFECTOS EN TODAS LAS IMAGENES

#iraf.gfreduce(inimages=imagen_in,outpref='default',fl_inter='no',fl_vardq='yes',fl_addmdf='yes',fl_over='yes',fl_trim='no',fl_gscrrej='no',fl_wavtran='no',fl_skysub='no',trace='no',recenter='no',bias=SBias,fl_bias='yes',fl_gsappwave='no',fl_fluxcal='no',fl_fulldq='no',fl_fixgap='no',line=100,nsum=10)

#iraf.gfreduce(inimages=ARC,outpref='default',fl_inter='no',fl_vardq='yes',fl_addmdf='yes',fl_over='yes',fl_trim='no',fl_gscrrej='no',fl_wavtran='no',fl_skysub='no',trace='no',recenter='no',bias=SBias,fl_bias='yes',fl_gsappwave='no',fl_fluxcal='no',fl_fulldq='no',fl_fixgap='no',line=100,nsum=10)

#iraf.gfreduce(inimages=FLAT,outpref='default',fl_inter='no',fl_vardq='yes',fl_addmdf='yes',fl_over='yes',fl_trim='no',fl_gscrrej='no',fl_wavtran='no',fl_skysub='no',trace='no',recenter='no',bias=SBias,fl_bias='yes',fl_gsappwave='no',fl_fluxcal='no',fl_fulldq='no',fl_fixgap='no',line=100,nsum=10)

 
#run script german!

####FLAT for reference
#
#   iraf.gfreduce(inimages =                
#   (outimages = "")          
#     (outpref = "default")   
#       (slits = "header")    
#     (exslits = "*")         
#(fl_nodshuffl = no)          
#    (fl_inter = yes)         
#    (fl_vardq = no)          
#   (fl_addmdf = yes)         
#     (fl_over = yes)         
#     (fl_trim = yes)         
#     (fl_bias = yes)         
#  (fl_scatsub = no)          
#   (fl_qecorr = no)          
#  (fl_gscrrej = yes)         
#   (fl_crspec = no)          
#(fl_gnsskysub = no)          
#  (fl_extract = yes)         #
#(fl_gsappwave = yes)        # 
#  (fl_wavtran = yes)         
#   (fl_skysub = yes)         
#  (fl_fluxcal = yes)         
#   (fl_fulldq = no)          
#    (dqthresh = 0.1)         
#     (rawpath = "")          
#     (key_mdf = "")          
#     (mdffile = "default")   
#      (mdfdir = "gmos$data/")
# (key_biassec = "BIASSEC")   
# (key_datasec = "DATASEC")   
#     (bpmfile = "gmos$data/chipgaps.dat") Bad pixel mask for column interpolation
#        (grow = 1.0)            
#        (bias = "")             
#   (reference = "")             
#   (sc_xorder = "1")            
#   (sc_yorder = "3")            
#    (sc_cross = no)             
#    (qe_refim = "")             
#(fl_keep_qeim = no)             
# (qe_corrpref = "qecorr")       
#(qe_corrimage = "")             
#(qe_data = "gmosQEfactors.dat") 
#  (qe_datadir = "gmos$data/")   
#    (response = "")             
#  (wavtraname = "")             
#   (sfunction = "")             
#  (extinction = "")             
#    (fl_fixnc = no)             
#  (fl_fixgaps = no)             
#   (fl_novlap = yes)            
#    (perovlap = 10.0)           
# (nbiascontam = "default")      
#    (biasrows = "default")      
#       (order = "default")      
#  (low_reject = 3.0)            
# (high_reject = 3.0)            
#    (niterate = 2)              
#   (cr_xorder = 9)              
#  (cr_sigclip = 4.5)            
#  (cr_sigfrac = 0.5)            
#   (cr_objlim = 1.0)            
#    (cr_niter = 4)              
#        (line = 100)            
#        (nsum = 10)             
#       (trace = yes)            
#    (recenter = yes)            
#      (thresh = 200.0)          
#    (function = "chebyshev")    
#     (t_order = 5)              
#      (t_nsum = 10)             
#     (weights = "variance")     
#   (gratingdb = "gmos$data/GMOSgratings.dat") Gratings database file
#    (filterdb = "gmos$data/GMOSfilters.dat") Filters database file
#     (xoffset = INDEF)          
#        (expr = "default")      
#    (sepslits = no)             
#          (w1 = INDEF)          
#          (w2 = INDEF)          
#          (dw = INDEF)          
#          (nw = INDEF)          
# (observatory = "default")      
#     (sci_ext = "SCI")          
#     (var_ext = "VAR")          
#      (dq_ext = "DQ")           
#     (logfile = "")             
#     (verbose = yes)            
#      (status = 0)              
#    (scanfile = )               
#        (mode = "al")           #



###############
# ##
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

#iraf.gfreduce(inimages=ARCg,outpref='ARC_',fl_inter='no',fl_vardq='no',fl_addmdf='no',fl_over='no',fl_trim='no',fl_gscrrej='no',fl_wavtran='no',fl_skysub='no',trace='no',recenter='no',fl_extract='yes',bias=SBias,fl_bias='no',orde='1',weight="none",fl_gsappwave='no',fl_fluxcal='no',ref='REF_FLAT_'+FLATg,fl_fulldq='yes')


#iraf.gswavelength(inimages='ARC_'+ARCg,coordlist='CuAr_GMOS.dat',fl_inter='yes')
 
##Apply QE to the FLATS
#################
#
#####iraf.gqecorr(inimages='b_rg'+FLAT,refimages='ARC_'+ARC[:-5],qecorr_data="gmosQEfactors.dat") 
#
#######################
FLATg='rrg'+FLAT[:-5]+'_tmp.fits'

#iraf.gqecorr(inimages=FLATg,outpref='q_',refimages='ARC_'+ARCg[:-5],qecorr_data="gmosQEfactors.dat") 


##############################################################################
#  ###Repitiendo la extraccion con el flat corregido por scattered qb_rg...
#
#  ####iraf.gfreduce(inimages='qtest_rg'+FLAT,outpref='FLAT_',fl_inter='no',fl_vardq='no',fl_addmdf='no',fl_over='no',fl_trim='no',fl_gscrrej='no',fl_wavtran='no',fl_skysub='no',trace='no',recenter='no',fl_extract='yes',bias=SBias,fl_bias='yes',fl_gsappwave='yes',fl_fluxcal='no',ref='REF_FLAT_'+FLAT,fl_fulldq='yes',fl_fixgap='no')
#################################################################################

#iraf.gfreduce(inimages='q_'+FLATg,outpref='FLAT_',fl_inter='no',fl_vardq='no',fl_addmdf='no',fl_over='no',fl_trim='no',fl_gscrrej='no',fl_wavtran='no',fl_skysub='no',trace='no',recenter='no',fl_extract='yes',bias=SBias,fl_bias='yes',fl_gsappwave='yes',fl_fluxcal='no',ref='REF_FLAT_'+FLATg,fl_fulldq='yes',fl_fixgap='no')

###haciendo el responso... normalizando el flat

####################################################
##  
##  iraf.gfresponse(inimage='FLAT_qb_rg'+FLAT,outimage='RESP_'+FLAT,skyimage="",fl_inter='yes',function='spline3',order='45',sample="*",wavtraname='ARC_'+ARC[:-5])
##  iraf.gfresponse(inimage='FLAT_q_rg'+FLAT,outimage='RESP_'+FLAT,skyimage="",fl_inter='yes',function='spline3',order='45',sample="*",wavtraname='ARC_'+ARC[:-5]) 
##
##########################################################

#OK

#iraf.gfresponse(inimage='FLAT_q_'+FLATg,outimage='RESP_'+FLATg,skyimage="",fl_inter='yes',function='spline3',order='45',sample="*") 
#iraf.gfresponse(inimage='FLAT_q_'+FLATg,outimage='RESP_'+FLATg,skyimage="",fl_inter='yes',function='spline3',order='45',sample="*",wavtraname='ARC_'+ARCg[:-5]) 



###############################################################
#
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
##correr cleanning de cosmicos.... aqui A ver como se agrega al script...

##Corregir los defectos del los CCDs

imagen_ing='rrg'+imagen_in[:-5]+'_tmp_CC2.fits'

 
##correr QE correction a la sciencia:

#iraf.gqecorr(inimages=imagen_ing[:-5],outpref='q_',refimages='ARC_'+ARCg[:-5],qecorr_data='gmosQEfactors.dat')  
 

#EXTRACCION TOTAL!!!

FLATg='rrg'+FLAT[:-5]+'_tmp.fits'
ARCg='rrg'+ARC[:-5]+'_tmp.fits'

#iraf.gfreduce(inimages='q_'+imagen_ing,outpref='SC_',fl_inter='no',fl_vardq='no',fl_addmdf='no',fl_over='no',fl_trim='no',fl_gscrrej='no',fl_extract='yes',fl_gsappwave='yes',fl_wavtran='yes',fl_skysub='no',trace='no',recenter='no',bias=SBias,fl_bias='yes',orde='1',weight="none",fl_fluxcal='no',ref='FLAT_q_'+FLATg,response='RESP_'+FLATg,wavtraname='ARC_'+ARCg[:-5],fl_fulldq='yes') 





##iraf.gfreduce(inimages='q_brg'+imagen_in,outpref='SC_',fl_inter='no',fl_vardq='no',fl_addmdf='no',fl_over='no',fl_trim='no',fl_gscrrej='no',fl_extract='yes',fl_gsappwave='yes',fl_wavtran='yes',fl_skysub='no',trace='no',recenter='no',bias=SBias,fl_bias='yes',orde='1',weight="none",fl_fluxcal='no',ref='FLAT_qb_rg'+FLAT,response='RESP_test'+FLAT,wavtraname='ARC_'+ARC[:-5],fl_fulldq='yes') 


##iraf.gfreduce(inimages='qbrg'+imagen_in,outpref='SC_',fl_inter='no',fl_vardq='no',fl_addmdf='no',fl_over='no',fl_trim='no',fl_gscrrej='no',fl_extract='yes',fl_gsappwave='yes',fl_wavtran='yes',fl_skysub='no',trace='no',recenter='no',bias=SBias,fl_bias='yes',orde='1',weight="none",fl_fluxcal='no',ref='FLAT_qb_rg'+FLAT,response='RESPONSE_TEST.fits',wavtraname='ARC_'+ARC[:-5],fl_fulldq='yes') 

##iraf.gfreduce(inimages=imagen_in,outpref='SC_',fl_inter='no',fl_vardq='no',fl_addmdf='no',fl_over='no',fl_trim='no',fl_gscrrej='no',fl_extract='yes',fl_gsappwave='yes',fl_wavtran='yes',fl_skysub='no',trace='no',recenter='no',bias=SBias,fl_bias='yes',orde='1',weight="none",fl_fluxcal='no',ref='FLAT_qb_rg'+FLAT,response='RESPONSE_TEST.fits',wavtraname='ARC_'+ARC[:-5],fl_fulldq='yes') 

##Sky subtraction:
