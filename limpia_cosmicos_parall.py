#!/usr/env python

####################
#
# Uses L.A. Cosmic from P. van Dokkum 
# Requires python version of la_cosmic
# http://www.astro.yale.edu/dokkum/lacosmic/
# it also uses Parallel python PP
#
#
##################


import sys, os, string

import astropy
from astropy.io import fits
import numpy as np
import pyfits
import scipy
from scipy import ndimage

# IMPORTING La Cosmic
import cosmics
import time

# Paralel python
import pp



if len(sys.argv) <=1 :
    print "\nExecute the program as follows:\n"
    print "python limpia_cosmicos_parall.py imagen_Gemini.fits Ncpus"
    print "where"
    print "imagen_Gemini.fits is gprepared with ONLY fl_addmdf=yes"
    print "Ncpus is the number of CPUs of the computer\n"
    exit()
else:
    print "Loading the image "

imagen_in = sys.argv[1]



def separa_Amps(imagen_a_corregir):
    #start_time = time.time()
    #header = pyfits.getheader(imagen_a_corregir)
    # Creating arrays
    imagen = []
    header_amp = []
    Amp = []

    for i in range(1, 13):
        print i
        # CC Cosmic Cleanned. Amp= Amplifiyer N (NOT CCed)
        header_amp.append(pyfits.getheader(imagen_a_corregir, 'sci', i))
        imagen.append(pyfits.getdata(imagen_a_corregir, 'sci', i))
        Amp.append(str(imagen_a_corregir[:-5] + "_Amp_" + str(i) + ".fits"))

    for k in range(0, 12):
        print "writing", Amp[k]
        pyfits.writeto(Amp[k], imagen[k], header_amp[k])

    return


def limpia_cosmicos(imagen_a_corregir, N_Amp):
    # imagen_a_corregir=imagen_in
    import cosmics
    # Amp.append(str(imagen_a_corregir[:-5]+"_Amp_"+str(i)+".fits"))
    # cleaned.append(str(imagen_a_corregir[:-5]+"_CC_"+str(i)+".fits"))

    array, header = cosmics.fromfits(str(imagen_a_corregir[:-5] + "_Amp_" + str(N_Amp) + ".fits"))
    # array, header = cosmics.fromfits(str(imagen_in[:-5]+"_Amp_"+str(N_Amp)+".fits"))
    gain_h = header["GAIN"]  # READING GAIN AND NOISE FROM HEADERS
    rnoise = header["RDNOISE"]
    print "READING GAIN AND NOISE ", gain_h, rnoise, str(imagen_a_corregir[:-5] + "_Amp_" + str(N_Amp) + ".fits")
    # print "READING GAIN AND NOISE ", gain_h,rnoise,str(imagen_in[:-5]+"_Amp_"+str(N_Amp)+".fits")

    # FOR IFU SCIENCE, NO STD STARS
    # c = cosmics.cosmicsimage(array, gain=gain_h, readnoise=rnoise, sigclip =3, sigfrac = 1, objlim = 5.0)
    # c = cosmics.cosmicsimage(array, gain=gain_h, readnoise=rnoise, sigclip=0.9, sigfrac=1, objlim=2.0)

    #FOR STD
    c = cosmics.cosmicsimage(array, gain=gain_h, readnoise=rnoise, sigclip=3, sigfrac=1, objlim=2.0)

    # There are other options, check the manual...

    # Run the full artillery 7 iterations is enough:
    c.run(maxiter=7)

    # Write the cleaned image into a new FITS file, conserving the original header :
    cleaned = str(imagen_a_corregir[:-5] + "_CC2_" + str(N_Amp) + ".fits")
    cosmics.tofits(cleaned, c.cleanarray, header)

    # yield (none)
    # If you want the mask, here it is :
    #    cosmics.tofits("mask.fits", c.mask, header)
    #


print "Separando los Amplificadores"
print "separa_Amps..."

start_time1 = time.time()

separa_Amps(imagen_in)

start_time2 = time.time()


print "7 iterations for cosmic rejection..."



ppservers = ()
#ppservers = ("10.0.0.1",)

if len(sys.argv) > 1:
    ncpus = int(sys.argv[2])
    # Creates jobserver with ncpus workers
    job_server = pp.Server(ncpus, ppservers=ppservers)
else:
    # Creates jobserver with automatically detected number of workers
    job_server = pp.Server(ppservers=ppservers)#

print "Starting pp with", job_server.get_ncpus(), "workers aka CPUs"


#FUERZA BRUTA....
Amplific= (1,2,3,4,5,6,7,8,9,10,11,12)
jobs=[(input,job_server.submit(limpia_cosmicos,(imagen_in,input,))) for input in Amplific]
for input, job in jobs:
    print "executing job N", input, "is", job()

print "Time elapsed parallel: ", time.time() - start_time2, "s"
job_server.print_stats()


print "Total  time elapsed: ", time.time() - start_time1, "s"
#job_server.print_stats()



###############################
#
# Rearmando la imagen en formato
# Gemini:
#
###############################

print "###############################"
print "#"
print "# Rearmando la imagen en formato"
print "# Gemini:"
print "#"
print "###############################"



imagen_end=imagen_in[:-5]
imagen_a_corregir=imagen_in
FILE_FIN=imagen_end+"_CC2.fits"    #imagen final
FILE_TEMP=imagen_end+"_tmp.fits"    #imagen Temporal donde se hacen los cambios
imagen_header_original=imagen_in   # Imagen Original con todos los headers,

import shutil

#   making a copy of the priginal file to a temporal one.
shutil.copyfile(imagen_header_original,FILE_TEMP)
imagen_header2=FILE_TEMP

for i in range (1,13):

    imagen_Corre=str(imagen_a_corregir[:-5]+"_CC2_"+str(i)+".fits")  # Output of correction "CC", NAME_CC.fits
    imagen_no_header=pyfits.getdata(imagen_Corre,header=False)
    header=pyfits.getheader(FILE_TEMP)
    header2=pyfits.getheader(FILE_TEMP,i)
    pyfits.update(imagen_header2,imagen_no_header,ext=i) #copia la imagen en si a la extension correspondiente
    pyfits.update(imagen_header2,imagen_no_header,header2,i) #copia el header original a la extension correspondiente


import shutil
shutil.copyfile(FILE_TEMP,FILE_FIN)##

##Final Sanity check

#pyfits.info(FILE_FIN)
#print "Deleting files"

print "Deleting _CC_ files and _Amp_ files"

for m in range (1,13):
    print " rm "+imagen_a_corregir[:-5]+"_CC2_"+str(m)+".fits"
    os.remove(imagen_a_corregir[:-5]+"_CC2_"+str(m)+".fits")
    print " rm "+imagen_a_corregir[:-5]+"_Amp_"+str(m)+".fits"
    os.remove(imagen_a_corregir[:-5]+"_Amp_"+str(m)+".fits")


print "Cosmic cleaning done"
print "Total execution time: %.1f minutes" % ((time.time() - start_time1)/60.)

print "\n Marcelo D. Mora, 07.09.2016 mmora@astro.puc.cl v.2.0 \n"

