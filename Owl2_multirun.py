#! /usr/bin/env python

###############################OWL 2.0 Multirun#######################################
##### Multiple runs of Owl 2 and analysis of the resulting recipies              #####
#####  of the Vibrational Spectrum extracted with Gaussian09.                    #####
######################################################################################

#F.Comitani 2015

#import shlex, subprocess, sys, os, math, random
import math
from Owl_core import *

print '#######################################################'
print '############ Welcome to OWL 2.0 Multirun! #############'
print '#######################################################\n'

#Create database
database, fullIngredients, nameIngredients = importDatabase('../Database')

#Import training recipe
database, fullTraining, nameTraining = importTraining('Recipe30B', database)

#Histogram parameters
maxi=4000.0
mini=0.0
nbin=200
#Fill the reference histogram
histoRef, histoRefNorm, allFreqRef = fillHistogram(mini, maxi, nbin, database, ['training'])

limit = 0.00190994735444
#Chisquare no Weights ~ 0.00190994735444
#Chisquare Weights ~ 0.000912924325706
#Use OwlError2.py to evaluate this limit
totals=[]
n=0
while n<1000:
	histoNew, histoNewNorm, added, error, errPrint, found = recipeSearch(limit, mini, maxi, nbin, histoRefNorm, database, nameIngredients, nameTraining)
	if found:	
		for w in added:
			totals.append(w)
		n+=1 
nameIngredients2=[]
for i in nameIngredients:
	if i!='training':
		nameIngredients2.append(int(i))
#Print occurrencies as percentage:
file=open("Occurrencies.list",'w')
for x in sorted(nameIngredients2):
	file.write(str(x)+' '+str(totals.count(str(x))*100.0/1000)+'\n')
file.close()

######################################################################################
#########                            THE END                               ###########
######################################################################################



