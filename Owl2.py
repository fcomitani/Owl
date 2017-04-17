#! /usr/bin/env python

##############################OWL 2.0#################################################
######################Olfactory White aLgorithm#######################################
##### Program to Elaborate White Odorants Recipies of Molecules as function      #####
#####  of the Vibrational Spectrum extracted with Gaussian09.                    #####
######################################################################################

#F.Comitani 2016

#import shlex, subprocess, sys, os, math, random
from Owl_core import *

print '#######################################################'
print '################ Welcome to OWL 2.0! ##################'
print '############# Olfactory White aLgorithm! ##############'
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
printHistogram(mini, maxi, nbin, histoRef, 'Ref')
printHistogram(mini, maxi, nbin, histoRefNorm, 'RefNorm')

limit = 0.00190994735444
#Chisquare no Weights ~ 0.00190994735444
#Chisquare Weights ~ 0.000912924325706
#Use Assaggiatore.py to evaluate this limit
histoNew, histoNewNorm, added, error, errPrint, found = recipeSearch(limit, mini, maxi, nbin, histoRefNorm, database, nameIngredients, nameTraining)
printIngredients(added, 'NewIngredients')
printErrorStepwise(errPrint, 'ErrorStepwise')
printHistogram(mini, maxi, nbin, histoNew, 'New')
printHistogram(mini, maxi, nbin, histoNewNorm, 'NewNorm')

#Comparison with the other white recipes
comparison(added)

######################################################################################
#########                            THE END                               ###########
######################################################################################


