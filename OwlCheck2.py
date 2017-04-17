#! /usr/bin/env python

#################################OWLCHECK2############################################
##### Program to Check Errors and Difference between recipies                    #####
#####  obtained with Owl 2                                                       #####
######################################################################################

#F.Comitani 2015

#import shlex, subprocess, sys, os, math, random
from Owl_core import *

print '#######################################################'
print '############## Welcome to OWLCHECK 2.0! ###############'
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


print 'Importing the New Recipe...'

Ing=[]
if os.path.isfile('./NewIngredients.list'):
	  IngTemp = readFile('./NewIngredients.list')
	  for n in range(len(IngTemp)):
		Ing.append(IngTemp[n][0])
	  histoNew = []
	  histoNewNorm = []
	  print 'Searching through Ingredients for an Alternative Recipe...'
	  #Define a new empty histogram and an error function as sum of the weighted sqared differences of 
	  #this histo normalised respect to the reference normalised one looking at each bin
	  for x in range(0,nbin):
		histoNew.append(0.0)
		histoNewNorm.append(0.0)
		#diff.append(histoNewNorm[x]-histoNorm[x])	
	  #Import the weights from a file if there is a file, otherwise all the frequencies below 2000 are weighted 1, those over 2000 are ignored.
	  if os.path.isfile('./Bin_weights.list'):
			w=importWeights('./Bin_weights.list')
	  else:
		print '\nWARNING:No Weights.list file found. No weights will be applied.\nRun OwlWeights.py if you want weights to be applied!'
		w=[]
 	  	for x in range(0,nbin):
			if x>nbin/2: ############### Weights OFF ###############
				w.append(1.0)
			else:
				w.append(1.0)
else:
	  print '\nERROR:No list of new ingredients found (I need the NewIngredients.list file)!'
          quit()

allFreq=[]
for n in Ing:
	histoNew, histoNewNorm, allFreq = addition(histoNew, histoNewNorm, mini, maxi, nbin, allFreq, database, n)
error=chisquare(histoNewNorm,histoRefNorm,w)

print "Done!\n"

#Now the comparison with 60A and 60B
comparison(Ing)

print   "Error from Training Recipe: "+str(error)+"\n"
 

######################################################################################
#########                            THE END                               ###########
######################################################################################

