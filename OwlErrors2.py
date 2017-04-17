#! /usr/bin/env python

###############################OWLERRORS 2############################################
#####           Program to Evaluate the error among the training recipes         #####
#####                                                                            #####
######################################################################################

#F.Comitani 2015

from Owl_core import *

print '#######################################################'
print '############ Welcome to OWLERRORS 2.0! ################'
print '#######################################################\n'

#Histogram parameters
maxi=4000.0
mini=0.0
nbin=200
os.chdir('../')

#Import training recipe1
print 'First Stage...'
database1={}
database1, fullTraining1, nameTraining1 = importTraining('Recipe30B', database1)
#Fill the reference histogram
histoRef1, histoRefNorm1, allFreqRef1 = fillHistogram(mini, maxi, nbin, database1, ['training'])
printHistogram(mini, maxi, nbin, histoRef1, 'Ref30B')
printHistogram(mini, maxi, nbin, histoRefNorm1, 'Ref30BNorm')
printSpline(mini, maxi, nbin, histoRefNorm1, 'Ref30BNorm')

#Import training recipe2
print 'Second Stage...'
database2={}
database2, fullTraining2, nameTraining2 = importTraining('Recipe60A', database2)
#Fill the reference histogram
histoRef2, histoRefNorm2, allFreqRef2 = fillHistogram(mini, maxi, nbin, database2, ['training'])
printHistogram(mini, maxi, nbin, histoRef2, 'Ref60A')
printHistogram(mini, maxi, nbin, histoRefNorm2, 'Ref60ANorm')
printSpline(mini, maxi, nbin, histoRefNorm2, 'Ref60ANorm')

#Import training recipe3
print 'Third Stage...'
database3={}
database3, fullTraining3, nameTraining3 = importTraining('Recipe60B', database3)
#Fill the reference histogram
histoRef3, histoRefNorm3, allFreqRef3 = fillHistogram(mini, maxi, nbin, database3, ['training'])
printHistogram(mini, maxi, nbin, histoRef3, 'Ref60B')
printHistogram(mini, maxi, nbin, histoRefNorm3, 'Ref60BNorm')
printSpline(mini, maxi, nbin, histoRefNorm3, 'Ref60BNorm')

print 'Evaluating Errors...'

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

error12 = chisquare(histoRefNorm1,histoRefNorm2,w)
error23 = chisquare(histoRefNorm2,histoRefNorm3,w)
error13 = chisquare(histoRefNorm1,histoRefNorm3,w)

print "\nError 30b-60a: " + str(error12) 
print "Error 60a-60b: " + str(error23) 
print "Error 30b-60b: " + str(error13) 

######################################################################################
#########                            THE END                               ###########
######################################################################################


