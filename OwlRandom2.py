#! /usr/bin/env python

#################################OWLRANDOM 2##########################################
##### Program to generate random results with Owl 2                              #####
##### and evaluate the typical values of the Error                               #####
######################################################################################

#F.Comitani 2015

import numpy
from Owl_core import *

print '#######################################################'
print '############## Welcome to OWLRANDOM 2.0! ##############'
print '#######################################################\n'

#Create database
database, fullIngredients, nameIngredients = importDatabase('../Database')

#Import training recipe
database, fullTraining, nameTraining = importTraining('Recipe30B', database)

#Histogram parameters
maxi=4000.0
mini=0.0
nbin=200
bin=(maxi-mini)/nbin
#Fill the reference histogram
histoRef, histoRefNorm, allFreqRef = fillHistogram(mini, maxi, nbin, database, ['training'])

print 'Building Random Recipes and Evaluating the Error...'

#Now the comparison with 60A and 60B and 30B
a60=[]
a60Full=[]
os.chdir("Recipe60A/")
for files in os.listdir("."):
    if files.endswith(".log"):
        a60Full.append(files)
for i in a60Full:
                end=i.find("_") 
                if end == -1:
                        end=i.find(".")
                a60.append(i[0:end])
b60=[]
b60Full=[]

os.chdir("../Recipe60B/")
for files in os.listdir("."):
    if files.endswith(".log"):
        b60Full.append(files)
for i in b60Full:
                end=i.find("_")
                if end == -1:
                        end=i.find(".")
                b60.append(i[0:end])

b30=[]
b30Full=[]
os.chdir("../Recipe30B/")
for files in os.listdir("."):
    if files.endswith(".log"):
        b30Full.append(files)
for i in b30Full:
                end=i.find("_")
                if end == -1:
                        end=i.find(".")
                b30.append(i[0:end])

os.chdir("../")

#Define weights
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

#Main program
n_histos=100
vecErrors = []
vecDiff60a = []
vecDiff60b = []
vecDiff30b = []
for a in range(0,n_histos):
	histoNew = []
	histoNormNew = []
	added=[]
	allFreq=[]
	nameIngUp=[]
	for x in range(0,nbin):
 			histoNew.append(0.0)
			histoNormNew.append(0.0)	
	randomMol=int(math.floor(random.random()*len(nameIngredients)))
	for n in nameIngredients:
		nameIngUp.append(n)
	for m in range(0,randomMol):
                # Addition
		pick=int(math.floor((len(nameIngUp))*random.random()))
	     	choice=nameIngUp[pick]	
		histoNew, histoNormNew, allFreq = addition(histoNew, histoNormNew, mini, maxi, nbin, allFreq, database, choice)
		nameIngUp.remove(choice)
		added.append(choice)
	error=chisquare(histoRefNorm,histoNormNew,w)
	vecErrors.append(error)
	#Compare the resulting recipe with these Whites
	counterA=0
	counterB=0
	counterC=0
	for i in added:
		if i in a60:
			counterA+=1
		if i in b60:
			counterB+=1
		if i in b30:
			counterC+=1
	unionA=len(added)
	unionB=len(added)
	unionC=len(added)
	for i in a60:
		if i not in added:
			unionA+=1
	for i in b60:
		if i not in added:
			unionB+=1
	for i in b30:
		if i not in added:
			unionC+=1
	ratioA=counterA*100.0/unionA
	ratioB=counterB*100.0/unionB
	ratioC=counterC*100.0/unionC
	vecDiff60a.append(ratioA)
	vecDiff60b.append(ratioB)
	vecDiff30b.append(ratioC)
print "Done!"

#Calculates means and standard deviations
meanE=numpy.mean(vecErrors)
mean60a=numpy.mean(vecDiff60a)
mean60b=numpy.mean(vecDiff60b)
mean30b=numpy.mean(vecDiff30b)
std=numpy.std(vecErrors)
std60a=numpy.std(vecDiff60a)
std60b=numpy.std(vecDiff60b)
std30b=numpy.std(vecDiff30b)

print	"Average Error = " + str(meanE) + " +- " + str(std)
print	"Average similarities: "+str(round(mean60a,2))+" +- "+str(round(std60a,2))+" % to White_60A, "+str(round(mean60b,2))+" +- "+str(round(std60b,2))+" % to White_60B and "+str(round(mean30b,2))+" +- "+str(round(std30b,2))+" % to White_30B."

######################################################################################
#########                            THE END                               ###########
######################################################################################







