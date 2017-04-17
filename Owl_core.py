#######################OWL CORE###############################
#####    Set of Functions used in Owl 2.0                #####
##############################################################

#F.Comitani 2016

import shlex, subprocess, sys, os, math, random

def readFile(nameFile):
    outList = []
    data =open(nameFile,"r")
    for line in data.readlines():
        outList.append([])
        for i in line.split():
            outList[-1].append(i)
    data.close()
    return outList

def splitDataset(dataset, n):
    x = []
    for j in range(len(dataset)):
        for k in range(len(dataset[j])):
            if k==n:
                x.append(float(dataset[j][k]))
    return x

############################################################################

def average_histo(h,w):
	summa=0
	zeros=0
	for f in range(len(h)):
		summa=summa+w[f]*float(h[f])
		if w[f]==0:
			zeros=zeros+1.0
	return summa/(len(h)-zeros)

def correlation(h1,h2,w):
	sigma=0
	sigma2=0
	av1 = average_histo(h1,w)
	av2 = average_histo(h2,w)
	for f in range(len(h1)):
		sigma=sigma+w[f]*(h1[f]-av1)*(h2[f]-av2)
		sigma2=sigma2+w[f]*((h1[f]-av1)**2)*((h2[f]-av2)**2)
	return sigma/math.sqrt(sigma2)

def chisquare(h1,h2,w):
	sigma=0
	for f in range(len(h1)):
		if w[f]!=0 and (h1[f]+h2[f])!=0:
			sigma=sigma+w[f]*((h1[f]-h2[f])**2)/(h1[f]+h2[f])
	return sigma

def intersection(h1,h2,w):
	sigma=0
	for f in range(len(h1)):
		sigma=sigma+w[f]*numpy.min(h1[f],h2[f])
	return sigma

def bhattacharyya(h1,h2,w):
	sigma=0
	zeros=0
	for f in range(len(h1)):
		sigma=sigma+w[f]*math.sqrt(h1[f]*h2[f])
		if w[f]==0:
			zeros=zeros+1.0
 	return math.sqrt(1-sigma/math.sqrt(average_histo(h1,w)*average_histo(h2,w)*(len(h2)-zeros)**2))

def spearman(h1,h2,w):
	sigma=0
	sigma1=0
	sigma2=0
	av1 = average_histo(h1,w)
	av2 = average_histo(h2,w)
	for f in range(len(h1)):
		if w[f]!=0:
			sigma+=w[f]*(h1[f]-av1)*(h2[f]-av2)
			sigma1+=(h1[f]-av1)**2
			sigma2+=(h2[f]-av2)**2
	if sigma1*sigma2==0:
		return 0
	else:
		return sigma/math.sqrt(sigma1*sigma2)

####################################################################################
def importDatabase(nameFolder):
	print 'Importing the Ingredients Database...'
	#Create a dictionary with a list of respective frequencies for each molecule
	#Take the log files
	fullIngredients=[]
	os.chdir(nameFolder+"/")
	for files in os.listdir("."):
    		if files.endswith(".log"):   
			fullIngredients.append(files)
	#Define a list with the numbers of the molecules
	nameIngredients=[]
	for i in fullIngredients:
		end=i.find("_") #Change it according to the names of the input files, or change the name of the files directly!
		if end == -1:
			end=i.find(".")
		nameIngredients.append(i[0:end])
	print 'Done!\n'
	print 'Assigning Frequencies to each Ingredient...'
	database={}
	for i in nameIngredients:
		database[i]=[]
	#Fill the dictionary
	for j in range(len(fullIngredients)):
		dataset = readFile(fullIngredients[j])
		for k in range(len(dataset)):
			if dataset[k] != []:
				if dataset[k][0]=='Frequencies':
					database[nameIngredients[j]].append(float(dataset[k][2]))
					database[nameIngredients[j]].append(float(dataset[k][3]))
					database[nameIngredients[j]].append(float(dataset[k][4]))
				if len(dataset[k])>1 and dataset[k][1]=='Thermochemistry':
					break
	os.chdir("../")	
	print 'Done!\n'
	return database, fullIngredients, nameIngredients

def importTraining(nameTrainingFolder, database):
	print 'Importing the Training Recipe...'
	#As before, the last entry of the dictionary is "training" to which all the frequencies of the training list are linked
	fullTraining=[]
	os.chdir(nameTrainingFolder+"/")
	for files in os.listdir("."):
    		if files.endswith(".log"):   
			fullTraining.append(files)
	nameTraining=[]
	for i in fullTraining:
		end=i.find("_") #Change it according to the names of the input files, or change the name of the files directly!
		if end == -1:
			end=i.find(".")
		nameTraining.append(i[0:end])
	database["training"]=[]
	for j in range(len(fullTraining)):
		datasetT = readFile(fullTraining[j])
		for k in range(len(datasetT)):
			if datasetT[k] != []:
				if datasetT[k][0]=='Frequencies':
					database["training"].append(float(datasetT[k][2]))
					database["training"].append(float(datasetT[k][3]))
					database["training"].append(float(datasetT[k][4]))
				if len(datasetT[k])>1 and datasetT[k][1]=='Thermochemistry':
						break
	os.chdir("../")	
	print 'Done!\n'
	return database, fullTraining, nameTraining

def fillHistogram(mini, maxi, nbin, database, ingredientsList):
	bin=(maxi-mini)/nbin
	histo = []
	histoNorm = []
	allFreq=[]
	for x in range(0,nbin):
		histo.append(0.0)
	for w in ingredientsList:
		for h in range(len(database[w])):
			paolinoX=float(database[w][h])	
			x=int((paolinoX-mini)/bin)
                	histo[x]= histo[x]+1.0   
			allFreq.append(database[w][h])
	for h in histo:
		histoNorm.append(h/(float(len(allFreq))*bin))   
	return histo, histoNorm, allFreq

def printHistogram(mini, maxi, nbin, histogram, name):
	bin=(maxi-mini)/nbin
	#Print it to file in a Gnuplot format
	file=open("Histo"+name+".plot",'w')
	for x in range(0,nbin):
		file.write(str(mini+bin*x) + '\t' + str(0)+ '\n')
		file.write(str(mini+bin*x) + '\t' + str(histogram[x])+ '\n')
		file.write(str(mini+bin*(x+1)) + '\t' + str(histogram[x])+ '\n')
	file.close()

def printSpline(mini, maxi, nbin, histogram, name):
	bin=(maxi-mini)/nbin
	#Print it to file in a Gnuplot format
	file=open("HistoSpline"+name+".plot",'w')
	for x in range(0,nbin):
		file.write(str(mini+bin*x) + '\t' + str(histogram[x])+ '\n')
	file.close()

def importWeights(namefile, minFreq, numBins):
		  mini=minFreq
		  nbin=numBins
		  weights=[]
		  weightsTemp = readFile('./Bin_weights.list')
		  for x in range(0,nbin):
			if ((x+0.5)*bin+mini)> float(weightsTemp[-1][1]):
				weights.append(0.0)
			else:
				for j in range(len(weightsTemp)):
					if ((x+0.5)*bin+mini)> float(weightsTemp[j][0]) and ((x+0.5)*bin+mini) <= float(weightsTemp[j][1]): 	
						weights.append(float(weightsTemp[j][2]))
		  return weights

def addition(histoNew, histoNormNew, mini, maxi, nbin, allFreq, database, choice):
	bin=(maxi-mini)/nbin
	#Initialise histograms etc.
	histoNewTemp=[]
	histoNormNewTemp=[]
	allFreqTemp=[]
	errTemp=0
	for y in range(len(histoNew)):	
		histoNewTemp.append(histoNew[y])
		histoNormNewTemp.append(histoNormNew[y])
	for z in allFreq:
		allFreqTemp.append(z)
	#Update the histograms 
	for f in range(len(database[choice])):
		allFreqTemp.append(database[choice][f])
                paolinoX=float(database[choice][f])
                x=int((paolinoX-mini)/bin)
                histoNewTemp[x]= histoNewTemp[x]+1.0
	for h in range(len(histoNewTemp)):
                histoNormNewTemp[h]=histoNewTemp[h]/(float(len(allFreqTemp))*bin)
		#diff[h]=histoNormNewTemp[h]-histoNorm[h]
		#diff[h]=histoNewTemp[h]-histoRef[h]
	return histoNewTemp, histoNormNewTemp, allFreqTemp

def substitution(histoNew, histoNormNew, mini, maxi, nbin, allFreq, database, choiceSub, choiceAdd):
	bin=(maxi-mini)/nbin
	histoNewTemp=[]
	histoNormNewTemp=[]
	allFreqTemp = []
	errTemp=0
	for y in range(len(histoNew)):	
		histoNewTemp.append(histoNew[y])
		histoNormNewTemp.append(histoNormNew[y])
	for z in allFreq:
		allFreqTemp.append(z)
	for f in range(len(database[choiceSub])):
	 	allFreqTemp.remove(database[choiceSub][f])
	        paolinoX=float(database[choiceSub][f])
       	        x=int((paolinoX-mini)/bin)
       		histoNewTemp[x]= histoNewTemp[x]-1.0
	for f in range(len(database[choiceAdd])):
	 	allFreqTemp.append(database[choiceAdd][f])
	        paolinoX=float(database[choiceAdd][f])
	        x=int((paolinoX-mini)/bin)
       		histoNewTemp[x]= histoNewTemp[x]+1.0
	for h in range(len(histoNewTemp)):
       	        histoNormNewTemp[h]=histoNewTemp[h]/(float(len(allFreqTemp))*bin)
		#diff[h]=histoNormNewTemp[h]-histoNorm[h]
		#diff[h]=histoNewTemp[h]-histoRef[h]
	return histoNewTemp, histoNormNewTemp, allFreqTemp	

		
def recipeSearch(limitValue, mini, maxi, nbin, histoRefNorm, database, nameIngredients, nameTraining):
	bin=(maxi-mini)/nbin
	histoNew = []
	histoNormNew = []
	print 'Searching through Ingredients for an Alternative Recipe...'
	#Define a new empty histogram and an error function as sum of the weighted sqared differences of 
	#this histo normalised respect to the reference normalised one looking at each bin
	for x in range(0,nbin):
		histoNew.append(0.0)
		histoNormNew.append(0.0)
		#diff.append(histoNormNew[x]-histoNorm[x])	
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
	error=chisquare(histoNormNew,histoRefNorm,w)
	#Cyclically pick a random molecule, add it and see if the error diminishes or not, if it does accept the move,
	#if it doesn't accept it with probability exp(-increase/limit)
	############################### WHILE loop #####################################################
	limit=limitValue
	added=[]
	allFreq=[]
	nameIngUp=[]
	errPrint=[]
	for n in nameIngredients:
		nameIngUp.append(n)
	found=True
	while (error>limit):
		#Control 
		if len(nameIngUp)==0:
			print "\nWARNING: Correct solution not found: increase the error limiting value!"
			found=False
			break
		########################## Addition ############################################################
		#Pick random element
		pick=int(math.floor((len(nameIngUp))*random.random()))
        	choice=nameIngUp[pick]	
        	histoNewTemp, histoNormNewTemp, allFreqTemp = addition(histoNew, histoNormNew, mini, maxi, nbin, allFreq, database, choice)
		#Recalculate the error
		errTemp=chisquare(histoNormNewTemp,histoRefNorm,w)
		#Accept or refuse the move and the update (or not!) histograms and error
		accept=False
		increase=errTemp-error
		test=random.random()
		exp=math.exp(-increase/limit)
		if test<=exp:
			accept=True
		if accept:
			added.append(choice)
			for y in range(len(histoNewTemp)):
        	        	histoNew[y]=histoNewTemp[y]
        	        	histoNormNew[y]=histoNormNewTemp[y]
			allFreq = []
			for z in allFreqTemp:
				allFreq.append(z)
			error=errTemp
			errPrint.append(error)
		nameIngUp.remove(choice)
		########################## Substitution #####################################################
		test2=random.random()
		if test2<1.0 and len(added)>0 and len(nameIngUp)>0:
			nameIngSub=[]
			for n in nameIngredients:
				if n not in added:
        				nameIngSub.append(n)
			#Pick random element from added and from the remaining ones
			pick2=int(math.floor((len(added))*random.random()))
			pick3=int(math.floor((len(nameIngSub))*random.random()))
			choiceSub=added[pick2]
			choiceAdd=nameIngSub[pick3]
			histoNewTemp, HistoNormNewTemp, allFreqTemp = substitution(histoNew, histoNormNew, mini, maxi, nbin, allFreq, database, choiceSub, choiceAdd)		
			#Ricalculate the error		
			errTemp=chisquare(histoNormNewTemp,histoRefNorm,w)
			increase2=errTemp-error
			accept2=False
			test3=random.random()
			exp2=math.exp(-increase2/limit)
			if test3<=exp2:
				accept2=True
			if accept2:
				added.remove(choiceSub)
				added.append(choiceAdd)
				for y in range(len(histoNewTemp)):
        	        		histoNew[y]=histoNewTemp[y]
        	        		histoNormNew[y]=histoNormNewTemp[y]
				allFreq = []
				for z in allFreqTemp:
					allFreq.append(z)
				error=errTemp
				errPrint.append(error)
				nameIngUp.append(choiceSub)
				if choiceAdd in nameIngUp:
					nameIngUp.remove(choiceAdd)
        	#If found a Recipe equal to the starting one, then restart
		counter = 0
		for i in added:
			if i in nameTraining:
				counter+=1
		if error == 0.0 and counter==len(nameTraining):
			print "Oops, I need to start again!"
			added=[]
			allFreq=[]
			nameIngUp=[]
			for n in nameIngredients:
       				nameIngUp.append(n)
			histoNew = []
			histoNormNew = []
			for x in range(0,nbin):
        	       		histoNew.append(0.0)
        	       		histoNormNew.append(0.0)
				error=chisquare(histoNormNew,histoRefNorm,w)			
				errPrint=[]
	####End of the WHILE loop##########################################################################
	print "Done!\n"
	if found:
		print "New Recipe Found!"
	else:
		print "This is the Best I can do!"
	print "Number of elements: " + str(len(added))+ ";"
	print   "\nError: "+str(error)
	return histoNew, histoNormNew, added, error, errPrint, found

def printOutputHistograms(histoNew, histoNewNorm, name): 
	#Print the new histogram normalised to file in a Gnuplot format
	file=open(name+"Norm.plot",'w')
	for x in range(0,nbin):
		file.write(str(mini+bin*x) + '\t' + str(0)+ '\n')
		file.write(str(mini+bin*x) + '\t' + str(histoNormNew[x])+ '\n')
		file.write(str(mini+bin*(x+1)) + '\t' + str(histoNormNew[x])+ '\n')
	file.close()
	#Print the new histogram to file in a Gnuplot format
	file=open(name+".plot",'w')
	for x in range(0,nbin):
		file.write(str(mini+bin*x) + '\t' + str(0)+ '\n')
		file.write(str(mini+bin*x) + '\t' + str(histoNormNew[x]*len(allFreq)*bin)+ '\n')
		file.write(str(mini+bin*(x+1)) + '\t' + str(histoNormNew[x]*len(allFreq)*bin)+ '\n')
	file.close()

def printIngredients(added, name):
	#Print the Ingredients to file
	file=open(name+".list",'w')
	for x in added:
		file.write(x+'\n')
	file.close()

def printErrorStepwise(errPrint, name):
	#Print the error stepwise
	file=open(name+".plot",'w')
	for x in range(len(errPrint)):
		file.write(str(x) + '\t' + str(errPrint[x])+ '\n')
	file.close()

def comparison(added):
	#Now the comparison with 60A and 60B
	#Take the log files
	a60=[]
	a60Full=[]
	os.chdir("Recipe60A/")
	for files in os.listdir("."):
	    if files.endswith(".log"):
	        a60Full.append(files)
	#Define a list with the numbers of the molecules
	for i in a60Full:
                end=i.find("_") #Change it according to the names of the input files, or change the name of the files directly!
                if end == -1:
                        end=i.find(".")
                a60.append(i[0:end])
	b60=[]
	b60Full=[]
	os.chdir("../Recipe60B/")
	for files in os.listdir("."):
	    if files.endswith(".log"):
	        b60Full.append(files)
	#Define a list with the numbers of the molecules
	for i in b60Full:
                end=i.find("_") #Change it according to the names of the input files, or change the name of the files directly!
                if end == -1:
                        end=i.find(".")
                b60.append(i[0:end])
	b30=[]
	b30Full=[]
	os.chdir("../Recipe30B/")
	for files in os.listdir("."):
	    if files.endswith(".log"):
	        b30Full.append(files)
	#Define a list with the numbers of the molecules
	for i in b30Full:
                end=i.find("_") #Change it according to the names of the input files, or change the name of the files directly!
                if end == -1:
                        end=i.find(".")
                b30.append(i[0:end])
	os.chdir("../")
	#compare the input list with those of the whites
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
	print	"Similarities: "+str(round(ratioA,2))+"% to White_60A, "+str(round(ratioB,2))+"% to White_60B and "+str(round(ratioC,2))+"% to White_30B.\n"


