import numpy as np
from numpy import *
from pylab import *
import matplotlib.pyplot as plt
import glob
import pyfits as pf
#from astropy.io import fits
from operator import add

#Gives plots of
	# the normalized fluxes and their mean
	# the residuals 
	# the standard deviatons together with the scaled mean normalized flux spectrum 

## Functions

def getData(direct):
	data = genfromtxt(direct,skip_header=0,usecols=(1),dtype=[('flux','float64')])
	fluxes=[]
	for i in range(len(data)):
		fluxes.append(data[i][0])
	return(fluxes[483:])

def getWavel(direct):
	wavel= genfromtxt(direct,skip_header=0,usecols=(0),dtype=[('wavelength','float64')])
	#Something weird, our datapoints are given as number, io just number, so remove the komma
	waves=[]
	for i in range(len(wavel)):
		waves.append(wavel[i][0])
	return (waves[483:])


def find_index(wavelength,waves):
	index=0
	for i in range(len(waves)):
		if waves[i]<wavelength:
			index=index+1
		else:
			break
	return(index)


def getRes(data,mean_flux):
	res=[]
	for i in range(len(data)):
		res.append((data[i]-mean_flux[i])/mean_flux[i])
	return(res)


def calc_mean_flux(night,number_of_datapoints,):
	#Calculation of the mean of the normalized fluxes
	mean_flux = []
	mean_flux_scaled=[]
	mean_first_part=[]
	mean_second_part=[]
	half=int(len(night)/2)

	for i in range(number_of_datapoints):
		mean=0
		mean1=0
		mean2=0
		for j in range(len(night)):
			mean=mean+night[j][i]
		mean=mean/len(night)

		for j in range(half):
			mean1=mean1+night[j][i]
		mean1=mean1/half	

		for j in range(half,len(night)):
			mean2=mean2+night[j][i]
		number_of_elm=len(night)-half
		mean2=mean2/number_of_elm


		mean_flux.append(mean)
		mean_first_part.append(mean1)
		mean_second_part.append(mean2)
		mean_flux_scaled.append(mean/10)

	return(mean_flux, mean_first_part, mean_second_part, mean_flux_scaled)

##Load in all files

files = glob.glob('ourData' + '/*normalized_MA.txt')
files.sort()
number_of_files = len(files)


##Get data and wavelengths
flux=[]
wavelengths=[]
for i in range(number_of_files):
	#fluxes
	f=getData(files[i])
	flux.append(f)
	#wavelengths
	w = getWavel(files[i])
	wavelengths.append(w)


waves=wavelengths[0]
number_of_datapoints=len(waves)

## Make sure all datasets have the same amount of datapoints (can be done with interpolate, but didnt work)
for i in range(number_of_files):
	if len(wavelengths[i])> number_of_datapoints :
		wavelengths[i]=wavelengths[i][:number_of_datapoints]
		flux[i]=flux[i][:number_of_datapoints]
	


#Interpolate to make sure for every dataset we have the same amount of datapoints: ERROR
#flux=[]
#for i in range(1,number_of_files):
#	flux[i]=np.interp(waves,wavelengths[i],data[i])


##Test extra normalisation by dividing all spectra by their maximum value
def extra_normalisation():
	for i in range(len(flux)):
		flux[i]=flux[i]/(max(flux[i]))

#extra_normalisation()


night1=flux[0:11]
night2=flux[11:26]
night3=flux[26:36]
night4a=flux[36:44]
night4b=flux[44:]


def all_plots(night): 

	print('start plots of one night')

	mean_flux, mean_first_part, mean_second_part, mean_flux_scaled= calc_mean_flux(night,number_of_datapoints)	

	def plotNormFl():
		#Plot of the normalized fluxes and their mean
		fig=plt.figure()
		plt.plot(waves,mean_flux,label = 'mean')
		for i in range(len(night)):
			plt.plot(waves,night[i], label=i)
		plt.legend()
		plt.xlabel('wavelength')
		plt.title('The normalized fluxes')
		plt.show()

	plotNormFl()




	def compare():
		plt.plot(waves,mean_flux,label = 'mean')
		plt.plot(waves,mean_first_part, label='first part')
		plt.plot(waves,mean_second_part, label='second part')
		
		plt.legend()
		plt.show()

	compare()

	#Calculation of the residuals of the normalized fluxes
	res=[]
	for i in range(len(night)):
		res.append(getRes(night[i],mean_flux))



	def plotResTimething():
		#Get residiual plot thing
		matrix=[]
		index1= find_index(4680,waves)
		index2= find_index(4695,waves)
		for i in range(len(night)):
			matrix.append(res[i][index1:index2])


		plt.imshow(matrix,cmap=plt.cm.BuPu_r)
		plt.title('Timeseries of the residiuals')
		plt.colorbar()
		plt.show()

	plotResTimething()



	def plotRes():
		#Plot of the residuals 
		plt.plot(waves,mean_flux,label = 'Mean Fluxes ')

		for i in range(len(night)):
			plt.plot(waves, res[i], label=i)
		
		plt.legend()
		plt.title('Residuals')

		plt.xlabel('wavelength')
		plt.ylim(-1,5)
		plt.show()

	plotRes()	
		

	#Calculation of the standard deviation of the normalized fluxes
	standardDeviation = []
	for i in range(number_of_datapoints):
		std=0
		for j in range(len(night)):
			std=std+(night[j][i]-mean_flux[i])**2
		std=(std/(len(night)-1))**(0.5)
		standardDeviation.append(std)
	

	
	def plotStd():
	#Plot of the standard deviatons together with the scaled mean normalized flux spectrum 

		plt.plot(waves,standardDeviation,label = 'std')
		plt.plot(waves,mean_flux_scaled, label='mean spectrum / 10')
		plt.legend()
		plt.title('Standard deviation of normalized spectra')
		plt.xlabel('wavelength')
		plt.xlim(4000,6740)
		plt.show()

	plotStd()

	print('End of this night')



#all_plots(night1)
#all_plots(night2)
#all_plots(night3)
#all_plots(night4a)
#all_plots(night4b)

def different_nights():
	#compare mean nights
	mf1,mf_first_part_1, mf_second_part_1, mf_scaled_1=calc_mean_flux(night1,number_of_datapoints)
	mf2,mf_first_part_2, mf_second_part_2, mf_scaled_2=calc_mean_flux(night2,number_of_datapoints)
	mf3,mf_first_part_3, mf_second_part_3, mf_scaled_3=calc_mean_flux(night3,number_of_datapoints)
	mf4a,mf_first_part_4a, mf_second_part_4a, mf_scaled_4a=calc_mean_flux(night4a,number_of_datapoints)
	mf4b,mf_first_part_4b, mf_second_part_4b, mf_scaled_4b=calc_mean_flux(night4b, number_of_datapoints)

	fig=plt.figure()
	#plt.plot(waves,mf2,label = 'mean night2')
	#plt.plot(waves,mf3,label = 'mean night3')
	plt.plot(waves,mf4a,label = 'mean night 4 part 1')
	plt.plot(waves,mf4b,label = 'mean night 4 part 2')
	plt.legend()
	plt.xlabel('wavelength')
	plt.title('The mean fluxes of different nights')
	plt.show()

different_nights()
