import numpy as np
from numpy import *
from pylab import *
import matplotlib.pyplot as plt
import glob
import pyfits as pf
from operator import add


def getData(direct):
	data = genfromtxt(direct,skip_header=0,usecols=(1),dtype=[('flux','float64')])
	fluxes=[]
	for i in range(len(data)):
		fluxes.append(data[i][0])
	return(fluxes[483:])

def getWavel(direct):
	wavel= genfromtxt(direct,skip_header=0,usecols=(0),dtype=[('wavelength','float64')])
	#Something weird, our datapoints are given as number, io just number, so remove the comma
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
		res.append((data[i]-mean_flux[i]))
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

for i in range(number_of_files):
	if len(wavelengths[i])> number_of_datapoints :
		wavelengths[i]=wavelengths[i][:number_of_datapoints]
		flux[i]=flux[i][:number_of_datapoints]
	


night1=flux[0:11]
night2=flux[11:26]
night3=flux[26:36]
night4a=flux[36:44]
night4b=flux[44:]


def all_spectra_of_night(night):
	for i in range(len(night)):
		for j in range(number_of_datapoints):
			night[i][j]= night[i][j] + i

	def plot_spectra():
		fig = plt.figure()
		for i in range(len(night)):
			plt.plot(waves,night[i], label=i)
		plt.legend()
		plt.ylabel('Relative intensity')
		plt.xlabel('wavelength')
		plt.title('Spectral Time series')
		plt.show()
	plot_spectra()

#all_spectra_of_night(night1)

################################ MY CODE     ############################################3

#they contain wavelength intervals
HeII4686 = []
HeII5411 = []
NIV4058  = []

index1 = find_index(4630,waves)
index2 = find_index(4740,waves)

index3 = find_index(5300,waves)
index4 = find_index(5500,waves)

index5 = find_index(4020,waves)
index6 = find_index(4070,waves)


def get_line(line,index_One,index_two):
	for i in range(index_One,index_two):
		line.append(waves[i])

get_line(HeII4686,index1,index2)
get_line(HeII5411,index3,index4)
get_line(NIV4058,index5,index6)

#For night1, night2, night3
def get_residuals_per_line(night,line):
	residuals = []
	mean_flux, mean_first_part, mean_second_part, mean_flux_scaled= calc_mean_flux(night,number_of_datapoints)	

	for i in range(len(night)):
		index1 = find_index(line[0],waves)
		index2 = find_index(line[-1],waves)
		res= []
		for j in range(index1,index2+1):
			res.append(abs(night[i][j]-mean_flux[j]))
		residuals.append(res)

	for i in range(len(night)):
		index1 = find_index(line[0],waves)
		index2 = find_index(line[-1],waves)
		each_residual_length = len(residuals[i])
		for j in range(each_residual_length):
			residuals[i][j]= residuals[i][j] + i*0.3

	def res_fig():
		fig=plt.figure()
		plt.axvline(4058)
		for i in range(len(night)):
			plt.plot(line,residuals[i])
		plt.legend()
		plt.xlabel('wavelength')
		plt.ylabel('Residuals + Offset')
		plt.title('Absolute Residuals')
		plt.show()

	res_fig()

#get_residuals_per_line(night1,NIV4058)
#get_residuals_per_line(night2,NIV4058)
#get_residuals_per_line(night3,NIV4058)



## Final night
def get_residuals_per_line_final_night(line):
	residuals = []
	mean_flux4a, mean_first_part4a, mean_second_part4a, mean_flux_scaled4a = calc_mean_flux(night4a,number_of_datapoints)	
	mean_flux4b, mean_first_part4b, mean_second_part4b, mean_flux_scaled4b = calc_mean_flux(night4b,number_of_datapoints)	

	for i in range(len(night4a)):
		index1 = find_index(line[0],waves)
		index2 = find_index(line[-1],waves)
		res= []
		for j in range(index1,index2+1):
			res.append(abs(night4a[i][j]-mean_flux4a[j]))
		residuals.append(res)

	for i in range(len(night4b)):
		index3 = find_index(line[0],waves)
		index4 = find_index(line[-1],waves)
		res2= []
		for j in range(index3,index4+1):
			res2.append(abs(night4b[i][j]-mean_flux4b[j]))
		residuals.append(res2)


	for i in range(len(residuals)):
		each_residual_length = len(residuals[i])
		for j in range(each_residual_length):
			residuals[i][j]= residuals[i][j] + i*0.3

	def res_fig():
		fig=plt.figure()
		plt.axvline(4058)
		for i in range(len(residuals)):
			plt.plot(line,residuals[i])
		plt.legend()
		plt.xlabel('wavelength')
		plt.ylabel('Residuals + Offest')
		plt.title('Absolute Residuals Final Night')
		plt.show()

	res_fig()


#get_residuals_per_line_final_night(NIV4058)

#################################################################################################################3


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
		plt.title('The normalized fluxes for the Final Night')
		plt.legend()
		plt.show()

	compare()

	#Calculation of the residuals of the normalized fluxes
	res=[]
	for i in range(len(night)):
		res.append(getRes(night[i],mean_flux))

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

	#plotRes()	
		

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

#different_nights()
