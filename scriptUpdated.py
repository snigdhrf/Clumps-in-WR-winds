## TO DO
# Plot flux-mean, on top of each other, time as y axis
# Variability study! (see a papers (for example: Wind inhomogeneities in WR stars 1, 1996), 2009 St-Louis is for WN star)
# NIII higher ionisation, formed deeper in the wind, than HeII

# Read of velocity, and get acceleration by comparing

# Greyscale instead of colorscale

# delta v of a clump (by its width of the line) is related to mass/density

# Do normalisation on part of the spectrum, around the line


import numpy as np
from numpy import *
from pylab import *
import matplotlib.pyplot as plt
import glob
import pyfits as pf
import scipy.optimize
#from astropy.io import fits
from operator import add
import datetime as dt
import scipy.ndimage.filters


#Gives plots of
	# the normalized fluxes and their mean
	# the residuals 
	# the standard deviatons together with the scaled mean normalized flux spectrum 

#LaTeX#
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


print('------------------START:', dt.datetime.now(),'---------------------')
print('Data is loading, please wait')

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

def smoothen(values,sigma):
	#Sigma chosen random, since trying to set sigma=standard_deviations did not work
	
	fit=scipy.ndimage.filters.gaussian_filter(values, sigma)
	return fit



def getRes(data,mean_flux):
	resi=[]
	for i in range(len(data)):
		resi.append(data[i]-mean_flux[i])

	#Uncomment if you dont want the residuals to be smoothened!!
	res=smoothen(resi,5)

	return(res)


def calc_mean_flux(night,number_of_datapoints):
	#Calculation of the mean of the normalized fluxes
	mf = []
	mf_scaled=[]
	mean_fp=[]
	mean_sp=[]
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


		mf.append(mean)
		mean_fp.append(mean1)
		mean_sp.append(mean2)
		mf_scaled.append(mean/10)

	mean_flux=smoothen(mf,10)
	mean_first_part=smoothen(mean_fp,10)
	mean_second_part=smoothen(mean_sp,10)
	mean_flux_scaled=smoothen(mf_scaled,10)



	return(mean_flux, mean_first_part, mean_second_part, mean_flux_scaled)



##Def to smoothen a part of the spectrum , now never used
def smoothen_line(flux,i1,i2):
	x=waves[i1:i2]
	flux_tr=np.transpose(flux)
	fluxSel_tr=flux_tr[i1:i2]
	fluxSel= np.transpose(fluxSel_tr)
	
	sigma=10

	fit=scipy.ndimage.filters.gaussian_filter(fluxSel, sigma)
	plt.plot(x,fit,label='fit')
	plt.plot(x,fluxSel, label= 'smoothend flux')
	plt.legend()
	plt.show()
	return(fit)




def plotResTimeseries(res,night):
	#Get residiual plot thing
	##HeII 
	matrix=[]
	index1= find_index(4675,waves)
	index2= find_index(4695,waves)


	for i in range(len(night)):
		matrix.append(res[i][index1:index2])
	plt.imshow(matrix,cmap=plt.cm.BuPu_r)
	plt.title('Timeseries of the residiuals for HeII line',fontsize=16)
	plt.colorbar()
	plt.show()

	##NIII
	matrix=[]
	index1= find_index(4622,waves)
	index2= find_index(4648,waves)


	for i in range(len(night)):
		matrix.append(res[i][index1:index2])
	plt.imshow(matrix,cmap=plt.cm.BuPu_r)
	plt.title('Timeseries of the residiuals for NIII line',fontsize=16)
	plt.colorbar()
	plt.show()


	##NIV
	matrix=[]
	index1= find_index(4044,waves)
	index2= find_index(4065,waves)

	for i in range(len(night)):
		matrix.append(res[i][index1:index2])
	plt.imshow(matrix,cmap=plt.cm.BuPu_r)
	plt.title('Timeseries of the residiuals for NIV line',fontsize=16)
	plt.colorbar()
	plt.show()


##Load in all files

files = glob.glob('ourData' + '/*normalized_MA.txt')
files.sort()
number_of_files = len(files)


flux=[]
wavelengths=[]
for i in range(number_of_files): 
	data= np.loadtxt(files[i], skiprows=0, unpack = True)
	flux.append(data[1][483:].tolist())
	wavelengths.append(data[0][483:].tolist())

# ##Get data and wavelengths
# flux=[]
# wavelengths=[]
# for i in range(number_of_files):
# 	#fluxes
# 	f=getData(files[i])
# 	flux.append(f)
# 	#wavelengths
# 	w = getWavel(files[i])
# 	wavelengths.append(w)


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


###Test extra normalisation by dividing all spectra by their maximum value
# def extra_normalisation():
# 	for i in range(len(flux)):
# 		flux[i]=flux[i]/(max(flux[i]))
# 
# extra_normalisation()



night1=flux[0:11]
night2=flux[11:26]
night3=flux[26:36]
night4a=flux[36:44]
night4b=flux[44:]

nights=[night1,night2,night3,night4a,night4b]

def all_plots(night): 

	print('----------start plotting----------')

	mean_flux, mean_first_part, mean_second_part, mean_flux_scaled= calc_mean_flux(night,number_of_datapoints)	

	

	def plotNormFl():
		#Plot of the normalized fluxes and their mean
		fig=plt.figure()
		plt.plot(waves,mean_flux,label = 'mean')
		for i in range(len(night)):
			plt.plot(waves,night[i], label=i)
		#plt.legend()

		plt.xlabel('wavelength', fontsize=16)
		plt.title('The normalized fluxes',fontsize=16)
		plt.show()

	
	# plotNormFl()


	def compare():
		plt.plot(waves,mean_flux,label = 'mean')
		plt.plot(waves,mean_first_part, label='first part')
		plt.plot(waves,mean_second_part, label='second part')
		plt.title('Comparison mean of first and second part of timeseries',fontsize=16)
		plt.legend()
		plt.show()

	# compare()

	#Calculation of the residuals of the normalized fluxes
	print('Calculating the residuals')

	res=[]
	for i in range(len(night)):
		res.append(getRes(night[i],mean_flux))



	def plotRes(res):
		#Plot of the residuals 
		fig,axs = plt.subplots(2)
		for i in range(len(night)):
			axs[0].plot(waves,res[i], label=i+1)
		
		axs[1].plot(waves, mean_flux, 'b' ,label = 'Mean Flux')
		axs[1].set_ylabel('Mean Flux',fontsize=12)
		axs[0].set_ylabel('residuals',fontsize=12)
		axs[1].set_xlabel('wavelength',fontsize=12)

		#plt.legend()
		#ax[1].ylim(-1,6)
		plt.show()

	plotRes(res)	
		
	plotResTimeseries(res,night)


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

		fig,axs = plt.subplots(2)
		axs[0].plot(waves,standardDeviation,'g',label = 'std. deviation')
		axs[1].plot(waves, mean_flux, 'b' ,label = 'Mean Flux')
		axs[1].set_xlabel(r'Wavelength ([$\AA$])',fontsize=16)
		axs[0].set_ylabel(r'Stan. deviation',fontsize = 12)
		axs[1].set_ylabel(r'Mean normalized Flux',fontsize = 12)
		#plt.xlim(4000,6740)
		plt.show()
	plotStd()



# 
# 	def plotResByGauss():

# 		plotResNIII(night)
# 			
# 		plotResNIV(night)		

# 		plotResHeII(night)

# 	plotResByGauss()


	print('End of this night')

print('Data has loaded, start plotting:', dt.datetime.now())
 

def doPlotting(nightNumber):
	if nightNumber=='1':
		all_plots(night1)

	elif nightNumber=='2':
		all_plots(night2)

	elif nightNumber=='3':
		all_plots(night3)

	elif nightNumber=='4a':
		all_plots(night4a)

	elif nightNumber=='4b':
		all_plots(night4b)

	else: 
		plot('Sorry, you gave wrong input')



nightNumber = input('For which night do you want to see the plots? If you want to quit, press q. Otherwise choose between 1, 2, 3, 4a and 4b:  ')

while (nightNumber!='q'):
	doPlotting(nightNumber)
	nightNumber = input('For which night do you want to see the plots? If you want to quit, press q. Otherwise choose between 1, 2, 3, 4a and 4b:   ')



#mean nights
mf1,mf_first_part_1, mf_second_part_1, mf_scaled_1=calc_mean_flux(night1,number_of_datapoints)
mf2,mf_first_part_2, mf_second_part_2, mf_scaled_2=calc_mean_flux(night2,number_of_datapoints)
mf3,mf_first_part_3, mf_second_part_3, mf_scaled_3=calc_mean_flux(night3,number_of_datapoints)
mf4a,mf_first_part_4a, mf_second_part_4a, mf_scaled_4a=calc_mean_flux(night4a,number_of_datapoints)
mf4b,mf_first_part_4b, mf_second_part_4b, mf_scaled_4b=calc_mean_flux(night4b, number_of_datapoints)



def doPlotting2(diff_nights):
	if diff_nights=='a':
		fig=plt.figure()
		plt.plot(waves,mf2,label = 'mean night2')
		plt.plot(waves,mf3,label = 'mean night3')
		plt.plot(waves,mf4a,label = 'mean night 4 part 1')
		plt.plot(waves,mf4b,label = 'mean night 4 part 2')
		plt.legend()
		plt.xlabel('wavelength', fontsize=16)
		plt.title('The mean fluxes of different nights', fontsize=16)
		plt.show()


	elif diff_nights=='b':
		fig=plt.figure()
		plt.plot(waves,mf4a,label = 'mean night 4 part 1')
		plt.plot(waves,mf4b,label = 'mean night 4 part 2')
		plt.legend()
		plt.xlabel('wavelength', fontsize=16)
		plt.title('The mean fluxes of night 4', fontsize=16)
		plt.show()

	else:
		print('Sorry you gave wrong input')


diff_nights= input('To compare night 2,3 and 4? press a. To see the plot to compare night 4 part a and b, press b. To quit, press q.     ')
while(diff_nights!= 'q'):
	doPlotting2(diff_nights)
	diff_nights= input('To compare night 2,3 and 4? Press a. To see the plot to compare night 4 part a and b, press b. To quit, press q.     ')





def calc_all_res(index1,index2):
	matrix=[]
	res=getRes(night2,mf2)
	for i in range(len(night2)):
	 	matrix.append(res[i][index1:index2])
 
	res=getRes(night3,mf3)
	for i in range(len(night3)):
	 	matrix.append(res[i][index1:index2])

	res=getRes(night4a,mf4a)
	for i in range(len(night4a)):
		matrix.append(res[i][index1:index2])

	res=getRes(night4b,mf4b)
	for i in range(len(night4b)):
		matrix.append(res[i][index1:index2])

	return(matrix)
	

def resiual_plotting(matrix):

	##NIV
	index1= find_index(4044,waves)
	index2= find_index(4065,waves)
	
	matrix=plot_all_res(index1,index2)
	
	cm = plt.cm.get_cmap('nipy_spectral') #the choise of the used colormap
	plt.imshow(matrix,cmap=plt.cm.BuPu_r)
	plt.title('Timeseries of the residiuals for NIV line',fontsize=16)
	plt.colorbar()
	plt.show()
	
	
	##HeII
	index1= find_index(4675,waves)
	index2= find_index(4695,waves)
	
	
	cm = plt.cm.get_cmap('nipy_spectral') #the choise of the used colormap
	plt.imshow(matrix,cmap=plt.cm.BuPu_r)
	plt.title('Timeseries of the residiuals for HeII line',fontsize=16)
	plt.colorbar()
	plt.show()
	
	
	##NIII
	index1= find_index(4622,waves)
	index2= find_index(4648,waves)
	
	matrix=plot_all_res(index1,index2)
	
	cm = plt.cm.get_cmap('nipy_spectral') #the choise of the used colormap
	plt.imshow(matrix,cmap=plt.cm.BuPu_r)
	plt.title('Timeseries of the residiuals for NIII line',fontsize=16)
	plt.colorbar()
	plt.show()
	
pr =input('Do you want to see the plots of the residuals for all nights together? y/n:   ')

if pr=='y':
	matrix=calc_all_res(index1,index2)
	resiual_plotting(matrix)

	
	
## TRASH
	
# Definitions to fit Gauss trough line to find residuals:
## Not a good idea, not a Gaussian profile

# def fitGaussian(index1,index2,flux):
# 
# 	x=np.array(waves[index1:index2])
# 	y=np.array(flux[index1:index2])
# 
# 	mean=sum(x*y)/sum(y)
# 	sigma=sum(y*(x-mean)**2)/sum(y)
# 
# 
# 	def gauss_function(x,a,x0,sigma):
# 		return a*np.exp(-(x-x0)**2/(2*sigma**2))
# 
# 	#p0 is first guess
# 
# 	popt,pcov= scipy.optimize.curve_fit(gauss_function, x, y, p0 = [1, mean, sigma])
# 	fit= gauss_function(x, *popt)
# 	res= abs(y-fit)/y
# 
# 	#plt.plot(x,fit, label='Gaussian fit')
# 	#plt.plot(x,y,label='original data')
# 	#plt.plot(x,res, label='residuals')
# 
# 	#plt.legend()
# 	#plt.show()
# 	return(res)	
# 
# 
# def plotResNIII(night):
# 	##NIII line
# 	resNIII=[]
# 	index1= find_index(4610,waves)
# 	index2= find_index(4665,waves)
# 	x=np.array(waves[index1:index2])
# 	for i in range(len(night)):
# 		resNIII.append(fitGaussian(index1,index2,night[i]))
# 	for i in range(len(night)):
# 		j=i+1
# 		plt.plot(x,resNIII[i],label=j)
# 
# 	plt.legend()
# 	plt.title('residualsNIII',fontsize=20)
# 	plt.show()
# 
# #NIV line
# def plotResNIV(night):
# 	resNIV=[]
# 	index1=find_index(4040,waves)
# 	index2=find_index(4075,waves)
# 	x=np.array(waves[index1:index2])
# 	for i in range(len(night)):
# 		resNIV.append(fitGaussian(index1,index2,night[i]))
# 	for i in range(len(night)):
# 		j=i+1
# 		plt.plot(x,resNIV[i],label=j)
# 
# 	plt.legend()	
# 	plt.title('residualsNIV', fontsize=20)
# 	plt.show()
# 
# 
# 
# #HeII line
# def plotResHeII(night):
# 	resHeII=[]	
# 	index1=find_index(4665,waves)
# 	index2=find_index(4710,waves)
# 	x=np.array(waves[index1:index2])
# 
# 	for i in range(len(night)):
# 		resHeII.append(fitGaussian(index1,index2,night[i]))
# 
# 	for i in range(len(night)):
# 		j=i+1
# 		plt.plot(x,resHeII[i],label=j)
# 
# 	plt.legend()
# 	plt.title('residualsHeII', fontsize=20)
# 	plt.show()


## Plot of residuals using Gaussian thing, not good idea
# def resNight4():
# 	spectraChosen=[]
# 	spectraChosen.append(night4a[0])
# 	spectraChosen.append(night4a[len(night4a)-1])
# 	spectraChosen.append(night4b[0])
# 	spectraChosen.append(night4b[len(night4b)-1])
# 
# 	plotResNIV(spectraChosen)
# 	plotResNIII(spectraChosen)
# 	plotResHeII(spectraChosen)
# 
# resNight4()


#all_plots(night1)
#all_plots(night2)
#all_plots(night3)
#all_plots(night4a)
#all_plots(night4b)

# def different_nights():
# 	#compare mean nights
# 	mf1,mf_first_part_1, mf_second_part_1, mf_scaled_1=calc_mean_flux(night1,number_of_datapoints)
# 	mf2,mf_first_part_2, mf_second_part_2, mf_scaled_2=calc_mean_flux(night2,number_of_datapoints)
# 	mf3,mf_first_part_3, mf_second_part_3, mf_scaled_3=calc_mean_flux(night3,number_of_datapoints)
# 	mf4a,mf_first_part_4a, mf_second_part_4a, mf_scaled_4a=calc_mean_flux(night4a,number_of_datapoints)
# 	mf4b,mf_first_part_4b, mf_second_part_4b, mf_scaled_4b=calc_mean_flux(night4b, number_of_datapoints)
# 
# 	fig=plt.figure()
# 	plt.plot(waves,mf2,label = 'mean night2')
# 	plt.plot(waves,mf3,label = 'mean night3')
# 	plt.plot(waves,mf4a,label = 'mean night 4 part 1')
# 	plt.plot(waves,mf4b,label = 'mean night 4 part 2')
# 	plt.legend()
# 	plt.xlabel('wavelength')
# 	plt.title('The mean fluxes of different nights')
# 	plt.show()
# 
	##plot of first and last spectrum of night 4
	# plt.plot(waves,night4a[0],label='first')
	# plt.plot(waves,night4b[len(night4b)-1],label='last')
	# plt.legend()
	# plt.xlabel('wavelength')
	# plt.title('The first and last spectrum of night 4')
	# plt.show()
	
#different_nights()