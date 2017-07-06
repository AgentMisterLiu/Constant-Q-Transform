'''
	Constant-Q Transform with FFT optimisation
	Developer: Zhenning Liu
'''
import random
import numpy as np
import math
from matplotlib import pyplot as plt

def w(nk,n):
	if (n>nk):
		return 0
	else:
		return 1

def CQT(ks,x,b,fs,f0):
	'''
	ks: integer, total number of filters
	x: numpy array, data in time series
	b: float, number of filters per octave
	fs: float, sampling rate
	f0: floar, fundamental frequency

	return powerspec: 2-d numpy array
	'''
	current=0
	powerspec=np.zeros((len(x),ks),dtype=complex)
	q=1.0/(2**(1.0/b)-1)
	n=int(q*fs/f0)
	y=np.zeros((ks,n),dtype=complex)
	yft=np.zeros((n,ks),dtype=complex)
	lenx=len(x)
	for k in range(ks):
		fk=f0*(2**(k/b))
		nk=nk=int(q*fs/fk)
		for i in range(n):
			y[k,i]=1.0/float(nk)*w(nk,i)*np.exp(2j*np.pi*fk*i/fs)
	for k in range(ks):
		ykft=np.fft.fft(y[k,:])	
		ykft.shape=(len(ykft),1)
		yft[:,k]=ykft.transpose()
	ystar=yft.conj()
	while (current+n) <= lenx :
		#print current,n,len(x)
		xn=x[current:current+n]
		xft=np.fft.fft(xn)
		xcq=np.dot(xft,ystar)
		xcq.shape=(len(xcq),1)
		powerspec[current]=xcq.transpose()
		current+=1
	return powerspec[0:current,:]

if __name__ == '__main__':
	ks=10
	x = np.random.rand(100000)
	b=3.0
	fs=1000.0
	f0=50.0
	ps=CQT(ks,x,b,fs,f0)
	print ps

