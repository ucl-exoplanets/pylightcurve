import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import pymc
import fcmodel
import time
##########################################################################################




##########################################################################################
def model_fit(data,iter,burn,ldcoeffs,RP,RPvar,P,A,Avar,E,I,Ivar,W,T0,WW=0):
	if isinstance(data, basestring):
		datax, datay = np.loadtxt(data,usecols=[0,1],unpack=True)
	else:
		datax, datay = data
	directory = 'mcmc_fit'
	if os.path.isdir(directory)==False:
		os.mkdir(directory)
		os.chdir(directory)
	else:
		i = 2
		while os.path.isdir(directory + '_' + str(i) ):
			i = i + 1
		directory = directory + '_' + str(i)
		os.mkdir(directory)
		os.chdir(directory)
	a1, a2, a3, a4 = ldcoeffs
	T0	=	T0 + P*int((datax[-1]-T0)/P)
	rp	=	pymc.Uniform( 'rp', (1-RPvar)*RP, (1+RPvar)*RP      , value = RP )
	a	=	pymc.Uniform( 'a' , (1-Avar)*A  , (1+Avar)*A        , value = A  )
	i	=	pymc.Uniform( 'i' , (1-Ivar)*I  , min(90,(1+Ivar)*I), value = I  )
	t0	=	pymc.Uniform( 't0', -0.05       , 0.05              , value = 0  )
	@pymc.deterministic
	def model_mcmc(rp=rp, a=a, i=i, t0=t0):
		return fcmodel.model((a1,a2,a3,a4),rp,P,a,E,i,W,T0+t0,datax,WW)
	y = pymc.Normal('y', mu = model_mcmc, tau = 10**8, observed = True, value = datay)
	count0 = time.time()
	M = pymc.MCMC([rp, a, i, t0, model_mcmc, y], db = 'pickle', dbname = 'model.pickle')
	M.db
	M.isample(iter, burn = burn, verbose = 1)
	M.db.close()
	count1 = time.time()
	print 'Sampling finished in {0} mins ( {1} points, {2} iterations ) \n'.format((count1-count0)/60,len(datax),iter)
	rpchain = M.trace('rp')[:]
	achain  = M.trace('a' )[:]
	ichain  = M.trace('i' )[:]
	t0chain = M.trace('t0')[:]
	pymc.Matplot.plot(M)
	plt.close('all')
	## save results
	result = [rpchain,achain,ichain,t0chain]
	w=open('result.txt','w')
	w.write('param \tinitial \tfinal \terror\n')
	w.write('RPORS\t{}\t{}\t{}\n'.format(RP,np.mean(result[0]),np.std(result[0])))
	mtp			=	T0 + np.mean(result[3])
	w.write('MTP\t{}\t{:.15}\t{}\n'.format(T0,mtp,np.std(result[3])))
	w.write('INCL\t{}\t{}\t{}\n'.format(I,np.mean(result[2]),np.std(result[2])))
	w.write('SMAORS\t{}\t{}\t{}\n'.format(A,np.mean(result[1]),np.std(result[1])))
	w.write('LDA1\t{}\n'.format(a1))
	w.write('LDA2\t{}\n'.format(a2))
	w.write('LDA3\t{}\n'.format(a3))
	w.write('LDA4\t{}\n'.format(a4))	
	w.close()
	phase = (datax-mtp)/P
	model_lc = fcmodel.model((a1,a2,a3,a4),np.mean(result[0]),P,np.mean(result[1]),E,np.mean(result[2]),W,mtp,datax,WW)
	model_rs = datay-model_lc
	np.savetxt('model_lc.txt',np.swapaxes([datax,phase,datay,model_lc,model_rs],0,1))
	np.savetxt('resultlists.txt',np.swapaxes(result,0,1))
	## plot
	print 'Plotting model' 
	plt.subplot2grid((4,1), (0,0), rowspan=3)
	plt.plot(phase,datay,'bo',mec='b')
	datax2 = np.arange(datax[0],datax[-1],(datax[1]-datax[0])/100)
	phase2 = (datax2-mtp)/P
	model_lc2 = fcmodel.model((a1,a2,a3,a4),np.mean(result[0]),P,np.mean(result[1]),E,np.mean(result[2]),W,mtp,datax2,WW)
	plt.plot(phase2,model_lc2,'r-')
	xlim = np.max(np.abs([plt.xlim()[0],plt.xlim()[1]]))
	plt.xlim(-xlim,xlim)
	plt.text(	plt.xlim()[0]+0.02*(plt.xlim()[-1]-plt.xlim()[0]),
				plt.ylim()[0]+0.02*(plt.ylim()[-1]-plt.ylim()[0]),
				r'$rms_{res} = %.3e$' % np.std(model_rs)
			)
	plt.yticks(plt.yticks()[0][1:])
	plt.ylabel(r'$relative \, flux$')
	plt.tick_params(labelbottom='off')
	plt.subplot(4,1,4)
	plt.axhline(0,color='k')
	plt.plot(phase,model_rs,'ko',ms=3)
	plt.xlim(-xlim,xlim)
	plt.xlabel(r'$phase$')
	plt.ylabel(r'$residuals$')
	plt.subplots_adjust(hspace=0.0)
	plt.savefig('fit.png',dpi=200)
	os.system("rm model.pickle")
	os.chdir('..')
##########################################################################################

