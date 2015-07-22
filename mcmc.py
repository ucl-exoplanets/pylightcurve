import numpy as np
import os
import matplotlib.pyplot as plt
import pymc
import models
import time


def transit(data, iterations, burn, binning, ldcoeffs, rprs, rprsvar, p, a, avar, e, i, ivar, w, ww, t0):
    if isinstance(data, basestring):
        datax, datay = np.loadtxt(data, usecols=[0, 1], unpack=True)
    else:
        datax, datay = data
    if binning > 0:
        start = np.mod(len(datax), binning)
        datax = np.mean(np.reshape(datax[start:], (len(datax) / binning, binning)), 1)
        datay = np.mean(np.reshape(datay[start:], (len(datay) / binning, binning)), 1)
    directory = 'mcmc_transit'
    if not os.path.isdir(directory):
        os.mkdir(directory)
        os.chdir(directory)
    else:
        fi = 2
        while os.path.isdir(directory + '_' + str(fi)):
            fi += 1
        directory = directory + '_' + str(fi)
        os.mkdir(directory)
        os.chdir(directory)
    t0 += p * int((datax[-1] - t0) / p)
    fit_rprs = pymc.Uniform('rp', (1 - rprsvar) * rprs, (1 + rprsvar) * rprs, value=rprs)
    fit_a = pymc.Uniform('a', (1 - avar) * a, (1 + avar) * a, value=a)
    fit_i = pymc.Uniform('i', i - ivar, min(90, i + ivar), value=i)
    fit_t0 = pymc.Uniform('t0', -0.05, 0.05, value=0)
    @pymc.deterministic
    def model_mcmc(model_rprs=fit_rprs, model_a=fit_a, model_i=fit_i, model_t0=fit_t0):
        return models.transit(ldcoeffs, model_rprs, p, model_a, e, model_i, w, ww, t0 + model_t0, datax)
    y = pymc.Normal('y', mu=model_mcmc, tau=np.std(datay[1:] - datay[:-1]) ** (-2), observed=True, value=datay)
    count0 = time.time()
    mcmc = pymc.MCMC([fit_rprs, fit_a, fit_i, fit_t0, model_mcmc, y], db='pickle', dbname='model.pickle')
    mcmc.isample(iterations, burn=burn, verbose=1)
    count1 = time.time()
    print 'Sampling finished in {0:.5}mins\n'.format((count1 - count0) / 60)
    pchain = mcmc.trace('rp')[:]
    achain = mcmc.trace('a')[:]
    ichain = mcmc.trace('i')[:]
    t0chain = mcmc.trace('t0')[:]
    pymc.Matplot.plot(mcmc)
    plt.close('all')
    # save results
    result = [pchain, achain, ichain, t0chain]
    www = open('result.txt', 'w')
    www.write('param \tinitial \tfinal \terror\n')
    www.write('RPORS\t{}\t{}\t{}\n'.format(rprs, np.mean(result[0]), np.std(result[0])))
    www.write('MTP\t{}\t{:.15}\t{}\n'.format(t0, t0 + np.mean(result[3]), np.std(result[3])))
    www.write('INCL\t{}\t{}\t{}\n'.format(i, np.mean(result[2]), np.std(result[2])))
    www.write('SMAORS\t{}\t{}\t{}\n'.format(a, np.mean(result[1]), np.std(result[1])))
    rprs = np.mean(result[0])
    a = np.mean(result[1])
    i = np.mean(result[2])
    t0 = t0 + np.mean(result[3])
    www.write('LDA1\t{}\n'.format(ldcoeffs[0]))
    www.write('LDA2\t{}\n'.format(ldcoeffs[1]))
    www.write('LDA3\t{}\n'.format(ldcoeffs[2]))
    www.write('LDA4\t{}\n'.format(ldcoeffs[3]))
    www.close()
    phase = (datax - t0) / p
    model_lc = models.transit(ldcoeffs, rprs, p, a, e, i, w, ww, t0, datax)
    model_rs = datay - model_lc
    np.savetxt('model_lc.txt', np.swapaxes([datax, phase, datay, model_lc, model_rs], 0, 1))
    np.savetxt('resultlists.txt', np.swapaxes(result, 0, 1))
    # plot
    print 'Plotting model' 
    plt.subplot2grid((4, 1), (0, 0), rowspan=3)
    plt.plot(phase, datay, 'bo', mec='b')
    datax2 = np.arange(datax[0], datax[-1], (datax[1] - datax[0]) / 100)
    phase2 = (datax2 - t0) / p
    model_lc2 = models.transit(ldcoeffs, rprs, p, a, e, i, w, ww, t0, datax2)
    plt.plot(phase2, model_lc2, 'r-')
    xlim = np.max(np.abs([plt.xlim()[0], plt.xlim()[1]]))
    plt.xlim(-xlim, xlim)
    plt.text(plt.xlim()[0] + 0.02 * (plt.xlim()[-1] - plt.xlim()[0]),
             plt.ylim()[0] + 0.02 * (plt.ylim()[-1] - plt.ylim()[0]),
             r'$rms_{res} = %.3e$' % np.std(model_rs))
    plt.yticks(plt.yticks()[0][1:])
    plt.ylabel(r'$relative \, flux$')
    plt.tick_params(labelbottom='off')
    plt.subplot(4, 1, 4)
    plt.axhline(0, color='k')
    plt.plot(phase, model_rs, 'ko', ms=3)
    plt.xlim(-xlim, xlim)
    plt.xlabel(r'$phase$')
    plt.ylabel(r'$residuals$')
    plt.subplots_adjust(hspace=0.0)
    plt.savefig('fit.png', dpi=200)
    os.system("rm model.pickle")
    os.chdir('..')


def transit_ld(data, iterations, burn, binning, ldcoeffs, rprs, rprsvar, p, a, avar, e, i, ivar, w, ww, t0):
    if isinstance(data, basestring):
        datax, datay = np.loadtxt(data, usecols=[0, 1], unpack=True)
    else:
        datax, datay = data
    if binning > 0:
        start = np.mod(len(datax), binning)
        datax = np.mean(np.reshape(datax[start:], (len(datax) / binning, binning)), 1)
        datay = np.mean(np.reshape(datay[start:], (len(datay) / binning, binning)), 1)
    directory = 'mcmc_transit'
    if not os.path.isdir(directory):
        os.mkdir(directory)
        os.chdir(directory)
    else:
        fi = 2
        while os.path.isdir(directory + '_' + str(fi)):
            fi += 1
        directory = directory + '_' + str(fi)
        os.mkdir(directory)
        os.chdir(directory)
    t0 += p * int((datax[-1] - t0) / p)
    fit_rprs = pymc.Uniform('rp', (1 - rprsvar) * rprs, (1 + rprsvar) * rprs, value=rprs)
    fit_a = pymc.Uniform('a', (1 - avar) * a, (1 + avar) * a, value=a)
    fit_i = pymc.Uniform('i', i - ivar, min(90, i + ivar), value=i)
    fit_t0 = pymc.Uniform('t0', -0.05, 0.05, value=0)
    ld1, ld2, ld3, ld4 = ldcoeffs
    fit_ld1 = pymc.Uniform('ld1', min(0.7 * ld1, 1.3 * ld1), max(0.7 * ld1, 1.3 * ld1), value=ld1)
    fit_ld2 = pymc.Uniform('ld2', min(0.7 * ld2, 1.3 * ld2), max(0.7 * ld2, 1.3 * ld2), value=ld2)
    fit_ld3 = pymc.Uniform('ld3', min(0.7 * ld3, 1.3 * ld3), max(0.7 * ld3, 1.3 * ld3), value=ld3)
    fit_ld4 = pymc.Uniform('ld4', min(0.7 * ld4, 1.3 * ld4), max(0.7 * ld4, 1.3 * ld4), value=ld4)
    @pymc.deterministic
    def model_mcmc(model_ld1=fit_ld1, model_ld2=fit_ld2, model_ld3=fit_ld3, model_ld4=fit_ld4,
                   model_rprs=fit_rprs, model_a=fit_a, model_i=fit_i, model_t0=fit_t0):
        model_ldcoeffs = (model_ld1, model_ld2, model_ld3, model_ld4)
        return models.transit(model_ldcoeffs, model_rprs, p, model_a, e, model_i, w, ww, t0 + model_t0, datax)
    y = pymc.Normal('y', mu=model_mcmc, tau=np.std(datay[1:] - datay[:-1]) ** (-2), observed=True, value=datay)
    count0 = time.time()
    mcmc = pymc.MCMC([fit_ld1, fit_ld2, fit_ld3, fit_ld4, fit_rprs, fit_a, fit_i, fit_t0, model_mcmc, y],
                     db='pickle', dbname='model.pickle')
    mcmc.isample(iterations, burn=burn, verbose=1)
    count1 = time.time()
    print 'Sampling finished in {0:.5}mins\n'.format((count1 - count0) / 60)
    pchain = mcmc.trace('rp')[:]
    achain = mcmc.trace('a')[:]
    ichain = mcmc.trace('i')[:]
    t0chain = mcmc.trace('t0')[:]
    ldcoeffs = (np.mean(mcmc.trace('ld1')[:]), np.mean(mcmc.trace('ld2')[:]),
                np.mean(mcmc.trace('ld3')[:]), np.mean(mcmc.trace('ld4')[:]))
    pymc.Matplot.plot(mcmc)
    plt.close('all')
    # save results
    result = [pchain, achain, ichain, t0chain]
    www = open('result.txt', 'w')
    www.write('param \tinitial \tfinal \terror\n')
    www.write('RPORS\t{}\t{}\t{}\n'.format(rprs, np.mean(result[0]), np.std(result[0])))
    www.write('MTP\t{}\t{:.15}\t{}\n'.format(t0, t0 + np.mean(result[3]), np.std(result[3])))
    www.write('INCL\t{}\t{}\t{}\n'.format(i, np.mean(result[2]), np.std(result[2])))
    www.write('SMAORS\t{}\t{}\t{}\n'.format(a, np.mean(result[1]), np.std(result[1])))
    rprs = np.mean(result[0])
    a = np.mean(result[1])
    i = np.mean(result[2])
    t0 = t0 + np.mean(result[3])
    www.write('LDA1\t{}\n'.format(ldcoeffs[0]))
    www.write('LDA2\t{}\n'.format(ldcoeffs[1]))
    www.write('LDA3\t{}\n'.format(ldcoeffs[2]))
    www.write('LDA4\t{}\n'.format(ldcoeffs[3]))
    www.close()
    phase = (datax - t0) / p
    model_lc = models.transit(ldcoeffs, rprs, p, a, e, i, w, ww, t0, datax)
    model_rs = datay - model_lc
    np.savetxt('model_lc.txt', np.swapaxes([datax, phase, datay, model_lc, model_rs], 0, 1))
    np.savetxt('resultlists.txt', np.swapaxes(result, 0, 1))
    # plot
    print 'Plotting model'
    plt.subplot2grid((4, 1), (0, 0), rowspan=3)
    plt.plot(phase, datay, 'bo', mec='b')
    datax2 = np.arange(datax[0], datax[-1], (datax[1] - datax[0]) / 100)
    phase2 = (datax2 - t0) / p
    model_lc2 = models.transit(ldcoeffs, rprs, p, a, e, i, w, ww, t0, datax2)
    plt.plot(phase2, model_lc2, 'r-')
    xlim = np.max(np.abs([plt.xlim()[0], plt.xlim()[1]]))
    plt.xlim(-xlim, xlim)
    plt.text(plt.xlim()[0] + 0.02 * (plt.xlim()[-1] - plt.xlim()[0]),
             plt.ylim()[0] + 0.02 * (plt.ylim()[-1] - plt.ylim()[0]),
             r'$rms_{res} = %.3e$' % np.std(model_rs))
    plt.yticks(plt.yticks()[0][1:])
    plt.ylabel(r'$relative \, flux$')
    plt.tick_params(labelbottom='off')
    plt.subplot(4, 1, 4)
    plt.axhline(0, color='k')
    plt.plot(phase, model_rs, 'ko', ms=3)
    plt.xlim(-xlim, xlim)
    plt.xlabel(r'$phase$')
    plt.ylabel(r'$residuals$')
    plt.subplots_adjust(hspace=0.0)
    plt.savefig('fit.png', dpi=200)
    os.system("rm model.pickle")
    os.chdir('..')
