import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from gatspy.periodic import LombScargleFast
import matplotlib
import os
import warnings
warnings.filterwarnings("ignore")


def mazeh_variability_amp(time,flux,period):
    time = time-np.nanmin(time)
    cycle_number = np.floor(time/period)
    amp_list =[]
    #print(np.unique(cycle_number))
    for i in np.unique(cycle_number):
        #print(np.nanstd(flux[cycle_number==i]))
        amp_list.append(np.nanpercentile(flux[cycle_number==i],[95])-np.nanpercentile(flux[cycle_number==i],[5]))
    #print(amp_list)
    return np.nanmedian(np.array(amp_list,'d'))

def bin(time, flux, nbins):
    time_mod = np.ones_like(time)
    t_min = np.min(time)
    length = np.max(time)-np.min(time)
    time_mod[:] = time[:]-t_min
    time_bin = np.zeros(nbins)
    flux_bin = np.zeros(nbins)
    unc_bin = np.zeros(nbins)
    for i in range(nbins):
        inside = np.where(abs(time_mod - length*(i+0.5)/(nbins*1.0)) < length/(nbins*2.0))[0]

        if len(flux[inside])>1:
            time_bin[i] = np.mean(time_mod[inside])#+t_min
            flux_bin[i] = np.mean(flux[inside])
            tmpp= np.sqrt(len(flux[inside])*1.)
            unc_bin[i] = np.std(flux)/tmpp
        elif len(flux[inside]) == 1:
            time_bin[i] = np.mean(time_mod[inside])
            flux_bin[i] = np.mean(flux[inside])
            unc_bin[i] = np.nanstd(flux)
        else:
            time_bin[i] = length*(i+0.5)/(nbins*1.0)
            flux_bin[i] = 0.0
            unc_bin[i] = np.nanstd(flux)
    good =np.where(flux_bin!=0.)[0]
    return time_bin[good]+t_min,flux_bin[good],unc_bin[good]

def plot_fold(time,flux,p_fold):

    fit_i = np.array([1],'i')
    from matplotlib import gridspec
    gs = gridspec.GridSpec(1,5 , width_ratios=[1.,1.,1.,1,1.])

    plt.rcParams.update({'font.size': 15})
    #fig = plt.figure(figsize=(6, 12))
    fig, ax = plt.subplots(1,5, sharey=True,figsize=(16,18))

    plt.subplots_adjust(left=0.1, right=0.9, top=0.85, bottom=0.1)
    fig.subplots_adjust(wspace=0)
    #ax = fig.add_subplot(1,1,1)
    j = 0
    i_tmp = 0
    for i in range(int((np.max(time)-np.min(time))/p_fold)):
        if len(time[(time<(np.min(time)+(i+1)*p_fold))&(time>(np.min(time)+i*p_fold))]) < 10: continue

        ##print(i, i_tmp, j)

        ccolor = 'gray'
        if np.any(fit_i == i): ccolor = 'Blue'
        ax[j].scatter(time[(time<(np.min(time)+(i+1)*p_fold))&(time>(np.min(time)+i*p_fold))] % p_fold, flux[(time<(np.min(time)+(i+1)*p_fold))&(time>(np.min(time)+i*p_fold))]+i_tmp*0.002/4, marker = '.',s = 3, color = ccolor)
        #ax[j].text(3,1+i_tmp*0.002/4,str(i))
        i_tmp += 1
        if i_tmp >= 19:
            j+=1
            i_tmp = 0
    ax[0].set_ylabel('Relative Flux',fontsize = 20)

    for j in range(5):

        ax[j].set_xticks([0,3,6,9,12])
        ax[j].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax[2].set_xlabel('Rotational Phase (days)',fontsize = 20)
    plt.ylim(0.9995,1.010)
    plt.savefig('plots/'+desig+'_fold.pdf')
    #plt.show()


def acf(t,x):
    x = x-np.median(x)

    delta_t = t[1:]-t[0:-1]
    cadence = np.nanmedian(delta_t)

    #print(cadence,delta_t,t[-1],t[0])
    finalnum = int(np.floor((t[-1]-t[0])/cadence)+1)
    ###print 'finalnum',finalnum

    time = []
    flux = []
    t_min = t[0]
    for i in range(finalnum):
        ###print i, t[abs(t-t_min-i*cadence) < 0.25*cadence], len(t[abs(t-t_min-i*cadence) < 0.25*cadence])

        time.append(t_min+i*cadence)
        if len(t[abs(t-t_min-i*cadence) < 0.25*cadence]) >= 1:
            flux.append(x[abs(t-t_min-i*cadence) < 0.25*cadence][0])
        if len(t[abs(t-t_min-i*cadence) < 0.25*cadence]) == 0:
            flux.append(0.)#flux[-1] two choices here either assign 0 or the last flux measurement
                ###print len(t),len(time),len(flux)
    plt.rcParams.update({'font.size': 14})
    plt.close()
    plt.figure(figsize=(8,4))
    plt.subplots_adjust(left=0.15, right=0.9, top=0.85, bottom=0.3)

    #plt.scatter(time,flux, color = 'blue', marker= '.', s =1)
    plt.scatter(t,x, color = 'black', marker= '.', s =1)

    plt.xlabel('BJD')
    plt.ylabel('Flux')
    plt.savefig('plots/'+desig+'_flux.png')


    result = np.correlate(flux, flux, mode='full')
    result = result[int(result.size/2):]
    per_grid = time-t_min


    result[0] = 0.
#let's undersample the acf when for sc data

    if cadence<=30.1/24./60.:
        average_width= 5
    if cadence<=2.1/24./60.:
        average_width= 200

    result_mean = result[0::average_width]
    per_grid_mean = per_grid[0::average_width]
    print('acf_width',average_width,len(result_mean))
    #print(average_width)
    for m in range(len(result_mean)):
        ##print(m, m*average_width,average_width*m+average_width+1,len(result_mean))
        result_mean[m] = np.nanmean(result[m*average_width:np.min([m*average_width+average_width+1,len(result)-1])])
        per_grid_mean[m] = np.nanmean(per_grid[m*average_width:np.min([m*average_width+average_width+1,len(result)-1])])

    result = result_mean
    per_grid = per_grid_mean
    #print( np.nanmin(result))
#plt.scatter(per_grid,result)
#plt.show()


    per_min = per_grid[result == np.nanmin(result)]
    re_max = np.nanmax(result[per_grid>per_min])
    per_max = per_grid[result == re_max]



    p_max_index = np.where(result == re_max)
    ##print('test',p_max_index, p_max_index[0],len(result)-1-p_max_index[0],int(len(result)-1-p_max_index[0]))
    p_higher = per_max
    half = 0.5*re_max
    #print(len(result),p_max_index[0])
    if len(p_max_index[0]) != 1:return per_max,0.,0.
    for k in range(int(len(result)-1-p_max_index[0])):
        dif_tmp = result[int(p_max_index[0]+k+1)]-half

        if dif_tmp >=0: continue
        if dif_tmp <0:
            p_higher = per_grid[int(p_max_index[0]+k+1)]
            break
    ##print('p_higher', p_higher)


    p_lower = per_max
    half = 0.5*re_max
    for k in range(int(-1+p_max_index[0])):
        dif_tmp = result[int(p_max_index[0]-k-1)]-half

        if dif_tmp >=0: continue
        if dif_tmp <0:
            p_lower = per_grid[int(p_max_index[0]-k-1)]
            break
    ##print('p_lower', p_lower)

    #identify local maximum
    from scipy.signal import argrelextrema
    smoothed_result = np.convolve(result,np.ones(10), 'same')
    local_max = argrelextrema(smoothed_result, np.greater)
    #local_max = np.r_[True, result[1:] < result[:-1]] & np.r_[result[:-1] < result[1:], True]#(np.diff(np.sign(np.diff(result))) < 0).nonzero()[0] + 1

    plt.rcParams.update({'font.size': 14})
    plt.close()
    plt.figure(figsize=(8,2))
    plt.subplots_adjust(left=0.15, right=0.9, top=0.85, bottom=0.3)

    plt.plot(per_grid,result)

    plt.xlabel('Period (days)')
    plt.ylabel('ACF Amplitude')

    plt.xlim(0.,np.min([100.,(t[-1]-t[0])*1.1]))
    height = np.max(result)
#plt.ylim([0,1.1*height])

    p_fix = [per_max,p_higher,p_lower]
#p_fix= [12.4,12.4-3.0,12.4+3.0]
#height = 0.000135
#half = 0.00015/2
    plt.plot([p_fix[0],p_fix[0]],[0,1.*height], '-', dashes = [3,1], color = 'red')
    plt.plot([p_fix[1],p_fix[1]],[0,half], '-', dashes = [3,1], color = 'orange')
    plt.plot([p_fix[2],p_fix[2]],[0,half], '-', dashes = [3,1], color = 'orange')

    '''
    newdata = pd.read_csv('local_maximum.txt', header = -1, sep='\s+', names=["order","period","height"])
    order = np.array(newdata['order'],'d')
    p_fix = np.array(newdata['period'],'d')
    power_fix = np.array(newdata['height'],'d')
    for k in range(len(order)):
        plt.plot([p_fix[k],p_fix[k]],[0,power_fix[k]], '-', dashes = [3,1], color = 'red')
    '''
    #plt.scatter(per_grid[local_max],result[local_max], color = 'red')

    plt.locator_params(nbins=4)
#plt.xscale('log')
    #plt.plot([5.1,5.1],[0,0.0004])
    plt.savefig('plots/'+desig+'_acf.pdf')

    plt.close('all')
    fig = plt.figure()
    plt.figure(figsize=(8,4))
    plt.subplots_adjust(left=0.15, right=0.9, top=0.85, bottom=0.3)
    ax = fig.add_subplot()

    t_bin, x_bin, unc_bin = bin(t % per_max,x,100)
    plt.scatter(t % per_max,x,marker = '.', s = 1, alpha = 0.5)
    plt.scatter(t_bin, x_bin,marker = '.', color = 'black')
    plt.ylabel('Relative Flux')
    plt.xlabel('Rotational Phase (days)')
    plt.text(0.9, 0.25, '%.1f' % per_max +'+'+'%.1f' % (p_higher-per_max)+'-'+'%.1f' % (per_max-p_lower)+'days', horizontalalignment='center',verticalalignment='center',transform=ax.transAxes)
    plt.savefig('plots/'+desig+'_fold_acf.png')
###print 'acf result'
#for k in range(2000): ##print per_grid[local_max][k], result[local_max][k]
    return per_max,p_higher,p_lower

def ls_bootstrap(t,x,dx,number):

    power = []
    for i in range(number):
        np.random.shuffle(t)
        ###print 't', t
        model = LombScargleFast().fit(t, x, dx)


        #model.optimizer.period_range=(1., 100.)
        periods = np.exp(np.linspace(np.log(0.1),np.log(100),3000))
        dist = model.periodogram(periods)

        power.append(np.max(dist))
    ###print np.max(dist)
    power = np.array(power)
    #plt.close()
    #plt.hist(power)
    #plt.show()
    return np.max(power)

def lsperi(t,x,dx,title, scramble):
    plt.rcParams.update({'font.size': 14})
    model = LombScargleFast().fit(t, x, dx)
    #plt.close()
    #print('test',np.min(t),np.min(x),np.min(dx))
    #plt.errorbar(t,x,yerr = [dx,dx])
    #plt.show()

    model.optimizer.period_range=(1., np.min([100.,(t[-1]-t[0])*0.99]))
    period = model.best_period
    ###print 'period', period


    #periods, power = model.periodogram_auto(nyquist_factor=20)
    periods = np.exp(np.linspace(np.log(1.),np.log(np.min([100.,(t[-1]-t[0])*1.1])),500))
    power = model.periodogram(periods)


    #bootstrap to find the FAP

    ###print 't_orig', t
    t_tmp = np.zeros_like(t)
    t_tmp[:] = t[:]
    if scramble==True: fap_level = ls_bootstrap(t_tmp,x,dx,100)


    ###print 't_final', t
    #for k in range(len(periods)): print(periods[k],power[k])
    #plt.close()
    #plt.plot(periods,power)
    #plt.show()

    power_dif = power[1:]-power[0:-1]

    min_dev = np.min(abs(np.max(power)-power))
    p_max = periods[abs(power-np.max(power)) == min_dev]
    p_max_index = np.where(abs(power-np.max(power)) == min_dev)[0]

    #if len(power[p_max_index][0]) == 0:
    half = power[p_max_index][0]*0.5
    #print('ok',half)
    '''
        half = np.max(power[power_dif > 0.])*0.5
        min_dev = np.min(abs(half-power))
        p_higher = periods[abs(power-half) == min_dev]
        ##print 'p_higher', p_higher
        '''
    p_higher = p_max
    ###print 'test', (len(power)-1-p_max_index[0])
    for k in range(int(len(power)-1-p_max_index[0])):
        dif_tmp = power[int(p_max_index[0]+k+1)]-half
        #print(k, power[int(p_max_index[0]+k+1)], half)
        #dif_tmp = power[p_max_index[0]+k+1]-power[p_max_index[0]+k]
        if dif_tmp >=0: continue
        if dif_tmp <0:
            p_higher = periods[int(p_max_index[0]+k+1)]
            break
    ##print 'p_higher', p_higher

    p_lower = p_max

    for k in range(int(p_max_index[0]-1)):
        dif_tmp = power[int(p_max_index[0]-k-1)]-half
        ###print k, power[p_max_index[0]-k-1], half
        #dif_tmp = power[p_max_index[0]+k+1]-power[p_max_index[0]+k]
        if dif_tmp >=0: continue
        if dif_tmp <0:
            p_lower = periods[int(p_max_index[0]-k-1)]
            break
    ##print 'p_lower', p_lower

    plt.close()
    plt.figure(figsize=(8,2))
    plt.subplots_adjust(left=0.15, right=0.9, top=0.85, bottom=0.3)

    plt.plot(periods,power, linewidth = 1)
    plt.plot([period,period],[0,np.max(power)])
    plt.plot([p_higher,p_higher],[0,0.5*np.max(power)])
    #plt.xlim([0,100])






    plt.xlabel('Period (days)')
    plt.ylabel('Power')







    plt.xlim(0.,np.min([100.,(t[-1]-t[0])*1.1]))




    height = np.max(power)

    #--------------

    if scramble==True:
        plt.plot([.1,100],[fap_level,fap_level],'-', dashes = [3,1],color = 'gray',alpha = 0.5)
        height = np.max([fap_level, np.max(power)])
        #plt.text(0.5,0.85*fap_level,'FAP = 0.001')
        ##print 'fap_level', fap_level, height,power
    #--------------

    #plt.ylim([-0.5*height,1.1*height])
    #plt.text(10,0.8*height,title)
    '''
    p_fix = [period,p_higher,p_lower]
    plt.plot([p_fix[0],p_fix[0]],[0,1.*height], '-', dashes = [3,1], color = 'red', zorder = 0, linewidth = 0.5)
    plt.plot([p_fix[1],p_fix[1]],[0,half], '-', dashes = [3,1], color = 'orange', zorder = 0, linewidth = 0.5)
    plt.plot([p_fix[2],p_fix[2]],[0,half], '-', dashes = [3,1], color = 'orange', zorder = 0, linewidth = 0.5)
    '''
#plt.plot([p_fix[2],p_fix[2]],[0,1.2*height], '-', dashes = [3,1])
#plt.plot([p_fix[3],p_fix[3]],[0,1.2*height], '-', dashes = [3,1])
#plt.plot([p_fix[4],p_fix[4]],[0,1.2*height])
    plt.title(title)
    plt.locator_params(nbins=4)
#plt.xscale('log')
    #plt.plot([5.1*1.6,5.1*1.6],[0,0.04])
    plt.savefig('plots/'+desig+title+'_periodogram.pdf')
#plt.show()

    plt.close('all')
    fig = plt.figure()
    plt.figure(figsize=(8,4))
    plt.subplots_adjust(left=0.15, right=0.9, top=0.85, bottom=0.3)
    ax = fig.add_subplot()

    t_bin, x_bin, unc_bin = bin(t % period,x,100)

    plt.scatter(t % period,x,marker = '.', s = 1, alpha = 0.5)
    plt.scatter(t_bin, x_bin,marker = '.', color = 'black')
    plt.ylabel('Relative Flux')
    plt.xlabel('Rotational Phase (days)')
    plt.text(0.9, 0.25, '%.1f' % period +'+'+'%.1f' % (p_higher-period)+'-'+'%.1f' % (period-p_lower)+'days', horizontalalignment='center',verticalalignment='center',transform=ax.transAxes)
    plt.savefig('plots/'+desig+'_fold_ls.png')


    return period,p_higher, p_lower


newdata = pd.read_csv('list.txt', header = None, sep='\s+', names=["epic"])#,"cam","can"
epic = np.array(newdata['epic'],'str')
#cam = np.array(newdata['cam'],'str')
#can_test = np.array(newdata['can'],'i')
#can_test = np.ones(len(epic),'i')

# Set can_test false for now to just run
can_test = np.zeros(len(epic),'i')

###print epic

txt_file = open('result.txt','w')
txt_file.write('TIC P_acf P_upper_acf P_lower_acf P_ls  P_upper_ls P_lower_ls Amplitude\n')
for i in range(0,len(epic)):#len(epic)
    print('Working on', i, epic[i])
    file = './lc/raw/'+epic[i]+'.txt'
    #file = '../../usp_rotation/data/'+epic[i]+'_raw.txt'
    print(file)
    if os.path.exists(file) == False:
        print('missing ', file)
        continue
    newdata = pd.read_csv('./lc/raw/'+epic[i]+'.txt', sep=',', header = 0, names=["Time", "Flux"])
    #newdata = pd.read_csv('../lc/raw/'+epic[i]+'.txt', header = None, names=["index","Time", "Flux"])
    #newdata = pd.read_table('data/hlsp_k2sff_k2_lightcurve_'+epic[i]+'-'+str(cam[i])+'_kepler_v1_llc-default-aper.txt', sep=',', header = 0, names=["Time", "Flux", "tmp"])
    ###print newdata
    time = np.array(newdata['Time'],'d')
    flux = np.array(newdata['Flux'],'d')
    #unc = np.array(newdata['unc'],'d')
    unc = np.ones_like(flux)*np.std(flux[0:20])
    #plt.scatter(time,flux)
    #plt.show()
    if can_test[i] ==1:
        width = 1.5
        tran = pd.read_csv('../input/'+epic[i]+"_input.txt", header = 0, sep='\s+', names = ['Rp/Rs','b','u1','u2','Rs/a','P','tc','duration','e','omega'])
        p = np.array(tran['P'], 'd')
        tc = np.array(tran['tc'], 'd')
        dura = np.array(tran['duration'], 'd')/24.0
        n_p = len(p)
        depth = np.array(tran['Rp/Rs'],'d')**2

        time_no_tran = time[:]
        flux_no_tran = flux[:]
        unc_no_tran = unc[:]
        for i_p in range(n_p):
            non_tran = np.where((((time_no_tran - tc[i_p]) % p[i_p]) >= width * dura[i_p]) & (((time_no_tran - tc[i_p]) % p[i_p]) <= (p[i_p]-width * dura[i_p])))
            flux_no_tran = flux_no_tran[non_tran]
            time_no_tran = time_no_tran[non_tran]
            unc_no_tran = unc_no_tran[non_tran]
        time = time_no_tran
        flux = flux_no_tran
        unc = unc_no_tran
    #plt.scatter(time,flux,color = 'green')
    #plt.show()
#good = [time<1000]
#time = time[good]
#flux = flux[good]
#unc = unc[good]

    global desig
    desig = epic[i]


    #remove outliers


    std = np.nanstd(flux)
    print('test',std)
    good = np.where(abs(flux-np.nanmedian(flux))< 5.*std)
    time = time[good]
    flux = flux[good]
    unc = unc[good]

    #plot_fold(time,flux,12.9553886426)#12.41245
    '''
    plot_fold(time[(time>000) & (time<200)],flux[(time>000) & (time<200)],12.41245)
    plot_fold(time[(time>200) & (time<400)],flux[(time>200) & (time<400)],12.41245)
    plot_fold(time[(time>400) & (time<600)],flux[(time>400) & (time<600)],12.41245)
    plot_fold(time[(time>600) & (time<800)],flux[(time>600) & (time<800)],12.41245)
    plot_fold(time[(time>800) & (time<1000)],flux[(time>800) & (time<1000)],12.41245)
    plot_fold(time[(time>1000) & (time<1200)],flux[(time>1000) & (time<1200)],12.41245)
    '''


    #fit away any linear trend
    n  = 2
    a = np.polyfit(time,flux,n)
    flux_poly = np.zeros_like(flux)#*a[n]
    for i in range(n+1):
        flux_poly +=a[i]*time**(n-i)
    flux /= flux_poly
    #plt.plot(time,flux_poly)
    #plt.scatter(time,flux)
    #print(len(time),len(flux))
    #plt.show()
    #print((np.nanmax(time)-np.nanmin(time))/0.2)

    nbins = int((np.nanmax(time)-np.nanmin(time))/0.05)
    time_bin, flux_bin, unc_bin = bin(time, flux, nbins)
    #for k in range(len(flux_bin)): print(flux_bin[k])
    #print(np.mean(flux_bin))
    #print('ok',np.min(time_bin),len(time_bin),len(flux_bin),len(unc_bin))
    #plt.close()
    #plt.scatter(time_bin, flux_bin)
    #plt.show()

    #print('no of data points', len(time),len(time_bin))
    print('amp',mazeh_variability_amp(time,flux,5.1),mazeh_variability_amp(time,flux,10.2),np.nanstd(flux))
    period_acf,p_higher_acf, p_lower_acf = acf(time,flux)
    period,p_higher, p_lower = lsperi(time_bin,flux_bin,unc_bin,'', False)


    txt_file.write(desig +'  ')
    txt_file.write(str(period_acf[0]) +'  ')
    txt_file.write(str(p_higher_acf-period_acf[0]) +'  ')
    txt_file.write(str(period_acf[0]-p_lower_acf) +'  ')
    txt_file.write(str(period) +'  ')
    txt_file.write(str(p_higher-period) +'  ')
    txt_file.write(str(period-p_lower) +'  ')
    txt_file.write(str(np.percentile(flux,90)-np.percentile(flux,10)) +'\n')

txt_file.close()
#steps 1) remove 5 sigma outlier; 2) remove a linear trend; 3) ls; 4) acf
