#!/usr/bin/env python

import numpy as np
import colorednoise as cn # https://github.com/jleute/colorednoise
import allantools as at # https://github.com/aewallin/allantools  or "pip install allantools"
import matplotlib.pyplot as plt
plt.rc('text', usetex=True) # for latex
plt.rc('font', family='serif')
import math

"""

colorednoise example with phase-PSD, frequency-PSD, and ADEV/MDEV
for different noise types.

Anders Wallin 2016-04-09, update 2018-10-28
"""


def main():

    nr = 2**14 # number of datapoints in time-series
    adev0 = 1.0e-11
    
    #
    # for each noise type generated with cn.noiseGen()
    # compute 
    #   - phase PSD and coefficient g_b
    #   - frequency PSD and coefficient h_a
    #   - ADEV
    #   - MDEV
    #
    tau0 = 0.1 # sample interval for simulated time-series
    sample_rate = 1.0/tau0
    
    qd0 = qd1 = qd2 = qd3 = qd4 = pow(adev0, 2) # discrete variance for noiseGen()
    # This scaling makes all g_i coefficients equal
    # -> PSDs cross at 1 Hz.
    qd1 = qd0*2*np.pi/sample_rate
    qd2 = qd1*2*np.pi/sample_rate
    qd3 = qd2*2*np.pi/sample_rate
    qd4 = qd3*2*np.pi/sample_rate

    
    print("generating timeseries..."),
    x0 = cn.noiseGen(nr, qd0, 0)   # white phase noise (WPM)
    x1 = cn.noiseGen(nr, qd1, -1)  # flicker phase noise (FPM)
    x2 = cn.noiseGen(nr, qd2, -2)  # white frequency noise (WFM)
    x3 = cn.noiseGen(nr, qd3 , -3) # flicker frequency noise (FFM)
    x4 = cn.noiseGen(nr, qd4 , -4) # random walk frequency noise (RWFM)
    
    
    # compute frequency time-series
    y0 = at.phase2frequency(x0, sample_rate)
    y1 = at.phase2frequency(x1, sample_rate)
    y2 = at.phase2frequency(x2, sample_rate)
    y3 = at.phase2frequency(x3, sample_rate)
    y4 = at.phase2frequency(x4, sample_rate)
    print("done.")
    
    print("computing PSDs..."),
    # compute phase PSD
    (f0, psd0) = at.noise.scipy_psd(x0, f_sample=sample_rate, nr_segments=4)
    (f1, psd1) = at.noise.scipy_psd(x1, f_sample=sample_rate, nr_segments=4)
    (f2, psd2) = at.noise.scipy_psd(x2, f_sample=sample_rate, nr_segments=4)
    (f3, psd3) = at.noise.scipy_psd(x3, f_sample=sample_rate, nr_segments=4)
    (f4, psd4) = at.noise.scipy_psd(x4, f_sample=sample_rate, nr_segments=4)
    # compute phase PSD prefactor g_b
    g0 = cn.phase_psd_from_qd( qd0, 0, tau0 )
    g1 = cn.phase_psd_from_qd( qd1, -1, tau0 )
    g2 = cn.phase_psd_from_qd( qd2, -2, tau0 )
    g3 = cn.phase_psd_from_qd( qd3, -3, tau0 )
    g4 = cn.phase_psd_from_qd( qd4, -4, tau0 )
    print g0
    print g1
    print g2
    print g3
    print g4
    

    # compute frequency PSD
    (ff0, fpsd0) = at.noise.scipy_psd(y0, f_sample=sample_rate, nr_segments=4)
    (ff1, fpsd1) = at.noise.scipy_psd(y1, f_sample=sample_rate, nr_segments=4)
    (ff2, fpsd2) = at.noise.scipy_psd(y2, f_sample=sample_rate, nr_segments=4)
    (ff3, fpsd3) = at.noise.scipy_psd(y3, f_sample=sample_rate, nr_segments=4)
    (ff4, fpsd4) = at.noise.scipy_psd(y4, f_sample=sample_rate, nr_segments=4)
    # compute frequency PSD prefactor h_a
    a0 = cn.frequency_psd_from_qd( qd0, 0, tau0 )
    a1 = cn.frequency_psd_from_qd( qd1, -1, tau0 )
    a2 = cn.frequency_psd_from_qd( qd2, -2, tau0 )
    a3 = cn.frequency_psd_from_qd( qd3, -3, tau0 )
    a4 = cn.frequency_psd_from_qd( qd4, -4, tau0 )
    print "a_i == g_i*(4pi^2)"
    print a0, g0*pow(2*np.pi,2)
    print a1, g1*pow(2*np.pi,2)
    print a2, g2*pow(2*np.pi,2)
    print a3, g3*pow(2*np.pi,2)
    print a4, g4*pow(2*np.pi,2)
    
    print("done.")

    print("computing ADEV/MDEV..."),
    # compute ADEV
    (t0,d0,e,n) = at.oadev(x0, rate=sample_rate)
    (t1,d1,e,n) = at.oadev(x1, rate=sample_rate)
    (t2,d2,e,n) = at.oadev(x2, rate=sample_rate)
    (t3,d3,e,n) = at.oadev(x3, rate=sample_rate)
    (t4,d4,e,n) = at.oadev(x4, rate=sample_rate)
    
    # compute MDEV
    (mt0,md0,e,n) = at.mdev(x0, rate=sample_rate)
    (mt1,md1,e,n) = at.mdev(x1, rate=sample_rate)
    (mt2,md2,e,n) = at.mdev(x2, rate=sample_rate)
    (mt3,md3,e,n) = at.mdev(x3, rate=sample_rate)
    (mt4,md4,e,n) = at.mdev(x4, rate=sample_rate)
    print("done.")
    

    plt.figure()
    # Phase PSD figure
    plt.subplot(2,2,1)
    plt.loglog(f0,[ g0*pow(xx, 0.0) for xx in f0],'--',label=r'$g_0f^0 = {h_2\over4\pi^2}f^0$', color='black')
    plt.loglog(f0,psd0,'.', color='black')
    
    plt.loglog(f1[1:],[ g1*pow(xx, -1.0) for xx in f1[1:]],'--',label=r'$g_{-1}f^{-1}  = {h_1\over4\pi^2}f^{-1} $', color='red')
    plt.loglog(f1,psd1,'.' ,color='red')
    
    plt.loglog(f2[1:],[ g2*pow(xx,-2.0) for xx in f2[1:]],'--',label=r'$g_{-2}f^{-2}   = {h_0\over4\pi^2}f^{-2} $', color='green')
    plt.loglog(f2,psd2,'.',  color='green')
    
    plt.loglog(f3[1:], [ g3*pow(xx,-3.0) for xx in f3[1:]],'--',label=r'$g_{-3}f^{-3}   = {h_{-1}\over4\pi^2}f^{-3} $', color='pink')
    plt.loglog(f3, psd3, '.', color='pink')
    
    plt.loglog(f4[1:], [ g4*pow(xx,-4.0) for xx in f4[1:]],'--',label=r'$g_{-4}f^{-4}   = {h_{-2}\over4\pi^2}f^{-4}$', color='blue')
    plt.loglog(f4, psd4, '.', color='blue')
    plt.grid()
    plt.legend(framealpha=0.9, fontsize=22)
    plt.title(r'Phase Power Spectral Density')
    plt.xlabel(r'Frequency (Hz)')
    plt.ylabel(r' $S_x(f)$ $(s^2/Hz)$')
    
    # frequency PSD figure
    plt.subplot(2,2,2)
    plt.loglog(ff0[1:], [ a0*pow(xx,2) for xx in ff0[1:]],'--',label=r'$h_{2}f^{2}$', color='black')
    plt.loglog(ff0,fpsd0,'.',  color='black')
    
    plt.loglog(ff1[1:], [ a1*pow(xx,1) for xx in ff1[1:]],'--',label=r'$h_{1}f^{1}$', color='red')
    plt.loglog(ff1,fpsd1,'.',  color='red')

    plt.loglog(ff2[1:], [ a2*pow(xx,0) for xx in ff2[1:]],'--',label=r'$h_{0}f^{0}$', color='green')
    plt.loglog(ff2,fpsd2,'.', color='green')
    
    plt.loglog(ff3[1:], [ a3*pow(xx,-1) for xx in ff3[1:]],'--',label=r'$h_{-1}f^{-1}$', color='pink')
    plt.loglog(ff3,fpsd3,'.', color='pink')

    plt.loglog(ff4[1:], [ a4*pow(xx,-2) for xx in ff4[1:]],'--',label=r'$h_{-2}f^{-2}$', color='blue')
    plt.loglog(ff4,fpsd4,'.', color='blue')
    plt.grid()
    plt.legend(framealpha=0.9, fontsize=22)
    plt.title(r'Frequency Power Spectral Density')
    plt.xlabel(r'Frequency (Hz)')
    plt.ylabel(r' $S_y(f)$ $(1/Hz)$')
    
    # ADEV figure
    plt.subplot(2,2,3)
    f_h = 0.5/tau0
    K0 = 3.0*f_h / (4.0*pow(np.pi,2))
    plt.loglog(t0, [ np.sqrt(K0*a0)/xx for xx in t0],'--', label=r'$\sqrt{3f_H\over4\pi^2} \sqrt{h_{2}}\tau^{-3/2}$', color='black')
    plt.loglog(t0,d0,'o', color='black')
    
    # IEEE 1139-2008, table B.2
    def K1(f_H, tau):
        gamma = 1.038/2
        return ((3.0/2.0)*np.log(2.0*np.pi*f_h*tau)+gamma)/(2.0*pow(np.pi,2))
    
    plt.loglog(t1, [ np.sqrt( K1(f_h, xx)*a1 )/xx for xx in t1],'--', label=r'$\sqrt{ 3\ln{(2 \pi f_H \tau)} + \gamma \over 4\pi^2 }\sqrt{h_{1}}\tau^{-1}$', color='red')
    plt.loglog(t1,d1,'o', color='red')
    
    K2 = 0.5
    plt.loglog(t2, [ np.sqrt(K2*a2)/math.sqrt(xx) for xx in t2],'--', label=r'$\sqrt{1\over2} \sqrt{h_{0}}\tau^{-1/2}$', color='green')
    plt.loglog(t2,d2,'o', color='green')

    K3 = 2*np.log(2)
    plt.loglog(t3, [ np.sqrt(K3*a3) for xx in t3],'--', label=r'$ \sqrt{2\ln{2}}\sqrt{h_{-1}}\tau^0$', color='pink')
    plt.loglog(t3,d3,'o', color='pink')

    K4 = (2*np.pi**2.0)/3.0
    plt.loglog(t4, [ np.sqrt(K4*a4)*math.sqrt(xx) for xx in t4],'--', label=r'$\sqrt{2\pi^2\over3}\sqrt{h_{-2}}\tau^{+1/2}$', color='blue')
    plt.loglog(t4,d4,'o', color='blue')

    plt.legend(framealpha=0.9, loc='lower left', fontsize=22)
    plt.grid()
    plt.title(r'Allan Deviation')
    plt.xlabel(r'$\tau$ (s)')
    plt.ylabel(r'ADEV')
    
    # MDEV
    plt.subplot(2, 2, 4)
    M0 = 3.0/(8.0*pow(np.pi,2))
    plt.loglog(t0, [ np.sqrt(M0*a0)/pow(xx,3.0/2.0) for xx in t0],'--', label=r'$\sqrt{ 3\over 8\pi^2 } \sqrt{h_{2}}\tau^{-3/2}$', color='black')
    plt.loglog(mt0,md0,'o', color='black')

    M1 = (24.0*np.log(2)-9.0*np.log(3)) / (8.0*pow(np.pi,2)) 
    
    plt.loglog(t1, [ np.sqrt(M1*a1) /xx for xx in t1],'--', label=r'$\sqrt{ 24\ln(2)-9\ln(3) \over 8\pi^2 } \sqrt{h_{1}}\tau^{-1}$', color='red')
    plt.loglog(mt1,md1,'o', color='red')

    M2 = 0.25
    plt.loglog(t2, [ np.sqrt( M2*a2 )/math.sqrt(xx) for xx in t2],'--', label=r'$ \sqrt{1\over4}\sqrt{h_{0}}\tau^{-1/2}$', color='green')
    plt.loglog(mt2,md2,'o', color='green')

    M3 = 27.0/20.0*np.log(2)
    plt.loglog(t3, [ np.sqrt( M3*a3 ) for xx in t3],'--', label=r'$\sqrt{ {27\over20}\ln(2) } \sqrt{h_{-1}}\tau^0$', color='pink')
    plt.loglog(mt3,md3,'o', color='pink')

    M4 = 11.0/20.0*pow(np.pi,2)
    plt.loglog(t4, [  np.sqrt( M4*a4 )*math.sqrt(xx) for xx in t4],'--', label=r'$ \sqrt{ {11\over20}\pi^2  } \sqrt{h_{-2}}\tau^{+1/2}$', color='blue')
    plt.loglog(mt4,md4,'o', color='blue')

    plt.legend(framealpha=0.9, loc='lower left', fontsize=22)
    plt.grid()
    plt.title(r'Modified Allan Deviation')
    plt.xlabel(r'$\tau$ (s)')
    plt.ylabel(r'MDEV')
    plt.show()

if __name__=="__main__":
    main()
