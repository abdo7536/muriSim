################# This script plots the spectra for uu, vv, ww, TT and epsilon for kx, ky and kz using the data in x_spec.xxxx, y_spec.xxxx and z_spec.xxxx files
################# Composer: Dr. Abhiram Doddi
################# Created on: 9th february 2022;
################# NOTE: Needs 'x_spec.xxxx', 'y_spec.xxxx' and 'z_spec.xxxx' files to be present in the script execution directory

################# Input(s):
### tau = Buoyancy period
tau = 1
### savDir = Link to the directory where the output file(s) are saved
savDir = '/Volumes/Public/Projects/Documents/Research/MURIDNS_2022_ANALYSIS_ONGOING/SHITDNS/run14/spectra/'
### figW = Figure width
figW = 8
### figSvOf = Figure filename offset for saving
figSvOf = 110
### figH = Figure Height
figH = 6
### figDPI = Figure DPI
figDPI = 300
### figFtSz = Figure text Fontsize
figFtSz = 8

################# All the dependencies to run inbuilt or custom functions that are called within the listed functions below need to be imported beforehand
### This section imports all the required function
try:
    import numpy as np
    print("numpy module is loaded...")
except ImportError:
    print("numpy module not installed...")
try:
    import math
    print("math module is loaded...")
except ImportError:
    print("math module not installed...")
try:
    import matplotlib.pyplot as plt
    print("matplotlib module is loaded...")
except ImportError:
    print("matplotlib module not installed...")
try:
    import time                             
    print("time module is loaded...")
except ImportError:
    print("time module not installed...")
try:
    import glob
    print("glob module is loaded...")
except ImportError:
    print("glob module not installed...")
try:
    import sys
    print("sys module is loaded...")
except ImportError:
    print("sys module not installed...")

################################## Script bit for processing X spectra 
################# List the files of a particular type - 'x_spec.xxxx'
arrNameX = []
arrNameX = sorted(glob.glob('x_spec.*'))

################# Loop over the number of files
for j in range(0,len(arrNameX)):
    ################# read data from 'x_spec.xxxx' file
    with open(arrNameX[j]) as e:
        Xspec = e.readlines()
    ################# Declare data array sizes
    tmpX = Xspec[0]
    tim = float(tmpX[18:30])/tau
    kx = np.zeros(len(Xspec)-3)
    uuX = np.zeros(len(Xspec)-3)
    vvX = np.zeros(len(Xspec)-3)
    wwX = np.zeros(len(Xspec)-3)
    epsX = np.zeros(len(Xspec)-3)
    ### dump data from file into data arrays
    for i in range(3,len(Xspec)):
        aX = Xspec[i]
        kx[i-3] = float(aX[1:5])
        uuX[i-3] = float(aX[6:17])
        vvX[i-3] = float(aX[18:29])
        wwX[i-3] = float(aX[30:41])
        epsX[i-3] = float(aX[42:53])
    ################# Plot the spectra and save them as image files
    ### for uu
    fig1x = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(kx,uuX)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim([min(kx)-1, max(kx)+1])
    plt.ylim([10**-12, 10**0])
    plt.grid(True)
    plt.xlabel(r'$k_{x}$', fontsize=figFtSz)
    plt.ylabel(r'uu', fontsize=figFtSz)
    fig1x.suptitle(r'Plot showing uu spectra in X at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam1x = savDir  + 'uuX_' + '{:04d}'.format(j+figSvOf) + '.jpg'
    fig1x.savefig(savNam1x, bbox_inches="tight", format='jpg')
    ### for vv
    fig2x = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(kx,vvX)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim([min(kx)-1, max(kx)+1])
    plt.ylim([10**-12, 10**0])
    plt.grid(True)
    plt.xlabel(r'$k_{x}$', fontsize=figFtSz)
    plt.ylabel(r'vv', fontsize=figFtSz)
    fig2x.suptitle(r'Plot showing vv spectra in X at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam2x = savDir  + 'vvX_' + '{:04d}'.format(j+figSvOf) + '.jpg'
    fig2x.savefig(savNam2x, bbox_inches="tight", format='jpg')
    ### for ww
    fig3x = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(kx,wwX)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim([min(kx)-1, max(kx)+1])
    plt.ylim([10**-12, 10**0])
    plt.grid(True)
    plt.xlabel(r'$k_{x}$', fontsize=figFtSz)
    plt.ylabel(r'ww', fontsize=figFtSz)
    fig3x.suptitle(r'Plot showing ww spectra in X at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam3x = savDir  + 'wwX_' + '{:04d}'.format(j+figSvOf) + '.jpg'
    fig3x.savefig(savNam3x, bbox_inches="tight", format='jpg')
    ### for eps
    fig5x = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(kx,epsX)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim([min(kx)-1, max(kx)+1])
    plt.ylim([10**-12, 10**0])
    plt.grid(True)
    plt.xlabel(r'$k_{x}$', fontsize=figFtSz)
    plt.ylabel(r'$\epsilon$', fontsize=figFtSz)
    fig5x.suptitle(r'Plot showing $\epsilon$ spectra in X at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam5x = savDir  + 'epsX_' + '{:04d}'.format(j+figSvOf) + '.jpg'
    fig5x.savefig(savNam5x, bbox_inches="tight", format='jpg')
################################## End of script bit for processing X spectra

################################## Script bit for processing Y spectra 
################# List the files of a particular type - 'y_spec.xxxx'
arrNameY = []
arrNameY = sorted(glob.glob('y_spec.*'))

################# Loop over the number of files
for j in range(0,len(arrNameY)):
    ################# read data from 'y_spec.xxxx' file
    with open(arrNameY[j]) as f:
        Yspec = f.readlines()
    ################# Declare data array sizes
    tmpY = Yspec[0]
    tim = float(tmpY[18:30])/tau
    ky = np.zeros(len(Yspec)-3)
    uuY = np.zeros(len(Yspec)-3)
    vvY = np.zeros(len(Yspec)-3)
    wwY = np.zeros(len(Yspec)-3)
    epsY = np.zeros(len(Yspec)-3)
    ### dump data from file into data arrays
    for i in range(3,len(Yspec)):
        aY = Yspec[i]
        ky[i-3] = float(aY[1:5])
        uuY[i-3] = float(aY[6:17])
        vvY[i-3] = float(aY[18:29])
        wwY[i-3] = float(aY[30:41])
        epsY[i-3] = float(aY[42:53])
    ################# Plot the spectra and save them as image files
    ### for uu
    fig1y = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(ky,uuY)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim([min(ky)-1, max(ky)+1])
    plt.ylim([10**-12, 10**0])
    plt.grid(True)
    plt.xlabel(r'$k_{y}$', fontsize=figFtSz)
    plt.ylabel(r'uu', fontsize=figFtSz)
    fig1y.suptitle(r'Plot showing uu spectra in Y at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam1y = savDir  + 'uuY_' + '{:04d}'.format(j+figSvOf) + '.jpg'
    fig1y.savefig(savNam1y, bbox_inches="tight", format='jpg')
    ### for vv
    fig2y = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(ky,vvY)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim([min(ky)-1, max(ky)+1])
    plt.ylim([10**-12, 10**0])
    plt.grid(True)
    plt.xlabel(r'$k_{y}$', fontsize=figFtSz)
    plt.ylabel(r'vv', fontsize=figFtSz)
    fig2y.suptitle(r'Plot showing vv spectra in Y at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam2y = savDir  + 'vvY_' + '{:04d}'.format(j+figSvOf) + '.jpg'
    fig2y.savefig(savNam2y, bbox_inches="tight", format='jpg')
    ### for ww
    fig3y = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(ky,wwY)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim([min(ky)-1, max(ky)+1])
    plt.ylim([10**-12, 10**0])
    plt.grid(True)
    plt.xlabel(r'$k_{y}$', fontsize=figFtSz)
    plt.ylabel(r'ww', fontsize=figFtSz)
    fig3y.suptitle(r'Plot showing ww spectra in Y at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam3y = savDir  + 'wwY_' + '{:04d}'.format(j+figSvOf) + '.jpg'
    fig3y.savefig(savNam3y, bbox_inches="tight", format='jpg')
    ### for eps
    fig5y = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(ky,epsY)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim([min(ky)-1, max(ky)+1])
    plt.ylim([10**-12, 10**0])
    plt.grid(True)
    plt.xlabel(r'$k_{y}$', fontsize=figFtSz)
    plt.ylabel(r'$\epsilon$', fontsize=figFtSz)
    fig5y.suptitle(r'Plot showing $\epsilon$ spectra in Y at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam5y = savDir  + 'epsY_' + '{:04d}'.format(j+figSvOf) + '.jpg'
    fig5y.savefig(savNam5y, bbox_inches="tight", format='jpg')
################################## End of script bit for processing Y spectra

################################## Script bit for processing Z spectra 
################# List the files of a particular type - 'z_spec.xxxx'
arrNameZ = []
arrNameZ = sorted(glob.glob('z_spec.*'))

################# Loop over the number of files
for j in range(0,len(arrNameZ)):
    ################# read data from 'z_spec.xxxx' file
    with open(arrNameZ[j]) as g:
        Zspec = g.readlines()
    ################# Declare data array sizes
    tmpZ = Zspec[0]
    tim = float(tmpZ[18:30])/tau
    kz = np.zeros(len(Zspec)-3)
    uuZ = np.zeros(len(Zspec)-3)
    vvZ = np.zeros(len(Zspec)-3)
    wwZ = np.zeros(len(Zspec)-3)
    epsZ = np.zeros(len(Zspec)-3)
    ### dump data from file into data arrays
    for i in range(3,len(Zspec)):
        aZ = Zspec[i]
        kz[i-3] = float(aZ[1:5])
        uuZ[i-3] = float(aZ[6:17])
        vvZ[i-3] = float(aZ[18:29])
        wwZ[i-3] = float(aZ[30:41])
        epsZ[i-3] = float(aZ[42:53])
    ################# Plot the spectra and save them as image files
    ### for uu
    fig1z = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(kz,uuZ)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim([min(kz)-1, max(kz)+1])
    plt.ylim([10**-12, 10**0])
    plt.grid(True)
    plt.xlabel(r'$k_{z}$', fontsize=figFtSz)
    plt.ylabel(r'uu', fontsize=figFtSz)
    fig1z.suptitle(r'Plot showing uu spectra in Z at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam1z = savDir  + 'uuZ_' + '{:04d}'.format(j+figSvOf) + '.jpg'
    fig1z.savefig(savNam1z, bbox_inches="tight", format='jpg')
    ### for vv
    fig2z = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(kz,vvZ)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim([min(kz)-1, max(kz)+1])
    plt.ylim([10**-12, 10**0])
    plt.grid(True)
    plt.xlabel(r'$k_{z}$', fontsize=figFtSz)
    plt.ylabel(r'vv', fontsize=figFtSz)
    fig2z.suptitle(r'Plot showing vv spectra in Z at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam2z = savDir  + 'vvZ_' + '{:04d}'.format(j+figSvOf) + '.jpg'
    fig2z.savefig(savNam2z, bbox_inches="tight", format='jpg')
    ### for ww
    fig3z = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(kz,wwZ)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim([min(kz)-1, max(kz)+1])
    plt.ylim([10**-12, 10**0])
    plt.grid(True)
    plt.xlabel(r'$k_{z}$', fontsize=figFtSz)
    plt.ylabel(r'ww', fontsize=figFtSz)
    fig3z.suptitle(r'Plot showing ww spectra in Z at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam3z = savDir  + 'wwZ_' + '{:04d}'.format(j+figSvOf) + '.jpg'
    fig3z.savefig(savNam3z, bbox_inches="tight", format='jpg')
    ### for eps
    fig5z = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(kz,epsZ)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim([min(kz)-1, max(kz)+1])
    plt.ylim([10**-12, 10**0])
    plt.grid(True)
    plt.xlabel(r'$k_{z}$', fontsize=figFtSz)
    plt.ylabel(r'$\epsilon$', fontsize=figFtSz)
    fig5z.suptitle(r'Plot showing $\epsilon$ spectra in Z at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam5z = savDir  + 'epsZ_' + '{:04d}'.format(j+figSvOf) + '.jpg'
    fig5z.savefig(savNam5z, bbox_inches="tight", format='jpg')
################################## End of script bit for processing Z spectra

################# Quit out of the program
sys.exit()
