################# This script plots the spectra for uu, vv, ww, TT and epsilon for kx, ky, kz and E, T, D, vt for kr using the data in x_spec.xxxx, y_spec.xxxx, z_spec.xxxx and r_spec.xxxx files
################# Composer: Dr. Abhiram Doddi
################# Created on: 9th february 2022;
################# Modified: Added radial spectra processing with E, T, D, vt variables
################# NOTE: Needs 'x_spec.xxxx', 'y_spec.xxxx', 'z_spec.xxxx' and 'r_spec.xxxx' files to be present in the script execution directory

################# Input(s):
### tau = Buoyancy period
tau = 1
### savDir = Link to the directory where the output file(s) are saved
savDir = './'  # Save to present working directory
### figW = Figure width
figW = 8
### figSvOf = Figure filename offset for saving
figSvOf = 0
### figH = Figure Height
figH = 6
### figDPI = Figure DPI
figDPI = 300
### figFtSz = Figure text Fontsize
figFtSz = 8
### yLimits = Y-axis limits [min, max] for all plots
yLimits = [10**-12, 10**0]

################# Bright Color Palette Definition
# Define bright colors for different variables
bright_colors = {
    'uu': '#FF1744',   # Bright Red
    'vv': '#2196F3',   # Bright Blue
    'ww': '#4CAF50',   # Bright Green
    'TT': '#FF9800',   # Bright Orange
    'eps': '#E91E63'   # Hot Pink
}

# Alternative bright color options (uncomment to use):
# bright_colors = {
#     'uu': '#E91E63',   # Hot Pink
#     'vv': '#00BCD4',   # Cyan
#     'ww': '#8BC34A',   # Light Green
#     'TT': '#FF5722'   # Deep Orange
# }

# Or use this vibrant palette:
# bright_colors = {
#     'uu': '#FF0000',   # Pure Red
#     'vv': '#0066FF',   # Electric Blue
#     'ww': '#00FF00',   # Pure Green
#     'TT': '#FF6600'   # Bright Orange
# }

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
try:
    import os
    print("os module is loaded...")
except ImportError:
    print("os module not installed...")

################# Create subdirectories for different spectra
# Create subdirectories if they don't exist
x_spectra_dir = os.path.join(savDir, 'x_spectra')
y_spectra_dir = os.path.join(savDir, 'y_spectra')
z_spectra_dir = os.path.join(savDir, 'z_spectra')
r_spectra_dir = os.path.join(savDir, 'r_spectra')

# Create directories
os.makedirs(x_spectra_dir, exist_ok=True)
os.makedirs(y_spectra_dir, exist_ok=True)
os.makedirs(z_spectra_dir, exist_ok=True)
os.makedirs(r_spectra_dir, exist_ok=True)

print(f"Created/verified subdirectories:")
print(f"  X spectra: {x_spectra_dir}")
print(f"  Y spectra: {y_spectra_dir}")
print(f"  Z spectra: {z_spectra_dir}")
print(f"  R spectra: {r_spectra_dir}")

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
    TTX = np.zeros(len(Xspec)-3)
    epsX = np.zeros(len(Xspec)-3) 
    ### dump data from file into data arrays
    for i in range(3,len(Xspec)):
        aX = Xspec[i]
        kx[i-3] = float(aX[1:5])
        uuX[i-3] = float(aX[6:17])
        vvX[i-3] = float(aX[18:29])
        wwX[i-3] = float(aX[30:41])
        TTX[i-3] = float(aX[42:53])
        epsX[i-3] = float(aX[54:65]) 

    ################# Automatically set x limits based on data range
    kx_min = np.min(kx[kx > 0])  # Use minimum positive value for log scale
    kx_max = np.max(kx)
    x_range = kx_max - kx_min
    x_margin = 0.1 * x_range  # 10% margin
    x_limits = [max(kx_min - x_margin, kx_min * 0.5), kx_max + x_margin]
    
    ################# Plot the spectra and save them as image files
    ### for uu
    fig1x = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(kx, uuX, color=bright_colors['uu'], linewidth=2)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(x_limits)
    plt.ylim(yLimits)
    plt.grid(True)
    plt.xlabel(r'$k_{x}$', fontsize=figFtSz)
    plt.ylabel(r'uu', fontsize=figFtSz)
    fig1x.suptitle(r'Plot showing uu spectra in X at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam1x = os.path.join(x_spectra_dir, 'uuX_' + '{:04d}'.format(j+figSvOf) + '.png')
    fig1x.savefig(savNam1x, bbox_inches="tight", format='png')
    ### for vv
    fig2x = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(kx, vvX, color=bright_colors['vv'], linewidth=2)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(x_limits)
    plt.ylim(yLimits)
    plt.grid(True)
    plt.xlabel(r'$k_{x}$', fontsize=figFtSz)
    plt.ylabel(r'vv', fontsize=figFtSz)
    fig2x.suptitle(r'Plot showing vv spectra in X at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam2x = os.path.join(x_spectra_dir, 'vvX_' + '{:04d}'.format(j+figSvOf) + '.png')
    fig2x.savefig(savNam2x, bbox_inches="tight", format='png')
    ### for ww
    fig3x = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(kx, wwX, color=bright_colors['ww'], linewidth=2)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(x_limits)
    plt.ylim(yLimits)
    plt.grid(True)
    plt.xlabel(r'$k_{x}$', fontsize=figFtSz)
    plt.ylabel(r'ww', fontsize=figFtSz)
    fig3x.suptitle(r'Plot showing ww spectra in X at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam3x = os.path.join(x_spectra_dir, 'wwX_' + '{:04d}'.format(j+figSvOf) + '.png')
    fig3x.savefig(savNam3x, bbox_inches="tight", format='png')
    ### for TT
    fig5x = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(kx, TTX, color=bright_colors['TT'], linewidth=2)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(x_limits)
    plt.ylim(yLimits)
    plt.grid(True)
    plt.xlabel(r'$k_{x}$', fontsize=figFtSz)
    plt.ylabel(r'TT', fontsize=figFtSz)
    fig5x.suptitle(r'Plot showing TT spectra in X at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam5x = os.path.join(x_spectra_dir, 'TTX_' + '{:04d}'.format(j+figSvOf) + '.png')
    fig5x.savefig(savNam5x, bbox_inches="tight", format='png')
    ### for eps
    fig4x = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(kx, epsX, color=bright_colors['eps'], linewidth=2)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(x_limits)
    plt.ylim(yLimits)
    plt.grid(True)
    plt.xlabel(r'$k_{x}$', fontsize=figFtSz)
    plt.ylabel(r'$\epsilon$', fontsize=figFtSz)
    fig4x.suptitle(r'Plot showing $\epsilon$ spectra in X at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam4x = os.path.join(x_spectra_dir, 'epsX_' + '{:04d}'.format(j+figSvOf) + '.png')
    fig4x.savefig(savNam4x, bbox_inches="tight", format='png')
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
    TTY = np.zeros(len(Yspec)-3)
    epsY = np.zeros(len(Yspec)-3)

    ### dump data from file into data arrays
    for i in range(3,len(Yspec)):
        aY = Yspec[i]
        ky[i-3] = float(aY[1:5])
        uuY[i-3] = float(aY[6:17])
        vvY[i-3] = float(aY[18:29])
        wwY[i-3] = float(aY[30:41])
        TTY[i-3] = float(aY[42:53])
        epsY[i-3] = float(aY[54:65]) 
    
    ################# Automatically set x limits based on data range
    ky_min = np.min(ky[ky > 0])  # Use minimum positive value for log scale
    ky_max = np.max(ky)
    y_range = ky_max - ky_min
    y_margin = 0.1 * y_range  # 10% margin
    y_limits = [max(ky_min - y_margin, ky_min * 0.5), ky_max + y_margin]
    
    ################# Plot the spectra and save them as image files
    ### for uu
    fig1y = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(ky, uuY, color=bright_colors['uu'], linewidth=2)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(y_limits)
    plt.ylim(yLimits)
    plt.grid(True)
    plt.xlabel(r'$k_{y}$', fontsize=figFtSz)
    plt.ylabel(r'uu', fontsize=figFtSz)
    fig1y.suptitle(r'Plot showing uu spectra in Y at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam1y = os.path.join(y_spectra_dir, 'uuY_' + '{:04d}'.format(j+figSvOf) + '.png')
    fig1y.savefig(savNam1y, bbox_inches="tight", format='png')
    ### for vv
    fig2y = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(ky, vvY, color=bright_colors['vv'], linewidth=2)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(y_limits)
    plt.ylim(yLimits)
    plt.grid(True)
    plt.xlabel(r'$k_{y}$', fontsize=figFtSz)
    plt.ylabel(r'vv', fontsize=figFtSz)
    fig2y.suptitle(r'Plot showing vv spectra in Y at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam2y = os.path.join(y_spectra_dir, 'vvY_' + '{:04d}'.format(j+figSvOf) + '.png')
    fig2y.savefig(savNam2y, bbox_inches="tight", format='png')
    ### for ww
    fig3y = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(ky, wwY, color=bright_colors['ww'], linewidth=2)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(y_limits)
    plt.ylim(yLimits)
    plt.grid(True)
    plt.xlabel(r'$k_{y}$', fontsize=figFtSz)
    plt.ylabel(r'ww', fontsize=figFtSz)
    fig3y.suptitle(r'Plot showing ww spectra in Y at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam3y = os.path.join(y_spectra_dir, 'wwY_' + '{:04d}'.format(j+figSvOf) + '.png')
    fig3y.savefig(savNam3y, bbox_inches="tight", format='png')
    ### for TT
    fig5y = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(ky, TTY, color=bright_colors['TT'], linewidth=2)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(y_limits)
    plt.ylim(yLimits)
    plt.grid(True)
    plt.xlabel(r'$k_{y}$', fontsize=figFtSz)
    plt.ylabel(r'TT', fontsize=figFtSz)
    fig5y.suptitle(r'Plot showing TT spectra in Y at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam5y = os.path.join(y_spectra_dir, 'TTY_' + '{:04d}'.format(j+figSvOf) + '.png')
    fig5y.savefig(savNam5y, bbox_inches="tight", format='png')
    ### for eps
    fig4y = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(ky, epsY, color=bright_colors['eps'], linewidth=2)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(y_limits)
    plt.ylim(yLimits)
    plt.grid(True)
    plt.xlabel(r'$k_{y}$', fontsize=figFtSz)
    plt.ylabel(r'$\epsilon$', fontsize=figFtSz)
    fig4y.suptitle(r'Plot showing $\epsilon$ spectra in Y at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam4y = os.path.join(y_spectra_dir, 'epsY_' + '{:04d}'.format(j+figSvOf) + '.png')
    fig4y.savefig(savNam4y, bbox_inches="tight", format='png')

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
    TTZ = np.zeros(len(Zspec)-3)
    epsZ = np.zeros(len(Zspec)-3)
    ### dump data from file into data arrays
    for i in range(3,len(Zspec)):
        aZ = Zspec[i]
        kz[i-3] = float(aZ[1:5])
        uuZ[i-3] = float(aZ[6:17])
        vvZ[i-3] = float(aZ[18:29])
        wwZ[i-3] = float(aZ[30:41])
        TTZ[i-3] = float(aZ[42:53])
        epsZ[i-3] = float(aZ[54:65])
    
    ################# Automatically set x limits based on data range
    kz_min = np.min(kz[kz > 0])  # Use minimum positive value for log scale
    kz_max = np.max(kz)
    z_range = kz_max - kz_min
    z_margin = 0.1 * z_range  # 10% margin
    z_limits = [max(kz_min - z_margin, kz_min * 0.5), kz_max + z_margin]
    
    ################# Plot the spectra and save them as image files
    ### for uu
    fig1z = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(kz, uuZ, color=bright_colors['uu'], linewidth=2)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(z_limits)
    plt.ylim(yLimits)
    plt.grid(True)
    plt.xlabel(r'$k_{z}$', fontsize=figFtSz)
    plt.ylabel(r'uu', fontsize=figFtSz)
    fig1z.suptitle(r'Plot showing uu spectra in Z at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam1z = os.path.join(z_spectra_dir, 'uuZ_' + '{:04d}'.format(j+figSvOf) + '.png')
    fig1z.savefig(savNam1z, bbox_inches="tight", format='png')
    ### for vv
    fig2z = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(kz, vvZ, color=bright_colors['vv'], linewidth=2)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(z_limits)
    plt.ylim(yLimits)
    plt.grid(True)
    plt.xlabel(r'$k_{z}$', fontsize=figFtSz)
    plt.ylabel(r'vv', fontsize=figFtSz)
    fig2z.suptitle(r'Plot showing vv spectra in Z at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam2z = os.path.join(z_spectra_dir, 'vvZ_' + '{:04d}'.format(j+figSvOf) + '.png')
    fig2z.savefig(savNam2z, bbox_inches="tight", format='png')
    ### for ww
    fig3z = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(kz, wwZ, color=bright_colors['ww'], linewidth=2)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(z_limits)
    plt.ylim(yLimits)
    plt.grid(True)
    plt.xlabel(r'$k_{z}$', fontsize=figFtSz)
    plt.ylabel(r'ww', fontsize=figFtSz)
    fig3z.suptitle(r'Plot showing ww spectra in Z at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam3z = os.path.join(z_spectra_dir, 'wwZ_' + '{:04d}'.format(j+figSvOf) + '.png')
    fig3z.savefig(savNam3z, bbox_inches="tight", format='png')
    ### for TT
    fig5z = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(kz, TTZ, color=bright_colors['TT'], linewidth=2)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(z_limits)
    plt.ylim(yLimits)
    plt.grid(True)
    plt.xlabel(r'$k_{z}$', fontsize=figFtSz)
    plt.ylabel(r'TT', fontsize=figFtSz)
    fig5z.suptitle(r'Plot showing TT spectra in Z at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam5z = os.path.join(z_spectra_dir, 'TTZ_' + '{:04d}'.format(j+figSvOf) + '.png')
    fig5z.savefig(savNam5z, bbox_inches="tight", format='png')
    ### for eps
    fig4z = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(kz, epsZ, color=bright_colors['eps'], linewidth=2)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(z_limits)
    plt.ylim(yLimits)
    plt.grid(True)
    plt.xlabel(r'$k_{z}$', fontsize=figFtSz)
    plt.ylabel(r'$\epsilon$', fontsize=figFtSz)
    fig4z.suptitle(r'Plot showing $\epsilon$ spectra in Z at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam4z = os.path.join(z_spectra_dir, 'epsZ_' + '{:04d}'.format(j+figSvOf) + '.png')
    fig4z.savefig(savNam4z, bbox_inches="tight", format='png')
################################## End of script bit for processing Z spectra

################################## Script bit for processing R spectra 
################# List the files of a particular type - 'r_spec.xxxx'
arrNameR = []
arrNameR = sorted(glob.glob('r_spec.*'))

################# Loop over the number of files
for j in range(0,len(arrNameR)):
    ################# read data from 'z_spec.xxxx' file
    with open(arrNameR[j]) as g:
        Rspec = g.readlines()
    ################# Declare data array sizes
    tmpR = Rspec[0]
    tim = float(tmpR[18:30])/tau
    kr = np.zeros(len(Rspec)-3)
    Er = np.zeros(len(Rspec)-3)
    Tr = np.zeros(len(Rspec)-3)
    Dr = np.zeros(len(Rspec)-3)
    vtr = np.zeros(len(Rspec)-3)
    
    ### dump data from file into data arrays
    for i in range(3,len(Rspec)):
        aR = Rspec[i]
        kr[i-3] = float(aR[1:5])
        Er[i-3] = float(aR[6:17])
        Tr[i-3] = float(aR[18:29])
        Dr[i-3] = float(aR[30:41])
        vtr[i-3] = float(aR[42:53])
    
    ################# Automatically set x limits based on data range
    kr_min = np.min(kr[kr > 0])  # Use minimum positive value for log scale
    kr_max = np.max(kr)
    r_range = kr_max - kr_min
    r_margin = 0.1 * r_range  # 10% margin
    r_limits = [max(kr_min - r_margin, kr_min * 0.5), kr_max + r_margin]
    
    ################# Plot the spectra and save them as image files
    ### for E
    fig1r = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(kr, Er, color=bright_colors['uu'], linewidth=2)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(r_limits)
    plt.ylim(yLimits)
    plt.grid(True)
    plt.xlabel(r'$k_{r}$', fontsize=figFtSz)
    plt.ylabel(r'E', fontsize=figFtSz)
    fig1r.suptitle(r'Plot showing E spectra in R at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam1r = os.path.join(r_spectra_dir, 'Er_' + '{:04d}'.format(j+figSvOf) + '.png')
    fig1r.savefig(savNam1r, bbox_inches="tight", format='png')
    ### for T
    fig2r = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(kr, Tr, color=bright_colors['vv'], linewidth=2)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(r_limits)
    plt.ylim(yLimits)
    plt.grid(True)
    plt.xlabel(r'$k_{r}$', fontsize=figFtSz)
    plt.ylabel(r'T', fontsize=figFtSz)
    fig2r.suptitle(r'Plot showing T spectra in R at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam2r = os.path.join(r_spectra_dir, 'Tr_' + '{:04d}'.format(j+figSvOf) + '.png')
    fig2r.savefig(savNam2r, bbox_inches="tight", format='png')
    ### for D
    fig3r = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(kr, Dr, color=bright_colors['ww'], linewidth=2)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(r_limits)
    plt.ylim(yLimits)
    plt.grid(True)
    plt.xlabel(r'$k_{r}$', fontsize=figFtSz)
    plt.ylabel(r'D', fontsize=figFtSz)
    fig3r.suptitle(r'Plot showing D spectra in R at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam3r = os.path.join(r_spectra_dir, 'Dr_' + '{:04d}'.format(j+figSvOf) + '.png')
    fig3r.savefig(savNam3r, bbox_inches="tight", format='png')
    ### for vt
    fig4r = plt.figure(figsize=(figW, figH), dpi=figDPI)
    plt.plot(kr, vtr, color=bright_colors['TT'], linewidth=2)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(r_limits)
    plt.ylim(yLimits)
    plt.grid(True)
    plt.xlabel(r'$k_{r}$', fontsize=figFtSz)
    plt.ylabel(r'vt', fontsize=figFtSz)
    fig4r.suptitle(r'Plot showing vt spectra in R at t/$\tau_{b}$ = ' + str(format(tim, '.4f')), fontsize=figFtSz, fontweight='bold')
    savNam4r = os.path.join(r_spectra_dir, 'vtr_' + '{:04d}'.format(j+figSvOf) + '.png')
    fig4r.savefig(savNam4r, bbox_inches="tight", format='png')
################################## End of script bit for processing Z spectra

################# Quit out of the program
sys.exit()