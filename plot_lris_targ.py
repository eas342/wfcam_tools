from astropy.io import ascii, fits
from matplotlib import pyplot as plt
from astropy.wcs import WCS
import numpy as np
import pdb
from matplotlib.patches import Rectangle

fileName = 'all_phot/images_and_conf/w20161229_01063_sf_st_ngc2506_j_band.fit'
ext = 1

coorFile = '../pan_starrs/pro/output/lris_targsNGC2506.csv'
## Tutorial stuff
# from astropy.utils.data import get_pkg_data_filename
# fileName=get_pkg_data_filename('tutorials/FITS-images/HorseHead.fits')
# ext = 0

#clusterCen = 

def plot_cluster():
    """ Makes a cluster plot with likely targets"""
    
    hdu = fits.open(fileName)[ext]
    wcs = WCS(hdu.header,fix=True)
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection=wcs)
    ax.imshow(hdu.data, origin='lower', cmap=plt.cm.viridis,
              vmin=390,vmax=700)
    ax.set_xlim(790,3350)
    ax.set_ylim(1150,3100)
    #ax.set_xlabel('x')
    #ax.set_ylabel('y')
    
    # Coordinate overlay
    #overlay = ax.get_coords_overlay('fk5')
    #overlay.grid(color='white', ls='dotted')
    #overlay[1].set_axislabel('RA (J2000)')
    #overlay[0].set_axislabel('Dec (J2000)')
    
    #ra = ax.coords[1]
    #dec = ax.coords[0]
    #ra.set_major_formatter('hh:mm:ss')
    #dec.set_major_formatter('dd:mm:ss')
    
    ## Show the targets of interest
    dat = ascii.read(coorFile)
    pts = dat['GROUP'] == 1
    ax.scatter(dat['RA'][pts],dat['DEC'][pts],
               facecolors='none',edgecolors='red',
               transform=ax.get_transform('world'))
    
    # Show the FOVs
    cenX, cenY = 120.004, -10.77
    
    xLRIS, yLRIS = cenX - 2.9/60., cenY - 3.9/60.
    r = Rectangle((xLRIS,yLRIS), 6./60., 7.8/60.,
                  edgecolor='cyan', facecolor='none',
                  transform=ax.get_transform('fk5'))
    ax.text(xLRIS - 0.001, yLRIS,'LRIS',color='cyan',transform=ax.get_transform('fk5'))
    
    widthX, widthY = 2.2/60., 2.2/60.
    
    offsetX = 2.2/60. + 43./3600.
    
    xA, yA = cenX - widthX, cenY - 0.5 * widthY + 0.5/60.
    xB, yB = xA + offsetX, yA
    
    rA = Rectangle((xA,yA), widthX, widthY,
                  edgecolor='orange',facecolor='none',
                  transform=ax.get_transform('fk5'))
    
    rB = Rectangle((xB, yB), widthX, widthY,
                  edgecolor='orange',facecolor='none',
                  transform=ax.get_transform('fk5'))
    ax.text(xA - 0.001, yA,'NIRCam',color='orange',transform=ax.get_transform('fk5'))
    
    ax.add_patch(r)
    ax.add_patch(rA)
    ax.add_patch(rB)
    
    fig.savefig('plots/targs_fovs.pdf',bbox_inches='tight')
    
    
if __name__ == '__main__':
    plot_cluster()
    