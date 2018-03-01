import numpy as np
import matplotlib.pyplot as plt
import glob
from astropy.table import Table
from astropy.io import fits, ascii
from astropy.coordinates import SkyCoord
from astropy import units as u
import os

## My goal here is to merge catalogs for each cluster

## I want a table like this, with averaged photometry
## Object   RA  Dec J   H   K

fileTable = ascii.read('ukirt_images_summary.csv')
## The NGC 6811 files are sometimes called "Base"
baseCols = (fileTable['Source'] == 'Base')
fileTable['Source'][baseCols] = 'NGC 6811'

## Columns of interest
## Mike Irwin suggests 1arsec radius aperture phot
## (ie Aper_flux_3)

useColumns = ['Aper']

class ukirt_catalog(object):
    """ Class to store information on a UKIRT catalog table
    """
    def __init__(self,catFile,chip=1):
        """ Initiate with a file name
        Parameters
        -------------
        catFile: str
            Catalog file name
        chip: int
            Chip number
        """
        HDUList = fits.open(catFile)
        self.catFile = catFile
        
        self.head0 = HDUList[0].header
        self.head1 = HDUList[chip].header
        self.dat1 = Table(HDUList[chip].data)
        self.coors = SkyCoord(ra=self.dat1['RA'] * u.radian,
                              dec=self.dat1['DEC'] * u.radian)
        self.filter = self.head1['FILTER']
        self.object = self.head1['OBJECT']
        
        basename = os.path.splitext(os.path.basename(catFile))[0]
        self.descrip = "{}_{}_{}".format(self.object,self.filter,basename)
        
        self.FOVcenter = SkyCoord(ra=self.head1['TELRA'] * u.hourangle,
                                  dec=self.head1['TELDEC'] * u.degree)
        
        HDUList.close()
    
    def distortion_cor(self,saveplot=False):
        """ Distortion correction function
        See http://adsabs.harvard.edu/abs/2009MNRAS.394..675H
        Parameters
        --------------
        r: Astropy SyCoord object
            Distance from the 
            Number of radians away from 
        """
        sep = self.FOVcenter.separation(self.coors)
        r = sep.radian
        k3 = self.head1['PV2_3']
        fcor = 1./((1. + 3. * k3 * r**2) * (1. + k3 * r**2))
        mi = -2.5 * np.log10(fcor)
        if saveplot == True:
            plt.plot(sep.arcmin,mi * (-1.),'o')
            plt.xlabel('Distance (arcmin)')
            plt.ylabel('$\Delta$ mag')
            plt.savefig('plots/'+self.descrip+'_distortion_cor.png')
            plt.close()
        
        return mi
    
    def get_mags(self,aper='3'):
        """ Return the magnitude for a given aperture
        Source: see http://wsa.roe.ac.uk//flatFiles.html for formula
        or http://adsabs.harvard.edu/abs/2009MNRAS.394..675H
        """
        fl = self.dat1['Aper_flux_{}'.format(aper)]
        flerr = self.dat1['Aper_flux_{}_err'.format(aper)]
        extCorr = (-1. * self.head1['EXTINCT'] * 
                  ((self.head1['AMSTART']+self.head1['AMEND'])/2.-1.))
        texpCorr = 2.5 * np.log10(self.head1['EXP_TIME'])
        ApCor = self.head1['APCOR{}'.format(aper)]
        MagZpt = self.head1['MAGZPT']
        magDistortion = self.distortion_cor()
        m = -2.5 * np.log10(fl) + MagZpt + extCorr + texpCorr - ApCor + magDistortion
        merr = 2.5 * np.log10(np.e) * flerr/fl
        return m, merr
    
    def make_simpleTable(self):
        """ Makes a simple table of results """
        t = Table()
        t['Ra (deg)'] = self.coors.ra.deg
        t['Dec (deg)'] = self.coors.dec.deg
        m, merr = self.get_mags()
        t["{} mag".format(self.filter)] = m
        t["{} mag err".format(self.filter)] = merr
        t["Classification"] = self.dat1['Classification']
        return t
    

def merge_files(cluster='NGC 2420'):
    """ Merges all files for a cluster
    Parameters
    -------------------
    cluster: str
        Cluster - either 'NGC 6811', 'NGC 2420' or 'NGC 2506'
    """
    
    print('Merging {} ...'.format(cluster))
    
    clusterInd = (fileTable['Source'] == cluster)
    cols = np.unique(fileTable['Filter'][clusterInd])
    firstTimeThrough = True
    for oneCol in cols:
        exposureInd = (fileTable['Filter'] == oneCol) & clusterInd
        exposures = fileTable['File'][exposureInd]
        for oneExp in exposures:

            subData = da
            if firstTimeThrough:
                t = Table()
                
                t['']