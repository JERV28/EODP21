
# LEVEL-1C
import sys

from l1c.src.initL1c import initL1c
from common.io.writeToa import writeToa, readToa
from common.io.readGeodetic import readGeodetic, getCorners
import mgrs
import numpy as np
from scipy.interpolate import bisplrep, bisplev
import matplotlib.pyplot as plt
from common.io.l1cProduct import writeL1c
from matplotlib import cm

class l1c(initL1c):

    def __init__(self, auxdir, indir, outdir):
        super().__init__(auxdir, indir, outdir)

    def processModule(self):

        self.logger.info("Start of the L1C Processing Module")

        for band in self.globalConfig.bands:

            self.logger.info("Start of BAND " + band)

            # Read TOA - output of the L1B in Radiances
            # -------------------------------------------------------------------------------
            toa = readToa(self.l1bdir, self.globalConfig.l1b_toa + band + '.nc')
            lat,lon = readGeodetic(self.gmdir, self.globalConfig.gm_geoloc)
            self.checkSize(lat,toa)

            # L1C reprojection onto the MGRS grid
            # -------------------------------------------------------------------------------
            lat_l1c, lon_l1c, toa_l1c = self.l1cProjtoa(lat, lon, toa, band)

            # Write output TOA
            # -------------------------------------------------------------------------------
            writeL1c(self.outdir, self.globalConfig.l1c_toa + band, lat_l1c, lon_l1c, toa_l1c)

            self.logger.info("End of BAND " + band)


            lat = getCorners(lat)
            lon = getCorners(lon)
            fig = plt.figure(figsize=(20,10))
            plt.plot(lon, lat, 'k', linewidth=2, label="L1B")
            plt.plot(lon_l1c, lat_l1c, 'r.', markersize=5, label="L1C MGRS")
            plt.title('Projection on ground', fontsize=20)
            plt.xlabel('Longitude [deg]', fontsize=16)
            plt.ylabel('Latitude [deg]', fontsize=16)
            plt.grid()
            plt.axis('equal')
            plt.legend()
            plt.savefig(self.outdir + 'footprint_' + band + '.png')
            plt.close(fig)


        self.logger.info("End of the L1C Module!")


    def l1cProjtoa(self, lat, lon, toa, band):
        '''
        This function reprojects the L1B radiances into the MGRS grid.

        The MGRS reference system
        https://www.bluemarblegeo.com/knowledgebase/calculator-2020/Military_Grid_Reference_System_(MGRS).htm
        MGRS: '31REQ4367374067'
        31 is the UTM zone, R is the UTM latitude band; EQ are the MGRS column and row band letters
        43673 is the MGRS Easting (5 dig); 74067 is the MGRS Northing (5dig)

        Python mgrs library documentation
        https://pypi.org/project/mgrs/

        :param lat: L1B latitudes [deg]
        :param lon: L1B longitudes [deg]
        :param toa: L1B radiances
        :param band: band
        :return: L1C radiances, L1C latitude and longitude in degrees
        '''
        #TODO

        tck = bisplrep(lat,lon,toa)

        Alat=np.array(lat)
        Alon=np.array(lon)

        m = mgrs.MGRS()
        mgrs_tiles = set([])

        for i in range (0,Alat.shape[0]):
            for j in range(0,Alat.shape[1]):
                mgrs_tiles.add(str(m.toMGRS(Alat[i,j],Alon[i,j],inDegrees=True,MGRSPrecision=self.l1cConfig.mgrs_tile_precision)))

        mgrs_tiles=list(mgrs_tiles)

        lat_l1c=[]
        lon_l1c=[]
        toa_l1c=[]

        for k in range (len(mgrs_tiles)):
            Alat, Alon = m.toLatLon(mgrs_tiles[k],inDegrees=True)
            lat_l1c.append(Alat)
            lon_l1c.append(Alon)
            toa_l1c.append(bisplev(Alat, Alon, tck))


        return lat_l1c, lon_l1c, toa_l1c

    def checkSize(self, lat,toa):
        '''
        Check the sizes of the input radiances and geodetic coordinates.
        If they don't match, exit.
        :param lat: Latitude 2D matrix
        :param toa: Radiance 2D matrix
        :return: NA
        '''
        #TODO

        if lat.shape != toa.shape:

            sys.exit()



