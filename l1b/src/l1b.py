
# LEVEL-1B MODULE

from l1b.src.initL1b import initL1b
from common.io.writeToa import writeToa, readToa
from common.src.auxFunc import getIndexBand
from common.io.readFactor import readFactor, EQ_MULT, EQ_ADD, NC_EXT
import numpy as np
import os
import matplotlib.pyplot as plt

class l1b(initL1b):

    def __init__(self, auxdir, indir, outdir):
        super().__init__(auxdir, indir, outdir)

    def processModule(self):

        self.logger.info("Start of the L1B Processing Module")

        for band in self.globalConfig.bands:

            self.logger.info("Start of BAND " + band)

            # Read TOA - output of the ISM in Digital Numbers
            # -------------------------------------------------------------------------------
            toa = readToa(self.indir, self.globalConfig.ism_toa + band + '.nc')

            toa_l1b_ref = readToa("/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-L1B/output/", self.globalConfig.l1b_toa + band + '.nc')



            # Equalization (radiometric correction)
            # -------------------------------------------------------------------------------
            if self.l1bConfig.do_equalization:
                self.logger.info("EODP-ALG-L1B-1010: Radiometric Correction (equalization)")

                # Read the multiplicative and additive factors from auxiliary/equalization/
                eq_mult = readFactor(os.path.join(self.auxdir,self.l1bConfig.eq_mult+band+NC_EXT),EQ_MULT)
                eq_add = readFactor(os.path.join(self.auxdir,self.l1bConfig.eq_add+band+NC_EXT),EQ_ADD)

                # Do the equalization and save to file
                toa = self.equalization(toa, eq_add, eq_mult)
                writeToa(self.outdir, self.globalConfig.l1b_toa_eq + band, toa)

            # Restitution (absolute radiometric gain)
            # -------------------------------------------------------------------------------
            self.logger.info("EODP-ALG-L1B-1020: Absolute radiometric gain application (restoration)")
            toa = self.restoration(toa, self.l1bConfig.gain[getIndexBand(band)])

            # Write output TOA
            # -------------------------------------------------------------------------------
            writeToa(self.outdir, self.globalConfig.l1b_toa + band, toa)
            #self.plotL1bToa(toa, self.outdir, band)

            Diftoa = self.differences(toa, toa_l1b_ref)  # Punto 1

            toa_isrf=readToa("/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-ISM/output/","ism_toa_isrf_" + band + '.nc' )

            figone= self.plottwo(toa, toa_isrf)
            figone.savefig("/home/luss/my_shared_folder/test_L1b/Figure_" + band + 'png')


            self.logger.info("End of BAND " + band)

        self.logger.info("End of the L1B Module!")



    def equalization(self, toa, eq_add, eq_mult):
        """
        Equlization. Apply an offset and a gain.
        :param toa: TOA in DN
        :param eq_add: Offset in DN
        :param eq_mult: Gain factor, adimensional
        :return: TOA in DN, equalized
        """
        #TODO
        toa_out=(toa-eq_add)/(eq_mult)

        return toa_out

    def restoration(self,toa,gain):
        """
        Absolute Radiometric Gain - restore back to radiances
        :param toa: TOA in DN
        :param gain: gain in [rad/DN]
        :return: TOA in radiances [mW/sr/m2]
        """
        #TODO
        toa = toa*gain


        self.logger.debug('Sanity check. TOA in radiances after gain application ' + str(toa[1,-1]) + ' [mW/m2/sr]')

        return toa

    #def plotL1bToa(self, toa_l1b, outputdir, band):
        #TODO


    def differences(self, toa, toa_l1b_ref):
        toaA = np.array(toa)
        toaB = np.array(toa_l1b_ref)
        diffe = toaB - toaA
        Multi = toaB*0.01
        C=0

        for i in range(toaA.shape[0]):

            for j in range(toaA.shape[1]):

                if diffe[i,j] > Multi[i,j]:

                    C=C+1

        if C != 0:
            print("Error")



        return C

    def plottwo(self, toa, toa_isrf):


        toaPA = np.array(toa)
        toaPB = np.array(toa_isrf)

        PA=toaPA[49]
        PB=toaPB[49]
        fig=plt.figure(figsize=(10, 10))
        plt.plot(range(150), PA, label="toa 1")
        plt.plot(range(150), PB, label="toa 2")

        plt.ylabel("TOA")
        plt.xlabel("pixels across track")
        plt.legend()

        return fig






