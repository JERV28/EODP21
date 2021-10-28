import os

from common.io.writeToa import writeToa, readToa
from ism.src.initIsm import initIsm
from math import pi
from ism.src.mtf import mtf
from numpy.fft import fftshift, ifft2, fft2
import numpy as np
from common.io.writeToa import writeToa
from common.io.readIsrf import readIsrf
from scipy.interpolate import interp1d, interp2d
from common.plot.plotMat2D import plotMat2D
from common.plot.plotF import plotF
from scipy.signal import convolve2d
from common.src.auxFunc import getIndexBand
import matplotlib.pyplot as plt

class opticalPhase(initIsm):

    def __init__(self, auxdir, indir, outdir):
        super().__init__(auxdir, indir, outdir)

    def compute(self, sgm_toa, sgm_wv, band):
        """
        The optical phase is in charge of simulating the radiance
        to irradiance conversion, the spatial filter (PSF)
        and the spectral filter (ISRF).
        :return: TOA image in irradiances [mW/m2/nm],
                    with spatial and spectral filter
        """
        self.logger.info("EODP-ALG-ISM-1000: Optical stage")

        # Calculation and application of the ISRF
        # -------------------------------------------------------------------------------
        self.logger.info("EODP-ALG-ISM-1010: Spectral modelling. ISRF")
        toa = self.spectralIntegration(sgm_toa, sgm_wv, band)
        toa_isrf=toa
        self.logger.debug("TOA [0,0] " +str(toa[0,0]) + " [e-]")

        if self.ismConfig.save_after_isrf:
            saveas_str = self.globalConfig.ism_toa_isrf + band
            writeToa(self.outdir, saveas_str, toa)

        # Radiance to Irradiance conversion
        # -------------------------------------------------------------------------------
        self.logger.info("EODP-ALG-ISM-1020: Radiances to Irradiances")
        toa = self.rad2Irrad(toa,
                             self.ismConfig.D,
                             self.ismConfig.f,
                             self.ismConfig.Tr)

        self.logger.debug("TOA [0,0] " +str(toa[0,0]) + " [e-]")



        # Spatial filter
        # -------------------------------------------------------------------------------
        # Calculation and application of the system MTF
        self.logger.info("EODP-ALG-ISM-1030: Spatial modelling. PSF/MTF")
        myMtf = mtf(self.logger, self.outdir)
        Hsys = myMtf.system_mtf(toa.shape[0], toa.shape[1],
                                self.ismConfig.D, self.ismConfig.wv[getIndexBand(band)], self.ismConfig.f, self.ismConfig.pix_size,
                                self.ismConfig.kLF, self.ismConfig.wLF, self.ismConfig.kHF, self.ismConfig.wHF,
                                self.ismConfig.defocus, self.ismConfig.ksmear, self.ismConfig.kmotion,
                                self.outdir, band)

        # Apply system MTF
        toa = self.applySysMtf(toa, Hsys) # always calculated
        self.logger.debug("TOA [0,0] " +str(toa[0,0]) + " [e-]")



        # Write output TOA & plots
        # -------------------------------------------------------------------------------
        if self.ismConfig.save_optical_stage:
            saveas_str = self.globalConfig.ism_toa_optical + band

            writeToa(self.outdir, saveas_str, toa)

            title_str = 'TOA after the optical phase [mW/sr/m2]'
            xlabel_str='ACT'
            ylabel_str='ALT'
            plotMat2D(toa, title_str, xlabel_str, ylabel_str, self.outdir, saveas_str)

            idalt = int(toa.shape[0]/2)
            saveas_str = saveas_str + '_alt' + str(idalt)
            plotF([], toa[idalt,:], title_str, xlabel_str, ylabel_str, self.outdir, saveas_str)



            #----------------------------------------------------------------------

            writeToa(self.outdir, self.globalConfig.l1b_toa + band, toa)
            #self.plotL1bToa(toa, self.outdir, band)




            toa_isrf_L=readToa("/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-ISM/output/", self.globalConfig.ism_toa_isrf + band + '.nc' )
            toa_optical=readToa("/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-ISM/output/","ism_toa_optical_" + band + '.nc' )

            Diftoa = self.differences(toa, toa_optical)  #
            DiftoaI = self.differencesTwo(toa, toa_isrf_L)  #

            figone= self.plottwo(toa_isrf, toa_isrf_L)
            figone.savefig("/home/luss/my_shared_folder/test_ISM/Figure_" + band + 'png')

            figtwo= self.plottwo2(toa, toa_optical)
            figtwo.savefig("/home/luss/my_shared_folder/test_ISM/Figure2_" + band + 'png')

            Tr = self.ismConfig.Tr
            print(Tr)






        return toa

    def rad2Irrad(self, toa, D, f, Tr):
        """
        Radiance to Irradiance conversion
        :param toa: Input TOA image in radiances [mW/sr/m2]
        :param D: Pupil diameter [m]
        :param f: Focal length [m]
        :param Tr: Optical transmittance [-]
        :return: TOA image in irradiances [mW/m2]
        """
        toa = Tr*toa*((pi/4)*(D/f)**2)


        return toa


    def applySysMtf(self, toa, Hsys):
        """
        Application of the system MTF to the TOA
        :param toa: Input TOA image in irradiances [mW/m2]
        :param Hsys: System MTF
        :return: TOA image in irradiances [mW/m2]
        """
        # TODO
        GE= fft2(toa)
        Hsys2= fftshift(Hsys)
        toa_ft=np.real(ifft2(GE*Hsys2))


        return toa_ft

    def spectralIntegration(self, sgm_toa, sgm_wv, band):
        """
        Integration with the ISRF to retrieve one band
        :param sgm_toa: Spectrally oversampled TOA cube 3D in irradiances [mW/m2]
        :param sgm_wv: wavelengths of the input TOA cube
        :param band: band
        :return: TOA image 2D in radiances [mW/m2]
        """
        # TODO

        isrf, wv_isrf = readIsrf(os.path.join(self.auxdir,self.ismConfig.isrffile),band)
        wv_isrf=wv_isrf*1000
        isrf_n=isrf/np.sum(isrf)

        sgm_toa=np.array(sgm_toa)
        fsgm=np.zeros([sgm_toa.shape[0],sgm_toa.shape[1]])
        for i in range(sgm_toa.shape[0]):
            for j in range(sgm_toa.shape[1]):
                cs=interp1d(sgm_wv,sgm_toa[i,j,:],fill_value=(0,0),bounds_error=False)
                toa_i=cs(wv_isrf)
                fsgm[i,j]=np.sum(isrf_n*toa_i)

        return fsgm




##-----------------------


    def differences(self, toa, toa_optical):
        toaA = np.array(toa)
        toaB = np.array(toa_optical)
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



    def differencesTwo(self, toa, toa_isrf_L):
        toaA = np.array(toa)
        toaB = np.array(toa_isrf_L)
        diffe = toaB - toaA
        Multi = toaB*0.01
        D=0

        for i in range(toaA.shape[0]):

            for j in range(toaA.shape[1]):

                if diffe[i,j] > Multi[i,j]:

                    D=D+1

        if D != 0:
            print("Error")



        return D













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


    def plottwo2(self, toa, toa_optical):


        toaPA = np.array(toa)
        toaPB = np.array(toa_optical)

        PA=toaPA[49]
        PB=toaPB[49]
        fig=plt.figure(figsize=(10, 10))
        plt.plot(range(150), PA, label="toa 1")
        plt.plot(range(150), PB, label="toa 2")

        plt.ylabel("TOA")
        plt.xlabel("pixels across track")
        plt.legend()

        return fig
