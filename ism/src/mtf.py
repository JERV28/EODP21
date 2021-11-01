from math import pi
from config.ismConfig import ismConfig
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.special import j1
from numpy.matlib import repmat
from common.io.readMat import writeMat
from common.plot.plotMat2D import plotMat2D
from scipy.interpolate import interp2d
from numpy.fft import fftshift, ifft2
import os




class mtf:
    """
    Class MTF. Collects the analytical modelling of the different contributions
    for the system MTF
    """
    def __init__(self, logger, outdir):
        self.ismConfig = ismConfig()
        self.logger = logger
        self.outdir = outdir

    def system_mtf(self, nlines, ncolumns, D, lambd, focal, pix_size,
                   kLF, wLF, kHF, wHF, defocus, ksmear, kmotion, directory, band):
        """
        System MTF
        :param nlines: Lines of the TOA
        :param ncolumns: Columns of the TOA
        :param D: Telescope diameter [m]
        :param lambd: central wavelength of the band [m]
        :param focal: focal length [m]
        :param pix_size: pixel size in meters [m]
        :param kLF: Empirical coefficient for the aberrations MTF for low-frequency wavefront errors [-]
        :param wLF: RMS of low-frequency wavefront errors [m]
        :param kHF: Empirical coefficient for the aberrations MTF for high-frequency wavefront errors [-]
        :param wHF: RMS of high-frequency wavefront errors [m]
        :param defocus: Defocus coefficient (defocus/(f/N)). 0-2 low defocusing
        :param ksmear: Amplitude of low-frequency component for the motion smear MTF in ALT [pixels]
        :param kmotion: Amplitude of high-frequency component for the motion smear MTF in ALT and ACT
        :param directory: output directory
        :return: mtf
        """

        self.logger.info("Calculation of the System MTF")

        # Calculate the 2D relative frequencies
        self.logger.debug("Calculation of 2D relative frequencies")
        fn2D, fr2D, fnAct, fnAlt = self.freq2d(nlines, ncolumns, D, lambd, focal, pix_size)

        # Diffraction MTF
        self.logger.debug("Calculation of the diffraction MTF")
        Hdiff = self.mtfDiffract(fr2D)

        # Defocus
        Hdefoc = self.mtfDefocus(fr2D, defocus, focal, D)

        # WFE Aberrations
        Hwfe = self.mtfWfeAberrations(fr2D, lambd, kLF, wLF, kHF, wHF)

        # Detector
        Hdet  = self. mtfDetector(fn2D)

        # Smearing MTF
        Hsmear = self.mtfSmearing(fnAlt, ncolumns, ksmear)

        # Motion blur MTF
        Hmotion = self.mtfMotion(fn2D, kmotion)

        # Calculate the System MTF
        self.logger.debug("Calculation of the Sysmtem MTF by multiplying the different contributors")
        Hsys = Hdiff*Hwfe*Hdefoc*Hdet*Hsmear*Hmotion

        # Plot cuts ACT/ALT of the MTF
        fig=self.plotMtf(Hdiff, Hdefoc, Hwfe, Hdet, Hsmear, Hmotion, Hsys, nlines, ncolumns, fnAct, fnAlt, directory, band)
        fig.savefig("/home/luss/my_shared_folder/test_ISM/MTF_ACT_" + band)

        figtwo=self.plotMtftwo(Hdiff, Hdefoc, Hwfe, Hdet, Hsmear, Hmotion, Hsys, nlines, ncolumns, fnAct, fnAlt, directory, band)
        figtwo.savefig("/home/luss/my_shared_folder/test_ISM/MTF_ALT_" + band)



        return Hsys

    def freq2d(self,nlines, ncolumns, D, lambd, focal, w):
        """
        Calculate the relative frequencies 2D (for the diffraction MTF)
        :param nlines: Lines of the TOA
        :param ncolumns: Columns of the TOA
        :param D: Telescope diameter [m]
        :param lambd: central wavelength of the band [m]
        :param focal: focal length [m]
        :param w: pixel size in meters [m]
        :return fn2D: normalised frequencies 2D (f/(1/w))
        :return fr2D: relative frequencies 2D (f/(1/fc))
        :return fnAct: 1D normalised frequencies 2D ACT (f/(1/w))
        :return fnAlt: 1D normalised frequencies 2D ALT (f/(1/w))
        """
        #TODO
        fstepAlt = 1/nlines/w
        fstepAct = 1/ncolumns/w

        eps=10e-6
        fAlt = np.arange(-1/(2*w),1/(2*w)-eps,fstepAlt)
        fAct = np.arange(-1/(2*w),1/(2*w)-eps,fstepAct)

        fnAlt=fAlt/(1/w)
        fnAct=fAct/(1/w)

        [fnAltxx,fnActxx] = np.meshgrid(fnAlt,fnAct,indexing='ij')
        fn2D=np.sqrt(fnAltxx*fnAltxx + fnActxx*fnActxx)

        ecutoff = D/(lambd*focal)

        fr2D=fn2D*((1/w)/ecutoff)

        return fn2D, fr2D, fnAct, fnAlt

    def mtfDiffract(self,fr2D):
        """
        Optics Diffraction MTF
        :param fr2D: 2D relative frequencies (f/fc), where fc is the optics cut-off frequency
        :return: diffraction MTF
        """
        #TODO
        Hdiff=(2/pi)*(np.arccos(fr2D)-fr2D*(1-fr2D**2)**(1/2))


        return Hdiff


    def mtfDefocus(self, fr2D, defocus, focal, D):
        """
        Defocus MTF
        :param fr2D: 2D relative frequencies (f/fc), where fc is the optics cut-off frequency
        :param defocus: Defocus coefficient (defocus/(f/N)). 0-2 low defocusing
        :param focal: focal length [m]
        :param D: Telescope diameter [m]
        :return: Defocus MTF
        """
        #TODO

        x=pi*defocus*fr2D*(1-fr2D)
        Hdefoc=(2*j1(x))/(x)

        return Hdefoc

    def mtfWfeAberrations(self, fr2D, lambd, kLF, wLF, kHF, wHF):
        """
        Wavefront Error Aberrations MTF
        :param fr2D: 2D relative frequencies (f/fc), where fc is the optics cut-off frequency
        :param lambd: central wavelength of the band [m]
        :param kLF: Empirical coefficient for the aberrations MTF for low-frequency wavefront errors [-]
        :param wLF: RMS of low-frequency wavefront errors [m]
        :param kHF: Empirical coefficient for the aberrations MTF for high-frequency wavefront errors [-]
        :param wHF: RMS of high-frequency wavefront errors [m]
        :return: WFE Aberrations MTF
        """
        #TODO
        Hwfe= np.exp(-fr2D*(1-fr2D)*(kLF*((wLF/lambd)**2)+kHF*((wHF/lambd)**2)))

        return Hwfe

    def mtfDetector(self,fn2D):
        """
        Detector MTF
        :param fnD: 2D normalised frequencies (f/(1/w))), where w is the pixel width
        :return: detector MTF
        """
        #TODO

        Hdet=abs(np.sinc(fn2D))

        return Hdet

    def mtfSmearing(self, fnAlt, ncolumns, ksmear):
        """
        Smearing MTF
        :param ncolumns: Size of the image ACT
        :param fnAlt: 1D normalised frequencies 2D ALT (f/(1/w))
        :param ksmear: Amplitude of low-frequency component for the motion smear MTF in ALT [pixels]
        :return: Smearing MTF
        """
        #TODO

        Hsmear= np.sinc(ksmear*fnAlt)
        Hsmear= np.transpose(repmat(Hsmear, ncolumns, 1))

        return Hsmear

    def mtfMotion(self, fn2D, kmotion):
        """
        Motion blur MTF
        :param fnD: 2D normalised frequencies (f/(1/w))), where w is the pixel width
        :param kmotion: Amplitude of high-frequency component for the motion smear MTF in ALT and ACT
        :return: detector MTF
        """
        #TODO

        Hmotion= np.sinc(kmotion*fn2D)


        return Hmotion

    def plotMtf(self,Hdiff, Hdefoc, Hwfe, Hdet, Hsmear, Hmotion, Hsys, nlines, ncolumns, fnAct, fnAlt, directory, band):
        """
        Plotting the system MTF and all of its contributors
        :param Hdiff: Diffraction MTF
        :param Hdefoc: Defocusing MTF
        :param Hwfe: Wavefront electronics MTF
        :param Hdet: Detector MTF
        :param Hsmear: Smearing MTF
        :param Hmotion: Motion blur MTF
        :param Hsys: System MTF
        :param nlines: Number of lines in the TOA
        :param ncolumns: Number of columns in the TOA
        :param fnAct: normalised frequencies in the ACT direction (f/(1/w))
        :param fnAlt: normalised frequencies in the ALT direction (f/(1/w))
        :param directory: output directory
        :param band: band
        :return: N/A
        """
        #TODO


        Hdiff_A=np.array(Hdiff)
        Hdefoc_A=np.array(Hdefoc)
        Hwfe_A=np.array(Hwfe)
        Hdet_A=np.array(Hdet)
        Hsmear_A=np.array(Hsmear)
        Hmotion_A=np.array(Hmotion)
        Hsys_A=np.array(Hsys)

        H1=Hdiff_A[49]
        H2=Hdefoc_A[49]
        H3=Hwfe_A[49]
        H4=Hdet_A[49]
        H5=Hsmear_A[49]
        H6=Hmotion_A[49]
        H7=Hsys_A[49]


        fig=plt.figure(figsize=(10, 10))
        plt.plot(fnAct, H1, label="Difraction MTF")
        plt.plot(fnAct, H2, label="Defocus MTF")
        plt.plot(fnAct, H3, label="WFE Aberrations MTF")
        plt.plot(fnAct, H4, label="Detector MTF")
        plt.plot(fnAct, H5, label="Smearing MTF")
        plt.plot(fnAct, H6, label="Motion blur MTF")
        plt.plot(fnAct, H7, label="System MTF")

        plt.ylabel("MTF")
        plt.xlabel("Spatial frequencies")
        plt.legend()

        #h=Hsys[:,149]
        #print("MTF=",h)


        return fig


    def plotMtftwo(self,Hdiff, Hdefoc, Hwfe, Hdet, Hsmear, Hmotion, Hsys, nlines, ncolumns, fnAct, fnAlt, directory, band):


        Hdiff_A=np.array(Hdiff)
        Hdefoc_A=np.array(Hdefoc)
        Hwfe_A=np.array(Hwfe)
        Hdet_A=np.array(Hdet)
        Hsmear_A=np.array(Hsmear)
        Hmotion_A=np.array(Hmotion)
        Hsys_A=np.array(Hsys)


        H11=Hdiff_A[:,74]
        H22=Hdefoc_A[:,74]
        H33=Hwfe_A[:,74]
        H44=Hdet_A[:,74]
        H55=Hsmear_A[:,74]
        H66=Hmotion_A[:,74]
        H77=Hsys_A[:,74]


        figtwo=plt.figure(figsize=(10, 10))
        plt.plot(fnAlt, H11, label="Difraction MTF")
        plt.plot(fnAlt, H22, label="Defocus MTF")
        plt.plot(fnAlt, H33, label="WFE Aberrations MTF")
        plt.plot(fnAlt, H44, label="Detector MTF")
        plt.plot(fnAlt, H55, label="Smearing MTF")
        plt.plot(fnAlt, H66, label="Motion blur MTF")
        plt.plot(fnAlt, H77, label="System MTF")

        plt.ylabel("MTF")
        plt.xlabel("Spatial frequencies")
        plt.legend()

        return figtwo
