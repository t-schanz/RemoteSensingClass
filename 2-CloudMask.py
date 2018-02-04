
import matplotlib

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from datetime import datetime as dt
from calendar import timegm
import argparse
import sys

from scipy.ndimage.morphology import binary_closing,binary_opening

from scipy import ndimage

# =========================================
# Radar Class
# =========================================
class RadarClass:
    """
    Class for working with Radar data. \n
    For initialisation it needs the path of a netCDF-file with Radar-data. \n
    functions: \n
        Internal:\n
        - __init__:          initializes the class. :param nc_FilePath: Path of a netCDF-file containing Radar-data\n
        - __time2obj:        converts epoche-time to datetime-objects\n
        - __getCloudMask:    calculates the cloud mask\n
        - __getRainMask:     calculates the rain mask\n
        - below:             just for testing and debugging\n
        - __dt_obj:          converts seconds since epoche to a datetime-object\n
        - __ut_obj:          converts a datetime-obejct to seconds since epoche\n
        External:\n
        - print_nc_infos:    prints information of the netCDF-file \n
        - time:              returns a 1D array containing the time as datetime-objects \n
        - range:             returns a 1D array containing the range-gates in m \n
        - Zf:                returns a 1D array containing the reflectivity in dBz \n
        - MeltHeight:        returns a 1D array containing the height of the melting layer in m \n
        - VEL:               returns a 1D array containing the doppler vertical velocity in m/s \n
        - cloudMask():       returns a 2D array containing a cloud-mask (see __getClodMask for more info) \n
        - rainMask():        returns a 2D array containing a rain-mask (see __getRainMask for more info) \n
        - rainRate():        retunrs a 2D array containing the rain intensity in mm/h \n
    """
    def __init__(self, nc_FilePath:str):
        """
        :param nc_FilePath: Path of a netCDF-file containing Radar-data
        """
        nc = Dataset(nc_FilePath)
        self._time = nc.variables['time'][:].copy()  # time in sec since 1970
        self._BCOtime = nc.variables["bco_day"][:].copy()
        self._range = nc.variables['range'][:].copy()  # range in m
        self._Zf = nc.variables['Zf'][:].copy()  # filtered reflectivity
        self._Zf[np.where(np.less(self._Zf,-50))] = np.nan # remove all Values < -50 dbz
        self._MeltHeight = nc.variables['MeltHei'][:].copy()  # Height of the meltinglayer in m
        self._VEL = nc.variables['VEL'][:].copy()  # vertical velocity of all hydrometeors
        self._VELg = nc.variables['VELg'][:].copy() # Doppler velocity of all targets
        self._LDR = nc.variables['LDR'][:].copy() # Linear depolarization rate of all Hydrometeors
        self._LDRg = nc.variables['LDRg'][:].copy() # Linear depolarization rate of all targets



        self.print_nc_infos(nc)
        nc.close()

        self._cloudMask = np.asarray(self._Zf).copy()
        self._rainMask = np.asarray(self._Zf).copy()
        self._rainRate = np.asarray(self._Zf).copy()
        self._notSphericMask = np.asarray(self._VELg).copy()

        self._obj_time = []
        self._cloudFraction = self.__getCloudFraction()

        # functions to be executed on initiation:
        self.__time2obj()  # creates datetime-objects from "time"

        self.__filter()

        self.__getCloudMask()  # creates a cloudmask
        self.__getRainMask()  # creates a rainmask
        self.__getNotSpheric() #gets a mask for all not spheric particles


    def __filter(self):
        self._Zf[self._Zf < -62] = np.nan

    def __time2obj(self):
        for element in self._time:
            self._obj_time.append(dt.replace(self.__dt_obj(element),tzinfo=None))
        self._obj_time = np.asarray(self._obj_time)

    def __getNotSpheric(self):
        self._notSphericMask[self._notSphericMask > -99999] = False
        # self._notSphericMask[]

    def __getCloudMask(self):
        """
        Creates a cloud mask.
        Values stand for:
        0 = cloud-beard
        1 = rain
        2 = cloud above melting layer
        3 = cloud below melting layer
        :return:
        """
        # TODO: Unterscheidung zwischen Cirrus und Cumulus via LDR and Melting Layer hight

        self._cloudMask[self._cloudMask > -99999] = np.nan # set complete array to nan
        # self._cloudMask[np.logical_and(self.VEL() < -1, self._Zf <= 0)] = 0# downdraft in cloud
        self._cloudMask[np.logical_and((self.Zf() > -50), (self.VEL() > -1))] = 30  # cloud
        self._cloudMask[self.Zf() <= -50] = 0  # cloud-beards

        self._cloudMask[self.VEL() <= -2] = 1  # rain

        __Below0CloudMask = np.asarray(self._cloudMask).copy()
        __Below0CloudMask = [__Below0CloudMask == -99999][0]
        __Below0CloudBeard = np.asarray(self._cloudMask).copy()
        __Below0CloudBeard = [__Below0CloudBeard == -99999][0]


        for i in range(len(self._time)):
            __Below0CloudMask[i][np.where(self._range < self._MeltHeight[i])] = True
            __Below0CloudBeard[i][np.where(self._range < self._MeltHeight[i])] = True

        self._cloudMask[np.logical_and(~__Below0CloudMask,self._cloudMask==30)] = 2 # everything above Melting layer height is cirrus
        self._cloudMask[np.logical_and(__Below0CloudMask,self._cloudMask==30)] = 3  # everything below Melting layer height is cumulus
        self._cloudMask[np.logical_and(~__Below0CloudBeard,self._cloudMask==0)] = np.nan  # cloud-beards just occur below melting layer height
        del __Below0CloudMask, __Below0CloudBeard


    def __getCloudFraction(self):
        sums = np.nansum(self.Zf(),axis=1)
        cloudy = 0

        for i in sums:
            if i!=0:
                cloudy += 1


        return np.divide(cloudy,len(sums))


    def __getRainMask(self):
        self._rainMask[self._rainMask > -99999] = np.nan # set complete array to nan
        self._rainRate[self._rainRate > -99999] = np.nan # set complete array to nan

        # rain where VEL < 1 and Zf > 0:
        self._rainMask[self.VEL() <= -2] = self._Zf[self.VEL() <= -2]  # rain

        # Precipitation only below melting-layer:
        __Below0C = np.asarray(self._rainMask)
        __Below0C = [__Below0C == -99999][0]
        for i in range(len(self._time)):
            __Below0C[i][np.where(self._range < self._MeltHeight[i])] = True

        self._rainMask[~__Below0C] = np.nan
        del __Below0C

        # Marshall-Palmer z-R relationship:
        # RR = 0.036 * 10^(0.0625 * dBZ)
        self._rainRate = np.multiply(0.036, np.power(10,np.multiply(0.0625 ,self._rainMask)))

        # throw away to small values
        self._rainRate[self._rainRate < 0.1] = np.nan

    def print_nc_infos(self, nc):
        print(nc)
        print('---------------------------------------')
        for key in nc.variables.keys():
            print(key)

    def __dt_obj(self, u):
        return dt.utcfromtimestamp(u)

    def __ut_obj(self, d):
        return timegm(d.timetuple())

    def time(self, shape:str="dt"):
        """
        :param shape: what shape the returned values should have. Allowed: \n
                "dt" = datetime-object \n
                "ut" = seconds since 1970
        :return: datetime-obj or seconds since 1970, depending on parameter-settings
        """

        if shape == "dt":
            return self._obj_time

        elif shape == "ut":
            return self._time

    def range(self):
        return self._range

    def Zf(self):
        return self._Zf

    def MeltHeight(self):
        return self._MeltHeight

    def VEL(self):
        return self._VEL

    def cloudMask(self):
        return self._cloudMask

    def rainMask(self):
        return self._rainMask

    def rainRate(self):
        return self._rainRate

    def LDR(self):
        return self._LDR

    def LDRg(self):
        return self._LDRg

    def BCOtime(self):
        return self._BCOtime



# ----------------------------------------------------------------------------------------------------------------------




def smooth(Radar):
    """
    This function uses binary opening and closing to get rid of weak signals which could be false interpreted as clouds. (binary opening)
    Furthermore it brings together parts of a cloud into one patch, so that it will later be counted as one cloud. (binary closing)
    :param Radar: Object of the Radar-class
    :return: smoothed cloudmask. Do only use this for counting clouds.
    """
    CumulusCloudMask = Radar.cloudMask().copy()
    CumulusCloudMask[CumulusCloudMask != 3] = np.nan
    CumulusCloudMask[np.isnan(CumulusCloudMask)] = 0
    CumulusCloudMask =  binary_opening(CumulusCloudMask, iterations=2).astype(int)
    CumulusCloudMask = binary_closing(CumulusCloudMask, iterations=10).astype(float)
    CumulusCloudMask[CumulusCloudMask == 0] = np.nan

    CirrusCloudMask = Radar.cloudMask().copy()
    CirrusCloudMask[CirrusCloudMask != 2] = np.nan
    CirrusCloudMask[np.isnan(CirrusCloudMask)] = 0
    CirrusCloudMask = binary_closing(CirrusCloudMask, iterations=30).astype(int)
    CirrusCloudMask =  binary_opening(CirrusCloudMask, iterations=5).astype(float)

    CirrusCloudMask[CirrusCloudMask == 0] = np.nan

    CM_smooth = CirrusCloudMask.copy()
    CM_smooth[np.logical_or((CirrusCloudMask == 1),(CumulusCloudMask==1))] = 1
    return CM_smooth

def countClouds(Radar):
    smoothed = smooth(Radar)
    smoothed[np.isnan(smoothed)] = 0
    mask = smoothed > np.mean(smoothed)
    structure_array = np.ones([3,3])
    label_im, nb_labels = ndimage.label(mask, structure=structure_array)
    print("Clouds in picture: %i" %nb_labels)
    label_im = label_im.astype(float)
    label_im[label_im == 0] = np.nan
    # plt.contourf(label_im.transpose())
    return label_im, nb_labels


def plotResults(Radar):
    fig, (ax1,ax2,ax3,ax4) = plt.subplots(nrows=4,sharex=True,figsize=(16,9))


    # a):
    ref = ax1.contourf(Radar.time(),Radar.range(),Radar.Zf().transpose(),cmap="jet")
    ax1.set_ylim(0,15000)
    ax1.set_ylabel("Height [m]")
    cb1 = fig.colorbar(ref,ax=ax1,label="Reflectivity [dBz]")

    # b):


    # c):
    easyRain = Radar.Zf()
    easyRain[~np.isnan(easyRain)] = 1
    easyRain[Radar.VEL() > -2] = np.nan
    ax3.contourf(Radar.time(),Radar.range(),easyRain.transpose(),cmap="Accent")
    ax3.set_ylabel("Height [m]")
    rain_patch = mpatches.Patch(facecolor="darkblue", label='Precipitation (where velocity < -2 m/s)')
    ax3.legend(handles=[rain_patch])

    # d)
    mp = ax4.contourf(Radar.time(),Radar.range(),Radar.rainRate().transpose(),cmap="jet")
    cb4 = fig.colorbar(mp, ax=ax4,label="Rain Rate [mm/h]")
    ax4.set_ylabel("Height [m]")





    plt.savefig("NR2.png")



if __name__ == "__main__":

    # Set nc-file depending on parsed argument:
    NC_FILE = "C:/Users/darkl/Downloads/Task2_CloudMaskAndCloudFraction/Task2_CloudMaskAndCloudFraction/MMCR__MBR__Spectral_Moments__10s__155m-25km__170923.nc"


    SAVE_PATH = ""
    SAVE_FILE = NC_FILE[:-3] + "_CloudMask.nc"

    # Initiate class:
    Radar = RadarClass(NC_FILE)


    # contourf_plot(Radar, Radar.rainRate
    # mask,count = countClouds(Radar)
    print(Radar._cloudFraction)
    # plotCloudmask(Radar, Radar.Zf(),count)
    plotResults(Radar)
