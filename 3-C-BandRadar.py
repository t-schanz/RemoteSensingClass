import matplotlib.pyplot as plt
import struct
from mpl_toolkits.basemap import Basemap
import numpy as np
import os

class RadarData:
    from datetime import datetime as dt
    from datetime import timedelta
    import numpy as np
    import sys
    import urllib.request

    def __init__(self,strftime="latest",forecast=0):
        """
        :param strftime: Datestring of the format YYYYMMDDhhmm. Example: "2017082150"
        :param forecast: minutes to forecast. 0=now. Only steps of 5 are allowed. Max=120
        """
        self.data_raw = []
        self.data = self.np.zeros([900,900]).astype(float)
        self.forecast_time = forecast # time to forecast (max 120 min). Just multiples of 5 allowed
        if strftime != "latest":
            year_str = strftime[0:4]
            month_str = strftime[4:6]
            day_str = strftime[6:8]
            hour_str = strftime[8:10]
            minute_str = strftime[10:12]

            try: # check for correct input
                self.time = self.dt(int(year_str), int(month_str), int(day_str), int(hour_str), int(minute_str))
                self.time_parameter = "datetime"

                if self.np.mod(self.time.minute,5) != 0: # if input minute is not devidable by 5 go to last 5 min. step
                    self.time = self.time - self.timedelta(minutes=self.np.mod(self.time.minute,5))

            except:
                print("Not a valid datestring")
                self.sys.exit(0)

        else:
            self.time = self.dt.now()
            self.time_parameter = "now"

        self.__downloadRadarData()

    def __downloadRadarData(self):
        PAGE = "http://opendata.dwd.de/weather/radar/composit/fx/"

        if self.time_parameter == "now":
            FILE = "FXLATEST_%03d_MF002"%self.forecast_time
        else:
            __time_str = self.time.strftime("%y%m%d%H%M")
            FILE = "FX" + __time_str + "_%03d_MF002"

        __radarFile = self.urllib.request.urlopen(PAGE + FILE)
        with open(FILE, 'wb') as output:
            output.write(__radarFile.read())

        with open(FILE,"rb") as f:
            text = f.read()
            self.header= text[:157]
            data = text[158:]
            counter = 0
            for x in range(0, len(data), 2):
                sample = data[x:x + 2]
                value = struct.unpack('<h', sample)[0]
                if value == 0:
                    dbz = self.np.nan
                elif value == 10692:
                    dbz = -999
                else:
                    dbz = (value / 10) / 2 - 32.5


                i = self.np.mod(counter,900)
                j = int(counter/900)
                self.data[i,j] = dbz

                counter += 1

        os.remove(FILE)

if __name__ == "__main__":
    Radar = RadarData()

    hamburg_coords = (53.551086, 9.993682)
    fig, ax1 = plt.subplots(nrows=1,ncols=1,figsize=(9,9))
    m = Basemap(llcrnrlon=4,llcrnrlat=46,urcrnrlon=16,urcrnrlat=56, epsg=31467)
    m.arcgisimage(service='ESRI_Imagery_World_2D', xpixels=1000, verbose=True)
    # m.drawcoastlines(linewidth=0.25)
    # m.drawcounties(linewidth=0.25)
    m.drawmeridians(np.arange(-180, 180, 2))
    m.drawparallels(np.arange(-90, 90, 2))
    lats,lons = m.makegrid(900,900)
    X,Y = m(lats,lons)
    img = m.contourf(x=X,
                     y=Y,
                     data=Radar.data.transpose(),cmap="jet")

    hamburg_x,hamburg_y = m(hamburg_coords[0],hamburg_coords[1])
    m.scatter(x=hamburg_x,y=hamburg_y, marker="o", markersize=2000, markerfacecolor="black",zorder=100)
    # plt.colorbar(img,cax=m)
    plt.savefig("Lates.png")