import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

def plot_map(service='World_Physical_Map', epsg=4269, xpixels=1000):
    # note, you need change the epsg for different region,
    # US is 4269, and you can google the region you want
    plt.figure(figsize=(8, 8))
    m = Basemap(projection='mill',llcrnrlon=4,llcrnrlat=46,urcrnrlon=16,urcrnrlat=56, resolution='l', epsg=epsg)

    # xpixels controls the pixels in x direction, and if you leave ypixels
    # None, it will choose ypixels based on the aspect ratio
    m.arcgisimage(service=service, xpixels=xpixels, verbose=False)

    plt.show()

plot_map(service = 'ESRI_Imagery_World_2D',epsg=31467)