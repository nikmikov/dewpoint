#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.basemap import Basemap
from sys import stdin

import _thread
import json
import traceback
import sys
import queue

"""
Plotting weather simulation data:

Input format (JSON):

HEADER
FRAME
FRAME
....

HEADER:
{
    "title" : TITLE,
    "lats"  : [LATS],
    "lons"  : [LONS]
}

FRAME:

{
    "frame" : "<DATETIME>",
    "temp"  : [t1,..],
    "uwind"  : [u1,..],
    "vwind"  : [v1,..]
}


"""

q = queue.Queue() #blocking queue

def create_basemaps(lats, lons):
    """ Setup global basemaps for eventual plotting """
    print("Creating basemaps for plotting", len(lats), len(lons))

    long, latg = np.meshgrid(lons, lats)
    print(len(long), len(latg))
    # Set up a global map
    bmap_globe = Basemap(projection='cea',llcrnrlat=-70, urcrnrlat=70,\
                         llcrnrlon=0,urcrnrlon=360,lat_ts=20,resolution='c')
    xg,yg = bmap_globe(long,latg)


    return {'global' : bmap_globe,
            'global_x' : xg,
            'global_y' : yg,
            }

def read_header():
    while True:
        try:
            userinput = stdin.readline()
            j = json.loads(userinput)
            if 'lats' in j and 'lons' in j:
                return j
        except Exception:
            traceback.print_exc()


def data_update_thread():
    while True:
        try:
            userinput = stdin.readline()
            j = json.loads(userinput)
            q.put(j)
        except Exception:
            traceback.print_exc()

def main(argv=None):
    print("Plot process started. Waiting for header")
    h = read_header()
    lats, lons = h['lats'], h['lons']
    bmaps = create_basemaps(lats, lons)
    xi = bmaps['global_x']
    yi = bmaps['global_y']
    a_shape = (len(lats),len(lons))

    fig, ax = plt.subplots()
    lev = [x for x in range(-40, 50, 1)]
    def update(data):
        Zc = [ x - 275 for x in data['temp'] ]
        Z = np.reshape(np.asarray(Zc), a_shape)
        U = np.reshape(np.asarray(data['uwind']), a_shape)
        V = np.reshape(np.asarray(data['vwind']), a_shape)
        ax.cla()
        bmaps['global'].drawcoastlines()
        cont = ax.contourf(xi, yi, Z,  extend="both",antialiasing=True, levels = lev)
        ax.hold(True)
        ax.quiver(xi,yi,U, V)
        print("plotted")



    def data_gen():
        while True:
            j = q.get()
            yield j

    ani = animation.FuncAnimation(fig, update, data_gen, interval=100)

    _thread.start_new_thread(data_update_thread, () )
    plt.show()


if __name__ == "__main__":
    sys.exit(main())
