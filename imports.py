import sys
import site
#sys.path.append(['', '/home/fdavilakurban/miniconda2/envs/orientations/lib/python27.zip', '/home/fdavilakurban/miniconda2/envs/orientations/lib/python2.7', '/home/fdavilakurban/miniconda2/envs/orientations/lib/python2.7/plat-linux2', '/home/fdavilakurban/miniconda2/envs/orientations/lib/python2.7/lib-tk', '/home/fdavilakurban/miniconda2/envs/orientations/lib/python2.7/lib-old', '/home/fdavilakurban/miniconda2/envs/orientations/lib/python2.7/lib-dynload', '/home/fdavilakurban/miniconda2/envs/orientations/lib/python2.7/site-packages','/home/fdavilakurban/illustris/','/home/fdavilakurban/illustris/illustris_python'])
#site.addsitedir(['', '/home/fdavilakurban/miniconda2/envs/orientations/lib/python27.zip', '/home/fdavilakurban/miniconda2/envs/orientations/lib/python2.7', '/home/fdavilakurban/miniconda2/envs/orientations/lib/python2.7/plat-linux2', '/home/fdavilakurban/miniconda2/envs/orientations/lib/python2.7/lib-tk', '/home/fdavilakurban/miniconda2/envs/orientations/lib/python2.7/lib-old', '/home/fdavilakurban/miniconda2/envs/orientations/lib/python2.7/lib-dynload', '/home/fdavilakurban/miniconda2/envs/orientations/lib/python2.7/site-packages','/home/fdavilakurban/illustris/','/home/fdavilakurban/illustris/illustris_python'])
site.addsitedir('/home/fdavilakurban/illustris/') 
site.addsitedir('/home/fdavilakurban/illustris/illustris_python') 
import illustris_python as il
import numpy as np

#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec

import random as ran
import numpy as np

from astropy.io import ascii
from astropy import table
from astropy.table import Table

#from astroML import correlation as cf

from collections import Counter

from scipy import spatial
from scipy.optimize import curve_fit
import scipy.stats

#from pyqt_fit import kde, kde_methods

#from periodic_kdtree import PeriodicCKDTree
#import multiprocessing
#from multiprocessing import Queue
import time 


#import pandas as pd
#import pyfof
import seaborn
seaborn.set_style('ticks')

from statsmodels.distributions.empirical_distribution import ECDF


rdr = ascii.get_reader(Reader=ascii.Basic)
rdr.header.splitter.delimiter = ' '
rdr.data.splitter.delimiter = ' '
rdr.header.start_line = 0
rdr.data.start_line = 1
rdr.data.end_line = None
rdr.header.comment = r'\s*#'
rdr.data.comment = r'\s*#'

plt.close()

sys.path.append('/home/fdavilakurban/illustris')
basePath = '/mnt/mirta3/illustris/Z0'

subgroups = il.groupcat.loadSubhalos(basePath,135)

