#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

__author__  = "Martin De Kauwe"
__version__ = "1.0 (07.03.2018)"
__email__   = "mdekauwe@gmail.com"

data = np.genfromtxt("x")

x = np.arange(len(data))
y = data[:,0]
yh = data[:,0] + data[:,1]
yl = data[:,0] - data[:,1]
plt.plot(x, y, "k-")
plt.fill_between(x, yl, yh)
plt.show()
