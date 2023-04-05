# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 12:26:49 2023

@author: Ignatio Denton
"""
#######################################################################

import matplotlib.pyplot as plt
from Device import device
from Params import params, pwem_params
from PWEM_model import PWEM_2D
from Plot import Plot_Results

#######################################################################

plt.rcParams["figure.figsize"] = (6.4, 4.8)
plt.rcParams.update({'font.size': 18,
                       "text.usetex": True});

#######################################################################

params = params();
device = device();  
device.Ellipse(params);
device.Plot_Device(params)

#######################################################################

pwem_params = pwem_params(params)
pwem_params.Symmetry();

#######################################################################

WN   = PWEM_2D.Run(params,device, pwem_params)

if pwem_params.BC == 0:
    Data = Plot_Results(WN, 1 , 10).Plot_Contours3D(params, pwem_params)
else:
    Data = Plot_Results(WN).Plot_Bands(params, pwem_params, 10, 1)

#######################################################################