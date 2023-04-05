# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 12:29:59 2023

@author: Ignatio Denton
"""
import numpy as np
import matplotlib.pyplot as plt
from PWEM_model import Conv_Mat

class device():
  
  # initialize class objects
  def __init__(self):
    self.ER  = []
    self.UR  = []
    self.ERC = []
    self.URC = []
  
  # Plot the unit cell
  def Plot_Device(self,params):
    plt.figure()
    plt.imshow(np.real(self.ER), extent = [-params.Lx/2, params.Lx/2, -params.Ly/2, params.Ly/2], 
               aspect = 'auto',  cmap = 'jet')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.title('Unit cell Re ${\epsilon(x,y)}$')
    plt.colorbar()
    plt.savefig("ER.pdf", format="pdf", bbox_inches="tight")

  # Plot convmat
  def Plot_ERC(self, params):
    plt.figure()
    Re = np.real(self.ERC); 
    Im = np.imag(self.ERC);
    plt.subplot(1, 2, 1)
    plt.imshow(Re, aspect = 'auto', cmap = 'cividis')
    plt.colorbar()
    plt.subplot(1, 2, 2)
    plt.imshow(Im, aspect = 'auto', cmap = 'cividis')
    plt.colorbar()
    plt.savefig("ERC.pdf", format="pdf", bbox_inches="tight")

  # Define Ellipse 
  def Ellipse(self, params):
    nx, ny   = params.dim
    ER       = np.zeros((nx,ny))
    X, Y     = np.meshgrid(params.x, params.y)                      # Meshgrid for R
    ER[:,:]  = (X/params.ax)**2 + (Y/params.ay)**2 >= params.r**2;  # Define the circle of ones
    self.ER  = ER*(params.er2 - params.er1) + params.er1
    self.ERC = Conv_Mat(self.ER, params)
    return self
