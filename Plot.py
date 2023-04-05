# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 18:27:46 2023

@author: Ignatio Denton
"""

import matplotlib.pyplot as plt
import numpy as np

class Plot_Results():
    
    """
    
        W     -> the eigenvalue matrix
        sheet -> the eigen-surface (iso-freq surf)
        bx,by -> wavevector components
        m     -> Bloch mode number
        lines -> contour lines
    """
    
######################################################################################################            
    
    def __init__(self, W, m = 1, lines = 10):
        
        self.W     = W
        self.sheet = []
        self.bx    = []
        self.by    = []
        self.m     = m
        self.lines = lines
    
    def Plot_Bands(self, params, pwem_params,k0_num, ylim):
        
######################################################################################################            
            
        if pwem_params.BC == 1:
            
            beta  = ["$\Gamma$", "$X$", "$M$", "$\Gamma$"]
            
######################################################################################################
        
        elif pwem_params.BC == 2:
            
            beta = ["$\Gamma$", "$X$", "$M$", "$\Delta$", "$\Gamma$"]
            
######################################################################################################

        plt.figure()
        for i in range(0,k0_num):
            plt.plot(self.W[i,:], 'blue')
            
        for i in range(1,len(beta) - 1):
            plt.axvline(x =  pwem_params.N_Points*i,   color = 'gray',ls='--')
            
######################################################################################################
        
        Ticks = [i * (pwem_params.N_Points) for i in range(0,len(beta))]
        plt.xticks(Ticks,beta)
        
######################################################################################################

        plt.xlim(0, pwem_params.N_Points*(len(beta)-1) )
        plt.ylim(0, ylim)

        plt.xlabel(r'$\vec{\beta}$')
        plt.ylabel(r'$a/\lambda_0$')
        
        plt.title(r'${}$ $mode$'.format(params.Mode))
        plt.savefig("Bands{}.pdf".format(params.Mode), format="pdf", bbox_inches="tight")
    
    def Plot_Contours3D(self,params, pwem_params):
        
######################################################################################################
        
        Harmonic = int(params.Harmonics[0] * params.Harmonics[1])
        sheet    = np.zeros((pwem_params.N_Points, pwem_params.N_Points, Harmonic))
        for i in range(0, Harmonic):
            sheet[:,:,i] = np.reshape(self.W[i,:], (pwem_params.N_Points, pwem_params.N_Points))
            
######################################################################################################
        
        x    = np.linspace(-np.pi,  np.pi, pwem_params.N_Points)
        y    = np.linspace( np.pi, -np.pi, pwem_params.N_Points)
        x, y = np.meshgrid(x,y)
        
######################################################################################################
        
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(10, 10))

        m, lines = self.m, self.lines

        ax.view_init(20, -30)
        
######################################################################################################
        
        surf = ax.plot_surface(x, y, sheet[:,:,m], cmap="jet", alpha  = 0.5)
        cp   = ax.contour(x, y, sheet[:,:,m], lines, cmap="jet", offset = 0)
        
######################################################################################################
        
        ax.clabel(cp, inline=True, 
                  fontsize=18, offset = 0)
        ax.set_zlim(0, np.max(sheet[:,:,m]))
        ax.set_xlim(-np.pi, np.pi)
        ax.set_ylim(-np.pi, np.pi)
        
######################################################################################################
        
        ax.set_xlabel(r'$\beta_x$')
        ax.set_ylabel(r'$\beta_y$')
        ax.set_zlabel(r'$a/\lambda_0$')
        
######################################################################################################
        
        ax.set_aspect("auto")
        plt.colorbar(surf, pad = 0.05, shrink=0.5, aspect=20, orientation = 'horizontal')
        
######################################################################################################
        
        if   m == 1:
            ax.set_title('${}$nd Bloch mode, ${}$ $mode$'.format(m+1, params.Mode))
        elif m == 0:
            ax.set_title('${}$st Bloch mode, ${}$ $mode$'.format(m+1, params.Mode))
        elif m == 2:
            ax.set_title('${}$rd Bloch mode, ${}$ $mode$'.format(m+1, params.Mode))
        else:
            ax.set_title('${}$th Bloch mode, ${}$ $mode$'.format(m+1, params.Mode))

        plt.savefig("Mode pol {}, Bloch {}.pdf".format(params.Mode, m+1), format="pdf")
        
######################################################################################################
        
        del ax, fig
        
        plt.figure()
        contours = plt.contour(x, y, sheet[:,:,m], lines, colors='black')
        plt.clabel(contours, inline=True, fontsize=8)
        
        plt.imshow(sheet[:,:,m], extent=[np.min(x), np.max(x), np.min(y), np.max(y)], origin='lower',
                   cmap='jet', alpha=0.3) # cmap RdGy or jet
        cbar = plt.colorbar();
        
        cbar.ax.set_ylabel(r'$a/\lambda_0$')

        plt.xlabel(r'$\beta_x$')
        plt.ylabel(r'$\beta_y$')
        
        if   m == 1:
            plt.title('${}$nd Bloch mode, ${}$ $mode$'.format(m+1, params.Mode))
        elif m == 0:
            plt.title('${}$st Bloch mode, ${}$ $mode$'.format(m+1, params.Mode))
        elif m == 2:
            plt.title('${}$rd Bloch mode, ${}$ $mode$'.format(m+1, params.Mode))
        else:
            plt.title('${}$th Bloch mode, ${}$ $mode$'.format(m+1, params.Mode))
        
        plt.savefig("Mode pol {}, Bloch {} 2D Contour.pdf".format(params.Mode, m+1), format="pdf")
        
######################################################################################################
        
        self.sheet = sheet
        self.bx    = x
        self.by    = y
        
######################################################################################################
        
        return self
        