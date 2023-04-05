# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 12:28:40 2023

@author: Ignatio Denton
"""
import numpy as np
from scipy.linalg import eigh

def Conv_Mat(A, Params):
    
    # Extract spatial harmonics (P, Q, R) of a general 3D unit cell
    
    P, Q    = Params.Harmonics;
    
    # Extract shape of an array
    
    Nx, Ny = np.shape(A);
    
    # Spatial harmonic indices
    
    NH = P * Q;                                 # Total num of spatial harmonics
    
    p  = np.array(np.arange(-np.floor(P/2), np.floor(P/2) + 1))  # Idx in x dir
    q  = np.array(np.arange(-np.floor(Q/2), np.floor(Q/2) + 1))  # Idx in y dir
    
    # Array indices for the zeroth harmonic
    
    p_0 = int(np.floor(Nx/2)); # add +1 in matlab
    
    q_0 = int(np.floor(Ny/2));
    
    # Fourier coefficients of A
    
    A  = np.fft.fftshift(np.fft.fftn(A) / (Nx * Ny)); # Ordered Fourier coeffs
    
    # Init Convolution matrix;
    
    C  = np.zeros((NH, NH), dtype = 'complex');
    
    # Looping
    for q_row in range(1,Q + 1):
      for p_row in range(1,P + 1):
        row = (q_row - 1) * P + p_row
        for q_col in range(1,Q + 1):
          for p_col in range(1,P + 1):
            col   = (q_col - 1) * P + p_col
            p_fft = int(p[p_row - 1] - p[p_col - 1]); # cut - 1 in matlab
            q_fft = int(q[q_row - 1] - q[q_col - 1]);
            C[row - 1, col - 1] = A[p_0 + p_fft, q_0 + q_fft];
    return C

class PWEM_2D():
    
    def Run(params, device, pwem_params):

        #######################################################################
        
        # Total number of spatial harmonics
        P,Q    = params.Harmonics
        M      = int(P*Q);
    
        # Bloch wave vectors
        bx     = pwem_params.beta[0,:];
        by     = pwem_params.beta[1,:];
        Nbx    = np.shape(bx);
    
        # Harmonic axes
        p      = np.arange( -np.floor(P/2), np.floor(P/2) + 1);
        q      = np.arange( -np.floor(Q/2), np.floor(Q/2) + 1);
    
        # Convolution matrices
        ERC    = device.ERC
        URC    = device.URC
    
        # Initialize normalized frequency arrays
        if params.Mode == 'E':
            W  = np.zeros((M,Nbx[0]));

        #######################################################################
    
        # Solve generalized eigen-value problem
        for nbeta in range(0,Nbx[0]):
            
            Kx     = bx[nbeta] - 2 * np.pi * p/params.Lx;
            Ky     = by[nbeta] - 2 * np.pi * q/params.Ly;
            Kx, Ky = np.meshgrid(Kx,Ky);
            Kx, Ky = Kx.flatten(), Ky.flatten();
            Kx, Ky = np.diag(Kx), np.diag(Ky);
    
          #######################################################################
    
            if params.Mode == 'E':
                if not params.is_magnetic:
                      A  = Kx**2 + Ky**2;
                else:
                  A  = Kx @ np.linalg.inv(URC) @ Kx + Ky @ np.linalg.inv(URC) @ Ky; # Operator for dielectric matrix
    
          #######################################################################
                
                k0           = eigh(A, ERC, eigvals_only = True);                   # Eigen values in general form (eigh (A,B))
                k0           = np.sort(k0)                                          # Sort eig vals (from lowest to highest)
                k0           = np.real(np.sqrt(k0)) / params.norm;                  # Normalize eig-vals
                W[:,nbeta] = k0;                                                    # Append eig-vals
        
          #######################################################################
          
            else:
                
                A = Kx @ np.linalg.inv(ERC) @ Kx + Ky @ np.linalg.inv(ERC) @ Ky;
                
                if not params.is_magnetic:
                    k0       = np.linalg.eigvals(A) 
                else: 
                    k0       = eigh(A, URC, eigvals_only = True); 
                k0           = np.sort(k0)
                k0           = np.real(np.sqrt(k0)) / params.norm;
                W[:,nbeta] = k0;
                  
          #######################################################################
    
        return W