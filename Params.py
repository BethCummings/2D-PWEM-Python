# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 12:29:22 2023

@author: Ignatio Denton
"""
import numpy as np
#######################################################################

class params():

  # initialize class objects
  def __init__(self):
    self.Lx          = 1;                  # period in x dir
    self.Ly          = 1;                  # period in y dir
    self.r           = 0.35;               # circle diameter
    self.ax          = 1;                # ellipse length (x axis)
    self.ay          = 1;               # ellipse width (y axis)
    self.er1         = 1;               # hole epsilon
    self.er2         = 2.25;                  # material epsilon
    self.Mode        = 'E';                # or H mode
    self.Harmonics   = [11, 11];        # harmonics in x, y, z dirs or P, Q, R
    self.dim         = [1024, 1024];    # Nx, Ny, Nz
    self.norm        = 2 * np.pi/self.Lx;  # normalization constant
    self.is_magnetic = 0;                  # magnetic constant
    self.x           = np.linspace(-self.Lx/2, self.Lx/2,self.dim[0])
    self.y           = np.linspace(-self.Ly/2, self.Ly/2,self.dim[1])

#######################################################################

class pwem_params():

    # initialize class objects
    def __init__(self, params):
        
      self.N_Points  = int(50);
      self.BC        = 0;
      self.beta      = []
      
      if self.BC == 1 or 2 or 3:
          self.T1    = 2*np.pi/params.Lx * np.array([[1],[0]])
          self.T2    = 2*np.pi/params.Ly * np.array([[0],[1]])
      else:
          component  = [1/params.Lx, -1/(params.Ly * np.sqrt(3))]
          self.T1    = 2*np.pi * np.array( ( component ) ) 
          self.T2    = 2*np.pi * np.array( ( np.abs(component) ) )  
    
    # Key points of symmetry
    def Symmetry(self):
        
        if self.BC == 1:
            
        #######################################################################
        
            Gamma =  np.array([[0],[0]])
            X     = self.T1/2
            M     = (self.T2 + self.T1)/2
            
        #######################################################################
        
            bx1   = np.linspace(Gamma[0], X[0], self.N_Points)
            bx2   = np.linspace(X[0],     X[0], self.N_Points)
            bx3   = np.linspace(M[0], Gamma[0], self.N_Points)
            
        #######################################################################
        
            by1   = np.linspace(Gamma[1], Gamma[1], self.N_Points)
            by2   = np.linspace(X[1],         M[1], self.N_Points)
            by3   = np.linspace(M[1],     Gamma[1], self.N_Points)
            
            bx    = np.concatenate((bx1,bx2, bx3))
            by    = np.concatenate((by1,by2, by3))
            
            self.beta = np.stack((bx,by), axis = 0) 
            
        elif self.BC == 2:
            
        #######################################################################
            Gamma =  np.array([[0],[0]])
            X     =  self.T1/2
            M     = (self.T1 + self.T2)/2
            Delta =  self.T2/2
        #######################################################################
        
            bx1   = np.linspace(Gamma[0], X[0], self.N_Points)
            bx2   = np.linspace(X[0],     X[0], self.N_Points)
            bx3   = np.linspace(M[0], Gamma[0], self.N_Points)
            bx4   = np.linspace(Gamma[0],Gamma[0],self.N_Points)
            
        #######################################################################
        
            by1   = np.linspace(Gamma[1], Gamma[1], self.N_Points)
            by2   = np.linspace(X[1],         M[1], self.N_Points)
            by3   = np.linspace(M[1],         M[1], self.N_Points)
            by4   = np.linspace(M[1],     Gamma[1], self.N_Points)
            
            bx    = np.concatenate((bx1,bx2, bx3, bx4))
            by    = np.concatenate((by1,by2, by3, by4))
            
            self.beta = np.stack((bx,by), axis = 0) 
        
        elif self.BC == 3:
            
            "Half Square lattice - under construction"
            Delta = np.array([0,0])
            Gamma = 1/2*np.array([0,-np.pi*2/params.Ly])
            Sigma = -(self.T1 + self.T2)/2
            M     = self.T1/2
            Z     = (self.T1 + self.T2)/2
            X     = self.T2/2
            
            
            
            pass
        elif self.BC == 4:
            "Hexagonal lattice with circle inside"
            
            pass

        #######################################################################
        

            
        else:
            
        #######################################################################
            
            bx    = np.linspace(-np.pi,  np.pi, self.N_Points)
            by    = np.linspace( np.pi, -np.pi, self.N_Points)
            beta  = np.zeros((2, self.N_Points**2))
            
            idx   = int(0);
            
            for nx in range(0, self.N_Points):
                for ny in range(0, self.N_Points):
                    beta[:,idx] = np.array(([bx[nx], by[ny]]))
                    idx = idx + 1

        #######################################################################
        
            self.beta = np.array(beta)
    
        return self