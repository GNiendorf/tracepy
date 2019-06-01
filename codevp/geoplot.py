import matplotlib.pyplot as plt
import numpy as np
from numpy import cos, sin, pi, sqrt

from .geometry import geometry
from .transforms import lab_frame

class geoplot:
    """ Class for plotting all surfaces and propagated rays in a optical system. """
    def __init__(self, geo_params, rays=[], pltparams={'c': 'red', 'alpha': 0.5 }):
        self.rays = np.array([rayiter for rayiter in rays if rayiter.P is not None])
        self.geo_params = geo_params
        self.surfaces = [geometry(surf) for surf in self.geo_params]
        self.surfpoints = []
        self.pltparams = pltparams
        self.gen_points()
        self.surfpoints = np.array(self.surfpoints)
        
    def gen_points(self):
        """ Generates the mesh points for each surface in the obj frame. """
        for surface in self.surfaces:
            bound = 10 if surface.r==0 else surface.r
            linspace = np.linspace(-bound, bound, 200) #General mesh points.
            x_mesh, y_mesh = np.meshgrid(linspace, -linspace)
            x_points = np.append(x_mesh, [linspace, np.zeros(200)]) #Used for cross sections.
            y_points = np.append(y_mesh, [np.zeros(200), linspace]) #Used for cross sections.
            meshpoints_2d = np.vstack((x_points.ravel(), y_points.ravel())).T
            z_points = surface.get_surface(meshpoints_2d, plot=True) #Function values in obj. frame.
            meshpoints = np.vstack((x_points.ravel(), y_points.ravel(), z_points)).T
            self.surfpoints.append(meshpoints) #Append to surfpoints which holds all surfaces' points.
    
    def plot_rays(self, axes):
        """ Plots 2d ray history points. Takes list axes to specify axes (0, 1, 2) to plot. """
        for ray in self.rays:
            for idx,_ in enumerate(ray.P_hist[:-1]):
                F, G = ray.P_hist[idx][axes]
                F_p, G_p = ray.P_hist[idx+1][axes]
                H_p, I_p = ray.D_hist[idx+1][axes] #Alpha, beta, gamma rotations.
                plt.plot([F, F_p], [G, G_p], **self.pltparams)
            plt.plot([F_p, F_p+2*H_p],[G_p, G_p+2*I_p], **self.pltparams) #Plot direction of ray after stop. 
    
#     def clip_lens(self, idx):
#         """ Clips points ouside of a lens intersection point. """
#         surf1, surf2 = self.surfaces[idx], self.surfaces[idx+1]
#         d = sqrt(np.sum(np.square(surf1.P - surf2.P)))
#         points1, points2 = np.nan_to_num(self.surfpoints[idx]), np.nan_to_num(self.surfpoints[idx+1])
#         s1_dis, s2_dis = sqrt(np.sum(np.square(points1), axis=1)), sqrt(np.sum(np.square(points2), axis=1))
#         return (s2_dis - s1_dis) >= 1e-4
            
    def plot_surfaces(self, axes):
        """ Plots 2d surface cross sections. Takes list axes to specify axes (0, 1, 2) to plot. """
        for idx, surf in enumerate(self.surfaces):
#             if (idx+1 < len(self.surfaces) and
#                 self.surfaces[idx].inter == 'refraction' and self.surfaces[idx+1].inter == 'refraction'):
#                 #Clip surface points past lens intersection
#                 clipped_idx = self.clip_lens(idx)
#                 points = self.surfpoints[idx][clipped_idx]
#                 print(points[:10000])
#                 #self.surfpoints[idx][:,2] *= clipped_idx
#                 #self.surfpoints[idx+1][:,2] *= clipped_idx      
            if np.any(np.mod(surf.D/pi, 1) != 0) and surf.shape == 'plane' and surf.r2 == 0:
                with np.errstate(invalid='ignore'):
                    cross_idx = abs(self.surfpoints[idx][:,axes[1]]) == 0 #Find cross section points.
            else:
                with np.errstate(invalid='ignore'):
                    cross_idx = abs(self.surfpoints[idx][:,1-axes[0]]) == 0 #Find cross section points.
            cross_points = self.surfpoints[idx][cross_idx]
            points = lab_frame(surf.R, surf, cross_points) #Transform to lab frame.
            F, G = points[:,axes[0]], points[:,axes[1]]
            plt.plot(F, G, 'k')
    
    def plotxz(self, both=None):
        """ Plots the xz coordinates of all rays and surface cross sections. """
        if both is None: #Override 1,1,1 subplot if displaying side-by-side.
            plt.subplot(1,1,1, aspect='equal') #Keeps aspect ratio equal.
        self.plot_rays(axes = [0, 2])
        self.plot_surfaces(axes = [0,2])
        plt.xlabel("X")
        plt.ylabel("Z")
        
    def plotyz(self, both=None):
        """ Plots the yz coordinates of all rays and surface cross sections. """
        if both is None: #Override 1,1,1 subplot if displaying side-by-side.
            plt.subplot(1,1,1, aspect='equal') #Keeps aspect ratio equal.
        self.plot_rays(axes = [1, 2])
        self.plot_surfaces(axes = [1, 2])
        plt.xlabel("Y")
        plt.ylabel("Z")
        
    def plot2d(self):
        """ Plots both xz and yz side-by-side. """
        plt.subplot(1,2,1, aspect='equal')
        self.plotxz(both=True)
        plt.subplot(1,2,2, aspect='equal') 
        self.plotyz(both=True)
        plt.tight_layout()
