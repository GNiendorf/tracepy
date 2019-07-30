import matplotlib.pyplot as plt
import numpy as np
from numpy import cos, sin, pi, sqrt

from .geometry import geometry
from .transforms import lab_frame
    
def _gen_points(surfaces):
    """ Generates the mesh points for each surface in the obj frame. """
    surfpoints = []
    for surface in surfaces:
        bound = surface.Diam/2.
        linspace = np.linspace(-bound, bound, 200) #General mesh points.
        x_mesh, y_mesh = np.meshgrid(linspace, -linspace)
        x_points = np.append(x_mesh, [linspace, np.zeros(200)]) #Used for cross sections.
        y_points = np.append(y_mesh, [np.zeros(200), linspace]) #Used for cross sections.
        meshpoints_2d = np.vstack((x_points.ravel(), y_points.ravel())).T
        z_points = surface.get_surface_plot(meshpoints_2d) #Function values in obj. frame.
        meshpoints = np.vstack((x_points.ravel(), y_points.ravel(), z_points)).T
        surfpoints.append(meshpoints) #Append to surfpoints which holds all surfaces' points.
    return np.array(surfpoints)

def _plot_rays(rays, axes, pltparams):
    """ Plots 2d ray history points. Takes list axes to specify axes (0, 1, 2) to plot. """
    for ray in rays:
        for idx,_ in enumerate(ray.P_hist[:-1]):
            F, G = ray.P_hist[idx][axes]
            F_p, G_p = ray.P_hist[idx+1][axes]
            H_p, I_p = ray.D_hist[idx+1][axes] #Alpha, beta, gamma rotations.
            plt.plot([G, G_p], [F, F_p], **pltparams)
        plt.plot([G_p, G_p+I_p],[F_p, F_p+H_p], **pltparams) #Plot direction of ray after stop. 

def _clip_lens(surfaces, surfpoints, idx):
    """ Clips points ouside of a lens intersection point. """
    surf1, surf2 = surfaces[idx], surfaces[idx+1]
    d = sqrt(np.sum(np.square(surf1.P - surf2.P)))
    points1, points2 = np.nan_to_num(surfpoints[idx]), np.nan_to_num(surfpoints[idx+1])
    points2[:,2] += d
    clipped_idx = (points2[:,2] - points1[:,2]) <= 0.
    surfpoints[idx][:,2][clipped_idx] = np.nan
    surfpoints[idx+1][:,2][clipped_idx] = np.nan  
    return surfpoints 
        
def _plot_surfaces(geo_params, axes):
    """ Plots 2d surface cross sections. Takes list axes to specify axes (0, 1, 2) to plot. """
    surfaces = [geometry(surf) for surf in geo_params]
    surfpoints = _gen_points(surfaces)
    lens_check = 0
    start = None
    for idx, surf in enumerate(surfaces):
        lens_condition = (idx+1 < len(surfaces) and
                            surfaces[idx].action == surfaces[idx+1].action == 'refraction')
        if lens_condition:
            surfpoints = _clip_lens(surfaces, surfpoints, idx)
        with np.errstate(invalid='ignore'):
            if np.any(np.mod(surf.D/pi, 1) != 0) and surf.c == 0 and surf.diam == 0:
                    cross_idx = abs(surfpoints[idx][:,axes[1]]) == 0 #Find cross section points.
            else:
                    cross_idx = abs(surfpoints[idx][:,1-axes[0]]) == 0 #Find cross section points.
        cross_points = surfpoints[idx][cross_idx]
        points = lab_frame(surf.R, surf, cross_points) #Transform to lab frame.
        F, G = points[:,axes[0]], points[:,axes[1]]
        #Connect the surfaces in a lens
        if surfaces[idx].action == surfaces[idx-1].action == 'refraction' and start is not None:
            lens_check = 1 - lens_check
            if lens_check == 1:
                start_2 = np.array([F[0], G[0]])
                end_2 = np.array([F[-1], G[-1]])
                dis1 = np.sqrt(np.sum(np.square(start - start_2)))
                dis2 = np.sqrt(np.sum(np.square(start - end_2)))
                if dis1 <= dis2:
                    idx = [0,-1]
                else:
                    idx = [-1,0]
                F = np.insert(F, idx, [start[0], end[0]])
                G = np.insert(G, idx, [start[1], end[1]])
        if lens_condition: #Store first and last point to connect surfaces.
            start = np.array([F[0], G[0]])
            end = np.array([F[-1], G[-1]])
        plt.plot(G, F, 'k')

def plotxz(geo_params, rays, pltparams={'c': 'red', 'alpha': 0.3 }, both=None):
    """ Plots the xz coordinates of all rays and surface cross sections. """
    rays = np.array([rayiter for rayiter in rays if rayiter.P is not None])
    if both is None: #Override 1,1,1 subplot if displaying side-by-side.
        plt.subplot(1,1,1, aspect='equal') #Keeps aspect ratio equal.
    _plot_rays(rays, [0, 2], pltparams)
    _plot_surfaces(geo_params, axes = [0,2])
    plt.xlabel("Z")
    plt.ylabel("X")
    
def plotyz(geo_params, rays, pltparams={'c': 'red', 'alpha': 0.3 }, both=None):
    """ Plots the yz coordinates of all rays and surface cross sections. """
    rays = np.array([rayiter for rayiter in rays if rayiter.P is not None])
    if both is None: #Override 1,1,1 subplot if displaying side-by-side.
        plt.subplot(1,1,1, aspect='equal') #Keeps aspect ratio equal.
    _plot_rays(rays, [1, 2], pltparams)
    _plot_surfaces(geo_params, axes = [1, 2])
    plt.xlabel("Z")
    plt.ylabel("Y")
    
def plot2d(geo_params, rays, pltparams={'c': 'red', 'alpha': 0.3 }):
    """ Plots both xz and yz side-by-side. """
    plt.subplot(2,1,1, aspect='equal')
    plotxz(geo_params, rays, pltparams, both=True)
    plt.subplot(2,1,2, aspect='equal') 
    plotyz(geo_params, rays, pltparams, both=True)
    plt.tight_layout()
