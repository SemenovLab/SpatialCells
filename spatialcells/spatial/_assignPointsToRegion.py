import numpy as np
import math
from shapely.geometry import Polygon, Point
from ._utils import *

def isInPolygon(polygons, point):
    """
    Check if a point is in any of the polygons
    """
    for polygon in polygons:
        if polygon.contains(Point(point)):
            return True
    return False

def assignPointsToRegion(anndata, boundaries_compoent, assigncolumn="region", target="target", donelist=[]):

    # edges
    boundary = getEdgesOnBoundaries(boundaries_compoent)
    # polygons
    polygons = getPolygons(boundaries_compoent)
    print("Assigning a region for each cell...")
    boundary_points = boundary[:, 0, :]
    x_min, x_max = math.floor(min(boundary_points[:,0]))-1, math.ceil(max(boundary_points[:,0]))
    y_min, y_max = math.floor(min(boundary_points[:,1]))-1, math.ceil(max(boundary_points[:,1]))
    nGrid = 10
    xv, yv, grid_label = getGrid(x_min, x_max, y_min, y_max, nGrid, boundary, polygons)
    assignRegion(xv, yv, grid_label, assigncolumn, target, boundary, polygons, anndata, donelist, nGrid = nGrid)
    print("Finish assigning.")
    return xv, yv, grid_label



def assignRegion(xv, yv, grid_label, assigncolumn, target, boundaries, polygons, adata, donelist, nGrid):
    '''
    Will change the "region" value as "target" if target is in the region given by boundaries
    @donelist: the values for region are ["0In", "1Bi", "2Bo", "3Out"]. If already 0In, we dont need to work on it.
    '''
    # print("assigning a region for each cell...", grid_label.shape, (xv).min(), (xv).max(), (yv).min(), (yv).max())
    for i in range(grid_label.shape[0]-1):
        for j in range(grid_label.shape[1]-1):
            
            label = sum([grid_label[i,j], grid_label[i, j+1], grid_label[i+1, j], grid_label[i+1, j+1]])
            # label = grid_label[i,j]
            if label == 0:# out
                pass
            elif label == 4: # in 
                #print("-"*nGrid,  label)
                cond = (adata.obs.X_centroid > xv[i,j]) & (adata.obs.X_centroid <= xv[i,j+1]) & (adata.obs.Y_centroid > yv[i,j]) & (adata.obs.Y_centroid <= yv[i+1,j])
                for item in donelist:
                    cond = cond & (adata.obs.region != item)
                adata.obs.loc[cond, assigncolumn] = target
                
            else: # partially in
                # find the points in the rectangle
                #print("-"*nGrid,  label, "nGrid:", nGrid)
                tmp = adata
                tmp = tmp[tmp.obs.X_centroid > xv[i,j]]
                tmp = tmp[tmp.obs.X_centroid <= xv[i,j+1]]

                tmp = tmp[tmp.obs.Y_centroid > yv[i,j]]
                tmp = tmp[tmp.obs.Y_centroid <= yv[i+1,j]]
                
                for item in donelist:
                    tmp = tmp[tmp.obs.region != item]

                IDs = tmp.obs["id"]
                tmp = tmp.obs[["X_centroid", "Y_centroid"]].to_numpy()
                newNGrid = 4 #(nGrid//2)+1
                # print("number of points:", tmp.shape, newNGrid)
                
                if (tmp.shape[0] < 100):
                    for idx in range(tmp.shape[0]):
                        # if (isInRegion(tmp[idx,:], boundaries)): #In
                        if (isInPolygon(polygons, tmp[idx,:])):
                            adata.obs.loc[adata.obs.id == IDs[idx], assigncolumn] = target
                            # adata.obs.loc[idx, assigncolumn] = target
                else:
                    xv1, yv1, grid_label1 = getGrid(xv[i,j], xv[i,j+1], yv[i,j], yv[i+1,j], newNGrid, boundaries, polygons)
                    assignRegion(xv1, yv1, grid_label1, assigncolumn, target, boundaries, polygons, adata, donelist, nGrid = newNGrid)



def getGrid(x_min, x_max, y_min, y_max, nGrids, edges, polygons):
    """
    Rectangle: x_min, x_max, y_min, y_max
    nGrids: the numbers of grids on the x direction
    Region: edges

    @output: the grid over the Region
    """
    step = (x_max - x_min) // nGrids
    # print("nGrid:", nGrid, "step:", step)
    points_x = [i for i in range(x_min, x_max + step + 1, step)]
    points_y = [i for i in range(y_min, y_max + step + 1, step)]
    points_x[-1] = x_max
    points_y[-1] = y_max
    xv, yv = np.meshgrid(points_x, points_y)

    grid_label = np.zeros(xv.shape)
    for i in range(xv.shape[0]):
        for j in range(xv.shape[1]):
            point = [xv[i, j], yv[i, j]]
            if hasEdge(point, step, polygons):
                grid_label[i, j] = 0.5
            # elif isInRegion(point, edges):
            elif isInPolygon(polygons, point):
                grid_label[i, j] = 1
    return xv, yv, grid_label
