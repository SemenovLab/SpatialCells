import numpy as np
import math
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial import Delaunay


def estimateInitialDistance(X):
    '''Use hierarchical clustering to get the checkpoints'''
    model = AgglomerativeClustering(distance_threshold=0, n_clusters=None,compute_distances=True)
    model = model.fit(X)
    print("Computing distances...")
    nsamples = X.shape[0]
    newnodes = []
    layer = 0
    distances = []
    for i, item in enumerate(model.children_):
        # print("i="+str(i)+":", item[0], item[1],"d:", model.distances_[i], "-->",  nsamples+i)
        newnodes.append(nsamples+i)
        if (item[0] in newnodes) and (item[1] in newnodes) :
            # print("-------",layer, model.distances_[i])
            newnodes = []
            layer += 1
            distances.append(model.distances_[i])
    # print(min(model.distances_), max(model.distances_), np.median(model.distances_), np.mean(model.distances_))
    return distances


def PointsInCircum(eachPoint, r, n=100):
    '''
    Return n points within r distance from eachPoint
    '''
    return [(eachPoint[0] + math.cos(2*math.pi/n*x)*r, eachPoint[1] + math.sin(2*math.pi/n*x)*r) for x in range(0,n+1)]

def bufferPoints (inPoints, stretchCoef, n):
    '''
    Return n*len(inPoints) points that are within r distance of each point.
    '''
    newPoints = []
    for eachPoint in inPoints:
        newPoints += PointsInCircum(eachPoint, stretchCoef, n)
    newPoints = np.array(newPoints)

    return newPoints

def isInRegion(point, edges, xy):
    '''
    Check whether "point" in in the region given by edges and xy.
    @edges are indexes of the region
    @xy are the data points
    '''
    count = 0
    counted = {}
    for i, j in edges:
        if (xy[i, 0]  >= point[0] or xy[j, 0] >= point[0]): # one of x values is bigger than the point.x
            cond1 = xy[i, 1]  >= point[1] and xy[j, 1] <= point[1]
            cond2 = xy[i, 1]  <= point[1] and xy[j, 1] >= point[1]
            if cond1 or cond2:
                if xy[i, 1] in counted.keys():
                    y = xy[i, 1]
                    y2 = xy[counted[xy[i, 1]], 1]
                    if y > min(xy[j, 1], y2) and y < max(xy[j, 1], y2):
                        continue
                if xy[j, 1] in counted.keys():
                    y = xy[j, 1]
                    y2 = xy[counted[xy[j, 1]], 1]
                    if y > min(xy[i, 1], y2) and y < max(xy[i, 1], y2):
                        continue
                count += 1
                counted[xy[i, 1]] =j
                counted[xy[j, 1]] =i
                
    return (count % 2) != 0



def getRegionArea(edges, xy):
    indexes = []
    for i, j in edges:
        indexes.append(i)
        indexes.append(j)
    points = xy[indexes, :]

    c = (int(max(points[:, 0]) + min(points[:, 0]))//2, 
     int(max(points[:, 1]) + min(points[:, 1]))//2)
    print("center:", c)

    from shapely.geometry import Polygon
    area = 0
    for i, j in edges:
        poly = Polygon([xy[i], xy[j], c])
        area += poly.area
    return int(area)

def getSorted(sub_li):
    ## [(i,j), (i,j), ...]
    ## sorted based on the minimum of the two
    sub_li = list(sub_li)
    sub_li.sort(key = lambda x: min(x))
    return sub_li

def getEdgeComponent(edges):
    '''
    The output is a dictionary. 
    '''
    edges = getSorted(edges)
    
    idx = 0
    indexes_visited = set()
    edge_components = {}
    # 1: [(i,j,k), [(i,j), (i, k)]]
    
    for i, j in edges:
        # new component
        if (i not in indexes_visited) and (j not in indexes_visited):
            edge_components[idx] = [set([i,j]), [(i,j)]]
            idx += 1
        else:
            for key, value in edge_components.items():
                if i in value[0] or j in value[0]:
                    value[0].add(i)
                    value[0].add(j)
                    value[1].append((i,j))
                    break
        indexes_visited.add(i)
        indexes_visited.add(j)
    
    return edge_components

def collapse(edge_components):
    for target_key in sorted(range(len(edge_components)),reverse=True):
        target_value = edge_components[target_key]
        if len(target_value[0]) == 0:
            continue

        collapse_set = set()
        collapse_set.add(target_key)
        for pi in target_value[0]:
            for idx, value in edge_components.items():
                if target_key != idx and pi in value[0]:
                    collapse_set.add(idx)
        
        if len(collapse_set) > 1:
            #print(target_key, ":", collapse_set)

            key_new = min(collapse_set)
            value_new = edge_components[key_new]
            for i in collapse_set:
                if i != key_new and len(edge_components[i]) >0:
                    #print(i, "--", len(value_new[0]), len(value_new[1]))

                    for j in edge_components[i][0]:
                        value_new[0].add(j)
                    value_new[1] += edge_components[i][1]

                    edge_components[i][0] = set()
                    edge_components[i][1] = []
                    #print(len(value_new[0]), len(value_new[1]))
    return edge_components

def removeEdgeComponents(edge_components, 
                         nedges_min = 50,
                         nedges_out_min = None,
                         xy = None):
    
    # only keep the components with >= nedges_min
    edges_new = []
    for key, value in edge_components.items():
        if len(value[1]) >= nedges_min: 
            edges_new += value[1]
            
    # if out side the main tumor region, we may want to keep
    if (nedges_out_min != None):
        for key, value in edge_components.items():
            if len(value[1]) < nedges_min and len(value[1]) >= nedges_out_min: 
                if isInRegion(xy[value[1][0][0],:], edges_new, xy):
                    pass
                else:
                    edges_new += value[1]

    return (edges_new)


def get_alphas(points):
    """
    Compute alpha candiates for the alpha shape (concave hull) of a set of points.
    :param points: np.array of shape (n,2) points.
    :return: a list of alpha candiates
    """
    assert points.shape[0] > 3, "Need at least four points"
    tri = Delaunay(points) #triangulation
    ret = set()
    # Loop over triangles:
    # ia, ib, ic = indices of corner points of the triangle
    for ia, ib, ic in tri.vertices:
        pa = points[ia]
        pb = points[ib]
        pc = points[ic]
        # Computing radius of triangle circumcircle
        # www.mathalino.com/reviewer/derivation-of-formulas/derivation-of-formula-for-radius-of-circumcircle
        a = np.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = np.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = np.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)
        s = (a + b + c) / 2.0
        area = np.sqrt(s * (s - a) * (s - b) * (s - c))
        circum_r = a * b * c / (4.0 * area)
        ret.add(int(circum_r))
    sorted(ret)
    return ret

def alpha_shape(points, alpha, only_outer=True):
    """
    Compute the alpha shape (concave hull) of a set of points.
    :param points: np.array of shape (n,2) points.
    :param alpha: alpha value.
    :param only_outer: boolean value to specify if we keep only the outer border
    or also inner edges.
    :return: set of (i,j) pairs representing edges of the alpha-shape. (i,j) are
    the indices in the points array.
    """
    assert points.shape[0] > 3, "Need at least four points"
    def add_edge(edges, i, j):
        """
        Add an edge between the i-th and j-th points,
        if not in the list already
        """
        if (i, j) in edges or (j, i) in edges:
            # already added
            assert (j, i) in edges, "Can't go twice over same directed edge right?"
            if only_outer:
                # if both neighboring triangles are in shape, it's not a boundary edge
                edges.remove((j, i))
            return
        edges.add((i, j))
    tri = Delaunay(points)
    edges = set()
    # Loop over triangles:
    # ia, ib, ic = indices of corner points of the triangle
    for ia, ib, ic in tri.vertices:
        pa = points[ia]
        pb = points[ib]
        pc = points[ic]
        # Computing radius of triangle circumcircle
        # www.mathalino.com/reviewer/derivation-of-formulas/derivation-of-formula-for-radius-of-circumcircle
        a = np.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = np.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = np.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)
        s = (a + b + c) / 2.0
        area = np.sqrt(s * (s - a) * (s - b) * (s - c))
        circum_r = a * b * c / (4.0 * area)
#         print(circum_r, alpha)
        if circum_r < alpha:
            add_edge(edges, ia, ib)
            add_edge(edges, ib, ic)
            add_edge(edges, ic, ia)
    return edges

def getGrid(x_min, x_max, y_min, y_max, nGrids, edges_tmp, points_tmp):
    '''
    Rectangle: x_min, x_max, y_min, y_max
    nGrids: the numbers of grids on the x direction
    Region: edges_tmp, points_tmp
    
    @output: the grid over the Region
    '''
    step = (x_max - x_min)//nGrids
    # print("nGrid:", nGrid, "step:", step)
    points_x = [i for i in range(x_min, x_max+step+1, step)]
    points_y = [i for i in range(y_min, y_max+step+1, step)]
    points_x[-1] = x_max
    points_y[-1] = y_max
    xv, yv = np.meshgrid(points_x, points_y)
    
    grid_label = np.zeros(xv.shape)
    for i in range(xv.shape[0]):
        for j in range(xv.shape[1]):
            point = [xv[i,j],yv[i,j]]
            if isInRegion(point, edges_tmp, points_tmp):
                grid_label[i,j] = 1
    return xv, yv, grid_label

def assignRegion(xv, yv, grid_label, assigncolumn, target, edges_tmp, points_tmp, adata, donelist, nGrid):
    '''
    Will change the "region" value as "target" if target is in the region (edges_tmp, points_tmp)
    @donelist: the values for region are ["0In", "1Bi", "2Bo", "3Out"]. If already 0In, we dont need to work on it.
    '''
    print("assigning a region for each cell...", grid_label.shape)
    for i in range(grid_label.shape[0]-1):
        for j in range(grid_label.shape[1]-1):
            
            label = sum([grid_label[i,j], grid_label[i, j+1], grid_label[i+1, j], grid_label[i+1, j+1]])
            
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
                newNGrid = nGrid #(nGrid//2)+1
                # print("number of points:", tmp.shape, newNGrid)
                
                if (tmp.shape[0] < 200) or (newNGrid < 2) :
                    for idx in range(tmp.shape[0]):
                        if (isInRegion(tmp[idx,:], edges_tmp, points_tmp)): #In
                            adata.obs.loc[adata.obs.id == IDs[idx], assigncolumn] = target
                else:
                    xv1, yv1, grid_label1 = getGrid(xv[i,j], xv[i,j+1], yv[i,j], yv[i+1,j], newNGrid, edges_tmp, points_tmp)
                    assignRegion(xv1, yv1, grid_label1, assigncolumn, target, edges_tmp, points_tmp, adata, donelist, nGrid = newNGrid)



