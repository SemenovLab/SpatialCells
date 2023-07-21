from sklearn.cluster import DBSCAN

from spacells_utils import *

def getCommunities(anndata, markerList):
    ad_temp = anndata[anndata.obs[markerList].all(axis=1)].copy()
    X = ad_temp.obs[["X_centroid", "Y_centroid"]].to_numpy()
    db = DBSCAN(eps = 200, min_samples = 20, algorithm="kd_tree").fit(X)
    labels = db.labels_
    ad_temp.obs["community"] = labels
    return ad_temp

def getBoundary(anndata, communityIndexList, alpha=400, nedges_min=50, nedges_out_min=None):
    """
    Get the boundary of the community
    @anndata: the anndata object
    @communityIndexList: the list of community indexes
    @alpha: the alpha parameter for alpha shape
    @nedges_min: the minimum number of edges for a community
    @nedges_out_min: the minimum number of edges for a community to be considered as "out"

    Return: the boundary of the community
    """
    xy = anndata.obs[anndata.obs["community"].isin(communityIndexList)][
        ["X_centroid", "Y_centroid"]
    ].to_numpy()
    print("xy shape:", xy.shape)
    edge_components = getAlphaShapes(xy, alpha)
    print("edge_components:", len(edge_components), [len(i) for i in edge_components], edge_components[0].shape)
    grouped_components = groupRemoveEdgeComponents(edge_components, nedges_min, nedges_out_min)
    return grouped_components

def getEdgesOnBoundary(boundary):
    points = []
    for comp in boundary:
        for ring in comp:
            print(ring.shape)
            points.append(ring)
    return np.concatenate(points)

def getRingsOnBoundary(boundary):
    rings = []
    for comp in boundary:
        for ring in comp:
            rings.append(ring[:,0,:])
    return rings

def getExtendedBoundary(boundary, offset=400, alpha=400):
    boundary_edges = getEdgesOnBoundary(boundary)
    boundary_points = boundary_edges[:,0,:]
    pointsStretched = bufferPoints(boundary_points, offset, n=100)
    edge_components = getAlphaShapes(pointsStretched, alpha)
    extended_boundaries = []
    for comp in edge_components:
        if not isInRegion(comp[0][0], boundary_edges):
            extended_boundaries.append(comp)
    return np.concatenate(extended_boundaries)

def getshrunkenBoundary(boundary, offset=400, alpha=400):
    boundary_edges = getEdgesOnBoundary(boundary)
    boundary_points = boundary_edges[:,0,:]
    pointsStretched = bufferPoints(boundary_points, offset, n=100)
    edge_components = getAlphaShapes(pointsStretched, alpha)
    extended_boundaries = []
    for comp in edge_components:
        if isInRegion(comp[0][0], boundary_edges):
            extended_boundaries.append(comp)
    return np.concatenate(extended_boundaries)

def assignPointsToRegion(anndata, boundary, assigncolumn="region", target="target", donelist=[]):
    # boundary_edges = getEdgesOnBoundary(boundary)
    boundary_edges = boundary
    boundary_points = boundary_edges[:, 0, :]
    x_min, x_max = math.floor(min(boundary_points[:,0]))-1, math.ceil(max(boundary_points[:,0]))
    y_min, y_max = math.floor(min(boundary_points[:,1]))-1, math.ceil(max(boundary_points[:,1]))
    nGrid = 10
    xv, yv, grid_label = getGrid(x_min, x_max, y_min, y_max, nGrid, boundary_edges)
    assignRegion(xv, yv, grid_label, assigncolumn, target, boundary_edges, anndata, donelist, nGrid = nGrid)


def isInRegion(point, boundary):
    """
    Check whether "point" in in the region given by edges and xy.
    @edges are indexes of the region
    @xy are the data points
    """
    count = 0
    larger_x2_count = 0
    smaller_x2_count = 0

    point = np.array(point)
    for pi, pj in boundary:
        # if point is on the edge, return True
        if (point[0] == pi[0] and point[1] == pi[1]) or (point[0] == pj[0] and point[1] == pj[1]):
            return True
        # handles edge case of touching line segments / segments on a straight line
        # if either x == point.x, let this point be x1.
        if pi[0] == point[0] or pj[0] == point[0]:
            if (pi[0] + pj[0]) / 2 > point[0]:  # x2 larger
                larger_x2_count += 1
            elif (pi[0] + pj[0]) / 2 < point[0]:  # x2 smaller
                smaller_x2_count += 1
        # check if line segment intersects with point, (0, point.y). if so count += 1
        # only check if the line segment is on the left side of the point and point between the two y values
        if (pi[1] < point[1]) != (pj[1] < point[1]) and (pi[0] < point[0] or pj[0] < point[0]):
            count += isIntersect(point, (0, point[1]), pi, pj)

    count += larger_x2_count % 2 and smaller_x2_count % 2
    return (count % 2) != 0


def assignRegion(xv, yv, grid_label, assigncolumn, target, boundaries, adata, donelist, nGrid):
    '''
    Will change the "region" value as "target" if target is in the region (edges_tmp, points_tmp)
    @donelist: the values for region are ["0In", "1Bi", "2Bo", "3Out"]. If already 0In, we dont need to work on it.
    '''
    print("assigning a region for each cell...", grid_label.shape, (xv).min(), (xv).max(), (yv).min(), (yv).max())
    for i in range(grid_label.shape[0]-1):
        for j in range(grid_label.shape[1]-1):
            
            # label = sum([grid_label[i,j], grid_label[i, j+1], grid_label[i+1, j], grid_label[i+1, j+1]])
            label = grid_label[i,j]
            if label == 0:# out
                pass
            elif label == 1: # in 
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

                # IDs = tmp.obs["id"]
                tmp = tmp.obs[["X_centroid", "Y_centroid"]].to_numpy()
                newNGrid = 4 #(nGrid//2)+1
                # print("number of points:", tmp.shape, newNGrid)
                
                if (tmp.shape[0] < 100):
                    for idx in range(tmp.shape[0]):
                        if (isInRegion1(tmp[idx,:], boundaries)): #In
                            adata.obs.loc[idx, assigncolumn] = target
                else:
                    xv1, yv1, grid_label1 = getGrid(xv[i,j], xv[i,j+1], yv[i,j], yv[i+1,j], newNGrid, boundaries)
                    assignRegion(xv1, yv1, grid_label1, assigncolumn, target, boundaries, adata, donelist, nGrid = newNGrid)






def getGrid(x_min, x_max, y_min, y_max, nGrids, edges):
    """
    Rectangle: x_min, x_max, y_min, y_max
    nGrids: the numbers of grids on the x direction
    Region: edges_tmp, points_tmp

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
            if hasEdge(point, step, edges):
                grid_label[i, j] = 0.5
            elif isInRegion1(point, edges):
                grid_label[i, j] = 1
    return xv, yv, grid_label