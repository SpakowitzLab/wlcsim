"""lines_in_box monte carlo simulator
pseudocode for col checker
poca - two points at which infinite lines get closest
bpoca - poca of two finite segments defining centers of cylinders
loca - line conencts poca
for periodization old line
  for periodization of new line
      if d(poca) > 2r
          return false
      if bpoca is not on tip of cylinder for either cylinder
          return naive formula
      elif one tip and one interior point (of closest approach)
          if loca's midpoint is inside box # already know: |loca| < 2 r
                  or check if intersections with box bdry collide
              return true
          else:
              return false
      else (if two tips)
          for each tip
              for each edge adjacent to the face tip is hitting
                  if |loca btween edge and cylinder| < r
                     mark wall opposite that edge for that tip
          for each wall that both tips share
              return if their ellipse outlines on that wall collide inside of the box wall

Yes I know this could have been done more prettily with object oriented shapes,
but I'm busy getting real results, do it yourself if you want.

"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import pandas as pd
from ..utils import path as wpath
import os
from numba import jit
import multiprocessing
import pscan
from scipy.optimize import curve_fit
import statsmodels.api as sm
from scipy.spatial import cKDTree


# our "floating precision" cutoff
SMALL_NUM = np.power(10.0, -8)

# in what follows, a cylinder's centerline's is used as a proxy for the
# cylinder itself. comments often only make sense with this in mind (e.g. the
# cylinder "hits" a wall means the intersection of its centerline with the
# plane defining the wall).

# if a cylinder hits one of the infinite planes defined by the ws, w0s below at
# a distance far enough away from the unit box's extents, then the outline of
# the intersection of the full cylinder with the unit box's wall will be
# effectively a rectange instead of an ellipse. this is the maximum magnitude
# of the collision (in units of box-widths away from the box's center) before the
# collision code simply uses the parallelogram approximation
max_boxes_away_for_ellipse_collision = 20

# now comes an absurd amount of precomputed constants that describe a unit box
# in the first quadrant

# normal vectors and intercepts defining the planes of the unit box in the
# first quadrant, used throughout collision detection
ws = [[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0],
      [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
ws = [np.array(w) for w in ws]
# three planes pass through (0,0,0), other three pass through (1,1,1)
w0s = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
# how to project a vector into any of the box walls
plane_projs = [np.array(w) == 0.0 for w in ws]
non_zero_indices = [[i for i,c in enumerate(plane_proj) if c]
                    for plane_proj in plane_projs]
plane_projs = [p.astype(float) for p in plane_projs]
# convenience variables for defining points on the box
ux = np.array([1.0, 0.0, 0.0])
uy = np.array([0.0, 1.0, 0.0])
uz = np.array([0.0, 0.0, 1.0])
uo = np.array([0.0, 0.0, 0.0])
edges = [(uo, ux), (ux, ux+uz), (ux+uz, uz), (uz, uo),
         (uo, uy), (uy, uy+uz), (uy+uz, uz),
         (ux, ux+uy), (ux+uy, uy),
         (ux+uy, ux+uy+uz), (ux+uz, ux+uy+uz), (uy+uz, ux+uy+uz)]
# I didn't choose the best edge labelling, but if you want to change it, just
# change all the relevant variables up here in the global scope, I didn't use
# the edge labeling explicitly anywhere in the main program, just through these
# variables

# each edge is adjacent to two faces
adjacent_faces = [(0, 2), (0, 4), (0, 5), (0, 1),
                  (1, 2), (1, 3), (1, 5), (2, 4),
                  (2, 3), (3, 4), (4, 5), (3, 5)]
# a map from faces to clockwise (from inside box) listing of adjacent edges to
# that face, list of lists of edge indexes.
adjacent_edges = [[0, 1, 2, 3],  # face 0
                  [4, 5, 6, 3],  # face 1
                  [0, 4, 8, 7],  # face 2
                  [11, 5, 8, 9], # face 3
                  [1, 10, 9, 7], # face 4
                  [2, 6, 11, 10]] # face 5
# face, edge -> other face
face_edge_to_face = {(faces[0], edge): faces[1] for edge, faces
        in enumerate(adjacent_faces)}
face_edge_to_face.update({(faces[1], edge): faces[0] for edge, faces
        in enumerate(adjacent_faces)})
# set of pairs of faces that share an edge, with both index orders for fast
# lookup
face_pairs = {face_pair for face_pair in adjacent_faces}
face_pairs.update({(face_pair[1], face_pair[0]) for face_pair in adjacent_faces})
# lookup the opposite of a face quickly as face_opposite[face]
face_opposite = [3, 4, 5, 0, 1, 2]
# periodization vectors to add to get possible periodizations about a
# particular face can be looked up easily here
# we quickly use the normals of the faces that share an edge with the face we
# want to keep constant during our periodization
periodizers_given_face = [
        [ws[face_edge_to_face[(face,edge)]] for edge in adjacent_edges[face]]
        for face in range(6)]
for plist in periodizers_given_face:
    plist.append(plist[0] + plist[1])
    plist.append(plist[1] + plist[2])
    plist.append(plist[2] + plist[3])
    plist.append(plist[3] + plist[0])

# semantically useful throughout code, create once
box_center = np.array([0.5, 0.5, 0.5])

def new_object(shape_params, shape):
    if 'sphere' in shape:
        r = shape_params
        return (uniform_point_from_unit_cube(), r)
    elif 'cylinder' in shape:
        return uniform_segment_from_unit_cube()
    else:
        raise ValueError('Invalid shape requested!')
    return None

def volume_of_shape(obj, shape_params, shape, periodic):
    if 'sphere' in shape:
        (c, r) = obj
        return volume_of_sphere(c, r, periodic)
    elif 'cylinder' in shape:
        xi, xf = obj
        d = shape_params
        return volume_of_cylinder(xi, xf, d, periodic)
    else:
        raise ValueError('Invalid shape requested!')
    return np.nan

def detect_collision(obj, obj_new, shape_params, shape, periodic=False):
    if 'sphere' in shape:
        (c, r) = obj
        (c_new, r_new) = obj_new
        return detect_collision_spheres(c, r, c_new, r_new, periodic)
    elif 'cylinder' in shape:
        xi, xf = obj
        xi_new, xf_new = obj_new
        d = shape_params
        return detect_collision_cylinders(xi, xf, xi_new, xf_new, d, periodic)
    else:
        raise ValueError('Invalid shape requested!')
    return False

def columns_for_shape(shape):
    if 'cylinder' in shape:
        return ['xi1', 'xi2', 'xi3', 'xf1', 'xf2', 'xf3', 'phi', 'num_objects',
                'did_succeed', 'is_removal']
    elif 'sphere' in shape:
        return ['c1', 'c2', 'c3', 'r', 'phi', 'num_objects',
                'did_succeed', 'is_removal']
    else:
        raise ValueError('Invalid shape requested!')
    return []

def numeric_columns_for_shape(shape):
    if 'sphere' in shape:
        return ['c1', 'c2', 'c3', 'r', 'phi', 'num_objects']
    elif 'cylinder' in shape:
        return ['xi1', 'xi2', 'xi3', 'xf1', 'xf2', 'xf3', 'phi', 'num_objects']
    else:
        raise ValueError('Invalid shape requested!')
    return []

def row_for_shape(obj, shape, phi, num_objects, did_succeed, is_removal):
    if 'sphere' in shape:
        if obj is None:
            c = [np.nan, np.nan, np.nan]
            r = np.nan
        else:
            c, r = obj
        return [c[0], c[1], c[2], r, phi, num_objects, did_succeed, is_removal]
    elif 'cylinder' in shape:
        if obj is None:
            xf = [np.nan, np.nan, np.nan]
            xi = [np.nan, np.nan, np.nan]
        else:
            xi, xf = obj
        return [xi[0], xi[1], xi[2], xf[0], xf[1], xf[2], phi, num_objects,
                did_succeed, is_removal]
    else:
        raise ValueError('Invalid shape requested!')
    return []

def get_removal_amount(obj, shape):
    """Get parameter that is used to scale probability of removal.
    e.g. for cylinders, this is their length. for spheres, this is their
    radius. etc."""
    if 'sphere' in shape:
        return 1 # all have same volume, so do per number
    elif 'cylinder' in shape:
        xi, xf = obj
        return np.linalg.norm(xf - xi)
    else:
        raise ValueError('Invalid shape requested!')
    return np.nan

@jit(nopython=True)
def line_hit_hyperplane(xi, xf, w, w0):
    """Get point of intersection of a lines in N-D with a hyperplane of
    dimension N-1. The line should be defined by two pionts and the hyperplane
    should be defined by its normal vector and a point on it."""
    w = w/np.linalg.norm(w)
    x = xf - xi
    # projection onto plane's normal of x
    x_proj_w = np.dot(w, x)
    if np.abs(x_proj_w) < SMALL_NUM: # effectively parallel to plane
        return None
    # intersection point as parameter of xi + (xf - xi)*t
    t = (w0 - np.dot(xi, w))/x_proj_w
    # intersection point itself
    return xi + x*t

@jit(nopython=True)
def get_face_given_point(x):
    # xz-plane
    if np.abs(x[1] - 0.0) < SMALL_NUM:
        return 0
    # yz-plane
    elif np.abs(x[0] - 0.0) < SMALL_NUM:
        return 1
    elif np.abs(x[2] - 0.0) < SMALL_NUM:
        return 2
    elif np.abs(x[1] - 1.0) < SMALL_NUM:
        return 3
    elif np.abs(x[0] - 1.0) < SMALL_NUM:
        return 4
    elif np.abs(x[2] - 1.0) < SMALL_NUM:
        return 5

@jit(nopython=True)
def d_line_line(x1, x2, y1, y2):
    """Find distance between two lines defined by two non-degenerate
    points on each of those lines."""
    w = x1 - y1
    u = x2 - x1
    v = y2 - y1
    uu = np.dot(u,u)
    uv = np.dot(u,v)
    vv = np.dot(v,v)
    uw = np.dot(u,w)
    vw = np.dot(v,w)
    D = uu*vv - uv*uv # "descrimant" / denominator ( always > 0 )

    if D < SMALL_NUM: # lines must be nearly parallel
        # subtract projection of w onto either line from w
        # same as distance formula from just one line to a point
        s = 0.0
        t = uw/uu if uu > vv else vw/vv
    else:
        s = (uv*vw - vv*uw)/D
        t = (uu*vw - uv*uw)/D
    # vector of closest approach
    voca = w + s*u - t*v
    return np.linalg.norm(voca)


@jit(nopython=True)
def d_segment_segment(xinit0, xfinal0, xinit1, xfinal1):
    """Find distance between two segments defined by initial
    and ending position vectors."""
    # see e.g.
    # http://geomalgorithms.com/a07-_distance.html#dist3D_Segment_to_Segment()
    # for details of algorithm
    u = xfinal0 - xinit0
    v = xfinal1 - xinit1
    w = xinit0 - xinit1
    a = np.dot(u,u)
    b = np.dot(u,v)
    c = np.dot(v,v)
    d = np.dot(u,w)
    e = np.dot(v,w)
    D = a*c - b*b
    sc = D; sN = D; sD = D;
    tc = D; tN = D; tD = D;

    if D < SMALL_NUM:
        sN = 0.0         # force using point P0 on segment S1
        sD = 1.0         # to prevent possible division by 0.0 later
        tN = e
        tD = c
    else:                 # get the closest points on the infinite lines
        sN = (b*e - c*d)
        tN = (a*e - b*d)
        if sN < 0.0:        # sc < 0 => the s=0 edge is visible
            sN = 0.0
            tN = e
            tD = c
        elif sN > sD:  # sc > 1  => the s=1 edge is visible
            sN = sD
            tN = e + b
            tD = c
    if tN < 0.0:            # tc < 0 => the t=0 edge is visible
        tN = 0.0
        # recompute sc for this edge
        if -d < 0.0:
            sN = 0.0
        elif -d > a:
            sN = sD
        else:
            sN = -d
            sD = a
    elif tN > tD:      # tc > 1  => the t=1 edge is visible
        tN = tD
        # recompute sc for this edge
        if (-d + b) < 0.0:
            sN = 0
        elif (-d + b) > a:
            sN = sD
        else:
            sN = (-d +  b)
            sD = a
    # finally do the division to get sc and tc
    sc = sN / sD if np.abs(sN) > SMALL_NUM else 0.0
    tc = tN / tD if np.abs(tN) > SMALL_NUM else 0.0

    # get the difference of the two closest points
    dP = w + (sc * u) - (tc * v)  # =  S1(sc) - S2(tc)

    return np.linalg.norm(dP)

@jit(nopython=True)
def norm3_squared(x):
    """Avoid extra operations when doing something like np.dot(f(x), f(x))"""
    # we don't need the generality of np.linalg.norm, so @jit of this simple function
    # will probably outperform np
    return x[0]*x[0] + x[1]+x[1] + x[2]*x[2]
    # return (x[0] - y[0])**2 + (x[1] - y[1])**2 + (x[2] - y[2])**2

def periodic_cylinder_collision(ui, uf, vi, vf, d):
    """Decide whether or not two cylinders inside of the box of side length 1
    are overlapping each other, when the box has "periodic" boundary conditions
    along all sides except for those on which the cylinders terminate."""

#     # renormalize to the equivalent problem on unit cube
#     if a != 1:
#         d = d/float(a)
#         ui = ui/a
#         vi = vi/a
#         uf = uf/a
#         vf = vf/a
#         a = 1.0
    r = d/2.0
    u = uf - ui
    v = vf - vi
    a = np.dot(u,u)
    b = np.dot(u,v)
    c = np.dot(v,v)
    d = np.dot(u,w)
    e = np.dot(v,w)
    D = a*c - b*b

    # if lines are nearly parallel, don't try to find point they're closest
    if D < SMALL_NUM:
        # any point works to check for collisions
        return norm3_squared(ui - vi) < d*d
    # calculate fractional path from xinit to xfinal that point of closest
    # approach between the two infinite lines defined by the input points is
    sN = (b*e - c*d) # numerator of "sC" on geomalgorithms.com
    tN = (a*e - b*d) # numerator of "tC"
    # if that point of closest approach is within the box (i.e. 0-1.0
    # fractionally), then collision detection is straightforward
    if 0.0 < sN and sN < D and 0.0 < tN and tN < D:
        sC = sN / sD if np.abs(sN) > SMALL_NUM else 0.0
        tC = tN / tD if np.abs(tN) > SMALL_NUM else 0.0
        return norm3_squared(ui + (sC*u) - (vi + (tC*v))) < d*d
    # if the point of closest approach is described by the tip of one
    # cylinder's center and an interior point of the other cylinder's center,
    # then there is a collision IFF (midpoint of line joining closest approach
    # points is inside of the box OR their extents on the box's wall overlap)
    # we take care of the first case here
    elif (sN <= 0.0 or D <= sN) and (0.0 < tN and tN < D):
        # uC = ui + (sC*u)
        # vC = vi + (tC*v)
        # wc = uC - vC
        # mp = vC + wc/2
        mp = ((ui + sC*u) + (vi + tC*v))/2
        # mp in unit box in first quadrant is same as ||R(pi/2)*(mp-0.5)||_1 < 1
        # but not sure if a rotation matrix is faster or slower than a bunch of
        # if statements
        if np.all((0.0 < mp) & (mp < 1.0)):
            return True
    # at this point, we know the point of closest approach is described by
    # either one tip and one interior point, and that we need to check if their
    # extents on the box's wall overlap, or the stick's two tips describe their
    # pointof closest approach, in which case we still only need to check if
    # their extents on the box overlap

    # then the sticks come closer and closer until they both hit their
    # respective ending walls. so we need only figure out what these
    # walls are and check if their outlines on the wall overlap or
    # not (accounting for fact that they might be fat enough for their
    # outlines to "spill" over onto a neighboring wall).
    # we number the walls 0-5, as in uniform_segment_from_unit_cube xz,
    # yz, xy planes and their opposites in the same order
    ui_walls = get_walls_with_outline(ui) #TODO periodize
    uf_walls = get_walls_with_outline(uf)
    vi_walls = get_walls_with_outline(vi)
    vi_walls = get_walls_with_outline(vf)
    isCollided = False
    for i in range(6):
        # if they cylinders' outlines on the box both have segment on a
        # particular wall, and these outlines overlap on the box's
        # wall, then the cylinders themselves overlap
        if (ui_walls[i] or uf_walls[i]) and (vi_walls[i] or vf_walls[i]) \
                and check_outline_overlap(ui, uf, vi, vf, wall=i):
            isCollided = True
            break
    return isCollided

def check_outline_overlap(ui, uf, vi, vf, wall):
    ellipse1 = get_ellipse_from_cylinder_wall(ui, uf, wall)
    ellipse2 = get_ellipse_from_cylinder_wall(vi, vf, wall)
    if ellipse1 and ellipse2:
        return check_ellipse_collisions(ellipse1, ellipse2)
    elif ellipse1:
        tube2 = get_tube2_from_cyclinder_wall(vi, vf, wall)
        return check_ellipse_tube2_collisions(ellipse=ellipse1, tube=tube2)
    elif ellipse2:
        tube1 = get_tube2_from_cyclinder_wall(ui, uf, wall)
        return check_ellipse_tube2_collisions(ellipse=ellipse2, tube=tube1)
    else:
        tube1 = get_tube2_from_cyclinder_wall(ui, uf, wall)
        tube2 = get_tube2_from_cyclinder_wall(vi, vf, wall)
        return check_tube2_tube2_collisions(tube1, tube2)

def get_ellipse_from_cylinder_wall(xi, xf, wall, r):
    w = ws[wall]
    w0 = w0s[wall]
    non_zero_index = non_zero_indices[wall]
    x = xf - xi
    Ip = line_hit_hyperplane(xi, xf, w, w0)
    # maybe they're not effectively parallell by numerical precision
    # (SMALL_NUM) test, but still effectively parallel in practice...say
    # collision happens very far from actual box
    if Ip is None \
    or norm3_squared(Ip - box_center) > MAX_BOXES_AWAY_FOR_ELLIPSE_COLLISION**2:
        return None
    # otherwise, we've determined the coordinates in the plane of interest
    h = Ip[non_zero_index[0]]
    k = Ip[non_zero_index[1]] # other index of Ip should be 1 or 0
    # angle from first coordinate to direction of long axis of ellipse
    theta = np.arctan(x[non_zero_index[1]]/x[non_zero_index[0]])
    # the minor axis is always the radius of the cylinder by definition
    b = r
    # the major axis is r*sec(phi), where phi is the (positive) angle between x and w
    # but we know phi = arccos(w.(x/||x||)), and sec(arccos(z)) = 1/z, so
    a = 1.0/np.dot(w, x/np.linalg.norm(x))
    return (a, b, theta, h, k)

def get_tube2_from_cyclinder_wall(xi, xf, wall, r):
    """Gets extents of 2D tube in which a (plane) wall, cylinder collision is
    effectively taking place for cylinders that are nearly perpendicular to the
    normal vector of the plane. This is returns as a start an end vector and a
    diameter in the plane."""
    w = ws[wall]
    w0 = w0s[wall]
    non_zero_index = non_zero_indices[wall]
    # the vector is nearly parallel to the plane, so we need only get the
    # distance from the plane to any point on the cylinder's center line to
    # determine its distance from the wall
    # also, w is already normalized in its definition
    dist = np.abs(d + np.dot(w, xi))
    # a plane parallel to a cylinder with radius r has a rectangular
    # intersection with the plane in the coordinates of the cylinder. the
    # height of this rectangle is just
    h = 2*np.sqrt(r*r - dist*dist)
    pi = xi*non_zero_index
    pf = xf*non_zero_index
    return (pi, pf, h)

def check_tube2_tube2_collisions(tube1, tube2):
    """Checks if the rectangular croos-sections of two cylinders parallel to
    the box's wall intersect in the the plane of the wall. The 2D tubular
    cross sections are passed in via two points along their center line and the
    the size of their extent as a distance from this center line. Thus, we need
    only check if the two corresponding infinite tube collide, then if they
    collide within the box's wall's boundaries. Two infinite tubes in 2D always
    collide unless their center lines are parallel. When they collide, their
    overlapping volume is a parallelogram. Either this parallelogram is
    contained by the box wall, or the two intersect. If they intersect, since
    they are both convex hulls of four points, it is enough to check that none
    of the vertices of the parallelgoram are contained in the box (and
    vice-versa) to conclude that they are not colliding."""
    xi, xf, h1 = tube1
    yi, yf, h2 = tube2
    x = xf - xi
    ui = x/np.linalg.norm(x)
    y = yf - yi
    vi = y/np.linalg.norm(y)
    # first check if tubes are parallell. if so, just return collision by
    # checking if the two lines as a whole are the correct distance from each
    # other
    if np.linalg.norm(ui - vi) < SMALL_NUM or np.linalg.norm(ui + vi) < SMALL_NUM:
        # distance between two parallel lines each defined by two points isn't
        # too hard to work out, just take the vector between any two points on
        # the line, then subtract off that vector's projection onto either line
        return np.linalg.norm((yi - xi) - np.dot(yi - xi, vi)) < h1 + h2
    # get the four vertices of the parallelogram that the two
    # tubes intersect in
    raise NotImplementedError('not done coding check_tube2_tube2_collisions')


def check_ellipse_collisions(ellipse1, ellipse2):
    a1,b1,theta1,h1,k1 = ellipse1
    a2,b2,theta2,h2,k2 = ellipse2
    raise NotImplementedError('not done coding check_ellipse_collisions')

@jit(nopython=True)
def uniform_point_from_unit_cube():
    return np.random.rand(3)

@jit(nopython=True)
def uniform_segment_from_unit_cube():
    """
    TODO: check if drawing a line through a random interior point with a random
    angle gives different distribution of stick density in box.
    Draws uniformly random segments from (unit) cube (in first quadrant)"""
    # sample a unit vector perpendicular to the xy-plane
    # using the "three normals" trick to generate a point on a sphere, then
    # ignoring the sign of the third component
    face_dx = np.random.randn(3)
    face_dx = face_dx/np.linalg.norm(face_dx)
    face_dx[2] = np.abs(face_dx[2]) # ignore sign
    # draw one of the six faces of the cube
    face = int(np.floor(np.random.rand(1)*6)[0]) # awkward casting to help jit
    # sample a point uniformly on that face
    face_x = np.random.rand(2)
    xinit = np.zeros(3)
    dxinit = np.zeros(3)
    # transform xinit and face_dx into face's coordinates
    if face == 0: # in xz-plane
        xinit[0] = face_x[0]
        xinit[2] = face_x[1]
        dxinit[1] = face_dx[2]
        dxinit[0] = face_dx[0]
        dxinit[2] = face_dx[1]
    elif face == 1: # in yz-plane
        xinit[1] = face_x[0]
        xinit[2] = face_x[1]
        dxinit[0] = face_dx[2]
        dxinit[1] = face_dx[0]
        dxinit[2] = face_dx[1]
    elif face == 2: # in xy-plane
        xinit[0:2] = face_x
        dxinit = face_dx
    elif face == 3: # opposite xz-plane
        xinit[0] = face_x[0]
        xinit[1] = 1.0
        xinit[2] = face_x[1]
        dxinit[1] = -face_dx[2]
        dxinit[0] = face_dx[0]
        dxinit[2] = face_dx[1]
    elif face == 4: # opposite yz-plane
        xinit[0] = 1.0
        xinit[1] = face_x[0]
        xinit[2] = face_x[1]
        dxinit[0] = -face_dx[2]
        dxinit[1] = face_dx[0]
        dxinit[2] = face_dx[1]
    elif face == 5: # opposite xy-plane
        xinit[0:2] = face_x
        xinit[2] = 1.0
        dxinit = face_dx
        dxinit[2] = -dxinit[2]
    #DEBUG
    # else:
    #     raise ValueError('internal error: face \\notin {0,1,2,3,4,5}')

    # now get collision point with other side of cube, by pretending our
    # segment is a path and asking at what time the path intersects the cube
    # if dx_i < 0, then x_i + dx_i*t exits the cube when x_i == 0
    # if dx_i > 0, then at x_i == 1
    tout = -xinit/dxinit # waste a couple cycles
    is_dx_pos = dxinit > 0
    if np.any(is_dx_pos):
        tout[is_dx_pos] = np.min((1.0 - xinit[is_dx_pos])/dxinit[is_dx_pos])
    xfinal = xinit + np.min(tout)*dxinit
    #DEBUG
    # if np.any(xfinal > 1.0):
    #     raise ValueError('internal error: xfinal should be inside unit cube')

    # get the segment initial and final locations
    return xinit, xfinal

@jit(nopython=True)
def uniform_segment_from_cube(a):
    """Draws uneformly random segments from cube (in first quadrant) defined by
    its side length."""
    xi,xf = uniform_segment_from_unit_cube()
    return a*xi, a*xf

@jit(nopython=True)
def volume_of_sphere(c, r, periodic=True):
    if periodic:
        return 4.0/3.0*np.pi*r*r*r
    else:
        raise NotImplementedError('volume of non-periodic sphere')
        return np.nan

# @jit(nopython=True)
def volume_of_cylinder(xi, xf, d, periodic=True):
    """Calculat the volume of a cylinder with center running from xi to xf
    inside a cube of side length a in the 1st quadrant. The formula used is the
    naive one, \pi r^2 l, but it is exactly correct if you're using periodic
    boundary conditions."""
    #TODO: ultra-low priority: implement exact formula for case when we are not
    # using periodic BCs
    # incorrect method: (see notes for full, correct method)
    # cases:
    #     1) hitting no edges - naive calculation works, same "in" as out
    #     2) hitting one edge - small modification
    #     3) hitting two+edges - ignored for now, hope these events are "rare"
    # case 1) d^2*pi*L is the correct formula. draw the "extra" volume you'd be
    # adding and subtracting from the "corners" that "jut out" of the box to
    # convince yourselves that they cancel out and you get back the naive
    # formula.
    # case 2) if you drew the first case out, you saw that the intersection of
    # the cylinder with the box edge is an ellipse, with major and minor axes
    # determined by the angle of intersection as (r|sec(\theta)|, r),
    # respectfully, where theta is angle from the normal to the wall to the
    # normal of the cylinder's natural cap. thus, you can consider case two as
    # the same as case 1, but subtracting the part ommited by the "collision"
    # with the wall not intersected by the center of the cylinder. since
    # finding the equation for this volume to be omitted seems hard in general,
    # since the ellipse's axes will in general not line up with the box edge,
    # we approximate by an intersection that is instead a circle with the
    # average of the major and minor axes as
    if periodic:
        xi, xf = get_true_xif_periodic(xi, xf)
    return np.pi*d*d*np.linalg.norm(xf - xi)

# @jit(nopython=True)
def get_true_xif_periodic(xi, xf):
    xi_face = get_face_given_point(xi)
    xf_face = get_face_given_point(xf)
    # if the faces are opposite of each other, then we are already fine
    if (xi_face, xf_face) not in face_pairs:
        return (xi, xf)
    # otherwise, we have to calculate which one is at a less extreme angle and
    # extend the other point beyond the box (recall ws are already normalized)
    mag_x = np.linalg.norm(xf - xi)
    thetai = np.arccos(np.dot(xi - xf, ws[xi_face])/mag_x)
    thetaf = np.arccos(np.dot(xf - xi, ws[xf_face])/mag_x)
    # if the angle is smaller, then the vector is more parallel to the normal
    # to the face, so we want to accept that face as the one we won't periodize
    # about
    if thetai > thetaf:
        # then we want new xi, found by continuing xf + t*(xi - xf) until it
        # hits the plane opposite to xf_face
        opposite_face = face_opposite[xf_face]
    else:
        opposite_face = face_opposite[xi_face]
    w = ws[opposite_face]
    w0 = w0s[opposite_face]
    # intersection point must always exist since we chose the plane of
    # intersection cleverly. order of xi,xf don't matter, the point will be
    # found regardless
    Ip = line_hit_hyperplane(xi, xf, w, w0)
    if thetai > thetaf:
        return (Ip, xf)
    else:
        return (xi, Ip)

def calculate_sphericity(lines):
    vectors = [line[1] - line[0] for line in lines]
    # vectors = [vector/np.linalg.norm(vector) for vector in vectors]
    outers = [np.outer(a, a) for a in vectors]
    inners = [np.inner(a, a) for a in vectors]
    S = np.zeros((3,3))
    Z = np.sum(inners)
    for s in outers:
        S += s
    S = S/Z
    eig = np.linalg.eig(S)
    evals = eig[0]
    evals.sort()
    return evals, (evals[0]+evals[1])*3/2

def plot_lines(lines):
    xinits = [line[0] for line in lines]
    xfinals = [line[1] for line in lines]
    return plot_segment(xinits, xfinals)

def plot_segment(xinits, xfinals):
    """Plot a list of uniformly generated segments."""
    # xfinals = np.stack(xfinals, axis=-1)
    # # 2-by-3-by-num_segments paths arrays
    # paths = np.stack(xinits, xfinals, axis=0)
    # easier than above, just plot one by one
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in range(len(xinits)):
        xinit = xinits[i]
        xfinal = xfinals[i]
        ax.plot([xinit[0], xfinal[0]], [xinit[1], xfinal[1]], [xinit[2], xfinal[2]])
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_zlim(0,1)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

def periodize_cylinders(xi, xf, d):
    """A generator for each of the copies of the cylinder that we have to check
    for collisions."""
    # of course the cylinder itself must be used, so we return that first in
    # case anybody only cares to do one, then they do the most important one
    yield xi, xf
    # now for each edge adjacent to the face that the tip of the line defining
    # the center of the cylinder is hitting, we can check if the cylinder
    # spills over onto the next face in that direction by just checking if the
    # points of closest approach from the cylinder to the edge are within r of
    # each other
    xip, xfp = get_true_xif_periodic(xi, xf)
    periodizers = periodizers_given_face[get_face_given_point(xip)]
    for periodizer in periodizers:
        yield periodizer+xip, periodizer+xfp
# old way to do periodization. unfinished. was going to "faster" by manually
# checking to make sure that we had to actually periodize by a face before we
# do. but logic for avoiding repeats,etc. seems hard and it costs a d_line_line
# to do the check anyway, which is the price we pay for just periodizing
# regardless and just making sure it's fine
    # periodization_faces = []
    # for edge_idx in adjacent_edges[xi_face]:
    #     edge = edges[edge_idx]
    #     if d_line_line(*edge, xi, xf):
    #         periodization_faces.append(face_edge_to_face[(xi_face, edge_idx)])

# @jit(nopython=True)
def periodize_point(p):
    periodizers = [-1.0, 0.0, 1.0]
    # points = []
    for offset_x in periodizers:
        for offset_y in periodizers:
            for offset_z in periodizers:
                # points.append(p + np.array([offset_x, offset_y, offset_z]))
                yield p + np.array([offset_x, offset_y, offset_z])
    # return points

# jitting this causes a memory leak
# @jit(nopython=True)
def detect_collision_spheres(c, r, c_new, r_new, periodic):
    if not periodic:
        raise NotImplementedError('collision detection for non-periodic spheres')
    for cp in periodize_point(c):
        for c_newp in periodize_point(c_new):
            if np.linalg.norm(cp - c_newp) < r + r_new:
                return True
    return False

def detect_collision_cylinders(xi, xf, xi_new, xf_new, d, periodic):
    """Are cylinders defined by start and end centerpoint sof xi,xf and
    xi_new,xf_new, each of diameter d, overlapping? use periodic boundary
    conditions in the box of side length a centered in the first quadrant if
    requested."""
    if not periodic:
        return d_segment_segment(xi, xf, xi_new, xf_new) < d
    # one for each of the 26 cubes adjacent to you that can affect you
    else:
        for xi_p, xf_p in periodize_cylinders(xi, xf, d):
            for xi_new_p, xf_new_p in periodize_cylinders(xi_new, xf_new, d):
                # not accurate if POCA is at tip of on or both cylinders
                if d_segment_segment(xi_p, xf_p, xi_new_p, xf_new_p) < d:
                    return True
        return False


# alternate collision detection, more accurate:
# first, if the point of closest approach of the cylinders' centers is interior
# to the cube, then the #TODO

def monte_carlo_with_removal(nsteps, mu, shape_params, shape, objects=None, phi=None,
                             periodic=True):
    """Perform nsteps MC steps. Half steps try to remove a uniformly random
    line with probability exp(-mu). Other half attempts (nsteps times)
    to insert a random line into a box of side length
    a. If the line overlaps with any other line (in the sense that they
    represent the centers of cylinders of radius d), reject the new line.

    Returns table of line insertion attempts and stats about the box at that
    time with a column specifying whether or not the move was accepted.

    TODO: calculate actual volume being added to box instead of just using the
    length*d^2*pi approximation, which only works for thin d."""
    if objects is None:
        objects = []
    num_objects = len(objects)
    if phi is None:
        phi = 0 # volume fraction of cylinders in box
    columns = columns_for_shape(shape)
    move_history = pd.DataFrame(index=np.arange(0, nsteps, 1), columns=columns)
    for i in range(nsteps):
        # two monte carlo moves available: insertion and removal
        if np.random.rand(1) < 0.5: # removal MC move
            # if no object exist, doens't make sense to try to remove them
            if num_objects == 0:
                move_history.loc[i] = row_for_shape(obj=None, shape=shape,
                        phi=phi, num_objects=num_objects, did_succeed=False,
                        is_removal=True)
                continue
            # choose a random object in the box
            remove_ind = int(np.random.rand(1)*num_objects)
            obj_rem = objects[remove_ind]
            # remove it with probability given by chemical potential "mu"
            # according to some parameter of the shape
            L = get_removal_amount(obj_rem, shape)
            did_succeed = np.random.rand(1) < np.exp(-mu*L)
            move_history.loc[i] = row_for_shape(obj_rem, shape, phi, num_objects,
                    did_succeed, is_removal=True)
            if did_succeed:
                objects.pop(remove_ind)
                num_objects -= 1
                phi -= volume_of_shape(obj_rem, shape_params, shape, periodic)
        else: # insertion MC move
            new_obj = new_object(shape_params, shape)
            did_succeed = True
            for existing_obj in objects:
                if detect_collision(existing_obj, new_obj, shape_params, shape, periodic):
                    did_succeed = False
                    break
            move_history.loc[i] = row_for_shape(new_obj, shape, phi, num_objects,
                    did_succeed, is_removal=False)
            if did_succeed:
                num_objects += 1
                phi += volume_of_shape(new_obj, shape_params, shape, periodic)
                objects += [new_obj]
    # the simulation is now done
    # hack in the correct types to the pandas dataframe holding the history fo
    # the moves we've done, since this can't be done at instantiation
    # (open bug in pandas as of 2017-02-24)
    numeric_columns = numeric_columns_for_shape(shape)
    move_history[numeric_columns] = move_history[numeric_columns].apply(pd.to_numeric)
    move_history['did_succeed'] = move_history.did_succeed.astype(bool)
    move_history['is_removal'] = move_history.is_removal.astype(bool)
    # return account of simulation as well as everything you need to continue
    # the simulation where it left off
    return objects, phi, move_history

#DEPRECATED!
# def monte_carlo(nsteps, d, a=1, lines=None, phi=None):
#     """Attempt (nsteps times) to put a random lines into a box of side length
#     a. If the line overlaps with any other line (in the sense that they
#     represent the centers of cylinders of radius d), reject the new line.

#     Returns table of line insertion attempts and stats about the box at that
#     time with a column specifying whether or not the move was accepted.

#     TODO: calculate actual volume being added to box instead of just using the
#     length*d^2*pi approximation, which only works for thin d."""
#     if lines is None:
#         lines = []
#     num_objects = len(lines)
#     if phi is None:
#         phi = 0 # volume fraction of cylinders in box
#     move_history = pd.DataFrame(index=np.arange(0, nsteps, 1),
#                                 columns=['xi1', 'xi2', 'xi3', 'xf1', 'xf2',
#                                          'xf3', 'phi', 'num_objects',
#                                          'did_succeed'])
#     for i in range(nsteps):
#         xi_new, xf_new = uniform_segment_from_cube(a)
#         did_succeed = True
#         for xi, xf in lines:
#             if d_segment_segment(xi, xf, xi_new, xf_new) < d:
#                 did_succeed = False
#                 break
#         move_history.loc[i] = [xi_new[0], xi_new[1], xi_new[2],
#                                xf_new[0], xf_new[1], xf_new[2],
#                                phi, num_objects, did_succeed]
#         if did_succeed:
#             num_objects += 1
#             phi += np.pi*d*d*np.linalg.norm(xf_new - xi_new)
#             lines += [(xi_new, xf_new)]
#     move_history = move_history.apply(pd.to_numeric)

#     return lines, phi, move_history

# this kind of refactoring is not nice because we can't save all the phi's in
# move_history this way
# def try_add_uniform_segment(lines, a=1):
#     xi_new, xf_new = uniform_segment_from_cube(a)
#     did_succeed = True
#     for xi, xf in lines:
#         if d_segment_segment(xi, xf, xi_new, xf_new) < d:
#             did_succeed = False
#             break
#     move_history.loc[i] = [xi_new[0], xi_new[1], xi_new[2],
#                             xf_new[0], xf_new[1], xf_new[2],
#                             phi, num_objects, did_succeed]
#     if did_succeed:
#         line = (xi_new, xf_new)
#         lines += [line]
#         return line
#     return None
def get_all_line_lengths(moves_df):
    return np.sqrt(
        np.power(moves_df['xf3'] - moves_df['xi3'], 2)
      + np.power(moves_df['xf2'] - moves_df['xi2'], 2)
      + np.power(moves_df['xf1'] - moves_df['xi1'], 2)
    )

# original function to run a bunch of sims to get the dependence of the
# acceptance probabilities on phi. we will instead be using the equilibrium
# values from an insertion/removal process to get the acceptance probabilities
# at equilibrium. if you want to come back and still do this, you should modify
# this function to work with the new monte_carlo functions
def get_accept_prob(outdir, nPhi, nL, steps_per_batch, mu,
                    shape_params, shape, num_sims=None):
    """Bin acceptance ratio dependence on phi into nP bins. Restart the
    MC sim every time that the rate that phi is increasing per nS steps is less
    than restart_thresh."""
    if steps_per_batch < 100:
        Warning('nS < 100 will result in bad averaging.')
    if nPhi < 2:
        raise ValueError('nPhi < 2 corresponds to no binning!')
    if nL < 2:
        raise ValueError('nL < 2 corresponds to no binning!')
    # get directory to save to
    outdir = wpath.get_unique_folder(outdir)
    # save index, currently not used
    num_restarts = 0
    save_ind = 0
    # file names to save to
    prob_file = "prob.txt"
    prob_file = os.path.join(outdir, prob_file)
    acc_file = "acc.txt"
    acc_file = os.path.join(outdir, acc_file)
    tot_file = "tot.txt"
    tot_file = os.path.join(outdir, tot_file)
    # initialize bin stuffs
    phiBins = np.linspace(0, 1, nPhi+1)
    phibins_file = os.path.join(outdir, "phibins.txt")
    np.savetxt(phibins_file, phiBins)
    LBins = np.linspace(0, 1*np.sqrt(3), nL+1)
    lbins_file = os.path.join(outdir, "lbins.txt")
    np.savetxt(lbins_file, LBins)
    accepted_moves = np.zeros((nPhi, nL))
    total_moves = np.zeros((nPhi, nL))
    # initialize inter-mc save state
    lines = None
    phi = None
    # initialize restart checker so that we always do a second simulation
    prev_phi_ave = -float('inf')
    while True:
        # perform MC
        lines, phi, new_moves = monte_carlo_with_removal(steps_per_batch,
                mu, shape_params, shape, lines, phi)
        new_moves['L'] = get_all_line_lengths(new_moves)
        # update ratios
        #TODO: check if np.digitize is faster
        for i in range(nPhi):
            for j in range(nL):
                tot_mask = ((phiBins[i] < new_moves['phi'])
                        & (new_moves['phi'] < phiBins[i+1])
                        & (LBins[j] < new_moves['L'])
                        & (new_moves['L'] < LBins[j+1])
                        & (~new_moves['is_removal']))
                acc_mask = tot_mask & (new_moves['did_succeed'])
                accepted_moves[i,j] += np.sum(acc_mask)
                total_moves[i,j] += np.sum(tot_mask)
        # check if should restart
        phi_ave = new_moves['phi'].mean()
        if phi_ave < prev_phi_ave:
            # reset mc state
            lines = None
            phi = None
            # reset to base case
            prev_phi_ave = -float('inf')
            num_restarts += 1
        else:
            # keep track fo previous phi for convergence criterion
            prev_phi_ave = phi_ave
        # write out current values for rate
        np.savetxt(prob_file, np.divide(accepted_moves, total_moves))
        np.savetxt(acc_file, accepted_moves)
        np.savetxt(tot_file, total_moves)
        save_ind += 1
        if num_sims is not None and num_restarts >= num_sims:
            break
    return accepted_moves, total_moves

def get_equilibrium(steps_per_batch, mu, shape_params, shape, lines=None,
                    phi=None, periodic=True):
    """For the values of d and mu provide, run until equilibrium. Here,
    equilibrium means run in batches of steps_per_batch MC steps, then if the
    average phi (packing fraction) in a batch drops below the average phi of
    the previous batch. One more batch is then run to calculate the equilibrium
    value.

    This works since the approach to equlibrium is an exponential (monotone)
    plus noise, so as long as steps_per_batch is chosen large enough that the
    reandom fluctuations during the rise to equilibrium get averaged out.

    For cylinder packign of radius 0.01 in a box of side length 1, a reasonable
    number for steps_per_batch is 1000. For all values tested that return
    reasonable answers, a reasonable number for steps_per_batch is 10000.
    """
    prev_phi_ave = -float('inf')
    lines = []; phi = 0;
    while True:
        lines, phi, moves = monte_carlo_with_removal(steps_per_batch, mu,
                shape_params, shape, lines, phi, periodic=periodic)
        moves = moves[~moves.is_removal]
        phi_ave = moves['phi'].mean()
        # if termination condition reached, run again so taht we decorrelate
        # from the artificial lower value of phi
        if phi_ave < prev_phi_ave:
            _, _, moves = monte_carlo_with_removal(steps_per_batch, mu,
                    shape_params, shape, lines, phi, periodic=periodic)
            moves = moves[~moves.is_removal]
            return moves['phi'].mean()
        prev_phi_ave = phi_ave


def equilibrium_mapper(p):
    return (p['s'], p['mu'], p['shape_param'],
            get_equilibrium(p['s'], p['mu'], p['shape_param'], p['shape'], periodic=p['p']))

def scan_equilibriums(steps_per_batch, mus, shape_params, shape, periodic=True, num_cores=8):
    """Use pscan to get a list of equilibrium phi levels as a function of the
    removal probability (i.e. mu a.k.a. chemical potential of the bath of
    "sticks" feeding our box filling process) and the diameter/box width ratio.

    steps_per_batch must be scalar, or same size as mus, they'll be
    jparams.
    """
    script_name = os.path.basename(__file__)
    #print(script_name + ': Running scan_equilibriums!')
    p = multiprocessing.Pool(num_cores)
    if type(steps_per_batch) is int:
        steps_per_batch = [steps_per_batch]
    if type(shape) is str:
        shape = [shape]
    if len(steps_per_batch) == 1:
        steps_per_batch = steps_per_batch*np.ones_like(mus, dtype=np.dtype(int))
    jparam = {'s': steps_per_batch, 'mu': mus}
    scan = pscan.Scan({'p': [periodic], 'shape_param': shape_params,
            'shape': shape})
    scan.add_jparam(jparam)
    print("mu\twidth_ratio\tequilibrium_phi\tnsteps")
    # for s,mu,sp,phi in map(equilibrium_mapper, scan.params()):
    for s,mu,sp,phi in p.imap_unordered(
            equilibrium_mapper, scan.params(), chunksize=10):
        print(str(mu) + '\t' + str(sp) + '\t' + str(phi) + '\t' + str(s))

def plot_equilibrium_scan(outfile):
    df = pd.read_table(outfile)
    hf = plt.figure()
    ha = hf.add_subplot(111, projection='3d')
    ha.scatter(df.width_ratio, df.equilibrium_phi, df.mu)
    ha.set_xlabel('$d$: Cylinder diameter/Box side length')
    ha.set_ylabel('$\phi$: Equilibrium packing density')
    ha.set_zlabel('$\mu$: Chemical potential/length')
    ha.set_title('Equilibrium Packing Densities from MC')
    return ha

# approx packing fract of cylinders is 0.73
def balls_approximation(phi, N, phi_full=0.73):
    phi_solvent = 1 - phi/phi_full
    return N*( phi_solvent*np.log(phi_solvent) - phi_solvent + 1 )

def fit_equilibrium_scan(outfile):
    """Try to get a fit for mu as a function of d and phi. So far, we tried
    "balls_approximation" above, which Quinn derived just based on
    inserting a comparable number of identifiable beads into a bin. but that
    function is concave up in planes of fixed d, and ours is concave down
    for the first long time before it slopes up rapidly (presumedly due to
    nematic effects, should check this)."""
    df = pd.read_table(outfile)
    ds = df.width_ratio.unique()
    ds.sort()
    for d in ds:
        dfd = df[df.width_ratio == d]
        phi = np.array(dfd.equilibrium_phi)
        mu = np.array(dfd.mu)
        popt, pcov = curve_fit(balls_approximation, phi, mu, p0=[d*d*d])
        import pdb; pdb.set_trace()



def plot_accept_prob(outdir):
    files = os.listdir(outdir)
    prob = None; acc = None; tot = None
    phiBins = None; Lbins = None
    for file in files:
        if file.startswith('prob'):
            prob = np.loadtxt(os.path.join(outdir, file))
        elif file.startswith('acc'):
            acc = np.loadtxt(os.path.join(outdir, file))
        elif file.startswith('tot'):
            tot = np.loadtxt(os.path.join(outdir, file))
        elif file.startswith('phibins'):
            phiBins = np.loadtxt(os.path.join(outdir, file))
        elif file.startswith('lbins'):
            Lbins = np.loadtxt(os.path.join(outdir, file))

    if prob is None or acc is None or tot is None \
            or phiBins is None or Lbins is None:
        raise FileNotFoundError('one of [phi/l]bins, prob, acc, tot not found in ' + outdir)
    # due to a bug, accepted rates not seem to be saved correctly
    acc = tot*prob
    # get bin centers for plotting
    dphiBin = phiBins[1] - phiBins[0]
    phiBinC = phiBins[:-1] + dphiBin/2
    dLbin = Lbins[1] - Lbins[0]
    LbinC = Lbins[:-1] + dphiBin/2
    phimesh = np.zeros((len(phiBinC), len(LbinC)))
    Lmesh = np.zeros((len(phiBinC), len(LbinC)))
    for i in range(len(phiBinC)):
        for j in range(len(LbinC)):
            phimesh[i,j] = phiBinC[i]
            Lmesh[i,j] = LbinC[j]
    X = phimesh; Y = Lmesh; Z = -np.log(prob)/Lmesh;
    Zave = -np.log(prob) # *Lmesh/Lmesh
    is_missing = np.isnan(Zave) | np.isinf(Zave)
    Zave[is_missing] = 0.0
    Lnorm = Lmesh;
    Lnorm[is_missing] = 0.0
    Zave = np.sum(Zave, axis=1)/np.sum(Lnorm, axis=1);
    hf = plt.figure()
    ha = hf.add_subplot(111, projection='3d')
    ha.scatter(X,Y,Z)
    ha.set_xlabel('$\phi$: Density Occupied')
    ha.set_ylabel('$\Delta{}L$: length of inserted stick')
    ha.set_zlabel('-ln[Prob(Accept)]/$\Delta{}L$')
    hf2 = plt.figure()
    ha2 = hf2.add_subplot(111)
    plt.scatter(X,Z,s=1,c=Y)
    plt.colorbar(label='$\Delta{}L')
    ha2.set_xlabel('$\phi$: Density Occupied')
    ha2.set_ylabel('-ln[Prob(Accept)]/$\Delta{}L$')
    plt.plot(X,Zave)
    return ha,ha2

def get_mu_of_d_phi(equilibrium_outfile, nphigrid=100):
    df = pd.read_table(equilibrium_outfile)
    dgrid = df.width_ratio.unique()
    dgrid.sort()
    phigrid = np.linspace(0, 1, nphigrid)
    mu = np.zeros((len(dgrid), nphigrid))
    mu[:,:] = np.nan
    zs = {}
    for i,d in enumerate(dgrid):
        lowess = sm.nonparametric.lowess
        dfd = df[df.width_ratio == d]
        x = dfd.equilibrium_phi
        y = dfd.mu
        z = lowess(y, x, frac=1.0/3.0)
        zs[d] = z
        xfit = z[:,0]
        yfit = z[:,1]
        # this should in principle be faster and generalize to arbitrary
        # dimensions, but DEFINITELY doesn't work as is. not sure what's
        # wrong. gives "almost" correct results
        # tree = cKDTree(np.expand_dims(xfit, axis=1))
        # dist, ind = tree.query(np.expand_dims(phigrid, axis=1), k=2)
        # d1, d2 = dist.T
        # x1, x2 = xfit[ind].T
        # v1, v2 = yfit[ind].T
        # mult = np.ones_like(v1)
        # mult[x1 > x2] = -1.0
        # mu[i,:] = mult*(d1)/(d1 + d2)*(v2 - v1) + v1
        # mu[i, phigrid < xfit[0]] = -float('inf')
        # mu[i, phigrid > xfit[-1]] = float('inf')
        for j,phi in enumerate(phigrid):
            inds = np.argwhere(xfit < phi)
            if inds.size == 0:
                mu[i,j] = -float('inf')
                continue
            ind_l = inds[-1]
            inds = np.argwhere(xfit > phi)
            if inds.size == 0:
                mu[i,j] = float('inf')
                continue
            ind_r = inds[0]
            x1 = xfit[ind_l]
            x2 = xfit[ind_r]
            y1 = yfit[ind_l]
            y2 = yfit[ind_r]
            mu[i,j] = y1 + (phi - x1)*(y2 - y1)/(x2 - x1)
    return dgrid, phigrid, mu

def plot_equilibrium_scan_with_interp(outfile, *args, **kwargs):
    ha = plot_equilibrium_scan(outfile)
    dgrid, phigrid, mu = get_mu_of_d_phi(outfile, *args, **kwargs)
    X, Y = np.meshgrid(dgrid, phigrid)
    ha.plot_surface(X, Y, mu.T, alpha=0.5)
    return ha

if __name__ == '__main__':

    # nsteps = np.linspace(1000, 5000, 101).astype(int)
    nsteps = 2000
    mus = np.flipud(np.linspace(0.01, 10, 101))
    ds = np.flipud(np.linspace(0.01, 0.3, 31))
    scan_equilibriums(nsteps, mus, ds, 'sphere', num_cores=32, periodic=True)

    # get_accept_prob('tmp/accept_probs', nPhi=100, nL=50,
    #                 steps_per_batch=1000, mu=5, d=0.01,
    #                 num_sims=1000)

    # lines, phi, move_history = monte_carlo_with_removal(10000, 1, 0.1, periodic=True)
    # print(move_history.phi.mean())



