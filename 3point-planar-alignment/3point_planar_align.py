#!/usr/bin/env python3
#------------------------------------------------------------------------------
# filename: 3point_planar_align.py
# author: Jon David
# date: Friday, March 9, 2018
# description:
#   Given an mx3 matrix M, where each row is a point in 3d space,
#   and three points p1, p2, and p3, returns M', a rotation of M
#   such that p1, p2, and p3 are co-planar with the XZ-plane,
#   where y=0; also returns T1, T2, T3, the three rotation
#   matrices used for this transformation.
#------------------------------------------------------------------------------

import numpy as np
import sys


USAGE_STR = """
ERROR: Not enough input arguments.
USAGE:
  python3 3point_planar_align.py <fp_name> <fm_name> <of_name>

  where,
    <fp_name> is the csv of three 3d points to align with XZ-plane, y=0;
    <fm_name> is the csv of all atom 3d points;
    <of_name> is the output csv of all the rotated 3d points;
"""


def read(fp_name, fm_name):
    """Reads fp_name, csv file of three 3d points; and fm_name, csv file of all atom 3d points. Returns (P,M), where P is a matrix of three 3d points, and M is a matrix of all atom 3d points."""
    P = np.loadtxt(fp_name, delimiter=',')
    M = np.loadtxt(fm_name, delimiter=',')
    return (P,M)

def rotate_all_points(T,X):
    """Rotates all points in X using transformation matrix T."""
    #return T.dot(X)
    #return X.dot(T)
    return T.dot(X.transpose())

def find_angle_between_vectors(v,w):
    """Finds the angle between vectors v,w."""
    v_magnitude = np.linalg.norm(v)
    w_magnitude = np.linalg.norm(w)

    print("--------------------")
    print("(DEBUG) v: {}, type(w): {}, shape: {}".format(v, type(v), v.shape))
    print("(DEBUG) w: {}, type(w): {}, shape: {}".format(w, type(w), w.shape))
    print("(DEBUG) v.dot(w): {}, type: {}".format(v.dot(w), type(v.dot(w))))
    print("(DEBUG) v_magnitude: {}, type: {}".format(v_magnitude, type(v_magnitude)))
    print("(DEBUG) w_magnitude: {}, type: {}".format(w_magnitude, type(w_magnitude)))
    
    theta = np.arccos( v.dot(w) / (v_magnitude*w_magnitude) )
    print("theta = {}".format( theta ))

    return theta


def align_vec1_to_XY_plane(v):
    """Align v to the XY-plane, z=0 by rotating about the y-axis"""
    v_proj_onto_XY_plane = np.array([v[0],v[1],0])  ## proj_XZ=Y(v), where z=0
    x_std_vec = np.array([1,0,0])
    #theta = find_angle_between_vectors(v,x_std_vec)
    theta = find_angle_between_vectors(v_proj_onto_XY_plane, x_std_vec)
    Ry = np.array([[np.cos(theta), 0, np.sin(theta)],
                   [0,1,0],
                   [-1*np.sin(theta), 0, np.cos(theta)]])
    return Ry


def align_vec1_to_x_axis(v):
    """Aligns v to the X-axis."""
    x_std_vec = np.array([1,0,0])
    theta = find_angle_between_vectors(v, x_std_vec)
    Rz = np.array([[np.cos(theta), -1*np.sin(theta),0],
                  [np.sin(theta), np.cos(theta),0],
                  [0,0,1]])
    return Rz


def align_vec2_to_YZ_plane(v):
    """Aligns v to the YZ-plane, y=0."""
    #z_std_vec = np.array([0,0,1])
    y_std_vec = np.array([0,1,0])
    v_proj_onto_YZ_plane = np.array([0,v[1],v[2]])  ## proj_YZ(v), where x=0
    #theta = find_angle_between_vectors(v_proj_onto_YZ_plane, z_std_vec)
    theta = find_angle_between_vectors(v_proj_onto_YZ_plane, y_std_vec)
    #theta = find_angle_between_vectors(v, z_std_vec)
    #theta = -1.0*theta
    Rx = np.array([[1,0,0],
                   [0,np.cos(theta),-1*np.sin(theta)],
                   [0,np.sin(theta),np.cos(theta)]])
    return Rx


##===== main ====
if len(sys.argv) < 4:
    print(USAGE_STR)
    exit()

fp_name = sys.argv[1]
fm_name = sys.argv[2]
of_name = sys.argv[3]
(P,M) = read(fp_name, fm_name)

p1 = P[0,:]
p2 = P[1,:]
p3 = P[2,:]

p1p2 = p2 - p1         #p1p2 and p1p3 lie on the same plane
p1p3 = p3 - p1
p1p2_orig = p1p2 - p1  #p1p2_orig, p1p3_orig have their base at origin
p1p3_orig = p1p3 - p1

print("p1p2: {}".format(p1p2))
print("p1p3: {}".format(p1p3))
print("p1p2_orig: {}".format(p1p2_orig))
print("p1p3_orig: {}".format(p1p3_orig))

##--- find series of transformations
T1 = align_vec1_to_XY_plane(p1p2_orig)       # rotate p1p2_orig about y-axis
T2 = align_vec1_to_x_axis( T1.dot(p1p2_orig) )      # now rotate about z-axis
T3 = align_vec2_to_YZ_plane( T2.dot(T1.dot(p1p3_orig)) ) # rotate p1p3_orig about x-axis
T = T3.dot(T2.dot(T1))

print("---- begin validation of transformations ----")
print("p1p2 as aligned to XY-plane: {}".format(T1.dot(p1p2_orig)))
print("p1p2 as aligned to X-plane: {}".format(T2.dot(T1.dot(p1p2_orig))))
print("p1p3 as aligned to XZ-plane: {}".format(T.dot(p1p3_orig)))
##--- rotate all points using those transformations
numPoints = max(M.shape)
X_orig = M - p1*np.ones(numPoints)  # move points towards origin
X_orig_rot = rotate_all_points(T,X_orig)          # rotate
X_rot = X_orig_rot + p1*np.ones(numPoints)   # move away from origin


print("\n---- begin debug points -----")
print("\n---- M ----")
print("{}".format(M))
print("\n---- M moved towards origin ----")
print("{}".format(X_orig))
print("\n---- M rotated ----")
print("{}".format(X_orig_rot))
print("\n---- M rotated, moved away from origin ----")
print("{}".format(X_rot))
np.savetxt(of_name, X_rot, delimiter=',')
