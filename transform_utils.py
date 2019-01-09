import numpy as np
from sklearn.neighbors import NearestNeighbors
import vtk
from vtk import *
import os
import sys
from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
from mpl_toolkits.mplot3d import Axes3D

def get_icp(spoints, tpoints, mode='Similarity'):
	#iterative closest point, using VTK.
	sourcePoints = vtk.vtkPoints()
	sourceVertices = vtk.vtkCellArray()

	for i in range(0, spoints.shape[0]):
		row = spoints[i,:]
		id = sourcePoints.InsertNextPoint(row[:])
		sourceVertices.InsertNextCell(1)
		sourceVertices.InsertCellPoint(id)

	source = vtk.vtkPolyData()
	source.SetPoints(sourcePoints)
	source.SetVerts(sourceVertices)
	if vtk.VTK_MAJOR_VERSION <= 5:
		source.Update()

	targetPoints = vtk.vtkPoints()
	targetVertices = vtk.vtkCellArray()

	for i in range(0, tpoints.shape[0]):
		row = tpoints[i,:]
		id = targetPoints.InsertNextPoint(row[:])
		targetVertices.InsertNextCell(1)
		targetVertices.InsertCellPoint(id)
	target = vtk.vtkPolyData()
	target.SetPoints(targetPoints)
	target.SetVerts(targetVertices)
	if vtk.VTK_MAJOR_VERSION <= 5:
		target.Update()
	icp = vtk.vtkIterativeClosestPointTransform()
	icp.SetSource(source)
	icp.SetTarget(target)
	if mode == 'Affine':
		icp.GetLandmarkTransform().SetModeToAffine()
	elif mode == 'Rigid':
		icp.GetLandmarkTransform().SetModeToRigidBody()
	icp.SetMaximumNumberOfIterations(1)	#iterations=1 for finding nearest neighbors only
	icp.StartByMatchingCentroidsOn()
	icp.SetMaximumMeanDistance(0.001)
	icp.Modified()
	icp.Update()
	return icp

def get_transformed(icp, spoints):
	#returns transformed points given the transformation matrix (icp) and the original source points (spoints).
	sourcePoints = vtk.vtkPoints()
	sourceVertices = vtk.vtkCellArray()

	for i in range(0, spoints.shape[0]):
		row = spoints[i,:]
		id = sourcePoints.InsertNextPoint(row[:])
		sourceVertices.InsertNextCell(1)
		sourceVertices.InsertCellPoint(id)

	source = vtk.vtkPolyData()
	source.SetPoints(sourcePoints)
	source.SetVerts(sourceVertices)
	if vtk.VTK_MAJOR_VERSION <= 5:
		source.Update()

	icpTransformFilter = vtk.vtkTransformPolyDataFilter()
	if vtk.VTK_MAJOR_VERSION <= 5:
		icpTransformFilter.SetInput(source)
	else:
		icpTransformFilter.SetInputData(source)

	icpTransformFilter.SetTransform(icp)
	icpTransformFilter.Update()

	transformedSource = icpTransformFilter.GetOutput()
	pointCount = spoints.shape[0]

	transformed_coords = []

	for index in range(pointCount):
		point = [0,0,0]
		transformedSource.GetPoint(index, point)
		if np.isnan(point[0]):
			ipdb.set_trace()
		transformed_coords.append(point)
	transformed_coords = np.array(transformed_coords)
	return transformed_coords


def create_coord_map(transformed_coords, tpoints, sID, tID):
	'''
	returns a dictionary mapping the cell number in the source points to a tuple containing
		(cell ID in the target points, distance between target point and transformed source point),
		as well as  a list of the mapping (to keep things in order).
	'''
	if isinstance(tID[0], str) or isinstance(tID[0], unicode):
		transformed_coords, tpoints, sID, tID = get_pr_only(transformed_coords, tpoints, sID, tID)
	coord_map = {}
	nn = NearestNeighbors(n_neighbors=min(tpoints.shape[0],transformed_coords.shape[0])).fit(transformed_coords)
	distances, indices = nn.kneighbors(tpoints)
	coord_map = get_closest_points(distances, indices, coord_map)
	new_coord_map = {}
	new_ind_lst = []
	for key in coord_map.keys():
		num, dist = coord_map[key]
		new_entry = (tpoints[num,:], dist)
		new_key = transformed_coords[key, :]
		new_ind_lst.append((new_key, new_entry))
		new_entry = (tID[num], dist)
		new_key = sID[key]
		new_coord_map[new_key] = new_entry
	return new_coord_map, new_ind_lst


def get_pr_only(coords, tpoints, sID, tID):
	#returns only photoreceptor cell data given data of all cells.
	i = 0
	while i < tID.size:
		if not tID[i].lower().startswith('pr'):
			tpoints = np.delete(tpoints, i, axis=0)
			tID = np.delete(tID, i)
		else:
			i += 1
	i = 0
	while i < sID.size:
		if not isinstance(sID[i], int):
			coords = np.delete(coords, i, axis=0)
			sID = np.delete(sID, i)
		else:
			i += 1
	return coords, tpoints, sID, tID


def get_relay_only(coords, tpoints, sID, tID):
	#returns only relay neuron data given data of all cells.
	i = 0
	while i < tID.size:
		if 'rn' not in tID[i].lower():
			tpoints = np.delete(tpoints, i, axis=0)
			tID = np.delete(tID, i)
		else:
			i += 1
	i = 0
	while i < sID.size:
		if not isinstance(sID[i], int):
			coords = np.delete(coords, i, axis=0)
			sID = np.delete(sID, i)
		else:
			i += 1
	return coords, tpoints, sID, tID


def get_closest_points(distances, indices, coord_map, method='simple'):
	#return mapping of indices to cell ID using varioous forms of nearest neighbors.
	if method=='minimax':
		order = np.argsort(distances[:,0])
		order = order[::-1]
		distances = distances[order]
		indices = indices[order]
	for i in range(indices.shape[0]):
		dist = distances[i,:]
		min_dist = np.min(dist)
		dist_ind = np.where(dist == min_dist)[0][0]
		nearest_ind = indices[i, dist_ind]
		if method=='minimin':
			if nearest_ind in coord_map.keys():
				prev_ind, prev_dist = coord_map[nearest_ind]
				if prev_ind == i:
					continue
				if min_dist < prev_dist:
					coord_map[nearest_ind] = (i, min_dist)
					distances[prev_ind, np.where(distances[prev_ind,:]==prev_dist)[0][0]] = np.inf
					return get_closest_points(distances, indices, coord_map)
				else:
					while (prev_dist <= min_dist) and (nearest_ind in coord_map.keys()):
						distances[i, dist_ind] = np.inf
						dist = distances[i,:]
						min_dist = np.min(distances[i,:])
						if np.isinf(min_dist):
							return coord_map
						dist_ind = np.where(dist == min_dist)
						nearest_ind = indices[i, dist_ind]
					if i >= indices.shape[1]:
						return coord_map
					else:
						return get_closest_points(distances, indices, coord_map)
			else:
				coord_map[nearest_ind] = (i, min_dist)
		elif method=='simple' or method=='minimax':
			if np.isinf(min_dist):
				continue
			elif min_dist > 15:
				continue
			coord_map[nearest_ind] = (i, min_dist)
			distances[np.where(indices == nearest_ind)] = np.inf
	return coord_map


def get_length(vector):
	# returns the length (2-norm) of a vector.
	return np.sqrt(np.sum(np.square(vector)))


def calc_rot(points_a, points_b):
	# returns rotation matrix, assumes points_a and points_b are 2x3 matrices.
	v1 = points_a[1,:] - points_a[0,:]
	v2 = points_b[1,:] - points_b[0,:]
	v = np.cross(v1, v2)
	v = v/(get_length(v))
	theta = np.arccos(np.dot(v1, v2) / (get_length(v1)*get_length(v2)))
	rotmat = create_rot(v, theta)
	return rotmat


def create_rot(v, theta):
	# returns theta in radians.
	cos = np.cos(theta)
	sin = np.sin(theta)
	oneminus = 1-cos

	rot = np.array([[cos+v[0]**2*oneminus,
						v[0]*v[1]*oneminus - v[2]*sin,
						v[0]*v[2]*oneminus + v[2]*sin],
					[v[1]*v[0]*oneminus + v[2]*sin,
						cos + v[1]**2*oneminus,
						v[1]*v[2]*oneminus - v[0]*sin],
					[v[2]*v[0]*oneminus - v[1]*sin,
						v[2]*v[1]*oneminus + v[0]*sin,
						cos + v[2]**2*oneminus]])
	return rot

def create_scatterplot(points1, points2):
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.scatter(points1[:,0],
			points1[:,1],
			points1[:,2], c='r')
	ax.scatter(points2[:,0],
			points2[:,1],
			points2[:,2], c='b')
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
	plt.show()


def visualize(iteration, error, X, Y, ax):
    #visualize coherent point drift algorithm.
	plt.cla()
	ax.scatter(X[:,0],	X[:,1], X[:,2], color='red', label='Target')
	ax.scatter(Y[:,0],	Y[:,1], Y[:,2], color='blue', label='Source')
	ax.text2D(0.87, 0.92, 'Iteration: {:d}\nError: {:06.4f}'.format(iteration, error), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize='x-large')
	ax.legend(loc='upper left', fontsize='x-large')
	plt.draw()
	plt.pause(0.001)
