import pandas as pd
import numpy as np
from scipy.io import savemat, loadmat
from scipy.ndimage.interpolation import affine_transform
from skimage import transform as tf
from sklearn.neighbors import NearestNeighbors
import vtk
from vtk import *
import os
import sys
from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
from mpl_toolkits.mplot3d import Axes3D

from transform_utils import *
from pycpd import rigid_registration, affine_registration, deformable_registration
from functools import partial
import time

'''
Calculates cell assignments based on an affine transformation between each dataset  and the EM data based on nuclei location in 3D.
	- for Relay Neuron cells
	- calculates initial rotation matrix based on anchor points
	- calculates affine transformation using the coherent point drift algorithm implemented by pycpd
	- calculates and shows confusion matrix between datasets and the neurotransmitter assignment frequencies for each cell

Author: Angela Zhang
First Written: 10/4/2018
Last Edited: 1/4/2019
'''

plot_process = False
plot_transformed = False
compare_metrics = True

''' ----------- start algo ----------- '''

home_folder = '/Users/Unimusic/Desktop/Smith_Lab'
label_folder = os.path.join(home_folder, 'flipped rns 10.4.18')
gt_file = os.path.join(home_folder, 'REformatted Nuclear coordinates.xlsx')
transformed_all = []
num = 0
for subfolder in os.listdir(label_folder):
	if not os.path.isdir(subfolder):
		continue

	print subfolder

	currfolder = os.path.join(label_folder, subfolder)
	reg_file = os.path.join(currfolder, 'All_cells.xls')
	affine_name = os.path.join(currfolder, 'Affine_Matrix.mat')

	df1 = pd.read_excel(gt_file)
	df2 = pd.read_excel(reg_file)
	df3 = pd.merge(df1, df2)

	# load cells, gather into groups
	names = np.unique(df2.values[:,0])
	cellinds = df2.index[df2['ID'].str.lower().str.contains('vgat|vacht')]
	cellinds2 = df2.index[df2['ID'].str.lower() == 'ddn']
	cells1 = df2.iloc[cellinds,:]
	cells2 = df2.iloc[cellinds2,:]

	df = df2

	cell_coords = df.values[:, 1:4].astype(float)
	cell_type = df.values[:, 0]
	cell_ID = df.values[:, 8]
	cell_ID[np.where(cell_type == 'ddN')] = 'ddN'
	neurotransmitters = df.values[:, 0]

	inds1 = df1.index[df1['ID'].str.lower().str.contains('rn')]
	inds2 = df1.index[df1['ID'].str.lower().str.startswith('ddn')]
	inds_1 = cellinds
	inds_2 = cellinds2

	inds = inds1.append(inds2)

	target_cells = df1.iloc[inds, :]
	target_coords = target_cells.values[:, 1:4].astype(float)
	tc_ID = target_cells.values[:, 0]

	target_a1 = df1.iloc[inds1, :]
	ta1_coords = target_a1.values[:, 1:4].astype(float)
	ta1_ID = target_a1.values[:, 0]

	target_a2 = df1.iloc[inds2, :]
	ta2_coords = target_a2.values[:, 1:4].astype(float)
	ta2_ID = target_a2.values[:, 0]

	curr_a1 = df2.iloc[inds_1]
	ca1_coords = curr_a1.values[:, 1:4].astype(float)
	ca1_ID = curr_a1.values[:, 8]

	curr_a2 = df2.iloc[inds_2]
	ca2_coords = curr_a2.values[:, 1:4].astype(float)
	ca2_ID = curr_a2.values[:, 8]

	target_anchor = np.array([np.mean(ta1_coords, axis=0), np.mean(ta2_coords, axis=0)])
	source_anchor = np.array([np.mean(ca1_coords, axis=0), np.mean(ca2_coords, axis=0)])

	num += 1

	# calculate intial rotation matrix based on anchor cells
	source_anchor -= np.mean(source_anchor, axis=0)
	target_anchor -= np.mean(target_anchor, axis=0)
	rot_mat = calc_rot(source_anchor, target_anchor)
	print rot_mat
	cell_coords -= np.mean(cell_coords, axis=0)
	transformed_init = np.matmul(rot_mat, cell_coords.transpose()).transpose()

	# calculate affine transform using coherent point drift
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	callback = partial(visualize, ax=ax)
	reg = affine_registration(**{ 'X': target_coords, 'Y': transformed_init, 'max_iterations': 10000})
	reg.register(callback)
	if plot_process:
		plt.show()
	transformed_coords = reg.TY

	new_transformed_coords, new_target_coords, new_tcID, new_tcID = get_relay_only(transformed_coords,
															target_coords,
															cell_ID,
															tc_ID)

	# get nearest neighbor assignments only (iterations=1), could also get affine transform using icp
	icp1 = get_icp(new_transformed_coords, new_target_coords, 'Similarity')
	transformed_coords = get_transformed(icp1, transformed_coords)

	if plot_process:
		create_scatterplot(cell_coords, target_coords)
		create_scatterplot(transformed_init, target_coords)
		create_scatterplot(transformed_coords, target_coords)

	cell_map, ind_map = create_coord_map(transformed_coords, target_coords, cell_ID, tc_ID)

	if plot_transformed:
		color=cm.rainbow(np.linspace(0,1,len(ind_map)))
		fig = plt.figure()
		ax1 = fig.add_subplot(211, projection='3d')
		ax2 = fig.add_subplot(212, projection='3d')
		color_ind = 0
		for pair in ind_map:
			target_coord = pair[0]
			transformed_coord = pair[1][0]
			ax1.scatter(transformed_coord[0],
					transformed_coord[1],
					transformed_coord[2], c=color[color_ind], marker='o')
			ax2.scatter(target_coord[0],
						target_coord[1],
						target_coord[2], c=color[color_ind], marker='^')
			color_ind += 1
		ax1.set_xlabel('X')
		ax1.set_ylabel('Y')
		ax1.set_zlabel('Z')
		ax2.set_xlabel('X')
		ax2.set_ylabel('Y')
		ax2.set_zlabel('Z')

		plt.show()

	neuro_map = {}
	neuro_record = {}
	for j in range(len(neurotransmitters)):
		transmitter = neurotransmitters[j]
		if isinstance(transmitter, float) or cell_ID[j] not in cell_map.keys():
			continue
		neuro_map[cell_map[cell_ID[j]][0]] = transmitter
		neuro_record[cell_ID[j]] = transmitter

	transformed_all.append((subfolder, transformed_coords, cell_ID, cell_map, neuro_map, neuro_record))

	final_df = pd.DataFrame.from_dict(cell_map, orient='index', columns = ['Assigned Cell', 'Distance'])
	final_df = final_df.reset_index()
	final_df = final_df.rename(index=str, columns={'index':'ID'})
	final_df = final_df[['ID', 'Assigned Cell', 'Distance']]
	final_df.to_csv(os.path.join(currfolder, 'registered_cells.csv'), index=False)

if compare_metrics:
	# calculate consistency
	inds = df1.index[df1['ID'].str.lower().str.contains('rn-')]
	all_keys = df1.iloc[inds]['ID'].values
	all_keys.sort()
	all_keys = all_keys.tolist()
	print all_keys
	confusion_matrix = np.zeros((len(all_keys),len(all_keys)))
	error_bar = np.zeros((len(all_keys),1))
	neuro_matrix = np.zeros((3,len(all_keys)))
	neuro_consistency = np.zeros((2,2))

	for i in range(len(transformed_all)):
		_,_,_,_,neuro_map,_ = transformed_all[i]
		for cell, transmitter in neuro_map.items():
			if transmitter == 'vgat':
				neuro_matrix[0, all_keys.index(cell)] += 1
			elif transmitter == 'vacht':
				neuro_matrix[1, all_keys.index(cell)] += 1
			else:
				neuro_matrix[2, all_keys.index(cell)] += 1

	for i in range(len(transformed_all)):
		for j in range(i+1, len(transformed_all)):
			subfolder1, rn1, id1, map1, _, nr1 = transformed_all[i]
			subfolder2, rn2, id2, map2, _, nr2 = transformed_all[j]

			print subfolder1
			print subfolder2

			inv_map1 = {v[0]: k for k, v in map1.iteritems()}
			inv_map2 = {v[0]: k for k, v in map2.iteritems()}

			curr_ind = 0
			for curr_id in id1:
				if curr_id not in map1.keys():
					rn1 = np.delete(rn1, curr_ind, axis=0)
					id1 = np.delete(id1, curr_ind, axis=0)
					curr_ind -= 1
				curr_ind += 1

			curr_ind = 0
			for curr_id in id2:
				if curr_id not in map2.keys():
					rn2 = np.delete(rn2, curr_ind, axis=0)
					id2 = np.delete(id2, curr_ind, axis=0)
					curr_ind -= 1
				curr_ind += 1

			fig = plt.figure()
			ax = fig.add_subplot(111, projection='3d')
			callback = partial(visualize, ax=ax)
			reg = affine_registration(**{ 'X': rn1, 'Y': rn2, 'max_iterations': 10000})
			reg.register(callback)
			if plot_process:
				plt.show()
			coords3 = reg.TY

			icp3 = get_icp(coords3, rn1, mode='Affine')
			coords3 = get_transformed(icp3, rn2)
			map3 = create_coord_map(coords3, rn1, id2, id1)[0]

			for key in map3.keys():
				t_key = map3[key][0]
				if t_key not in map1.keys() or key not in map2.keys():
					continue
				tkey_cellID = map1[t_key][0]
				key_cellID = map2[key][0]
				if 'rn' not in tkey_cellID.lower() or 'rn' not in key_cellID.lower():
					continue
				x_axis = all_keys.index(tkey_cellID)
				y_axis = all_keys.index(key_cellID)
				if 'rn' in tkey_cellID.lower() and 'rn' in key_cellID.lower():
					error_bar[np.abs(x_axis-y_axis)] += 1
				confusion_matrix[x_axis, y_axis] += 1
				confusion_matrix[y_axis, x_axis] += 1
				if key not in nr2.keys() or t_key not in nr1.keys():
					continue
				if nr2[key] == 'vgat':
					if nr1[t_key] == 'vgat':
						neuro_consistency[0,0] += 1
					elif nr1[t_key] == 'vacht':
						neuro_consistency[0,1] += 1
				elif nr2[key] == 'vacht':
					if nr1[t_key] == 'vacht':
						neuro_consistency[1,1] += 1
					elif nr1[t_key] == 'vgat':
						neuro_consistency[1,0] += 1

	plt.matshow(confusion_matrix, cmap='hot')
	plt.xticks(range(error_bar.size),all_keys, rotation='vertical')
	plt.yticks(range(error_bar.size),all_keys)
	plt.colorbar()
	plt.show()

	print neuro_consistency

	plt.matshow(neuro_consistency, cmap='YlOrRd_r')
	plt.colorbar()
	plt.show()

	plt.matshow(neuro_matrix[0:2,:], cmap='YlOrRd_r')
	plt.xticks(range(error_bar.size),all_keys, rotation='vertical')
	plt.yticks(range(3), ['vgat', 'vacht'])
	plt.colorbar()
	plt.show()
