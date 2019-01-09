import pandas as pd
import numpy as np
from scipy.io import savemat, loadmat
from scipy.ndimage.interpolation import affine_transform
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
Calculates transformation between each dataset in a folder and the EM data based on nuclei location in 3D.
	- for photoreceptor cells only
	- calculates initial rotation matrix based on anchor points (pr body, ant/ddn cells
	- calculates affine transformation using coherent point drift algorithm implemented by pycpd
	- aligns points based on nearest neighbors post-transformation
	- plots confusion matrix and neurotransmitter type assignment frequency

Author: Angela Zhang
First Written: 10/3/2018
Last Edited; 1/4/2019
'''

plot_process = False
plot_transformed = False
compare_metrics = True

''' ----------- start algo ----------- '''

home_folder = '/Users/Unimusic/Desktop/Smith_Lab'
label_folder = os.path.join(home_folder, 'flipped prs 10.3.18')
gt_file = os.path.join(home_folder, 'REformatted Nuclear coordinates.xlsx')
transformed_prs = []
for subfolder in os.listdir(label_folder):
	print subfolder
	use_ant = False
	if not os.path.isdir(subfolder):
		continue
	currfolder = os.path.join(label_folder, subfolder)
	reg_file = os.path.join(currfolder, 'All_cells.xlsx')
	if not os.path.exists(reg_file):
		reg_file=reg_file[:-1]
	affine_name = os.path.join(currfolder, 'Affine_Matrix.mat')

	df1 = pd.read_excel(gt_file)
	df2 = pd.read_excel(reg_file)
	if any(df2['ID'].str.lower().str.contains('ant')):
		use_ant = True

	use_pr2 = False

	# load pr cells, assign groups
	pr_names = np.unique(df2.values[:,0])
	if (not 'ddN' in pr_names) and (not use_ant):
		use_pr2 = True
		indspr1 = df2.index[df2['ID'].str.lower() == 'pr']
		indspr2 = df2.index[df2['ID'].str.lower().str.startswith('pr (ii)')]
		pr1 = df2.iloc[indspr1,:]
		pr2 = df2.iloc[indspr2,:]
		pr_df = [df2]
	else:
		pr_df = [df2]

	df = pr_df[0]

	pr_coords = df.values[:, 1:4].astype(float)
	neurotransmitters = df.values[:, 9]
	pr_type = df.values[:, 0]
	pr_ID = df.values[:, 8]
	for prind in range(pr_type.size):
		if not pr_type[prind].lower().startswith('pr'):
			pr_ID[prind] = pr_type[prind].lower()

	if use_ant:
		inds1 = df1.index[df1['ID'].str.lower().str.startswith('pr-')]
		inds2 = df1.index[df1['ID'].str.lower().str.startswith('ant')]
		inds_1 = df2.index[df2['ID'].str.lower().str.startswith('pr')]
		inds_2 = df2.index[df2['ID'].str.lower().str.startswith('ant')]
		inds = inds1.append(inds2)

	elif not use_pr2:
		inds1 = df1.index[df1['ID'].str.lower().str.startswith('pr-')]
		inds2 = df1.index[df1['ID'].str.lower().str.startswith('ddn')]
		inds_1 = df2.index[df2['ID'].str.lower().str.startswith('pr')]
		inds_2 = df2.index[df2['ID'].str.lower().str.startswith('ddn')]
		inds = inds1.append(inds2)

	else:
		inds1 = df1.index[df1['ID'].str.lower().str.startswith('pr-') & df1['ID'].str.lower().str.split('-').str[1].str.isnumeric()]
		inds2 = df1.index[df1['ID'].str.lower().str.startswith('pr-') & df1['ID'].str.lower().str.split('-').str[1].str.isalpha()]
		inds_1 = df2.index[df2['ID'].str.lower()=='pr']
		inds = inds1.append(inds2)

	target_pr = df1.iloc[inds, :]
	tp_coords = target_pr.values[:, 1:4].astype(float)
	tp_ID = target_pr.values[:, 0]

	if not use_pr2:

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

		source_anchor -= np.mean(source_anchor, axis=0)
		target_anchor -= np.mean(target_anchor, axis=0)
		rot_mat = calc_rot(source_anchor, target_anchor)
		print rot_mat
		pr_coords -= np.mean(pr_coords, axis=0)
		transformed_init = np.matmul(rot_mat, pr_coords.transpose()).transpose()
	else:
		transformed_init = pr_coords

	# calculate affine transform using coherent point drift
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	callback = partial(visualize, ax=ax)
	reg = affine_registration(**{ 'X': tp_coords, 'Y': transformed_init, 'max_iterations': 10000})
	reg.register(callback)
	if plot_process:
		plt.show()
	transformed_coords = reg.TY

	new_transformed_coords, new_tp_coords, new_pr_ID, new_tp_ID = get_pr_only(transformed_coords,
															tp_coords,
															pr_ID,
															tp_ID)
	# This is just to get nearest neighbors, but could be used to calculate the similarity transformation using icp
	icp1 = get_icp(new_transformed_coords, new_tp_coords, 'Similarity')
	transformed_coords = get_transformed(icp1, transformed_coords)

	pr_map, ind_map = create_coord_map(transformed_coords, tp_coords, pr_ID, tp_ID)

	if plot_transformed:
		color=cm.rainbow(np.linspace(0,1,len(ind_map)))
		fig = plt.figure()
		ax1 = fig.add_subplot(1,2,1, projection='3d')
		ax2 = fig.add_subplot(1,2,2, projection='3d')
		color_ind = 0
		for pair in ind_map:
			tp_coord = pair[0]
			transformed_coord = pair[1][0]
			ax1.scatter(transformed_coord[0],
					transformed_coord[1],
					transformed_coord[2], c=color[color_ind], marker='o')
			ax2.scatter(tp_coord[0],
						tp_coord[1],
						tp_coord[2], c=color[color_ind], marker='^')
			color_ind += 1
		ax1.set_xlabel('X')
		ax1.set_ylabel('Y')
		ax1.set_zlabel('Z')
		ax1.set_xlabel('X')
		ax1.set_ylabel('Y')
		ax1.set_zlabel('Z')

		plt.show()


	neuro_map = {}
	neuro_record = {}
	for j in range(len(neurotransmitters)):
		transmitter = neurotransmitters[j]
		if isinstance(transmitter, float) or pr_ID[j] not in pr_map.keys():
			continue
		neuro_map[pr_map[pr_ID[j]][0]] = transmitter
		neuro_record[pr_ID[j]] = transmitter

	transformed_prs.append((subfolder, transformed_coords, pr_ID, pr_map, neuro_map, neuro_record))

	final_pr_df = pd.DataFrame.from_dict(pr_map, orient='index', columns = ['Assigned Cell', 'Distance'])
	final_pr_df = final_pr_df.reset_index()
	final_pr_df = final_pr_df.rename(index=str, columns={'index':'ID'})
	final_pr_df = final_pr_df[['ID', 'Assigned Cell', 'Distance']]
	final_pr_df.to_csv(os.path.join(currfolder, 'registered_pr.csv'), index=False)

if compare_metrics:
	# calculate consistency
	inds = df1.index[df1['ID'].str.lower().str.startswith('pr-')]
	all_keys = df1.iloc[inds]['ID'].values
	all_keys.sort()
	all_keys = all_keys.tolist()
	print all_keys
	confusion_matrix = np.zeros((len(all_keys),len(all_keys)))
	error_bar = np.zeros((len(all_keys),1))
	neuro_matrix = np.zeros((2,len(all_keys)))
	neuro_consistency = np.zeros((2,2))
	cell_freq = {}
	for i in range(len(transformed_prs)):
		_,_,_,_,neuro_map,_ = transformed_prs[i]
		for cell, transmitter in neuro_map.items():
			if transmitter == 'vgat':
				neuro_matrix[0, all_keys.index(cell)] += 1
			elif transmitter == 'vglut':
				neuro_matrix[1, all_keys.index(cell)] += 1
			else:
				neuro_matrix[0, all_keys.index(cell)] += 1
				neuro_matrix[1, all_keys.index(cell)] += 1
			if cell in cell_freq.keys():
				cell_freq[cell] += 1
			else:
				cell_freq[cell] = 1
	print cell_freq
	neuro_matrix[0,:] /= float(7)
	neuro_matrix[1,:] /= float(4)

	for i in range(len(transformed_prs)):
		for j in range(i+1, len(transformed_prs)):
			subfolder1, pr1, id1, map1, _, nr1 = transformed_prs[i]
			subfolder2, pr2, id2, map2, _, nr2 = transformed_prs[j]

			print subfolder1
			print subfolder2

			inv_map1 = {v[0]: k for k, v in map1.iteritems()}
			inv_map2 = {v[0]: k for k, v in map2.iteritems()}

			curr_ind = 0
			for curr_id in id1:
				if curr_id not in map1.keys():
					pr1 = np.delete(pr1, curr_ind, axis=0)
					id1 = np.delete(id1, curr_ind, axis=0)
					curr_ind -= 1
				curr_ind += 1

			curr_ind = 0
			for curr_id in id2:
				if curr_id not in map2.keys():
					pr2 = np.delete(pr2, curr_ind, axis=0)
					id2 = np.delete(id2, curr_ind, axis=0)
					curr_ind -= 1
				curr_ind += 1

			icp3 = get_icp(pr2, pr1, mode='Affine')
			coords3 = get_transformed(icp3, pr2)
			map3 = create_coord_map(coords3, pr1, id2, id1)[0]

			for key in map3.keys():
				t_key = map3[key][0]
				if t_key not in map1.keys() or key not in map2.keys():
					continue
				tkey_cellID = map1[t_key][0]
				key_cellID = map2[key][0]
				if not tkey_cellID.lower().startswith('pr') or not key_cellID.lower().startswith('pr'):
					continue
				x_axis = all_keys.index(tkey_cellID)
				y_axis = all_keys.index(key_cellID)
				if tkey_cellID.lower().startswith('pr') and key_cellID.lower().startswith('pr'):
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
				elif nr2[key] == 'vglut':
					if nr1[t_key] == 'vglut':
						neuro_consistency[1,1] += 1
					elif nr1[t_key] == 'vgat':
						neuro_consistency[1,0] += 1

	plt.matshow(confusion_matrix, cmap='hot')
	plt.xticks(range(error_bar.size),all_keys, rotation='vertical')
	plt.yticks(range(error_bar.size),all_keys)
	plt.colorbar()

	plt.matshow(confusion_matrix[:-7,:-7], cmap='hot')
	plt.xticks(range(error_bar.size-7),all_keys[:-7], rotation='vertical')
	plt.yticks(range(error_bar.size-7),all_keys[:-7])
	plt.colorbar()

	df = pd.DataFrame(data=confusion_matrix, index=all_keys, columns=all_keys)
	df.to_csv('confusion_matrix.csv')

	plt.matshow(neuro_matrix, cmap='YlOrRd_r')
	plt.xticks(range(error_bar.size),all_keys, rotation='vertical')
	plt.yticks(range(3), ['vgat', 'vglut'])
	plt.colorbar()

	plt.matshow(neuro_matrix[:,:-7], cmap='YlOrRd_r')
	plt.xticks(range(error_bar.size-7),all_keys[:-7], rotation='vertical')
	plt.yticks(range(3), ['vgat', 'vglut'])
	plt.colorbar()

	plt.show()
	print 'neuro_stats:'
	print neuro_matrix[0,:]
	print neuro_matrix[1,:]
