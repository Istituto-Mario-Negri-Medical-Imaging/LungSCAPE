import os
import nibabel as nib
import numpy as np

folder_path = '/mnt/sdb1/lungs_longitudinale/Lobes_Segmentations_TS'
output_path = '/mnt/sdb1/lungs_longitudinale/Lobes_Segmentations_TS_unique'

# Iterate over subfolders in the specified folder
for subfolder in os.listdir(folder_path):
	subfolder_path = os.path.join(folder_path, subfolder)

	# Check if it's a directory
	if os.path.isdir(subfolder_path):

		file1 = 'lung_upper_lobe_left.nii.gz'
		idx1 = 1
		path1 = os.path.join(subfolder_path, file1)
		
		file2 = 'lung_lower_lobe_left.nii.gz'
		idx2 = 2
		path2 = os.path.join(subfolder_path, file2)
		
		file3 = 'lung_upper_lobe_right.nii.gz'
		idx3 = 3
		path3 = os.path.join(subfolder_path, file3)
		
		file4 = 'lung_middle_lobe_right.nii.gz'
		idx4 = 4
		path4 = os.path.join(subfolder_path, file4)

		file5 = 'lung_lower_lobe_right.nii.gz'
		idx5 = 5
		path5 = os.path.join(subfolder_path, file5)
	
		# Save the aggregated image
		output_file = os.path.join(output_path, f"{subfolder}_lobesTS.nii.gz")

		affine = nib.load(path1).affine  # Use affine from the first image
		img1 = nib.load(path1).get_fdata().astype('uint8')
		img2 = nib.load(path2).get_fdata().astype('uint8')
		img3 = nib.load(path3).get_fdata().astype('uint8')
		img4 = nib.load(path4).get_fdata().astype('uint8')	
		img5 = nib.load(path5).get_fdata().astype('uint8')

		img1[img2!=0]=2
		img1[img3!=0]=3
		img1[img4!=0]=4
		img1[img5!=0]=5

		aggregated_nifti = nib.Nifti1Image(img1, affine)
		nib.save(aggregated_nifti, output_file)

		print(f"Aggregated image for subfolder '{subfolder}' saved to '{output_file}'.")


