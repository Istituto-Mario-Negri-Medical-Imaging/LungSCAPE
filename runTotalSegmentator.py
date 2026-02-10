# Training_Lobes_preproc

import os
import nibabel as nib
import numpy as np

in_folder_path = '...'
out_folder_path = '...'

all_files = sorted(os.listdir(in_folder_path))
for i in range(0, len(all_files), 1):
	ct_file = all_files[i]

	# Extract the patient ID from the CT file name (assuming consistent naming convention)
	patient_id = ct_file.split('.')[0]  # Removes the file extension to get patient code,
	new_name = patient_id

	# Construct full file paths
	ct_path = os.path.join(in_folder_path, ct_file)
	output = os.path.join(out_folder_path, new_name)
	
	os.system(' '.join(['TotalSegmentator -i', ct_path, '-o', output, '--roi_subset lung_upper_lobe_left lung_lower_lobe_left lung_upper_lobe_right lung_middle_lobe_right	lung_lower_lobe_right']))
	
	os.system(' '.join(['TotalSegmentator -i', ct_path, '-o', output, '--ta lung_vessels']))


