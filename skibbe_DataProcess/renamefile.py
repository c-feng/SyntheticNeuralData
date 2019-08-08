import os

files_path = "/media/fcheng/skibbe-dense_axon_sim/data/synthetic/chaos/300_0.5_2_4_6/"

files_name = os.listdir(files_path)
files_name.sort()

# root = os.getcwd()
for file_name in files_name:
    suffix = file_name.split('.')[-1]
    name = file_name.replace('.'+suffix, '')

    new_name = name.replace('.', '-') + '.' + suffix

    os.rename(files_path+file_name, files_path+new_name)
