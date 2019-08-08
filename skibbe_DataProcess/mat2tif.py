import os
import numpy as np
import scipy.io
from queue import Queue
from skimage.external import tifffile

def readmat(file_path):
    """ .mat return img, gt, connection
        return:
            img: (512, 512, 512)
            connection: (6, M)
            coords: (3, N)
            labels: (N, )
    """
    mat = scipy.io.loadmat(file_path)
    
    img = mat['img']
    imgD = mat['imgD']
    GT = mat['gt']  # 0,1,2 行: position
                    # 20: pathid
    connection = mat['con']  # 0,1: connectionnection
    
    coords = np.minimum(GT[:3, :], 299)
    labels = GT[20, :].astype(np.int)
    return img, imgD, connection, coords, labels

def decision_end_cross(connection):
    """ 有前后连接关系 分析出 端点 和 交叉点
        return: (8, M), 0,1 is the same, 6,7中-1端点, -2交叉点
    """
    connection_e_c = []
    for i in range(connection.shape[1]):
        head = np.where(connection[0, :]==connection[0, i])[0]
        if len(head) > 1:
            stateH = -2
        elif len(np.where(connection[1, :]==connection[0, i])[0]) == 0:
            stateH = -1
        else:
            stateH = 0
        
        end = np.where(connection[0, :]==connection[1, i])[0]
        if len(end) == 0:
            stateE = -1
        elif len(end) > 1:
            stateE = -2
        else:
            stateE = 0
        temp = connection[:, i].tolist() + [stateH, stateE]
        connection_e_c.append(temp)
    return np.array(connection_e_c).T

def gene_SwcsList(coords, con_e_c, labels):
    swcs_list = []
    cnt = 0

    curr_id = 1
    label = labels[cnt]
    coord = coords[:, 0].tolist()[::-1]
    pointstate = -1
    last_id = -1

    temp = []
    temp += [curr_id, label]
    temp += coord
    temp += [pointstate, last_id]

    swc_list = []
    swc_list.append(temp)
    cnt += 1

    # 转为.swc 格式
    cross = Queue()
    for i, C in enumerate(con_e_c.T[0:, :], 0):
        
    #     if C[0] != con[i-1, 1] and len(np.where(con[1, :]==C[0])[0]) > 0:
    #         continue
    #     elif C[0] != con[i-1, 1] and len(np.where(con[1, :]==C[0])[0]) == 0:
    #         swcs_list.append(swc_list)
    #         swc_list = []
        temp = []
        if i != 0 and C[0] != con_e_c[1, i-1]:
            if len(np.where(con_e_c[1, :]==C[0])[0]) == 0:
                swcs_list.append(swc_list)
            
                swc_list = []
                curr_id = 1
                last_id = -1
                label = labels[cnt]
                coord = coords[:, C[0]].tolist()[::-1]
                pointstate = C[6]  # -1
                
                temp += [curr_id, label]
                temp += coord
                temp += [pointstate, last_id]
                
                swc_list.append(temp)
                cnt += 1
                temp = []
            else:
                if C[0] == con_e_c[0, i-1]:  # 分支第一个点
                    cross.put(last_id)
                    continue
                else:
                    if len(np.where(con_e_c[0, :i-1]==con_e_c[0, i-1])[0]) == 0:  # 紧接着交叉点的位置
                        #C[1]
                        last_id = cross.get()
                        curr_id += 1
                        label = labels[i]
                        coord = coords[:, C[0]].tolist()[::-1]
                        pointstate = C[6]  # 0

                        temp += [curr_id, label]
                        temp += coord
                        temp += [pointstate, last_id]
                        swc_list.append(temp)
                        cnt += 1
                        temp = []


        last_id = curr_id
        curr_id += 1
        label = labels[cnt]
        coord = coords[:, C[1]].tolist()[::-1]
        pointstate = C[7]  # 0
        
        temp += [curr_id, label]
        temp += coord
        temp += [pointstate, last_id]
        swc_list.append(temp)
        cnt += 1

    swcs_list.append(swc_list)
    return swcs_list

def gene_SwcsList_v1(coords, con_e_c, labels):
    swcs_list = []

    curr_id = 1
    label = labels[0]
    coord = coords[:, 0].tolist()[::-1]
    pointstate = -1
    last_id = -1

    temp = []
    temp += [curr_id, label]
    temp += coord
    temp += [pointstate, last_id]

    swc_list = []
    swc_list.append(temp)

    # 转为.swc 格式
    cross = Queue()
    for i, C in enumerate(con_e_c.T[0:, :], 0):
        
    #     if C[0] != con[i-1, 1] and len(np.where(con[1, :]==C[0])[0]) > 0:
    #         continue
    #     elif C[0] != con[i-1, 1] and len(np.where(con[1, :]==C[0])[0]) == 0:
    #         swcs_list.append(swc_list)
    #         swc_list = []
        temp = []
        if i != 0 and C[0] != con_e_c[1, i-1]:
            if len(np.where(con_e_c[1, :]==C[0])[0]) == 0:
                swcs_list.append(swc_list)
            
                swc_list = []
                curr_id = 1
                last_id = -1
                label = labels[C[0]]
                coord = coords[:, C[0]].tolist()[::-1]
                pointstate = C[6]  # -1
                
                temp += [curr_id, label]
                temp += coord
                temp += [pointstate, last_id]
                
                swc_list.append(temp)
                temp = []
                cross = Queue()

            else:
                if C[0] == con_e_c[0, i-1]:  # 分支第一个点
                    if C[1] != con_e_c[0, i+1]:
                        cross.put(last_id)
                        continue
                    else:
                        last_id = last_id - len(np.where(con_e_c[0, :i]==C[0])[0]) + 1
                        curr_id += 1
                        label = labels[C[1]]
                        coord = coords[:, C[1]].tolist()[::-1]
                        pointstate = C[7]
                        
                        temp += [curr_id, label]
                        temp += coord
                        temp += [pointstate, last_id]
                        swc_list.append(temp)
                        temp = []
                else:
                    if len(np.where(con_e_c[0, :i-1]==con_e_c[0, i-1])[0]) == 0:  # 紧接着交叉点的位置
                        #C[1]
                        last_id = cross.get()
                        curr_id += 1
                        label = labels[C[0]]
                        coord = coords[:, C[0]].tolist()[::-1]
                        pointstate = C[6]  # 0

                        temp += [curr_id, label]
                        temp += coord
                        temp += [pointstate, last_id]
                        swc_list.append(temp)
                        temp = []


        last_id = curr_id
        curr_id += 1
        label = labels[C[1]]
        coord = coords[:, C[1]].tolist()[::-1]
        pointstate = C[7]  # 0
        
        temp += [curr_id, label]
        temp += coord
        temp += [pointstate, last_id]
        swc_list.append(temp)

    swcs_list.append(swc_list)
    return swcs_list

def swcslist2txt(swcs_str_list, file_dir, name):
    for i, swc_str_list in enumerate(swcs_str_list, 1):
        file_name = name.split('.')[0] + "_{:d}.swc".format(i)
        file = open(os.path.join(file_dir, file_name), 'w')
        for j, swc_str in enumerate(swc_str_list):
            file.write(swc_str + '\n')
        file.close()

def create_gt(coords, shape, labels):
    gt = np.zeros(shape)
    coords_int = np.round(coords).astype(np.int)
    gt[coords_int[0, :], coords_int[1, :], coords_int[2, :]] = labels
    return gt

def gene_end_cross_mask(c_e_c, coords, shape):
    """ 将endpoint: 1, cross: 2 记录为mask 的形式
    """
    mask = np.zeros(shape)
    coords_int = np.round(coords).astype(np.int)
    for C in c_e_c.T:
        if C[6] == -1:
            mask[coords_int[0, C[0]], coords_int[1, C[0]], coords_int[2, C[0]]] = 1
        elif C[6] == -2:
            mask[coords_int[0, C[0]], coords_int[1, C[0]], coords_int[2, C[0]]] = 2
        if C[7] == -1:
            mask[coords_int[0, C[1]], coords_int[1, C[1]], coords_int[2, C[1]]] = 1
        elif C[7] == -2:
            mask[coords_int[0, C[1]], coords_int[1, C[1]], coords_int[2, C[1]]] = 2
    return mask


# file_path = r"C:\Users\Administrator\Desktop\test_neural\testSyn\chaos\chaos.mat"
# path, name = os.path.split(file_path)

# # 读取.mat
# img, con, coords, labels = readmat(file_path)

# # 交叉点, 端点
# con_e_c = decision_end_cross(con)

# # 得到swc格式的list()
# swcslist = gene_SwcsList(coords, con_e_c, labels)
# swcs_str_list = [[str(l).strip('[').strip(']') for l in i] for i in swcslist]

# # 读写.swc
# swcslist2txt(swcs_str_list, path, name)

def process_single_mat(file_path, out_dir):
    path, name = os.path.split(file_path)

    # 读取.mat
    img, imgD, con, coords, labels = readmat(file_path)

    # 创建.tif 
    tifpath = os.path.join(out_dir, "tiffs")
    if not os.path.exists(tifpath):
        os.makedirs(tifpath)
    tifffile.imsave(os.path.join(tifpath, name.split('.')[0]+'.tif'), img)

    noisepath = os.path.join(out_dir, "noises")
    if not os.path.exists(noisepath):
        os.makedirs(noisepath)
    tifffile.imsave(os.path.join(noisepath, name.split('.')[0]+'.tif'), imgD)


    gt = create_gt(coords, img.shape, labels)
    gtpath = os.path.join(out_dir, "gts")
    if not os.path.exists(gtpath):
        os.makedirs(gtpath)
    tifffile.imsave(os.path.join(gtpath, name.split('.')[0]+'.tif'), gt.astype(np.int16))

    # 交叉点, 端点
    con_e_c = decision_end_cross(con)

    # 生成end_cross_msk
    ec_mask = gene_end_cross_mask(con_e_c, coords, img.shape)
    ec_maskpath = os.path.join(out_dir, "end_cross_mask")
    if not os.path.exists(ec_maskpath):
        os.makedirs(ec_maskpath)
    tifffile.imsave(os.path.join(ec_maskpath, name.split('.')[0]+'.tif'), ec_mask.astype(np.int16))

    # 得到swc格式的list()
    swcslist = gene_SwcsList_v1(coords, con_e_c, labels)
    swcs_str_list = [[str(l).strip('[').strip(']').replace(', ', ' ') for l in i] for i in swcslist]

    # 读写.swc
    swcpath = os.path.join(out_dir, "swcs")
    if not os.path.exists(swcpath):
        os.makedirs(swcpath)
    swcslist2txt(swcs_str_list, swcpath, name)

    print("{} have benn processed".format(name))


mats_dir = "/media/fcheng/skibbe-dense_axon_sim/data/synthetic/chaos/300_0-5_2_4_6/"
swcs_name = os.listdir(mats_dir)
swcs_name.sort()
# print(swcs_name)

output_dir = "/media/fcheng/synthetic/SyntheticData/chaos/300_0-5_2_4_6/"
for swc_name in swcs_name:
    out_dir = os.path.join(output_dir, swc_name.split('.')[0])  # .../Syn_chaos_300_25_2_000/
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    process_single_mat(os.path.join(mats_dir, swc_name), out_dir)