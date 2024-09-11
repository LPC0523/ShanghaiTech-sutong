from reconstruction import Paramaters, recon_fan_param
import matplotlib.pyplot as plt
from scipy.ndimage import shift
from os.path import join as ospj
import tifffile as tif
import numpy as np


# 读取BMP图片并取得投影图像
def read_bmp_file(file_path: str):
    """
    Reads a BMP file and displays a preview of the image.

    Args:
        file_path (str): Path to the BMP file.

    Returns:
        None
    """
    try:
        # Read the binary file
        with open(file_path, 'rb') as fid:
            a = np.fromfile(fid, dtype=np.uint16)

        # Remove the BMP header info (first 54 bytes)
        data_img = a[27:]

        # Reshape data into a 2D array
        img = data_img.reshape((480, 640))

        return img
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")


# 从投影图像构建弦图
def get_sinogram(Parm, data_path, ref_path):
    sinogram = np.zeros((Parm.param['nProj'], Parm.param['dect_count']))
    Blank_row = np.zeros((1, Parm.param['dect_count']))
    num_projections = Parm.param['nProj']

    # 将所有亮场投影求和取平均，得到预处理用的亮场行
    for i in range(num_projections):
        filename = '{:03d}.bmp'.format(i)
        blank_image = read_bmp_file(ref_path + filename)
        blank_row = blank_image[blank_image.shape[0] // 2]
        Blank_row += blank_row
    Blank_row = Blank_row / num_projections + 1e-8

    # 按顺序将投影图取中间行，预处理并拼接为弦图
    for i in range(num_projections):
        filename = '{:03d}.bmp'.format(i)
        projection_image = read_bmp_file(data_path + filename)
        projection_row = projection_image[projection_image.shape[0]//2]

        # Blank field correction
        corrected_row =  np.abs(np.log(Blank_row)- np.log(projection_row))

        # Put corrected data into the corresponding column of sinogram
        sinogram[i, :] = corrected_row

    # 对探测器偏移的影响进行校正
    corrected_sinogram = shift(sinogram, [8, -5], cval=0)
    plt.figure()
    plt.imshow(corrected_sinogram, cmap='gray')
    return  corrected_sinogram



if __name__ == "__main__":
    Parm = Paramaters()

    # 设置重建参数，其中dso、dsd和detector_width参数的单位是图像像素个数，因此要除以图像像素大小
    Parm.param['nProj'] = 800
    Parm.param['dect_count'] = 640
    Parm.param['dsd'] = (476.90+45)/(0.006*476.90/(476.90+45))
    Parm.param['dso'] = 476.90/(0.006*476.90/(476.90+45))
    Parm.param['detector_width'] = (476.90+45)/476.90

    # 数据路径，根据自己的数据存放位置进行修改
    for i in range(3):
        data_path = './Linzihao'+str(i+1)+'/ScanData/'
        ref_path = './Linzihao'+str(i+1)+'/ScanRef/'
        sinogram = get_sinogram(Parm, data_path, ref_path)
        image = recon_fan_param(sinogram, Parm.param)
        file_name_sino = ('Sinogram_' + str(i + 1) + '.tif')
        file_name_im = ('Image_' + str(i+1) + '.tif')
        tif.imwrite(ospj('result/', file_name_sino), sinogram)
        tif.imwrite(ospj('result/', file_name_im), image)

        # 可视化结果
        plt.figure()
        plt.imshow(image, cmap='gray')
        plt.colorbar()
        plt.show()