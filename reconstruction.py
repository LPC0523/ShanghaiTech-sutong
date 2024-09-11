import astra
import numpy as np
from skimage.transform import radon, rescale

class Paramaters:
    def __init__(self):
        self.param = {}
        self.param['nx'] = 640  # image width
        self.param['ny'] = 640  # image height
        self.param['dect_count'] = 512  # number of detectors
        self.param['dsd'] = 1500  # distance from source to detector
        self.param['dso'] = 1000  # distance from source to object
        self.param['nProj'] = 720  # number of projection views
        self.param['startangle'] = 0  # start angle
        self.param['endangle'] = 2*np.pi  # end angle
        self.param['detector_width'] = 1  # detector spacing
        self.param['algorithm'] = 'FBP_CUDA'  # reconstruction algorithm SIRT
        self.param['interation'] = -1  # interations, only used in interative reconstruction algorithms
        self.param['short_scan'] = False  # park compensation, only used in fanbeam and conebeam
        self.param['noise'] = False
        self.reuse = False


def add_noise(fp, photon=10**5):
    sino = fp/100
    np.random.seed(42)
    proj_exp = np.maximum(photon * np.exp(-sino), 1)
    guass = np.random.normal(0, 100, size=np.shape(proj_exp))
    proj_exp = np.random.poisson(proj_exp) + guass
    proj_data = -np.log(np.minimum(proj_exp / photon, 1))
    return proj_data*100


def create_sinogram_parallel(image, param):
    angles = np.linspace(param['startangle'], param['endangle'],
                         param['nProj'], endpoint=False)
    sino = get_parallel_sino(image, dect_w=param['detector_width'],
                        dect_count=param['dect_count'], vol_geom_size=(param['nx'], param['ny']), angles=angles)
    return sino


def create_sinogram_parallel2(image, angles):
    sinogram = radon(image, theta=angles, circle=False).T
    return sinogram


def get_parallel_sino(image, dect_w, dect_count, vol_geom_size, angles, ):
    vol_geom_par = astra.create_vol_geom(vol_geom_size)
    proj_geom_par = astra.create_proj_geom('parallel', dect_w, dect_count,
                                           angles)

    proj_par_id = astra.create_projector('cuda', proj_geom_par, vol_geom_par)
    sino_par_id, sino_gram = astra.create_sino(image, proj_par_id)

    astra.projector.delete(proj_par_id)
    astra.projector.delete(sino_par_id)
    astra.clear()
    return sino_gram


def recon_parallel_param(sino, param):
    angles = np.linspace(param['startangle'], param['endangle'],
                         param['nProj'], endpoint=False)
    image = recon_parallel(alg=param['algorithm'], sino=sino,
                      dect_w=param['detector_width'], vol_geom_size=(param['nx'], param['ny']), angles=angles,
                      interations=param['interation'])
    return image


def recon_parallel(alg, sino, dect_w, vol_geom_size,
                   angles, interations=-1, filter="ram-lak"):
    astra.algorithm.clear()
    vol_geom_parallel = astra.create_vol_geom(vol_geom_size)
    proj_geom_parallel = astra.create_proj_geom('parallel', dect_w, sino.shape[1],
                                                angles)
    proj_parallel_id = astra.create_projector('cuda', proj_geom_parallel, vol_geom_parallel)
    sinogram_parallel_id = astra.data2d.create('-sino', proj_geom_parallel, sino)
    rec_parallel_id = astra.data2d.create('-vol', vol_geom_parallel)
    cfg_parallel = astra.astra_dict(alg)
    cfg_parallel['ReconstructionDataId'] = rec_parallel_id
    cfg_parallel['ProjectionDataId'] = sinogram_parallel_id
    cfg_parallel['ProjectorId'] = proj_parallel_id
    if alg.startswith("FBP"):
        cfg_parallel["FilterType"] = filter
    alg_parallel_id = astra.algorithm.create(cfg_parallel)
    if interations != -1:
        astra.algorithm.run(alg_parallel_id, interations)
    else:
        astra.algorithm.run(alg_parallel_id)
    rec_parallel = astra.data2d.get(rec_parallel_id)
    astra.algorithm.delete(alg_parallel_id)
    astra.data2d.delete(rec_parallel_id)
    astra.data2d.delete(sinogram_parallel_id)
    astra.projector.delete(proj_parallel_id)
    # rec_fan[rec_fan < 0] = 0
    # rec_fan = np.flipud(rec_fan)
    astra.algorithm.clear()
    astra.clear()
    return rec_parallel


def get_fan_sino_param(image, param):
    angles = np.linspace(param['startangle'], param['endangle'], param['nProj'], endpoint=False)
    sino = get_fan_sino(image, source_ori=param['dso'], ori_detector=param['dsd'] - param['dso'],
                        dect_w=param['detector_width'],
                        dect_count=param['dect_count'], vol_geom_size=(param['nx'], param['ny']), angles=angles, )
    return sino


def recon_fan_param(sino, param):
    angles = np.linspace(param['startangle'], param['endangle'],
                         param['nProj'], endpoint=False)
    image = recon_fan(alg=param['algorithm'], sino=sino, source_ori=param['dso'],
                      ori_detector=param['dsd'] - param['dso'],
                      dect_w=param['detector_width'], vol_geom_size=(param['nx'], param['ny']), angles=angles,
                      interations=param['interation'], short_scan=param['short_scan'])
    return image


def get_fan_sino(image, source_ori, ori_detector, dect_w, dect_count, vol_geom_size, angles, ):
    vol_geom_fan = astra.create_vol_geom(vol_geom_size)
    proj_geom_fan = astra.create_proj_geom('fanflat', dect_w, dect_count,
                                           angles,
                                           source_ori, ori_detector)

    proj_fan_id = astra.create_projector('cuda', proj_geom_fan, vol_geom_fan)
    sino_fan_id, sino_gram = astra.create_sino(image, proj_fan_id)

    astra.projector.delete(proj_fan_id)
    astra.projector.delete(sino_fan_id)
    astra.clear()
    return sino_gram


def recon_fan(alg, sino, source_ori, ori_detector, dect_w, vol_geom_size,
              angles, interations=-1, short_scan=False, filter="ram-lak"):
    astra.algorithm.clear()
    vol_geom_fan = astra.create_vol_geom(vol_geom_size)
    proj_geom_fan = astra.create_proj_geom('fanflat', dect_w, sino.shape[1],
                                           angles,
                                           source_ori, ori_detector)
    proj_fan_id = astra.create_projector('cuda', proj_geom_fan, vol_geom_fan)
    sinogram_fan_id = astra.data2d.create('-sino', proj_geom_fan, sino)
    rec_fan_id = astra.data2d.create('-vol', vol_geom_fan)
    cfg_fan = astra.astra_dict(alg)
    cfg_fan['ReconstructionDataId'] = rec_fan_id
    cfg_fan['ProjectionDataId'] = sinogram_fan_id
    # cfg_fan['ProjectorId'] = proj_fan_id
    cfg_fan['option'] = {'ShortScan': short_scan}
    if alg.startswith("FBP"):
        cfg_fan["FilterType"] = filter
    alg_fan_id = astra.algorithm.create(cfg_fan)
    if interations != -1:
        astra.algorithm.run(alg_fan_id, interations)
    else:
        astra.algorithm.run(alg_fan_id)
    rec_fan = astra.data2d.get(rec_fan_id)
    astra.algorithm.delete(alg_fan_id)
    astra.data2d.delete(rec_fan_id)
    astra.data2d.delete(sinogram_fan_id)
    astra.projector.delete(proj_fan_id)
    # rec_fan[rec_fan < 0] = 0
    # rec_fan = np.flipud(rec_fan)
    astra.algorithm.clear()
    astra.clear()
    return rec_fan


def recon_cone_param(projections, param):
    angles = np.linspace(param['startangle'], param['endangle'],
                         param['nProj'], endpoint=False)
    image = recon_cone(alg=param['algorithm'], projections=projections, source_ori=param['dso'],
                       ori_detector=param['dsd'] - param['dso'], nh=param['nh'], nw=param['nw'],
                       det_spacing=param['det_spacing'], vol_geom_size=(param['ny'], param['nx'], param['nz']),
                       angles=angles,
                       interations=param['interation'], short_scan=param['short_scan'])
    return image


def recon_cone(alg, projections, source_ori, ori_detector, det_spacing, vol_geom_size,
               nh, nw,
               angles, interations=-1, short_scan=False):
    vol_geom = astra.creators.create_vol_geom(vol_geom_size)

    proj_geom = astra.create_proj_geom('cone', det_spacing[0], det_spacing[1], nh, nw,
                                       angles, source_ori, ori_detector)
    projections_id = astra.data3d.create('-proj3d', proj_geom, projections)
    reconstruction_id = astra.data3d.create('-vol', vol_geom, data=0)
    cfg_cone = astra.astra_dict(alg)
    cfg_cone['ProjectionDataId'] = projections_id
    cfg_cone['ReconstructionDataId'] = reconstruction_id
    cfg_cone['option'] = {'ShortScan': short_scan}

    algorithm_id = astra.algorithm.create(cfg_cone)
    astra.algorithm.run(algorithm_id)

    # 重建volume舍负值
    result = astra.data3d.get(reconstruction_id)
    # result = np.flipud(result)
    # reconstruction[reconstruction < 0] = 0
    astra.algorithm.clear()
    astra.clear()

    return result
