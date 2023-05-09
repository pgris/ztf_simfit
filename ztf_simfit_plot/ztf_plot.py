from ztf_pipeutils.ztf_hdf5 import Read_LightCurve
from ztf_simfit.ztf_fit import SN_fit
import sncosmo


class VisuLC:
    def __init__(self, metaFileInput, metaDirInput, SNR=5):

        meta = Read_LightCurve(file_name=metaFileInput, inputDir=metaDirInput)
        metaTable = meta.get_table(path='meta')

        metadata = metaTable.meta
        # get lc
        lcDir = metadata['directory']
        lcName = metadata['file_name']

        self.lcs = Read_LightCurve(file_name=lcName, inputDir=lcDir)
        self.SNR = SNR
        
        # print SNIDS

        print(metaTable['path'])

    def plot(self, lcpath):

        lc = self.lcs.get_table(lcpath)

        # trying to fit here
        fit = SN_fit(lc)
        result, fitted_model = fit('sncosmo', plot=True)
