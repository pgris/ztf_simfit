from ztf_pipeutils.ztf_hdf5 import Read_LightCurve
import statistics as stat
import matplotlib.pylab as plt
import numpy as np
import operator as op


class Efficiency:
    "class which calculate the c_err per bins of z"

    def __init__(self, meta_Input='Meta_fit.hdf5', inputDir='dataLC'):

        self.cl = Read_LightCurve(file_name=meta_Input, inputDir=inputDir)
        self.metadata = self.cl.get_table(path='meta')

    def mask_to_apply(self, metadata, var, op, lim):
        for i, vr in enumerate(var):
            mask = op[i](metadata[vr], lim[i])
            metadata = metadata[mask]
        return metadata

    def z_efficiency(self):
        eff_rg, Z_rg = [], []
        eff_rgi, Z_rgi = [], []
        eff_, Z_ = [], []
        z = np.arange(np.min(self.metadata['z']), np.max(
            self.metadata['z']), 0.01)

        for i in range(0, len(z)-1):
            zmin = np.round(z[i], 2)
            zmax = np.round(z[i+1], 2)

            md_rg = self.mask_to_apply(self.metadata, var=['z', 'z', 'n_i_band', 'sel'], op=[
                                       op.ge, op.lt, op.eq, op.eq], lim=[zmin, zmax, 0, 1])
            md_rgi = self.mask_to_apply(self.metadata, var=['z', 'z', 'n_i_band', 'sel'], op=[
                                        op.ge, op.lt, op.ne, op.eq], lim=[zmin, zmax, 0, 1])

            # sans distinction
            md_ = self.mask_to_apply(self.metadata, var=['z', 'z', 'sel'], op=[
                                     op.ge, op.lt, op.eq], lim=[zmin, zmax, 1])
            md_tot = self.mask_to_apply(self.metadata, var=['z', 'z'], op=[
                                        op.ge, op.lt], lim=[zmin, zmax])

            # rg
            E_rg = len(md_rg)/len(md_tot)
            eff_rg.append(E_rg)
            Z_rg.append(np.round(z[i], 2))

            # rgi
            E_rgi = len(md_rgi)/len(md_tot)
            eff_rgi.append(E_rgi)
            Z_rgi.append(np.round(z[i], 2))

            # sans distinction
            E_ = len(md_)/len(md_tot)
            eff_.append(E_)
            Z_.append(np.round(z[i], 2))

        plt.figure(figsize=(8, 6))
        plt.plot(Z_rg, eff_rg, label='lc with r,g bands')
        plt.plot(Z_rgi, eff_rgi, label='lc with r,g,i bands')
        plt.plot(Z_, eff_, label='for all lc')
        plt.xlabel('$z$', fontsize=15)
        plt.ylabel('Efficiency', fontsize=15)

        plt.legend()
