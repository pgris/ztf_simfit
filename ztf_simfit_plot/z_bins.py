from ztf_pipeutils.ztf_pipeutils.ztf_hdf5 import Read_LightCurve
import statistics as stat
import numpy as np
import operator
from scipy import interpolate
import os
import astropy


class Apply_mask:
    "class which affect mask on meta data"

    def __init__(self, metaFitInput='Meta_fit.hdf5', inputDir='dataLC',
                 var={'sel': 'sel', 'c_err': 'c_err',
                      'zmin': 'z', 'zmax': 'z', 'chi2': 'chi2'},
                 op={'op sel': 'operator.eq', 'op c_err': 'operator.ne', 'op zmin': 'operator.ge',
                     'op zmax': 'operator.lt', 'op chi2': 'operator.le'},
                 lim={'lim sel': 1, 'lim c_err': -1, 'lim zmin': 0.01, 'lim zmax': 0.1,
                      'lim chi2': 2}):
        """
        Parameters
        --------------
        metaFitInput: str, opt
          input file (default: 'Meta_fit.hdf5')
        inputDir : str, opt
            input directory to find the file (default='dataLC').
        var : dict, opt
            dictionary with list of variable that you want to apply mask on it
        op : dict, opt
            dictionary with list of operator for the mask
        lim : dict, opt
            dictionary with lim

        Examples
        --------
        You want to apply a mask to your metadata : mk = metadata['z']<0.5
            So, you put 'z' on var dictionnary.
                you put 'operator.lt' on op dictionary.
                you put '0.5' on lim dictionary.
        """

        self.cl = Read_LightCurve(file_name=metaFitInput, inputDir=inputDir)
        self.metadata = self.cl.get_table(path='meta')
        self.dico = {'variables': var, 'operators': op, 'limites': lim}

    def meta_add_chi2(self):
        md = self.metadata
        md['chi2'] = self.metadata['chisq']/self.metadata['ndof']
        return md

    def mask_to_apply(self):
        md_ = self.meta_add_chi2()
        lim = [lim for lim in self.dico['limites'].values()]
        op = [op for op in self.dico['operators'].values()]
        var = [var for var in self.dico['variables'].values()]

        for i, vr in enumerate(var):
            op_ = eval(op[i])
            mask = op_(md_[vr], lim[i])
            md_ = md_[mask]
        return md_


class Z_bins:
    "class which calculate the c_err per bins of z"

    def __init__(self, metaFitInput='Meta_fit.hdf5', inputDir='dataLC', dico={'variables': {'sel': 'sel',
                                                                                            'c_err': 'c_err',
                                                                                            'zmin': 'z',
                                                                                            'zmax': 'z',
                                                                                            'chi2': 'chi2'},
                                                                              'operators': {'op sel': 'operator.eq',
                                                                                            'op c_err': 'operator.ne',
                                                                                            'op zmin': 'operator.ge',
                                                                                            'op zmax': 'operator.lt',
                                                                                            'op_chi2': 'operator.lt'},
                                                                              'limites': {'lim sel': 1,
                                                                                          'lim c_err': -1,
                                                                                          'lim zmin': 0.01,
                                                                                          'lim zmax': 0.1,
                                                                                          'lim_chi2': 2}}, zmin=0.01, zmax=0.2):

        self.zmin = zmin
        self.zmax = zmax
        self.cl = Apply_mask(metaFitInput=metaFitInput, inputDir=inputDir, var=dico['variables'], op=dico['operators'],
                             lim=dico['limites'])
        self.dico = self.cl.dico
        self.inputDir = inputDir
        self.metaFitInput = metaFitInput

    def moy_z_bin(self):

        z = np.arange(self.zmin, self.zmax, 0.01)
        z_moy, z_RMS, c_moy, c_RMS = [], [], [], []

        for i in range(0, len(z)-1):
            zmin_, zmax_ = np.round(z[i], 2), np.round(z[i+1], 2)
            self.dico['limites']['lim zmin'] = zmin_
            self.dico['limites']['lim zmax'] = zmax_

            md = self.cl.mask_to_apply()

            z_moy.append(np.mean(md['z']))
            z_RMS.append(np.std(md['z']))
            c_moy.append(np.mean(md['c_err']))
            c_RMS.append(1/(np.sqrt(np.sum(1/md['c_err']**2))))

        return z_moy, z_RMS, c_moy, c_RMS

    def interpolate1d(self, x, y):
        f = interpolate.interp1d(x, y, bounds_error=False, fill_value=0.)
        x_new = x
        y_new = f(x_new)
        z_comp = f(0.04)
        return x_new, y_new, z_comp

    def get_z(self):

        z_moy, z_RMS, c_moy, c_RMS = self.moy_z_bin()
        sup_list, inf_list = [], []
        for i in range(0, len(z_moy)):
            sup_list.append(c_moy[i] + c_RMS[i])
            inf_list.append(c_moy[i] - c_RMS[i])

        x_new, y_new, z_comp = self.interpolate1d(c_moy, z_moy)
        x_sup, y_sup, z_comp_sup = self.interpolate1d(sup_list, z_moy)
        x_inf, y_inf, z_comp_inf = self.interpolate1d(inf_list, z_moy)

        return z_comp, z_comp_sup, z_comp_inf
    
    def plot_z_comp(self):
        
        import matplotlib.pylab as plt
        z_moy, z_RMS, c_moy, c_RMS = self.moy_z_bin()
        z_comp, z_comp_sup, z_comp_inf = self.get_z()
        sup_list, inf_list = [], []
        for i in range(0, len(z_moy)):
            sup_list.append(c_moy[i] + c_RMS[i])
            inf_list.append(c_moy[i] - c_RMS[i])
        
        plt.figure(figsize=(9,7))
        #plt.plot(z_moy, c_moy, 'o', y_new, x_new, '-', color='orange')
        plt.plot(z_moy, c_moy, 'o', color='orange')
        plt.xlabel('$z$', fontsize=13)
        plt.ylabel('$\sigma_{c}$', fontsize=13)

        plt.text(0.06, 0.05, r'$z_{completeness}$ ='+'{}'.format(
                np.round(z_comp, 3)) + r'$\pm$' + '{}'.format(np.round(z_comp-z_comp_sup, 3)))
        plt.text(0.015, 0.042, r'$\sigma_{c}=0.04$', color='red')
        plt.axhline(y=0.04, color='red', linestyle='-')
        plt.axvline(x=z_comp, color='red', linestyle='-')
        plt.fill_between(z_moy, sup_list, inf_list,
                         color='gold', alpha=0.3)
        #plt.errorbar(z_moy, c_moy, yerr=c_RMS,
                     #fmt='none', capsize=2, elinewidth=1, color='red')
        plt.show()
        
    def plot_err_c_z(self, axhline=False, fontsize=15, error_bar=False, color='red'):
        
        import matplotlib.pylab as plt
        z_moy, z_RMS, c_moy, c_RMS = self.moy_z_bin()
        plt.scatter(z_moy, c_moy)
        plt.xlabel('$z$', fontsize=fontsize)
        plt.ylabel('$\sigma_c$', fontsize=fontsize)
        if axhline:
            plt.axhline(y=0.04, color='r', linestyle='-')
        if error_bar:
            plt.errorbar(z_moy, c_moy, yerr=c_RMS,
                         fmt='none', capsize=2, elinewidth=1, color=color)

    def add_data(self, metaDirOutput, metaFileOutput):

        cl = Read_LightCurve(file_name=self.metaFitInput,
                             inputDir=self.inputDir)
        metadata = cl.get_table(path='meta')
        new_meta = metadata.copy()

        z_comp, z_comp_sup, z_comp_inf = self.get_z(
            add_column=True)

        new_meta['z_comp'] = z_comp
        new_meta['z_comp_sup'] = z_comp_sup
        new_meta['z_comp_inf'] = z_comp_inf

        if not os.path.exists(metaDirOutput):
            os.makedirs(metaDirOutput)

        fOut = '{}/{}'.format(metaDirOutput, metaFileOutput)
        astropy.io.misc.hdf5.write_table_hdf5(
            new_meta, fOut, path='meta', overwrite=True, serialize_meta=False)
