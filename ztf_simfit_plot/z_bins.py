from ztf_pipeutils.ztf_hdf5 import Read_LightCurve
import statistics as stat
import matplotlib.pylab as plt
import numpy as np
import operator


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
            c_RMS.append(np.std(md['c_err']))

        return z_moy, z_RMS, c_moy, c_RMS

    def plot_err_c_z(self, axhline=False, fontsize=15, error_bar=False, color='red'):
        z_moy, z_RMS, c_moy, c_RMS = self.moy_z_bin()
        plt.scatter(z_moy, c_moy)
        plt.xlabel('$z$', fontsize=fontsize)
        plt.ylabel('$\sigma_c$', fontsize=fontsize)
        if axhline:
            plt.axhline(y=0.04, color='r', linestyle='-')
        if error_bar:
            plt.errorbar(z_moy, c_moy, xerr=z_RMS, yerr=c_RMS,
                         fmt='none', capsize=2, elinewidth=1, color=color)
