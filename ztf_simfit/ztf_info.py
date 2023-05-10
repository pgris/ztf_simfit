import operator
from astropy.table import Table, hstack, vstack
from ztf_pipeutils.ztf_hdf5 import Read_LightCurve
import numpy as np


def get_info(data, info, name_info='name', col_info='col', col_sel='colsel', type_sel='typesel', thresh_info_a='thresha', op_info_a='opa', thresh_info_b='threshb', op_info_b='opb'):
    """"
    Function to estimate infos from data

    Parameters
    --------------
    data: astropy table
      data to process (lc, ...)
    info: astropy table
      values to estimate
    name_info: str, opt
      name of the selection column in info (default: name)
    col_info: str, opt
       name of the col column in info (default: col)
    thresh_info: str, opt
       name of the lim column in info (default: lim)
    type_info: str, opt
        name of the type column in info (default: type)
    op_info: str, opt
         name of the op column in info (default: op)

    """

    res = Table()
    for row in info:

        col = row[col_sel]
        type_ = eval(row[type_sel])
        opa = eval(row[op_info_a])
        thresha = type_(row[thresh_info_a])
        opb = eval(row[op_info_b])
        threshb = type_(row[thresh_info_b])
        maska = opa(data[col], thresha)
        maskb = opb(data[col], threshb)
        new_tab = data[maska & maskb]

        nn = len(new_tab[row[col_info]])
        ccol = row[name_info]
        if 'epochs' in ccol or 'bands' in ccol:
            nn = len(np.unique(new_tab[row[col_info]]))
        res[row[name_info]] = [nn]

    return res


def get_selec(tab, selec_tab, name_selec='name', thresh_selec='thresh', op_selec='op', type_selec='type', sel_tag_col='sel'):
    """
    fonction to set a select flag (column) to tab data

    Parameters
    ---------------
    tab: astropy table
      data to process
    selec_tab: astropy table
       tab of selection criteria
    name_selec: str, opt
      name of the name column in selec_tab (default: name)
    thresh_selec: str, opt
      name of the thresh column in selec_tab (default: thresh)
     op_selec: str, opt
      name of the op column in selec_tab (default: op)
    type_selec: str, opt
      name of the type column in selec_tab (default: type)
    sel_tag_col: str, opt
      name of the column for the selection in the tab (default: sel)
    """

    idx = True
    tab[sel_tag_col] = 0
    for row in selec_tab:
        col = row[name_selec]
        op = eval(row[op_selec])
        type_ = eval(row[type_selec])
        lim = type_(row[thresh_selec])
        idx &= op(tab[col], lim)

    tab[sel_tag_col][idx] = 1

    return tab


def complete_lc(lc, snr):

    lc['SNR'] = lc['flux'] / lc['fluxerr']
    lc['phase'] = (lc['time'] - lc.meta['t0']) / (1-lc.meta['z'])
    idx = lc['SNR'] >= snr

    return lc[idx]


class Info:
    """
    class to estimate infos on data

    Parameters
    --------------
    metaName: str
      name of the meta file to process
    metaDir: str
       dir where metaFile is located
    info_tab: astropy table
       info definition
    snr: float, opt
      snr LC points min val (default: 1.)
    """

    def __init__(self, metaName, metaDir, info_tab, snr=1.0):

        # read metadata
        read_meta = Read_LightCurve(file_name=metaName, inputDir=metaDir)
        self.metadata = read_meta.get_table(path='meta')

        # get lc
        lcDir = self.metadata.meta['directory']
        lcName = self.metadata.meta['file_name']

        self.read_lc = Read_LightCurve(file_name=lcName, inputDir=lcDir)

        # infos
        self.info_tab = info_tab

        # snr
        self.snr = snr

    def __call__(self, bad_prefix='bad'):
        """
        Main method for processing

        Returns
        ----------
        astropy table with the list of info for data

        """
        restab = Table()
        restab.meta = self.metadata.meta

        for meta in self.metadata:
            tt = Table(meta)
            path = meta['path']
            if bad_prefix not in path:
                lc = self.read_lc.get_table(path)
                lc = complete_lc(lc, self.snr)
                lc['mjd_min'] = np.min(lc['time'])
                lc['night'] = lc['time']-lc['mjd_min']+1
                lc['night'] = lc['night'].astype(int)
                lc.sort('night')
                res = get_info(lc, self.info_tab)
                tt = hstack([tt, res])
            else:
                names = self.info_tab['name'].tolist()
                vals = [[-1]]*len(self.info_tab)
                tb = Table(vals, names=names)
                tt = hstack([tt, tb])
            restab = vstack([restab, tt])

        return restab
