import simsurvey_tools as sst
import pandas as pd
import simsurvey
import os


class Simul_lc:
    "Definition of a class that simule light curve"

    def __init__(self, folder_dir, sfd98File, rcidFile, csvFile, ztf_fields, z_range=(0.01, 0.1), dec_range=(-30, 90), n_det=1,
                 ntransient=11, seed=70, threshold=1, **kwargs):
        """
        Parameters
        ----------
        folder_dir : str
            Name of the folder directory to find sfd98 file.
        z_range : (int,int)
            redshift range (default=(0.01, 0.1)).
        dec_range : (int,int)
            declinaison range for observation (default=(-30, 90)).
        n_det : int
            required number of detections (default=1).
        n_transient : int
            we can change the number of transientor the rate (default=11).
        seed : int
            default=70
        threshold : int
            S/N requirement for detection (default=1).
        """

        self.sfd98_dir = os.path.join(folder_dir, sfd98File)
        self.rcid_dir = os.path.join(folder_dir, rcidFile)
        self.csv_dir = os.path.join(folder_dir, csvFile)
        self.ztf_fields_dir = os.path.join(folder_dir, ztf_fields)

        self.fields = sst.load_ztf_fields(filename=self.ztf_fields_dir)
        self.ccds = sst.load_ztf_ccds(
            filename=self.rcid_dir, num_segs=64)  # it's rcid

        self.obs = pd.read_csv(self.csv_dir)
        self.simul = self.simul_lc(
            z_range, dec_range, ntransient, seed, n_det, threshold, **kwargs)

    def __call__(self):
        """
        Return
        ------
        lc : LightcurveCollection
            Collection of simulated light curve, to cheak the firt lc : lc[0]
        """
        survey = self.simul
        lc = survey.get_lightcurves()

        print('light curves', lc.meta)
        return lc

    def simul_lc(self, z_range, dec_range, ntransient, seed, n_det, threshold, **kwargs):

        self.obs['field'] = self.obs['field'].astype('int64')
        self.obs['time'] = self.obs['time'] - 2400000.5

        plan = simsurvey.SurveyPlan(time=self.obs['time'], band=self.obs['band'], zp=self.obs['zp'], obs_field=self.obs['field'],
                                    obs_ccd=self.obs['rcid'], skynoise=self.obs['skynoise'],
                                    fields={k: v for k, v in self.fields.items(
                                    ) if k in ['ra', 'dec', 'field_id', 'width', 'height']},
                                    ccds=self.ccds)

        mjd_range = (plan.cadence['time'].min() - 30,
                     plan.cadence['time'].max() + 30)

        transientprop = {}
        transientprop['lcsimul_func'] = 'basic'
        transientprop['lcsimul_prop'] = {}
        transientprop['lcsimul_prop']['color_mean'] = kwargs.pop('color_mean')
        transientprop['lcsimul_prop']['color_sigma'] = kwargs.pop(
            'color_sigma')
        transientprop['lcsimul_prop']['stretch_mean'] = kwargs.pop(
            'stretch_mean')
        transientprop['lcsimul_prop']['stretch_sigma'] = kwargs.pop(
            'stretch_sigma')

        tr = simsurvey.get_transient_generator(zrange=z_range, transient='Ia', template='salt2',
                                               dec_range=dec_range, mjd_range=mjd_range, sfd98_dir=self.sfd98_dir,
                                               ntransient=ntransient, seed=seed, transientprop=transientprop)

        survey = simsurvey.SimulSurvey(
            generator=tr, plan=plan, n_det=n_det, threshold=threshold)
        return survey
