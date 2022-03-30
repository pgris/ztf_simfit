import ztf_pipeutils.simsurvey_tools as sst
import pandas as pd
import simsurvey
import os


class Simul_lc:
    "Definition of a class that simule light curve"

    def __init__(self, ztfdataDir, logDir, dustmapDir, rcidFile, csvFile, ztf_fields, z_range=(0.01, 0.1), dec_range=(-30, 90), n_det=1,
                 ntransient=11, seed=70, threshold=1, **kwargs):
        """
        Parameters
        ----------
        ztfdataDir : str
            directory with ztf data
        logDir:
           wdir with cadence logs
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

        self.dustmapDir = dustmapDir
        # grab dustmaps if necessary
        if not os.path.isdir(self.dustmapDir):
            os.makedirs(self.dustmapDir)
            self.dustmaps()

        self.dustmapDir += '/sfd'

        rcids = os.path.join(ztfdataDir, rcidFile)
        logObs = os.path.join(logDir, csvFile)
        ztf_fields = os.path.join(ztfdataDir, ztf_fields)

        self.fields = sst.load_ztf_fields(filename=ztf_fields)
        self.ccds = sst.load_ztf_ccds(
            filename=rcids, num_segs=64)  # it's rcid

        self.obs = pd.read_csv(logObs)

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
                                               dec_range=dec_range, mjd_range=mjd_range, sfd98_dir=self.dustmapDir,
                                               ntransient=ntransient, seed=seed, transientprop=transientprop)

        survey = simsurvey.SimulSurvey(
            generator=tr, plan=plan, n_det=n_det, threshold=threshold)
        return survey

    def dustmaps(self):
        """
        method to grab dustmaps
        Dust maps will be placed in self.dustmapDir/sfd

        """
        from dustmaps.config import config
        config['data_dir'] = self.dustmapDir
        import dustmaps.sfd
        dustmaps.sfd.fetch()
