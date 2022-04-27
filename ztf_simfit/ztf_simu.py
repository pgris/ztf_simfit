import ztf_pipeutils.simsurvey_tools as sst
import pandas as pd
import simsurvey
import os


class Simul_lc:
    "Definition of a class that simule light curve"

    def __init__(self, ztfdataDir, dustmapDir, rcidFile, ztf_fields, z_range=(0.01, 0.1), n_det=1,
                 ntransient=11, seed=70, threshold=1, **kwargs):
        """
        Parameters
        ----------
        ztfdataDir : str
            directory with ztf data
        dustmapDir: str
          dir with dust maps
        rcidFile: str
          rcid files (ZTF internal)
        ztf_fields: str
           file of ztf_fields (ZTF internal)
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

        self.z_range = z_range
        self.ntransient = ntransient
        self.seed = seed
        self.n_det = n_det
        self.threshold = threshold
        self.kwargs = kwargs
        self.dustmapDir = dustmapDir
        # grab dustmaps if necessary
        if not os.path.isdir(self.dustmapDir):
            os.makedirs(self.dustmapDir)
            self.dustmaps()

        self.dustmapDir += '/sfd'

        rcids = os.path.join(ztfdataDir, rcidFile)
        #obsPath = os.path.join(obsDir, obsFile)
        ztf_fields = os.path.join(ztfdataDir, ztf_fields)

        self.fields = sst.load_ztf_fields(filename=ztf_fields)
        self.ccds = sst.load_ztf_ccds(
            filename=rcids, num_segs=64)  # it's rcid

    def __call__(self, obs, ra_range, dec_range):
        """
        Parameters
        --------------
        obs: array
          array of observations
        Return
        ------
        lc : LightcurveCollection
            Collection of simulated light curve, to cheak the firt lc : lc[0]
        """
        #survey = self.simul
        survey = self.simul_lc(obs,  ra_range, dec_range, self.z_range,
                               self.ntransient, self.seed, self.n_det, self.threshold, **self.kwargs)

        lc = survey.get_lightcurves()

        return lc

    def simul_lc(self, obs, ra_range, dec_range, z_range, ntransient, seed, n_det, threshold, **kwargs):

        obs['field'] = obs['field'].astype('int64')
        #obs['time'] = obs['time'] - 2400000.5

        plan = simsurvey.SurveyPlan(time=obs['time'], band=obs['band'], zp=obs['zp'], obs_field=obs['field'],
                                    obs_ccd=obs['rcid'], skynoise=obs['skynoise'],
                                    fields={k: v for k, v in self.fields.items(
                                    ) if k in ['ra', 'dec', 'field_id', 'width', 'height']},
                                    ccds=self.ccds)

        mjd_range = (plan.cadence['time'].min(),
                     plan.cadence['time'].max())

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
                                               ra_range=ra_range,
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
