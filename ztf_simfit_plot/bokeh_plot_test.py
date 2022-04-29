from bokeh.layouts import column
from bokeh.plotting import Figure, output_file, show, save, output_notebook, curdoc
from bokeh.models import MathML
from bokeh.util.compiler import TypeScript
from bokeh.transform import factor_cmap, factor_mark
from bokeh.sampledata.penguins import data
from bokeh.plotting import figure, show
from ztf_pipeutils.ztf_hdf5 import Read_LightCurve

import bokeh
import bokeh.plotting
bokeh.plotting.output_notebook()


class Bokeh_plot:
    "class which plot with bokeh"

    def __init__(self, file_name='Meta.hdf5', inputDir='dataLC'):

        self.class_SN = Read_LightCurve(file_name=file_name, inputDir=inputDir)
        self.meta_SN = self.class_SN.get_table(path='meta')
        self.meta_SN['chi'] = self.meta_SN['chisq']/self.meta_SN['ndof']

    def plot1(self, xlim=[0], ylim=[0.04], tooltips=[('SN path', '@path'), ('SN chi', '@chi')], y_range=(0, 0.1),
              plot_width=800, plot_height=600):

        selec = self.meta_SN['sel'] == 1
        msk = self.meta_SN['c_err'] != -1.0
        M = selec & msk

        data = bokeh.models.ColumnDataSource({
            'redshift':    self.meta_SN["z"][M],
            'err_c':       self.meta_SN["c_err"][M],
            'path':        self.meta_SN["path"][M],
            'chi':         self.meta_SN['chi'][M],
        })

        fig = bokeh.plotting.figure(
            x_axis_label='$$z$$',
            y_axis_label='$$\sigma_{c}$$',
            y_range=y_range,
            plot_width=plot_width,
            plot_height=plot_height,
            background_fill_color='whitesmoke',
            background_fill_alpha=0.8
        )

        fig.add_tools(bokeh.models.HoverTool(tooltips=tooltips, mode='mouse',))
        fig.ray(x=xlim, y=ylim, length=300, angle=0,
                color='red', legend_label='0.04')
        fig.legend.location = "top_left"
        fig.legend.click_policy = "hide"

        fig.scatter(x='redshift', y='err_c',
                    source=data, fill_alpha=0.2, size=10)

        bokeh.plotting.show(fig)

    def plot2(self):

        selec = self.meta_SN['sel'] == 1
        msk = self.meta_SN['c_err'] != -1.0

        MASK = self.meta_SN['n_i_band'] == 0
        MASK2 = self.meta_SN['n_i_band'] != 0
        M = selec & msk & MASK
        M2 = selec & msk & MASK2
        M3 = selec & msk

        m1 = self.meta_SN.copy()[M]
        m2 = self.meta_SN.copy()[M2]

        data = bokeh.models.ColumnDataSource({
            'redshift_rg':    m1["z"],
            'err_c_rg':       m1["c_err"],
            'path_rg':        m1["path"], })
        data2 = bokeh.models.ColumnDataSource({
            'redshift_rgi':   m2["z"],
            'err_c_rgi':      m2["c_err"],
            'path_rgi':       m2["path"], })

        fig = bokeh.plotting.figure(
            x_axis_label='$$z$$',
            y_axis_label='$$\sigma_{c}$$',
            y_range=(0, 0.1),
            plot_width=800,
            plot_height=600,
            background_fill_color='whitesmoke',
            background_fill_alpha=0.8
        )

        fig.scatter(x='redshift_rg', y='err_c_rg', source=data, fill_alpha=0.2,
                    size=10, legend_label='r,g band',     color='LightSeaGreen')
        fig.scatter(x='redshift_rgi', y='err_c_rgi', source=data2,
                    fill_alpha=0.2, size=10, legend_label='r,g,i band', color='FireBrick')

        fig.add_tools(bokeh.models.HoverTool(tooltips=[
                      ('SN (rg) path ', '@path_rg'), ('SN (rgi) path ', '@path_rgi')], mode='mouse',))

        fig.ray(x=[0], y=[0.04], length=300, angle=0,
                color='red', legend_label='0.04')
        fig.legend.location = "top_left"
        fig.legend.click_policy = "hide"

        bokeh.plotting.show(fig)
