from ztf_pipeutils.ztf_hdf5 import Read_LightCurve
import statistics as stat
import matplotlib.pylab as plt
import numpy as np
import operator
import bokeh
import bokeh.plotting
from bokeh.plotting import figure
from bokeh.io import show, output_file
bokeh.plotting.output_notebook()


class Histo:
    "class which affect mask on meta data"

    def __init__(self, metaFitInput='Meta_fit.hdf5', inputDir='dataLC'):

        self.cl = Read_LightCurve(file_name=metaFitInput, inputDir=inputDir)
        self.metadata = self.cl.get_table(path='meta')
        self.metadata['chi2'] = self.metadata['chisq']/self.metadata['ndof']
        mk = self.metadata['sel'] != 0
        self.metadata = self.metadata[mk]

    def histo_plt(self, x='z', bins=30, range=(0, 0.2)):
        plt.hist(self.metadata[x], bins=bins, range=range)

    def histo_bokeh(self, x='z', bins=50, density=True, plot_width=800, plot_height=600, x_range=(0, 0.2)):
        hist, edges = np.histogram(
            self.metadata[x], density=density, bins=bins)

        data = bokeh.models.ColumnDataSource({'top': hist,
                                              'left': edges[:-1],
                                             'right': edges[1:]})

        p = figure(x_axis_label=x,
                   y_axis_label='count',
                   x_range=x_range,
                   plot_width=plot_width,
                   plot_height=plot_height,
                   background_fill_color='whitesmoke',
                   background_fill_alpha=0.8)

        p.quad(top='top', bottom=0, left='left',
               right='right', line_color="white", source=data)
        show(p)
