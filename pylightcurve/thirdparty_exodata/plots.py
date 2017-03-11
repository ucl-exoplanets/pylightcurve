# """ This module contains some plotting functions and plot types for easy plot creation
# """
# import os
# import sys
# import math
# from collections import OrderedDict
# 
# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib import rcParams
# import seaborn
# import time
# 
# from . import astroquantities as aq
# from . import astroclasses as ac
# 
# seaborn.set_style('whitegrid')
# rcParams.update({'figure.autolayout': True})
# 
# 
# class _GlobalFigure(object):
#     """ sets up the figure and subfigure object with all the global parameters.
#     """
# 
#     def __init__(self, size='small'):
#         self.setup_fig(size)
# 
#     def setup_fig(self, size='small'):
#         self.set_size(size)
# 
#     def set_size(self, size):
#         """ choose a preset size for the plot
#         :param size: 'small' for documents or 'large' for presentations
#         """
#         if size == 'small':
#             self._set_size_small()
#         elif size == 'large':
#             self._set_size_large()
#         else:
#             raise ValueError('Size must be large or small')
# 
#     def _set_size_small(self):
#         self.fig = plt.figure(figsize=(5, 4))
#         self.ax = self.fig.add_subplot(1, 1, 1)
#         self.set_title_size(10)
#         self.set_axis_label_size(12)
#         self.set_axis_tick_label_size(12)
# 
#     def _set_size_large(self):
#         self.fig = plt.figure(figsize=(10, 7.5))
#         self.ax = self.fig.add_subplot(1, 1, 1)
#         self.set_title_size(20)
#         self.set_axis_label_size(20)
#         self.set_axis_tick_label_size(20)
# 
#     def set_global_font_size(self, fontsize):
#         ax = self.ax
#         for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
#             item.set_fontsize(fontsize)
#         plt.draw()
# 
#     def set_title_size(self, fontsize):
#         self.ax.title.set_fontsize(fontsize)
#         plt.draw()
# 
#     def set_axis_label_size(self, fontsize):
#         for axis in (self.ax.xaxis.label, self.ax.yaxis.label):
#             axis.set_fontsize(fontsize)
#         plt.draw()
# 
#     def set_axis_tick_label_size(self, fontsize):
#         for axis in self.ax.get_xticklabels() + self.ax.get_yticklabels():
#             axis.set_fontsize(fontsize)
#         plt.draw()
# 
#     # set_foregroundcolor / set_backgroundcolor from Jasonmc https://gist.github.com/jasonmc/1160951
#     def set_foregroundcolor(self, color):
#          '''For the specified axes, sets the color of the frame, major ticks,
#              tick labels, axis labels, title and legend
#          '''
# 
#          ax = self.ax
# 
#          for tl in ax.get_xticklines() + ax.get_yticklines():
#              tl.set_color(color)
#          for spine in ax.spines:
#              ax.spines[spine].set_edgecolor(color)
#          for tick in ax.xaxis.get_major_ticks():
#              tick.label1.set_color(color)
#          for tick in ax.yaxis.get_major_ticks():
#              tick.label1.set_color(color)
#          ax.axes.xaxis.label.set_color(color)
#          ax.axes.yaxis.label.set_color(color)
#          ax.axes.xaxis.get_offset_text().set_color(color)
#          ax.axes.yaxis.get_offset_text().set_color(color)
#          ax.axes.title.set_color(color)
#          lh = ax.get_legend()
#          if lh != None:
#              lh.get_title().set_color(color)
#              lh.legendPatch.set_edgecolor('none')
#              labels = lh.get_texts()
#              for lab in labels:
#                  lab.set_color(color)
#          for tl in ax.get_xticklabels():
#              tl.set_color(color)
#          for tl in ax.get_yticklabels():
#              tl.set_color(color)
#          plt.draw()
# 
#     def set_backgroundcolor(self, color):
#          '''Sets the background color of the current axes (and legend).
#              Use 'None' (with quotes) for transparent. To get transparent
#              background on saved figures, use:
#              pp.savefig("fig1.svg", transparent=True)
#          '''
#          ax = self.ax
#          ax.patch.set_facecolor(color)
#          lh = ax.get_legend()
#          if lh != None:
#              lh.legendPatch.set_facecolor(color)
# 
#          plt.draw()
# 
#     def set_y_axis_log(self, logscale=True):
#         if logscale:
#             self.ax.set_yscale('log')
#         else:
#             self.ax.set_yscale('linear')
# 
#     def set_x_axis_log(self, logscale=True):
#         if logscale:
#             self.ax.set_xscale('log')
#         else:
#             self.ax.set_xscale('linear')
# 
# 
# class _AstroObjectFigs(_GlobalFigure):
#     """ contains extra functions for dealing with input of astro objects
#     """
# 
#     def __init__(self, objectList, size='small'):
#         _GlobalFigure.__init__(self, size)
#         self.objectList = objectList  # list of planets, stars etc
#         self._objectType = self._getInputObjectTypes()  # are we dealing with planets, stars etc
# 
#     def _getInputObjectTypes(self):
#         # get the type of the first object
#         firstObject = self.objectList[0]
#         firstObjectType = type(firstObject)
# 
#         for astroObject in self.objectList:
#             if not firstObjectType == type(astroObject):
#                 raise TypeError('Input object list contains mixed types ({0} and {1})'.format(firstObjectType,
#                                                                                               type(astroObject)))
#         return firstObjectType
# 
#     def _getParLabelAndUnit(self, param):
#         """ checks param to see if it contains a parent link (ie star.) then returns the correct unit and label for the
#          job from the parDicts
#         :return:
#         """
# 
#         firstObject = self.objectList[0]
# 
#         if isinstance(firstObject, ac.Planet):
#             if 'star.' in param:
#                 return _starPars[param[5:]]  # cut off star. part
#             else:
#                 return _planetPars[param]
#         elif isinstance(firstObject, ac.Star):
#             return _starPars[param]
#         else:
#             raise TypeError('Only Planets and Star object are currently supported, you gave {0}'.format(type(firstObject)))
# 
#     # OECPy specific functions
#     def _gen_label(self, param, unit):
#         # TODO could be star or other in future
#         parValues = self._getParLabelAndUnit(param)
# 
#         if unit is None:
#             return parValues[0]
#         else:
#             unitsymbol = self._get_unit_symbol(unit)
#             return '{0} ({1})'.format(parValues[0], unitsymbol)
# 
#     def _get_unit_symbol(self, unit):
#         """ Every quantities object has a symbol, but i have added `latex_symbol` for some types. This first checks for
#         a latex symbol and then if it fails will just use symbol
#         :return:
#         """
#         try:
#             return '${0}$'.format(unit.latex_symbol)
#         except AttributeError:
#             return unit.symbol
# 
# 
# class _BaseDataPerClass(_AstroObjectFigs):  # hangover from ETLOS's multiple child plot classes
#     """ Base class for plots counting the results by a attribute. Child classes must modify
#     * _classVariables (self._allowedKeys, )
#     * _getSortKey (take the planet, turn it into a key)
#     """
# 
#     def __init__(self, astroObjectList, unit=None, size='small'):  # added unit here as class will break without it anyway
#         _AstroObjectFigs.__init__(self, astroObjectList, size)
# 
#         self._classVariables()  # add info from child classes
#         self.unit = unit
#         self.resultsByClass = self._processResults()
# 
#     def _set_size_small(self):
#         _AstroObjectFigs._set_size_small(self)
#         self.xticksize = 8
# 
#     def _set_size_large(self):
#         _AstroObjectFigs._set_size_large(self)
#         self.xticksize = 15
# 
#     def _classVariables(self):
#         """ Variables to be loaded in init by child classes
# 
#         must set
#         *self._allowedKeys = (tuple of keys)
#         """
#         self._allowedKeys = ('Class')
# 
#     def _getSortKey(self, planet):
#         """ Takes a planet and turns it into a key to be sorted by
#         :param planet:
#         :return:
#         """
# 
#         return 'Class'
# 
#     def _processResults(self):
#         """ Checks each result can meet SNR requirments, adds to count
#         :return:
#         """
# 
#         resultsByClass = self._genEmptyResults()
# 
#         for astroObject in self.objectList:
#             sortKey = self._getSortKey(astroObject)
#             resultsByClass[sortKey] += 1
# 
#         return resultsByClass
# 
#     def _getPlotData(self):
#         """ Turns the resultsByClass Dict into a list of bin groups skipping the uncertain group if empty
# 
#         return: (label list, ydata list)
#         :rtype: tuple(list(str), list(float))
#         """
#         resultsByClass = self.resultsByClass
# 
#         try:
#             if resultsByClass['Uncertain'] == 0:  # remove uncertain tag if present and = 0
#                 resultsByClass.pop('Uncertain', None)
#         except KeyError:
#             pass
# 
#         plotData = list(zip(*resultsByClass.items()))  # (labels, ydata)
# 
#         return plotData
# 
#     def plotBarChart(self, title='', xlabel=None, c='#3ea0e4', label_rotation=False):
#         ax = self.ax
# 
#         labels, ydata = self._getPlotData()
# 
#         numItems = float(len(ydata))
# 
#         ind = np.arange(numItems)  # the x locations for the groups
#         ind /= numItems  # between 0 and 1
# 
#         spacePerBar = 1./numItems
#         gapratio = 0.5  # gap to bar ratio, 0.5 is even
#         width = spacePerBar*gapratio
#         gap = spacePerBar - width
# 
#         ax.bar(ind, ydata, width, color=c)
#         ax.set_xticklabels(labels)
#         ax.set_xticks(ind+(width/2.))
# 
#         for axis in self.ax.get_xticklabels():
#             axis.set_fontsize(self.xticksize)
# 
#         if self.unit is None:  # NOTE this is hacked in so it only works with DataPerParameterClass
#             self.unit = self._getParLabelAndUnit(self._planetProperty)[1]  # use the default unit defined in this class
#         self._yaxis_unit = self.unit
# 
#         if xlabel is None:
#             plt.xlabel(self._gen_label(self._planetProperty, self.unit))
#         else:
#             plt.xlabel(xlabel)
# 
#         if label_rotation:
#             plt.xticks(rotation=label_rotation)
#         plt.ylabel('Number of Planets')  # TODO could be stars, binaries etc
#         plt.title(title)
#         plt.xlim([min(ind)-gap, max(ind)+(gap*2)])
#         plt.draw()
# 
#     def plotPieChart(self, title=None, cmap_name='Pastel2'):
#         if sys.hexversion >= 0x02700000:
#             self.fig.set_tight_layout(True)
# 
#         # Generate plot data
#         labels, ydata = self._getPlotData()
#         totalobjects = float(np.sum(ydata))
#         fracs = [ynum / totalobjects for ynum in ydata]
# 
#         explode = np.zeros(len(fracs))  # non zero makes the slices come out of the pie
# 
#         # plot pie chart
#         cmap = plt.cm.get_cmap(cmap_name)
#         colors = cmap(np.linspace(0., 0.9, len(fracs)))
#         plt.pie(fracs, explode=explode, labels=labels, autopct='%1.1f%%', shadow=False, startangle=90,
#                 colors=colors)
# 
#         # generate plot / axis labels
#         for axis in self.ax.get_xticklabels():
#             axis.set_fontsize(self.xticksize)
# 
#         if self.unit is None:  # NOTE this is hacked in so it only works with DataPerParameterClass
#             self.unit = self._getParLabelAndUnit(self._planetProperty)[1]  # use the default unit defined in this class
#         self._yaxis_unit = self.unit
# 
#         if title is None:
#             title = 'Planet {0} Bins'.format(self._gen_label(self._planetProperty, self.unit))  # NOTE not always planet
# 
#         plt.title(title)
#         plt.xlim(-1.5, 1.5)
# 
#     def saveAllBarChart(self, filepath, *args, **kwargs):
#         self.plotBarChart(*args, **kwargs)
#         plt.savefig(os.path.join(filepath))
# 
#     def _genEmptyResults(self):
#         """ Uses allowed keys to generate a empty dict to start counting from
#         :return:
#         """
# 
#         allowedKeys = self._allowedKeys
# 
#         keysDict = OrderedDict()  # Note: list comprehension take 0 then 2 then 1 then 3 etc for some reason. we want strict order
#         for k in allowedKeys:
#             keysDict[k] = 0
# 
# 
#         resultsByClass = keysDict
# 
#         return resultsByClass
# 
# 
# class DataPerParameterBin(_BaseDataPerClass):
#     """ Generates Data for planets per parameter bin"""
# 
#     def __init__(self, results, planetProperty, binLimits, unit=None, size='small'):
#         """
#         :param planetProperty: property of planet to bin. IE 'e' for eccentricity, 'star.magV' for magV
#         :param binLimits: list of bin limits (lower limit, upper, upper, maximum) (note you can have maximum +)
#         :param unit: unit to scale param to (see general plotter)
#         :return:
#         """
# 
#         self._binlimits = binLimits
#         self._planetProperty = planetProperty
# 
#         self._genKeysBins()  # Generate the bin keys/labels (must do before base class processes results)
#         _BaseDataPerClass.__init__(self, results, unit, size)
# 
#     def _getSortKey(self, planet):
#         """ Takes a planet and turns it into a key to be sorted by
#         :param planet:
#         :return:
#         """
# 
#         value = eval('planet.'+self._planetProperty)
# 
#         # TODO some sort of data validation, either before or using try except
# 
#         if self.unit is not None:
#             try:
#                 value = value.rescale(self.unit)
#             except AttributeError:  # either nan or unitless
#                 pass
# 
#         return _sortValueIntoGroup(self._allowedKeys[:-1], self._binlimits, value)
# 
#     def _classVariables(self):
#         pass  # Overload as we dont want it to set anything in this class
# 
#     def _genKeysBins(self):
#         """ Generates keys from bins, sets self._allowedKeys normally set in _classVariables
#         """
#         binlimits = self._binlimits
# 
#         allowedKeys = []
#         midbinlimits = binlimits
# 
#         if binlimits[0] == -float('inf'):
#             midbinlimits = binlimits[1:]  # remove the bottom limit
#             allowedKeys.append('<{0}'.format(midbinlimits[0]))
# 
#         if binlimits[-1] == float('inf'):
#             midbinlimits = midbinlimits[:-1]
# 
#         lastbin = midbinlimits[0]
# 
#         for binlimit in midbinlimits[1:]:
#             if lastbin == binlimit:
#                 allowedKeys.append('{0}'.format(binlimit))
#             else:
#                 allowedKeys.append('{0} to {1}'.format(lastbin, binlimit))
#             lastbin = binlimit
# 
#         if binlimits[-1] == float('inf'):
#             allowedKeys.append('{0}+'.format(binlimits[-2]))
# 
#         allowedKeys.append('Uncertain')
#         self._allowedKeys = allowedKeys
# 
# 
# class GeneralPlotter(_AstroObjectFigs):
#     """ This class should be able to create a plot with lots of options like the online visual plots. In future it
#     should be turned into a GUI
#     """
# 
#     def __init__(self, objectList, xaxis=None, yaxis=None, xunit=None, yunit=None, xaxislog=False, yaxislog=False, size='small'):
#         """
#         :param objectList: list of astro objects to use in plot ie planets, stars etc
#         :param xaxis: value to use on the xaxis, should be a variable or function of the objects in objectList. ie 'R'
#             for the radius variable and 'calcDensity()' for the calcDensity function
#         :param yaxis: value to use on the yaxis, should be a variable or function of the objects in objectList. ie 'R'
#             for the radius variable and 'calcDensity()' for the calcDensity function
# 
#         :type objectList: list, tuple
#         :type xaxis: str
#         :type yaxis: str
#         """
#         _AstroObjectFigs.__init__(self, objectList, size)
# 
#         # setup vars - to be replaced by themes
#         self.set_marker_color()
#         self.set_marker_size()
# 
#         # params that can be modified by user
#         self.marker = 'o'
# 
#         # set later
#         self.xlabel = None
#         self.ylabel = None
# 
#         # Handle given parameters
#         if xaxis:
#             self.set_xaxis(xaxis, xunit)
#         else:
#             self._xaxis = None
# 
#         if yaxis:
#             self.set_yaxis(yaxis, yunit)
#         else:
#             self._yaxis = None
# 
#         self.set_x_axis_log(xaxislog)
#         self.set_y_axis_log(yaxislog)
# 
#     def _set_size_large(self):
#         _AstroObjectFigs._set_size_large(self)
#         self.set_marker_size(60)
# 
#     def plot(self):
#         xaxis = [float(x) for x in self._xaxis]
#         yaxis = [float(y) for y in self._yaxis]
# 
#         assert(len(xaxis) == len(yaxis))
# 
#         plt.scatter(xaxis, yaxis, marker=self.marker, facecolor=self._marker_color, edgecolor=self._edge_color,
#                     s=self._marker_size)
# 
#         plt.xlabel(self.xlabel)
#         plt.ylabel(self.ylabel)
# 
#     def set_xaxis(self, param, unit=None, label=None):
#         """ Sets the value of use on the x axis
#         :param param: value to use on the xaxis, should be a variable or function of the objects in objectList. ie 'R'
#         for the radius variable and 'calcDensity()' for the calcDensity function
# 
#         :param unit: the unit to scale the values to, None will use the default
#         :type unit: quantities unit or None
# 
#         :param label: axis label to use, if None "Parameter (Unit)" is generated here and used
#         :type label: str
#         """
# 
#         if unit is None:
#             unit = self._getParLabelAndUnit(param)[1]  # use the default unit defined in this class
#         self._xaxis_unit = unit
# 
#         self._xaxis = self._set_axis(param, unit)
#         if label is None:
#             self.xlabel = self._gen_label(param, unit)
#         else:
#             self.xlabel = label
# 
#     def set_yaxis(self, param, unit=None, label=None):
#         """ Sets the value of use on the yaxis
#         :param param: value to use on the yaxis, should be a variable or function of the objects in objectList. ie 'R'
#         for the radius variable and 'calcDensity()' for the calcDensity function
# 
#         :param unit: the unit to scale the values to
#         :type unit: quantities unit or None
# 
#         :param label: axis label to use, if None "Parameter (Unit)" is generated here and used
#         :type label: str
#         """
#         if unit is None:
#             unit = self._getParLabelAndUnit(param)[1]  # use the default unit defined in this class
#         self._yaxis_unit = unit
# 
#         self._yaxis = self._set_axis(param, unit)
#         if label is None:
#             self.ylabel = self._gen_label(param, unit)
#         else:
#             self.ylabel = label
# 
#     def _set_axis(self, param, unit):
#         """ this should take a variable or a function and turn it into a list by evaluating on each planet
#         """
#         axisValues = []
#         for astroObject in self.objectList:
#             try:
#                 value = eval('astroObject.{0}'.format(param))
#             except ac.HierarchyError:  # ie trying to call planet.star and one planet is a lone ranger
#                 value = np.nan
# 
#             if unit is None:  # no unit to rescale (a aq.unitless quanitity would otherwise fail with ValueError)
#                 axisValues.append(value)
#             else:
#                 try:
#                     axisValues.append(value.rescale(unit))
#                 except AttributeError:  # either nan or unitless
#                     axisValues.append(value)
# 
#         return axisValues
# 
#     def set_marker_color(self, color='#3ea0e4', edgecolor='k'):
#         """ set the marker color used in the plot
#         :param color: matplotlib color (ie 'r', '#000000')
#         """
#         # TODO allow a colour set per another variable
#         self._marker_color = color
#         self._edge_color = edgecolor
# 
#     def set_marker_size(self, size=30):
#         """ set the marker size used in the plot
#         :param size: matplotlib size
#         """
#         self._marker_size = size
# 
# 
# class DiscoveryMethodByYear(object):
# 
#     def __init__(self, planet_list, methods_to_plot=('RV', 'transit', 'Other'), skip_solar_system_planets=True):
#         """ Produces a labelled stacked bar chart of planet discovery methods
# 
# 
#         :param planet_list: list of planet objects
#         :param methods_to_plot: list of discovery methods to plot (OEC syntax)
#         , use 'other' to group non specified methods
#         :param skip_solar_system_planets: don't include solar system planets
#         (including them will really stretch the xscale)
#         """
#         self.planet_list = planet_list
#         self.skip_solar_system_planets = skip_solar_system_planets
#         self.methods_to_plot = methods_to_plot
# 
#     def setup_keys(self):
#         """ Build the initial data dictionary to store the values
#         """
# 
#         discovery_methods = {}
#         discovery_years = {}
#         nan_list = []
# 
#         # Initial Loop to get keys
#         for planet in self.planet_list:
#             if 'Solar System' in planet.params['list'] and self.skip_solar_system_planets:
#                 continue
#             try:
#                 discovery_methods[planet.discoveryMethod] += 1
#             except KeyError:
#                 discovery_methods[planet.discoveryMethod] = 1
# 
#             try:
#                 discovery_years[planet.discoveryYear] += 1
#             except KeyError:
#                 discovery_years[planet.discoveryYear] = 1
# 
#             if planet.discoveryMethod is np.nan:
#                 nan_list.append(planet)
# 
#         self.nan_list = nan_list
# 
#         return discovery_years
# 
#     def generate_data(self):
# 
#         discovery_years = self.setup_keys()
# 
#         _year_list = np.array([year for year in discovery_years.keys() if year is not np.nan])
#         year_list = range(_year_list.min(), _year_list.max()+1)  # +1 needed as its < not <=
# 
#         discovery_methods = {k:0 for k in self.methods_to_plot}
#         plot_matrix = {k:discovery_methods.copy() for k in year_list}
# 
#         # TODO (ryan) the same sort of  loop as setup_keys again? really?
#         for planet in self.planet_list:
#             if 'Solar System' in planet.params['list'] and self.skip_solar_system_planets:
#                 continue
# 
#             if planet.discoveryMethod in self.methods_to_plot:
#                 discovery_method = planet.discoveryMethod
#             else:
#                 discovery_method = 'Other'
# 
#             try:
#                 plot_matrix[planet.discoveryYear][discovery_method] += 1
#             except KeyError:  # nan year
#                 pass
# 
#         return plot_matrix, year_list, discovery_years
# 
#     def plot(self, method_labels=None, exodata_credit=True):
#         """
#         :param method_labels: list of legend labels for the methods, i.e. give
#         'Radial Velocity' rather than 'RV'
# 
#         :return:
#         """
# 
#         if method_labels is None:
#             method_labels = self.methods_to_plot
# 
#         plot_matrix, year_list, discovery_years = self.generate_data()
# 
#         fig = plt.figure()
#         ax = fig.add_subplot(1, 1, 1)
# 
#         ind = np.arange(len(year_list))
#         colors = seaborn.color_palette("Set1", n_colors=8)
#         width = 0.9
#         bottom = np.zeros_like(year_list)
#         for i, method in enumerate(self.methods_to_plot):
#             bar_data = [plot_matrix[year][method] for year in year_list]
#             plt.bar(ind, bar_data, width, color=colors[i], bottom=bottom, label=method_labels[i], linewidth=0)
# 
#             bottom += np.array(bar_data) # plot the next bar stacked
# 
#         plt.ylabel('Number of planets', fontsize=14)
#         plt.xlabel('Year', fontsize=14)
#         plt.xticks(ind+width/2., year_list, rotation=70, fontsize=13)
#         plt.xlim(0, ind[-1]+1)
#         plt.legend(loc=2, fontsize=13)
# 
#         if exodata_credit is True:
#             ax.text(2, 420, 'Generated by ExoData on {}'.format(time.strftime("%d/%m/%Y")), fontsize=11)
# 
#         for i, year in enumerate(year_list):
#             try:
#                 total_planets = discovery_years[year]
#             except KeyError:  # year with no planets discovered
#                 continue  # Dont label
#             y_pos = bottom[i]+10
#             x_pos = i
# 
#             ax.text(x_pos, y_pos, '{}'.format(total_planets), fontsize=9)
# 
#         plt.grid(b=True, which='both', color='0.9', linestyle='-')
# 
#         return fig
# 
# 
# def _sortValueIntoGroup(groupKeys, groupLimits, value):
#     """ returns the Key of the group a value belongs to
#     :param groupKeys: a list/tuple of keys ie ['1-3', '3-5', '5-8', '8-10', '10+']
#     :param groupLimits: a list of the limits for the group [1,3,5,8,10,float('inf')] note the first value is an absolute
#     minimum and the last an absolute maximum. You can therefore use float('inf')
#     :param value:
#     :return:
#     """
# 
#     if not len(groupKeys) == len(groupLimits)-1:
#         raise ValueError('len(groupKeys) must equal len(grouplimits)-1 got \nkeys:{0} \nlimits:{1}'.format(groupKeys,
#                                                                                                          groupLimits))
# 
#     if math.isnan(value):
#         return 'Uncertain'
# 
#     # TODO add to other if bad value or outside limits
#     keyIndex = None
# 
#     if value == groupLimits[0]:  # if value is == minimum skip the comparison
#         keyIndex = 1
#     elif value == groupLimits[-1]:  # if value is == minimum skip the comparison
#         keyIndex = len(groupLimits)-1
#     else:
#         for i, limit in enumerate(groupLimits):
#             if value < limit:
#                 keyIndex = i
#                 break
# 
#     if keyIndex == 0:  # below the minimum
#         raise BelowLimitsError('Value {0} below limit {1}'.format(value, groupLimits[0]))
# 
#     if keyIndex is None:
#         raise AboveLimitsError('Value {0} above limit {1}'.format(value, groupLimits[-1]))
# 
#     return groupKeys[keyIndex-1]
# 
# 
# class OutOfLimitsError(Exception):
#     pass
# 
# 
# class BelowLimitsError(Exception):
#     pass
# 
# 
# class AboveLimitsError(Exception):
#     pass
# 
# # A parameter dict based on the current astroclass list and the parameter selected. Here the standard
# # units, full name for labels and short name as stored. All values that could be used are listed here and many are
# # commented because i dont think they should be used or will not work without further code to interpret
# # These dicts allow axis labels for the plots
# _planetPars = {
#     # paramKey    : (axis label,      unit)
#     #'isTransiting': ('Is Transiting', bool), TODO add bool and grouped values
#     'calcTransitDuration()': ('Transit Duration', aq.day),
#     'calcTransitDepth()': ('Transit Depth', None),
#     #'type()': ('Planet Class', None),
#     # 'massType()': ('Planet Mass Class', None),
#     # 'radiusType()': ('Planet Radius Class', None),
#     # 'temptype()': ('Planet Temp Class', None),
#     'mu': ('Mean Molecular Weight', None),
#     'albedo': ('Planet Albedo', None),
#     #'calcTemperature()': ('Mean Planet Temperature (Calculated)', aq.K),
#     #'estimateMass()': ('Planet Mass Estimated From Radius', aq.M_j),
#     # 'calcSMA()': ('Semi-Major Axis', aq.au),
#     # 'calcSMAfromT()': ('Semi-Major Axis (calculated from T)', aq.au),
#     # 'calcPeriod': ('Planet Radius', aq.R_j),
#     # 'discoveryMethod': ('Planet Radius', aq.R_j),
#     'e': ('Orbital Eccentricity', None),
#     'discoveryYear': ('Discovery Year', None),
#     # 'lastUpdate': ('Last Updated', Date),
#     'age': ('Planet Age', aq.Gyear),
#     # 'ra': ('RA', None),
#     # 'dec': ('DEC', None),
#     'd': ('Distance to System', aq.pc),
#     'R': ('Planet Radius', aq.R_j),
#     'T': ('Planet Temperature', aq.K),
#     'M': ('Planet Mass', aq.M_j),
#     'calcSurfaceGravity()': ('Planet Surface Gravity', aq.ms2),
#     'calcLogg()': ('Planet logg', None),
#     'calcDensity()': ('Planet Density', aq.gcm3),
#     'i': ('Orbital Inclination', aq.deg),
#     'P': ('Period', aq.day),
#     'a': ('Semi-Major Axis', aq.au),
#     # 'transittime': ('Transit Time', aq.JulianDay),
#     # TODO Catalogue Unit issue with theese
#     # 'periastron': ('Orbit Periastron', aq.deg),
#     # 'longitude': ('Orbit Longitude', aq.deg),
#     # 'ascendingnode': ('Orbit Ascending Node', aq.deg),
# }
# 
# _starPars = {
#     # paramKey    : (axis label,      unit)
#     'magV': ('V Magnitude', None),
#     'magB': ('B Magnitude', None),
#     'magH': ('H Magnitude', None),
#     'magI': ('I Magnitude', None),
#     'magJ': ('J Magnitude', None),
#     'magK': ('K Magnitude', None),
#     'd': ('Distance to System', aq.pc),
#     'calcLuminosity()': ('Stellar Luminosity', aq.L_s),
#     #'calcTemperature()': ('Stellar Temperature (Calculated)', aq.K),
#     'Z': ('Stellar Metallicity', None),
#     # 'Spectral Type': ('Spectral Type', str),
#     # 'estimateAbsoluteMagnitude': ('Stellar Age', aq.Gyear),
#     # 'estimateDistance': ('Stellar Age', aq.Gyear),
#     'age': ('Stellar Age', aq.Gyear),
#     # 'ra': ('RA', None),
#     # 'dec': ('DEC', None),
#     'R': ('Stellar Radius', aq.R_s),
#     'T': ('Stellar Temperature', aq.K),
#     'M': ('Stellar Mass', aq.M_s),
#     'calcSurfaceGravity()': ('Planet Surface Gravity', aq.ms2),
#     'calcLogg()': ('Stellar logg', None),
#     'calcDensity()': ('Stellar Density', aq.gcm3),
# }
# 
# # TODO Binary and System support for plots