rake
=====

Handling rake in kinematic inversions
-------------------------------------

PYEQ core algorithm handles rake as constant through time and provided at every subfault.
However, a few options allows to search for the 'best' rake as well.

Handling time-constant rake in PYEQ is managed by two options:

``rake_type`` specifies how rake is calculated at subfaults

``rake_value`` specifies some values to calculate the rake at every subfault

Whatever ``rake_type``, the

rake_type='fixed'
-----------------

rake is a constant value provided by ``rake_value`` at all subfaults and all model time step.
rake_value is a float in the range [-180.,180.]

rake_type='euler'
-----------------

rake is calculated at every subfault from an euler pole provided as a string ``lon/lat/w/style``.
lon, lat are euler pole coordinates in decimal degrees.
w is the angular velocity in decimal degrees per Myr
style is among ['inverse','normal','left-lateral','right-lateral'] is used used to impose the non-negativity constraint.

rake_delta
----------

If rake_delta is provided (default is None), itallows some search for the best rake.

rake_delta = '1/-20/20/3/4' will search for a single value correcting the rake provided by the rake_type and rake_value.
The search will be performed over the [-20,+20] degrees, with a prior constraint of 3 degrees and a weight of 4 for
spatially smoothing rake change.

rake_delta = '2/-20/20/3/4' is similar as above except that the rake corrected is estimated at each subfault.




The `pandas
<https://https://pandas.pydata.org/>`_ library offers a versatile way to build dates arrays, through `pandas date_range documentation
<https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.date_range.html/>`_ method.
pandas.date_range command are recognized in PYEQ. Here are a few examples.

Most of the time, you want to model daily time series with daily estimates of slip.
To generate dates at 12:00 from April 17, 2016 to April 21, every day:

::

	-dates "pandas.date_range(start='2016/4/17 12:00:00', end='2016/4/21 12:00:00', freq='1D')"

You can check the results, running in a jupyter notebook or ipython:

::

	import pandas
	pandas.date_range(start='2016/4/17 12:00:00', end='2016/4/21 12:00:00', freq='1D')

To generate dates at 12:00 from April 17, 2016 to April 28, by step of 2 days:

::

	import pandas
	pandas.date_range(start='2016/4/17 12:00:00', end='2016/4/28 12:00:00', freq='2D')

freq can be 'W','M','Y' for week, month year. Add closed='left' (closed='right') to control inclusion of start and end dates respectively. See `pandas date_range documentation
<https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.date_range.html/>`_ for many other options.

To generate dates every 3 days at 12:00 from April 17 with 5 dates:

::

	import pandas
	pandas.date_range(start='2016/4/17 12:00:00', periods=5, freq='3D')

To generate 11 dates (that is model will have 10 time steps) between 2016/01/01 and 2016/08/01:

::

	import pandas
	pandas.date_range(start='2016/4/17 12:00:00', end='2016/8/1 12:00:00' , periods=11)

Note: in this case, dates are not anymore expressed at 12:00. In PYEQ the rounding option will handle this. Setting -rounding 'day' will round all observation and model dates at 12:00. -rounding 'minute' will preserve the original time but also interpret observation dates at their original minute.

Dates can ve comnbined using 'union'

::


	dates = pandas.date_range(start='2016/4/17 12:00:00', end='2016/5/17 12:00:00', freq='1D').union(pandas.date_range(start='2016/5/18 12:00:00', end='2019/1/30 12:00:00', freq='30D'))


-dates as a numpy npy file
---------------------------------

A numpy npy file can be provided. It is assumed to include a numpy 1D array of datetime objects.

-dates as a text file
----------------------------

The format assummed for a text file is:

::

	2016-04-17 12:00:00
	2016-04-18 12:00:00
	2016-04-19 12:00:00
	2016-04-20 12:00:00
	2016-04-21 12:00:00
	2016-04-22 12:00:00
	2016-04-23 12:00:00
	2016-04-24 12:00:00

-dates all
----------------------------

With this option, PYEQ will use all dates read from the input time series.

