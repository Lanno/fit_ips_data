# fit_ips_data
At CASS I assisted with a project that would analyze the remotely measured data of solar wind to determine its speed. The data in question is taken from interplanetary scintillation (IPS) radio arrays around the world and is in the form of an amplitude over time, a timeseries. This data is checked for excessive noise or interference and fourier transformed into the frequency space, then it is fit to a model that will yield a set of parameters including the speed of the solar wind along a particular axis.

This effort is unique in that presently, the most successful means of measuring the speed of the solar wind requires three radio arrays. Cross correlation is used on the data from these arrays to identify patterns in the shadow of the solar wind that are produced when it passes in front of stellar radio sources. These shadow patterns travel across the surface of the Earth as the solar wind moves outward from the Sun, the arrays measure these patterns when they pass over them and compare the times of measurement and the distance between the arrays to determine the speed. These arrays are in Japan and are managed by STELab. The project that I was assisting attempted to use a model that requires only one radio array to obtain a measurement, such as the ones in India, Korea, and Mexico.
