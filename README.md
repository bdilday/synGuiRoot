synGuiRoot
==========

a gui (using ROOT) to estimate the ionization temperature and chemical abundances for supernova spectra. 

the "syn" refers to synow and/or synapps. The model spectra are built from a collection of single element spectra produced by synow, and one of the applications is to get a good set of initial parameter guesses to pass to synapps.

the APPDATA directory contains spectra of SN 1991bg and SN 2004dt as examples.

USAGE:

* you must install root. you should be able to do this with a package manager. for full details consult the root website, http://root.cern.ch/drupal/

* cd to the synGuiRoot dir

* untar the APPDATA.tgz file

* start an interactive root session by typing root at the command line prompt

* compile the synGuiRoot code to a shared library. to do this type
  .x synGuiRoot.C++ 
  at the interactive-root prompt 

* this should launch the application. from here it is a gui. The ions are labeled by their atomic number followed by a 2-digit integer giving the ionization state. for example, 1401 is singly-ionized Si, SiII. 2601 is FeII, 0800 is OI, 0100 is HI, 0200 is neutral He, 0201 is HeII, etc...

* to include an ion in the fit, check the box, and to exclude an ion leave the box unchecked. 

* the slider bars control the abundance. 

* the ionization temperature (e_temp) can be changed to one of three discrete values. 

* the blackbody temperature of the photosphere (t_phot) can be changed to one of three discrete values. 

* the left-most upper slider bar controls the redshift

* the next slider over controls the normalization of the spectrum

* the next two text boxes set the wavelength range to consider in the fit.

* there are five buttons;

* GetChi2 returns the chi-squared for the current set of parameters.

* DoFit uses the current parameters as inputs to minuit and runs migrad to get the maximum likelihood abundance, given the ions that are enabled, the electron temperature, and the photospheric temperature.

* ClearAll resets all parameters

* LinearBBFit fixes the abundances and just fits the normalization constant using linear least squares.

* SaveFigure saves the figure to a file.

* The drop-down menu to the right of the SaveFigure button sets the taumin parameter of synow; a large value means consider only "strong" lines.

* to fit a different SN, you must edit the synGuiRoot.C file and change the hard-coded file name, here,
SetSNFile("./APPDATA/1991bg_19911214_3640_9282_00.dat");

  then recompile the code,
 .x synGuiRoot.C++
  
* to add new SNe, prepare an input ascii file with columns
wavelength-(in angstroms) flux uncertainty-on-flux

  if you don't have uncertainties, you can set this to -1 for every wavelength.

