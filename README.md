# Real-time-GHF
A real-time extension of the generalized Hartree Fock (GHF) method, based on the algorithm presented in _Ding, F.; Goings, J. J.; Frisch, M. J.; Li, X. Ab initio non-relativistic spin dynamics. J. Chem. Phys. 2014, 141, 214111. DOI: 10.1063/1.4902884._

`rt_ghf.py` is the main RT-GHF class. `h_atom.py` and `li_atom.py` are example scripts, both of which are from Ding et al. `h_mag_05.png` and `li_mag.png` are the output graphs for the scripts. 

The hydrogen atom example exactly matches Ding et al. The lithium atom example agrees with the results from Ding et al., with the exception that the initial magnetization of the electron is along the _x_ direction rather than the _y_ direction. 

_Note:_ As of Jan 8, 2024, the file `rt_ghf.py` has been updated to print data to a file and to include a separate function to plot data. The example files were created using an older version of the code that appended observables to numpy arrays and that included plotting within the dynamics definition, and so do not include the `var.plot()` call. As of Jan 17, 2024, energy is also calculated as an observable and a separate plot function has been included to plot energy as a function of time.
