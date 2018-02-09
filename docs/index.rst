****************************************
Core Cosmology Library
****************************************

The Core Cosmology Library (CCL) is a standardized library of routines to calculate basic observables used in cosmology. It will be the standard analysis package used by the LSST Dark Energy Science Collaboration (DESC).

The core functions of this package include:

 * Matter power spectra :math:`P(k)` from numerous models including CLASS, the Mira-Titan Emulator  and halofit
 * Hubble constant :math:`H(z)` as well as comoving distances :math:`\chi(z)` and distance moduli :math:`\mu(z)`
 * Growth of structure :math:`D(z)` and :math:`f`
 * Correlation functions :math:`C_\ell` for arbitrary combinations of tracers including galaxies, shear and number counts
 * Halo mass function :math:`dn/dM` and halo bias :math:`b(M)`
 * Approximate baryonic modifications to the matter power spectra :math:`\Delta^2_{\rm baryons}`
 * Simple modified gravity extensions :math:`\Delta f(z)`

The source code is available on github at https://github.com/LSSTDESC/CCL.

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   source/installation
   source/cosmology_object
   source/reporting_bugs

.. toctree::
   :maxdepth: 1
   :caption: For Developers

   source/navigating_the_code
	     
.. toctree::
   :maxdepth: 1
   :caption: Reference

   CCL Reference/API<api/modules>
