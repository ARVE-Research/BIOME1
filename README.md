# BIOME1
This is an implementation of the BIOME1 global vegetation model described in Prentice et al. (1992), with some additional modifications and updates for calculating radiation including accounting for changes in Earth's orbit around the sun in the past and surface incident shortwave and longwave radiation. This Fortran 90 code was originally written around 2007-2009 by Jed O. Kaplan.

Prentice, I. C., Cramer, W., Harrison, S. P., Leemans, R., Monserud, R. A., & Solomon, A. M. (1992). A Global Biome Model Based on Plant Physiology and Dominance, Soil Properties and Climate. *J Biogeogr*, 19(2), 117-134. [doi:10.2307/2845499](https://doi.org/10.2307/2845499)

The following software needs to be available in the path when building:

- Autotools (including Autoconf, Automake, M4, and libtool)
- pkg-config
- a Fortran compiler
- HDF5 (only C library required)
- netCDF compiled against the above HDF5, including C and Fortran libraries
