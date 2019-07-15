# Paraffin #

A collection of pre and post processing tools for Puffin.

Includes packages for matching the beam to a periodic lattice, scaling particle data from SU to Puffin notation (and vice-versa), and up-sampling the sparse beam from accelerator codes to a finer distribution more suitable for a FEL simulation.

The supplied Python script `fullConv.py` utilizes the main components in this package as an example. It takes a particle file in [SU format](https://github.com/UKFELs/FXFEL), matches the beam to a FEL undulator line employing a FODO focusing channel, upsamples the beam for use in a FEL code, and then scales the data for input into [Puffin](https://github.com/UKFELs/Puffin).

## To Install

Do

```
pip install .
```

in this directory to install as a Python package. The Python package `paraffin` will then be available for use in your own Python scripts or functions (see the supplied exemplar `fullConv.py` script for an example of use).

## Beam matching and focusing lattice routines:

### Top-level routine:

The routine `SU2Matched` will read in a file, in SU format, and produce a second file, in SU format, which is matched to the periodic lattice specified by the user.

The periodic lattice is assumed to be either a FODO lattice or a single undulator. The FODO lattice is assumed to begin with a half strength (that is, a half length) quad, focusing in x and defocusing in y, followed by a free space drift, then an undulator module, then a free space drift, then a quad (defocusing in x), then a drift, then an undulator module, then a drift, then a final half-strength quad (focusing in x). The lattice is specified by passing the subroutine an instance of the undulator module class, the quad focusing factor (the magnitude of which is assumed equal in x and y, see below), and the separation distance between undulator modules (NOTE: **not** the distance from undulator to quad!!).

If a quad focusing factor of zero is passed to the function, then the beam will be matched to the undulator only. The undulator may have a natural focusing channel and an optional **strong** intra-undulator focusing channel. The beam will be matched to take both into account.

For example, the routine can be called with:
```
fnameout = SU2Matched(fnamein, puffVars, undmod, qf, DL)
```
where `fnamein` is the filename of the particle file to be matched, in SU format, and `fnameout` is the name of the output file which is written by the function, containing the matched beam, in SU format.

### Other utilities:

These other routines are at a deeper part of the API, and are utilized by the top-level `SU2Matched`. They are exposed to the user directly for convenience, as they may be useful for more complex operations.

```
beta, alpha = getMatchTwiss(fullMat)
```
This returns the matched Twiss parameters for a given transform matrix. The matrix is '1D', so only for one transverse direction. The matrix should be a 2x2 matrix. The output Twiss parameters are labelled meaningfully, corresponding to the commonly used beta and alpha Twiss parameters.

```
twissx, twissy = getFODOTwiss(puffVars, undmod, f, DL, emitx, emity)
```
This returns the Twiss parameters required for matching to a FODO lattice. Output is in the form of 2 3-element arrays. Each array contains the emittance (geometric), beta, and alpha required for matching to the lattice. `twissx` contains the matched Twiss in the x direction, `twissy` the same for the y direction. Input emittance is the geometric emittance.

```
twissx, twissy = getUndTwiss(puffVars, undmod, emitx, emity)
```
Same as `getFODOTwiss`, but for a single undulator rather than a periodic lattice.


```
x2, px2, y2, py2 = matchTwiss(x, px, y, py, TX1, TX2, TY1, TY2)
```
This performs a transform on the given beam transverse coordinates, to match *from* Twiss parameters `TX1` and `TY1` to `TX2` and `TY2`. In this case, Twiss parameters take the form of a 2-element list or array, containing beta, and alpha. Particle coordinates are input as arrays, where `x` and `y` are the particle coordinates in meters, and `px` and `py` are the transverse phase space coordinates `dx/dz` and `dy/dz`, respectively.

```
twxt, twyt = getBeamTwiss(MPs)
```
This calculates and returns the Twiss parameters, in x and y respectively, of the beam in the input array MPs.

### Setting up the undulator beam line elements

The undulator is setup with reference to the Puffin scaling class, and is defined using the number of undulator periods, the type of undulator (a string, which is blank, indicating the Puffin variably polarized undulator, or one of `helical`, `planepole` or `curved`), the polarization (given by ux and uy, which is ignored if any but the Puffin variably polarized undulator is used), and the in-undulator *strong* focusing wavenumber in x and y. A simple setup would be:

```
undmod = undulator(puffVars, undtype = '', Nw = 26, ux = 1., uy = 0.)
```

The focusing factor of the quads, given by F = (rho B) / (g L), where (rho B) is the magnetic rigidity, g is the gradient of the quad field (i.e. dBy/dx) and L is the length of the quad (in meters). (thin quad approximation)

The drift is defined in meters.

We gratefully acknowledge the support of STFC's ASTeC department for HPC access, using the STFC Hartree Centre,
and the John von Neumann Institute for Computing (NIC) on JUROPA at Julich Supercomputing Centre (JSC), under project HHH20
