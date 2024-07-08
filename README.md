
#Here is macros for alignment of pixel-local-coordinates on the timing counter structure.

##Here we say:

- CAD coordinate : designed values (x,y,z) for pixels and reference points where the COBRA center is set as (0,0,0).
- FARO coordinate : measured points data for each pixel by FARO ScanArm.
- COBRA 20XX coordinate : measured values (x,y,z) of reference points' center for 20XX installation into COBRA magnet structure.

###1. Align FARO data and designed values in COBRA.
- make it easy to get the initial value of position and rotation angle of pixels.
`src/FARO_CAD_calib.cpp`

###2. Separate points into 512 groups. Each of them belongs to one pixel.
`Autokiridashi_FARO_CAD.cpp`

###3. Fit the pixels by cuboid (L x W x T = 120 x 40 (or 50) x 5 mm^3) or only by the top surface (L x W = 120 x 5 mm^3)
- parameters (dx, dy, dz, dphi, dtheta, dpsi).
- former one is essentially correct.
- latter one gives instant and not bad estimation.
`fit.cpp` (`optfit.cpp` is being developed) or `Initial_value.cpp`

###4. Combine the fit results above and reference points in the 20XX operation.
- treated as just a parallel shift: (dx2, dy2, dz2)
`Ana_VS_COBRA.cpp`

