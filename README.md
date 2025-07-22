In the examples folder you find `simulate_collagen.jl` with a few simulations of the SHG intensity and phase of the amplitude under various illumination and polarized detection conditions.
The coefficients were taken from 
P. Stoller, K.M. Reiser, P.M. Celliers, and A.M. Rubenchik* “Polarization-Modulated Second Harmonic Generation in Collagen”, Biophys. J. 82, 3330–3342 (2002).
with gamma = - 0.7.

The 3rd order tensor is obtained by calling `get_chi2_stoller(adim)`, where `adim` corresponds to the dimension the collagen 1 fiber is oriented towards.

The anglular dependence is obtained by variing the polarization angles of illumination and assuming a corresponding detection.
This can be done with the function `get_intensity` and `get_angle`.
