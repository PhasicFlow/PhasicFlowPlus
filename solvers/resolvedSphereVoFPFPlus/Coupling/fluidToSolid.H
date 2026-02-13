Foam::volSymmTensorField devRhoReff = -rho*mixture.nu()*dev(twoSymm(fvc::grad(U)));

// Calculating the fluid to solid interaction
coupling.calculateFluidInteraction
(
	p,
	devRhoReff,
	particleSTLs
);

//Sending fluid forces and torques to DEM
coupling.sendFluidForceToDEM();
coupling.sendFluidTorqueToDEM();

