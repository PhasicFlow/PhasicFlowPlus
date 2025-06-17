Foam::volSymmTensorField devRhoReff = -rho*turbulence->nu()*dev(twoSymm(fvc::grad(U)));

// Calculating the fluid to solid interaction
coupling.calculateFluidInteraction
(
	p*rho, // Correcting kinematic pressure
	devRhoReff,
	particleSTLs
);

//Sending fluid forces and torques to DEM
coupling.sendFluidForceToDEM();
coupling.sendFluidTorqueToDEM();

