// Solve the Momentum equation
phi = fvc::interpolate(U * coupling.alpha()) & mesh.Sf();

volScalarField alphaNuEff(coupling.alpha() * turbulence->nuEff());

fvVectorMatrix UEqn
(
    fvm::ddt(coupling.alpha(), U) 
  + fvm::div(phi, U)
  - fvc::div(alphaNuEff * dev(fvc::grad(U)().T()))
  + fvm::Sp(coupling.Sp()/rho, U)
 ==
  - coupling.Su()/rho
  + fvOptions(U)
);

UEqn.relax();

fvOptions.constrain(UEqn);

if (piso.momentumPredictor())
{
    solve(UEqn == - coupling.alpha() * fvc::grad(p));

    fvOptions.correct(U);
}
