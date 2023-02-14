// Solve the Momentum equation
phi = fvc::interpolate(U * coupling.alpha()) & mesh.Sf();

volScalarField alphaNuEff(coupling.alpha() * turbulence->nuEff());

MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    fvm::ddt(coupling.alpha(), U) 
  + fvm::div(phi, U)
  + MRF.DDt(coupling.alpha(), U)
  - fvc::div(alphaNuEff * dev(fvc::grad(U)().T()))
  - fvm::Sp(coupling.Sp()/rho, U)
 ==
    coupling.Su()/rho
  + fvModels.source(U)
);

UEqn.relax();

fvConstraints.constrain(UEqn);

if (piso.momentumPredictor())
{
    solve(UEqn == - coupling.alpha() * fvc::grad(p));

    fvConstraints.constrain(U);
}
