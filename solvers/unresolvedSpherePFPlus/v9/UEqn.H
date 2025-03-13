// Solve the Momentum equation
volScalarField alphaNuEff(coupling.alpha() * turbulence->nuEff());

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(coupling.alpha(), U)
  + fvm::div(phi, U)
  - fvc::div(alphaNuEff * dev(fvc::grad(U)().T()))
  + fvm::Sp(coupling.Sp()/rho, U)
 ==
  - coupling.Su()/rho  
  + fvModels.source(U)
);

fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvConstraints.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn ==  - coupling.alpha() * fvc::grad(p));

    fvConstraints.constrain(U);
}
