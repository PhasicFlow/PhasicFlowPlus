// Solve the Momentum equation
volScalarField alphaNuEff(coupling.alpha() * turbulence->nuEff());

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(coupling.alpha(), U)
  + fvm::div(phi, U)
  //- fvc::div(alphaNuEff * dev(fvc::grad(U)().T()))
  + turbulence->divDevReff(U)
  + fvm::Sp(coupling.Sp()/rho, U)
 ==
  - coupling.Su()/rho  
  + fvOptions(U)
);

fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn ==  - coupling.alpha() * fvc::grad(p));

    fvOptions.correct(U);
}
