{
    rho = thermo.rho();

    // 1st part of the density Correction
   	thermo.rho() -= psi*p;

    volScalarField rAU(1.0/UEqn.A());
    U = rAU*UEqn.H();

    if (pimple.transonic())
    {
        surfaceScalarField phiv
        (
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rAU, rho, U, phi)
        );

        phi = fvc::interpolate(rho)*phiv;

        surfaceScalarField phid
        (
            "phid",
            fvc::interpolate(thermo.psi())*phiv
        );

        fvScalarMatrix pDDtEqn
        (
            fvc::ddt(rho) + fvc::div(phi)
          + correction(fvm::ddt(psi, p) + fvm::div(phid, p))
        );

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                pDDtEqn
              - fvm::laplacian(rho*rAU, p)
            );

            pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

            if (pimple.finalNonOrthogonalIter())
            {
                phi += pEqn.flux();
            }
        }
    }
    else
    {
        phi =
            fvc::interpolate(rho)
           *(
                (fvc::interpolate(U) & mesh.Sf())
              + fvc::ddtPhiCorr(rAU, rho, U, phi)
            );

        fvScalarMatrix pDDtEqn
        (
            fvc::ddt(rho) + psi*correction(fvm::ddt(p))
          + fvc::div(phi)
        );

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                pDDtEqn
              - fvm::laplacian(rho*rAU, p)
            );

            pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

            if (pimple.finalNonOrthogonalIter())
            {
                phi += pEqn.flux();
            }
        }
    }

    // 2nd part of the density Correction
  	thermo.rho() += psi*p;

    #include "rhoEqn.H"
    #include "compressibleContinuityErrs.H"

    U -= rAU*fvc::grad(p);
    U.correctBoundaryConditions();
    K = 0.5*magSqr(U);

    dpdt = fvc::ddt(p);
}
