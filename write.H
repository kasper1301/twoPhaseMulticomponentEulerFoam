    if (runTime.outputTime())
    {
        volVectorField Ur
        (
            IOobject
            (
                "Ur",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            U1 - U2
        );
        surfaceScalarField massFlux1
        (
            IOobject
            (
                IOobject::groupName("massFlux", phase1.name()),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(rho1)*phi1
        );
        surfaceScalarField massFlux2
        (
            IOobject
            (
                IOobject::groupName("massFlux", phase2.name()),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(rho2)*phi2
        );
        
        runTime.write();
    }
