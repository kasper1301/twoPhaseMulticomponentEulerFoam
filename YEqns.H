{
    label inertIndex1 = -1;
    volScalarField Yt1(0.0*Y1[0]);

    volScalarField rightHandSide(0.0*alpha1);

    forAll(Y1, i)
    {
        if (Y1[i].name() != inertSpecie1)
        {
            volScalarField& Yi = Y1[i];

            fvScalarMatrix YiEqn1
            (
                fvm::ddt(alpha1, Yi)
              + fvm::div(alphaPhi1, Yi)
              - fvm::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), Yi)
              - fvm::laplacian(alpha1*phase1.turbulence().nuEff(), Yi)
             ==
                rightHandSide
            );

            YiEqn1.relax();
            YiEqn1.solve();

            Yi.max(0.0);
            Yt1 += Yi;
        }
        else
        {
            inertIndex1 = i;
        }
    }

    Y1[inertIndex1] = scalar(1) - Yt1;
    Y1[inertIndex1].max(0.0);
    
    label inertIndex2 = -1;
    volScalarField Yt2(0.0*Y2[0]);

    forAll(Y2, i)
    {
        if (Y2[i].name() != inertSpecie2)
        {
            volScalarField& Yi = Y2[i];

            fvScalarMatrix YiEqn2
            (
                fvm::ddt(alpha2, Yi)
              + fvm::div(alphaPhi2, Yi)
              - fvm::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), Yi)
              - fvm::laplacian(alpha2*phase2.turbulence().nuEff(), Yi)
             ==
                rightHandSide
            );

            YiEqn2.relax();
            YiEqn2.solve();

            Yi.max(0.0);
            Yt2 += Yi;
        }
        else
        {
            inertIndex2 = i;
        }
    }

    Y2[inertIndex2] = scalar(1) - Yt2;
    Y2[inertIndex2].max(0.0);
}
