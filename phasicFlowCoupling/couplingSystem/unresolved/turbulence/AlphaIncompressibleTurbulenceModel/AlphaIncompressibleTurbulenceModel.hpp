/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
      O        C enter of
     O O       E ngineering and
    O   O      M ultiscale modeling of
   OOOOOOO     F luid flow       
------------------------------------------------------------------------------
License
    This file is partially part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This file is partially part of phasicFlow. 

    It is a free software for simulating 
    granular and multiphase flows. You can redistribute it and/or modify it under
    the terms of GNU General Public License v3 or any other later versions. 
 
    phasicFlow is distributed to help others in their research in the field of 
    granular and multiphase flows, but WITHOUT ANY WARRANTY; without even the
    implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/


#ifndef __IncompressibleTurbulenceModel_hpp__
#define __IncompressibleTurbulenceModel_hpp__

#include "TurbulenceModel.H"
#include "incompressibleTurbulenceModel.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class IncompressibleTurbulenceModel Declaration
\*---------------------------------------------------------------------------*/

template<class TransportModel>
class AlphaIncompressibleTurbulenceModel
:
    public TurbulenceModel
    <
        volScalarField,
        geometricOneField,
        incompressibleTurbulenceModel,
        TransportModel
    >
{

public:

    typedef volScalarField      alphaField;
    typedef geometricOneField   rhoField;
    typedef TransportModel      transportModel;


    // Constructors

        //- Construct
        AlphaIncompressibleTurbulenceModel
        (
            const word& type,
            const volScalarField& alpha,
            const geometricOneField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const TransportModel& transport,
            const word& propertiesName
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<AlphaIncompressibleTurbulenceModel> New
        (
            const volScalarField& alpha,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const TransportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName
        );


    //- Destructor
    virtual ~AlphaIncompressibleTurbulenceModel() = default;


    // Member Functions

        //- Return the effective stress tensor
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the effective stress tensor based on a given velocity field
        virtual tmp<volSymmTensorField> devReff
        (
            const volVectorField& U
        ) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Return the effective stress tensor
        virtual tmp<volSymmTensorField> devRhoReff() const;

        //- Return the effective stress tensor based on a given velocity field
        virtual tmp<volSymmTensorField> devRhoReff
        (
            const volVectorField& U
        ) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "AlphaIncompressibleTurbulenceModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
