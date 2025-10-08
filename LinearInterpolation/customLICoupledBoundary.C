/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "customLICoupledBoundary.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mui.h"
#include "../muiconfig.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

customLICoupledBoundary::
customLICoupledBoundary
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch())
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 1.0;
}


customLICoupledBoundary::
customLICoupledBoundary
(
    const customLICoupledBoundary& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    temperatureCoupledBase(patch(), ptf)
{}


customLICoupledBoundary::
customLICoupledBoundary
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict)
{
    this->readValueEntry(dict, IOobjectOption::MUST_READ);

    //get rank from MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // intialise MUI interface
    // TODO: Generalise for arbitray postition not only left->right
    std::vector<std::string> ifsName;
    ifsName.emplace_back("ifs" + std::to_string(rank+1));
    mui_ifs=mui::create_uniface<mui::mui_config>( "OpenFOAM", ifsName );	    
    Info <<  "OF solver Finsihed creating MUI interface " << endl; 


    //set dt to minimum value among all solvers
    Time& runTime = const_cast<Time&>(this->db().time());
    scalar dt = runTime.deltaTValue();
    Info << "DT before: " << dt << endl;

    MPI_Allreduce(&dt, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    runTime.setDeltaT(dt, false);
    Info << "DT after: " << runTime.deltaTValue() << endl;

    

    if (this->readMixedEntries(dict))
    {
        // Full restart
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = Zero;
        valueFraction() = 1.0;
    }

// This blocks (crashes) with more than two worlds!
//
///    // Store patch value as initial guess when running in database mode
///    mappedPatchFieldBase<scalar>::initRetrieveField
///    (
///        this->internalField().name(),
///        *this
///    );
///    mappedPatchFieldBase<scalar>::initRetrieveField
///    (
///        this->internalField().name() + "_weights",
///        this->patch().deltaCoeffs()
///    );
}


customLICoupledBoundary::
customLICoupledBoundary
(
    const customLICoupledBoundary& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    temperatureCoupledBase(patch(), wtcsf)
{}


customLICoupledBoundary::
customLICoupledBoundary
(
    const customLICoupledBoundary& wtcsf
)
:
    mixedFvPatchScalarField(wtcsf),
    temperatureCoupledBase(patch(), wtcsf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void customLICoupledBoundary::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    mixedFvPatchScalarField::autoMap(mapper);
    temperatureCoupledBase::autoMap(mapper);
    //mappedPatchFieldBase<scalar>::autoMap(mapper);
}


void customLICoupledBoundary::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const customLICoupledBoundary& tiptf =
        refCast
        <
            const customLICoupledBoundary
        >(ptf);

    temperatureCoupledBase::rmap(tiptf, addr);
    //mappedPatchFieldBase<scalar>::rmap(ptf, addr);
}


tmp<Foam::scalarField>
customLICoupledBoundary::kappa
(
    const scalarField& Tp
) const
{
    // Get kappa from relevant thermo
    tmp<scalarField> tk(temperatureCoupledBase::kappa(Tp));
    return tk;
}


void customLICoupledBoundary::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    const int oldTag = UPstream::incrMsgType();

    // get current time for MUI
    scalar time =  this->db().time().value();  

    // calculate kP/deltaP and push to neighbour
    const scalarField& Tp = *this;
    // deltaCoeffs = 1 / (distance between internalField and patch) for all patch faces
    const scalarField myKDelta = kappa(Tp)*patch().deltaCoeffs();

    //TODO: this should be moved inside loop since kappa is not always constant along boundary
    mui::point1d my_point;my_point[0]=0;
    mui_ifs[0]->push( "weight", my_point, myKDelta[0]);

    //push internal field to neighbour
    scalarField internal = patchInternalField();
    forAll(internal, faceI){
        // height of current cell along patch
        my_point[0] = patch().Cf()[faceI].y();
        mui_ifs[0]->push("temp", my_point, internal[faceI]);
    }  
    mui_ifs[0]->commit(time);

    // fetch kN/deltaN from neighbour and create scalarField
    my_point[0] = 0;
    scalar nbrKDelta = mui_ifs[0]->fetch("weight", my_point, time, spatial_sampler, temporal_sampler);
    scalarField nbrKDeltaField = scalarField(this->size(), nbrKDelta);

    this->valueFraction() = nbrKDeltaField/(nbrKDeltaField + myKDelta);


    //fetch internal field from neighbour
    scalarField nbrIntFld = scalarField(this->size());
    std::vector<double> fetch_vals = mui_ifs[0]->fetch_values<mui::mui_config::REAL>( "temp", time, temporal_sampler );
    
    //create scalarField from fetched values
    // Possible optimisation here using MPI/MUI datatypes, might be able to fetch to a scalarField
    forAll(nbrIntFld, faceI){
        nbrIntFld[faceI] = fetch_vals[faceI];
    }
    this->refValue() = nbrIntFld;

    

    this->refGrad() = Zero;
    mixedFvPatchScalarField::updateCoeffs();
    
    UPstream::msgType(oldTag);  // Restore tag
}

void customLICoupledBoundary::write
(
    Ostream& os
) const
{
    mixedFvPatchField<scalar>::write(os);
    temperatureCoupledBase::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    customLICoupledBoundary
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam


// ************************************************************************* //
