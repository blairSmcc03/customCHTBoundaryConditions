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

#include "customNeumannDirichletCoupledBoundary.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mui.h"
#include "../muiconfig.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

customNeumannDirichletCoupledBoundary::
customNeumannDirichletCoupledBoundary
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


customNeumannDirichletCoupledBoundary::
customNeumannDirichletCoupledBoundary
(
    const customNeumannDirichletCoupledBoundary& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    temperatureCoupledBase(patch(), ptf)
{}


customNeumannDirichletCoupledBoundary::
customNeumannDirichletCoupledBoundary
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
}


customNeumannDirichletCoupledBoundary::
customNeumannDirichletCoupledBoundary
(
    const customNeumannDirichletCoupledBoundary& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    temperatureCoupledBase(patch(), wtcsf)
{}


customNeumannDirichletCoupledBoundary::
customNeumannDirichletCoupledBoundary
(
    const customNeumannDirichletCoupledBoundary& wtcsf
)
:
    mixedFvPatchScalarField(wtcsf),
    temperatureCoupledBase(patch(), wtcsf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void customNeumannDirichletCoupledBoundary::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    mixedFvPatchScalarField::autoMap(mapper);
    temperatureCoupledBase::autoMap(mapper);
    //mappedPatchFieldBase<scalar>::autoMap(mapper);
}


void customNeumannDirichletCoupledBoundary::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const customNeumannDirichletCoupledBoundary& tiptf =
        refCast
        <
            const customNeumannDirichletCoupledBoundary
        >(ptf);

    temperatureCoupledBase::rmap(tiptf, addr);
    //mappedPatchFieldBase<scalar>::rmap(ptf, addr);
}


tmp<Foam::scalarField>
customNeumannDirichletCoupledBoundary::kappa
(
    const scalarField& Tp
) const
{
    // Get kappa from relevant thermo
    tmp<scalarField> tk(temperatureCoupledBase::kappa(Tp));
    return tk;
}


void customNeumannDirichletCoupledBoundary::updateCoeffs()
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

    //get kappa from thermo object
    const scalarField& Tp = *this;
    const scalarField kappaTp = kappa(Tp);
    mui::point1d my_point;

    if(!boundaryInitialised){
        //Decide if this is neumann or dirichlet side, 
        // no point doing this every iteration thus global bool
        // NOTE: this cannot go in the constructor since the thermo object is not initialised
        scalarField myConvergenceCoefficient = patch().magSf()*kappaTp*patch().deltaCoeffs()*this->size();
        // Push convergence coefficient to neighbour

        //TODO: this should really be in a for loop since we could get blow up in individual cells, 
        //requires rewrite of Python solver since each cell would have to impose a unique condition
        my_point[0]=0;
        mui_ifs[0]->push("coupling", my_point, myConvergenceCoefficient[0]);
        mui_ifs[0]->commit( 0 );

        scalar totalConvergenceCoefficient = myConvergenceCoefficient[0] / mui_ifs[0]->fetch("coupling", my_point, 0, spatial_sampler, temporal_sampler);
        
        if(totalConvergenceCoefficient <= 1){
            dirichBoundary = true;
        }else{
            dirichBoundary = false;
        }
        boundaryInitialised = true;
        
    } 
    if(dirichBoundary){
        //Calculate heatflux
        scalarField Q = kappaTp*patch().magSf()*snGrad();

        // push flux to neighbour
        forAll(Q, faceI){
            my_point[0] = patch().Cf()[faceI].y();
            //multiply by number of nodes on patch to ensure sum of heat flux is consistent for non-conforming meshes
            mui_ifs[0]->push("flux", my_point, Q[faceI]*this->size());
        }
        mui_ifs[0]->commit( time );

        //fetch temperature from neighbour
        scalarField nbrIntFld = scalarField(this->size());
        std::vector<double> fetch_vals = mui_ifs[0]->fetch_values<mui::mui_config::REAL>( "temp", time, temporal_sampler );
        
        // create scalarField from fetched values
        // Possible optimisation here using MPI/MUI datatypes, might be able to fetch to a scalarField
        forAll(nbrIntFld, faceI){
            nbrIntFld[faceI] = fetch_vals[faceI];
        }
        this->refValue() = nbrIntFld;
        
    }else{
        //push internal field to neighbour
        forAll(Tp, faceI){
            // height of current cell along patch
            my_point[0] = patch().Cf()[faceI].y();
            mui_ifs[0]->push("temp", my_point, Tp[faceI]);
        }  
        mui_ifs[0]->commit(time);

        //fetch heat flux from neighbour
        scalarField internal = patchInternalField();
        scalarField nbrFluxFld = scalarField(this->size());
        std::vector<double> fetch_vals = mui_ifs[0]->fetch_values<mui::mui_config::REAL>( "flux", time, temporal_sampler );
        forAll(nbrFluxFld, faceI){
            nbrFluxFld[faceI] = (fetch_vals[faceI])/(kappaTp[faceI]*this->size());
        }
        this->refValue() = internal - nbrFluxFld/(patch().deltaCoeffs()*patch().magSf());
    }

    this->valueFraction() = 1.0;
    this->refGrad() = Zero;

    Info<<"Average T at rightWall: " <<  average(Tp) << endl;

    mixedFvPatchScalarField::updateCoeffs();
    UPstream::msgType(oldTag);  // Restore tag
}

void customNeumannDirichletCoupledBoundary::write
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
    customNeumannDirichletCoupledBoundary
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam


// ************************************************************************* //
