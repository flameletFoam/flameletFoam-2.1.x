/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

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

Class
    Foam::combustionModels::linearInterpolation

Description
    Linear Table Interpolation

Contributors/Copyright
    2014 Hagen Müller <hagen.mueller@unibw.de>

\*---------------------------------------------------------------------------*/

#ifndef LINEARINTERPOLATION_H_
#define LINEARINTERPOLATION_H_

#include "typeInfo.H"
#include "flameletTable.H"

namespace Foam
{

class fvMesh;

class linearInterpolation
:
    public flameletTable
{

private:

	//- Disallow default bitwise copy construct
    linearInterpolation(const linearInterpolation&);

    //- Disallow default bitwise assignment
    void operator=(const linearInterpolation&);

    //- Table Values
    List<List<scalarList> > tableValues_ ;

public:

        static const List<List<scalarList> > defaultList;

        //- Runtime type information
        TypeName("linearInterpolation");

        // Constructors

        //- Construct from components
        linearInterpolation(const fvMesh& mesh, const word& tableName);

	    //- Destructor
        virtual ~linearInterpolation();

        //- Member functions

        //- Interpolate tableValues_ in a 3D-Table
        virtual inline scalar interpolate(const List<int>& ub, const scalarList& pos) const;
};

} // End Foam namespace

#endif /* LINEARINTERPOLATION_H_ */
