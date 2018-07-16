/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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

Class
    Foam::fastInterpolation2DTable

Description
    2D table interpolation. The data must be in ascending order in both
    dimensions x and y.

SourceFiles
    fastInterpolation2DTable.C

\*---------------------------------------------------------------------------*/

#ifndef fastInterpolation2DTable_H
#define fastInterpolation2DTable_H

#include "List.H"
#include "Tuple2.H"
#include "tableReader.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type> class fastInterpolation2DTable;

template<class Type>
inline fastInterpolation2DTable<Type> operator+
(
    const fastInterpolation2DTable<Type>&,
    const fastInterpolation2DTable<Type>&
);


template<class Type>
inline fastInterpolation2DTable<Type> operator-
(
    const fastInterpolation2DTable<Type>&,
    const fastInterpolation2DTable<Type>&
);


template<class Type>
inline fastInterpolation2DTable<Type> operator*
(
    const scalar,
    const fastInterpolation2DTable<Type>&
);


/*---------------------------------------------------------------------------*\
                    Class fastInterpolation2DTable Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fastInterpolation2DTable
:
    public List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>
{
public:

    // Public data types

        //- Enumeration for handling out-of-bound values
        enum boundsHandling
        {
            ERROR,          //!< Exit with a FatalError
            WARN,           //!< Issue warning and clamp value (default)
            CLAMP           //!< Clamp value to the start/end value
        };

        //- Cconvenience typedef
        typedef List<Tuple2<scalar, List<Tuple2<scalar, Type>>>> table;


private:

    // Private data

        //- Enumeration for handling out-of-bound values
        boundsHandling boundsHandling_;

        //- File name
        fileName fileName_;

        //- The actual reader
        autoPtr<tableReader<Type>> reader_;
        
        //- If data is dirstributed uniformly
        bool uniform_;
        
        //- true or yes if all the values are 0
        Switch isNull_;


    // Private Member Functions

        //- Read the table of data from file
        void readTable();

        //- Return interpolated value in List
        Type interpolateValue
        (
            const List<Tuple2<scalar, Type>>& data,
            const scalar
        ) const;

        //- Return an X index from the matrix
        template<class L>
        label Xi
        (
            const List<Tuple2<scalar, L> >& t,
            const scalar valueX,
            const bool highVal
        ) const;
        
        //- Check if data is uniform 
        bool checkUniform() const;


public:

    // Constructors

        //- Construct null
        fastInterpolation2DTable();

        //- Construct from components
        fastInterpolation2DTable
        (
            const List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>& values,
            const boundsHandling bounds,
            const fileName& fName,
            const bool uniform,
            const Switch isNull
        );

        //- Construct given the name of the file containing the table of data
        fastInterpolation2DTable(const fileName& fName);

        //- Construct by reading the fileName and boundsHandling from dictionary
        fastInterpolation2DTable(const dictionary& dict);

        //- Construct copy
        fastInterpolation2DTable(const fastInterpolation2DTable& interpTable);


    // Member Functions

        //- Return the out-of-bounds handling as a word
        word boundsHandlingToWord(const boundsHandling& bound) const;

        //- Return the out-of-bounds handling as an enumeration
        boundsHandling wordToBoundsHandling(const word& bound) const;

        //- Set the out-of-bounds handling from enum, return previous setting
        boundsHandling outOfBounds(const boundsHandling& bound);

        //- Check that list is monotonically increasing
        //  Exit with a FatalError if there is a problem
        void checkOrder() const;
        
        //- Derivation of the first variable
        Type Xderivative(const scalar valueX,const scalar valueY) const;
        
        //- Derivate a tabulated value.
        //  Simple derivation for testing purpose. Must be improved.
        Type Yderivative(const scalar valueX,const scalar valueY) const;
        

        //- Write
        void write(Ostream& os) const;


    // Member Operators

        //- Return an element of constant Tuple2<scalar, Type>
        const List<Tuple2<scalar, Type>>& operator[](const label) const;

        //- Return an interpolated value
        Type operator()(const scalar, const scalar) const;
        
// Friend operators

        friend fastInterpolation2DTable operator+ <Type>
        (
            const fastInterpolation2DTable&,
            const fastInterpolation2DTable&
        );

        friend fastInterpolation2DTable operator- <Type>
        (
            const fastInterpolation2DTable&,
            const fastInterpolation2DTable&
        );

        friend fastInterpolation2DTable operator* <Type>
        (
            const scalar,
            const fastInterpolation2DTable&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fastInterpolation2DTable.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
