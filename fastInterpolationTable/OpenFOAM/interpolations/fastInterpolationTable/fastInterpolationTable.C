/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "fastInterpolationTable.H"
#include "openFoamTableReader.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Type>
void Foam::fastInterpolationTable<Type>::readTable()
{
    // preserve the original (unexpanded) fileName to avoid absolute paths
    // appearing subsequently in the write() method
    fileName fName(fileName_);

    fName.expand();

    // Read data from file
    reader_()(fName, *this);

    if (this->empty())
    {
        FatalErrorInFunction
            << "table read from " << fName << " is empty" << nl
            << exit(FatalError);
    }

    // Check that the data are okay
    check();
    uniform_ = checkUniform();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fastInterpolationTable<Type>::fastInterpolationTable()
:
    List<Tuple2<scalar, Type>>(),
    boundsHandling_(fastInterpolationTable::WARN),
    fileName_("fileNameIsUndefined"),
    reader_(nullptr),
    uniform_(false),
    isNull_(true)
{}


template<class Type>
Foam::fastInterpolationTable<Type>::fastInterpolationTable
(
    const List<Tuple2<scalar, Type>>& values,
    const boundsHandling bounds,
    const fileName& fName,
    const Switch isNull
)
:
    List<Tuple2<scalar, Type>>(values),
    boundsHandling_(bounds),
    fileName_(fName),
    reader_(nullptr),
    uniform_(false),
    isNull_(isNull)
{}


template<class Type>
Foam::fastInterpolationTable<Type>::fastInterpolationTable(const fileName& fName)
:
    List<Tuple2<scalar, Type>>(),
    boundsHandling_(fastInterpolationTable::WARN),
    fileName_(fName),
    reader_(new openFoamTableReader<Type>(dictionary())),
    uniform_(false),
    isNull_(false)
{
    readTable();
}


template<class Type>
Foam::fastInterpolationTable<Type>::fastInterpolationTable(const dictionary& dict)
:
    List<Tuple2<scalar, Type>>(),
    boundsHandling_(wordToBoundsHandling(dict.lookup("outOfBounds"))),
    fileName_(dict.lookup("fileName")),
    reader_(tableReader<Type>::New(dict)),
    uniform_(false),
    isNull_(false)
{
    readTable();
}


template<class Type>
Foam::fastInterpolationTable<Type>::fastInterpolationTable
(
     const fastInterpolationTable& interpTable
)
:
    List<Tuple2<scalar, Type>>(interpTable),
    boundsHandling_(interpTable.boundsHandling_),
    fileName_(interpTable.fileName_),
    reader_(interpTable.reader_),    // note: steals reader. Used in write().
    uniform_(interpTable.uniform_),
    isNull_(interpTable.isNull_)
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::word Foam::fastInterpolationTable<Type>::boundsHandlingToWord
(
     const boundsHandling& bound
) const
{
    word enumName("warn");

    switch (bound)
    {
        case fastInterpolationTable::ERROR:
        {
            enumName = "error";
            break;
        }
        case fastInterpolationTable::WARN:
        {
            enumName = "warn";
            break;
        }
        case fastInterpolationTable::CLAMP:
        {
            enumName = "clamp";
            break;
        }
        case fastInterpolationTable::REPEAT:
        {
            enumName = "repeat";
            break;
        }
    }

    return enumName;
}


template<class Type>
typename Foam::fastInterpolationTable<Type>::boundsHandling
Foam::fastInterpolationTable<Type>::wordToBoundsHandling
(
    const word& bound
) const
{
    if (bound == "error")
    {
        return fastInterpolationTable::ERROR;
    }
    else if (bound == "warn")
    {
        return fastInterpolationTable::WARN;
    }
    else if (bound == "clamp")
    {
        return fastInterpolationTable::CLAMP;
    }
    else if (bound == "repeat")
    {
        return fastInterpolationTable::REPEAT;
    }
    else
    {
        WarningInFunction
            << "bad outOfBounds specifier " << bound << " using 'warn'" << endl;

        return fastInterpolationTable::WARN;
    }
}


template<class Type>
typename Foam::fastInterpolationTable<Type>::boundsHandling
Foam::fastInterpolationTable<Type>::outOfBounds
(
    const boundsHandling& bound
)
{
    boundsHandling prev = boundsHandling_;
    boundsHandling_ = bound;
    return prev;
}


template<class Type>
void Foam::fastInterpolationTable<Type>::check() const
{
    label n = this->size();
    scalar prevValue = List<Tuple2<scalar, Type>>::operator[](0).first();

    for (label i=1; i<n; ++i)
    {
        const scalar currValue =
            List<Tuple2<scalar, Type>>::operator[](i).first();

        // avoid duplicate values (divide-by-zero error)
        if (currValue <= prevValue)
        {
            FatalErrorInFunction
                << "out-of-order value: "
                << currValue << " at index " << i << nl
                << "    Table "<<fileName_<<nl
                << exit(FatalError);
        }
        prevValue = currValue;
    }
}


template<class Type>
bool Foam::fastInterpolationTable<Type>::checkUniform()
{
    label n = this->size();
    
    const List<Tuple2<scalar, Type>>& table = *this;
    
    scalar deltaX = (table.last().first()-table.first().first())/n;
    
    
    for (label i=0;i<n-1;i++)
    {
        scalar diffX = std::abs(((table[i+1].first()-table[i].first())/deltaX)-1);
        if (diffX > 0.02)
        {
            return false;
        }
    }
    
    return true;
}

template<class Type>
void Foam::fastInterpolationTable<Type>::write(Ostream& os) const
{
    os.writeKeyword("file")
        << fileName_ << token::END_STATEMENT << nl;
    os.writeKeyword("outOfBounds")
        << boundsHandlingToWord(boundsHandling_) << token::END_STATEMENT << nl;
    if (reader_.valid())
    {
        reader_->write(os);
    }
}


template<class Type>
Type Foam::fastInterpolationTable<Type>::rateOfChange(const scalar value) const
{
    label n = this->size();

    if (n <= 1)
    {
        // There are not enough entries to provide a rate of change
        return 0;
    }

    scalar minLimit = List<Tuple2<scalar, Type>>::operator[](0).first();
    scalar maxLimit = List<Tuple2<scalar, Type>>::operator[](n-1).first();
    scalar lookupValue = value;

    if (lookupValue < minLimit)
    {
        switch (boundsHandling_)
        {
            case fastInterpolationTable::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") underflow" << nl
                    << "    Table "<<fileName_<<nl
                    << exit(FatalError);
                break;
            }
            case fastInterpolationTable::WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") underflow" << nl
                    << "    Table "<<fileName_<<nl
                    << "    Zero rate of change."
                    << endl;
                // fall-through to 'CLAMP'
            }
            case fastInterpolationTable::CLAMP:
            {
                return 0;
                break;
            }
            case fastInterpolationTable::REPEAT:
            {
                // adjust lookupValue to >= minLimit
                scalar span = maxLimit-minLimit;
                lookupValue = fmod(lookupValue-minLimit, span) + minLimit;
                break;
            }
        }
    }
    else if (lookupValue >= maxLimit)
    {
        switch (boundsHandling_)
        {
            case fastInterpolationTable::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") overflow" << nl
                    << "    Table "<<fileName_<<nl
                    << exit(FatalError);
                break;
            }
            case fastInterpolationTable::WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") overflow" << nl
                    << "    Table "<<fileName_<<nl
                    << "    Zero rate of change."
                    << endl;
                // fall-through to 'CLAMP'
            }
            case fastInterpolationTable::CLAMP:
            {
                return 0;
                break;
            }
            case fastInterpolationTable::REPEAT:
            {
                // adjust lookupValue <= maxLimit
                scalar span = maxLimit-minLimit;
                lookupValue = fmod(lookupValue-minLimit, span) + minLimit;
                break;
            }
        }
    }

    label lo = 0;
    label hi = 0;
    
    List<Tuple2<scalar, Type>>& table = *this;
    
    // look for the correct range
    if (uniform_)
    {
        // For equidistant tables get range in table
        scalar psi = (lookupValue-table.first().first())/(table.last().first() - table.first().first());
        
        // Bound psi
        psi =std::max(std::min(psi,1.0),0.0);
        
        lo = std::floor(psi*(n-1));
        hi = std::ceil(psi*(n-1));
    }
    else
    {
        for (label i = 0; i < n; ++i)
        {
            if (lookupValue >= List<Tuple2<scalar, Type>>::operator[](i).first())
            {
                lo = hi = i;
            }
            else
            {
                hi = i;
                break;
            }
        }
    }

    if (lo == hi)
    {
        // we are at the end of the table - or there is only a single entry
        return 0;
    }
    else if (hi == 0)
    {
        // this treatment should should only occur under these conditions:
        //  -> the 'REPEAT' treatment
        //  -> (0 <= value <= minLimit)
        //  -> minLimit > 0
        // Use the value at maxLimit as the value for value=0
        lo = n - 1;

        return
        (
            (
                List<Tuple2<scalar, Type>>::operator[](hi).second()
              - List<Tuple2<scalar, Type>>::operator[](lo).second()
            )
           /(
               List<Tuple2<scalar, Type>>::operator[](hi).first()
             + minLimit
             - List<Tuple2<scalar, Type>>::operator[](lo).first()
            )
        );
    }
    else
    {
        // normal rate of change
        return
        (
            (
                List<Tuple2<scalar, Type>>::operator[](hi).second()
              - List<Tuple2<scalar, Type>>::operator[](lo).second()
            )
           /(
                List<Tuple2<scalar, Type>>::operator[](hi).first()
              - List<Tuple2<scalar, Type>>::operator[](lo).first()
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
const Foam::Tuple2<Foam::scalar, Type>&
Foam::fastInterpolationTable<Type>::operator[](const label i) const
{
    label ii = i;
    label n  = this->size();

    if (n <= 1)
    {
        ii = 0;
    }
    else if (ii < 0)
    {
        switch (boundsHandling_)
        {
            case fastInterpolationTable::ERROR:
            {
                FatalErrorInFunction
                    << "index (" << ii << ") underflow" << nl
                    << "    Table "<<fileName_<<nl
                    << exit(FatalError);
                break;
            }
            case fastInterpolationTable::WARN:
            {
                WarningInFunction
                    << "index (" << ii << ") underflow" << nl
                    << "    Continuing with the first entry"
                    << "    Table "<<fileName_<<nl
                    << endl;
                // fall-through to 'CLAMP'
            }
            case fastInterpolationTable::CLAMP:
            {
                ii = 0;
                break;
            }
            case fastInterpolationTable::REPEAT:
            {
                while (ii < 0)
                {
                    ii += n;
                }
                break;
            }
        }
    }
    else if (ii >= n)
    {
        switch (boundsHandling_)
        {
            case fastInterpolationTable::ERROR:
            {
                FatalErrorInFunction
                    << "index (" << ii << ") overflow" << nl
                    << "    Table "<<fileName_<<nl
                    << exit(FatalError);
                break;
            }
            case fastInterpolationTable::WARN:
            {
                WarningInFunction
                    << "index (" << ii << ") overflow" << nl
                    << "    Continuing with the last entry"
                    << "    Table "<<fileName_<<nl
                    << endl;
                // fall-through to 'CLAMP'
            }
            case fastInterpolationTable::CLAMP:
            {
                ii = n - 1;
                break;
            }
            case fastInterpolationTable::REPEAT:
            {
                while (ii >= n)
                {
                    ii -= n;
                }
                break;
            }
        }
    }

    return List<Tuple2<scalar, Type>>::operator[](ii);
}


template<class Type>
Type Foam::fastInterpolationTable<Type>::operator()(const scalar value) const
{
    label n = this->size();

    if (n <= 1)
    {
        return List<Tuple2<scalar, Type>>::operator[](0).second();
    }

    scalar minLimit = List<Tuple2<scalar, Type>>::operator[](0).first();
    scalar maxLimit = List<Tuple2<scalar, Type>>::operator[](n-1).first();
    scalar lookupValue = value;

    if (lookupValue < minLimit)
    {
        switch (boundsHandling_)
        {
            case fastInterpolationTable::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") underflow" << nl
                    << "    Table "<<fileName_<<nl
                    << exit(FatalError);
                break;
            }
            case fastInterpolationTable::WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") underflow" << nl
                    << "    Continuing with the first entry"
                    << "    Table "<<fileName_<<nl
                    << endl;
                // fall-through to 'CLAMP'
            }
            case fastInterpolationTable::CLAMP:
            {
                return List<Tuple2<scalar, Type>>::operator[](0).second();
                break;
            }
            case fastInterpolationTable::REPEAT:
            {
                // adjust lookupValue to >= minLimit
                scalar span = maxLimit-minLimit;
                lookupValue = fmod(lookupValue-minLimit, span) + minLimit;
                break;
            }
        }
    }
    else if (lookupValue >= maxLimit)
    {
        switch (boundsHandling_)
        {
            case fastInterpolationTable::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") overflow" << nl
                    << "    Table "<<fileName_<<nl
                    << exit(FatalError);
                break;
            }
            case fastInterpolationTable::WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") overflow" << nl
                    << "    Continuing with the last entry"
                    << "    Table "<<fileName_<<nl
                    << endl;
                // fall-through to 'CLAMP'
            }
            case fastInterpolationTable::CLAMP:
            {
                return List<Tuple2<scalar, Type>>::operator[](n-1).second();
                break;
            }
            case fastInterpolationTable::REPEAT:
            {
                // adjust lookupValue <= maxLimit
                scalar span = maxLimit-minLimit;
                lookupValue = fmod(lookupValue-minLimit, span) + minLimit;
                break;
            }
        }
    }

    label lo = 0;
    label hi = 0;

    // look for the correct range
    for (label i = 0; i < n; ++i)
    {
        if (lookupValue >= List<Tuple2<scalar, Type>>::operator[](i).first())
        {
            lo = hi = i;
        }
        else
        {
            hi = i;
            break;
        }
    }

    if (lo == hi)
    {
        // we are at the end of the table - or there is only a single entry
        return List<Tuple2<scalar, Type>>::operator[](hi).second();
    }
    else if (hi == 0)
    {
        // this treatment should should only occur under these conditions:
        //  -> the 'REPEAT' treatment
        //  -> (0 <= value <= minLimit)
        //  -> minLimit > 0
        // Use the value at maxLimit as the value for value=0
        lo = n - 1;

        return
        (
            List<Tuple2<scalar, Type>>::operator[](lo).second()
          + (
                List<Tuple2<scalar, Type>>::operator[](hi).second()
              - List<Tuple2<scalar, Type>>::operator[](lo).second()
            )
           *(lookupValue / minLimit)
        );
    }
    else
    {
        // normal interpolation
        return
        (
            List<Tuple2<scalar, Type>>::operator[](lo).second()
          + (
                List<Tuple2<scalar, Type>>::operator[](hi).second()
              - List<Tuple2<scalar, Type>>::operator[](lo).second()
            )
           *(
                lookupValue
              - List<Tuple2<scalar, Type>>::operator[](lo).first()
            )
           /(
                List<Tuple2<scalar, Type>>::operator[](hi).first()
              - List<Tuple2<scalar, Type>>::operator[](lo).first()
            )
        );
    }
}


// ************************************************************************* //
