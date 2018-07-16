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

#include "IFstream.H"
#include "openFoamTableReader.H"
#include "Vector.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Type>
void Foam::fastInterpolation2DTable<Type>::readTable()
{
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

    // Check that the data are in ascending order
    checkOrder();
    // Check if the data is uniformly distributed
    uniform_=checkUniform();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fastInterpolation2DTable<Type>::fastInterpolation2DTable()
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>(),
    boundsHandling_(fastInterpolation2DTable::WARN),
    fileName_("fileNameIsUndefined"),
    reader_(nullptr),
    uniform_(false),
    isNull_(true)
{}


template<class Type>
Foam::fastInterpolation2DTable<Type>::fastInterpolation2DTable
(
    const List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>& values,
    const boundsHandling bounds,
    const fileName& fName,
    const bool uniform,
    const Switch isNull
)
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>(values),
    boundsHandling_(bounds),
    fileName_(fName),
    reader_(nullptr),
    uniform_(uniform),
    isNull_(isNull)
{}


template<class Type>
Foam::fastInterpolation2DTable<Type>::fastInterpolation2DTable(const fileName& fName)
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>(),
    boundsHandling_(fastInterpolation2DTable::WARN),
    fileName_(fName),
    reader_(new openFoamTableReader<Type>(dictionary())),
    uniform_(false),
    isNull_(false)
{
    readTable();
}


template<class Type>
Foam::fastInterpolation2DTable<Type>::fastInterpolation2DTable(const dictionary& dict)
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>(),
    boundsHandling_(wordToBoundsHandling(dict.lookup("outOfBounds"))),
    fileName_(dict.lookup("fileName")),
    reader_(tableReader<Type>::New(dict)),
    uniform_(false),
    isNull_(false)
{
    readTable();
}


template<class Type>
Foam::fastInterpolation2DTable<Type>::fastInterpolation2DTable
(
     const fastInterpolation2DTable& interpTable
)
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>(interpTable),
    boundsHandling_(interpTable.boundsHandling_),
    fileName_(interpTable.fileName_),
    reader_(interpTable.reader_),    // note: steals reader. Used in write().
    uniform_(false)
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::fastInterpolation2DTable<Type>::interpolateValue
(
    const List<Tuple2<scalar, Type>>& data,
    const scalar lookupValue
) const
{
    label n = data.size();

    scalar minLimit = data.first().first();
    scalar maxLimit = data.last().first();

    if (lookupValue < minLimit)
    {
        switch (boundsHandling_)
        {
            case fastInterpolation2DTable::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") less than lower "
                    << "bound (" << minLimit << ")" <<nl
                    << "    Table "<<fileName_<<nl
                    << exit(FatalError);
                break;
            }
            case fastInterpolation2DTable::WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") less than lower "
                    << "bound (" << minLimit << ")" << nl
                    << "    Table "<<fileName_<<nl
                    << "    Continuing with the first entry"
                    << endl;
                // fall-through to 'CLAMP'
            }
            case fastInterpolation2DTable::CLAMP:
            {
                return data.first().second();
                break;
            }
        }
    }
    else if (lookupValue >= maxLimit)
    {
        switch (boundsHandling_)
        {
            case fastInterpolation2DTable::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") greater than upper "
                    << "bound (" << maxLimit << ")" << nl
                    << "    Table "<<fileName_<<nl
                    << exit(FatalError);
                break;
            }
            case fastInterpolation2DTable::WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") greater than upper "
                    << "bound (" << maxLimit << ")" << nl
                    << "    Table "<<fileName_<<nl
                    << "    Continuing with the last entry"
                    << endl;
                // fall-through to 'CLAMP'
            }
            case fastInterpolation2DTable::CLAMP:
            {
                return data.last().second();
                break;
            }
        }
    }

    // look for the correct range in X
    label lo = 0;
    label hi = 0;

    if (uniform_)
    {
        // For equidistant tables get range in table
        double psi = (lookupValue-data.first().first())/(data.last().first() - data.first().first());
        
        // Bound psi
        psi =std::max(std::min(psi,1.0),0.0);
        
        lo = std::floor(psi*(n-1));
        hi = std::ceil(psi*(n-1));
    }
    else
    {
        for (label i = 0; i < n; ++i)
        {
            if (lookupValue >= data[i].first())
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
        return data[lo].second();
    }
    else
    {
        Type m =
            (data[hi].second() - data[lo].second())
           /(data[hi].first() - data[lo].first());

        // normal interpolation
        return data[lo].second() + m*(lookupValue - data[lo].first());
    }
}


template<class Type>
template<class L>
Foam::label Foam::fastInterpolation2DTable<Type>::Xi
(
    const List<Tuple2<scalar, L> >& t,
    const scalar valueX,
    const bool highVal
) const
{

    label limitI = 0;

    if (highVal)
    {
        limitI = t.size() - 1;
    }
    
    if (   
        ((valueX>t[limitI].first()) && highVal)
        || ((valueX<t[limitI].first()) && !highVal)
       )
    {
        switch (boundsHandling_)
        {
            case fastInterpolation2DTable::ERROR:
            {
                FatalErrorInFunction
                    << "value (" << valueX << ") out of bounds. "
                    << "Table "<<fileName_<<nl
                    << exit(FatalError);
                break;
            }
            case fastInterpolation2DTable::WARN:
            {
                WarningInFunction
                    << "value (" << valueX << ") out of bounds. "
                    << "Table "<<fileName_<<nl
                    << endl;
                // fall-through to 'CLAMP'
            }
            case fastInterpolation2DTable::CLAMP:
            {
                return limitI;
            }
            default:
            {
                FatalErrorInFunction
                    << "Un-handled enumeration " << boundsHandling_
                    << abort(FatalError);
            }
        }
    }

    label i = 0;
    
    if (uniform_)
    {
        // For equidistant tables get range in table
        label nX = t.size();
        double psi = (valueX-t.first().first())/(t.last().first() - t.first().first());
        // Bound psi
        psi =std::max(std::min(psi,1.0),0.0);
        
        // Get high value
        if (highVal)
        {
            i =  std::ceil(psi*(nX-1));
        }
        else
        {
            // Get low value
            i = std::floor(psi*(nX-1));
        }
    }
    else
    {
        if (highVal)
        {
            label nX = t.size();
            i = 0;
            while ((i < nX) && (valueX > t[i].first()))
            {
                i++;
            }
        }
        else
        {
            i = t.size() - 1;
            while ((i > 0) && (valueX < t[i].first()))
            {
                i--;
            }
        }
    }

    return i;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
Type Foam::fastInterpolation2DTable<Type>::operator()
(
    const scalar valueX,
    const scalar valueY
) const
{
    // Considers all of the list in Y being equal
    label nX = this->size();

    const table& t = *this;

    if (nX == 0)
    {
        WarningInFunction
            << "cannot interpolate a zero-sized table - returning zero" << endl;

        return Zero;
    }
    else if (nX == 1)
    {
        // only 1 column (in X) - interpolate to find Y value
        return interpolateValue(t.first().second(), valueY);
    }
    else
    {
        // have 2-D data, interpolate

        // find low and high indices in the X range that bound valueX
        label x0i = Xi(t, valueX, false);
        label x1i = Xi(t, valueX, true);

        if (x0i == x1i)
        {
            return interpolateValue(t[x0i].second(), valueY);
        }
        else
        {
            Type y0(interpolateValue(t[x0i].second(), valueY));
            Type y1(interpolateValue(t[x1i].second(), valueY));

            // gradient in X
            scalar x0 = t[x0i].first();
            scalar x1 = t[x1i].first();
            Type mX = (y1 - y0)/(x1 - x0);

            // interpolate
            return y0 + mX*(valueX - x0);
        }
    }
}


template<class Type>
Foam::word Foam::fastInterpolation2DTable<Type>::boundsHandlingToWord
(
     const boundsHandling& bound
) const
{
    word enumName("warn");

    switch (bound)
    {
        case fastInterpolation2DTable::ERROR:
        {
            enumName = "error";
            break;
        }
        case fastInterpolation2DTable::WARN:
        {
            enumName = "warn";
            break;
        }
        case fastInterpolation2DTable::CLAMP:
        {
            enumName = "clamp";
            break;
        }
    }

    return enumName;
}


template<class Type>
typename Foam::fastInterpolation2DTable<Type>::boundsHandling
Foam::fastInterpolation2DTable<Type>::wordToBoundsHandling
(
    const word& bound
) const
{
    if (bound == "error")
    {
        return fastInterpolation2DTable::ERROR;
    }
    else if (bound == "warn")
    {
        return fastInterpolation2DTable::WARN;
    }
    else if (bound == "clamp")
    {
        return fastInterpolation2DTable::CLAMP;
    }
    else
    {
        WarningInFunction
            << "bad outOfBounds specifier " << bound << " using 'warn'" << endl;

        return fastInterpolation2DTable::WARN;
    }
}


template<class Type>
typename Foam::fastInterpolation2DTable<Type>::boundsHandling
Foam::fastInterpolation2DTable<Type>::outOfBounds
(
    const boundsHandling& bound
)
{
    boundsHandling prev = boundsHandling_;
    boundsHandling_ = bound;
    return prev;
}


template<class Type>
void Foam::fastInterpolation2DTable<Type>::checkOrder() const
{
    label n = this->size();
    const table& t = *this;

    scalar prevValue = t[0].first();

    for (label i=1; i<n; ++i)
    {
        const scalar currValue = t[i].first();

        // avoid duplicate values (divide-by-zero error)
        if (currValue <= prevValue)
        {
            FatalErrorInFunction
                << "out-of-order value: "
                << currValue << " at index " << i << nl
                << exit(FatalError);
        }
        prevValue = currValue;
    }
}

template<class Type>
bool Foam::fastInterpolation2DTable<Type>::checkUniform() const
{
    label n = this->size();
    const table& t = *this;
    
    // Check first dimension
    scalar deltaX = (t.last().first()-t.first().first())/n;
    
    // check if it is in the range of 2%
    for (int i=0;i<n-1;i++)
    {
        scalar diff = std::abs((t[i+1].first()-t[i].first())/deltaX-1);
        if  (diff > 0.02)
        {
            return false;
        }
    }
    
    // check second dimension
    List<Tuple2<scalar, Type> > column = t[0].second();
    scalar nCol = column.size();
    scalar deltaY = column.last().first()-column.first().first();
    
    for (int i=0; i<nCol-1; i++)
    {
        scalar diff = std::abs(column[i+1].first()-column[i].first())/deltaY -1;
        if (diff > 0.02)
        {
            return false;
        }
    }
    
    return true;
}


template<class Type>
Type Foam::fastInterpolation2DTable<Type>::Xderivative
(
    const scalar valueX,
    const scalar valueY
) const
{
    // Derivation of the value to the second input parameter.
    // ==> For Tables generated with CoolProp the format is (p,T)
    //     Meaning that this is the derivation to the temperature!

    // This works only well for uniform distributed data
    if (!uniform_)
    {
        WarningIn
        (
            "Type Foam::fastInterpolation2DTable<Type>::Xderivative"
            "("
            "const scalar, "
            "const scalar"
            ") const"
        )
            << "No uniform distributed data. Method not implemented yet"
        << endl;
        return pTraits<Type>::zero;
    }
    const table& t = *this;

    if (t.size() <= 1)
    {
        WarningIn
        (
            "Type Foam::fastInterpolation2DTable<Type>::Xderivative"
            "("
            "const scalar, "
            "const scalar"
            ") const"
        )
            << "cannot derivate a zero- or one-sized table - returning zero"
        << endl;
        return pTraits<Type>::zero;
    }
    
    List<Tuple2<scalar, Type> > row0;
    List<Tuple2<scalar, Type> > row1;
    // find low and high indices in the X range that bound valueX
    label x0i = Xi(t, valueX, false);
    label x1i = Xi(t, valueX, true);

    if (x0i == x1i)
    {
        // check if y1i is at max 
        if (x1i==t.size()-1)
        {
            x0i = x1i-1;
        }
        else if (x0i==0)
        {
            x1i = 1;
        }
    }
    
    // Get the bounding value of X
    Type x0 = t[x0i].first();
    Type x1 = t[x1i].first();
    
    row0 = t[x0i].second();
    row1 = t[x1i].second();
    
    // Find corresponding low and high Y index
    label y0i = Xi(row0,valueY,false);
    label y1i = Xi(row0,valueY,true);
    
    // Get bounding Y value
    Type y0 = t[x0i].second()[y0i].first();
    Type y1 = t[x0i].second()[y1i].first(); 
    

    
    // Get the 4 result values
    Type xy00 = t[x0i].second()[y0i].second();
    Type xy10 = t[x1i].second()[y0i].second();
    Type xy01 = t[x0i].second()[y1i].second();
    Type xy11 = t[x1i].second()[y1i].second();
    
    Type v1, v2;
    
    // Interpolate values
    if (y0i==y1i)
    {
        v1 = xy00;
        v2 = xy10;
    }
    else
    {
        v1 = (xy01-xy00)/(y1-y0) * (valueY-y0) + xy00;
        v2 = (xy11-xy10)/(y1-y0) * (valueY-y0) + xy10;
    }
    
    return (v2-v1)/(x1-x0);
}


template<class Type>
Type Foam::fastInterpolation2DTable<Type>::Yderivative
(
    const scalar valueX,
    const scalar valueY
) const
{
    // Derivation of the value to the second input parameter.
    // ==> For Tables generated with CoolProp the format is (p,T)
    //     Meaning that this is the derivation to the temperature!

    // This works only well for uniform distributed data
    if (!uniform_)
    {
        WarningIn
        (
            "Type Foam::fastInterpolation2DTable<Type>::Yderivative"
            "("
            "const scalar, "
            "const scalar"
            ") const"
        )
            << "No uniform distributed data. Method not implemented yet"
        << endl;
        return pTraits<Type>::zero;
    }
    const table& t = *this;

    if (t.size() <= 1)
    {
        WarningIn
        (
            "Type Foam::fastInterpolation2DTable<Type>::Yderivative"
            "("
            "const scalar, "
            "const scalar"
            ") const"
        )
            << "cannot derivate a zero- or one-sized table - returning zero"
        << endl;
        return pTraits<Type>::zero;
    }
    
    List<Tuple2<scalar, Type> > row0;
    List<Tuple2<scalar, Type> > row1;
    // find low and high indices in the X range that bound valueX
    label x0i = Xi(t, valueX, false);
    label x1i = Xi(t, valueX, true);
    
    // Get the bounding value of X
    Type x0 = t[x0i].first();
    Type x1 = t[x1i].first();
    
    row0 = t[x0i].second();
    row1 = t[x1i].second();
    // Find corresponding low and high Y index
    
    label y0i = Xi(row0,valueY,false);
    label y1i = Xi(row0,valueY,true);
    
    if (y0i == y1i)
    {
        // check if y1i is at max 
        if (y1i==t[x0i].second().size()-1)
        {
            y0i = y1i-1;
        }
        else if (y0i==0)
        {
            y1i = 1;
        }
    }
    
    // Get the 4 result values
    Type xy00 = t[x0i].second()[y0i].second();
    Type xy10 = t[x1i].second()[y0i].second();
    Type xy01 = t[x0i].second()[y1i].second();
    Type xy11 = t[x1i].second()[y1i].second();
    
    Type v1, v2;
    // Interpolate values
    if (x0i==x1i)
    {
        v1 = xy00;
        v2 = xy01;
    }
    else
    {
        v1 = (xy10-xy00)/(x1-x0) * (valueX-x0) + xy00;
        v2 = (xy11-xy01)/(x1-x0) * (valueX-x0) + xy01;
    }

    // Get bounding Y value
    Type y0 = t[x0i].second()[y0i].first();
    Type y1 = t[x0i].second()[y1i].first();    

    return (v2-v1)/(y1-y0);
}



template<class Type>
void Foam::fastInterpolation2DTable<Type>::write(Ostream& os) const
{
    os.writeKeyword("file")
        << fileName_ << token::END_STATEMENT << nl;
    os.writeKeyword("outOfBounds")
        << boundsHandlingToWord(boundsHandling_) << token::END_STATEMENT << nl;

    *this >> os;
}

// * * * * * * * * * * * * * * Operator Overload  * * * * * * * * * * * * * //
template<class Type>
inline Foam::fastInterpolation2DTable<Type> Foam::operator+
(
    const fastInterpolation2DTable<Type>& et1,
    const fastInterpolation2DTable<Type>& et2
)
{
    List<Tuple2<scalar, List<Tuple2<scalar, Type> > > > etn = et1, ett = et2;

    if (et1.size() != et2.size())
    {
	FatalErrorInFunction
            << "attempt to sum list with different sizes."
            << abort(FatalError);
    }

    if (et1.isNull_)
    {
    	return fastInterpolation2DTable<Type>
        (
            ett,
    	    et2.boundsHandling_,
    	    et2.fileName_,
            false,
    	    et2.isNull_
        );
    }
    else if (et2.isNull_)
    {
    	return fastInterpolation2DTable<Type>
        (
            etn,
    	    et1.boundsHandling_,
    	    et1.fileName_,
            false,
    	    et1.isNull_
        );
    }

    for (int i = 0 ; i < etn.size() ; ++i)
    {
	for (int j = 0 ; j < etn[i].second().size() ; ++j)
	{
	    etn[i].second()[j].second() += ett[i].second()[j].second();
	}
    }

    return fastInterpolation2DTable<Type>
    (
	etn,
	et1.boundsHandling_,
	et1.fileName_,
	false,
    false
    );
}


template<class Type>
inline Foam::fastInterpolation2DTable<Type> Foam::operator-
(
    const fastInterpolation2DTable<Type>& et1,
    const fastInterpolation2DTable<Type>& et2
)
{
    List<Tuple2<scalar, List<Tuple2<scalar, Type> > > > etn = et1, ett = et2;

    if (et1.size() != et2.size())
    {
	FatalErrorInFunction
            << "attempt to sum list with different sizes."
            << abort(FatalError);
    }

    for (int i = 0 ; i < etn.size() ; ++i)
    {
	for (int j = 0 ; j < etn[i].second().size() ; ++j)
	{
	    etn[i].second()[j].second() -= ett[i].second()[j].second();
	}
    }

    return fastInterpolation2DTable<Type>
    (
	 etn,
	 et1.boundsHandling_,
	 et1.fileName_,
     false,
	 et1.isNull_
    );
}


template<class Type>
inline Foam::fastInterpolation2DTable<Type> Foam::operator*
(
    const scalar s,
    const fastInterpolation2DTable<Type>& et
)
{
    List<Tuple2<scalar, List<Tuple2<scalar, Type> > > > etn = et;

    // I'm not sure that it saves time but I try it.
    if (s == 1 || et.isNull_)
    {
    	return fastInterpolation2DTable<Type>
        (
    	    etn,
    	    et.boundsHandling_,
    	    et.fileName_,
            false,
            et.isNull_
        );
    }
    else
    {
    	for (int i = 0 ; i < etn.size() ; ++i)
    	{
    	    for (int j = 0 ; j < etn[i].second().size() ; ++j)
    	    {
    		etn[i].second()[j].second() *= s;
    	    }
    	}
    }

    if (s == 0)
    {
    	return fastInterpolation2DTable<Type>
    	(
    	    etn,
    	    et.boundsHandling_,
    	    et.fileName_,
    	    false,
            true
    	);
    }
    else
    {
	return fastInterpolation2DTable<Type>
	(
	    etn,
	    et.boundsHandling_,
	    et.fileName_,
        false,
	    et.isNull_
	);
    }
}




// ************************************************************************* //
