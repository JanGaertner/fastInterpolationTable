# fastInterpolationTable

Interpolation tables for OpenFOAM based on the interpolationTable.H and interpolation2DTable.H of the standard OpenFOAM 5.0 distribution. 

Code was also adopted from the tabulatedThermophysicalProperties directory of Yuusha0


## Direct Lookup

If the data for the tables is generated with a uniform spacing the values inside the table can be accesed directly by index, instead of looping through the complete table. 

## Changes to the extrapolation2DTable from Yuusha0/tabulatedThermophysicalProperties

 * If the distance between neighbouring points is nearly equal (0.002 variance) then it is assumed to be of uniform spacing.
 * Derivative to the X or Y direction of the 2D table are computed with that assumption of uniform spacing. Originally the Y derivation was called Tderivative in extrapolation2DTable 
 * the ```isNull_``` switch is still included but the sense of this is not checked and also not made use of. 


