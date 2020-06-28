Those normalized DSMs are produced in batch mode using lastools and some own development.

The files (8bit) encode the height above ground in dm, all heights above 25.5m are set to 25.5 (i.e. grey value 255)

There was no manual check, so no quantitative analysis, nor check for gross errors.


==> We do not guarantee error free data here, this is just for researchers to help using height data, other than the absolute DSM.

Per tile you find two files ..._lastools.jpg: produced directly by lastools (see shell script). There, the ground points get tesselated to have a closed ground surface, then per above-ground point the difference is taken. This leads to some artefacts where the ground interpolation was not correct, eg when many above ground points are labelled as ground.

the ..._ownapproach.jpg was computed by avoiding tesselation: per above-ground-point simply the closest on-ground-point is searched and the respective height substracted. There are of course also artefacts, but a quantitative comparision between the two was not conducted.

the *laz files are the original DSMs converted to laz. So everyone interested could perform the ground filtering, as well.



