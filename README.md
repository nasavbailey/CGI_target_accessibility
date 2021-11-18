# Roman Coronagraph Instrument (CGI) target accessibility calculator

Simple tool to show approximately which days of the year Roman may point to a given target. It is based only on the Sun-angle constraints shown in "figures/Roman-pointing-constraints" (based on a figure from Holler et al). This tool only provides approximate accessibility dates because it has several simplifying assumptions:
 * the observatory is at geocenter, not the true L2 halo orbit
 * keep-out zones for the Earth, Moon, and other solar system bodies are not included. The Earth and Moon keep-outs have only a small effect for the majority of targets. The other solar system objects' keep-out zones are so small that their impacts are negligible.

Overplots a shaded region for the 2 seasons per year when the Galactic Bulge is observable, as the Galactic Bulge Time Domain Survey (GBTDS) will typically take priority during these times
 
Required packages: astropy, matplotlib, pandas, numpy

Optional packages: astroquery; this package is used to query simbad to retrieve target coordinates; without this package the user must provide a csv file of target coordinates.

Input: A text file of either target names (see example.txt) or target names and coordinates (see "coords_example_file.csv")


# References
B. Holler, et al., "Solar system science with the Wide-Field Infrared Survey Telescope," _JATIS_, 4, 034003 (2018)
[ads entry](http://adsabs.harvard.edu/abs/2018JATIS...4c4003H)

## License and acknowledgements
yp.py was largely written by Brian Kern (Jet Propulsion Laboratory, California Institute of Technology)

Government sponsorship acknowledged. This work was carried out in part at the Jet Propulsion Laboratory, California Institute of Technology, under a contract with the National Aeronautics and Space Administration.

Copyright 2021.

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License in the LICENSE file or at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.




