# A JAVA implementation of the Soil Water Balance Model SIMPEL
## General description
This repository includes a JAVA implementation of the Soil Water Balance Model SIMPEL, which has been originally developed by Georg Hörmann (University of Kiel) in a spreadsheet format (Hörmann 1997, Hörmann et al., 2007). The new implementation includes new features, while preserving relevant design criteria of the minimal original model, which makes it a versatile modelling tool:

* Flexible time step, also enabling sub-daily temporal integration of the water balance
* New interception model (Rutter et al., 1971)
* Snow model added (temperature index)
* External evapotranspiration calculation

## How to run the model
Compile the source and call java on the command line with the arguments thud.simpel.SimpelModel (i.e., the main class) and the example path (see run.sh script).

## References
Hörmann, G. (1997): SIMPEL - ein einfaches, benutzerfreundliches Bodenwassermodell zum Einsatz in der Ausbildung. Deutsche Gewässerkundliche Mitteilungen 41(2):67-72
Hörmann, G., X. Zhang and N. Fohrer (2007): Comparison of a simple and a spatially distributed hydrologic model for the simulation of a lowland catchment in Northern Germany. Ecological Modelling, 209 (1): 21-28. 
Rutter, A. J., Kershaw, K. A., Robins, P. C., and Morton, A. J. (1971). A predictive model of rainfall interception in forests, 1. Derivation of the model from observations in a plantation of Corsican pine. Agr. Meteorol. 9, 367–384.

