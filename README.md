# SubtypeDrug: Prioritization of Candidate Cancer Subtype Specific Drugs
A systematic biology tool was developed to prioritize cancer subtype-specific drugs by integrating genetic perturbation, drug action, biological pathway, and cancer subtype. The capabilities of this tool include inferring patient-specific subpathway activity profiles in the context of gene expression profiles with subtype labels, calculating differentially expressed subpathways based on cultured human cells treated with drugs in the 'cMap' (connectivity map) database, prioritizing cancer subtype specific drugs according to drug-disease reverse association score based on subpathway, and visualization of results (Castelo (2013) <doi:10.1186/1471-2105-14-7>; Han et al (2019) <doi:10.1093/bioinformatics/btz894>; Lamb and Justin (2006) <doi:10.1126/science.1132939>).   


## Installation method：
```R
#1. 
library(devtools); 
install_github("hanjunwei-lab/SubtypeDrug")

#2.
install.packages("SubtypeDrug")

#Use：
library(SubtypeDrug)
```

`SubtypeDrug` is published in Bioinformatics. Please cite the following article when using `SubtypeDrug`: https://doi.org/10.1093/bioinformatics/btab011.
