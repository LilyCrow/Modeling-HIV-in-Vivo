# HIV_ICR

This study aims to model the within host dynamics of HIV. Inside this repository you will find the following four scripts: 2017_InVivo.R, PI.R, RTI.R, and RTI_PI.R.

Each script is based upon the equations outlined in "The In Vivo Dynamics of HIV Infection with the Influence of Cytotoxic T Lymphocyte Cells" which can be found here: https://www.hindawi.com/journals/isrn/2017/2124789/. 2017_InVivo.R replicates the plots outlined in "The In Vivo Dynamics of HIV". The other three scripts-- PI.R, RTI.R, and PI_RTI.R-- account for two types of drug therapies in different variations (individually and combined). 

In order to run each script you will need to install an R interpretor (unless you already have one). To do so type
```
sudo apt install r-base
```
into the terminal.

Open the R interpretor by typing
```
R
```
into the terminal.


To run 2017_InVivo.R type
```
source("2017_InVivo.R")
```
into the R interpretor.

If the file ran correctly
[1] "wrote output to file  result_hiv.png"
[1] "wrote output to file  result_hiv.pdf"
will be printed out.

Type
```
q()
```
to exit the R interpretor, or open another terminal.

Type
```
open result_hiv.png
```
(or reuslt_hiv.pdf) to view the plot.


To run the model with a protease inhibitor type
```
source("PI.R")
```
into the R interpretor.
Type
```
open result_pi.png
```
(or .pdf) into the terminal.

To run the model with a reverse transcriptase inhibitor type
```
source("RTI.R")
```
into the R interpretor.
Type
```
open result_rti.png
```
(or .pdf) into the terminal.

To run the model with combination therapy (a PI and an RTI) type
```
source("PI_RTI.R")
```
into the R interpretor.
Type
```
open result_cARV.png
```
(or .pdf) into the terminal.


