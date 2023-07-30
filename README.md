# HIV_ICR

This study aims to model the within host dynamics of HIV with a focus on CD4+ T-cells and viral load. Inside this repository you will find the following scripts: In_Vivo.R and RTI.R.

In_Vivo.R is based upon the equations outlined in "An Optimal Control Approach to HIV Immunology" which can be found here: https://pdfs.semanticscholar.org/022f/3909c2f6a319345ef2d62f3526cdd5dac0ca.pdf. In_Vivo.R replicates equation (1) in "An Optimal Approach" using the parameters from tables (2) and (3). The model replicates the concentration of CD4+ T-cells, infected CD4+ T-cells, viruses, CD8+ T-cells, and activated CD8+ T-cells per mm3 inside an HIV-1 infected person over the first 365 days of infection.

RTI.R builds upon In_Vivo.R by accounting for a reverse transcriptase inhibitor (RTI), a type of antiretroviral therapy used to treat HIV. RTI's interfere with the viruses ability to productively infect a cell, thereby preventing the creation of new viruses. The effectiveness of the drug was calculated based on an industry form of measurement called the EC50 value and the minimum and maximum concentration as a function of time. The RTI used in this model is Emtricitabine. The EC50 values of Emtricitabine can be found here: https://www.accessdata.fda.gov/drugsatfda_docs/label/2008/021500s010,021896s004lbl.pdf along with the Cmin and Cmax values here: https://accp1.onlinelibrary.wiley.com/doi/abs/10.1177/0091270007300951 , the Hill Coefficient used in this model is 1, calculated by this paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3236558/.

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


To run In_Vivo.R type
```
source("In_Vivo.R")
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
source("RTI.R")
```
into the R interpretor.
Type
```
open result_rti.png
```
(or .pdf) into the terminal.


