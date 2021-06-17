# Semi-structured-KineticModel - Model framework to predict the coculture growth under a wide range of growth conditions 

![alt text](https://github.com/AU-Wang-He-Group/Semi-structured-KineticModel/blob/main/Output/Graphical%20abstract.jpg?raw=true)


This repository contains code and experimental data for running the Semi-structured kinetic modeling for methanotroph-photoautotroph (M-P) cocultures (2021) 

This model has been implemented and tested with a series of designed experiments using *Methylomicrobium buryatense* 5GB1- *Arthrospira platensis*,
and it is expected to be applicable to other coculture systems with different species or varieties.


# 1. Folder structure

- __Code__: codes of the model for different conditions
- __Experimental Data__: all the data from the experiments used in the modeling (all analyzed data are embeded in codes as well)
- __Output__: figures generated by the model for some conditions 

![warning-icon](https://img.icons8.com/emoji/48/000000/warning-emoji.png)
__WARNING__: functionV7 needs to be downloaded to run the codes at different conditions- all codes are similar with the same variables and just initial conditions at each code are different.


# 2. Experimental Data folder structure

The Experimental Data folder includes the following files:

- __Exp.Data__: contains the individual biomass concentration and also the individual substrate consumption and product excretion rates.
   - all this set of data has been analyzed using an experimental-computational (E-C) protocol developed by our group: (URL: https://onlinelibrary.wiley.com/doi/full/10.1002/bit.27603).

   - _Different conditions_: five different series of experiments has been carried out in this work.
   1-Coculture 2-single culture 3-Gas composition 4-Light intensity 5-Inoculum ratio
     
   The discription of each file and the corresponding experiment are included in the following table:

   |File name                                 |Description                           |
   |------------------------------------------|--------------------------------------|
   |_Exp.Data_Coculture-Singleculture_        |the coculture and single cultures growth experiments  |
   |_Exp.Data_GasComposition_                 |the effect of the gas phase composition on the coculture growth     |
   |_Exp.Data_lightIntensity_                 |the effect of different light intensities on the coculture growth   |
   |_Exp.Data_InoculumRatio_                  |the effect of different inoculum ratio on the coculture growth     |
 

![info-icon](https://img.icons8.com/flat_round/48/000000/info.png)
__NOTE__: photoautotroph has been specified with Green color and methanotroph with Orange color 
   to make reading the excel file easier.


- __Raw.Data__: contains measurements of each sample points for all experiments.
   - these data will not be used in the code directly and it is just has been provided as additional information.
   - raw data of the measurements included:
      - Optical density (OD 670)
      - Inorganic carbon concentration (IC) (mg/L)
      - Gas area and concentrations using Gas Chromatography (GC) for -Oxygen (O2) -Methane (CH4) -Carbon dioxide (CO2)


# 3. Code folder structure

- User can produce all figures available in the paper: 

- In order to make it clear, we use the following table to introduce each file:

![Table1](https://user-images.githubusercontent.com/67964457/122443411-2f786d00-cf65-11eb-9ea3-8d64c39a4b36.jpg)

- __Matlab codes__: The discription of each file and the corresponding simulation are included in the following table:

   |File name                                 |Description                           |
   |------------------------------------------|--------------------------------------|
   |_CocultureVsSingle_        | To reproduce the analysis of Case A and Case E |
   |_GasComposition_                 |To reproduce the analysis of Case C     |
   |_LightIntensity_                 |To reproduce the analysis of Case B   |
   |_InoculumRatio_                  |To reproduce the analysis of Case D     |
   |_functionV7_                  |the function used for runnig all above codes     |
   |_WithoutSelfShad_                  |A code for investigating the absence of self-shading effect     |
   |_functionWO_                  |the function used with _WithoutSelfShad_ file     |
   
![warning-icon](https://img.icons8.com/emoji/48/000000/warning-emoji.png)
__WARNING__: Function file needs to be downloaded to run the codes at different conditions.


## 3.1. How to read the codes (preparation)
- Flowchart of semi-structure modeling of the M-P coculture: 

 ![image](https://user-images.githubusercontent.com/67964457/122436706-a827fb00-cf5e-11eb-9c9e-f56282463082.png)


- The structure of the codes has been shown in the following subsections:

### 3.1.1. Experimental data

each code has been started with the experimental section, which is copied from analyzed experimental data (__Exp.Data__).

![warning-icon](https://img.icons8.com/emoji/48/000000/warning-emoji.png)
__WARNING__: If you wish to use your own data start by cleaning each one of the variables in _Experimental data_ section in the code: 
Further explanation about -Name- of each variable has been provided inside the codes.


### 3.1.2. Parameters

- Insert Gas composition:
 initial percentage of each gas component in the system


- Insert light intensity:
light intensity needs to be entered for photoautotroph growth based on (umol photon/m2/s)

- Insert volume of liquid and gas phase (L):
volume of liquid and gas phase on the bioreactor

- photoautotroph Yields:
needed for calculating amount of biomass and _in situ_ gas consumption production rate

- methanotroph yields:
needed for calculating amount of biomass and _in situ_ gas consumption production rate

- Monod parameters:
for calculation of the species growth rate 

- time(hours):
all sample points in the experiment. to give the running time.

### 3.1.3. Initialization 

- Insert initial individual biomass (inoculation):
concentration of each species in the system (gDCW/L)

- Insert total dissoved inorganic carbon in the liquid (mmol/L)(if any): 
amount of bicarbonate or carbonate salt in the medium

- Henry's constant (M/atm):
of all gases used in the system


## 3.2. Solution of the ODE equations

The codes use each time steps to solve the equations.

![info-icon](https://img.icons8.com/flat_round/48/000000/info.png)
__NOTE__: Time steps need to be defined in this section.

## 3.2. Output (plots)

Model outputs will be genrated in different plots such as 

(i) biomass prediction
         
 (ii) population ratio
 
 (iii) Change in Gas Phase
 
 etc. further explanation about the used variable are provided inside the codes.


# 4. Contacts

If you need assistance using these analysis scripts or to adjust them to your specific aims, 
do contact us at:

_Kiumars Badr_: kzb0054 [at] auburn.edu

_Jin Wang_: wang [at] auburn.edu

_Wang group Jun 2021_
