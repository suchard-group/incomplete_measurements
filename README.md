# Inferring phenotypic trait evolution on large trees with many incomplete measurements

This repository contains the necessary data and scripts to reproduce the analyses and figures in the [manuscript](https://arxiv.org/abs/1906.03222) entitled _Inferring phenotypic trait evolution on large trees with many incomplete measurements_ by Gabriel Hassler, Max R. Tolkoff, William L. Allen, Lam Si Tung Ho, Philippe Lemey, and Marc A. Suchard.

The instructions below should be sufficient to perform all analyses on OSX, Windows, and Linux, but please raise an issue on this GitHub repository if you encounter any problems.

## File Structure

- __beast__
    * __beast.jar__ - Compressed [BEAST](http://beast.community/) source code. This file was compiled from the version of BEAST located [here](https://github.com/beast-dev/beast-mcmc/tree/0d143129454605510b9112ff9abf29ef7a7bf670).
- __data__
    * __ele_1307_sm_sa1.tre__ - Mammalian phylogeny from [Fritz et al. (2009)](https://doi.org/10.1111/j.1461-0248.2009.01307.x). Downloaded from [here](https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fj.1461-0248.2009.01307.x&file=ELE_1307_sm_SA1.tre).
    * __hiv_dates.csv__ - Sampling dates associated with HIV data.
    * __hiv_newick.txt__ - Newick HIV tree. Tree originated from ML tree in [Blanquart et al. (2017)](https://doi.org/10.1371/journal.pbio.1002608) (processing described in manuscript).
    * __hiv_processed_data.csv__ - Processed HIV data used in analysis (processing described in manuscript).
    * __mammals_log_data.csv__ - CSV file that stores logged mammals data. Missing values are prepresented by NaN. File was created by running Julia script `./scripts/mammals_data.jl`.
    * __mammals_newick.txt__ - The __mammalST_MSW05_bestDates__ tree manually copied from  __ele_1307_sm_sa1.tre__.
    * __mammals_trimmed_newick.txt__ - Mammals newick file with taxa not present in the data set removed. File was created by running Julia script `./scripts/mammals_data.jl`.
    * __PanTHERIA_1-0_WR05_Aug2008.txt__ - Mammalian phenotype data from PanTHERIA database. Downloaded from [here](http://esapubs.org/archive/ecol/E090/184/).
    * __prokaryotes_newick.txt__ - Obtained tree log file after running `./xml/prokaryotes.xml` in BEAST, resampled tree every 100,000 trees to reduce the file size using the BEAST software [LogCombiner](https://beast.community/logcombiner), determined the maximum clade credibility tree using [TreeAnnotator](https://beast.community/treeannotator), and output a newick file using [FigTree](https://beast.community/figtree). This newick had a single negative branch length, which was manually edited to 0.0 with the child banch lengths shortened accordingly.
    * __prokaryotes_processed_data.csv__ - Processed prokaryote data used in analysis (processing described in manuscript).
- __logs__
    * __hiv_prediction__ - directory for storing BEAST log files in HIV posterior predictive power analysis.
    * __PCMBase_timing__ - directory for storing results of speed comparison between [PCMBaseCpp](https://github.com/venelin/PCMBaseCpp) and BEAST.
    * __simulation_study__ - directory for storing BEAST log files from simulation studies.
    * __timing__ - directory for storing BEAST log files and timing files for comparing computational efficiency with other Bayesian inference regimes.
    * __hiv.log__ - BEAST log file from running `./xml/hiv.xml` in BEAST.
    * __mammals.log__ - BEAST log file from running `./xml/mammals.xml` in BEAST.
    * __prokaryotes.log__ - BEAST log file from running `./xml/prokaryotes.xml` in BEAST. The log file present in this repository has been subsampled from the original to meet GitHub file size limits.
- __scripts__
    * __plots__
        * __HIV_MSE.R__ - Generates Figure 7.
        * __bacteriaCorrelation.csv__ - Stores values for generating Figure 5. Created by running Julia script `./scripts/plots/bacteria_correlation.jl`.
        * __bacteriaCorrelationLabels.csv__ - Stores the order of the traits in Figure 5. Created by running Julia script `./scripts/plots/bacteria_correlation.jl`.
        * __bacteria_correlation.jl__ - Processes `./logs/prokatyrotes.log` and prepares csv files for Figure 5.
        * __bacteria_plot.r__ - Generates Figure 5.
        * __correlation_plot.r__ - Functions for making Figures 3 and 5.
        * __correlation_prep.jl__ - Functions for processing BEAST log files into correlation csv files.
        * __custom_boxplot.R__ - Functions for customized box plots used by `./scripts/plots/HIV_MSE.R`. Modified from [here](https://github.com/tidyverse/ggplot2/blob/master/R/stat-boxplot.r).
        * __hiv_prediction.csv__ - Stored values for producing Figure 7. Created by `/scripts/hiv_prediction_analysis.jl`.
        * __mammalsCorrelation.csv__ - Stored values for producint Figure 3. Created by `./scripts/plots/mammals_correlation.jl`.
        * __mammalsCorrelationLabels.csv__ - Order of the traits for Figure 3. Created by `./scripts/plots/mammals_correlation.jl`.
        * __mammals_correlation.jl__ - Processes `./logs/mammals.log` and prepares csv files for Figure 3.
        * __mammals_plot.r__ - Generates Figure 3.
        * __simulation_plotting.r__ - Generates Figure 2, SI Figures 1, 2, and 3.
    * __storage__ - Directory for storing intermediate files produced by various scripts.
    * __dependencies.jl__ - Installs all Julia packages necessary to run the other scripts.
    * __dependencies.r__ - Installs all R packages necessary to run the other scripts. Note that on Windows machines you will also need to install [RTools](https://cran.r-project.org/bin/windows/Rtools/).
    - __hiv_pcmTiming.r__ - Times how quickly PCMBaseCpp evaluates the likelihood for the HIV example.
    - __hiv_prediction_analysis.jl__ - Analyzes results of XML created by __hiv_prediction_setup.jl__ and outputs a summary to `./scripts/plots/hiv_prediciton.csv`.
    - __hiv_prediction_setup.jl__ - Creates XML to test posterior predictive performance of various models.
    - __mammals_data.jl__ - Inputs the raw data files for the mammals example and processes them into formats used for XML files and timing.
    - __mammals_pcmTiming.r__ - Times how quickly PCMBaseCpp evaluates the likelihood for the mammals example.
    - __pcm_test_analysis.jl__ - Analyzes results of speed tests and outputs them to `./scripts/storage/speed_results.csv`.
    - __pcm_test_setup.jl__ - Creates XML and R readable data files for speed comparison between BEAST and PCMBaseCpp
    - __pcm_timing.r__ - Function for timing how quickly PCMBaseCpp evaluates likelihoods.
    - __prok_pcmTiming.r__ - Times how quickly PCMBaseCpp evaluates the likelihood for the prokaryotes example.
    - __run_timing_sim.sh__ - Linux shell script to evaluate speed of BEAST and PCMBaseCpp likelihood calculations on simulated data sets.
    - __run_timing.sh__ - Linux shell scripts to evaluate speed of BEAST and PCMBaseCpp likelihood calculations on real data sets.
    - __sim_pcmTiming.r__ - R script for timing PCMBaseCpp with simulated data.
    - __simulation_analysis.jl__ - Analyzes log files from simulation study. The log files should be in the `./logs/simulation_study` directory. Produces the __trait_simulation.csv__ and __matrix_simulation.csv__ files in the `./scripts/storage` directory.
    - __simulation_setup.jl__ - Creates XML files form simulation study.
    - __timing_analysis.jl__ - Processes results of speed comparison against earlier methods and outputs `./scripts/storage/timing_summary.csv`.
    - __timing_setup.jl__ - Creates BEAST xml files for speed tests for Table 1.
    - __xml_setup.jl__ - Creates __mammals.xml__ and __hiv.xml__ files.
- __xml__
    * __PCMBase_comparison__ - Directory for storing BEAST xml files for comparision with [PCMBaseCpp](https://github.com/venelin/PCMBaseCpp).
    * __hiv_precision__ - Directory for storing BEAST xml files for the posterior predictive performance analysis on the HIV data set. All internal files are templates used by `./scripts/hiv_prediction_setup.jl` for generating new xml and are not intdended to be run directly.
    * __simulation_study__ - Directory for storing BEAST xml files for the simulation study.
    * __timing__ - Directory for storing xml files for comparing the computational efficiency of our methods agains other Bayesian inference regimes. All internal files are templates used by `./scripts/timing_setup.jl` for generating new xml and are not intdended to be run directly.
    * __hiv.xml__ - BEAST xml file for the HIV heritability analysis in Section 7.3. Generated by `./scripts/xml_setup.jl` Julia script.
    * __mammals.xml__ - BEAST xml file for the mammals life history analysis in Section 7.1. Generated by `./scripts/xml_setup.jl` Julia script.
    * __prokaryotes.xml__ - BEAST xml file for prokaryotes analysis in Section 7.2. Created with [BEAUti](http://beast.community/beauti), the graphical user interface for generating BEAST XML files, and manually edited.

## Computing Environment

### Java

Java is required to run BEAST.
All code was compiled and run using Java 1.8. To check if you have Java installed, enter `java -version` into the command line.
If Java is installed, this command will return something like `java version "1.8.0_241"`.
If Java is not installed, it will likely return an error.

If Java is not installed, please follow the instructions [here](https://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html).
Note that if a version of Java other than 1.8 is installed, you will likely be able to run all analyses.
If you encounter problems running BEAST with another version, download Java 1.8 as a troubleshooting step.
Alternatively, follow the insctructions __Alternative BEAST Setup__ at the end of this document.

### BEAST

There are two ways to run an xml file in BEAST

1. BEAST Graphical User Interface (GUI)
    * From the command line, navigate to the `./beast` directory.
    * Enter into the command line `java -jar beast.jar`. The BEAST GUI should launch.
    * Click "Choose File..." next to the "BEAST XML File:" option, and choose the BEAST xml file you want to run.
    * At the bottom of the GUI, click "Run".
    * BEAST should start running on the command terminal.
2. Command line
    * From the command line, navigate to the `./beast` directory.
    * Enter into the command line `java -jar beast.jar <path/to/your/xml.xml>`. For example, if you wanted to run the __hiv.xml__ file in the located at `./xml/hiv.xml`, you would enter `java -jar beast.jar ../xml/hiv.xml` on Linux or OSX, and `java -jar beast.jar ..\xml\hiv.xml` on Windows.
    * BEAST should start running on the command terminal.

Regardless of which method you use to run BEAST, any log files will be located in the directory from which you run BEAST (in this case the `./beast` directory).
These log files can be investigated using [Tracer](http://beast.community/tracer).

Note that to run the `./xml/prokaryotes.xml` file, you will have to install [BEAGLE](https://github.com/beagle-dev/beagle-lib).
Instructions for installing BEAGLE with BEAST can be found [here](https://beast.community/beagle). BEAGLE is not required to run any other BEAST xml files.



### R
    
Some scripts (those ending in `.r` or `.R`) related to plotting and the comparison with [PCMBase](https://venelin.github.io/PCMBase/) are written in R. All scripts were run with R-3.6.2, but will likely work on any R version 3.4 or above.
Instructions for installing R can be found [here](https://www.r-project.org/).

If you use a Windows machine, you will also need to install [RTools](https://cran.r-project.org/bin/windows/Rtools/) to install all relevant packages.

__To ensure all necessary packages are installed, run the `./scripts/dependencies.r` R script.__

__All R scripts in this repository should be run from the directory they are located in__.
To run an R script, you can use [RStudio](https://rstudio.com/), or simply navigate on the comand line to the appropriate directory and enter `Rscript <file.r>` into the command line, replacing `<file.r>` with the desired script.

<!--

2. Install [R](https://www.r-project.org/) (and Rscript).
    * All scripts were tested on R-3.6.2
3. Install [RTools](https://cran.r-project.org/bin/windows/Rtools/).
4. Install relevant R packages by running  __dependencies.r__ script (see note under __Work Flows__ for instructions on how to run R scripts).

-->

### Julia

Most scripts (those ending in `.jl`) related to data pre-processing, simulating data, and processing BEAST log files are written in [Julia](https://julialang.org/).
All scripts were run with Julia v1.4, but will likely work on any version 1.0 or above.
Instructions for installing Julia can be found [here](https://julialang.org/downloads/).

__To ensure all necessary packages are installed, run the `./scripts/dependencies.jl` Julia script.__

To run a Julia script, use the command line to navigate to the directory containing the script of interest and enter `julia <file.jl>` into the command line, replacing `<file.jl>` with the script you want to run.







<!--

## Setting up BEAST

1. Ensure you have Java version 1.8 installed.
    * From the command line, type `java -version`.
    * The output should start with `java version "1.8.0_<other numbers>"`.
    * If java is not installed or the wrong version of java is installed, please follow the instructions [here](https://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html).
    * Note that if a different version of Java is already installed, you will probably be able to run all of the xml. If you encounter problems, you can update your version of Java later.
2. Ensure you have `git` installed.
    * From the command line, type `git --version`. The version of git should be returned.
    * If `git` is not installed, you can find directions to install it [here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).
3. Ensure you have `ant` installed.
    * From the command line, type `ant -version`. The version of ant should be returned.
    * If `ant` is not installed, you can find directions to install in [here](https://ant.apache.org/manual/install.html).
5. Download and build BEAGLE (this is only necessary to run the `prokaryote.xml` example).
    * Follow the instructions [here](https://github.com/beagle-dev/beagle-lib) to download and install BEAGLE.
6. Download and build BEAST
    * Enter the following code on the command line. 
    ```
    git clone https://github.com/beast-dev/beast-mcmc.git
    cd beast-mcmc
    git checkout repeated_measures
    ant
    ```
    
## Running an XML file on BEAST 
1. Navigate to the directory with the `beast.jar` file, which can be found in `./beast` directory.
    * From the `beast-mcmc` directory (which you set up above), enter the following into the command line:
    ```
    cd build
    cd dist
    ```
2. Run an `.xml` file on BEAST
    * Option 1: launch the BEAST gui
        * Enter `java -jar beast.jar`
        * A gui will pop up.
        * Click __Choose File__ and select the `.xml` file you wish to run.
        * Click __Run__.
        * BEAST should begin running in a gui window.
    * Option 2: run BEAST from the command line
        * Enter `java -jar beast.jar <path/to/your/file.xml>`.
        * BEAST should begin running on the command line.
        
    * The output `.log` file will be in the `././beast` directory.
3. Use Tracer to view the `.log` file.
    * Tracer installation and usage instructions can be found [here](http://beast.community/tracer).
    



## Scripts
- __dependencies.jl__ - Installs all Julia packages necessary to run the other scripts.
_ __dependencies.r__ - Installs all R packages necessary to run the other scripts. Note that on Windows machines you will also need to install [RTools](https://cran.r-project.org/bin/windows/Rtools/).
- __hiv_pcmTiming.r__ - Times how quickly PCMBaseCpp evaluates the likelihood for the HIV example.
- __mammals_data.jl__ - Inputs the raw data files for the mammals example and processes them into formats used for XML files and timing.
- __hiv_prediction_analysis.jl__ - Analyzes results of XML created by __hiv_prediction_setup.jl__ and outputs a summary to `./scripts/plots/hiv_prediciton.csv`.
- __hiv_prediction_setup.jl__ - Creates XML to test posterior predictive performance of various models.
- __mammals_pcmTiming.r__ - Times how quickly PCMBaseCpp evaluates the likelihood for the mammals example.
- __pcm_test_analysis.jl__ - Analyzes results of speed tests and outputs them to `./scripts/storage/speed_results.csv`.
- __pcm_test_setup.jl__ - Creates XML and R readable data files for speed comparison between BEAST and PCMBaseCpp
- __pcm_timing.r__ - Function for timing how quickly PCMBaseCpp evaluates likelihoods.
- __prok_pcmTiming.r__ - Times how quickly PCMBaseCpp evaluates the likelihood for the prokaryotes example.
- __run_timing_sim.sh__ and __run_timing.sh__ - Linux shell scripts to evaluate speed of BEAST and PCMBaseCpp likelihood calculations on simulated data sets. Note that the directories are hard-coded and you will have to manually modify these files to fit your own file structure. Specifically, `repo_dir` should be the top directory of this git repo, and `beast_jar` should be set to the location of your BEAST .jar file.
- __sim_pcmTiming.r__ - R script for timing PCMBaseCpp with simulated data.
- __simulation_analysis.jl__ - Analyzes log files from simulation study. The log files should be in the `./logs/simulation_study` directory. Produces the __trait_simulation.csv__ and __matrix_simulation.csv__ files in the `./scripts/storage` directory.
- __simulation_setup.jl__ - Creates XML files form simulation study.
- __timing_analysis.jl__ - Processes results of speed comparison against earlier methods and outputs `./scripts/storage/timing_summary.csv`.
- __timing_setup.jl__ - Creates XML files for speed tests for Table 1.
- __xml_setup.jl__ - Creates __mammals.xml__ and __hiv.xml__ files.

### Plotting Scripts

- __bacteria_plot.r__ - Generates Figure ?.
- __HIV_MSE.R__ - Generates Figure 7.

-->

## Work Flows

Unless otherwise noted, all scripts are in the `./scripts` directory of this repository.

<!--
Unless otherwise noted, all command line arguments should be run from the `./scripts` directory of this repository.
All scripts are written in either R (file ending with `.r` or `.R`) for Julia (files ending with `.jl`).
To run an R script, simpyly enter `Rscript <script_name.r>` into the command line, replacing `<script_name.r>` with the name of the relevant script.
For example to run the __dependencies.r__ script, you would enter `Rscript dependencies.r` from the `./scripts` directory.
To run a Julia script, simply enter `julia <script_name.jl>` into the command line, again replacing `<script_name.jl>` with the name of the relevant Julia script.

All references to the `beast-mcmc` directory point to the location of the BEAST installation from the __Setting up BEAST__ section above. All other file or directory locations reference this repository.


### Computing Environment
1. Set up Java and BEAST as described above.
2. Install [R](https://www.r-project.org/) (and Rscript).
    * All scripts were tested on R-3.6.2
3. Install [RTools](https://cran.r-project.org/bin/windows/Rtools/).
4. Install relevant R packages by running  __dependencies.r__ script (see note under __Work Flows__ for instructions on how to run R scripts).
5. Install [Julia](https://julialang.org/).
    * All scripts were tested on Julia v1.4.0, but will likely work for any v1.x
6. Install relevant Julia packages by running __dependencies.jl__ script.
-->

### Mammals Correlation and Figure 3
1. Run __mammals_data.jl__ script.
2. Run __xml_setup.jl__ script.
    * Creates __mammals.xml__ file in the `./xml` directory.
3. Run the `./xml/mammals.xml` file according to the instructions in the __BEAST__ section above.
4. Move the log file, located at `./beast/mammals.log`, to the `./logs` directory of this repository.
    * You can examine this file using [Tracer](http://beast.community/tracer).
5. Run the __mammals_correlation.jl__ script in the `./scripts/plots` directory to extract relevant information from log file.
6. Run the __mammals_plot.r__ script in the `./scripts/plots` directory to generate Figure 3.

### Prokaryotes Correlation and Figure 5
1. Run the `./xml/prokatyotes.xml` file according to the instructions in the __BEAST__ section above.
    * The __prokaryotes.xml__ file in the `./xml` directory was created by [BEAUti](https://beast.community/beauti) and manually edited.
2. Move the __prokaroytes.log__ file (located in the `./beast` directory) to the `./logs` directory.
    * You can examine this file using [Tracer](http://beast.community/tracer).
3. Run the __bacteria_correlation.jl__ script in the `./scripts/plots` directory to extract relevant information from log file.
4. Run the __bacteria_plot.r__ script in the `./scripts/plots` directory to generate Figure 5.

### HIV Heritability
1. Run __xml_setup.jl__ script.
    * Creates __hiv.xml__ file in the `./xml` directory.
2. Run the __hiv.xml__ file in the `./xml` directory according to the instructions in the __BEAST__ section above.
3. Move the __hiv.log__ file (located in the `./beast` directory) to the `./logs` directory.
4. You can examine this file using [Tracer](http://beast.community/tracer). The relevant parameters are:
    * varianceProportionStatistic11 $\rightarrow$ GSVL
    * varianceProportionStatistic22 $\rightarrow$ SPVL
    * varianceProportionStatistic33 $\rightarrow$ CD4 slope

### HIV Posterior Predictive Power and Figure 7
1. Run __hiv_prediction_setup.jl__ script to generate xml files.
    * All xml are output to the `./xml/hiv_prediction` directory.
2. Run all xml files ending in a number in the `./xml/hiv_prediction` directory according to the instructions in the __BEAST__ section above.
    * Note that the four xml files not ending in a number are templates and have not had any data removed.
3. Copy all log files to the `./logs/hiv_prediction` directory.
4. Run the __hiv_prediction_analysis.jl__ script.
5. Run the __HIV_MSE.R__ script in the `./scripts/plots` directory to generate Figure 7.

### Computational Efficiency and Table 1
1. Run __timing_setup.jl__ script to generate xml for timing.
    * All xml are in the `./xml/timing` directory.
2. Run all xml files ending in a number in the `./xml/timing` directory according to the instructions in the __BEAST__ section above.
    * Note that you will need to save the screen output to a `.txt` file with the same name as the original xml file. To do this, using the command line, use the following command `java -jar beast.jar <path/to/your/file.xml> > <file.txt>`, replacing `<path/to/your/file.xml>` and `<file.txt>` with the relevant paths and filenames.
3. Move all `.log` files and saved `.txt` files generated to the `./logs/timing` directory of this repository.
4. Run the __timing_analysis.jl__ script. The output will be written to `./scripts/storage/timing_summary.csv` file.

### Simulation Study; Figure 2; and SI Figures 1, 2, and 3
1. Ensure the __mammals.log__, __hiv.log__, and __prokaryotes.log__ files are in the `./logs` directory of this repository.
2. Run the __simulation_setup.jl__ script to generate simulated data sets and xml files.
3. Run all xml files ending in the `./xml/simulation_study` directory according to the instructions in the __BEAST__ section above.
4. Copy all `.log` files produced from running the xml files above to the `./logs/simulation_study` directory of this repository.
5. Run the __simulation_analysis.jl__ script.
6. Run the __simulation_plotting.r__ script in the `./scripts/plots` directory to generate Figure 2 and SI Figures 1, 2, and 3.

### Comparison with PCMBaseCpp and SI Tables 1 and 2
1. Ensure the __mammals.log__, __hiv.log__, and __prokaryotes.log__ files are in the `./logs` directory of this repository.
2. Run the __mammals_data.jl__ script to generate the relevant newick files.
3. Run the __pcm_test_setup.jl__ script to generate simulated data sets and xml files.
4. Execute the __run_timing.sh__ and __run_timing_sim.sh__ linux shell scripts.
    * The __run_timing_sim.sh__ script will probably take a long time (several hours) to complete.
5. Run the __pcm_test_analysis.jl__ script. The results are located in the `./scripts/storage/speed_results.csv` file.

