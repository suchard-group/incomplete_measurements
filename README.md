# missing_traits_paper

## Setting up BEAST (Windows, Mac, and Linux)
1. Ensure you have Java version 1.8 installed.
    * From the command line, type `java -version`.
    * The output should start with `java version "1.8.0_<other numbers>"`.
    * If java is not installed or the wrong version of java is installed, please follow the instructions [here](https://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html).
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
    git clone --branch repeated_measures --depth 1 https://github.com/beast-dev/beast-mcmc.git
    cd beast-mcmc
    ant
    ```
## Running an XML file on BEAST 
1. Navigate to the directory with the `beast.jar` file.
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
        
    * The output `.log` file will be in the `./beast-mcmc/build/dist` directory.
3. Use Tracer to view the `.log` file.
    * Tracer installation and usage instructions can be found [here](http://beast.community/tracer).
        
        
    
