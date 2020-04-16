#!/bin/bash

repo_dir=/media/ghassler/SD_Storage/missing_traits_paper/
script_dir=$repo_dir/scripts/
xml_dir=$repo_dir/xml/PCMBase_comparison
time_dir=$repo_dir/logs/PCMBase_timing
storage_dir=$script_dir/storage/PCMBase_comparison

beast_jar=/media/ghassler/SD_Storage/beast-mcmc/build/dist/beast.jar

cd $script_dir

for i in $(cat sim_timing_files.txt); do
 
    echo "Starting $i"
    echo "Beast"
    java -jar $beast_jar -overwrite $xml_dir/$i.xml #> /dev/null 2>&1
    echo "PCMBase"
    Rscript sim_pcmTiming.r $i $storage_dir # > /dev/null 2>&1
#    echo "PhylogeneticEM"
#    Rscript "$j"_phyloEMTiming.r #> /dev/null

    mv "$i"Timer_beast.txt $time_dir
    mv "$i"Timer_pcm.txt $time_dir
#    mv "$j"Timer_phy.txt $time_dir
    echo "Done"
done
        
