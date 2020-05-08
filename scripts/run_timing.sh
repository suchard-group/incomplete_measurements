#!/bin/bash

repo_dir=..
script_dir=$repo_dir/scripts/
xml_dir=$repo_dir/xml/PCMBase_comparison
time_dir=$repo_dir/logs/PCMBase_timing

beast_jar=$repo_dir/beast/beast.jar

for ((i=1; i <= $1; ++i)); do
    for j in hiv prok mammals; do
        echo "Starting $j $i"
        echo "Beast"
        java -jar $beast_jar -overwrite $xml_dir/"$j"PCMComparison.xml #> /dev/null 2>&1
        echo "PCMBase"
        Rscript "$j"_pcmTiming.r #> /dev/null


        mv "$j"PCMComparisonTimer_beast.txt $time_dir/"$j"PCMComparisonTimer_beast_$i.txt
        mv "$j"PCMComparisonTimer_pcm.txt $time_dir/"$j"PCMComparisonTimer_pcm_$i.txt
        echo "Done"
        
    done
done
        
