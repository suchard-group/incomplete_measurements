@echo off
setlocal enabledelayedexpansion

set script_dir=D:\missing_traits_paper\scripts
set beast_jar=D:\beast-mcmc\build\dist\beast.jar
set xml_dir=D:\missing_traits_paper\xml\PCMBase_comparison
set time_dir=D:\missing_traits_paper\logs\PCMBase_timing
set storage_dir=%script_dir%\storage\PCMBase_comparison

cd %script_dir%



echo %storage_dir%

for /F %%i in (sim_timing_files.txt) do (

	echo %%i
	
	java -jar %beast_jar% -overwrite %xml_dir%\%%i.xml > nul 2>&1
	Rscript sim_pcmTiming.r %%i %storage_dir% > nul 2>&1
	
	move %%iTimer_pcm.txt %time_dir%\
	move %%iTimer_beast.txt %time_dir%\

)
