@echo off
setlocal enabledelayedexpansion

set script_dir=D:\missing_traits\paper\scripts
set beast_jar=D:\beast-mcmc\build\dist\beast.jar
set xml_dir=D:\missing_traits_paper\xml\PCMBase_comparison
set time_dir=D:\missing_traits_paper\logs\PCMBase_timing

echo %xml_dir%

for /L %%i in (1, 1, %1) do (
	echo "Starting for loop"
	echo Test %%i

	for %%j in (hiv, prok, mammals) do (
		echo %%j

		java -jar %beast_jar% -overwrite %xml_dir%\%%j%PCMComparison.xml > nul 2>&1
		Rscript %%j%_pcmTiming.r > nul 2>&1

		set old_name_pcm=%%j%PCMComparisonTimer_pcm.txt
		set new_name_pcm=%%jPCMComparisonTimer_pcm_%%i.txt
		set old_name_beast=%%j%PCMComparisonTimer_beast.txt
		set new_name_beast=%%jPCMComparisonTimer_beast_%%i.txt


		ren !old_name_pcm! !new_name_pcm!
		ren !old_name_beast! !new_name_beast!
		
		move !new_name_pcm! !time_dir!\
		move !new_name_beast! !time_dir!\



	)

	echo =============================
)
