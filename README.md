# PYTEST TESTING TOOL (PTT)

## What is a Pytest

Simply put, a Pytest is a pass or fail Python test. For instance, with the WCS step, we 
have Python scripts (which we are calling auxiliary code within the frame of the testing 
tool) that compare the pipeline product with the ESA corresponding intermediary file, and 
calculates a difference. The Pytest is to asses if that difference is less than or equal 
to an X threshold value. Hence, a failed test means that the condition was not met. If an
error should occur with the pipeline, the test will be flagged as an error.



## Quick Start Guide


NOTE: This guide assumes that Conda has been installed. If you have not yet done 
so, please follow the instructions at:
https://astroconda.readthedocs.io/en/latest/
Please use python 3.5


1. Create the conda environment for testing and get the configuration files.  

a. Conda environment for this testing campaign:
- Current testing version is 7.1 for release candidate 0.7.9.0. However, this version has
a broken version of the assign WCS step. Hence, we are currently testing with the 
development version of the pipeline (please see the NOTE below).
- In case you need to test some specific step of release candidate 0.7.9.0, please
follow these instructions to create the build7.1 environment. In a terminal type:
```bash
conda create -n jwst_b7.1 --file url_depending_on_your_system python=3.5
```
for the current release candidate, the ulr options are:
- Linux: http://ssb.stsci.edu/releases/jwstdp/0.9.x/latest-linux
- OS X: http://ssb.stsci.edu/releases/jwstdp/0.9.x/latest-osx


NOTE:
If you need to use the development version of the pipeline then do the following:
```bash
conda create -n jwst_dev -c http://ssb.stsci.edu/conda-dev jwst python=3.5
```
Then, to update the development environment, activate the environment and then type:
```bash
conda update --override-channels -c http://ssb.stsci.edu/conda-dev -c defaults --all
```

b. Configuration files corresponding to this build. Create a directory (e.g. 
```b7.1cfg_files```) somewhere in your testing working space, and ```cd``` into it. Now 
type the following command within the conda environment you just created (see step 2).
```bash
collect_pipeline_cfgs .
```


2. Activate the conda environment for testing the pipeline, e.g. type:
```bash
source activate jwst_b7.1
```
From here on, every step of this guide should happen within the conda testig environment.

NOTE: 
- If you forget what did you name your new environment type:
```bash
conda env list
```
this will list all environments you have created.
- If you want to remove an old testing environment type:
```bash
conda env remove -n name_of_your_old_environment
```


3. Install the pytest html plugin. Within the conda testing environment, type:
```bash
pip install pytest-html
```
NOTE: Every time you create a new conda environment you need to install html plugin.


4. Clone the repository so you have the PTT. To do this click at the top right 
of this page, in the dropdown button that says clone or download, then copy the ulr that
appears there. Now, within the conda testing environment, go to or create the directory 
where you want the PTT to "live" in. However, make sure that the configuration files 
directory is at least at the top level of the directory tree where the PTT will live, e.g. 
the ```b7.1cfg_files``` directory and the ```nirspec_testing_tool``` directory can be at 
the same level, but the ```b7.1cfg_files``` directory cannot be inside the 
```nirspec_testing_tool``` directory because otherwise the .cfg files
will be picked up by Git and it will try to put them in your space.
```bash
git clone the_ulr_you_copied
```
After this is done, you should see a full copy of the PTT in your directory.

NOTE: 
- If you have already cloned the repository, in the terminal go to where you placed 
the ```nirspec_pipe_testing_tool``` directory. Then, use the following command to update 
the code:
```git pull```
- If, however, you had written script(s) in the tool's directory tree, git will not let 
you continue until you move the script(s) to another directory. 


5. Prepare the data to run through the pipeline. To do this:

a. Copy the test data you will use from the NIRSpec vault directory. Go to the directory 
where you want to have the testing data, and from there type:
```bash
cp /grp/jwst/wit4/nirspec_vault/prelaunch_data/testing_sets/b7.1_pipeline_testing/
                                             test_data_suite/the_data_you_want_to_copy .
```

b. In the directory where you copied the test data, you will need to run a script PER
fits file you want to test. Do not worry, this will only be done once. This script will
create a new subdirectory with the necessary input file to run the SSB script 
that converts raw data into uncal type files. You can choose to either keep this 
subdirectory, or tell the script to remove it after the operation is done. In the 
terminal type:
```bash
python /path_to_the_testing_tool/nirspec_pipe_testing_tool/utils/
                                        prepare_data2run.py fits_file.fits MODE -u
```
where the MODE is expected to be one of: FS, MOS, IFU. This command will update the uncal 
keyword header without creating a new file, and will also keep the subdirectory. To 
remove it, simply add ```-rm``` at the end. To save the keyword changes in a new fits 
file (instead of updating), remove the ```-u```. The new uncal fits file is now ready 
for pipeline ingest.

c. Optional. If you want to see the header of any file, you can use the another script
in the ```utils``` directory of the PTT. If you just want to see on-screen the 
header, go where your fits file "lives" and type:
```bash
python /path_to_the_testing_tool/nirspec_pipe_testing_tool/utils/
                                                        read_hdr.py fits_file.fits
```
This command will show the main header. To save the header to a text file add a ```-s``` 
at the end. If you want to see/save a different extension add at the end ```-e=1``` for 
extension 1, and so on.

d. Now, the data is ready to be ran through cal_detector1. Please go ahead with step 6
of this guide to do that.


6. Set the PTT configuration file. In a text file editor, you are going to modify the 
configuration file named ```cwspec2_config.cfg```, which lives at 
```/nirspec_pipetesting_tool/calwebb_spec2_pytests/```. 
This is the file that controls all the input that the tool needs. Please open it and 
make sure that:
- All the paths point to the right places. The files can be located anywhere, but both,
the pipeline and the tool, will run faster if the files are local on your computer.
- The input file for the PTT is the final output file from calwebb_detector1.
- The adequate mode for the data to be tested is set correctly, choices are: FS, IFU,
or MOS.
- The steps that you want to be ran or not are set to True or False.
- In the bottom part of the file, all the additional arguments for the PTT are 
correct, e.g. threshold values, figure switches, and additional fits files.


7. Run the ```calwebb_detector1``` pipeline. The final output of this is the level 2 data
required to run the PTT. In a terminal, go into the directory where the testing tool lives 
(i.e. at the level of the ```calwebb_spec2_pytests``` directory), and make sure that the 
testing conda environment is on. Now type:
```bash
python ../utils/run_cal_detector1.py /path_where_the_uncal_file_lives/uncal_file.fits -sbs
```
This script will execute the calwebb detector 1 pipeline step by step. If you want to run
it in a single run, using the configuration file that you have for it in the ```utils``` 
directory, simply remove the ```-sbs``` from previous command. 
If everything went well, you will see a text file called 
```cal_detector1_outputs_and_times.txt```, which contains the steps ran, the name of the 
output fits file, and the time each step took to run. However, if you chose to run the 
calwebb detector1 in a single run, you will only see the total running time. This text file,
along with the intermediary products will be located in the path you set for the 
```working_directory``` variable in the configuration file of the PTT.

You now are able to run the MESA calwebb 
detector 1 testing tool. Steps to obtain the MESA testing tool:

a. At the same level as the top directory of the PTT (i.e. the ```nirspec_pipe_testing_tool``` 
directory), create a new directory called ```MESA_cal_detector1```.

b. Inside ```MESA_cal_detector1```, you will clone their Git repository. Please follow the 
directions for this at:
http://calibration-pipeline-testing-tool.readthedocs.io/en/latest/

c. Now, in the ```utils``` directory of ```nirspec_pipe_testing_tool```, you will find a 
sample json file that you can modify in order to use as input for the MESA calwebb
detector 1 testing tool. Notice that it differs from the example given in the MESA 
documentation. In our example, the steps are in the order of the pipeline, to make it 
easier to determine what file is input to which step. Please copy our sample json file 
from the ```utils``` directory into the ```MESA_cal_detector1``` directory you created,
and or modify the json file with the name and path of your intermediary fits products 
from running the calwebb_detector1 pipeline.

d. Run the MESA testing tool with the following command:
```bash
test_pipeline --config ./cal_detector1_input.json
```
Please record your progress. Make a PDF export of the html report, and place it in 
the corresponding column of Table 2 of our testing campaign Confluence page:
https://confluence.stsci.edu/display/JWST/NIRSpec+Pipeline+Testing+Build+7.1+part+1

 
8. Ready to run PTT. Go back to the directory where PTT lives and into the 
```calwebb_spec2_pytests``` directory, copy final output file from calwebb detector1 into 
the working directory you indicated in the ```cwspec2_config.cfg``` file, and make sure 
that the input file for the PTT matches the file you just copied into the working 
directory. Now, to ensure that everything is in order, and to see what pytests will be 
executed and in which order type:
```bash
pytest --collect-only
```

9. Do the first PTT run. As an output of the testing tool you will see an html 
file, ```report.html```, and an intermediary product file name map will appear in the 
```calwebb_spec2_pytests``` directory. The fits files of intermediary products will 
be saved in the path you indicated at the ```cwspec2_config.cfg``` file with the variable
 ```working_directory```. In the terminal type:
```bash
pytest -s --config_file=cwspec2_config.cfg --html=report.html
```
The ```-s``` will capture all the print statements in the code on screen.


10. Report your findings. If all went well, you should have the html report in the 
```calwebb_spec2_pytests``` directory, along with the file name map of intermediary 
products that appeared in the working directory. Each section of the report specifies 
what tests passed or failed, and, if there were any errors, the report will tell you 
(with the name of the script where the error occurred) if it is a pipeline error or a 
PTT error. Keep updating your testing progress in our testing campaign Confluence page:
https://confluence.stsci.edu/display/JWST/NIRSpec+Pipeline+Testing+Build+7.1+part+1

Please follow these actions for reporting your progress, depending on the outcome:
- Do not place the html file in the Confluence page because it will get corrupted. 
Instead, please create a PDF export from the report and use that file for the 
Confluence page and for any sharing purposes.
- If there are no testing tool errors and some of the pytests failed, you can check  
off the steps and please put that report in Table 2 of the Confluence page.
- If there are testing tool errors please place your report in the Confluence page, 
add a small comment in the "Testing Tool" column of Table 2 in the Confluence page, and 
send an email to Maria Pena-Guerrero (pena@stsci.edu), so that the 
PTT error can be addressed as soon as possible.
- If there are errors with the pipeline, please create a Basecamp message 
(https://stsci-ins.basecamphq.com/projects/11477312-jwst-pipeline/posts), following the 
NIRSpec guidelines for doing so, which you can find at:
https://confluence.stsci.edu/display/JWST/NIRSpec+Guidelines+for+Pipeline+Issue+Reporting
and copy the link in the corresponding column of Table 2 of the Confluence page.
- Keep updating the html report as you continue with the testing, and including the final
report.


11. Put final results in the NIRSpec vault. At the end of the testing campaign, please do 
the following within your working directory: 

a. Create a new directory called ```yourname_MODE_Name_of_the_raw_data_file_results```, 
e.g. ```maria_IFU_NRSSMOS-MOD-G1M-17-5344175105_1_491_SE_2015-12-10T18h00m06_results```.
(If you already copied a previous version of your results and you want to keep it then add
a ```_v2``` -or the version number that corresponds- after the word results, e.g. 
```maria_IFU_NRSSMOS-MOD-G1M-17-5344175105_1_491_SE_2015-12-10T18h00m06_results_v2```) 

b. Inside that new direcoty place all the intermediary fits products as well as the html
report

c. Also inside that directory, create a text file named ```YourNameMODEresults.txt```, 
e.g. ```mariaIFU_NRSSMOS-MOD-G1M-17-5344175105_1_491_SE_2015-12-10T18h00m06_results.txt```. 
In this text file you will only type the full path where you obtained the testing data, 
e.g. 
```
/grp/jwst/wit4/nirspec_vault/prelaunch_data/testing_sets/b7.1_pipeline_testing/
                                                              test_data_suite/IFU_CV3
```

d. Finally, place the results directory you created in the staging directory of the 
NIRSpec vault, and send Gray an email (gkanarek@stsci.edu). The path of the staging 
directory is:
```/grp/jwst/wit4/nirspec_vault/staging```
	


## TO KEEP IN MIND

- A text file containing an intermediary product name map will be created in the pytests 
directory.
- If any of the central store directory calls do not respond (e.g. when looking at the 
flats), the pytest will be skipped even if the step is set to True in the config file. 
The failing message will say that the step was set to False (this is a known bug). To 
force the tests run, you will have to download the files the tool is calling, and change 
the corresponding paths in the configuration file.
- The output in the terminal can be a bit overwhelming if there was a failed test or an 
error, since it shows both, the pipeline messages and the PTT messages. In the html 
report is much clearer to understand what happened.
- As part of the testing campaign, it is important that you run the pipeline from the 
command line as well, and that you make sure that the outcome intermediary files are 
consistent with the ones ran with scripts, i.e. the PTT. This sanity check is minor 
but important to verify. 
In the ```utils``` directory there are two text files named 
```terminal_commands_calwebb_detector1_steps.txt``` and
```terminal_commands_calwebb_spec2_steps.txt```. These files contain all the commands 
you can use from the terminal to run the calwebb_detector1 and calwebb spec2 steps,
respectively.
- Finally, remember that:

a. Whenever you need to read either the main or science headers of a file,
you can always use the ```read_hdr.py``` script located in the ```utils``` directory of
the testing tool.

b. If you need to change/add a keyword value to a specific extension of a file, you can
use the ```change_keywd.py``` script, also located at the ```utils``` directory.



## POSSIBLE OUTCOMES OF THE PYTESTS

- Passed = the assertion was true, so the test condition was met.
- Failed = the assertion was false, the test condition was NOT met (there will be an 
AssertionError on-screen and in the html file, with a clear PTT customized message of 
what happened).
- Skipped = the test was skipped (there will be a message on the html report and on-screen
explaining why the test was skipped).
- Error = this is a coding error, or a bug in either the pipeline or the PTT. Please 
follow the appropriate steps described in the "For Testers" section of our testing 
campaign Confluence page:
https://confluence.stsci.edu/display/JWST/NIRSpec+Pipeline+Testing+Build+7.1+part+1



## ADDING TESTING ROUTINES

You can add additional testing routines. We prefer if these are written in python 3.5, 
however, due to time constraints, we devised a method for including scripts in other 
programing languages. For this testing campaign, the guidelines for auxiliary code are 
the following: 
1. Programing language = preferably Python 3.5, however, other languages (e.g. IDL, C) 
are OK.

2. Please comment as much as you can so that it is easy for other testers to follow 
your code.

3. At the end of your script you must create a pass or fail test, e.g. after you calculate 
a root-mean-square check if that value is within some range. If instead you have a 
threshold value, please make the threshold an input variable with a default value.

4. Your script must produce a text file named ```scriptname_result.txt``` at
the end, which will only contain one word: True or False.

5. Create a directory named ```scriptname``` with the following contents:
a. A text file named ```scriptname_init.txt```, and in there type in line 1 
the pipeline step after which your code should be run, in line 2 the language you used, 
and in line 3 the exact command you would use to run the script. If this is an IDL script, 
line 2 should be the exact call from within the IDL session.

b. The script you wrote.

c. One level up from where you created the script directory, create a new directory named 
```YourName_PTT_auxiliary_code```. In this directory, please place the directory 
```scriptname```, and a create a text file named 
```YourName_PTT_auxiliary_code.txt```. 
This text file will only contain the following path:
```/grp/jwst/wit4/nirspec_vault/pipe_testing_tool/auxiliary_code```

6. Now, copy the ```YourName_PTT_auxiliary_code``` directory in the NIRSpec vault, and 
please send an email to Gray Kanarek (gkanarek@stsci.edu). The vault's path is:
```/grp/jwst/wit4/nirspec_vault/staging```.

7. Your code will be implemented for the next testing campaign and if there is a threshold 
associated with it, it will appear as a variable in the PTT configuration file. For the 
current testing campaign, please place it in the "PTT Comments and/or errors" column of  
Table 2 of our progress keeping Confluence page:
https://confluence.stsci.edu/display/JWST/NIRSpec+Pipeline+Testing+Build+7.1+part+1
You will be asked to write a description of your code for the final testing report 
of the current testing campaign.



## Enjoy your pipeline testing!

If you have any question of what a specific step does, you can get a description at:
http://jwst-pipeline.readthedocs.io/en/latest/jwst/pipeline/description.html#stage-2-spectroscopic-pipeline-step-flow-calwebb-spec2
