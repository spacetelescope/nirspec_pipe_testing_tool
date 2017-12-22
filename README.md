# PYTEST TESTING TOOL

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


1. Create the conda environment for testing the correct version of the pipeline
* Current testing version is 7.1
* Specific instructions to create the build7.1 environment.
For example, in a terminal type:
```bash
conda create -n jwst_b7.1 --file url_depending_on_your_system python=3.5
```
for the current release candidate, the ulr options are:
- Linux: http://ssb.stsci.edu/releases/jwstdp/0.7.8/latest-linux
- O5S X: http://ssb.stsci.edu/releases/jwstdp/0.7.8/latest-osx


2. Activate the conda environment for testing the pipeline, e.g. type:
```bash
source activate jwst_b7.1
```
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


3. Install the pytest html plugin. Within the conda environment type:
```bash
pip install pytest-html
```


4. Clone the repository so you have the testing tool. To do this click at the top right 
of this page, in the dropdown button that says clone or download, then copy the ulr that
appears there. Now, go/create to the directory where you want the testing tool to "live" 
in, and type:
```bash
git clone the_ulr_you_copied
```
After this is done you should see a full copy of the testing tool in your directory.

NOTE: If you have already cloned the repository, in the terminal go to where you placed 
the ```nirspec_pipe_testing_tool``` directory. Then use this command to update the code:
```git pull```
However, if you had written script(s) in the tool's directory tree, git will not let you
until you move the script(s) to another directory. 


5. Prepare the data to run through the pipeline. To do this:
a) Copy the test data you will use from the NIRSpec vault directory. Go to the directory 
where you want to have the testing data:
```bash
cp /grp/jwst/wit4/nirspec_vault/prelaunch_data/testing_sets/b7.1_pipeline_testing/
                                             test_data_suite/the_data_you_want_to_copy .
```

b) In the directory where you copied the test data, you will need to run a script PER
fits file you want to test. Do not worry, this will only be done once. This script will
create a new subdirectory with the necessary input file to run the SSB script 
that converts raw data into uncal type files. You can choose to either keep this subfolder,
or ask to remove it after the operation is done. In the terminal type:
```bash
python /path_to_the_testing_tool/nirspec_pipe_testing_tool/utils/
                                                   prepare_data2run.py fits_file.fits -u
```
This command will update the uncal keyword header without creating a new file, and
will also keep the subdirectory. To remove it, simply add ```-rm``` at the end. To save
the keyword changes in a new fits file (instead of updating), remove the ```-u```.
The new uncal fits file is now ready for pipeline ingest.

c) Optional. If you want to see the header of any file, you can use the another script
in the ```utils``` directory of the testing tool. If you just want to see on-screen the 
header, go where your fits file "lives" and type"
```bash
python /path_to_the_testing_tool/nirspec_pipe_testing_tool/utils/
                                                              read_hdr.py fits_file.fits
```
This command will show the main header. To save the header to a text file add a ```-s``` 
at the end. If you want to see/save a different extension add at the end ```-e=1``` for 
extension 1, and so on.

d) Now, the data is ready to be ran through cal_detector1. Please go ahead on to step 6
of this guide to do that.


6. In a text file editor, you are going to modify the configuration file 
```cwspec2_config.cfg```, which lives at 
```/nirspec_pipetesting_tool/calwebb_spec2_pytests/```. 
This is the file that controls all the input that the tool needs. Please open it and 
make sure that:
- All the paths point to the right places (the files can be located anywhere, but the both
the pipeline and the tool will run faster if the files are local to your computer).
- Set the adequate mode for the data to be tested.
- Set to True or False the steps that you want to be ran or not.
- All the additional arguments for the testing tool are correct, e.g. the threshold values
and figure switches in the bottom part of the file.


7. You are now ready to run ```calwebb_detector1``` to obtain the level 2 data required to 
run the testing tool. In a terminal, go into the directory where the testing tool lives 
(i.e. at the level of the ```/calwebb_spec2/``` directory), and make sure the testing
conda environment is on. Now type:
```bash
python ../utils/runcal_detector1.py /path_where_the_uncal_file_lives/uncal_file.fits
```
This will execute calwebb detector 1 in a single run, using the configuration file that 
you have for it in the ```utils``` directory. You can also add a ```-sbs``` in order to 
execute the processing of level 1 step by step. If everything went well, you now are able
to run the calwebb detector 1 testing tool, and you should have the right input for 
running the calwebb spec 2 testing tool. 
To run the calwebb detector 1 testing tool please follow the directions at:
http://calibration-pipeline-testing-tool.readthedocs.io/en/latest/
In the ```utils``` directory you will find a sample json file that you can modify to run 
that tool.
 

8. See what pytests will be executed and in which order. Within the active conda 
environment type:
```bash
pytest --collect-only
```

9. Do the first testing tool run. As an output of the testing tool you will see an html 
file, ```report.html```, and an intermediary product file name map will appear in the 
```/calwebb_spec2/``` directory, along with fits files of intermediary products in the 
path you pointed with the ```cwspec2_config.cfg``` file with the variable
 ```working_directory```. In the terminal type:
```bash
pytest -s --config_file=cwspec2_config.cfg --html=report.html
```
The ```-s``` will capture all the print statements in the code on screen.


10. Report your findings. If all went well you should have the html report in the 
```/calwebb_spec2/``` directory, along with the file name map of intermediary products. 
Each section of the report specifies what tests passed or failed, and, if there were any 
errors, the report will tell you (with the name of the file where the error occured) if 
it is a pipeline error or a testing tool error. We keep the testing progress and the 
reports in this Confluence page:
https://confluence.stsci.edu/display/JWST/NIRSpec+Pipeline+Testing+Build+7.1+part+1

Please follow these actions for reporting depending on the outcome:
- If there are no testing tool errors and some of the pytests failed, you can check  
off the steps and please put that report in Table 2 of the Confluence page.
- If there are testing tool errors please place your report in the Confluence page, 
add a small comment in the Testing Tool column, and send Maria Pena-Guerrero 
(pena@stsci.edu) an email that you have done this, so the error can be addressed as 
soon as possible.
- If there are errors in the pipeline, please create a Basecamp ticket 
(https://stsci-ins.basecamphq.com/projects/11477312-jwst-pipeline/posts), following the 
guidelines for doing this at:
https://confluence.stsci.edu/display/JWST/NIRSpec+Guidelines+for+Pipeline+Reporting
and copy the link in the corresponding column of Table 2 of the Confluence table.
- Keep updating the html report as you continue with the testing, and including the final
report.


11. Put final results in the NIRSpec vault. At the end of the testing campaign, please do 
the following within your working directory: 
a) Create a directory ```YourNameMODEresults```, e.g. ```MariaIFUresults```
b) Inside that new direcoty place all the intermediary fits products as well as the html
report
c) Also inside that directory, create a text file named ```YourNameMODEresults.txt```, 
e.g. ```MariaIFUresults.txt```. In this text file you will only type the full path where 
you obtained the testing data, e.g. 
```/grp/jwst/wit4/nirspec_vault/prelaunch_data/testing_sets/b7.1_pipeline_testing/
                                                              test_data_suite/IFU_CV3```
d) Finally, place the ```YourNameMODEresults``` directory in the staging directory of the 
NIRSpec vault, and send Gray an email that you have done so at: gkanarek@stsci.edu
The path of the staging directory is:
```/grp/jwst/wit4/nirspec_vault/staging```
	


## TO KEEP IN MIND

- A text file containing an intermediary product name map will be created in the pytests 
directory.
- If any of the central store directory calls do not respond, the pytest will be skipped 
even if the step is set to True in the config file. To make the tests run, you will have 
to download the files the tool is calling, and change the corresponding paths in the 
configuration file.
- The output in the terminal can be a bit overwhelming if there was a failed test or an 
error. In the html report is much clearer to understand what happened.
- In the ```utils``` directory there is a text file named 
```terminal_commands_calwebb_spec2_steps.txt```.
This file contains all the commands you can use from the terminal for running 
calwebb spec2 steps. As part of the testing campaign, it is important that you run the 
pipeline from the command line and that you make sure that the outcome intermediary 
files are the same as those ran with the testing tool. This sanity check is minor but 
important to verify.
- Remember that:
a) Whenever you need to read either the main or science headers of a file,
you can always use the ```read_hdr.py``` script located in the ```utils``` directory of
the testing tool.
b) If you need to change/add a keyword value to a specific extension of a file, you can
use the ```change_keywd.py``` script, also located at the ```utils``` directory.



## POSSIBLE OUTCOMES OF THE PYTESTS

- Passed = the assertion was true, so the test condition was met.
- Failed = the assertion was false, the test condition was NOT met (there will be an 
AssertionError on-screen and in the html file, with a message of what happened).
- Skipped = the test was skipped (there will be a message on the html report and on-screen
explaining why the test was skipped).
- Error = this is a coding error (a bug), please send me an email with the html report. 



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

4. Your script must produce a text file named ```same_name_as_the_script_result.txt``` at
the end, which will only contain one word: True or False.

5. Create a directory named ```same_name_as_the_script`` with the following contents:
a) A text file named ```same_name_as_the_script_init.txt```, and in there type in line 1 
the pipeline step after which your code should be run, in line 2 the language you used, 
and in line 3 the exact command you would use to run the script. If this is an IDL script, 
line 2 should be the exact call from within the IDL session.
b) The script you wrote.
c) One level up from where you created the script directory, create a new directory named 
```YourName_PTT_auxiliary_code```. In this directory, please place the directory 
```same_name_as_the_script```, and a create a text file named 
```YourName_PTT_auxiliary_code.txt```. 
This file will only contain the following path:
```/grp/jwst/wit4/nirspec_vault/pipe_testing_tool/auxiliary_code```

6. Now, copy the ```YourName_PTT_auxiliary_code``` directory in the NIRSpec vault, and 
please send an email to Gray Kanarek (gkanarek@stsci.edu). The vault's path is:
```/grp/jwst/wit4/nirspec_vault/staging```.

7. Your code will be implemented for the next testing campaign and if there is a threshold 
associated with it, it will appear as a variable in the PTT configuration file. For the 
current testing campaign, please place it in the "PTT Comments and/or errors" column of  
Table 2 of our progress keeping confluence page:
https://confluence.stsci.edu/display/JWST/NIRSpec+Pipeline+Testing+Build+7.1+part+1
You will be asked to write a description of your code for the final testing report 
of the current testing campaign.



## Enjoy your pipeline testing!

If you have any question of what a specific step does, you can get a description at:
http://jwst-pipeline.readthedocs.io/en/latest/jwst/pipeline/description.html#stage-2-spectroscopic-pipeline-step-flow-calwebb-spec2
