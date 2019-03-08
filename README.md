# PYTEST TESTING TOOL (PTT)


## What is a Pytest

Simply put, a Pytest is a pass or fail Python test. For instance, with the WCS step, we 
have Python scripts (which we are calling auxiliary code within the frame of the testing 
tool) that compare the pipeline product with the ESA corresponding intermediary file, and 
calculates a difference. The Pytest is to assert if that difference is less than or equal 
to an X threshold value. Hence, a failed test means that the condition was not met. If an
error should occur with the pipeline, the test will be flagged as an error.



## Useful links

- PTT Documentation: 
https://innerspace.stsci.edu/pages/viewpage.action?pageId=123011558

- Testing data sets:
https://innerspace.stsci.edu/display/JWST/NIRSpec+Pipeline+Testing+Build+7.1+part+2

- SCSB GitHub repository: 
https://github.com/spacetelescope/jwst

- Pipeline Documentation:
http://jwst-pipeline.readthedocs.io 



## Quick Start Guide


NOTE: This guide assumes that Conda has been installed. If you have not yet done 
so, please follow the instructions at:
https://astroconda.readthedocs.io/en/latest/
Please use python 3.6 (though python 3.5 is supported as well)


THREE THINGS BEFORE STARTING

I.- You may want to clean your PYTHONPATH so that you do not get mysterious failures. To 
do this simply type the following command in the terminal:
```bash
unset PYTHONPATH
```
You can do this every time you run the pytests, or when you start getting strange 
failures. Another option is not to define PYTHONPATH at all in the .profile (or 
equivalent: .bashrc, .bash_profile, .cshrc, .login, etc.) file.

II.- If you work outside the internal network, i.e. in the visitors network or at home, 
you also want to set the following environment variables in the terminal or add them to 
your .profile (or equivalent) file:
```bash
export CRDS_SERVER_URL=https://jwst-crds.stsci.edu
export CRDS_PATH=${HOME}/crds_cache
```
These changes will not affect your work while working with the internal network at ST.

III.- A brief description of what each pipeline step does, as well as a brief description
of all the pytests implemented in the tool, the tests that are still in concept phase, and 
the tests that are still needed, can be found in the Confluence space for PTT. You can 
get to it from the main page of NIRSpec/NIRSpec JWST Pipeline/NIRSpec Calibration
Pipeline Testing Tool (PTT), or by clicking in the following link:
https://confluence.stsci.edu/pages/viewpage.action?pageId=123011558


QUICK START GUIDE

1. Create the conda environment for testing and get the configuration files.  

a. Conda environment for this testing campaign:
- Current testing version is 7.1.3 Please install this version by typing the following 
command in a terminal:
```bash
conda create -n jwst_b7131 --file url_depending_on_your_system 
```
for the current release candidate, the ulr options are:
- Linux: http://ssb.stsci.edu/releases/jwstdp/0.9.6/latest-linux
- OS X: http://ssb.stsci.edu/releases/jwstdp/0.9.6/latest-osx

As bugfixes are announced your current pipeline software may be updated by issuing the 
command:
```bash
conda update --file http://ssb.stsci.edu/releases/jwstdp/0.9.6/latest-osx
```

Please note that this is the most stable release candidate since 06/08/2018, and until
further notice from SCSB. To be completely sure that this is the latest stable version 
you can go to the GitHub repo of SCSB: https://github.com/spacetelescope/jwst. In the 
README file you will find a table called "Software vs DMS build version map", there you 
will see at the top which one is the latest stable release candidate.

NOTE:
If you need to use the development version of the pipeline then do the following:
```bash
conda create -n jwst_dev -c http://ssb.stsci.edu/conda-dev jwst python=3.6
```
Then, to update the development environment, activate the environment and then type:
```bash
conda update --override-channels -c http://ssb.stsci.edu/conda-dev -c defaults --all
```
To get a specific version of the pipeline you need to do 2 steps:
i) Find the version you want with (newest at the bottom)
```bash
conda search -c http://ssb.stsci.edu/astroconda-dev jwst
```
ii) Install the version you want
```bash
conda install --override-channels -c https://ssb.stsci.edu/astroconda-dev -c defaults jwst=version_you_want
```
where ```jwst=version_you_want``` is usually something like ```jwst=1.2.3dev456```.


b. Configuration files corresponding to this build. Create a directory (e.g. 
```b713cfg_files```) somewhere in your testing working space, and ```cd``` into it. Now 
type the following command within the conda environment you just created (see step 2).
```bash
collect_pipeline_cfgs .
```


2. Activate the conda environment for testing the pipeline, e.g. type:
```bash
source activate your_newly_created_environment
```
If the above command does not work try:
```bash
conda activate your_newly_created_environment
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


3. Install the following plugins for the pytests to work properly. Within the conda 
testing environment, type:
```bash
pip install msgpack
pip install pytest-html
```
NOTE: Every time you create a new conda environment you need to install these plugins.


4. Clone the repository so you have the PTT. To do this click at the top right 
of this page, in the dropdown button that says clone or download, then copy the ulr that
appears there. Now, within the conda testing environment, go to or create the directory 
where you want the PTT to "live" in. However, make sure that the configuration files 
directory is at least at the top level of the directory tree where the PTT will live, e.g. 
the ```b713cfg_files``` directory and the ```nirspec_testing_tool``` directory can be at 
the same level, but the ```b713cfg_files``` directory cannot be inside the 
```nirspec_testing_tool``` directory because otherwise the .cfg files
will be picked up by Git and it will try to put them in your space.
Remember you are in the GitHub repository page so go all the way to the top of the page,
click on the green button and copy the ulr that appears there.
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
cp -r /grp/jwst/wit4/nirspec_vault/prelaunch_data/testing_sets/b7.1_pipeline_testing/test_data_suite/the_data_you_want_to_copy .
```

NOTE:
You can start with the FS benchmark data to make sure you are doing the right thing. To 
get the data go to 
```bash
/grp/jwst/wit4/nirspec_vault/pipe_testing_tool/PTT_FS_benchmark_run
```
There you will find a FS raw file, a ```PTT_config.cfg``` file, and a directory called
```results_491```, which contains all the intermediary fits products obtained from running
```calwebb_detector1```, the output text files from running the corresponding script 
(which include the ```cal_detector1_outputs_and_times_DETECTOR.txt``` and the added 
keywords to the ```_uncal``` file), the all the intermediary fits products obtained from 
running ```calwebb_spec2```, all the plots created with the PTT, the ```report.html```, 
and the PTT text file outputs (```screen_output_DETECTOR.txt``` and 
```True_steps_suffix_map_DETECTOR.txt```). You can use the ```PTT_config.cfg```
(changing the paths appropriately) provided to make sure you obtain the same results 
from the PTT run.


b. In the directory where you copied the test data, you will need to run a script PER
fits file you want to test. Do not worry, this will only be done once. This script will
create a new subdirectory with the necessary input file to run the SSB script 
that converts raw data into uncal type files. You can choose to either keep this 
subdirectory, or tell the script to remove it after the operation is done. In the 
terminal type:
```bash
python /path_to_the_testing_tool/nirspec_pipe_testing_tool/utils/prepare_data2run.py fits_file.fits MODE -u
```
where the MODE is expected to be one of: FS, MOS, IFU, BOTS, or MOS_sim (use this last one
only for MOS simulations, simulations for other modes should use the corresponding mode). 
This command will update the uncal keyword header without creating a new file, and will 
also keep the subdirectory. To remove it, simply add ```-rm``` at the end. To save the 
keyword changes in a new fits file (instead of updating), remove the ```-u```. The new 
uncal fits file is now ready for pipeline ingest.

For now, this script sets the FILTER keyword to default to the corresponding science 
configuration. This is because the pipeline does not know yet what steps to do with
FILTER=OPAQUE.

c. Optional. If you want to see the header of any file, you can use the another script
in the ```utils``` directory of the PTT. If you just want to see on-screen the 
header, go where your fits file "lives" and type:
```bash
python /path_to_the_testing_tool/nirspec_pipe_testing_tool/utils/read_hdr.py fits_file.fits
```
This command will show the main header. To save the header to a text file add a ```-s``` 
at the end. If you want to see/save a different extension add at the end ```-e=1``` for 
extension 1, and so on.

d. Now, the data is ready to be ran through cal_detector1. Please go ahead with step 6
of this guide to do that.


6. Set the PTT configuration file. In a text file editor, you are going to modify the 
configuration file named ```PTT_config.cfg```, which lives at 
```/nirspec_pipetesting_tool/calwebb_spec2_pytests/```. 
This is the file that controls all the input that the tool needs. Please open it and 
make sure that:
- All the paths point to the right places. The files can be located anywhere, but both,
the pipeline and the tool, will run faster if the files are local on your computer.
- The input file for the PTT is the final output file from ```calwebb_detector1```.
- The adequate mode for the data to be tested is set correctly, choices are: FS, IFU,
MOS, BOTS, or MOS_sim.
- The variable ```change_filter_opaque``` should be set to False unless you want to change
the FILTER keyword back to OPAQUE.
- The variable ```raw_data_root_file``` should be the name of the raw file you downloaded
from the NIRSpec vault; for ground observations it starts with NRS. If you are running 
simulations then you can look into the ```ESA_Int_products``` directory and see what is
the name of the directory that corresponds to your data, copy that name and add .fits to
the end, e.g. for my simulation file ```F170LP-G235M_MOS_observation-6-c0e0_001_DN_NRS1.fits```
go into ```/grp/jwst/wit4/nirspec_vault/prelaunch_data/testing_sets/b7.1_pipeline_testing/test_data_suite/simulations/ESA_Int_products```,
then set ```raw_data_root_file = F170LP-G235M_MOS_observation-6-c0e0_001.fits```
- The steps that you want to be ran or not are set to True or False.
- In the bottom part of the file, all the additional arguments for the PTT are 
correct, e.g. threshold values, figure switches, and additional fits files.


7. Run the ```calwebb_detector1``` pipeline. The final output of this is the level 2 data
required to run the PTT. In a terminal, go into the directory where the testing tool lives 
(i.e. at the level of the ```calwebb_spec2_pytests``` directory), and make sure that the 
testing conda environment is on. Please check the name used in the ```DETECTOR``` 
keyword in the fits file (you can use the ```read_hdr.py``` script in the ```utils``` 
directory to dump the header into the screen and/or a text file). Now type the following 
command making sure that the appropriate detector name is used in the name of the output 
text file:
```bash
script -a -t 0 caldetector1_screenout_DETECTOR.txt python ../utils/run_cal_detector1.py /path_where_the_uncal_file_lives/uncal_file.fits
```
The first part of the command will create the text file ```caldetector1_screenout_DETECTOR.txt```
which will contain all the screen output of the pipeline run. This text file will be used 
to determine the time that each step took to run. The next part of the command is what 
actually runs the ```run_cal_detector1.py``` script, which executes the calwebb detector 1
pipeline in a single run, using the configuration file that you have for it in 
the ```utils``` directory. 

WARNING: Do not change the name of the output text file, as it is referenced with that 
name in the PTT code. 

To run calwebb detector 1 step-by-step, simply add ```-sbs``` at the end of  
the previous command. Note that running the pipeline in full for calwebb detector 1 takes 
about half the time as it does running it step-by-step, due to IO processing time.

NOTE ON SCRIPT COMMAND:
- The script command will run on Mac machines, if you are using Linux, please follow these
steps:

i) Make sure you are in ```bash```.

ii) In the terminal type the following command making sure you have the appropriate 
detector name in the output text file:
```bash
python ../utils/run_cal_detector1.py /path_where_the_uncal_file_lives/uncal_file.fits 2>&1 | tee caldetector1_screenout_DETECTOR.txt
```

The command means the following:
- 2> redirect standard error
- &1 to standard output
- | tee pipe all output to the given text file and to console at the same time


If everything went well, you will see a text file called 
```cal_detector1_outputs_and_times_DETECTOR.txt```, which contains the steps ran, the 
name of the output fits file, and the time each step took to run. This text file, along 
with the intermediary products will be located in the path you set for the 
```working_directory``` variable in the configuration file of the PTT.

NOTE ON ```calwebb_detector1``` ERRORS:
- If you were not able to get the file to run though cal detector1 due to an error saying 
that the pipeline was not able to find a best reference for dark or superbias, it is 
possible this is due to the filter keyword  in the main header set to OPAQUE. 
- In this case, you can use an auxiliary script named ```change_filter_opaque2science```, 
which lives at ```/nirspec_pipetesting_tool/calwebb_spec2_pytests/auxiliary_code```. 
- You can run the script from a terminal in any directory by typing:
```bash
python ../calwebb_spec2_pytests/auxiliary_code/change_filter_opaque2science.py file.fits 
```

If all went well and you have a ```gain_scale_DETECTOR.fits``` file, you now are able to 
run the MESA calwebb detector 1 testing tool. Steps to obtain the MESA testing tool:

a. At the same level as the top directory of the PTT (i.e. the ```nirspec_pipe_testing_tool``` 
directory), create a new directory called ```MESA_cal_detector1```.

b. Inside ```MESA_cal_detector1```, you will clone their Git repository. Please follow the 
directions for this at:
http://calibration-pipeline-testing-tool.readthedocs.io/en/latest/

c. Now, in the ```utils/data``` directory of ```nirspec_pipe_testing_tool```, you will 
find a sample json file that you can modify in order to use as input for the MESA calwebb
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
Please record your progress in our testing campaign Confluence page:
https://confluence.stsci.edu/display/JWST/NIRSpec+Pipeline+Testing+Build+7.1+part+2. 

 
8. Ready to run PTT. Go back to the directory where PTT lives and into the 
```calwebb_spec2_pytests``` directory, copy final output file from calwebb detector1 into 
the working directory you indicated in the ```PTT_config.cfg``` file, and make sure 
that the input file for the PTT matches the file you just copied into the working 
directory. Now, to ensure that everything is in order, and to see what pytests will be 
executed and in which order type:
```bash
pytest --collect-only
```

9. Do the first PTT run. As an output of the testing tool you will see an html 
file, ```report.html```, and an intermediary product file name map will appear in the 
```calwebb_spec2_pytests``` directory, and at the end it will be moved to the path you 
indicated at the ```PTT_config.cfg``` file with the variable ```working_directory```. The 
output fits files of intermediary products will also be saved in the working directory. 
Please check the name used in the ```DETECTOR``` keyword in the fits file (you can use 
the ```read_hdr.py``` script in the ```utils``` directory to dump the header into the 
screen and/or a text file). In the terminal type:
```bash
script -a -t 0 calspec2_screenout_DETECTOR.txt pytest -s --config_file=PTT_config.cfg --html=report.html --self-contained-html
```
The first part of the command will create the text file ```calspec2_screenout_DETECTOR.txt``` 
in the ```calweb_spec2_pytests``` directory. The text file will contain all the screen 
output of the pipeline run, and it will be moved to the working directory at the end. The 
contents of the text file will be used to determine the time that each step took to run. 

The next part of the above command runs the pytest scripts; the ```-s``` will capture all 
the print statements in the code on screen. For now, the report will always be saved in 
the ```calweb_spec2_pytests``` directory, and you will have to manually move it to where
you want it to live, unless you give the full path of the destination in the ```--html=``` 
flag.

NOTE: 
- If you are running the tool in Linux, the script command will not work. Please see the 
NOTE ON SCRIPT COMMAND of step 7 of this guide.

- If you are doing a re-run of the pytest only, be careful not to overwrite the 
```calspec2_screenout_DETECTOR.txt``` file, in this case you should use the following 
command to run the pytest:
```bash
pytest -s --config_file=PTT_config.cfg --html=report.html --self-contained-html > pytest_only_DETECTOR.txt
```

-> To only run a few pipeline steps you need to:
a) Make sure that the variable ```run_calwebb_spec2``` in the ```PTT_config.cfg``` file
is set to False (otherwise the pipeline will run in full and we have no control of
individual steps).
b) Turn off (i.e. set to False) the steps you do not want to run in the ```PTT_config.cfg``` 
file, which are located in the section ```run_pipe_steps``` of the file.

-> To run a few pytest you need to select which pytest to run (i.e. set to True) in the 
```PTT_config.cfg``` file, which are located in the section ```run_pytest``` of the file.

-> To only run pytest and skip running the pipeline entirely:
a) Make sure that the variable ```run_calwebb_spec2``` in the ```PTT_config.cfg``` file
is set to False.
b) Set to False all the pipeline steps in the ```PTT_config.cfg``` file. The steps are 
located in the section ```run_pipe_steps``` of the file.
c) Set to True all pytest you want to run in the ```PTT_config.cfg``` file. These are 
located in the section ```run_pytest``` of the file.

NOTE: Setting a pytest to False in the ```PTT_config.cfg``` file will skip running all 
the pytest for that step (i.e. completion, reference file, and validation pytest).


10. Report your findings. If all went well, you should have the html report in the 
```calwebb_spec2_pytest``` directory. Please move it manually to the working directory. 
The file name map of intermediary products, the screen_output text file, and all the 
intermediary products should already be in the working directory. 

Each section of the html report specifies what tests passed or failed, and, if there were 
any errors, the report will tell you (with the name of the script where the error 
occurred) if it is a pipeline error or a PTT error. Keep updating your testing progress 
in the Confluence page where you obtained the name of your testing data sets, in the 
column named Test reports (if there is none, please go ahead and add it). In this column 
place a copy of your html report and the text file of the screen output. 

Please follow these actions for reporting your progress, depending on the outcome:
* The html file is portable, i.e. it will not get corrupted if you move it to another
directory. 
- Place the html file in the Confluence page. You can also create a PDF export from the 
report and use that file for the Confluence page and for any sharing purposes.
- If there are no testing tool errors and some of the pytests failed, you can check  
off the steps and please put that report in Table 2 of the Confluence page.
- If there are testing tool errors please place your report in the Confluence page, 
add a small comment in the "Testing Tool" column of Table 2 in the Confluence page, and 
send an email to Maria Pena-Guerrero (pena@stsci.edu), so that the 
PTT error can be addressed as soon as possible.
- If there are errors with the pipeline, please do the following:
i) Look up if there are any JIRA issues already filed related to the error.
ii) If there are JIRA issues open you can add a comment that you also encountered the
problem for your data.
iii) If there are closed JIRA issues you can re-open them and report in the comments 
section.
iv) If there is no JIRA issue for the problem please go ahead and create one. You can find
the NIRSpec guidelines for doing so at:
https://confluence.stsci.edu/display/JWST/NIRSpec+Guidelines+for+Pipeline+Issue+Reporting
and copy the link (or JIRA issue number) next to your report.
- Keep updating the html report as you continue with the testing, and including the final
report.


11. Put final results in the NIRSpec vault. There are two stages of putting results in the
vault: (i) After running ```calwebb_detector1``` and (ii) after running ```calwebb_spec2```
and at the end of the testing campaign. For either stage please do the following within 
your working directory:

a. Create a new directory called ```yourname_MODE_Name_of_the_raw_data_file_caldetector1```
for stage (i), and ```yourname_MODE_Name_of_the_raw_data_file_calwebbspec2``` for stage 
(ii), e.g. ```maria_IFU_NRSSMOS-MOD-G1M-17-5344175105_1_491_SE_2015-12-10T18h00m06_calwebbspec2```.
(If you already copied a previous version of your results and you want to keep it then add
a ```_v2``` -or the version number that corresponds- after the word results, e.g. 
```maria_IFU_NRSSMOS-MOD-G1M-17-5344175105_1_491_SE_2015-12-10T18h00m06_calwebbspec2_v2```) 

b. Inside that new direcoty place all the intermediary fits products as well as the html
report, and other report text files.

c. Also inside that directory, create a text file named ```path2results.txt```. 
In this text file you will only type the full path where you obtained the testing data, 
e.g. 
```
/grp/jwst/wit4/nirspec_vault/prelaunch_data/testing_sets/b7.1_pipeline_testing/test_data_suite/IFU_CV3
```

d. Place the results directory you created in the staging directory of the 
NIRSpec vault, and send Gray an email (gkanarek@stsci.edu). The path of the staging 
directory is:
```/grp/jwst/wit4/nirspec_vault/staging```

e. Finally, change the permissions recursively on all files in your newly created 
directory, so that the NIRSpec curators team can move the data. To do this type the 
following command in a terminal:
```chmod -R 777 /grp/jwst/wit4/nirspec_vault/staging/your_newly_created_directory```
	


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
In the ```utils/data``` directory there are two text files named 
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
