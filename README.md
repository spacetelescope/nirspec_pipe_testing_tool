# PYTEST TESTING TOOL

## What is a Pytest

Simply put, a Pytest is a pass or fail Python test. For instance, with the WCS step, we 
have Python scripts (which we are calling auxiliary code within the frame of the testing 
tool) that compare the pipeline product with the ESA corresponding intermediary file, and 
calculates a difference. The Pytest is to asses if that difference is less than or equal 
to an X threshold value. Hence, a failed test means that the condition was not met. If an
error should occur with the pipeline, the test will be flagged as an error.



## Quick Start Guide


*** Note *** : This guide assumes that Conda has been installed. If you have not yet done 
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
the ulr options are:
- Linux: http://ssb.stsci.edu/releases/jwstdp/0.7.8/latest-linux
- O5S X: http://ssb.stsci.edu/releases/jwstdp/0.7.8/latest-osx


2. Activate the conda environment for testing the pipeline, e.g. type:
```bash
source activate jwst_b7.1
```
* NOTE: 
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

4. In a text file editor, you are going to modify the configuration file 
"cwspec2_config.cfg" (but everything should be ok for the first run but better to check).
- This is where you will give all the arguments to the testing tool.
- Make sure all the paths are what you need.
- Set to True or False the steps that you want to be ran or not.
- Modify the threshold values and figure switches in the bottom part of the file.

5. See what pytests will be executed and in which order. Within the active conda 
environment type:
```bash
pytest --collect-only
```

6. Do the first run and create a "report.html" file as an output of the testing tool. In 
the terminal type:
```bash
pytest -s --config_file=cwspec2_config.cfg --html=report.html
```
The "-s" will capture all the print statements in the code on screen.


If all went well you should have a report.html in your directory.


NOTE THAT:
- A text file containing a suffix name map will be created in the pytests directory.
- If any of the central store directory calls do not respond, the pytest will be skipped 
even if the step is set to True in the config file. To make the tests run, you will have 
to download the files the tool is calling, and change the corresponding paths in the 
configuration file.
- The output in the terminal can be a bit overwhelming if there was a failed test or an 
error. The html report is much clearer to understand what happened.



## POSSIBLE OUTCOMES OF THE PYTESTS

- Passed = the assertion was true, so the test condition was met.
- Failed = the assertion was false, the test condition was NOT met (there will be an 
AssertionError on-screen and in the html file, with a message of what happened).
- Skipped = the test was skipped (there will be a message on the html report and on-screen
explaining why the test was skipped).
- Error = this is a coding error (a bug), please send me an email with the html report. 



If you have any question of what a specific step does, you can get a description at:
http://ssb.stsci.edu/doc/jwst_dev/jwst/pipeline/description.html#stage2-imaging-flow


Enjoy your pipeline testing!
