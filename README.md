# PYTEST TESTING TOOL

## What is a Pytest

Simply put, a Pytest is a pass or fail Python test. For instance, with the WCS step, we have Python scripts (which we are calling auxiliary code within the frame of the testing tool) that compare the pipeline product with the ESA corresponding intermediary file, and calculates a difference. The Pytest is to asses if that difference is less than or equal to an X threshold value. 



## Quick Start Guide


** Note ** : This guide assumes that Conda has been installed. If you have not yet done so, please follow the instructions at:
https://astroconda.readthedocs.io/en/latest/


1. Create the conda environment for testing the correct version of the pipeline
* Current testing version is 7.1
* Specific instructions to create the build7.1 environment.
For example, in a terminal type:
```bash
conda create -n jwst_b7.1 --file url_to_be_given_to_us
```

2. Activate the conda environment for testing the pipeline, e.g. type:
```bash
source activate jwst_b7.1
```

3. Install the pytest html plugin. Within the conda environment type:
```bash
pip install pytest-html
```

4. In a text file editor, you are going to modify the configuration file "cwspec2_config.cfg" (but everything should be ok for the first run but better to check).
- This is where you will give all the arguments to the testing tool.
- Make sure all the paths are what you need.
- Set to True or False the steps that you want to be ran or not.
- Modify the threshold values and figure switches in the bottom part of the file.

5. See what pytests will be executed and in which order. Within the active conda environment type:
```bash
pytest --collect-only
```

6. Do the first run and create a "report.html" file as an output of the testing tool. In the terminal type:
```bash
pytest -s --config_file=cwspec2_config.cfg --html=report.html
```


If all went well you should have a report.html in your directory.


NOTE THAT:
- A text file containing a suffix name map will be created in the pytests directory.
- If any of the central store directory calls do not respond, the pytest will be sikipped even if the step is set to True in the config file. To make the tests run, you will have to download the files the tool is calling, and change the corresponding paths in the configuration file.
- The output in the terminal can be a bit overwhelming if there was a failed test or an error. The html report is much clearer to understand what happened.


If you have any question of what a specific step does, you can get a description at:
http://ssb.stsci.edu/doc/jwst_dev/jwst/pipeline/description.html#stage2-imaging-flow


Enjoy your pipeline testing!
