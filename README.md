This pipeline should be able to automatically calibrate phase reference data from the VLBA, but I make no garauntees. I know it runs in CASA 6.4, but I think it should also run in CASA 5.7. Again, no promises.

You can run each script (make sure you go in number order) individually by typing

> casa --nologger -c /path/to/the/phasecal_pipeline_0X_the rest of the name.py

Or you can run the whole pipeline at once using:

> casa --nologger -c /path/to/the/phasecal_pipeline.py

**impotant note**

This pipeline depends on the UF_PIPELINE_UTILS file. You either need to have this file in the directory you are running the script, or you need to add the file to your .casa directory. Instructions for that are below. 

- go to ~/.casa
- make a place to store the UF_PIPELINE_UTILS.py file
    - for example I have /Users/lucas.hunt/.casa/USNO_VLBA_PIPELINE on my mac and /home/lucas.hunt/USNO_VLBA_PIPELINE on fridadev
- if it doesn't already exist, make an init.py file (for casa5.7) or config.py file (for casa6.4)
- in the file add the lines
    - import sys
    - sys.path.append('/path/to/directory/with/UF_PIPELINE_UTILS.py)
    - for example on my Mac I have sys.path.append('/Users/lucas.hunt/.casa/USNO_VLBA_PIPELINE')

If you follow the above instructions correctly you should be able to start casa and import UF_PIPELINE_UTILS

