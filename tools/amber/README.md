# Installation of Amber for structural relaxation

We recommend downloading and installing the appropriate version of **Amber** from (http://ambermd.org). The results reported in the manuscript were obtained using **Amber version 14 (Amber14)**.

After installation, you need to update the Amber path in the Amber-based refinement scripts.

Specifically, modifying the following line in the script `/tools/amber/3amberRefine1.sh`:
```
source /share/amber14/amber.sh
```
Replace it with the path to your own Amber installation, like this:
```
source /path/to/your/amber/amber.sh
```

Similarly, update the following lines in the script `/tools/amber/fixamberinput.sh`:
```
protreslib="/share/amber14/dat/reslib/leap/protein.amberua.prepin"
nareslib="/share/amber14/dat/reslib/leap/na.amberua.prepin"
```
Replace them with:
```
protreslib="/path/to/your/amber/dat/reslib/leap/protein.amberua.prepin"
nareslib="/path/to/your/amber/dat/reslib/leap/na.amberua.prepin"
```
