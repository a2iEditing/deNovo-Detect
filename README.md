# Orshai sites analyst pipeline
A tool for analyze editing sites from RNA-seq data
## Installation and Requirements
### Dependencies
- _[SAMtools](http://samtools.sourceforge.net/)_ - version 1.9 or higher (tested on 1.9)
- _[bedtools](https://bedtools.readthedocs.io/en/latest)_ - version 2.30 or higher (tested on 2.30)
- _[GAWK](https://www.gnu.org/software/gawk/)_ - version 5.0.0 or higher (tested on 5.0.0)

- _[Perl](https://www.perl.org/get.html)_ - version 5.1.6 or higher (tested on 5.1.6)
- _[Python 3.7](https://www.python.org/downloads/source/)_ (a clean installation is sufficient)
#### Required Python Libraries
- os
- shutil
- logging
- argparse
- logomaker
- pandas
- glob
- multiprocessing
- subprocess
- matplotlib.pyplot
- multiprocessing.Pool
- pybedtools.bedtool


### OS Requirements
**Currently, the program supports only GNU/Linux operating systems** (and probably any other POSIX OS)

##### CPU and Memory
The program has low demand for system resources (CPU and memory) - only the default resource requirements of SAMtools and bedtools are needed (run with default CPU and memory parameters to generate the CMPileups). For the rest of the processing, the program demands very little. However **the default thread number is high** (and can be easily changed using command line parameters)

##### Disk Space
The installation requires a bit more than 7G of free disk space, almost all of which is for the built-in resources (built-in genomes and tables which are not mandatory for running, see further details below for installation without downloading and running)

## Running
Run _AnalyzeTissue.py -h_  to see full help.

### An example for a simple run:
```
_pythonInterpeter scriptLocation/AnalyzeTissue.py -a _alignmentDirectory_ -i _regionFile_ -o _outputDirectory_
```

### Inputs

#### ALignment 
The input directory contains alignment (BAM) files. Therefore, the program looks for the BAM files within all alignment directories.  
**Note: alignment should be unique.** (non-unique alignment may create unpredicted, algorithm dependent, and biases)

#### Genome
Currently, only hg38 fit for this script

#### Region file
A 6 columns bed file needs to include: chromosome, start position, end position, name, zero, strand bias.
*Currently, only the CDS region file is fit to the script

### Output

#### Final list
A list of the editing site for each cluster size of: 0, 50, 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200, 102400.
In addition, a list of all the editing sites filtered by minimum editing rate of 0.02 and minimum read per site of 100. the user can change those values.

#### Logo graphs
Logo graph of all the mismaches type orgenaize in cluster size of: 0, 50, 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200, 102400

#### Log file
The log file will create in the path that the user inserts. If the user didn't insert the log file path, the log file would generate inside the Log directory.























