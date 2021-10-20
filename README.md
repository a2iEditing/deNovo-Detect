# deNovo Detect
A tool for the *de novo* detection of editing sites from RNA-seq data


## Installation and Requirements
### Dependencies
- _[SAMtools](http://samtools.sourceforge.net/)_
- _[bedtools](https://bedtools.readthedocs.io/en/latest)_ - version 2.27 or higher
- _[GAWK](https://www.gnu.org/software/gawk/)_ - version 5.0.0 or higher (tested on 5.0.0)

- _[Perl](https://www.perl.org/get.html)_ - version 5.1.6 or higher (tested on 5.1.6)
- _[Python 3](https://www.python.org/downloads/source/)_ - version 3.7 and up

#### Additional Python Libraries Required 
- logomaker
- pandas
- matplotlib
- pybedtools


### OS Requirements
**Currently, the program supports only GNU/Linux operating systems** (and probably any other POSIX OS)

##### CPU and Memory
The program has low demand for system resources (CPU and memory) - only the default resource requirements of SAMtools and bedtools are needed (run with default CPU and memory parameters to generate the CMPileups). For the rest of the processing, the program demands very little. However **the default thread number is high** (and can be easily changed using command line parameters)

##### Disk Space
The installation requires a bit more than 7G of free disk space, almost all of which is for the built-in resources
## Running
Run _AnalyzeTissue.py -h_  to see full help.

### An example for a simple run:
```
_pythonInterpeter scriptLocation/AnalyzeTissue.py -a _alignmentDirectory_ -i _regionFile_ -o _outputDirectory_
```

### Inputs

#### ALignment 
An input directory containing alignment (BAM) files (it may be nested).   
**Note: alignment should be unique.** (non-unique alignment may create unpredicted, algorithm dependent, and biases)

#### Region file
A BED6 file, with the columns: chromosome, start position, end position, name, empty field (padded with 0), strand bias (+\-).
**Currently, only the CDS region file is avialble for the script by default*

#### Genomic files
To run this script, some genomic files needed to be downloaded and filtered. for guidelines, please download and follow the _[Guidelines for creating genome files](https://github.com/zivtzur6/Orshai_sites_analist_pipline/blob/main/Guidelines%20for%20creating%20genome%20files.docx)_

### Output

#### Final list
A list of the detected editing site for cluster sizes of: 0, 50, 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200, 102400.
In addition, a list of the aforementioned editing sites filtered by minimum editing rate of 0.02% and minimum coverage per site of 100 reads by default (the user may change these).

#### Logo graphs
For each pf the aforementioned cluster sizes, a logo graph (i.e. motif) for all mismatchs types. 

#### Log file
A log file created at the predefined location. If the user didn't specify a path, the log file would be generated inside the "Log" directory.


© 2021 Bar-Ilan University (Erez Y. Levanon, Erez.Levanon@biu.ac.il; Eli Eisenberg, elieis@post.tau.ac.il; Ziv Tzur, tzur.zivh@live.biu.ac.il).























