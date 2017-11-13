Dan Shea  
2017.06.07  

# IntroMap
Introgression screening software IntroMap

Please contact me (Dan Shea) if you have any questions regarding the use of IntroMap.  

I have tried to document the use of the software in the Jupyter notebook, as well as in a paper that we have recently
submitted to BMC Genetics for peer review. Our goal is to create software that will assist with (or possibly replace) marker
based assays for the screening of interspecific hybrids. If the paper is accepted I'll provide the information for citation
should you decide to employ this software in your own research.

If you find the software useful, or you have any ideas on how we can improve the methodology, please feel
free to drop us a line. One of the functions we would like to add is automation of the parameter tuning step.

Don't use IntroMap 30 minutes after eating, your mileage may vary, etc., etc., ad infinitum.  

Copyright (c) 2016, Daniel J. Shea, Niigata University, Dept. of Agriculture, Laboratory of Plant Breeding All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Anaconda, conda and the runtime environment dependencies
One of the things provided in this repository is a `YAML` file for conda, so that the environment in which IntroMap was developed can easily be duplicated. I'd highly recommend installing anaconda by following the instructions provided here: https://docs.continuum.io/anaconda/ and getting yourself a nice clean environment, as opposed to messing about with an OS supplied copy of python.

Once you have `conda` installedon your system, you can utilize the `YAML` file to replicate the environment.
You will need to make minor modifications to the file prior to setting up the new environment on your system.

The conda documentation can be found here https://conda.io/docs/index.html. But in particular, you will want to look at the following https://conda.io/docs/using/envs.html#id11. I've duplicated the instructions here for brevity, but I always recommend to RTFM.

(Note: You may need to edit the prefix in the environment file first to point to where your conda installation stores environments on your system.)

To create an environment using a supplied file, run:  
`conda env create -f IntroMap_conda_environment.yml`

To activate the new environment:  
`source activate IntroMap`

## Once you have an environment installed and activated, you can either run the Notebook or try the cli version
The cli version is `IntroMap.py`.  
First give it proper executable permissions i.e. - `chmod 755 IntroMap.py`  
Then, invoke it as follows to see the command line options i.e. - `./IntroMap.py -h`  
You should see the following output:  

<pre>
(IntroMap) $ ./IntroMap.py -h                                                                                                                                                                                                                    
usage: IntroMap.py [-h] -i INFILE -r REFERENCE -o OUTFILE [-t THRESHOLD]
                   [-b {True,False}] [-f FRAC]

IntroMap

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        The name of the BAM file to be examined.
  -r REFERENCE, --reference REFERENCE
                        The reference genome used for the aligned BAM file.
  -o OUTFILE, --outfile OUTFILE
                        The output filename format. Example: -o lowess would
                        create chr1.lowess to chrN.lowess outputs.
  -t THRESHOLD, --threshold THRESHOLD
                        The threshold value used by the binary classifier.
                        Default is 0.90
  -b {True,False}, --below {True,False}
                        Boolean flag for the threshold. True = report regions
                        below the threshold. False = report regions at or
                        above threshold. Default is True
  -f FRAC, --frac FRAC  The windowsize used when performing the locally
                        weighted linear regression, given as a fraction of the
                        chromosome length. Default is 0.05

If you use IntroMap, please consider citing our paper.
</pre>

