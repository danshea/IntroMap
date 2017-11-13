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

Here is a sample invocation of `IntroMap.py` using our \#174-12-26 genome data and the accompanying output.

<pre>
(IntroMap) $ ./IntroMap.py -i ../bowtie2_alignment/against_bra/174-12-26-1-28280253.bra.bowtie2.sorted.nodups.bam -r ../bowtie2_alignment/against_bra/bra.fa -o 174-12-26-1-28280253.out
63536100 Alignment records were processed from ../bowtie2_alignment/against_bra/174-12-26-1-28280253.bra.bowtie2.sorted.nodups.bam using Reference ../bowtie2_alignment/against_bra/bra.fa
20825947 records were processed for A08
38883802 records were processed for A09
16404182 records were processed for A10
26938828 records were processed for A02
31764690 records were processed for A03
26790030 records were processed for A01
25209370 records were processed for A06
25875098 records were processed for A07
19268591 records were processed for A04
25302534 records were processed for A05
A02     228453  6213700
A09     1       3703909
</pre>

And here are the files that are created as a result of running `IntroMap.py`.

<pre>
-rw-rw-r-- 1 dshea dshea 398873704 11月 13 18:42 A01.174-12-26-1-28280253.out
-rw-rw-r-- 1 dshea dshea         0 11月 13 18:44 A01.174-12-26-1-28280253.out.predicted.out
-rw-rw-r-- 1 dshea dshea    108824 11月 13 18:44 A01_174-12-26-1-28280253_out_predicted.png
-rw-rw-r-- 1 dshea dshea 401089206 11月 13 18:41 A02.174-12-26-1-28280253.out
-rw-rw-r-- 1 dshea dshea        19 11月 13 19:04 A02.174-12-26-1-28280253.out.predicted.out
-rw-rw-r-- 1 dshea dshea     96351 11月 13 19:04 A02_174-12-26-1-28280253_out_predicted.png
-rw-rw-r-- 1 dshea dshea 472940725 11月 13 18:42 A03.174-12-26-1-28280253.out
-rw-rw-r-- 1 dshea dshea         0 11月 13 19:08 A03.174-12-26-1-28280253.out.predicted.out
-rw-rw-r-- 1 dshea dshea    125827 11月 13 19:08 A03_174-12-26-1-28280253_out_predicted.png
-rw-rw-r-- 1 dshea dshea 286888025 11月 13 18:44 A04.174-12-26-1-28280253.out
-rw-rw-r-- 1 dshea dshea         0 11月 13 19:17 A04.174-12-26-1-28280253.out.predicted.out
-rw-rw-r-- 1 dshea dshea    116349 11月 13 19:17 A04_174-12-26-1-28280253_out_predicted.png
-rw-rw-r-- 1 dshea dshea 376726291 11月 13 18:44 A05.174-12-26-1-28280253.out
-rw-rw-r-- 1 dshea dshea         0 11月 13 19:42 A05.174-12-26-1-28280253.out.predicted.out
-rw-rw-r-- 1 dshea dshea    118286 11月 13 19:42 A05_174-12-26-1-28280253_out_predicted.png
-rw-rw-r-- 1 dshea dshea 375339187 11月 13 18:43 A06.174-12-26-1-28280253.out
-rw-rw-r-- 1 dshea dshea         0 11月 13 19:43 A06.174-12-26-1-28280253.out.predicted.out
-rw-rw-r-- 1 dshea dshea    124054 11月 13 19:43 A06_174-12-26-1-28280253_out_predicted.png
-rw-rw-r-- 1 dshea dshea 385251390 11月 13 18:43 A07.174-12-26-1-28280253.out
-rw-rw-r-- 1 dshea dshea         0 11月 13 19:51 A07.174-12-26-1-28280253.out.predicted.out
-rw-rw-r-- 1 dshea dshea    117447 11月 13 19:51 A07_174-12-26-1-28280253_out_predicted.png
-rw-rw-r-- 1 dshea dshea 310074562 11月 13 18:40 A08.174-12-26-1-28280253.out
-rw-rw-r-- 1 dshea dshea         0 11月 13 19:51 A08.174-12-26-1-28280253.out.predicted.out
-rw-rw-r-- 1 dshea dshea    116213 11月 13 19:52 A08_174-12-26-1-28280253_out_predicted.png
-rw-rw-r-- 1 dshea dshea 578936471 11月 13 18:41 A09.174-12-26-1-28280253.out
-rw-rw-r-- 1 dshea dshea        14 11月 13 19:55 A09.174-12-26-1-28280253.out.predicted.out
-rw-rw-r-- 1 dshea dshea     93137 11月 13 19:55 A09_174-12-26-1-28280253_out_predicted.png
-rw-rw-r-- 1 dshea dshea 244240188 11月 13 18:41 A10.174-12-26-1-28280253.out
-rw-rw-r-- 1 dshea dshea         0 11月 13 19:55 A10.174-12-26-1-28280253.out.predicted.out
-rw-rw-r-- 1 dshea dshea    120184 11月 13 19:55 A10_174-12-26-1-28280253_out_predicted.png
</pre>

The `*.predicted.out` files hold identified introgressed loci coordinates for a chromosome.  
<pre>
(IntroMap) $ cat A02.174-12-26-1-28280253.out.predicted.out
A02     228453  6213700
</pre>
And the _chrN_._outpostfix_.out files hold the homology scores at each position for the locally weighted linear regressions. Currently, the cli version re-generates these each time, but because this is a costly computational operation, we plan to add the option to re-plot from previously generated runs out files. Currently, this can be done using the Jupyter notebook. Instructions on how to do so are provided in the notebook cell that generates the plots.
