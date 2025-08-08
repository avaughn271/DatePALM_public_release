# DatePALM

## Purpose
DatePALM is a program to infer patterns of polygenic selection, including time-varying selective pressures and selection on correlated traits, using both modern and ancient DNA. 

## Running DatePALM

DatePALM has 3 steps:

## (1) Format GWAS summary statistics (or use the pre-specified summary statistics available in the XXX Directory)


There are 2 acceptable formats for GWAS metadata, one for analyzing 1 trait at a time and one for jointly analyzing 2 traits.

### 1 Trait
An example metadata file would look like:
```bash
LD_block	variant	derived_allele	minor_allele	minor_AF	beta	se	pval
8	1:11438910:T:C	C	C	0.0589645	0.00759057	0.00222401	0.00019
10	1:14007558:T:C	C	T	0.292772	-0.00415705	0.00112472	0.00038
11	1:16407987:G:T	G	G	0.281201	-0.00465798	0.0011484	4.6e-06
13	1:18907068:G:C	C	C	0.105726	-0.00809486	0.00176444	4.3e-05
17	1:25251326:G:A	A	G	0.425563	-0.00669086	0.00103425	6.9e-13
18	1:26987646:A:G	G	G	0.023938	-0.0126736	0.00343447	1.8e-05
```
It contains the following columns

* LD_block - the LD block of the human genome in which this variant is located (see for example [Berisa and Pickrell (2016)](https://doi.org/10.1093/bioinformatics/btv546)).
  
* variant - the ID of this variant, formatted as (CHR:BP:REF:ALT). It is assumed that the effect size *beta* measures the effect of the ALT allele (represented by the last character in the variant ID)
  
* derived_allele - the base corresponding to the derived alelle. It is assumed that the fitted SNP likelihoods (as described in Step 2) are for the derived allele

* minor_allele - the base corresponding to the minor alelle.

* minor_AF - the frequency of the minor allele

* beta - the effect size of the ALT allele on the trait (as defined above)

* se - the standard error of the effect size *beta*

* pval - the p-value of the GWAS for this SNP

### 2 Traits
An example metadata file would look like:
```bash
LD_block	variant	derived_allele	minor_allele	minor_AF	beta@TRAIT1	se@TRAIT1	pval@TRAIT1	beta@TRAIT2	se@TRAIT2	pval@TRAIT2
2	1:2452357:G:A	A	A	0.22085	-0.00712539	0.00174441	5.2e-05	-0.003	0.003	0.3685
2	1:1940071:C:T	T	T	0.154144	0.0012429	0.00200304	0.51	-0.016	0.004	0.0004858
3	1:4315204:G:T	T	G	0.46067	-0.00492573	0.00145255	0.00056	-0.003	0.002	0.3012
3	1:3830041:A:G	A	G	0.330162	0.000751367	0.0015385	0.65	-0.009	0.003	0.001125
4	1:5701154:T:C	C	C	0.157624	-0.00926937	0.0019892	3.2e-06	-0.006	0.003	0.09926
4	1:4888739:C:G	G	G	0.0630181	0.00188014	0.00297756	0.5	-0.023	0.007	0.001168
5	1:6995522:C:T	T	T	0.0409137	-0.0140566	0.00365428	0.00011	0.008	0.006	0.1772
```

* The columns are defined identically as in the 1 trait case, except that there are now 2 beta, se, and pval columns, one for each of the two traits. We denote which trait each value corresponds to by appending @TRAIT to the end of each column name where TRAIT is the name of the trait we are considering.


## (2) Fit the likelihood function of selection for a set of SNPs.

In this step, we calculate the likelihood function of a certain SNP, L<sup>SNP</sup>(s).

We do this with the command

```bash
$ python PATH/snp_lik.py
```

The arguments supplied to this script are almost identical to the ones for the *inference.py* script of CLUES2. A look at the [CLUES2 GitHub](https://github.com/avaughn271/CLUES2/tree/main) might be helpful for understanding how these arguments work.

## This step takes as input:

**--times** The input file of coalescent times. Identical to that described in CLUES2. You can supply any combination of a **--times** file, an **--ancientSamps** file, and/or an **--ancientHaps** file. See the [CLUES2 GitHub](https://github.com/avaughn271/CLUES2/tree/main) for more information on this format.

**--ancientSamps** The input file of ancient genotype probabilities. Identical to that described in CLUES2. You can supply any combination of a **--times** file, an **--ancientSamps** file, and/or an **--ancientHaps** file. See the [CLUES2 GitHub](https://github.com/avaughn271/CLUES2/tree/main) for more information on this format.

**--ancientHaps** The input file of ancient haplotype probabilities. Identical to that described in CLUES2. You can supply any combination of a **--times** file, an **--ancientSamps** file, and/or an **--ancientHaps** file. See the [CLUES2 GitHub](https://github.com/avaughn271/CLUES2/tree/main) for more information on this format.

**--popFreq** The modern derived allele frequency.

**--N** The effective population size (Ne). This is the HAPLOID effective population size, denoting the number of haplotypes in a given population, NOT the number of diploid individuals. 100 diploid humans corresponds to an N value of 200. Either this or a coal file must be supplied, but not both.

**--coal** The population size file denoting different population sizes through time. Identical to the file format used by [RELATE](https://myersgroup.github.io/relate/). The population sizes considered here are the HAPLOID population size, denoting the number of haplotypes in a given population, NOT the number of diploid individuals. 100 diploid humans corresponds to an N value of 200. Either this or a value of N must be supplied, but not both.

**--tCutoff** The maximum time (in generations) to be considered in the analysis.

**--df** (optional) This is the number of discretization points that is used to bin the allele frequencies. A higher number will result in less error due to rounding of allele frequencies but at increased computational cost. We find that having a finer discretization (higher df) is more important when N is large. This is because large population sizes result in smaller allele frequency fluctuations from generation to generation and only a fine discretization grid will be able to model them accurately. The most rigorous way to set df is to steadily increase df until the results appear to have converged, but we find that practically the default value of 450 is sufficient for nearly all cases.

**--timeBins** (optional) A list of epoch breakpoints, sorted in increasing order. These, along with 0 and tCutoff, give the endpoints of the disjoint time intervals in which independent selection coefficients will be inferred. For example, if --timeBins 200 300 is used, 3 separate selection coefficients will be inferred, for each of the intervals [0,200), [200,300), and [300,tCutoff). If this argument is not supplied, one selection coefficient will be inferred for the interval [0,tCutoff). The youngest breakpoint must be greater than 0 and the oldest breakpoint must be less than tCutoff.

**--out** The prefix of the output file. This should be "*DIR*/*chr*_*bp*" where *chr* is the chromosome the SNP is on, *bp* is the basepair position of that SNP on that chromosome, and *DIR* is the directory where this file should be stored (possibly the current working directory if you wish; directory must already exist). For example "SNPs/7_123456", would be the proper prefix for a SNP at position 123456 on chromosome 7, which would be stored in the directory SNPs. "7_123456" would also be acceptable, in which case the file would be stored in the current working directory. 

## This step will produce (in the current working directory)

***out.txt***  A .txt file representing the fitted form of the likelihood function. A (possibly multivariate) Gaussian distribution is fit to the likelihood function. If *n* epochs are considered, then this file will have *2n* lines. The first *n* lines list the mean of the Gaussian distribution. The following *n* lines each denote a row of the covariance matrix of the Gaussian distribution.

An example output where 1 epoch is considered (which will be the case when the **--timeBins** argument is not used):
```bash
0.0065
0.00001
```
This denotes that the likelihood function that is fit to this SNP is a Gaussian distribution with mean 0.0065 and variance 0.00001.

An example output where 3 epochs are considered (for example if  **--timeBins 100 200** is used):
```bash
0.0288910757054824
0.0015543206761774
-0.0003582844929024
0.0003224682663718 -0.0001247559211745 -0.0000819910274343
-0.0001247559211745 0.0016714404724434 -0.0005475150222164
-0.0000819910274343 -0.0005475150222164 0.0002322068500070
```
This denotes that the likelihood function that is fit to this SNP is a multivariate Gaussian distribution with mean 

[0.0288910757054824, 0.0015543206761774, -0.0003582844929024] 

and a covariance matrix of 

[[0.0003224682663718, -0.0001247559211745, -0.0000819910274343], 

[-0.0001247559211745, 0.0016714404724434, -0.0005475150222164],

[-0.0000819910274343, -0.0005475150222164, 0.0002322068500070]].

The selection coefficients are indexed beginning at the present and moving backwards in time. This means that when considering multivariate Gaussian distributions, the 1st entry corresponds to the most recent epoch, the 2nd entry corresponds to the next most recent epoch, etc.


## (3) Estimate selection gradients

After running the previous step for all the requisite SNPs, we can then estimate the selection gradient with the following function

```bash
$ python PATH/palm.py
```

## This step takes as input:

**--metadata** The name of a .tsv containing GWAS summary statistic information for each SNP. Can be formatted for either analyzing one trait or two traits.

**--snpDir** The location of the directory containing the SNP likelihood files from the previous step. For example, if the previous step was run with the argument --out "SNPs/7_123456", you should pass "SNPs" to this argument.

**--out** Prefix of the output file.

**--B** Number of boostrap replicates to use for the estimation of standard errors. Must be an integer greater than 1. A larger number will result in less noise in the estimates, at increased computational cost. Default value is 1000. 

**--maxp** A positive number between 0 and 1. Only SNPs with a p-value below this value (for either one or both traits) will be included in the analysis. Default value is 1.

**--makePlot** Can only be used if two traits are being analyzed jointly. If this flag is used, a plot (or plots) illustrating the joint and marginal likelihood functions will be produced.

## This step will produce (in the current working directory)

***out.pdf***  A plot showing contour plots of the joint and marginal selection gradient likelihood surfaces. The contours of the marginal selection gradient likelihood surface are shown in grayscale, with the marginal gradient MLE shown as a black hexagon. The contours of the joint selection gradient likelihood surface are shown in color, with the joint selection MLE shown as a yellow circle. 7 contours are shown for each likelihood surface, corresponding to the contours within which 12.5%, 25%, 37.5%, 50%, 62.5%, 75%, and 87.5% of the probability lies.

Example plot:

<img src="https://github.com/avaughn271/DatePALM/blob/main/output_joint.png">

***out.txt***  A .txt file giving the results of the analysis. The exact format depends on the number of traits being analyzed and the number of different epochs considered.

### 1 Trait - 1 Epoch
```bash
sgrad_mle	sgrad_se	sgrad_z
11.5453	128.1632	0.0901
```
A table with 3 entries is produced, describing one statistic, the selection gradient (abbreviated sgrad). The first entry is the MLE of the selection gradient. The second entry is the standard error of the selection gradient, as estimated by block bootstrapping. If *K* linkage blocks are present in the data and *B* bootstrap samples are taken, then for each of the *B* bootstrap samples, *K* blocks are sampled with replacement from the *K* possible blocks. This new dataset consisting of these *K* sampled blocks appended together is then used to generate a new estimate of the selection gradient. This is done *B* times and the standard error reported here is the standard deviation of this set of selection gradient estimates. The third entry is the z-score, which is simply the estimated selection gradient divided by the standard error. This z-score can be used in a Wald test to generate a (two-sided) p-value to test whether the selection gradient differs significantly from 0 by the simple relation p=2(1-F(|z|)) where F is the CDF of the standard normal distribution.

### 1 Trait - Multiple Epochs
```bash
sgrad_mle	sgrad_se	sgrad_z	delta_mle	delta_se	delta_z
3.6515	55.4455	0.0659	2.0295	5.3768	0.3775
1.6220	50.1745	0.0323	13.5456	1.2497	10.8392
-11.9236	51.4239	-0.2319	-1	-1	-1
```
A table with one row for each epoch is produced (the first row is the most recent epoch, the second row is the next most recent epoch, etc.). Two statistics are analyzed.

* The first column is the MLE of the selection gradient, which is now a vector. The second column is the standard error of each entry of the selection gradient, as estimated by block bootstrapping, as described previously. The third column is the z-score of each entry of the selection gradient.

* The fourth column is a new statistic, called delta. If *n* epochs are analyzed, delta is a vector of length *n-1*, where each entry is simply the difference between the selection gradient in adjacent epochs. For example, for our estimated selection gradient of [3.6515, 1.6220, -11.9236], delta is [3.6515-1.6220, 1.6220-(-11.9236)] = [2.0295, 13.5456]. We also report the standard error of delta, as estimated by block bootstrapping, as described previously, as the fifth column. The sixth column is the z-score of each entry of delta.

As before, these z-scores can be converted to p-values by the simple relation p=2(1-F(|z|)) where F is the CDF of the standard normal distribution. The test being run here is against the null hypothesis that the entries of the selection gradient are the same in adjacent epochs. As delta differs in length from the number of epochs, all delta-related values are set to -1 for the oldest epoch.

### 2 Traits - 1 Epoch
```bash
Trait	jgrad_mle	jgrad_se	jgrad_z	mgrad_mle	mgrad_se	mgrad_z	R_mle	R_se	R_z
MENOPAUSE_AGE	-2.6729	5.1069	-0.5234	-2.9933	1.2111	-2.4716	0.3204	4.1242	0.0777
ProstateCancer	-0.7787	5.1692	-0.1506	-3.2095	3.0028	-1.0688	2.4309	7.0038	0.3471
```

A table with two rows is produced, with each row corresponding to a trait. The first column is the name of each trait. Then, the MLE of 3 different statistics are computed, along with their standard errors (via block bootstrapping) and corresponding z-scores.

* The first statistic analyzed is the joint selection gradient (abbreviated jgrad), which is computed by considering selection acting jointly on both traits.

* The second statistic analyzed is the marginal selection gradient (abbreviated mgrad), which is computed by considering selection acting only on one trait at a time.

* The third statistic analyzed is called R, which is simply the difference between the joint selection gradient and the marginal selection gradient.

With the z-scores computed, you can then test whether an entry of the joint selection gradient is significantly different from 0, an entry of the marginal selection gradient is significantly different from 0, or whether the difference between corresponding entries in the joint and marginal selection gradients is significantly different from 0.

### 2 Traits - Multiple Epochs
```bash
Trait_epoch	jgrad_mle	jgrad_se	jgrad_z	mgrad_mle	mgrad_se	mgrad_z	R_mle	R_se	R_z	delta_mle	delta_se	delta_z	deltaR_mle	deltaR_se	deltaR_z
MENOPAUSE_AGE_1	-2.3350	3.3043	-0.7066	-2.1933	1.1891	-1.8445	-0.1416	2.1555	-0.0657	-4.5906	5.1924	-0.8841	-3.3102	4.4833	-0.7383
MENOPAUSE_AGE_2	2.2557	2.2721	0.9928	-0.9129	0.4732	-1.9293	3.1686	2.7162	1.1665	1.4641	1.3666	1.0713	4.2392	3.1353	1.3521
MENOPAUSE_AGE_3	0.7916	1.7687	0.4476	1.8622	1.5655	1.1895	-1.0706	1.9323	-0.5540	-1	-1	-1	-1	-1	-1
ProstateCancer_1	0.1867	3.0081	0.0621	-0.6201	3.4191	-0.1814	0.8068	5.1952	0.1553	4.0499	4.4015	0.9201	3.7538	6.2655	0.5991
ProstateCancer_2	-3.8632	2.4118	-1.6018	-0.9162	2.2654	-0.4044	-2.9471	1.3780	-2.1386	-5.6693	5.5744	-1.0170	-1.7302	1.6578	-1.0437
ProstateCancer_3	1.8061	4.0297	0.4482	3.0230	4.8738	0.6203	-1.2169	2.2465	-0.5417	-1	-1	-1	-1	-1	-1
```
A table with two rows for each epoch is produced, with each row corresponding to a trait-epoch pair. The first column is the name of each trait-epoch pair. Trait-epoch pairs are sorted first by trait, then by epoch, with the epochs indexed beginning at the most recent epoch (epoch 1) and increasing in number as the epochs get older. Then, the MLE of 5 different statistics are computed, along with their standard errors (via block bootstrapping) and corresponding z-scores.

* The first statistic analyzed is the joint selection gradient (abbreviated jgrad), which is computed by considering selection acting jointly on both traits in each of the different epochs.

* The second statistic analyzed is the marginal selection gradient (abbreviated mgrad), which is computed by considering selection acting only on one trait at a time across each of the different epochs.

* The third statistic analyzed is R, which is the same as defined before.

* The fourth statistic analyzed is delta, which is the difference between the values of the joint selection gradient for a particular trait in adjacent epochs.

* The fifth statistic analyzed is called deltaR, which is the difference between the values of R for a particular trait in adjacent epochs.

All delta and deltaR-related values are set to -1 for the oldest epoch.
