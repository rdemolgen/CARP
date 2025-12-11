# :tropical_fish:	 CARP :tropical_fish:	
## Chromosome Abmormalities Represented in Python

### Respository for plotting:
- Copy number changes
- B-allele frequency
- Homozygosity
- UPD

#### Dependencies
Dependencies to run CARP tools can be found in the `requirements.txt` file. `pip install -r requirements`.

### Ideograms
Control script `generate_plots.py` for Allele_Fraction and Sample_Dosage classes

There are 3 modes for generating plots:
- baf: directly controls `allele_fraction.py`. Use for generating BAF plots
- ideogram: to generate combined BAF/Dosage ideograms for the whole genome (set `-l` to 'all') or a single chromosome
- dosage: directly controls `savvycnv_dosage.py`. Use for generating coverage plots

**Required files**
- VCF file (bgzipped) and associated TBI index
- SavvyCNV CNV calls for each individual. File must be named in the format `cnvs_<sampleID>.20000`
- SavvyCNV data files containing per-bin metrics. File must be named in the format `<sampleID>.coverageBinner.20000.data`
- ISCA regions

The software looks for these files in a sub-directory named `test_files/`

**Required arguments**
- `-l`: Genomic region to plot `all` will plot across the entire genome. Alternatively, use `<chr:start-end>` for a specific region or `<chr>` for a single chromosome.
- `-fam`: Family number
- `-p`: Proband ID
- `-m`: Mode. Valid values are `baf`, `ideogram`, `dosage`

**Optional arguments**
- `-s`: List of samples to plot, first sample should be proband, e.g. `"WGS_EX1234567 WGS_EX1234568 WGS_EX1234569"`. Sample ids will be retrived from VCF is this option is not provided, assuming proband is the first sample id.
- `-g`: List of genotypes in sample order, e.g. `"0/1 0/0 1/1"`. BAF will automatically generate plots based on pre-defined genotypes if genotype option is not given. 
- `-o`: Give output directory for plots.
- `-f`: No-filtering. By default BAF script will ignore any variant that is not a PASS. This option will accept all variant quality filters flags.
- `-vq`: Variant quality score threshold, default=30.
- `-dp`: Variant read depth threshold, default=10.
- `-gq`: Variant genotype quality threshold, default=20.
- `-mq`: Variant mapping quality threshold, default=40.
- `-qd`: Variant qual-vy-depth threshold, default=2.

**Examples**
#### ideogram mode
1. Generate PDF report containing whole genome ideogram for proband, plus per-chromosome ideograms for all family members, using default parameters
`python generate_plots.py -l all -fam F09999 -p WGS_EX4440000 -m ideogram`
2. Same as 1., but with changes to the SNP filtering
`python generate_plots.py -l all -fam F09999 -vq 60 -gq 20 -mq 60 -dp 10 -qd 10 -p WGS_EX4440000 -m ideogram`
3. Generate PNG ideogram for chromosome 16 only
`python generate_plots.py -l 16 -fam F09999 -vq 60 -gq 20 -mq 60 -dp 10 -qd 10 -p WGS_EX4440000 -m ideogram`
4. Generate PDF report for a single sample from a joint-called family
`python generate_plots.py -l 16 -s WGS_EX4440000 -fam F09999 -vq 60 -gq 20 -mq 60 -dp 10 -qd 10 -p WGS_EX4440000 -m ideogram`

#### baf mode
1. Generate BAF plots for a single sample for a specific genomic region
`python generate_plots.py -l 12:110000000-112000000 -fam F09999  -s WGS_EX4440000 --proband WGS_EX4440000 --mode baf`
2. Generate BAF plots for all samples for a single chromosome
`python generate_plots.py -l 12:110000000-112000000 -fam F09999 --proband WGS_EX4440000 --mode baf`

#### dosage mode
Watch this space

#### Unit tests
Run tests `python -m unittest`

```                                                                                                                                                                                  
                                                   ++                                               
                                                  ++++**                                            
                                                  =+++=+**                                          
                                               ++++++=++++++                                        
   **                                 =+=+==+++++++++++++++*+**                                     
   *******                      +==+++++++++++++++*++******##***##*#*##                             
   ***######*                +++++++++*++*++*+****#*##****++++*+**+*+**####%%%                      
    ##########*            ++++++++****#*#######****+*+++++=++++++*+++++**+#*#%%%                   
    *##########*#*        =+****######**#************+*++=++++*++**=+==+=**+***##%#*                
     **###*#####*####      +++*###************+*+***+*=+========+==*==+=+=*++**#######              
      *#*#*#######**###########***************++++++++++-==-====-==-==-+-+-+=**##%%##%###           
       %#***####****##*##**********++++++*+*+*++++=+=+======-====--=-==--===+**+=*####*#####        
        ***##****************+***+*++++++*++++=+==++========+-==-==-=-===-+=+**===+****##+*#**#     
      ##***##**#*+**++++++++*+++++++++++=+*=+++=++=+==++=+==++=====-=====+-+#=++====*+=*#*****###   
     ##****##*****====-=---=====++==++=++=++==++=+====-==========-==-=-====**==-=++=++=-=+++++++*   
   *********#*****+=-===-=-==-----=====--=--=-===-==-==-==-=--=--=:----=--==-===----==---------+*   
  **********#*#**        #*+==----:::------:--:-=----:-------:--::-:-:-:-----:-=====-===--====      
 **********##            *++++++===------:-----:---:-:--:-::-::---=------=----==++                  
 *******                 ++*++++++=-             ======---==+== +==-----+                           
                         +*****+--              =====+=--       =======                             
                         **+*==-               =-=++==--        ==++=+                              
                                               -=====                                               
                                                ---                                                 
``` 