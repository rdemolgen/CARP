# :tropical_fish:	 CARP :tropical_fish:	
## Chromosome Abmormalities Represented in Python

### Respository for plotting:
- Copy number changes
- B-allele frequency
- Homozygosity
- UPD

#### Dependencies
Dependencies to run CARP tools can be found in the `requirements.txt` file. `pip install -r requirements`.

#### B-allele frequency (BAF)
The BAF script uses a VCF file to generate individual sample BAF plots and/or joint sample BAF plots.

**Usage** \
`python src/allele_fraction.py -v <vcfFile> -l <genomic_location> [options]`

**Required arguments**
- `-v`: Path to vcf file (requires `vcf.gz.tbi` file)
- `-l`: Genomic region to plot `<chr:start-end>` or `chr` for entire chromosome.

**Optional arguments**
- `-s`: List of samples to plot, first sample should be proband, e.g. `"WGS_EX1234567 WGS_EX1234568 WGS_EX1234569"`. Sample ids will be retrived from VCF is this option is not provided, assuming proband is the first sample id.
- `-g`: List of genotypes in sample order, e.g. `"0/1 0/0 1/1"`. BAF will automatically generate plots based on pre-defined genotypes if genotype option is not given. 
- `-o`: Give output directory for plots.
- `-f`: No-filtering. By default BAF script will ignore any variant that is not a PASS. This option will accept all variant quality filters flags.
- `-vq`: Variant quality score threshold, default=30.
- `-dp`: Variant read depth threshold, default=10.
- `-gq`: Variant genotype quality threshold, default=20.

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