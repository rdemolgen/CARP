 
![alt text](resources/carp_banner_25pct.png "Chromosomal Anomaly Reports Banner")

## Respository for plotting:
- Copy number changes
- B-allele frequency
- Homozygosity
- UPD

## Dependencies
Dependencies to run CARP tools can be found in the `requirements.txt` file. `pip install -r requirements`.

## CARP Operation
Control script `carp.py` requires the following arguments:
1. `--mode`: `ideogram` or `baf`
2. `--location`: Enter `chr`, `chr:start-end` or `all`
3. `--proband_id`: The index patient sample id
3. `--prefix`: File prefix to output files

Optional arguments use `python -m src.carp --help`

**input/output** \
`-i --iDir`: Carp will look in this directory for required files \
`-o --outDir`: Directory that Carp will output files to \
`-v --vcf`: Provide vcf and index file path \
`--cyto`: Provide cytoband file, default resources/hg38_cytoBand.txt \
`--isca`: Provide ClinGen curated dosage sensitivy file, default resources/

**samples/genotypes** \
`-s --samples`: List of samples to plot, first sample should be proband, e.g. `"WGS_EX1234567 WGS_EX1234568 WGS_EX1234569"`. Sample ids will be retrived from VCF is this option is not provided, assuming proband is the first sample id. \
`-g --genotypes`: List of genotypes in sample order, e.g. `"0/1 0/0 1/1"`. BAF will automatically generate plots based on pre-defined genotypes if genotype option is not given. \
web_ClinGen_region_curation_list_GRCh38_20250425.tsv \

**filtering** \
`-f`: No-filtering. By default BAF script will ignore any variant that is not a PASS. This option will accept all variant quality filters flags. \
`-vq`: Variant quality score threshold, default=30. \
`-dp`: Variant read depth threshold, default=10. \
`-gq`: Variant genotype quality threshold, default=20. \
`-mq`: Variant mapping quality threshold, default=40. \
`-qd`: Variant qual-vy-depth threshold, default=2. 


### DNAnexus
**Build applet** \
To use on DNAnexus build the Carp applet by compiling `carp.wdl` using `java -jar /mnt/data1/software/dxCompiler/dxCompiler-2.11.6.jar compile carp.wdl -extras extras.json -project <project-id>> -folder <folderPath>`

**Applet use**
- The DNAnexus Carp applet will find required files automatically.
- Enter the required arguments for the applet.
- Any optional arguments can be passed using the `options` input field.

### Local
**Direct**
- Create and activate virutal environment and install requirements `/usr/bin/python3 -m venv venv && source venv/bin/activate && pip install --upgrade pip && pip install -r requirements.txt`
- Run Carp: `python3 -m src.carp --mode <mode> --proband_id <sampleId> --location <location> --prefix <prefix> --inDir /data --outDir /data`

**Docker**
Carp can be run through interactive use of a Docker image.
- Build Docker: `docker build --network=host -f Dockerfile -t swglh/carp:<version> .`
- Run Docker: `docker run --rm -itv $(pwd):/data -w /usr/carp swglh/carp:<version>`
- Run Carp: `python3 -m src.carp --mode <mode> --proband_id <sampleId> --location <location> --prefix <prefix> --inDir /data --outDir /data`

### Required files
- VCF file (bgzipped) and associated TBI index
- SavvyCNV CNV calls for each individual. File must be named in the format `cnvs_<sampleID>.20000`
- SavvyCNV data files containing per-bin metrics. File must be named in the format `<sampleID>.coverageBinner.20000.data`
- ISCA regions `https://search.clinicalgenome.org/kb/downloads` *ClinGen_Dosage_Sensitivity.tsv*
- GRCh38 cytoband regions `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz`

### Unit tests
Run tests `python -m unittest`
Coverage
- `coverage run -m unittest -v`
- `coverage report -m` or `coverage html`

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