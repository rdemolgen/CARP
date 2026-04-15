# :tropical_fish:	 CARP :tropical_fish:	
# Chromosome Abmormalities Represented in Python

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

For optional arguments use `python -m src.carp --help`

### DNAnexus
**Build applet**
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