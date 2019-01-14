# FLORA(II): Functional Long-noncoding RNA Assembly Workflow (II)

FLORA(II) provides easy-to-use command line tools for fast lncRNA transcriptome assembly from RNA-seq BAM files.


## Prerequisites

Install bedtools and CPAT.
Please make sure bedtools and cpat.py are on the system path.

## Installing

Download/Clone the repository and change the working directory to the FLORA directory

```
python setup.py install
```


## Run

### generateFilteredBams.py

Generate BAM files with certain regions removed.

```
usage: generateFilteredBams [-h] -g INPUTGTF [-t TYPES [TYPES ...]]
                            [-o OUTPUTDIR] [-n NTHREAD]
                            inputBams

positional arguments:
  inputBams             The path to a txt file that contains paths to input
                        bam files separated by the newline character.

optional arguments:
  -h, --help            show this help message and exit
  -g INPUTGTF           The path to reference gene annotation file in GTF
                        format. It can be downloaded from gencode website.
  -t TYPES [TYPES ...]  The gene types to be removed from BAM files.
  -o OUTPUTDIR          Output directory for the output files
  -n NTHREAD            Number of threads to be used for bedtools intersect.
```

### filterTranscripts.py

Identify lncRNA transcripts from transcriptome assembled by Stringtie or Cufflinks.

```
usage: filterTranscript [-h] -r REFERENCE [-e EXON] [-l LENGTH] [-c CPAT] -x
                        HEXAMER -m LOGIT [-o OUTPUTGTF]
                        inputGTF

positional arguments:
  inputGTF      Input transcriptome file in GTF format. Can be output from
                stringtie or cufflinks

optional arguments:
  -h, --help    show this help message and exit
  -r REFERENCE  The path to reference genome sequences in FASTA format. It
                will be indexed automatically by CPAT if .fai file is not
                present.
  -e EXON       The least number of exons a transcipt should have in order to
                be kept in the final transcriptome
  -l LENGTH     The shortest transcript to be kept in the final transcriptome
  -c CPAT       CPAT cutoff used to filter transcripts
  -x HEXAMER    The path to hexamer table required by CPAT. Can be downloaded
                from CPAT website.
  -m LOGIT      The path to logit model required by CPAT. Can be downloaded
                CPAT website.
  -o OUTPUTGTF  Output prefix for the final transcriptome GTF file
```

### annotateTranscripts.py

Find overlapped and nearby genes of lncRNAs in reference annotation (RefSeq or GENCODE annotation in GFF format).

```
usage: annotateTranscripts [-h] -r REFERENCE -f PATH [-i IDENTIFIER]
                           [-d DISTANCE] [-n NUMBER] [-o OUTPUT]
                           inputGTF

positional arguments:
  inputGTF       Input transcriptome file in GTF format. The file could be
                 output from Cufflinks or StringTie

optional arguments:
  -h, --help     show this help message and exit
  -r REFERENCE   Identify the source of reference annotation. Currently, it
                 should be either refseq or gencode
  -f PATH        The path to the reference annotation file.
  -i IDENTIFIER  The path to the RefSeq assembly report. The file is used to
                 replace RefSeq sequence identifiers with USCS identifiers.
                 Required if you use RefSeq annotation.
  -d DISTANCE    Specify the distance within which the nearby genes should
                 locate. Default: 10,000 bp. Unit: bp.
  -n NUMBER      Specify the number of nearby genes obtained. Default: 1.
  -o OUTPUT      Output file name for the final annotation. Default:
                 FLORA_annotation.txt
```

### functionalPrediction


Gene regulatory netork is constructed via "ARACNe-AP" (https://github.com/califano-lab/ARACNe-AP) based on expression data.
[Command Line]
```
for i in {1..100}
do
  (
    java -Xmx120G -jar .../aracne.jar -e data/expression.matrix.txt  -o output/ --tfs data/lnc_tf.txt --pvalue 1E-8 --seed $i
  ) &
  wait
done

java -Xmx120G -jar .../aracne.jar -o output/ --consolidate

```

LncRNA's function was predicted based on gene regulatory network.
[R]
```
getnetwork()
    usage: getnetwork(lnc.info, coding.info,         #lncRNA and coding genes' information (id and name)
                      network,                       #gene regulatory network (predicted by ARACNe-AP)
                      lnc.name)                      #name of target lncRNA   
    output: coding - target lncRNA  regulatory network

functionalPrediction()
    usage: functionalPrediction(lnc.name,           #name of lncRNA 
                                lnc.coding,         #fcoding - target lncRNA  regulatory network, generated in getnetwork() function
                                gotype)             #["regulator","target","all"] use regulator/target/all genes in the netowrk to do the prediction
    output: list of GO terms


You need to install the R package “gProfileR” first.

```

[Example]
```

lnc.info <- read.table("data/lnc.info.txt", sep = "\t", head=T, row.names=1)
coding.info <- read.table("data/coding.info.txt", sep = "\t", head=T, row.names=1)
network <- read.table("data/network.txt", sep = "\t", head=T)

lnc.name <- "LINC01614"

lnc.coding <- getnetwork(lnc.info, coding.info, network, lnc.name)

results <- functionalPrediction(lnc.name, lnc.coding)
gene.use <- results$gene.use
GO <- subset(results$GO, domain %in% c("BP", "CC", "MF"))

library(ggplot2)
GO$term.name = factor(GO$term.name, levels = GO$term.name)
ggplot(aes(x = term.name, y = -log10(GO$p.value), fill = domain), data = GO) + 
  geom_bar(stat="identity",position="dodge") +
  labs(x = "", y = "-log10_Pvalue", title= lnc.name )  + coord_flip() +
  theme_classic() + scale_fill_manual(values = c('#8470FF','#87CEFA','#FFC125')) +
  theme(legend.title=element_blank(), axis.title.x = element_text(size=9) )

```
![image](https://github.com/shuangat/FLORA/blob/master/data/LINC01614.png?raw=true)

