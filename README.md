# A series of genome-editing plasmids for <I>Saccharomyces cerevisiae</I>.
Features of the plasmid series are:
- A Cas gene and the cognate sgRNA/crRNA gene are encoded on a single plasmid (<I>CEN</I>/<I>URA3</I>).
- The target sequence can be inserted by the Golden Gate Assembly.
- The expression of the Cas gene and/or the sgRNA/crRNA are under the control of the <I>GAL1</I> promoter.
- The plasmid can be eliminated by 5-FOA counter-selection after the completion of genome-editing.

# Available plasmids
|  Number  |  Cas type    |  Promoter for Cas gene  |  Promoter for sgRNA/crRNA  |  Terminator for sgRNA/crRNA  |  Plasmid type  |  Marker gene  |
| -------- | ------------ | ----------------------- | -------------------------- | ---------------------------- | -------------- | ------------- |
|  15-13   |  SpCas9      |  <I> pGAL1 </I>         | <I> pSNR52 </I>            | <I> tSUP4 </I>               | Centromeric    | <I> URA3 </I> |
|  16-15   |  SpCas9      |  <I> pGAL1 </I>         | <I> pGAL1 </I>             | <I> tCYC1 </I>               | Centromeric    | <I> URA3 </I> |
|  16-16   |  enAsCas12a  |  <I> pGAL1 </I>         | <I> pGAL1 </I>             | <I> tCYC1 </I>               | Centromeric    | <I> URA3 </I> |
|  17-31   |  SaCas9      |  <I> pGAL1 </I>         | <I> pGAL1 </I>             | <I> tCYC1 </I>               | Centromeric    | <I> URA3 </I> |

# Scripts for oligo DNA design
Oligo DNA sequences for the Golden Gate Assembly can be designed by a python script for each plasmid.
## How to use the scripts
The input file should be made in TSV (Tab Spaced Value) format. An example shown below:
```
YFG1a [tab] GCTAGTCGATCGATCGTACG
YFG1b [tab] CGTGGTCCCACGCGCGCACC
```
An appropriate python script should be run in a terminal. An example shown below:
```
python sgRNA_oligo_designer_for_(16-15).py input_file.tsv
```
The code above will make an output file (input_file.tsv.output.txt) below:
```
Target name	Target seq	HH + Target seq	Fwd seq for GGA (16-15)	Rev seq for GGA (16-15)
YFG1a	GCTAGTCGATCGATCGTACG	ACTAGCCTGATGAGTCCGTGAGGACGAAACGAGTAAGCTCGTCGCTAGTCGATCGATCGTACG	GGAGACTAGCCTGATGAGTCCGTGAGGACGAAACGAGTAAGCTCGTCGCTAGTCGATCGATCGTACG	AAACCGTACGATCGATCGACTAGCGACGAGCTTACTCGTTTCGTCCTCACGGACTCATCAGGCTAGT
YFG1b	CGTGGTCCCACGCGCGCACC	ACCACGCTGATGAGTCCGTGAGGACGAAACGAGTAAGCTCGTCCGTGGTCCCACGCGCGCACC	GGAGACCACGCTGATGAGTCCGTGAGGACGAAACGAGTAAGCTCGTCCGTGGTCCCACGCGCGCACC	AAACGGTGCGCGCGTGGGACCACGGACGAGCTTACTCGTTTCGTCCTCACGGACTCATCAGCGTGGT
```
You can order the custom oligo DNA synthesis of the two sequences ("Fwd seq for GGA (16-15)" and "Rev seq for GGA (16-15)") for the Golden Gate Assembly.
