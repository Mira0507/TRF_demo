## Utilization of Tandem Repeats Finder (TRF) by Mira Sohn

### Aims

#### - This workflow is aimed at finding and counting repetitive elements based on ensembl gene id using the Tandem Repeats Finder (Unix/Linux ver.)

#### - References: [TRF Web](https://tandem.bu.edu/trf/trf.html), [TRF linux](https://tandem.bu.edu/trf/trf.unix.help.html), [TRF paper](https://pubmed.ncbi.nlm.nih.gov/9862982/), [biomaRt Doc](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)

### Conda environment

#### - This analysis was performed under [conda](https://conda.io/projects/conda/en/latest/index.html) environment (see the [recipe](https://github.com/Mira0507/TRF_demo/blob/master/conda_r.yml))

### Loading packages

    library(data.table)
    library(tidyverse)
    library(biomaRt)

### Importing gene ids

#### - This demonstration uses ensembl human gene id

    # Importing my genes to a vector 
    geneid <- scan("geneid.txt", character(), quote="")

    # Exploring the imported object 
    geneid[1:10]

    ##  [1] "ENSG00000125148" "ENSG00000101361" "ENSG00000162676" "ENSG00000164687"
    ##  [5] "ENSG00000127951" "ENSG00000108691" "ENSG00000171791" "ENSG00000079459"
    ##  [9] "ENSG00000069011" "ENSG00000130164"

    class(geneid)

    ## [1] "character"

    length(geneid)

    ## [1] 11

### Retrieving cDNA sequences by connecting to Ensembl with biomaRt

#### - This step uses useMart() and getSequence() functions from biomaRt

    # listEnsembl() <------ shows the BioMart databases hosted by Ensembl (e.g. snp)

    # listDatasets(useMart('ensembl')) <--- shows available datasets 
    # (e.g.hsapiense_gene_ensembl, hsapiense_snp)  





    # Creating a connection to ensembl human gene db
    ens.mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

    # Exploring Ensembl filters attributes 

    # att <- listAttributes(ens.mart)
    # fil <- listFilters(ens.mart)   # Returns a data frame


    # att[1:20,]
    # fil[1:20,]



    # Creating an empty data frame
    seqTable <- data.frame()


    # Getting the longest cDNA sequences per gene id
    # (This execution can take up to several mimutes!)
    for (i in 1:length(geneid)) {

        # Assigning a gene id
        id <- geneid[i]

        # Getting sequences along with the length
        df <- getSequence(id=id, 
                          type="ensembl_gene_id", 
                          seqType="cdna", 
                          mart=ens.mart, 
                          verbose=F) %>%
                mutate(seqsize=nchar(cdna)) # This new column stores the size of each cDNA

        # Filtering the longest cDNA per GENEID
        # (Multiple cDNA (transcripts) are found in one gene id)
        df <- df[df$seqsize == max(df$seqsize), ] 

        seqTable <- rbind(seqTable, df)
    }

    # Explore the output
    head(seqTable)

    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          cdna
    ## 4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     TCCAAGTCCCAGCGAACCCGCGTGCAACCTGTCCCGACTCTAGCCGCCTCTTCAGCTCGCCATGGATCCCAACTGCTCCTGCGCCGCCGGTAAGAGGCTGGGGATGCCCAGTGTAGACTGTAGCGCTAGAGAAGCAATTTCTGACCCCTCTTTCTTTCTCTGGTCACTCAATTTCAGGACAGGAGTTGCTCCTTCCCAAAGAGTTTTGGGGTATCTTTCTCTCCATTCTAGGTTATTCGGAGCCCCCTTTTTACCGTTAAGGAGATCTGAGTTAATGGCTTGCTCAAGTTCCCAGGAATCGGTTGTGGACTGAGGAACTCGGCCCCGGGCTCTTAGTACGCCGTCCCTTGTTCAGGTATCCAGGGACGGTTCTCACCTCTGTCTTTTCTCCTTGCAGGTGACTCCTGCACCTGCGCCGGCTCCTGCAAATGCAAAGAGTGCAAATGCACCTCCTGCAAGAAAAGTAAGTGGGATCCTCTCTTTCCTCTACCCCTTCCCTGTCCTCCAGCCTGTCCCCTCTCCACCATCCTCAGGGGAATTAAAGCAGTCTGGGGATGCCCCATTGCGCGGAAATTGTTGCCTCCTCAGTGATCCTTATCAGGGAGAGCAGGAATCCTTATTCCCGGTGTCGCTAGTACTCATCTCTGCCGCCTCCTGTCTGCCCCCAGGCTGCTGCTCCTGCTGCCCTGTGGGCTGTGCCAAGTGTGCCCAGGGCTGCATCTGCAAAGGGGCGTCGGACAAGTGCAGCTGCTGCGCCTGATGCTGGGACAGCCCCGCTCCCAGATGTAAAGAACGCGACTTCCACAAACCTGGATTTTTTATGTACAACCCTGACCGTGACCGTTTGCTATATTCCTTTTTCTATGAAATAATGTGAATGATAATAAAACAGCTTTGACTTGA
    ## 6                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        CGAAGGGAACTGAATGAGGACAAGCTGGAGAAGCTGGAGGAGCTGACAATGGATGGGGCCAAGGCTAAGGCTATTCTGGATGCCTCACGGTCCTCCATGGGTCAGTGCAGAGCCTGGCAACCTGCATAAGGTATGGGGCTCTAAATAGCTGGCCTCTTGCATTCACACTTGGTTTTTCCTAGGCATGGACATATCTGCCATTGACTTGATAAACATCGAGAGCTTCTCCAGTCGTGTGGTGTCTTTATCTGAATACCGCCAGAGCCTACACACTTACCTGCGCTCCAAGATGAGCCAAGTAGCCCCCAGCCTGTCAGCCCTAATTGGGGAAGCGGTGCGTCACAGGGGACTCAAAAATGGGAGAATAAGGACTGTTGCCATGTGCACCTGCACTGCTGTATTTCGTGACCCACCATGTCTTCCCTAGTTGTGCTTGATGGGGAGGTGGGGAGCAGGGCTGTCGTGCAACTGGGCAGGTCAGCAGTTCATTTCTCTGACTGCTTCCTTGACTCTCTCTCCAGGTAGGTGCACGTCTCATCGCACATGCTGGCAGCCTCACCAACCTGGCCAAGTATCCAGCATCCACAGTGCAGATCCTTGGGGCTGAAAAGGCCCTGTTCAGGTACCAGTGAGGGCACCTGCCCACAATCAGGTGCCACTTCTGGTGCCCACTGCTTGTTGGGGGATCACGGTGATGGCTGACCAGGGCTCCCTGACCTATACAGGCCTCTGCTATGGGGGTGATGGCCAGTCCTGGTGTCTGAGTGATTCCCAGGGCCCAGCAAAGGGACCAAGTTTCCAGGTCAGCGACATTGGATGCCTTCCCTCTGCCTCTGGGAGCTATGGGTTGGCATGCATTGGGGTAGAGATCCAATCTGGCCTGAGGCTCACTCAGGACTTCGGGGTGAGAGGAGGGGAGGAGCTGAGCTGCCTTGGCTAATGGGGTTGAAATTTCTGATCTTAAACTCTCCACTGAATATTCTCTCAGAGCCCTGAAGACAAGGGGTAACACTCCAAAATATGGACTCATTTTCCACTCCACCTTCATTGGCCGAGCAGCTGCCAAGAACAAAGGCCGCATCTCCCGATACCTGGCAAACAAATGCAGTATTGCCTCACGAATCGATTGCTTCTCTGGTATGGGTGGGGGGGCGTTGGCAGGTGTGAGAAGGGGCTGGGTGGCTGGGTGGGGAGGCTTGCAACCATAGCTTCCACAATGATGGCAATATTTTTCGTCAACAGCAGTTCACCTAGTGAGTGTTGAGACTCTGGGTCTGAGTGAAGCTGAGGGTAGAGGGAACACAGGGTTGGGGTAGTTTCTCTCTTTGGGCTGACAGGCTTTGTCACCCACACACATCCAGAGGTGCCCACGAGTGTATTCGGGGAGAAGCTTCGAGAACAAGTTGAAGAGCGACTGTCCTTCTATGAGACTGGAGAGATACCACGAAAGAATCTGGATGTCATGAAGGAAGCAATGGTTCAGGTCAGTTGGGCTTTGCTGGGTGTGGAGTGGCATAGCTAGCTGTTGGAGGTGATGAACTGTCTGAGCCTGACCTTGTAGAATGGAGGCAAAAAAACTGATTTAATGAGCCTGATCCAATAAAGCCAGAAAGGAGTCCTCAGAGCACCAGAAGTCTTCAGGCCCTTTTAGCACTTTTCTTTGACCAGGCAGAGGAAGCGGCTGCTGAGATTACTAGGAAGCTGGAGAAACAGGAGAAGAAACGCTTAAAGAAGGAAAAGAAACGGCTGGCTGCACTTGCCCTCGCGTCTTCAGAAAACAGCAGTAGTACTCCAGAGGAGTGTGAGGTCAGTAGGCAGCACGGCCCTGGCAGAGATCCTAGGTTGTAGGATTTTCAACAGCAGAACAAAGGATATGCTGCATCAAGCTGTGGTCTTGAGTCCAGGCTTTTGGACTGAAACAAGGACCTGAAACATCTAAAACTACCTCTTGATTCTATAGGAAGGAGATAGGTGCTGAACTTGCTCAAGAGCCCAGAGAGCTGGTTGTAGCTCACACCCGTTCCCTGGGCATGTGTGTTCTGTCCTCGGCTGCCTCCCAGGAGTCCTCAACCTGGGGTAGTGTAAATTCCTGCTCTGCTTATTATCAGACGTGTGTCCGGAGGTGGTCGTGTTTCACAGTGGGGATGGGGGCAGGGAGGTCCCCAATGTGCTAAGCTACAATCATTCTCCCTGAGATTTTCATTTAGCACCCAGTTTCTTAAACAGTGTTTCAGGGCCCTGTCTGGAACTTGGCATGATGGTTCTGTTGCGACCAGCATGGTGGGTGTTTTTTAGGTTTTTTTTTTTAATGGGCTGAGGTAATTTCTCATGACATGTTTTCCTTCTAATTTGGGACAGCCTTTGGGGTGGATTTCTAAAGTTATACCCACACAATTAAACTATCCCAGAAACACTGGGCAATGTTAACGACACGCGTTCCCCTGCCTTGGCTACTTAATTGCTGAAGATGTAATGAGCACTGTTCTCACAGCCTGTTCCCCTGTCCTTCCCTTTAGGAGATGAGTGAAAAACCCAAAAAGAAGAAAAAGCAAAAGCCCCAGGAGGTTCCTCAGGAGAATGGAATGGAAGACCCATCTATCTCTTTCTCCAAACCCAAGAAAAAGAAATCTTTTTCCAAGGAGGAGTTGATGAGTAGCGATCTTGAAGAGACCGCTGGCAGCACCAGTATTCCCAAGAGGAAGAAGTCTACACCCAAGGAGGAAACAGTTAATGACCCTGAGGAGGCAGGCCACAGAAGTGGCTCCAAGAAAAAGAGGAAATTCTCCAAAGAGGAGCCGGTCAGCAGTGGGCCTGAAGAGGCGGTTGGCAAGAGCAGCTCCAAGAAGAAGAAAAAGTTCCATAAAGCATCCCAGGAAGATTAGAATGCAAATGGACATTCTCTGGGAGGTGGGGCATACCATAGCCCAAGGTGACATTTCCCACCCTGTGCCGTGTTCCCCAATAAAAACAAATTCACAAG
    ## 3  AGTCGAAGGCTGAGCGGACCCGGTCGAGAGCCGGCGAGCTCCGCGCAGGGAGGGTGCGCCCACCGGTCCCGCCGGGCGCCCGCGGGACGCGCCGCCAGGGCCCTCTCCGCCGGGGGCTCGGCGCTCGCCCACCTCTTCCAAATTTAACCATTACCTAAATCCGAAGGGAAATGAGCAAACCTCTCGGATTGGGTGTCAAGGTCTCCTCCGGGCTGGGGCTGAGCAAGCCCTCGGAGTGACCGTGGGTGACAGCGGCTCCAGGGACTCTTGGGGCGCAGTGGGGAAAGTGCCGGACCACCATGCCGCGCTCATTTCTCGTCAAAAGCAAGAAGGCTCACAGCTACCACCAGCCGCGCTCCCCAGGACCAGACTATTCCCTCCGTTTAGAGAATGTACCGGCGCCTAGCCGAGCAGACAGCACTTCAAATGCAGGCGGGGCGAAGGCGGAGCCCCGGGACCGTTTGTCCCCCGAATCGCAGCTGACCGAAGCCCCAGACAGAGCCTCCGCATCCCCAGACAGCTGCGAAGGCAGCGTCTGCGAACGGAGCTCGGAGTTTGAGGACTTCTGGAGGCCCCCGTCACCCTCCGCGTCTCCAGCCTCGGAGAAGTCAATGTGCCCATCGCTGGACGAAGCCCAGCCCTTCCCCCTGCCTTTCAAACCGTACTCATGGAGCGGCCTGGCGGGTTCTGACCTGCGGCACCTGGTGCAGAGCTACCGACCGTGTGGGGCCCTGGAGCGTGGCGCTGGCCTGGGCCTCTTCTGCGAACCCGCCCCGGAGCCTGGCCACCCGGCCGCGCTGTACGGCCCGAAGCGGGCTGCCGGCGGCGCGGGGGCCGGGGCGCCAGGGAGCTGCAGCGCAGGGGCCGGTGCCACCGCTGGCCCTGGCCTAGGGCTCTACGGCGACTTCGGGTCTGCGGCAGCCGGGCTGTATGAGAGGCCCACGGCAGCGGCGGGCTTGCTGTACCCCGAGCGTGGCCACGGGCTGCACGCAGACAAGGGCGCTGGCGTCAAGGTGGAGTCGGAGCTGCTGTGCACCCGCCTGCTGCTGGGCGGCGGCTCCTACAAGTGCATCAAGTGCAGCAAGGTGTTCTCCACGCCGCACGGGCTCGAGGTGCACGTGCGCAGGTCCCACAGCGGTACCAGACCCTTTGCCTGCGAGATGTGCGGCAAGACCTTCGGGCACGCGGTGAGCCTGGAGCAGCACAAAGCCGTGCACTCGCAGGAACGGAGCTTTGACTGTAAGATCTGTGGGAAGAGCTTCAAGAGGTCATCCACACTGTCCACACACCTGCTTATCCACTCAGACACTCGGCCCTACCCCTGTCAGTACTGTGGCAAGAGGTTCCACCAGAAGTCAGACATGAAGAAACACACTTTCATCCACACTGGTGAGAAGCCTCACAAGTGCCAGGTGTGCGGCAAGGCATTCAGCCAGAGCTCCAACCTCATCACCCACAGCCGCAAACACACAGGCTTCAAGCCCTTCGGCTGCGACCTCTGTGGGAAGGGTTTCCAGAGGAAGGTGGACCTCCGAAGGCACCGGGAGACGCAGCATGGGCTCAAATGAGCACCCTGGCTGGCTGCAAGCAGCAGCTACACAACACTACAGAGGGCAGCCTCCCTGCTTGCCACCACTCTGCTCCCTGCTTGCCTCCACTCCCTTCTGACTTTCCAGACCCCAGGTCCAGTCTGCAGATCCTACCAGGTTGCTCCTCCTTCGCCTTACCTCCTGGAGCTGCCAGAAGAAATGAGGTACCTTTTCAAAGTGCAGCCGAGAGTGAGAACCAAGTGACTCTCTAGGCTTCGGACACAAATAGGCTCCTCTACACCTGAAGACAAAGGCAAAGTCAAATGGGGACCAGAATAAATCTTAGACCCCACAGTCCTTCCCATTTCCAGCCCTAATCTACAGACAGGAATGCCCTTCAGGTTTCTTCCCTCCCCCCTCTTGACCTACCCCAGATATTTGTGTGGAAGAGGAGGAATCACCATTTACAAGGTGGACAAATGCTAATATTTTTATCTAGAAAGAAGAGTGAGTGTTAACTTTTATTTTTTTCCTTCTGGGGGGTCTGTTGACTCCTTTCTTTTGGGTGCTGCCTATAAATCTTGGAGGAATCATTTCTCCTCCTCAAAAACTGATTCAGAAACTGACTTGGGGAAGGAATTTAATACTTTGAAGTCATGAGATGCACCATCGAGGCTACCCCCAAGAAGAAGCAGAAGAGAAGTTGGTAATGAGAGGGGATTAGAGGTCCTCCCTTCAGTAGGGCTGTGAAAACCTCATCACTGGAGGTAAAAGCACAAGCAATGCCTGTGGACAAGATGTCATTCATTCACTCAGCAAATGTTCATGGATCACCGGCTACCAAGGTACCAGGCACCATGCTAGGTATTGGGGAAGAGAGACTGAAGTCACAACCCCTGACTGCTCCTCAAAAGCTAACGGTTGCACCTCCAAGTGGCTGGGTCTGTTCTTACTCTTGGAGGGAATTCTGAGAAGACAGCACAGAATTGTAAACCTTCCCTTTTGACCCTTTTGGATTTTATCAGGTGTAAACAAAAAGCTGAACAGTTACTTCAAAGATATGTGTGTATATTCAGTTTTTTATTGTTAAGCTGATATTTTAAAGATTTCTGAGCTAGCAGGCATGTGGGAAGGAAGGCTCTGTCTTCAACTCTTTGACCCTCCATGTGTACCATAGAGGGGGGAAAGGTGGTATTTTCACTTTGATGAGGTTGGTAAATGTTTTTAGATCTTCTGGTAAGCATTATGTTTGTTAATACATATTTATTAGAGTGATGTTTTAAGTTAATAAAGTATTAAGAGTATTACAGATTGCCTTTCCTTTAGGCAGGGAATTTGGGGAAGACCTTTGTTTTAGTCTGTAGGGGGAAAGCTGGCACTGAAGGAAACTGCAATGGAGTATAAAAGTGAAGCACACAAGCTAATATCCAGAAGAGTAAAATCCAATTGCTCGGGTTTCCATATACTTTCTCCTTTGCTGCTACACCTCAGAGATACTTCTTTCATCTTCTTTCATAAGAGAAAAACTAGGAGGTGTACAAAACTAAGCTACTACCTTGTAGAGCCCCTGATGGTTCCACAGGCACGTAAGTAAGCTCTTTCCTCTGAAAATGTCAACATATTTTTTTAGGGGCCTGAAGGCTCCCAGATGAATAAAGAAATCACAGGTCTTAGGTCTAGGTTTTAGTTTTAGTCACAGCTCTGTCACCAACACATCATGTGACTTGTGTAAATCTTTTGACTTCCCCGATACTGGTTTCCTTGATCAAGTAATACGGAGACAGCTCTGTCATTCCTGTCCAACCCTGACCGCTTTGATTCTGTGGAGACTGGAAAGGAAAAAAAGTGTGAGATTGGTTAGCTGCCTCCCCTAATGCAAGTCACAAATCTAAAGGCCCATGAGGATCTGGGTTATTTCTGGTTCTGTCCCACCTGTGGGAGTTTTTTGCATAAGTCAATCCTTTAAGAGAATGAAACCAGAATAAAATCCTGGCCAAAATCCCTACTACTAGTGCTTTAACACCTCCCAGCCAGACATTCAGGGTATGCATGCCCCATGGTGAGAACCAGCAGATTGCCCATTTCCAACTCAAATTGATTTCACACATATAATTGAGCACCCAGTATGTGCTGGGCATTCTTCTAGGGGGTGGAGATATAGCAGTGAGGACAATGTCTTACGAAGCTTTATTTATTCATTGTTGAATCCTCAGTTCCCAGAATAGTGCCTGGCACACAGTAGGAGATTAATAAATATTTACTGAGATGTTAAGTGAATAAACAGAAAATAAGCTCCTGCTTCATCCTAGGGCCTGGCTCTGCTCTCTACAGCACAAAGGGAAGAGAGATTGCTTGACAGCCTGGGAAACAAATAGAATCCACACGACTTTCTCTCCAGGTTGGCAGGAAAGGGGGGCTTGAATGCTGCTTTTCTCTGGTTTCTCCAAATCTCAGTCCTTGGCTGTCCTCACCTGATTTCCAGCACCATGTCTTTCCATGTCAGTTTACTCTCATTCAGAATTTAAATTCGTTTACTCTGTAATGGAGTTTAAGATAGTTTTTCAGAACTGATGGATGATGAACAGTTGCAGCGTCCAGGAACCACTTCCTAGTGTCCTGGGAGACAATCCTCCCAGTGAATGAAAATGACAAAAACATCCTTCCCACAGAAGTCACCCTGTGTACTGGTTGGTTATTTCCAACCTGGTTTGAAAAGTAAAGAAGTTGCACTGAGGGGTAAGTAGAACTGATTGAAATCTAGGTGCCTGTGCTGTCCTGCTCTTCCCTTATGACGTTTAGAACATTTGGTCCATTTTTGTTATCCCTGTTTATGCTCTAAGTGAGGAATTGTAATTTTGATGTTAACTGTATGCCTTTTTTTTTTTTTTTTGGTAAAAAATGTATCTTGATATCCCCAGTTCAGAATCAAAACAAATTACTGTGCCTGTTACCGTGTGAAATTTGAAGATTGCAAATATCTCTCTTAATGTTATTACTTGACTCCATTAAAAATAGATTGTGAATA
    ## 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  GGGAGCGTGCGCGCCTCTTGCCCGCCCGCGGGCCGCAAGATTCCGGAGGTGGTCCACCCCGTGGTCCACCTACCTCTGCTTCTTTCCCTTCGCCCAAACGGCAGCCCTTCCGCATTGCTTCCTGTCCTTTAGCGCGCGCACCCGTCACCTCACGCTGCACTTCTTTCGACCCCCTCCAGGCGACCCTGTATTTCCCTTTTTTCCCCCTTTACTCATCCTTCCCCTTCCTGCCACCATTCTGTCTCTCCCATTTACCCCATCGTGGCAGCACCCACTCCCCTTTCCCTCCCTGTCGCATCTTTTGTTCCCTGTGGCGCGCAGGTCACCCTCCGTTTTCTTCATGGGGACGCGGTGCTGGCGCGCAGTTTCCCGCAGAAATCCTGTGGAGGGGTCTGAGTGGCGCTACAGCTTCAGCTTGCATTCCTCTGTCAGCCCCGTCTCCCCCAGGAGTGGGAATAGCTTTGCGAAAAATGGGCGCAATGGCCAAGCCAGATTGTATCATCACTTGTGATGGTAAAAACCTCACCATAAAAACTGAGAGCACTTTGAAAACAACACAGTTTTCTTGTACCCTGGGAGAGAAGTTTGAAGAAACCACAGCTGATGGCAGAAAAACTCAGACTGTCTGCAACTTTACAGATGGTGCATTGGTTCAGCATCAGGAGTGGGATGGGAAGGAAAGCACAATAACAAGAAAATTGAAAGATGGGAAATTAGTGGTGGAGTGTGTCATGAACAATGTCACCTGTACTCGGATCTATGAAAAAGTAGAATAAAAATTCCATCATCACTTTGGACAGGAGTTAATTAAGAGAATGACCAAGCTCAGTTCAATGAGCAAATCTCCATACTGTTTCTTTCTTTTTTTTTTCATTACTGTGTTCAATTATCTTTATCATAAACATTTTACATGCAGCTATTTCAAAGTGTGTTGGATTAATTAGGATCATCCCTTTGGTTAATAAATAAATGTGTTTGTGCTAATA
    ## 11                                                                                                                                                                                                                                                                                                           CTGCTGGGGTGAGCAGCACTGTAAAGATGAAGCTGGCTAACTGGTACTGGCTGAGCTCAGCTGTTCTTGCCACTTACGGTTTTTTGGTTGTGGCAAACAATGAAACAGAGGAAATTAAAGATGAAAGAGCAAAGGATGTCTGCCCAGTGAGACTAGAAAGCAGAGGGAAATGCGAAGAGGCAGGGGAGTGCCCCTACCAGGTAAGCCTGCCCCCCTTGACTATTCAGCTCCCGAAGCAATTCAGCAGGATCGAGGAGGTGTTCAAAGAAGTCCAAAACCTCAAGGAAATCGTAAATAGTCTAAAGAAATCTTGCCAAGACTGCAAGCTGCAGGCTGATGACAACGGAGACCCAGGCAGAAACGGACTGTTGTTACCCAGTACAGGAGCCCCGGGAGAGGTTGGTGATAACAGAGTTAGAGAATTAGAGAGTGAGGTTAACAAGCTGTCCTCTGAGCTAAAGAATGCCAAAGAGGAGATCAATGTACTTCATGGTCGCCTGGAGAAGCTGAATCTTGTAAATATGAACAACATAGAAAATTATGTTGACAGCAAAGTGGCAAATCTAACATTTGTTGTCAATAGTTTGGATGGCAAATGTTCAAAGTGTCCCAGCCAAGAACAAATACAGTCACGTCCAGTTCAACATCTAATATATAAAGATTGCTCTGACTACTACGCAATAGGCAAAAGAAGCAGTGAGACCTACAGAGTTACACCTGATCCCAAAAATAGTAGCTTTGAAGTTTACTGTGACATGGAGACCATGGGGGGAGGCTGGACAGTGCTGCAGGCACGTCTCGATGGGAGCACCAACTTCACCAGAACATGGCAAGACTACAAAGCAGGCTTTGGAAACCTCAGAAGGGAATTTTGGCTGGGGAACGATAAAATTCATCTTCTGACCAAGAGTAAGGAAATGATTCTGAGAATAGATCTTGAAGACTTTAATGGTGTCGAACTATATGCCTTGTATGATCAGTTTTATGTGGCTAATGAGTTTCTCAAATATCGTTTACACGTTGGTAACTATAATGGCACAGCTGGAGATGCATTACGTTTCAACAAACATTACAACCACGATCTGAAGTTTTTCACCACTCCAGATAAAGACAATGATCGATATCCTTCTGGGAACTGTGGGCTGTACTACAGTTCAGGCTGGTGGTTTGATGCATGTCTTTCTGCAAACTTAAATGGCAAATATTATCACCAAAAATACAGAGGTGTCCGTAATGGGATTTTCTGGGGTACCTGGCCTGGTGTAAGTGAGGCACACCCTGGTGGCTACAAGTCCTCCTTCAAAGAGGCTAAGATGATGATCAGACCCAAGCACTTTAAGCCATAAATCACTCTGTTCATTCCTCCAGGTATTCGTTATCTAATAGGGCAATTAATTCCTTCAGCACTTTAGAATATGCCTTGTTTCATATTTTTCATAGCTAAAAAATGATGTCTGACGGCTAGGTTCTTATGCTACACAGCATTTGAAATAAAGCTGAAAAACAATGCATTTTAAAGGAGTCCTTTGTTGTTATGCTGTTATCCAATGAACACTTGCAAGCAATTAGCAATATTGAGAATTATACATTAGATTTACAATTCTTTTAATTTCTATTGAAACTTTTTCTATTGCTTGTATTACTTGCTGTATTTAAAAAATAATTGTTGGCTGGGTGTGGTAGCTCACGCCTGTAATCCCAGCACTTTGGAATGTCAAGGCAGGCAGATCACTTGAGGTCAGGAGTTTGAGACCAGCCTGGCCAAACATGTGAAACGCTGTCTCTATTAAAAATACAAAAATTAGCCGGGCATGGTGGTACATGCCTGTAATCCTAGCTACTTGGGAGGCTGAGGCAGGAGAATCGCTTGAACCTGAGAGGAAGAGGTTGCAGTGAGCCAAGACTGAGCCACTGCACTCCAGCATGGGTGACAGAGAAAACTCTGTCTCAAACAAAAAAATAATAAAATTTATTCAGTAGGCTGGATTCTACACAAAGTAATCTGTATTTGGGCCATGATTTAAGCACATCTGAAGGTATATCACTCTTTTCAGGCTATAATTATTTGGGTAATCTTCATTCTGAGACAAACTTAATCTATATCATTTACTTTGCAACAGAACAACCCTACAGCATTTTGGTTCCCAGACTAAGGGAACTAATATCTATATAATTAAACTTGTTCATTTATCATTCATGAAATATAAAATACTTGTCATTTAAACCGTTTAAAAATGTGGTAGCATAATGTCACCCCAAAAAGCATTCAGAAAGCAATGTAACTGTGAAGACCAGGGTTTAAAGGTAATTCATTTATAGTTTATAACTCCTTAGATGTTTGATGTTGAAAACTGCTTTAACATGAAAATTATCTTCCTCTGCTCTGTGTGAACAATAGCTTTTAATTTAAGATTGCTCACTACTGTACTAGACTACTGGTAGGTTTTTTTGGGGGGGGTGGGTAGGGATATGTGGGTAATGAAGCATTTACTTACAGGCTATCATACTCTGAGGCCAATTTTATCTCCAAAGCAATAATATCATTAAGTGATTCACTTCATAGAAGGCTAAGTTTCTCTAGGACAGATAGAAAACATGAATTTTGAAATATATAGAACAGTAGTTAAAATACTATATATTTCAACCCTGGCTGGTAGATTGCTTATTTTACTATCAGAAACTAAAAGATAGATTTTTACCCAAACAGAAGTATCTGTAATTTTTATAATTCATCAATTCTGGAATGCTATATATAATATTTAAAAGACTTTTTAAATGTGTTTAATTTCATCATCGTAAAAAGGGATCATCTCAGAGAGAACAGCAGTATTCTGCGTATTTTTAAAAATGCTCTAGAGTAACATTTGAAGTAATTCACTGTAGTGTATGCCAGTCCTAGAAATAATTTTTTTAATTTCTGGTGTCTGTTTCTAATACACTAACCAAGTTTTCAAAATATATTTACAAAGATGCATCTTTACCCATTATTTTAAAATGATTAAGGAGGATAGTTGCTTCAGGTAACAAGCTAATTTTTCAAATATTAGGCCCTTACAGAACTATTTAGTCAAAAAGTAAGATATTCCTTTAAAATATATAACCCAAAGCTTTCAGTTAAACATGATATATCACAAATACTATTAAAATGTTAAAGAGAAATGCAAATAGCATTAAATGATGACCAAAATGTAAAATATTGTAGATTTCAAAAGCTGTGTCTCTATTAGGTGGGATACCAAATGTAAATGATGTAACTGACGTTGTTTTTTACTTTTTACTTTTTAAAAAAGACTAAAAACGTTTTGATATTATACAATGTATTTGTTTCAGATAAGGTCATTGTCATTTAGTATATATAATTAATATATGTACAAGTTTAAGTAAATTCCTGTGAGTAAAAATGGACTTATCACAAAACATAGTTCTAAAGAAAGGTATATGCTCATATACACGGTGTCCATTAATTTAATGGGAACTAGGTATAACTTCAGGAGAATTTGGCAAATAATTCATTAATCCATGTAAATATTCAAAAGCTTGTTCTATCCACATTATTTCAAGGGATCACTTTATTTTTCATTATACTTTCACAGCACTTTTCTAGTAAATTCTGTAACACAGAAATTCCATTTTGGAATCATTTCATGTTACCAATAATTCTAGACTTTTATAACATTTAACATGTTGATGGAAATAGATTACATCTGCTAGAACCTTTTGCCTTAACTATTCACCAATATATGCTAATATTCATAAATATGGATTGACTGTTTACAAACATTAGAATCTTGTCTTGGTTCCATTTTGATGGCTAATATTTGTTATCTTAATTAAGACTATTTCTGAGGTCATGATTACTTGAAAATATTGACTAAAACTGGGTCCTTAGAAATTCCAGGTGGAGCTGATTTACCTATGACTGAGGGGAAAAAAAAATCAAATTTTACTGATAATAGTAATGCTCCAAATGAATTAATGACACATCTGTTCAATAAATAAAGAGCTTAAATATACAAAACATAAGAAATCTGGGCAACAAAACTTGTGGTCTTTACTTTTGAATAGCTACCCAAGAAAAGGTTTTAAAGGTAAAAGTTATGAGTAATGTCATCACAATAAGCTCTTGTTTAAAATTCTTTTCTTTTATGTATAATTAGGTTTATGTTTCATGTCTTTTTAAAACCTTATAAAAGATTTAATTATCACATCTATTCTTCAATGTGGAAATATTAAATATTGTTGGTTGTAAAATAATA
    ## 12                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       GGTTTATTTTCCAGATGCAATCAATGCCCCAGTCACCTGCTGTTATAACTTCACCAATAGGAAGATCTCAGTGCAGAGGCTCGCGAGCTATAGAAGAATCACCAGCAGCAAGTGTCCCAAAGAAGCTGTGATGTGAGTTCAGCACACCAACCTTCCCTGGCCTGAAGTTCTTCCTTGTGGAGCAAGGGACAAGCCTCATAAACCTAGAGTCAGAGAGTGCACTATTTAACTTAATGTACAAAGGTTCCCAATGGGAAAACTGAGGCACCAAGGGAAAAAGTGAACCCCAACATCACTCTCCACCTGGGTGCCTATTCAGAACACCCCAATTTCTTTAGCTTGAAGTCAGGATGGCTCCACCTGGACACCTATAGGAGCAGTTTGCCCTGGGTTCCCTCCTTCCACCTGCGTTCCTCCTCTAGCTCCCATGGCAGCCCTTTGGTGCAGAATGGGCTGCACTTCTAGACCAAAACTGCAAAGGAACTTCATCTAACTCTGTCCTCCCTCCCCACAGCTTCAAGACCATTGTGGCCAAGGAGATCTGTGCTGACCCCAAGCAGAAGTGGGTTCAGGATTCCATGGACCACCTGGACAAGCAAACCCAAACTCCGAAGACTTGAACACTCACTCCACAACCCAAGAATCTGCAGCTAACTTATTTTCCCCTAGCTTTCCCCAGACACCCTGTTTTATTTTATTATAATGAATTTTGTTTGTTGATGTGAAACATTATGCCTTAAGTAATGTTAATTCTTATTTAAGTTATTGATGTTTTAAGTTTATCTTTCATGGTACTAGTGTTTTTTAGATACAGAGACTTGGGGAAATTGCTTTTCCTCTTGAACCACAGTTCTACCCCTGGGATGTTTTGAGGGTCTTTGCAAGAATCATTAATACAAAGAATTTTTTTTAACATTCCAATGCATTGCTAAAATATTATTGTGGAAATGAATATTTTGTAACTATTACACCAAATAAATATATTTTTGTACAAAA
    ##    ensembl_gene_id seqsize
    ## 4  ENSG00000125148     903
    ## 6  ENSG00000101361    2980
    ## 3  ENSG00000162676    4554
    ## 1  ENSG00000164687     986
    ## 11 ENSG00000127951    4256
    ## 12 ENSG00000108691     996

    dim(seqTable)

    ## [1] 11  3

### Running TRF

#### - This workflow uses Linux TRF

#### - The input file format has to be **FASTA** (see below example)

#### - The output files are either html and/or .dat files depending on flags set

#### - This workflow uses recommended parameters (see [Unix/linux TRF](https://tandem.bu.edu/trf/trf.unix.help.html))

    # Below is an example of FASTA format
    #
    #
    # >ENSG00000125148
    # GGTTTATTTTCCAGATGCAATCAATGCCCC....
    # >ENSG00000101361
    # CACTTCTAGACCAAAACTGCAAAGGAAC......

    # Creating an input fasta file storing gene ids and cDNA sequences (as fasta format!!)
    for (i in 1:length(geneid)) {

        # Adding a line for my gene id (e.g. >ENSG00000125148)
        line1 <- paste0(">", geneid[i]) 

        # Adding a line for cDNA sequences 
        line2 <- seqTable$cdna[i]

        write(line1, file="trf_input.fa", append=T)
        write(line2, file="trf_input.fa", append=T) # trf_input.fa is my input for trf

    }


    # Running tandem repeats finder (trf)
    # Command: 
    # trf <File> <Match> <Mismatch> <Delta> <PM> <PI> <Minscore> <MaxPeriod>
    # -h: avoids to create html output
    # -d: creates a .dat file 
    system("trf trf_input.fa 2 7 7 80 10 50 500 -d -h")  # Returns a .dat file 

    # Done.> pops up when it's done  
    # "trf_input.fa.2.7.7.80.10.50.500.dat" (=output) is created in my working directory 

### Exploring analyzed repetitive elements

    # Converting the TRF results from the .dat file to a data frame
    # 
    # 
    # Isolating every single line of info to the trflines object
    trflines <- readLines("trf_input.fa.2.7.7.80.10.50.500.dat")  

    # Exploring the trflines object
    class(trflines)

    ## [1] "character"

    length(trflines)

    ## [1] 109

    trflines[1:20]

    ##  [1] "Tandem Repeats Finder Program written by:"
    ##  [2] ""                                         
    ##  [3] "Gary Benson"                              
    ##  [4] "Program in Bioinformatics"                
    ##  [5] "Boston University"                        
    ##  [6] "Version 4.09"                             
    ##  [7] ""                                         
    ##  [8] ""                                         
    ##  [9] "Sequence: ENSG00000125148"                
    ## [10] ""                                         
    ## [11] ""                                         
    ## [12] ""                                         
    ## [13] "Parameters: 2 7 7 80 10 50 500"           
    ## [14] ""                                         
    ## [15] ""                                         
    ## [16] ""                                         
    ## [17] ""                                         
    ## [18] "Sequence: ENSG00000101361"                
    ## [19] ""                                         
    ## [20] ""

    # Creating an empty data frame
    trfTable <- data.frame()

    # Storing every single line to the data frame except empty ("") and 
    # unnecessary (1 to 6) lines 
    # (The line 1 to 6 of trflines have title & description about this run)
    for (i in 7:length(trflines)) {

        if (trflines[i] != "" & !grepl("Parameters", trflines[i], fixed=T)) {

            trflines[i] <- str_replace(trflines[i], "Sequence: ", "")

            df <- data.frame(Line=trflines[i])

            trfTable <- rbind(trfTable, df)

        }

    }



    # Exploring the output data frame
    head(trfTable)

    ##                                                                                                            Line
    ## 1                                                                                               ENSG00000125148
    ## 2                                                                                               ENSG00000101361
    ## 3                                                                                               ENSG00000162676
    ## 4 1619 1660 22 1.9 22 95 0 75 7 52 11 28 1.64 CTCCCTGCTTGCCACCACTCTG CTCCCTGCTTGCCACCACTCTGCTCCCTGCTTGCCTCCACTC
    ## 5                                                                                               ENSG00000164687
    ## 6                                                                                               ENSG00000127951

    dim(trfTable)

    ## [1] 15  1

    class(trfTable)

    ## [1] "data.frame"

    # Clean the TRF data frame 
    # (removing rows containing parameters and adding line codes)
    trfTable <- trfTable %>%
        mutate(code=ifelse(grepl("ENSG", Line, fixed=T), 
                           "ID", 
                           "CDNA")) 


    # Exploring the output data frame
    head(trfTable)

    ##                                                                                                            Line
    ## 1                                                                                               ENSG00000125148
    ## 2                                                                                               ENSG00000101361
    ## 3                                                                                               ENSG00000162676
    ## 4 1619 1660 22 1.9 22 95 0 75 7 52 11 28 1.64 CTCCCTGCTTGCCACCACTCTG CTCCCTGCTTGCCACCACTCTGCTCCCTGCTTGCCTCCACTC
    ## 5                                                                                               ENSG00000164687
    ## 6                                                                                               ENSG00000127951
    ##   code
    ## 1   ID
    ## 2   ID
    ## 3   ID
    ## 4 CDNA
    ## 5   ID
    ## 6   ID

    dim(trfTable)

    ## [1] 15  2

    # Extracting the number of repetitive elements per GENEID
    #
    # Creating an empty data frame 
    trfnumTable <- data.frame()
    code.vector <- trfTable$code

    GENE <- c()
    COUNT <- c()

    index.id <- 0
    index.cdna <- 0
    for (i in 1:nrow(trfTable)) {

        if (code.vector[i] == "ID") {

            index.id <- index.id + 1
            index.cdna <- 0

        } else {

            index.id <- 0
            index.cdna <- index.cdna + 1
        }

        GENE[i] <- index.id # Indicates GENEID coded "ID"
        COUNT[i] <- index.cdna  # Indicates repetitive elements coded "CDNA" 


    }
    # Exploring the GENE and Count vectors
    GENE

    ##  [1] 1 2 3 0 1 2 3 4 0 0 1 2 0 1 2

    COUNT

    ##  [1] 0 0 0 1 0 0 0 0 1 2 0 0 1 0 0

    # GENE[i] = 0: repetitive element line 
    # GENE[i] > 0: gene id line
    #   1  2  3  0  1  2  3  4  0  0  1
    #         |  |           |  |  |
    #     geneid |       geneid |  |
    #           repeat1     repeat1|
    #                            repeat2
    #
    #
    #
    # COUNT[i] = 0: gene id line
    # COUNT[i] > 0: repetitive element line 
    #   0  0  0  1  0  0  0  0  1  2  0 
    #         |  |           |  |  |  
    #     geneid |       geneid |  |
    #          repeat1      repeat1|
    #                           repeat2


    # Assign vectors storing line indices of 
    # gene id (repetitive.geneid), the first repetitive elements (repetitive.first), 
    # and the last repetitive elements (repetitive.last)
    repetitive.geneid <- which(COUNT == 1) - 1  # right before the first repetitive elements
    repetitive.first <- which(COUNT == 1)  # the first repetitive elements 
    repetitive.last <- which(GENE == 1) - 1  # the last repetitive elements per gene id


    # Exploring the indices
    repetitive.geneid

    ## [1]  3  8 12

    repetitive.first

    ## [1]  4  9 13

    repetitive.last

    ## [1]  0  4 10 13

    # In case that the first and/or the last gene had one repetitive element: 
    if (repetitive.last[1] < repetitive.first[1]) {

        repetitive.last <- repetitive.last[-1]

        if (repetitive.last[length(repetitive.last)] < 
            repetitive.first[length(repetitive.first)]) {

            repetitive.last <- c(repetitive.last, 
                                 repetitive.first[length(repetitive.first)])

        }

        

    }

    # Explore the indices (all three vectors have to have the same length!)
    repetitive.geneid

    ## [1]  3  8 12

    repetitive.first

    ## [1]  4  9 13

    repetitive.last

    ## [1]  4 10 13

    # Summarize the gene id (GENEID), indices of the first (Repeats_First) and 
    # the last (Repeats_Last) repetitive elements in each gene id
    trfnumTable <- data.frame(GENEID=trfTable$Line[repetitive.geneid],  
                              Repeats_First=repetitive.first, 
                              Repeats_Last=repetitive.last) %>% mutate(Repeats_Count=Repeats_Last - Repeats_First + 1,
                                       GENEID=str_replace(GENEID, "Sequence: ", "")) 


    # Exploring the output data frame
    head(trfnumTable)

    ##            GENEID Repeats_First Repeats_Last Repeats_Count
    ## 1 ENSG00000162676             4            4             1
    ## 2 ENSG00000171791             9           10             2
    ## 3 ENSG00000069011            13           13             1

    dim(trfnumTable)

    ## [1] 3 4

    # Explore the output data frame
    glimpse(trfnumTable)

    ## Rows: 3
    ## Columns: 4
    ## $ GENEID        <chr> "ENSG00000162676", "ENSG00000171791", "ENSG00000069011"
    ## $ Repeats_First <int> 4, 9, 13
    ## $ Repeats_Last  <dbl> 4, 10, 13
    ## $ Repeats_Count <dbl> 1, 2, 1

    # Export as a csv file
    write.csv(trfnumTable, "TRF_demo.csv")

### Session Info

    sessionInfo()

    ## R version 4.0.3 (2020-10-10)
    ## Platform: x86_64-conda-linux-gnu (64-bit)
    ## Running under: Ubuntu 20.04.2 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /home/mira/miniconda3/envs/snakemake_r/lib/libopenblasp-r0.3.15.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] biomaRt_2.46.3    forcats_0.5.1     stringr_1.4.0     dplyr_1.0.5      
    ##  [5] purrr_0.3.4       readr_1.4.0       tidyr_1.1.3       tibble_3.1.1     
    ##  [9] ggplot2_3.3.3     tidyverse_1.3.1   data.table_1.14.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Biobase_2.50.0       httr_1.4.2           bit64_4.0.5         
    ##  [4] jsonlite_1.7.2       modelr_0.1.8         assertthat_0.2.1    
    ##  [7] askpass_1.1          BiocFileCache_1.14.0 stats4_4.0.3        
    ## [10] blob_1.2.1           cellranger_1.1.0     yaml_2.2.1          
    ## [13] progress_1.2.2       pillar_1.6.0         RSQLite_2.2.5       
    ## [16] backports_1.2.1      glue_1.4.2           digest_0.6.27       
    ## [19] rvest_1.0.0          colorspace_2.0-1     htmltools_0.5.1.1   
    ## [22] XML_3.99-0.6         pkgconfig_2.0.3      broom_0.7.6         
    ## [25] haven_2.4.1          scales_1.1.1         openssl_1.4.4       
    ## [28] generics_0.1.0       IRanges_2.24.1       ellipsis_0.3.2      
    ## [31] cachem_1.0.4         withr_2.4.2          BiocGenerics_0.36.0 
    ## [34] cli_2.5.0            magrittr_2.0.1       crayon_1.4.1        
    ## [37] readxl_1.3.1         memoise_2.0.0        evaluate_0.14       
    ## [40] ps_1.6.0             fs_1.5.0             fansi_0.4.2         
    ## [43] xml2_1.3.2           tools_4.0.3          prettyunits_1.1.1   
    ## [46] hms_1.0.0            lifecycle_1.0.0      S4Vectors_0.28.1    
    ## [49] munsell_0.5.0        reprex_2.0.0         AnnotationDbi_1.52.0
    ## [52] compiler_4.0.3       rlang_0.4.11         grid_4.0.3          
    ## [55] rstudioapi_0.13      rappdirs_0.3.3       rmarkdown_2.7       
    ## [58] gtable_0.3.0         curl_4.3.1           DBI_1.1.1           
    ## [61] R6_2.5.0             lubridate_1.7.10     knitr_1.31          
    ## [64] fastmap_1.1.0        bit_4.0.4            utf8_1.2.1          
    ## [67] stringi_1.5.3        parallel_4.0.3       Rcpp_1.0.6          
    ## [70] vctrs_0.3.8          dbplyr_2.1.1         tidyselect_1.1.1    
    ## [73] xfun_0.20
