## Utilization of Tandem Repeats Finder (TRF) by Mira Sohn

### Aims

#### - This workflow is aimed at finding and counting repetitive elements based on ensembl gene id using the Tandem Repeats Finder (Unix/Linux ver.)

#### - References: [TRF Web](https://tandem.bu.edu/trf/trf.html), [TRF linux](https://tandem.bu.edu/trf/trf.unix.help.html), [TRF paper](https://pubmed.ncbi.nlm.nih.gov/9862982/), [biomaRt Doc](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)

### Conda environment

#### - This analysis is used under [conda](https://conda.io/projects/conda/en/latest/index.html) environment (see the [recipe](https://github.com/Mira0507/TRF_demo/blob/master/conda_r.yml))

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

    ## [1] 216

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
    ## 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     TCCAAGTCCCAGCGAACCCGCGTGCAACCTGTCCCGACTCTAGCCGCCTCTTCAGCTCGCCATGGATCCCAACTGCTCCTGCGCCGCCGGTAAGAGGCTGGGGATGCCCAGTGTAGACTGTAGCGCTAGAGAAGCAATTTCTGACCCCTCTTTCTTTCTCTGGTCACTCAATTTCAGGACAGGAGTTGCTCCTTCCCAAAGAGTTTTGGGGTATCTTTCTCTCCATTCTAGGTTATTCGGAGCCCCCTTTTTACCGTTAAGGAGATCTGAGTTAATGGCTTGCTCAAGTTCCCAGGAATCGGTTGTGGACTGAGGAACTCGGCCCCGGGCTCTTAGTACGCCGTCCCTTGTTCAGGTATCCAGGGACGGTTCTCACCTCTGTCTTTTCTCCTTGCAGGTGACTCCTGCACCTGCGCCGGCTCCTGCAAATGCAAAGAGTGCAAATGCACCTCCTGCAAGAAAAGTAAGTGGGATCCTCTCTTTCCTCTACCCCTTCCCTGTCCTCCAGCCTGTCCCCTCTCCACCATCCTCAGGGGAATTAAAGCAGTCTGGGGATGCCCCATTGCGCGGAAATTGTTGCCTCCTCAGTGATCCTTATCAGGGAGAGCAGGAATCCTTATTCCCGGTGTCGCTAGTACTCATCTCTGCCGCCTCCTGTCTGCCCCCAGGCTGCTGCTCCTGCTGCCCTGTGGGCTGTGCCAAGTGTGCCCAGGGCTGCATCTGCAAAGGGGCGTCGGACAAGTGCAGCTGCTGCGCCTGATGCTGGGACAGCCCCGCTCCCAGATGTAAAGAACGCGACTTCCACAAACCTGGATTTTTTATGTACAACCCTGACCGTGACCGTTTGCTATATTCCTTTTTCTATGAAATAATGTGAATGATAATAAAACAGCTTTGACTTGA
    ## 9                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        CGAAGGGAACTGAATGAGGACAAGCTGGAGAAGCTGGAGGAGCTGACAATGGATGGGGCCAAGGCTAAGGCTATTCTGGATGCCTCACGGTCCTCCATGGGTCAGTGCAGAGCCTGGCAACCTGCATAAGGTATGGGGCTCTAAATAGCTGGCCTCTTGCATTCACACTTGGTTTTTCCTAGGCATGGACATATCTGCCATTGACTTGATAAACATCGAGAGCTTCTCCAGTCGTGTGGTGTCTTTATCTGAATACCGCCAGAGCCTACACACTTACCTGCGCTCCAAGATGAGCCAAGTAGCCCCCAGCCTGTCAGCCCTAATTGGGGAAGCGGTGCGTCACAGGGGACTCAAAAATGGGAGAATAAGGACTGTTGCCATGTGCACCTGCACTGCTGTATTTCGTGACCCACCATGTCTTCCCTAGTTGTGCTTGATGGGGAGGTGGGGAGCAGGGCTGTCGTGCAACTGGGCAGGTCAGCAGTTCATTTCTCTGACTGCTTCCTTGACTCTCTCTCCAGGTAGGTGCACGTCTCATCGCACATGCTGGCAGCCTCACCAACCTGGCCAAGTATCCAGCATCCACAGTGCAGATCCTTGGGGCTGAAAAGGCCCTGTTCAGGTACCAGTGAGGGCACCTGCCCACAATCAGGTGCCACTTCTGGTGCCCACTGCTTGTTGGGGGATCACGGTGATGGCTGACCAGGGCTCCCTGACCTATACAGGCCTCTGCTATGGGGGTGATGGCCAGTCCTGGTGTCTGAGTGATTCCCAGGGCCCAGCAAAGGGACCAAGTTTCCAGGTCAGCGACATTGGATGCCTTCCCTCTGCCTCTGGGAGCTATGGGTTGGCATGCATTGGGGTAGAGATCCAATCTGGCCTGAGGCTCACTCAGGACTTCGGGGTGAGAGGAGGGGAGGAGCTGAGCTGCCTTGGCTAATGGGGTTGAAATTTCTGATCTTAAACTCTCCACTGAATATTCTCTCAGAGCCCTGAAGACAAGGGGTAACACTCCAAAATATGGACTCATTTTCCACTCCACCTTCATTGGCCGAGCAGCTGCCAAGAACAAAGGCCGCATCTCCCGATACCTGGCAAACAAATGCAGTATTGCCTCACGAATCGATTGCTTCTCTGGTATGGGTGGGGGGGCGTTGGCAGGTGTGAGAAGGGGCTGGGTGGCTGGGTGGGGAGGCTTGCAACCATAGCTTCCACAATGATGGCAATATTTTTCGTCAACAGCAGTTCACCTAGTGAGTGTTGAGACTCTGGGTCTGAGTGAAGCTGAGGGTAGAGGGAACACAGGGTTGGGGTAGTTTCTCTCTTTGGGCTGACAGGCTTTGTCACCCACACACATCCAGAGGTGCCCACGAGTGTATTCGGGGAGAAGCTTCGAGAACAAGTTGAAGAGCGACTGTCCTTCTATGAGACTGGAGAGATACCACGAAAGAATCTGGATGTCATGAAGGAAGCAATGGTTCAGGTCAGTTGGGCTTTGCTGGGTGTGGAGTGGCATAGCTAGCTGTTGGAGGTGATGAACTGTCTGAGCCTGACCTTGTAGAATGGAGGCAAAAAAACTGATTTAATGAGCCTGATCCAATAAAGCCAGAAAGGAGTCCTCAGAGCACCAGAAGTCTTCAGGCCCTTTTAGCACTTTTCTTTGACCAGGCAGAGGAAGCGGCTGCTGAGATTACTAGGAAGCTGGAGAAACAGGAGAAGAAACGCTTAAAGAAGGAAAAGAAACGGCTGGCTGCACTTGCCCTCGCGTCTTCAGAAAACAGCAGTAGTACTCCAGAGGAGTGTGAGGTCAGTAGGCAGCACGGCCCTGGCAGAGATCCTAGGTTGTAGGATTTTCAACAGCAGAACAAAGGATATGCTGCATCAAGCTGTGGTCTTGAGTCCAGGCTTTTGGACTGAAACAAGGACCTGAAACATCTAAAACTACCTCTTGATTCTATAGGAAGGAGATAGGTGCTGAACTTGCTCAAGAGCCCAGAGAGCTGGTTGTAGCTCACACCCGTTCCCTGGGCATGTGTGTTCTGTCCTCGGCTGCCTCCCAGGAGTCCTCAACCTGGGGTAGTGTAAATTCCTGCTCTGCTTATTATCAGACGTGTGTCCGGAGGTGGTCGTGTTTCACAGTGGGGATGGGGGCAGGGAGGTCCCCAATGTGCTAAGCTACAATCATTCTCCCTGAGATTTTCATTTAGCACCCAGTTTCTTAAACAGTGTTTCAGGGCCCTGTCTGGAACTTGGCATGATGGTTCTGTTGCGACCAGCATGGTGGGTGTTTTTTAGGTTTTTTTTTTTAATGGGCTGAGGTAATTTCTCATGACATGTTTTCCTTCTAATTTGGGACAGCCTTTGGGGTGGATTTCTAAAGTTATACCCACACAATTAAACTATCCCAGAAACACTGGGCAATGTTAACGACACGCGTTCCCCTGCCTTGGCTACTTAATTGCTGAAGATGTAATGAGCACTGTTCTCACAGCCTGTTCCCCTGTCCTTCCCTTTAGGAGATGAGTGAAAAACCCAAAAAGAAGAAAAAGCAAAAGCCCCAGGAGGTTCCTCAGGAGAATGGAATGGAAGACCCATCTATCTCTTTCTCCAAACCCAAGAAAAAGAAATCTTTTTCCAAGGAGGAGTTGATGAGTAGCGATCTTGAAGAGACCGCTGGCAGCACCAGTATTCCCAAGAGGAAGAAGTCTACACCCAAGGAGGAAACAGTTAATGACCCTGAGGAGGCAGGCCACAGAAGTGGCTCCAAGAAAAAGAGGAAATTCTCCAAAGAGGAGCCGGTCAGCAGTGGGCCTGAAGAGGCGGTTGGCAAGAGCAGCTCCAAGAAGAAGAAAAAGTTCCATAAAGCATCCCAGGAAGATTAGAATGCAAATGGACATTCTCTGGGAGGTGGGGCATACCATAGCCCAAGGTGACATTTCCCACCCTGTGCCGTGTTCCCCAATAAAAACAAATTCACAAG
    ## 2  AGTCGAAGGCTGAGCGGACCCGGTCGAGAGCCGGCGAGCTCCGCGCAGGGAGGGTGCGCCCACCGGTCCCGCCGGGCGCCCGCGGGACGCGCCGCCAGGGCCCTCTCCGCCGGGGGCTCGGCGCTCGCCCACCTCTTCCAAATTTAACCATTACCTAAATCCGAAGGGAAATGAGCAAACCTCTCGGATTGGGTGTCAAGGTCTCCTCCGGGCTGGGGCTGAGCAAGCCCTCGGAGTGACCGTGGGTGACAGCGGCTCCAGGGACTCTTGGGGCGCAGTGGGGAAAGTGCCGGACCACCATGCCGCGCTCATTTCTCGTCAAAAGCAAGAAGGCTCACAGCTACCACCAGCCGCGCTCCCCAGGACCAGACTATTCCCTCCGTTTAGAGAATGTACCGGCGCCTAGCCGAGCAGACAGCACTTCAAATGCAGGCGGGGCGAAGGCGGAGCCCCGGGACCGTTTGTCCCCCGAATCGCAGCTGACCGAAGCCCCAGACAGAGCCTCCGCATCCCCAGACAGCTGCGAAGGCAGCGTCTGCGAACGGAGCTCGGAGTTTGAGGACTTCTGGAGGCCCCCGTCACCCTCCGCGTCTCCAGCCTCGGAGAAGTCAATGTGCCCATCGCTGGACGAAGCCCAGCCCTTCCCCCTGCCTTTCAAACCGTACTCATGGAGCGGCCTGGCGGGTTCTGACCTGCGGCACCTGGTGCAGAGCTACCGACCGTGTGGGGCCCTGGAGCGTGGCGCTGGCCTGGGCCTCTTCTGCGAACCCGCCCCGGAGCCTGGCCACCCGGCCGCGCTGTACGGCCCGAAGCGGGCTGCCGGCGGCGCGGGGGCCGGGGCGCCAGGGAGCTGCAGCGCAGGGGCCGGTGCCACCGCTGGCCCTGGCCTAGGGCTCTACGGCGACTTCGGGTCTGCGGCAGCCGGGCTGTATGAGAGGCCCACGGCAGCGGCGGGCTTGCTGTACCCCGAGCGTGGCCACGGGCTGCACGCAGACAAGGGCGCTGGCGTCAAGGTGGAGTCGGAGCTGCTGTGCACCCGCCTGCTGCTGGGCGGCGGCTCCTACAAGTGCATCAAGTGCAGCAAGGTGTTCTCCACGCCGCACGGGCTCGAGGTGCACGTGCGCAGGTCCCACAGCGGTACCAGACCCTTTGCCTGCGAGATGTGCGGCAAGACCTTCGGGCACGCGGTGAGCCTGGAGCAGCACAAAGCCGTGCACTCGCAGGAACGGAGCTTTGACTGTAAGATCTGTGGGAAGAGCTTCAAGAGGTCATCCACACTGTCCACACACCTGCTTATCCACTCAGACACTCGGCCCTACCCCTGTCAGTACTGTGGCAAGAGGTTCCACCAGAAGTCAGACATGAAGAAACACACTTTCATCCACACTGGTGAGAAGCCTCACAAGTGCCAGGTGTGCGGCAAGGCATTCAGCCAGAGCTCCAACCTCATCACCCACAGCCGCAAACACACAGGCTTCAAGCCCTTCGGCTGCGACCTCTGTGGGAAGGGTTTCCAGAGGAAGGTGGACCTCCGAAGGCACCGGGAGACGCAGCATGGGCTCAAATGAGCACCCTGGCTGGCTGCAAGCAGCAGCTACACAACACTACAGAGGGCAGCCTCCCTGCTTGCCACCACTCTGCTCCCTGCTTGCCTCCACTCCCTTCTGACTTTCCAGACCCCAGGTCCAGTCTGCAGATCCTACCAGGTTGCTCCTCCTTCGCCTTACCTCCTGGAGCTGCCAGAAGAAATGAGGTACCTTTTCAAAGTGCAGCCGAGAGTGAGAACCAAGTGACTCTCTAGGCTTCGGACACAAATAGGCTCCTCTACACCTGAAGACAAAGGCAAAGTCAAATGGGGACCAGAATAAATCTTAGACCCCACAGTCCTTCCCATTTCCAGCCCTAATCTACAGACAGGAATGCCCTTCAGGTTTCTTCCCTCCCCCCTCTTGACCTACCCCAGATATTTGTGTGGAAGAGGAGGAATCACCATTTACAAGGTGGACAAATGCTAATATTTTTATCTAGAAAGAAGAGTGAGTGTTAACTTTTATTTTTTTCCTTCTGGGGGGTCTGTTGACTCCTTTCTTTTGGGTGCTGCCTATAAATCTTGGAGGAATCATTTCTCCTCCTCAAAAACTGATTCAGAAACTGACTTGGGGAAGGAATTTAATACTTTGAAGTCATGAGATGCACCATCGAGGCTACCCCCAAGAAGAAGCAGAAGAGAAGTTGGTAATGAGAGGGGATTAGAGGTCCTCCCTTCAGTAGGGCTGTGAAAACCTCATCACTGGAGGTAAAAGCACAAGCAATGCCTGTGGACAAGATGTCATTCATTCACTCAGCAAATGTTCATGGATCACCGGCTACCAAGGTACCAGGCACCATGCTAGGTATTGGGGAAGAGAGACTGAAGTCACAACCCCTGACTGCTCCTCAAAAGCTAACGGTTGCACCTCCAAGTGGCTGGGTCTGTTCTTACTCTTGGAGGGAATTCTGAGAAGACAGCACAGAATTGTAAACCTTCCCTTTTGACCCTTTTGGATTTTATCAGGTGTAAACAAAAAGCTGAACAGTTACTTCAAAGATATGTGTGTATATTCAGTTTTTTATTGTTAAGCTGATATTTTAAAGATTTCTGAGCTAGCAGGCATGTGGGAAGGAAGGCTCTGTCTTCAACTCTTTGACCCTCCATGTGTACCATAGAGGGGGGAAAGGTGGTATTTTCACTTTGATGAGGTTGGTAAATGTTTTTAGATCTTCTGGTAAGCATTATGTTTGTTAATACATATTTATTAGAGTGATGTTTTAAGTTAATAAAGTATTAAGAGTATTACAGATTGCCTTTCCTTTAGGCAGGGAATTTGGGGAAGACCTTTGTTTTAGTCTGTAGGGGGAAAGCTGGCACTGAAGGAAACTGCAATGGAGTATAAAAGTGAAGCACACAAGCTAATATCCAGAAGAGTAAAATCCAATTGCTCGGGTTTCCATATACTTTCTCCTTTGCTGCTACACCTCAGAGATACTTCTTTCATCTTCTTTCATAAGAGAAAAACTAGGAGGTGTACAAAACTAAGCTACTACCTTGTAGAGCCCCTGATGGTTCCACAGGCACGTAAGTAAGCTCTTTCCTCTGAAAATGTCAACATATTTTTTTAGGGGCCTGAAGGCTCCCAGATGAATAAAGAAATCACAGGTCTTAGGTCTAGGTTTTAGTTTTAGTCACAGCTCTGTCACCAACACATCATGTGACTTGTGTAAATCTTTTGACTTCCCCGATACTGGTTTCCTTGATCAAGTAATACGGAGACAGCTCTGTCATTCCTGTCCAACCCTGACCGCTTTGATTCTGTGGAGACTGGAAAGGAAAAAAAGTGTGAGATTGGTTAGCTGCCTCCCCTAATGCAAGTCACAAATCTAAAGGCCCATGAGGATCTGGGTTATTTCTGGTTCTGTCCCACCTGTGGGAGTTTTTTGCATAAGTCAATCCTTTAAGAGAATGAAACCAGAATAAAATCCTGGCCAAAATCCCTACTACTAGTGCTTTAACACCTCCCAGCCAGACATTCAGGGTATGCATGCCCCATGGTGAGAACCAGCAGATTGCCCATTTCCAACTCAAATTGATTTCACACATATAATTGAGCACCCAGTATGTGCTGGGCATTCTTCTAGGGGGTGGAGATATAGCAGTGAGGACAATGTCTTACGAAGCTTTATTTATTCATTGTTGAATCCTCAGTTCCCAGAATAGTGCCTGGCACACAGTAGGAGATTAATAAATATTTACTGAGATGTTAAGTGAATAAACAGAAAATAAGCTCCTGCTTCATCCTAGGGCCTGGCTCTGCTCTCTACAGCACAAAGGGAAGAGAGATTGCTTGACAGCCTGGGAAACAAATAGAATCCACACGACTTTCTCTCCAGGTTGGCAGGAAAGGGGGGCTTGAATGCTGCTTTTCTCTGGTTTCTCCAAATCTCAGTCCTTGGCTGTCCTCACCTGATTTCCAGCACCATGTCTTTCCATGTCAGTTTACTCTCATTCAGAATTTAAATTCGTTTACTCTGTAATGGAGTTTAAGATAGTTTTTCAGAACTGATGGATGATGAACAGTTGCAGCGTCCAGGAACCACTTCCTAGTGTCCTGGGAGACAATCCTCCCAGTGAATGAAAATGACAAAAACATCCTTCCCACAGAAGTCACCCTGTGTACTGGTTGGTTATTTCCAACCTGGTTTGAAAAGTAAAGAAGTTGCACTGAGGGGTAAGTAGAACTGATTGAAATCTAGGTGCCTGTGCTGTCCTGCTCTTCCCTTATGACGTTTAGAACATTTGGTCCATTTTTGTTATCCCTGTTTATGCTCTAAGTGAGGAATTGTAATTTTGATGTTAACTGTATGCCTTTTTTTTTTTTTTTTGGTAAAAAATGTATCTTGATATCCCCAGTTCAGAATCAAAACAAATTACTGTGCCTGTTACCGTGTGAAATTTGAAGATTGCAAATATCTCTCTTAATGTTATTACTTGACTCCATTAAAAATAGATTGTGAATA
    ## 21                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 GGGAGCGTGCGCGCCTCTTGCCCGCCCGCGGGCCGCAAGATTCCGGAGGTGGTCCACCCCGTGGTCCACCTACCTCTGCTTCTTTCCCTTCGCCCAAACGGCAGCCCTTCCGCATTGCTTCCTGTCCTTTAGCGCGCGCACCCGTCACCTCACGCTGCACTTCTTTCGACCCCCTCCAGGCGACCCTGTATTTCCCTTTTTTCCCCCTTTACTCATCCTTCCCCTTCCTGCCACCATTCTGTCTCTCCCATTTACCCCATCGTGGCAGCACCCACTCCCCTTTCCCTCCCTGTCGCATCTTTTGTTCCCTGTGGCGCGCAGGTCACCCTCCGTTTTCTTCATGGGGACGCGGTGCTGGCGCGCAGTTTCCCGCAGAAATCCTGTGGAGGGGTCTGAGTGGCGCTACAGCTTCAGCTTGCATTCCTCTGTCAGCCCCGTCTCCCCCAGGAGTGGGAATAGCTTTGCGAAAAATGGGCGCAATGGCCAAGCCAGATTGTATCATCACTTGTGATGGTAAAAACCTCACCATAAAAACTGAGAGCACTTTGAAAACAACACAGTTTTCTTGTACCCTGGGAGAGAAGTTTGAAGAAACCACAGCTGATGGCAGAAAAACTCAGACTGTCTGCAACTTTACAGATGGTGCATTGGTTCAGCATCAGGAGTGGGATGGGAAGGAAAGCACAATAACAAGAAAATTGAAAGATGGGAAATTAGTGGTGGAGTGTGTCATGAACAATGTCACCTGTACTCGGATCTATGAAAAAGTAGAATAAAAATTCCATCATCACTTTGGACAGGAGTTAATTAAGAGAATGACCAAGCTCAGTTCAATGAGCAAATCTCCATACTGTTTCTTTCTTTTTTTTTTCATTACTGTGTTCAATTATCTTTATCATAAACATTTTACATGCAGCTATTTCAAAGTGTGTTGGATTAATTAGGATCATCCCTTTGGTTAATAAATAAATGTGTTTGTGCTAATA
    ## 1                                                                                                                                                                                                                                                                                                            CTGCTGGGGTGAGCAGCACTGTAAAGATGAAGCTGGCTAACTGGTACTGGCTGAGCTCAGCTGTTCTTGCCACTTACGGTTTTTTGGTTGTGGCAAACAATGAAACAGAGGAAATTAAAGATGAAAGAGCAAAGGATGTCTGCCCAGTGAGACTAGAAAGCAGAGGGAAATGCGAAGAGGCAGGGGAGTGCCCCTACCAGGTAAGCCTGCCCCCCTTGACTATTCAGCTCCCGAAGCAATTCAGCAGGATCGAGGAGGTGTTCAAAGAAGTCCAAAACCTCAAGGAAATCGTAAATAGTCTAAAGAAATCTTGCCAAGACTGCAAGCTGCAGGCTGATGACAACGGAGACCCAGGCAGAAACGGACTGTTGTTACCCAGTACAGGAGCCCCGGGAGAGGTTGGTGATAACAGAGTTAGAGAATTAGAGAGTGAGGTTAACAAGCTGTCCTCTGAGCTAAAGAATGCCAAAGAGGAGATCAATGTACTTCATGGTCGCCTGGAGAAGCTGAATCTTGTAAATATGAACAACATAGAAAATTATGTTGACAGCAAAGTGGCAAATCTAACATTTGTTGTCAATAGTTTGGATGGCAAATGTTCAAAGTGTCCCAGCCAAGAACAAATACAGTCACGTCCAGTTCAACATCTAATATATAAAGATTGCTCTGACTACTACGCAATAGGCAAAAGAAGCAGTGAGACCTACAGAGTTACACCTGATCCCAAAAATAGTAGCTTTGAAGTTTACTGTGACATGGAGACCATGGGGGGAGGCTGGACAGTGCTGCAGGCACGTCTCGATGGGAGCACCAACTTCACCAGAACATGGCAAGACTACAAAGCAGGCTTTGGAAACCTCAGAAGGGAATTTTGGCTGGGGAACGATAAAATTCATCTTCTGACCAAGAGTAAGGAAATGATTCTGAGAATAGATCTTGAAGACTTTAATGGTGTCGAACTATATGCCTTGTATGATCAGTTTTATGTGGCTAATGAGTTTCTCAAATATCGTTTACACGTTGGTAACTATAATGGCACAGCTGGAGATGCATTACGTTTCAACAAACATTACAACCACGATCTGAAGTTTTTCACCACTCCAGATAAAGACAATGATCGATATCCTTCTGGGAACTGTGGGCTGTACTACAGTTCAGGCTGGTGGTTTGATGCATGTCTTTCTGCAAACTTAAATGGCAAATATTATCACCAAAAATACAGAGGTGTCCGTAATGGGATTTTCTGGGGTACCTGGCCTGGTGTAAGTGAGGCACACCCTGGTGGCTACAAGTCCTCCTTCAAAGAGGCTAAGATGATGATCAGACCCAAGCACTTTAAGCCATAAATCACTCTGTTCATTCCTCCAGGTATTCGTTATCTAATAGGGCAATTAATTCCTTCAGCACTTTAGAATATGCCTTGTTTCATATTTTTCATAGCTAAAAAATGATGTCTGACGGCTAGGTTCTTATGCTACACAGCATTTGAAATAAAGCTGAAAAACAATGCATTTTAAAGGAGTCCTTTGTTGTTATGCTGTTATCCAATGAACACTTGCAAGCAATTAGCAATATTGAGAATTATACATTAGATTTACAATTCTTTTAATTTCTATTGAAACTTTTTCTATTGCTTGTATTACTTGCTGTATTTAAAAAATAATTGTTGGCTGGGTGTGGTAGCTCACGCCTGTAATCCCAGCACTTTGGAATGTCAAGGCAGGCAGATCACTTGAGGTCAGGAGTTTGAGACCAGCCTGGCCAAACATGTGAAACGCTGTCTCTATTAAAAATACAAAAATTAGCCGGGCATGGTGGTACATGCCTGTAATCCTAGCTACTTGGGAGGCTGAGGCAGGAGAATCGCTTGAACCTGAGAGGAAGAGGTTGCAGTGAGCCAAGACTGAGCCACTGCACTCCAGCATGGGTGACAGAGAAAACTCTGTCTCAAACAAAAAAATAATAAAATTTATTCAGTAGGCTGGATTCTACACAAAGTAATCTGTATTTGGGCCATGATTTAAGCACATCTGAAGGTATATCACTCTTTTCAGGCTATAATTATTTGGGTAATCTTCATTCTGAGACAAACTTAATCTATATCATTTACTTTGCAACAGAACAACCCTACAGCATTTTGGTTCCCAGACTAAGGGAACTAATATCTATATAATTAAACTTGTTCATTTATCATTCATGAAATATAAAATACTTGTCATTTAAACCGTTTAAAAATGTGGTAGCATAATGTCACCCCAAAAAGCATTCAGAAAGCAATGTAACTGTGAAGACCAGGGTTTAAAGGTAATTCATTTATAGTTTATAACTCCTTAGATGTTTGATGTTGAAAACTGCTTTAACATGAAAATTATCTTCCTCTGCTCTGTGTGAACAATAGCTTTTAATTTAAGATTGCTCACTACTGTACTAGACTACTGGTAGGTTTTTTTGGGGGGGGTGGGTAGGGATATGTGGGTAATGAAGCATTTACTTACAGGCTATCATACTCTGAGGCCAATTTTATCTCCAAAGCAATAATATCATTAAGTGATTCACTTCATAGAAGGCTAAGTTTCTCTAGGACAGATAGAAAACATGAATTTTGAAATATATAGAACAGTAGTTAAAATACTATATATTTCAACCCTGGCTGGTAGATTGCTTATTTTACTATCAGAAACTAAAAGATAGATTTTTACCCAAACAGAAGTATCTGTAATTTTTATAATTCATCAATTCTGGAATGCTATATATAATATTTAAAAGACTTTTTAAATGTGTTTAATTTCATCATCGTAAAAAGGGATCATCTCAGAGAGAACAGCAGTATTCTGCGTATTTTTAAAAATGCTCTAGAGTAACATTTGAAGTAATTCACTGTAGTGTATGCCAGTCCTAGAAATAATTTTTTTAATTTCTGGTGTCTGTTTCTAATACACTAACCAAGTTTTCAAAATATATTTACAAAGATGCATCTTTACCCATTATTTTAAAATGATTAAGGAGGATAGTTGCTTCAGGTAACAAGCTAATTTTTCAAATATTAGGCCCTTACAGAACTATTTAGTCAAAAAGTAAGATATTCCTTTAAAATATATAACCCAAAGCTTTCAGTTAAACATGATATATCACAAATACTATTAAAATGTTAAAGAGAAATGCAAATAGCATTAAATGATGACCAAAATGTAAAATATTGTAGATTTCAAAAGCTGTGTCTCTATTAGGTGGGATACCAAATGTAAATGATGTAACTGACGTTGTTTTTTACTTTTTACTTTTTAAAAAAGACTAAAAACGTTTTGATATTATACAATGTATTTGTTTCAGATAAGGTCATTGTCATTTAGTATATATAATTAATATATGTACAAGTTTAAGTAAATTCCTGTGAGTAAAAATGGACTTATCACAAAACATAGTTCTAAAGAAAGGTATATGCTCATATACACGGTGTCCATTAATTTAATGGGAACTAGGTATAACTTCAGGAGAATTTGGCAAATAATTCATTAATCCATGTAAATATTCAAAAGCTTGTTCTATCCACATTATTTCAAGGGATCACTTTATTTTTCATTATACTTTCACAGCACTTTTCTAGTAAATTCTGTAACACAGAAATTCCATTTTGGAATCATTTCATGTTACCAATAATTCTAGACTTTTATAACATTTAACATGTTGATGGAAATAGATTACATCTGCTAGAACCTTTTGCCTTAACTATTCACCAATATATGCTAATATTCATAAATATGGATTGACTGTTTACAAACATTAGAATCTTGTCTTGGTTCCATTTTGATGGCTAATATTTGTTATCTTAATTAAGACTATTTCTGAGGTCATGATTACTTGAAAATATTGACTAAAACTGGGTCCTTAGAAATTCCAGGTGGAGCTGATTTACCTATGACTGAGGGGAAAAAAAAATCAAATTTTACTGATAATAGTAATGCTCCAAATGAATTAATGACACATCTGTTCAATAAATAAAGAGCTTAAATATACAAAACATAAGAAATCTGGGCAACAAAACTTGTGGTCTTTACTTTTGAATAGCTACCCAAGAAAAGGTTTTAAAGGTAAAAGTTATGAGTAATGTCATCACAATAAGCTCTTGTTTAAAATTCTTTTCTTTTATGTATAATTAGGTTTATGTTTCATGTCTTTTTAAAACCTTATAAAAGATTTAATTATCACATCTATTCTTCAATGTGGAAATATTAAATATTGTTGGTTGTAAAATAATA
    ## 31                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       GGTTTATTTTCCAGATGCAATCAATGCCCCAGTCACCTGCTGTTATAACTTCACCAATAGGAAGATCTCAGTGCAGAGGCTCGCGAGCTATAGAAGAATCACCAGCAGCAAGTGTCCCAAAGAAGCTGTGATGTGAGTTCAGCACACCAACCTTCCCTGGCCTGAAGTTCTTCCTTGTGGAGCAAGGGACAAGCCTCATAAACCTAGAGTCAGAGAGTGCACTATTTAACTTAATGTACAAAGGTTCCCAATGGGAAAACTGAGGCACCAAGGGAAAAAGTGAACCCCAACATCACTCTCCACCTGGGTGCCTATTCAGAACACCCCAATTTCTTTAGCTTGAAGTCAGGATGGCTCCACCTGGACACCTATAGGAGCAGTTTGCCCTGGGTTCCCTCCTTCCACCTGCGTTCCTCCTCTAGCTCCCATGGCAGCCCTTTGGTGCAGAATGGGCTGCACTTCTAGACCAAAACTGCAAAGGAACTTCATCTAACTCTGTCCTCCCTCCCCACAGCTTCAAGACCATTGTGGCCAAGGAGATCTGTGCTGACCCCAAGCAGAAGTGGGTTCAGGATTCCATGGACCACCTGGACAAGCAAACCCAAACTCCGAAGACTTGAACACTCACTCCACAACCCAAGAATCTGCAGCTAACTTATTTTCCCCTAGCTTTCCCCAGACACCCTGTTTTATTTTATTATAATGAATTTTGTTTGTTGATGTGAAACATTATGCCTTAAGTAATGTTAATTCTTATTTAAGTTATTGATGTTTTAAGTTTATCTTTCATGGTACTAGTGTTTTTTAGATACAGAGACTTGGGGAAATTGCTTTTCCTCTTGAACCACAGTTCTACCCCTGGGATGTTTTGAGGGTCTTTGCAAGAATCATTAATACAAAGAATTTTTTTTAACATTCCAATGCATTGCTAAAATATTATTGTGGAAATGAATATTTTGTAACTATTACACCAAATAAATATATTTTTGTACAAAA
    ##    ensembl_gene_id seqsize
    ## 3  ENSG00000125148     903
    ## 9  ENSG00000101361    2980
    ## 2  ENSG00000162676    4554
    ## 21 ENSG00000164687     986
    ## 1  ENSG00000127951    4256
    ## 31 ENSG00000108691     996

    dim(seqTable)

    ## [1] 217   3

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

    ## [1] 10986

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

        if (trflines[i] != "") {

            df <- data.frame(Line=trflines[i])

            trfTable <- rbind(trfTable, df)

        }



    }


    # Exploring the output data frame
    head(trfTable)

    ##                             Line
    ## 1      Sequence: ENSG00000125148
    ## 2 Parameters: 2 7 7 80 10 50 500
    ## 3      Sequence: ENSG00000101361
    ## 4 Parameters: 2 7 7 80 10 50 500
    ## 5      Sequence: ENSG00000162676
    ## 6 Parameters: 2 7 7 80 10 50 500

    dim(trfTable)

    ## [1] 3420    1

    class(trfTable)

    ## [1] "data.frame"

    # Clean the TRF data frame 
    # (removing rows containing parameters and adding line codes)
    trfTable <- trfTable %>%
        mutate(code=ifelse(grepl("Sequence", Line, fixed=T), "S", 
                           ifelse(grepl("Parameters", Line, fixed=T), "P",
                           "CDNA"))) %>% dplyr::filter(code != "P")


    # Exploring the output data frame
    head(trfTable)

    ##                                                                                                            Line
    ## 1                                                                                     Sequence: ENSG00000125148
    ## 2                                                                                     Sequence: ENSG00000101361
    ## 3                                                                                     Sequence: ENSG00000162676
    ## 4 1619 1660 22 1.9 22 95 0 75 7 52 11 28 1.64 CTCCCTGCTTGCCACCACTCTG CTCCCTGCTTGCCACCACTCTGCTCCCTGCTTGCCTCCACTC
    ## 5                                                                                     Sequence: ENSG00000164687
    ## 6                                                                                     Sequence: ENSG00000127951
    ##   code
    ## 1    S
    ## 2    S
    ## 3    S
    ## 4 CDNA
    ## 5    S
    ## 6    S

    dim(trfTable)

    ## [1] 2340    2

    # Extracting the number of repetitive elements per GENEID
    #
    # Creating an empty data frame 
    trfnumTable <- data.frame()

    # Converting the line code to a character vector
    code.vector <- trfTable$code


    # Creating empty vectors 
    GENE <- c()    
    COUNT <- c()

    # Creating variables for storing counts 
    num.s <- 0
    num.c <- 0

    for (i in 1:nrow(trfTable)) {

        # If the line code is "S" (=Sequence: gene id)
        if (code.vector[i] == "S") {

            num.s <- num.s + 1
            num.c <- 0

        # If the line code is "CDNA" (=Repeats) 
        } else {

            num.s <- 0 
            num.c <- num.c + 1
        }

        GENE[i] <- num.s # Indicates gene id coded "S"
        COUNT[i] <- num.c  # Indicates repetitive elements coded "CDNA" 


    }


    # Exploring the GENE and Count vectors
    GENE

    ##    [1]  1  2  3  0  1  2  3  4  0  0  1  2  0  1  2  3  4  5  6  7  0  0  0  0
    ##   [25]  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  2  0  0
    ##   [49]  1  2  0  0  1  2  0  1  2  3  0  0  0  0  0  0  1  2  0  1  2  3  0  0
    ##   [73]  1  0  0  0  1  0  1  0  0  1  0  1  0  0  0  0  1  0  0  0  1  2  3  4
    ##   [97]  5  6  7  8  9 10 11 12  0  0  1  2  3  4  5  6  0  0  1  2  0  1  2  3
    ##  [121]  0  0  0  1  0  1  0  1  0  0  0  1  2  0  1  2  0  0  0  0  1  0  0  0
    ##  [145]  0  0  1  0  0  1  2  0  0  1  2  3  4  0  0  1  2  3  4  0  1  0  1  0
    ##  [169]  1  0  0  0  0  1  2  3  0  0  0  0  0  0  1  0  0  0  0  0  1  2  3  0
    ##  [193]  1  0  0  0  0  1  2  0  0  0  0  0  1  2  0  0  1  0  0  1  2  3  4  5
    ##  [217]  6  7  8  0  1  0  1  0  0  1  2  3  4  5  6  7  0  0  1  2  0  0  0  1
    ##  [241]  0  0  0  0  1  0  0  0  0  1  2  0  0  1  2  0  1  0  0  0  0  0  0  0
    ##  [265]  0  0  1  0  1  2  3  0  0  1  0  0  0  0  0  0  0  1  2  3  4  0  0  0
    ##  [289]  0  0  0  0  0  0  0  1  2  0  1  2  0  0  0  0  1  0  1  2  3  0  1  2
    ##  [313]  3  4  5  6  0  1  2  0  0  1  2  3  4  0  0  0  0  1  0  0  1  2  0  1
    ##  [337]  2  3  4  5  6  7  0  0  0  0  1  0  0  0  0  0  1  0  0  0  1  2  0  0
    ##  [361]  1  0  1  2  3  4  0  1  0  1  0  1  0  1  2  0  0  0  0  0  0  1  2  3
    ##  [385]  0  1  0  0  1  2  0  1  2  0  1  2  0  1  2  3  4  5  0  0  0  0  1  2
    ##  [409]  3  4  0  1  0  1  2  0  0  1  0  0  1  0  0  1  0  1  0  1  0  0  1  0
    ##  [433]  1  0  1  0  0  1  0  1  0  0  1  0  0  0  0  0  0  0  0  0  1  2  0  0
    ##  [457]  0  1  2  0  1  2  3  4  5  0  1  0  1  2  3  0  1  2  3  4  0  0  1  2
    ##  [481]  0  1  2  3  4  5  6  7  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
    ##  [505]  0  0  0  0  0  0  0  0  1  2  0  0  1  2  0  0  1  2  0  1  2  3  0  0
    ##  [529]  0  0  0  0  1  2  0  1  2  3  0  0  1  0  0  0  1  0  1  0  0  1  0  1
    ##  [553]  0  0  0  0  1  0  0  0  1  2  3  4  5  6  7  8  9 10 11 12  0  0  1  2
    ##  [577]  3  4  5  6  0  0  1  2  0  1  2  3  0  0  0  1  0  1  0  1  0  0  0  1
    ##  [601]  2  0  1  2  0  0  0  0  1  0  0  0  0  0  1  0  0  1  2  0  0  1  2  3
    ##  [625]  4  0  0  1  2  3  4  0  1  0  1  0  1  0  0  0  0  1  2  3  0  0  0  0
    ##  [649]  0  0  1  0  0  0  0  0  1  2  3  0  1  0  0  0  0  1  2  0  0  0  0  0
    ##  [673]  1  2  0  0  1  0  0  1  2  3  4  5  6  7  8  0  1  0  1  0  0  1  2  3
    ##  [697]  4  5  6  7  0  0  1  2  0  0  0  1  0  0  0  0  1  0  0  0  0  1  2  0
    ##  [721]  0  1  2  0  1  0  0  0  0  0  0  0  0  0  1  0  1  2  3  0  0  1  0  0
    ##  [745]  0  0  0  0  0  1  2  3  4  0  0  0  0  0  0  0  0  0  0  1  2  0  1  2
    ##  [769]  0  0  0  0  1  0  1  2  3  0  1  2  3  4  5  6  0  1  2  0  0  1  2  3
    ##  [793]  4  0  0  0  0  1  0  0  1  2  0  1  2  3  4  5  6  7  0  0  0  0  1  0
    ##  [817]  0  0  0  0  1  0  0  0  1  2  0  0  1  0  1  2  3  4  0  1  0  1  0  1
    ##  [841]  0  1  2  0  0  0  0  0  0  1  2  3  0  1  0  0  1  2  0  1  2  0  1  2
    ##  [865]  0  1  2  3  4  5  0  0  0  0  1  2  3  4  0  1  0  1  2  0  0  1  0  0
    ##  [889]  1  0  0  1  0  1  0  1  0  0  1  0  1  0  1  0  0  1  0  1  0  0  1  0
    ##  [913]  0  0  0  0  0  0  0  0  1  2  0  0  0  1  2  0  1  2  3  4  5  0  1  0
    ##  [937]  1  2  3  0  1  2  3  4  0  0  1  2  0  1  2  3  4  5  6  7  0  0  0  0
    ##  [961]  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  2  0  0
    ##  [985]  1  2  0  0  1  2  0  1  2  3  0  0  0  0  0  0  1  2  0  1  2  3  0  0
    ## [1009]  1  0  0  0  1  0  1  0  0  1  0  1  0  0  0  0  1  0  0  0  1  2  3  4
    ## [1033]  5  6  7  8  9 10 11 12  0  0  1  2  3  4  5  6  0  0  1  2  0  1  2  3
    ## [1057]  0  0  0  1  0  1  0  1  0  0  0  1  2  0  1  2  0  0  0  0  1  0  0  0
    ## [1081]  0  0  1  0  0  1  2  0  0  1  2  3  4  0  0  1  2  3  4  0  1  0  1  0
    ## [1105]  1  0  0  0  0  1  2  3  0  0  0  0  0  0  1  0  0  0  0  0  1  2  3  0
    ## [1129]  1  0  0  0  0  1  2  0  0  0  0  0  1  2  0  0  1  0  0  1  2  3  4  5
    ## [1153]  6  7  8  0  1  0  1  0  0  1  2  3  4  5  6  7  0  0  1  2  0  0  0  1
    ## [1177]  0  0  0  0  1  0  0  0  0  1  2  0  0  1  2  0  1  0  0  0  0  0  0  0
    ## [1201]  0  0  1  0  1  2  3  0  0  1  0  0  0  0  0  0  0  1  2  3  4  0  0  0
    ## [1225]  0  0  0  0  0  0  0  1  2  0  1  2  0  0  0  0  1  0  1  2  3  0  1  2
    ## [1249]  3  4  5  6  0  1  2  0  0  1  2  3  4  0  0  0  0  1  0  0  1  2  0  1
    ## [1273]  2  3  4  5  6  7  0  0  0  0  1  0  0  0  0  0  1  0  0  0  1  2  0  0
    ## [1297]  1  0  1  2  3  4  0  1  0  1  0  1  0  1  2  0  0  0  0  0  0  1  2  3
    ## [1321]  0  1  0  0  1  2  0  1  2  0  1  2  0  1  2  3  4  5  0  0  0  0  1  2
    ## [1345]  3  4  0  1  0  1  2  0  0  1  0  0  1  0  0  1  0  1  0  1  0  0  1  0
    ## [1369]  1  0  1  0  0  1  0  1  0  0  1  0  0  0  0  0  0  0  0  0  1  2  0  0
    ## [1393]  0  1  2  0  1  2  3  4  5  0  1  0  1  2  3  0  1  2  3  4  0  0  1  2
    ## [1417]  0  1  2  3  4  5  6  7  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
    ## [1441]  0  0  0  0  0  0  0  0  1  2  0  0  1  2  0  0  1  2  0  1  2  3  0  0
    ## [1465]  0  0  0  0  1  2  0  1  2  3  0  0  1  0  0  0  1  0  1  0  0  1  0  1
    ## [1489]  0  0  0  0  1  0  0  0  1  2  3  4  5  6  7  8  9 10 11 12  0  0  1  2
    ## [1513]  3  4  5  6  0  0  1  2  0  1  2  3  0  0  0  1  0  1  0  1  0  0  0  1
    ## [1537]  2  0  1  2  0  0  0  0  1  0  0  0  0  0  1  0  0  1  2  0  0  1  2  3
    ## [1561]  4  0  0  1  2  3  4  0  1  0  1  0  1  0  0  0  0  1  2  3  0  0  0  0
    ## [1585]  0  0  1  0  0  0  0  0  1  2  3  0  1  0  0  0  0  1  2  0  0  0  0  0
    ## [1609]  1  2  0  0  1  0  0  1  2  3  4  5  6  7  8  0  1  0  1  0  0  1  2  3
    ## [1633]  4  5  6  7  0  0  1  2  0  0  0  1  0  0  0  0  1  0  0  0  0  1  2  0
    ## [1657]  0  1  2  0  1  0  0  0  0  0  0  0  0  0  1  0  1  2  3  0  0  1  0  0
    ## [1681]  0  0  0  0  0  1  2  3  4  0  0  0  0  0  0  0  0  0  0  1  2  0  1  2
    ## [1705]  0  0  0  0  1  0  1  2  3  0  1  2  3  4  5  6  0  1  2  0  0  1  2  3
    ## [1729]  4  0  0  0  0  1  0  0  1  2  0  1  2  3  4  5  6  7  0  0  0  0  1  0
    ## [1753]  0  0  0  0  1  0  0  0  1  2  0  0  1  0  1  2  3  4  0  1  0  1  0  1
    ## [1777]  0  1  2  0  0  0  0  0  0  1  2  3  0  1  0  0  1  2  0  1  2  0  1  2
    ## [1801]  0  1  2  3  4  5  0  0  0  0  1  2  3  4  0  1  0  1  2  0  0  1  0  0
    ## [1825]  1  0  0  1  0  1  0  1  0  0  1  0  1  0  1  0  0  1  0  1  0  0  1  0
    ## [1849]  0  0  0  0  0  0  0  0  1  2  0  0  0  1  2  0  1  2  3  4  5  0  1  0
    ## [1873]  1  2  3  0  1  2  3  4  0  0  1  2  0  1  2  3  4  5  6  7  0  0  0  0
    ## [1897]  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  2  0  0
    ## [1921]  1  2  0  0  1  2  0  1  2  3  0  0  0  0  0  0  1  2  0  1  2  3  0  0
    ## [1945]  1  0  0  0  1  0  1  0  0  1  0  1  0  0  0  0  1  0  0  0  1  2  3  4
    ## [1969]  5  6  7  8  9 10 11 12  0  0  1  2  3  4  5  6  0  0  1  2  0  1  2  3
    ## [1993]  0  0  0  1  0  1  0  1  0  0  0  1  2  0  1  2  0  0  0  0  1  0  0  0
    ## [2017]  0  0  1  0  0  1  2  0  0  1  2  3  4  0  0  1  2  3  4  0  1  0  1  0
    ## [2041]  1  0  0  0  0  1  2  3  0  0  0  0  0  0  1  0  0  0  0  0  1  2  3  0
    ## [2065]  1  0  0  0  0  1  2  0  0  0  0  0  1  2  0  0  1  0  0  1  2  3  4  5
    ## [2089]  6  7  8  0  1  0  1  0  0  1  2  3  4  5  6  7  0  0  1  2  0  0  0  1
    ## [2113]  0  0  0  0  1  0  0  0  0  1  2  0  0  1  2  0  1  0  0  0  0  0  0  0
    ## [2137]  0  0  1  0  1  2  3  0  0  1  0  0  0  0  0  0  0  1  2  3  4  0  0  0
    ## [2161]  0  0  0  0  0  0  0  1  2  0  1  2  0  0  0  0  1  0  1  2  3  0  1  2
    ## [2185]  3  4  5  6  0  1  2  0  0  1  2  3  4  0  0  0  0  1  0  0  1  2  0  1
    ## [2209]  2  3  4  5  6  7  0  0  0  0  1  0  0  0  0  0  1  0  0  0  1  2  0  0
    ## [2233]  1  0  1  2  3  4  0  1  0  1  0  1  0  1  2  0  0  0  0  0  0  1  2  3
    ## [2257]  0  1  0  0  1  2  0  1  2  0  1  2  0  1  2  3  4  5  0  0  0  0  1  2
    ## [2281]  3  4  0  1  0  1  2  0  0  1  0  0  1  0  0  1  0  1  0  1  0  0  1  0
    ## [2305]  1  0  1  0  0  1  0  1  0  0  1  0  0  0  0  0  0  0  0  0  1  2  0  0
    ## [2329]  0  1  2  0  1  2  3  4  5  0  1  0

    COUNT

    ##    [1]  0  0  0  1  0  0  0  0  1  2  0  0  1  0  0  0  0  0  0  0  1  2  3  4
    ##   [25]  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24  0  0  1  2
    ##   [49]  0  0  1  2  0  0  1  0  0  0  1  2  3  4  5  6  0  0  1  0  0  0  1  2
    ##   [73]  0  1  2  3  0  1  0  1  2  0  1  0  1  2  3  4  0  1  2  3  0  0  0  0
    ##   [97]  0  0  0  0  0  0  0  0  1  2  0  0  0  0  0  0  1  2  0  0  1  0  0  0
    ##  [121]  1  2  3  0  1  0  1  0  1  2  3  0  0  1  0  0  1  2  3  4  0  1  2  3
    ##  [145]  4  5  0  1  2  0  0  1  2  0  0  0  0  1  2  0  0  0  0  1  0  1  0  1
    ##  [169]  0  1  2  3  4  0  0  0  1  2  3  4  5  6  0  1  2  3  4  5  0  0  0  1
    ##  [193]  0  1  2  3  4  0  0  1  2  3  4  5  0  0  1  2  0  1  2  0  0  0  0  0
    ##  [217]  0  0  0  1  0  1  0  1  2  0  0  0  0  0  0  0  1  2  0  0  1  2  3  0
    ##  [241]  1  2  3  4  0  1  2  3  4  0  0  1  2  0  0  1  0  1  2  3  4  5  6  7
    ##  [265]  8  9  0  1  0  0  0  1  2  0  1  2  3  4  5  6  7  0  0  0  0  1  2  3
    ##  [289]  4  5  6  7  8  9 10  0  0  1  0  0  1  2  3  4  0  1  0  0  0  1  0  0
    ##  [313]  0  0  0  0  1  0  0  1  2  0  0  0  0  1  2  3  4  0  1  2  0  0  1  0
    ##  [337]  0  0  0  0  0  0  1  2  3  4  0  1  2  3  4  5  0  1  2  3  0  0  1  2
    ##  [361]  0  1  0  0  0  0  1  0  1  0  1  0  1  0  0  1  2  3  4  5  6  0  0  0
    ##  [385]  1  0  1  2  0  0  1  0  0  1  0  0  1  0  0  0  0  0  1  2  3  4  0  0
    ##  [409]  0  0  1  0  1  0  0  1  2  0  1  2  0  1  2  0  1  0  1  0  1  2  0  1
    ##  [433]  0  1  0  1  2  0  1  0  1  2  0  1  2  3  4  5  6  7  8  9  0  0  1  2
    ##  [457]  3  0  0  1  0  0  0  0  0  1  0  1  0  0  0  1  0  0  0  0  1  2  0  0
    ##  [481]  1  0  0  0  0  0  0  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
    ##  [505] 17 18 19 20 21 22 23 24  0  0  1  2  0  0  1  2  0  0  1  0  0  0  1  2
    ##  [529]  3  4  5  6  0  0  1  0  0  0  1  2  0  1  2  3  0  1  0  1  2  0  1  0
    ##  [553]  1  2  3  4  0  1  2  3  0  0  0  0  0  0  0  0  0  0  0  0  1  2  0  0
    ##  [577]  0  0  0  0  1  2  0  0  1  0  0  0  1  2  3  0  1  0  1  0  1  2  3  0
    ##  [601]  0  1  0  0  1  2  3  4  0  1  2  3  4  5  0  1  2  0  0  1  2  0  0  0
    ##  [625]  0  1  2  0  0  0  0  1  0  1  0  1  0  1  2  3  4  0  0  0  1  2  3  4
    ##  [649]  5  6  0  1  2  3  4  5  0  0  0  1  0  1  2  3  4  0  0  1  2  3  4  5
    ##  [673]  0  0  1  2  0  1  2  0  0  0  0  0  0  0  0  1  0  1  0  1  2  0  0  0
    ##  [697]  0  0  0  0  1  2  0  0  1  2  3  0  1  2  3  4  0  1  2  3  4  0  0  1
    ##  [721]  2  0  0  1  0  1  2  3  4  5  6  7  8  9  0  1  0  0  0  1  2  0  1  2
    ##  [745]  3  4  5  6  7  0  0  0  0  1  2  3  4  5  6  7  8  9 10  0  0  1  0  0
    ##  [769]  1  2  3  4  0  1  0  0  0  1  0  0  0  0  0  0  1  0  0  1  2  0  0  0
    ##  [793]  0  1  2  3  4  0  1  2  0  0  1  0  0  0  0  0  0  0  1  2  3  4  0  1
    ##  [817]  2  3  4  5  0  1  2  3  0  0  1  2  0  1  0  0  0  0  1  0  1  0  1  0
    ##  [841]  1  0  0  1  2  3  4  5  6  0  0  0  1  0  1  2  0  0  1  0  0  1  0  0
    ##  [865]  1  0  0  0  0  0  1  2  3  4  0  0  0  0  1  0  1  0  0  1  2  0  1  2
    ##  [889]  0  1  2  0  1  0  1  0  1  2  0  1  0  1  0  1  2  0  1  0  1  2  0  1
    ##  [913]  2  3  4  5  6  7  8  9  0  0  1  2  3  0  0  1  0  0  0  0  0  1  0  1
    ##  [937]  0  0  0  1  0  0  0  0  1  2  0  0  1  0  0  0  0  0  0  0  1  2  3  4
    ##  [961]  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24  0  0  1  2
    ##  [985]  0  0  1  2  0  0  1  0  0  0  1  2  3  4  5  6  0  0  1  0  0  0  1  2
    ## [1009]  0  1  2  3  0  1  0  1  2  0  1  0  1  2  3  4  0  1  2  3  0  0  0  0
    ## [1033]  0  0  0  0  0  0  0  0  1  2  0  0  0  0  0  0  1  2  0  0  1  0  0  0
    ## [1057]  1  2  3  0  1  0  1  0  1  2  3  0  0  1  0  0  1  2  3  4  0  1  2  3
    ## [1081]  4  5  0  1  2  0  0  1  2  0  0  0  0  1  2  0  0  0  0  1  0  1  0  1
    ## [1105]  0  1  2  3  4  0  0  0  1  2  3  4  5  6  0  1  2  3  4  5  0  0  0  1
    ## [1129]  0  1  2  3  4  0  0  1  2  3  4  5  0  0  1  2  0  1  2  0  0  0  0  0
    ## [1153]  0  0  0  1  0  1  0  1  2  0  0  0  0  0  0  0  1  2  0  0  1  2  3  0
    ## [1177]  1  2  3  4  0  1  2  3  4  0  0  1  2  0  0  1  0  1  2  3  4  5  6  7
    ## [1201]  8  9  0  1  0  0  0  1  2  0  1  2  3  4  5  6  7  0  0  0  0  1  2  3
    ## [1225]  4  5  6  7  8  9 10  0  0  1  0  0  1  2  3  4  0  1  0  0  0  1  0  0
    ## [1249]  0  0  0  0  1  0  0  1  2  0  0  0  0  1  2  3  4  0  1  2  0  0  1  0
    ## [1273]  0  0  0  0  0  0  1  2  3  4  0  1  2  3  4  5  0  1  2  3  0  0  1  2
    ## [1297]  0  1  0  0  0  0  1  0  1  0  1  0  1  0  0  1  2  3  4  5  6  0  0  0
    ## [1321]  1  0  1  2  0  0  1  0  0  1  0  0  1  0  0  0  0  0  1  2  3  4  0  0
    ## [1345]  0  0  1  0  1  0  0  1  2  0  1  2  0  1  2  0  1  0  1  0  1  2  0  1
    ## [1369]  0  1  0  1  2  0  1  0  1  2  0  1  2  3  4  5  6  7  8  9  0  0  1  2
    ## [1393]  3  0  0  1  0  0  0  0  0  1  0  1  0  0  0  1  0  0  0  0  1  2  0  0
    ## [1417]  1  0  0  0  0  0  0  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
    ## [1441] 17 18 19 20 21 22 23 24  0  0  1  2  0  0  1  2  0  0  1  0  0  0  1  2
    ## [1465]  3  4  5  6  0  0  1  0  0  0  1  2  0  1  2  3  0  1  0  1  2  0  1  0
    ## [1489]  1  2  3  4  0  1  2  3  0  0  0  0  0  0  0  0  0  0  0  0  1  2  0  0
    ## [1513]  0  0  0  0  1  2  0  0  1  0  0  0  1  2  3  0  1  0  1  0  1  2  3  0
    ## [1537]  0  1  0  0  1  2  3  4  0  1  2  3  4  5  0  1  2  0  0  1  2  0  0  0
    ## [1561]  0  1  2  0  0  0  0  1  0  1  0  1  0  1  2  3  4  0  0  0  1  2  3  4
    ## [1585]  5  6  0  1  2  3  4  5  0  0  0  1  0  1  2  3  4  0  0  1  2  3  4  5
    ## [1609]  0  0  1  2  0  1  2  0  0  0  0  0  0  0  0  1  0  1  0  1  2  0  0  0
    ## [1633]  0  0  0  0  1  2  0  0  1  2  3  0  1  2  3  4  0  1  2  3  4  0  0  1
    ## [1657]  2  0  0  1  0  1  2  3  4  5  6  7  8  9  0  1  0  0  0  1  2  0  1  2
    ## [1681]  3  4  5  6  7  0  0  0  0  1  2  3  4  5  6  7  8  9 10  0  0  1  0  0
    ## [1705]  1  2  3  4  0  1  0  0  0  1  0  0  0  0  0  0  1  0  0  1  2  0  0  0
    ## [1729]  0  1  2  3  4  0  1  2  0  0  1  0  0  0  0  0  0  0  1  2  3  4  0  1
    ## [1753]  2  3  4  5  0  1  2  3  0  0  1  2  0  1  0  0  0  0  1  0  1  0  1  0
    ## [1777]  1  0  0  1  2  3  4  5  6  0  0  0  1  0  1  2  0  0  1  0  0  1  0  0
    ## [1801]  1  0  0  0  0  0  1  2  3  4  0  0  0  0  1  0  1  0  0  1  2  0  1  2
    ## [1825]  0  1  2  0  1  0  1  0  1  2  0  1  0  1  0  1  2  0  1  0  1  2  0  1
    ## [1849]  2  3  4  5  6  7  8  9  0  0  1  2  3  0  0  1  0  0  0  0  0  1  0  1
    ## [1873]  0  0  0  1  0  0  0  0  1  2  0  0  1  0  0  0  0  0  0  0  1  2  3  4
    ## [1897]  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24  0  0  1  2
    ## [1921]  0  0  1  2  0  0  1  0  0  0  1  2  3  4  5  6  0  0  1  0  0  0  1  2
    ## [1945]  0  1  2  3  0  1  0  1  2  0  1  0  1  2  3  4  0  1  2  3  0  0  0  0
    ## [1969]  0  0  0  0  0  0  0  0  1  2  0  0  0  0  0  0  1  2  0  0  1  0  0  0
    ## [1993]  1  2  3  0  1  0  1  0  1  2  3  0  0  1  0  0  1  2  3  4  0  1  2  3
    ## [2017]  4  5  0  1  2  0  0  1  2  0  0  0  0  1  2  0  0  0  0  1  0  1  0  1
    ## [2041]  0  1  2  3  4  0  0  0  1  2  3  4  5  6  0  1  2  3  4  5  0  0  0  1
    ## [2065]  0  1  2  3  4  0  0  1  2  3  4  5  0  0  1  2  0  1  2  0  0  0  0  0
    ## [2089]  0  0  0  1  0  1  0  1  2  0  0  0  0  0  0  0  1  2  0  0  1  2  3  0
    ## [2113]  1  2  3  4  0  1  2  3  4  0  0  1  2  0  0  1  0  1  2  3  4  5  6  7
    ## [2137]  8  9  0  1  0  0  0  1  2  0  1  2  3  4  5  6  7  0  0  0  0  1  2  3
    ## [2161]  4  5  6  7  8  9 10  0  0  1  0  0  1  2  3  4  0  1  0  0  0  1  0  0
    ## [2185]  0  0  0  0  1  0  0  1  2  0  0  0  0  1  2  3  4  0  1  2  0  0  1  0
    ## [2209]  0  0  0  0  0  0  1  2  3  4  0  1  2  3  4  5  0  1  2  3  0  0  1  2
    ## [2233]  0  1  0  0  0  0  1  0  1  0  1  0  1  0  0  1  2  3  4  5  6  0  0  0
    ## [2257]  1  0  1  2  0  0  1  0  0  1  0  0  1  0  0  0  0  0  1  2  3  4  0  0
    ## [2281]  0  0  1  0  1  0  0  1  2  0  1  2  0  1  2  0  1  0  1  0  1  2  0  1
    ## [2305]  0  1  0  1  2  0  1  0  1  2  0  1  2  3  4  5  6  7  8  9  0  0  1  2
    ## [2329]  3  0  0  1  0  0  0  0  0  1  0  1

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

    ##   [1]    3    8   12   20   46   50   54   58   66   70   73   77   79   82   84
    ##  [16]   89  104  112  116  120  124  126  128  133  136  141  147  151  157  163
    ##  [31]  165  167  169  176  183  191  193  199  206  209  219  221  223  232  236
    ##  [46]  240  245  251  255  257  267  271  274  285  297  300  305  309  316  319
    ##  [61]  325  330  334  342  347  353  358  361  366  368  370  372  375  384  386
    ##  [76]  390  393  396  402  410  412  415  418  421  424  426  428  431  433  435
    ##  [91]  438  440  443  454  459  465  467  471  476  480  488  514  518  522  526
    ## [106]  534  538  541  545  547  550  552  557  572  580  584  588  592  594  596
    ## [121]  601  604  609  615  619  625  631  633  635  637  644  651  659  661  667
    ## [136]  674  677  687  689  691  700  704  708  713  719  723  725  735  739  742
    ## [151]  753  765  768  773  777  784  787  793  798  802  810  815  821  826  829
    ## [166]  834  836  838  840  843  852  854  858  861  864  870  878  880  883  886
    ## [181]  889  892  894  896  899  901  903  906  908  911  922  927  933  935  939
    ## [196]  944  948  956  982  986  990  994 1002 1006 1009 1013 1015 1018 1020 1025
    ## [211] 1040 1048 1052 1056 1060 1062 1064 1069 1072 1077 1083 1087 1093 1099 1101
    ## [226] 1103 1105 1112 1119 1127 1129 1135 1142 1145 1155 1157 1159 1168 1172 1176
    ## [241] 1181 1187 1191 1193 1203 1207 1210 1221 1233 1236 1241 1245 1252 1255 1261
    ## [256] 1266 1270 1278 1283 1289 1294 1297 1302 1304 1306 1308 1311 1320 1322 1326
    ## [271] 1329 1332 1338 1346 1348 1351 1354 1357 1360 1362 1364 1367 1369 1371 1374
    ## [286] 1376 1379 1390 1395 1401 1403 1407 1412 1416 1424 1450 1454 1458 1462 1470
    ## [301] 1474 1477 1481 1483 1486 1488 1493 1508 1516 1520 1524 1528 1530 1532 1537
    ## [316] 1540 1545 1551 1555 1561 1567 1569 1571 1573 1580 1587 1595 1597 1603 1610
    ## [331] 1613 1623 1625 1627 1636 1640 1644 1649 1655 1659 1661 1671 1675 1678 1689
    ## [346] 1701 1704 1709 1713 1720 1723 1729 1734 1738 1746 1751 1757 1762 1765 1770
    ## [361] 1772 1774 1776 1779 1788 1790 1794 1797 1800 1806 1814 1816 1819 1822 1825
    ## [376] 1828 1830 1832 1835 1837 1839 1842 1844 1847 1858 1863 1869 1871 1875 1880
    ## [391] 1884 1892 1918 1922 1926 1930 1938 1942 1945 1949 1951 1954 1956 1961 1976
    ## [406] 1984 1988 1992 1996 1998 2000 2005 2008 2013 2019 2023 2029 2035 2037 2039
    ## [421] 2041 2048 2055 2063 2065 2071 2078 2081 2091 2093 2095 2104 2108 2112 2117
    ## [436] 2123 2127 2129 2139 2143 2146 2157 2169 2172 2177 2181 2188 2191 2197 2202
    ## [451] 2206 2214 2219 2225 2230 2233 2238 2240 2242 2244 2247 2256 2258 2262 2265
    ## [466] 2268 2274 2282 2284 2287 2290 2293 2296 2298 2300 2303 2305 2307 2310 2312
    ## [481] 2315 2326 2331 2337 2339

    repetitive.first

    ##   [1]    4    9   13   21   47   51   55   59   67   71   74   78   80   83   85
    ##  [16]   90  105  113  117  121  125  127  129  134  137  142  148  152  158  164
    ##  [31]  166  168  170  177  184  192  194  200  207  210  220  222  224  233  237
    ##  [46]  241  246  252  256  258  268  272  275  286  298  301  306  310  317  320
    ##  [61]  326  331  335  343  348  354  359  362  367  369  371  373  376  385  387
    ##  [76]  391  394  397  403  411  413  416  419  422  425  427  429  432  434  436
    ##  [91]  439  441  444  455  460  466  468  472  477  481  489  515  519  523  527
    ## [106]  535  539  542  546  548  551  553  558  573  581  585  589  593  595  597
    ## [121]  602  605  610  616  620  626  632  634  636  638  645  652  660  662  668
    ## [136]  675  678  688  690  692  701  705  709  714  720  724  726  736  740  743
    ## [151]  754  766  769  774  778  785  788  794  799  803  811  816  822  827  830
    ## [166]  835  837  839  841  844  853  855  859  862  865  871  879  881  884  887
    ## [181]  890  893  895  897  900  902  904  907  909  912  923  928  934  936  940
    ## [196]  945  949  957  983  987  991  995 1003 1007 1010 1014 1016 1019 1021 1026
    ## [211] 1041 1049 1053 1057 1061 1063 1065 1070 1073 1078 1084 1088 1094 1100 1102
    ## [226] 1104 1106 1113 1120 1128 1130 1136 1143 1146 1156 1158 1160 1169 1173 1177
    ## [241] 1182 1188 1192 1194 1204 1208 1211 1222 1234 1237 1242 1246 1253 1256 1262
    ## [256] 1267 1271 1279 1284 1290 1295 1298 1303 1305 1307 1309 1312 1321 1323 1327
    ## [271] 1330 1333 1339 1347 1349 1352 1355 1358 1361 1363 1365 1368 1370 1372 1375
    ## [286] 1377 1380 1391 1396 1402 1404 1408 1413 1417 1425 1451 1455 1459 1463 1471
    ## [301] 1475 1478 1482 1484 1487 1489 1494 1509 1517 1521 1525 1529 1531 1533 1538
    ## [316] 1541 1546 1552 1556 1562 1568 1570 1572 1574 1581 1588 1596 1598 1604 1611
    ## [331] 1614 1624 1626 1628 1637 1641 1645 1650 1656 1660 1662 1672 1676 1679 1690
    ## [346] 1702 1705 1710 1714 1721 1724 1730 1735 1739 1747 1752 1758 1763 1766 1771
    ## [361] 1773 1775 1777 1780 1789 1791 1795 1798 1801 1807 1815 1817 1820 1823 1826
    ## [376] 1829 1831 1833 1836 1838 1840 1843 1845 1848 1859 1864 1870 1872 1876 1881
    ## [391] 1885 1893 1919 1923 1927 1931 1939 1943 1946 1950 1952 1955 1957 1962 1977
    ## [406] 1985 1989 1993 1997 1999 2001 2006 2009 2014 2020 2024 2030 2036 2038 2040
    ## [421] 2042 2049 2056 2064 2066 2072 2079 2082 2092 2094 2096 2105 2109 2113 2118
    ## [436] 2124 2128 2130 2140 2144 2147 2158 2170 2173 2178 2182 2189 2192 2198 2203
    ## [451] 2207 2215 2220 2226 2231 2234 2239 2241 2243 2245 2248 2257 2259 2263 2266
    ## [466] 2269 2275 2283 2285 2288 2291 2294 2297 2299 2301 2304 2306 2308 2311 2313
    ## [481] 2316 2327 2332 2338 2340

    repetitive.last

    ##   [1]    0    4   10   13   44   48   52   55   64   67   72   76   78   81   83
    ##  [16]   88   92  106  114  117  123  125  127  131  134  140  146  149  153  159
    ##  [31]  164  166  168  173  182  188  192  197  204  208  211  220  222  225  234
    ##  [46]  239  244  249  253  256  266  268  273  281  295  298  304  306  310  317
    ##  [61]  321  329  332  335  346  352  356  360  362  367  369  371  373  381  385
    ##  [76]  388  391  394  397  406  411  413  417  420  423  425  427  430  432  434
    ##  [91]  437  439  442  452  457  460  466  468  472  478  481  512  516  520  523
    ## [106]  532  535  540  544  546  549  551  556  560  574  582  585  591  593  595
    ## [121]  599  602  608  614  617  621  627  632  634  636  641  650  656  660  665
    ## [136]  672  676  679  688  690  693  702  707  712  717  721  724  734  736  741
    ## [151]  749  763  766  772  774  778  785  789  797  800  803  814  820  824  828
    ## [166]  830  835  837  839  841  849  853  856  859  862  865  874  879  881  885
    ## [181]  888  891  893  895  898  900  902  905  907  910  920  925  928  934  936
    ## [196]  940  946  949  980  984  988  991 1000 1003 1008 1012 1014 1017 1019 1024
    ## [211] 1028 1042 1050 1053 1059 1061 1063 1067 1070 1076 1082 1085 1089 1095 1100
    ## [226] 1102 1104 1109 1118 1124 1128 1133 1140 1144 1147 1156 1158 1161 1170 1175
    ## [241] 1180 1185 1189 1192 1202 1204 1209 1217 1231 1234 1240 1242 1246 1253 1257
    ## [256] 1265 1268 1271 1282 1288 1292 1296 1298 1303 1305 1307 1309 1317 1321 1324
    ## [271] 1327 1330 1333 1342 1347 1349 1353 1356 1359 1361 1363 1366 1368 1370 1373
    ## [286] 1375 1378 1388 1393 1396 1402 1404 1408 1414 1417 1448 1452 1456 1459 1468
    ## [301] 1471 1476 1480 1482 1485 1487 1492 1496 1510 1518 1521 1527 1529 1531 1535
    ## [316] 1538 1544 1550 1553 1557 1563 1568 1570 1572 1577 1586 1592 1596 1601 1608
    ## [331] 1612 1615 1624 1626 1629 1638 1643 1648 1653 1657 1660 1670 1672 1677 1685
    ## [346] 1699 1702 1708 1710 1714 1721 1725 1733 1736 1739 1750 1756 1760 1764 1766
    ## [361] 1771 1773 1775 1777 1785 1789 1792 1795 1798 1801 1810 1815 1817 1821 1824
    ## [376] 1827 1829 1831 1834 1836 1838 1841 1843 1846 1856 1861 1864 1870 1872 1876
    ## [391] 1882 1885 1916 1920 1924 1927 1936 1939 1944 1948 1950 1953 1955 1960 1964
    ## [406] 1978 1986 1989 1995 1997 1999 2003 2006 2012 2018 2021 2025 2031 2036 2038
    ## [421] 2040 2045 2054 2060 2064 2069 2076 2080 2083 2092 2094 2097 2106 2111 2116
    ## [436] 2121 2125 2128 2138 2140 2145 2153 2167 2170 2176 2178 2182 2189 2193 2201
    ## [451] 2204 2207 2218 2224 2228 2232 2234 2239 2241 2243 2245 2253 2257 2260 2263
    ## [466] 2266 2269 2278 2283 2285 2289 2292 2295 2297 2299 2302 2304 2306 2309 2311
    ## [481] 2314 2324 2329 2332 2338

    # In case that the first and/or the last gene had one repetitive element: 
    if (repetitive.last[1] < repetitive.first[1]) {

        repetitive.last <- repetitive.last[-1]

        if (repetitive.last[length(repetitive.last)] < 
            repetitive.first[length(repetitive.first)]) {

            repetitive.last <- c(repetitive.last, 
                                 repetitive.first[length(repetitive.last)])

        }

        

    }

    # Explore the indices (all three vectors have to have the same length!)
    repetitive.geneid

    ##   [1]    3    8   12   20   46   50   54   58   66   70   73   77   79   82   84
    ##  [16]   89  104  112  116  120  124  126  128  133  136  141  147  151  157  163
    ##  [31]  165  167  169  176  183  191  193  199  206  209  219  221  223  232  236
    ##  [46]  240  245  251  255  257  267  271  274  285  297  300  305  309  316  319
    ##  [61]  325  330  334  342  347  353  358  361  366  368  370  372  375  384  386
    ##  [76]  390  393  396  402  410  412  415  418  421  424  426  428  431  433  435
    ##  [91]  438  440  443  454  459  465  467  471  476  480  488  514  518  522  526
    ## [106]  534  538  541  545  547  550  552  557  572  580  584  588  592  594  596
    ## [121]  601  604  609  615  619  625  631  633  635  637  644  651  659  661  667
    ## [136]  674  677  687  689  691  700  704  708  713  719  723  725  735  739  742
    ## [151]  753  765  768  773  777  784  787  793  798  802  810  815  821  826  829
    ## [166]  834  836  838  840  843  852  854  858  861  864  870  878  880  883  886
    ## [181]  889  892  894  896  899  901  903  906  908  911  922  927  933  935  939
    ## [196]  944  948  956  982  986  990  994 1002 1006 1009 1013 1015 1018 1020 1025
    ## [211] 1040 1048 1052 1056 1060 1062 1064 1069 1072 1077 1083 1087 1093 1099 1101
    ## [226] 1103 1105 1112 1119 1127 1129 1135 1142 1145 1155 1157 1159 1168 1172 1176
    ## [241] 1181 1187 1191 1193 1203 1207 1210 1221 1233 1236 1241 1245 1252 1255 1261
    ## [256] 1266 1270 1278 1283 1289 1294 1297 1302 1304 1306 1308 1311 1320 1322 1326
    ## [271] 1329 1332 1338 1346 1348 1351 1354 1357 1360 1362 1364 1367 1369 1371 1374
    ## [286] 1376 1379 1390 1395 1401 1403 1407 1412 1416 1424 1450 1454 1458 1462 1470
    ## [301] 1474 1477 1481 1483 1486 1488 1493 1508 1516 1520 1524 1528 1530 1532 1537
    ## [316] 1540 1545 1551 1555 1561 1567 1569 1571 1573 1580 1587 1595 1597 1603 1610
    ## [331] 1613 1623 1625 1627 1636 1640 1644 1649 1655 1659 1661 1671 1675 1678 1689
    ## [346] 1701 1704 1709 1713 1720 1723 1729 1734 1738 1746 1751 1757 1762 1765 1770
    ## [361] 1772 1774 1776 1779 1788 1790 1794 1797 1800 1806 1814 1816 1819 1822 1825
    ## [376] 1828 1830 1832 1835 1837 1839 1842 1844 1847 1858 1863 1869 1871 1875 1880
    ## [391] 1884 1892 1918 1922 1926 1930 1938 1942 1945 1949 1951 1954 1956 1961 1976
    ## [406] 1984 1988 1992 1996 1998 2000 2005 2008 2013 2019 2023 2029 2035 2037 2039
    ## [421] 2041 2048 2055 2063 2065 2071 2078 2081 2091 2093 2095 2104 2108 2112 2117
    ## [436] 2123 2127 2129 2139 2143 2146 2157 2169 2172 2177 2181 2188 2191 2197 2202
    ## [451] 2206 2214 2219 2225 2230 2233 2238 2240 2242 2244 2247 2256 2258 2262 2265
    ## [466] 2268 2274 2282 2284 2287 2290 2293 2296 2298 2300 2303 2305 2307 2310 2312
    ## [481] 2315 2326 2331 2337 2339

    repetitive.first

    ##   [1]    4    9   13   21   47   51   55   59   67   71   74   78   80   83   85
    ##  [16]   90  105  113  117  121  125  127  129  134  137  142  148  152  158  164
    ##  [31]  166  168  170  177  184  192  194  200  207  210  220  222  224  233  237
    ##  [46]  241  246  252  256  258  268  272  275  286  298  301  306  310  317  320
    ##  [61]  326  331  335  343  348  354  359  362  367  369  371  373  376  385  387
    ##  [76]  391  394  397  403  411  413  416  419  422  425  427  429  432  434  436
    ##  [91]  439  441  444  455  460  466  468  472  477  481  489  515  519  523  527
    ## [106]  535  539  542  546  548  551  553  558  573  581  585  589  593  595  597
    ## [121]  602  605  610  616  620  626  632  634  636  638  645  652  660  662  668
    ## [136]  675  678  688  690  692  701  705  709  714  720  724  726  736  740  743
    ## [151]  754  766  769  774  778  785  788  794  799  803  811  816  822  827  830
    ## [166]  835  837  839  841  844  853  855  859  862  865  871  879  881  884  887
    ## [181]  890  893  895  897  900  902  904  907  909  912  923  928  934  936  940
    ## [196]  945  949  957  983  987  991  995 1003 1007 1010 1014 1016 1019 1021 1026
    ## [211] 1041 1049 1053 1057 1061 1063 1065 1070 1073 1078 1084 1088 1094 1100 1102
    ## [226] 1104 1106 1113 1120 1128 1130 1136 1143 1146 1156 1158 1160 1169 1173 1177
    ## [241] 1182 1188 1192 1194 1204 1208 1211 1222 1234 1237 1242 1246 1253 1256 1262
    ## [256] 1267 1271 1279 1284 1290 1295 1298 1303 1305 1307 1309 1312 1321 1323 1327
    ## [271] 1330 1333 1339 1347 1349 1352 1355 1358 1361 1363 1365 1368 1370 1372 1375
    ## [286] 1377 1380 1391 1396 1402 1404 1408 1413 1417 1425 1451 1455 1459 1463 1471
    ## [301] 1475 1478 1482 1484 1487 1489 1494 1509 1517 1521 1525 1529 1531 1533 1538
    ## [316] 1541 1546 1552 1556 1562 1568 1570 1572 1574 1581 1588 1596 1598 1604 1611
    ## [331] 1614 1624 1626 1628 1637 1641 1645 1650 1656 1660 1662 1672 1676 1679 1690
    ## [346] 1702 1705 1710 1714 1721 1724 1730 1735 1739 1747 1752 1758 1763 1766 1771
    ## [361] 1773 1775 1777 1780 1789 1791 1795 1798 1801 1807 1815 1817 1820 1823 1826
    ## [376] 1829 1831 1833 1836 1838 1840 1843 1845 1848 1859 1864 1870 1872 1876 1881
    ## [391] 1885 1893 1919 1923 1927 1931 1939 1943 1946 1950 1952 1955 1957 1962 1977
    ## [406] 1985 1989 1993 1997 1999 2001 2006 2009 2014 2020 2024 2030 2036 2038 2040
    ## [421] 2042 2049 2056 2064 2066 2072 2079 2082 2092 2094 2096 2105 2109 2113 2118
    ## [436] 2124 2128 2130 2140 2144 2147 2158 2170 2173 2178 2182 2189 2192 2198 2203
    ## [451] 2207 2215 2220 2226 2231 2234 2239 2241 2243 2245 2248 2257 2259 2263 2266
    ## [466] 2269 2275 2283 2285 2288 2291 2294 2297 2299 2301 2304 2306 2308 2311 2313
    ## [481] 2316 2327 2332 2338 2340

    repetitive.last

    ##   [1]    4   10   13   44   48   52   55   64   67   72   76   78   81   83   88
    ##  [16]   92  106  114  117  123  125  127  131  134  140  146  149  153  159  164
    ##  [31]  166  168  173  182  188  192  197  204  208  211  220  222  225  234  239
    ##  [46]  244  249  253  256  266  268  273  281  295  298  304  306  310  317  321
    ##  [61]  329  332  335  346  352  356  360  362  367  369  371  373  381  385  388
    ##  [76]  391  394  397  406  411  413  417  420  423  425  427  430  432  434  437
    ##  [91]  439  442  452  457  460  466  468  472  478  481  512  516  520  523  532
    ## [106]  535  540  544  546  549  551  556  560  574  582  585  591  593  595  599
    ## [121]  602  608  614  617  621  627  632  634  636  641  650  656  660  665  672
    ## [136]  676  679  688  690  693  702  707  712  717  721  724  734  736  741  749
    ## [151]  763  766  772  774  778  785  789  797  800  803  814  820  824  828  830
    ## [166]  835  837  839  841  849  853  856  859  862  865  874  879  881  885  888
    ## [181]  891  893  895  898  900  902  905  907  910  920  925  928  934  936  940
    ## [196]  946  949  980  984  988  991 1000 1003 1008 1012 1014 1017 1019 1024 1028
    ## [211] 1042 1050 1053 1059 1061 1063 1067 1070 1076 1082 1085 1089 1095 1100 1102
    ## [226] 1104 1109 1118 1124 1128 1133 1140 1144 1147 1156 1158 1161 1170 1175 1180
    ## [241] 1185 1189 1192 1202 1204 1209 1217 1231 1234 1240 1242 1246 1253 1257 1265
    ## [256] 1268 1271 1282 1288 1292 1296 1298 1303 1305 1307 1309 1317 1321 1324 1327
    ## [271] 1330 1333 1342 1347 1349 1353 1356 1359 1361 1363 1366 1368 1370 1373 1375
    ## [286] 1378 1388 1393 1396 1402 1404 1408 1414 1417 1448 1452 1456 1459 1468 1471
    ## [301] 1476 1480 1482 1485 1487 1492 1496 1510 1518 1521 1527 1529 1531 1535 1538
    ## [316] 1544 1550 1553 1557 1563 1568 1570 1572 1577 1586 1592 1596 1601 1608 1612
    ## [331] 1615 1624 1626 1629 1638 1643 1648 1653 1657 1660 1670 1672 1677 1685 1699
    ## [346] 1702 1708 1710 1714 1721 1725 1733 1736 1739 1750 1756 1760 1764 1766 1771
    ## [361] 1773 1775 1777 1785 1789 1792 1795 1798 1801 1810 1815 1817 1821 1824 1827
    ## [376] 1829 1831 1834 1836 1838 1841 1843 1846 1856 1861 1864 1870 1872 1876 1882
    ## [391] 1885 1916 1920 1924 1927 1936 1939 1944 1948 1950 1953 1955 1960 1964 1978
    ## [406] 1986 1989 1995 1997 1999 2003 2006 2012 2018 2021 2025 2031 2036 2038 2040
    ## [421] 2045 2054 2060 2064 2069 2076 2080 2083 2092 2094 2097 2106 2111 2116 2121
    ## [436] 2125 2128 2138 2140 2145 2153 2167 2170 2176 2178 2182 2189 2193 2201 2204
    ## [451] 2207 2218 2224 2228 2232 2234 2239 2241 2243 2245 2253 2257 2260 2263 2266
    ## [466] 2269 2278 2283 2285 2289 2292 2295 2297 2299 2302 2304 2306 2309 2311 2314
    ## [481] 2324 2329 2332 2338 2338

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
    ## 4 ENSG00000159140            21           44            24
    ## 5 ENSG00000182568            47           48             2
    ## 6 ENSG00000081913            51           52             2

    dim(trfnumTable)

    ## [1] 485   4

    # Explore the output data frame
    glimpse(trfnumTable)

    ## Rows: 485
    ## Columns: 4
    ## $ GENEID        <chr> "ENSG00000162676", "ENSG00000171791", "ENSG00000069011",…
    ## $ Repeats_First <int> 4, 9, 13, 21, 47, 51, 55, 59, 67, 71, 74, 78, 80, 83, 85…
    ## $ Repeats_Last  <dbl> 4, 10, 13, 44, 48, 52, 55, 64, 67, 72, 76, 78, 81, 83, 8…
    ## $ Repeats_Count <dbl> 1, 2, 1, 24, 2, 2, 1, 6, 1, 2, 3, 1, 2, 1, 4, 3, 2, 2, 1…

    # Export as a csv file
    write.csv(trfnumTable, "TRF_demo.csv")

### Session Info

    sessionInfo()

    ## R version 4.0.3 (2020-10-10)
    ## Platform: x86_64-conda-linux-gnu (64-bit)
    ## Running under: Ubuntu 20.04.2 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /home/mira/miniconda3/envs/snakemake_r/lib/libopenblasp-r0.3.12.so
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
    ## [19] rvest_1.0.0          colorspace_2.0-0     htmltools_0.5.1.1   
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
    ## [52] compiler_4.0.3       rlang_0.4.10         grid_4.0.3          
    ## [55] rstudioapi_0.13      rappdirs_0.3.3       rmarkdown_2.7       
    ## [58] gtable_0.3.0         curl_4.3             DBI_1.1.1           
    ## [61] R6_2.5.0             lubridate_1.7.10     knitr_1.31          
    ## [64] fastmap_1.1.0        bit_4.0.4            utf8_1.2.1          
    ## [67] stringi_1.5.3        parallel_4.0.3       Rcpp_1.0.6          
    ## [70] vctrs_0.3.8          dbplyr_2.1.1         tidyselect_1.1.1    
    ## [73] xfun_0.20
