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

    ## [1] 19770

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

    ## [1] 6156    1

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

    ## [1] 4212    2

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
    ## [2329]  0  1  2  0  1  2  3  4  5  0  1  0  1  2  3  0  1  2  3  4  0  0  1  2
    ## [2353]  0  1  2  3  4  5  6  7  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
    ## [2377]  0  0  0  0  0  0  0  0  1  2  0  0  1  2  0  0  1  2  0  1  2  3  0  0
    ## [2401]  0  0  0  0  1  2  0  1  2  3  0  0  1  0  0  0  1  0  1  0  0  1  0  1
    ## [2425]  0  0  0  0  1  0  0  0  1  2  3  4  5  6  7  8  9 10 11 12  0  0  1  2
    ## [2449]  3  4  5  6  0  0  1  2  0  1  2  3  0  0  0  1  0  1  0  1  0  0  0  1
    ## [2473]  2  0  1  2  0  0  0  0  1  0  0  0  0  0  1  0  0  1  2  0  0  1  2  3
    ## [2497]  4  0  0  1  2  3  4  0  1  0  1  0  1  0  0  0  0  1  2  3  0  0  0  0
    ## [2521]  0  0  1  0  0  0  0  0  1  2  3  0  1  0  0  0  0  1  2  0  0  0  0  0
    ## [2545]  1  2  0  0  1  0  0  1  2  3  4  5  6  7  8  0  1  0  1  0  0  1  2  3
    ## [2569]  4  5  6  7  0  0  1  2  0  0  0  1  0  0  0  0  1  0  0  0  0  1  2  0
    ## [2593]  0  1  2  0  1  0  0  0  0  0  0  0  0  0  1  0  1  2  3  0  0  1  0  0
    ## [2617]  0  0  0  0  0  1  2  3  4  0  0  0  0  0  0  0  0  0  0  1  2  0  1  2
    ## [2641]  0  0  0  0  1  0  1  2  3  0  1  2  3  4  5  6  0  1  2  0  0  1  2  3
    ## [2665]  4  0  0  0  0  1  0  0  1  2  0  1  2  3  4  5  6  7  0  0  0  0  1  0
    ## [2689]  0  0  0  0  1  0  0  0  1  2  0  0  1  0  1  2  3  4  0  1  0  1  0  1
    ## [2713]  0  1  2  0  0  0  0  0  0  1  2  3  0  1  0  0  1  2  0  1  2  0  1  2
    ## [2737]  0  1  2  3  4  5  0  0  0  0  1  2  3  4  0  1  0  1  2  0  0  1  0  0
    ## [2761]  1  0  0  1  0  1  0  1  0  0  1  0  1  0  1  0  0  1  0  1  0  0  1  0
    ## [2785]  0  0  0  0  0  0  0  0  1  2  0  0  0  1  2  0  1  2  3  4  5  0  1  0
    ## [2809]  1  2  3  0  1  2  3  4  0  0  1  2  0  1  2  3  4  5  6  7  0  0  0  0
    ## [2833]  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  2  0  0
    ## [2857]  1  2  0  0  1  2  0  1  2  3  0  0  0  0  0  0  1  2  0  1  2  3  0  0
    ## [2881]  1  0  0  0  1  0  1  0  0  1  0  1  0  0  0  0  1  0  0  0  1  2  3  4
    ## [2905]  5  6  7  8  9 10 11 12  0  0  1  2  3  4  5  6  0  0  1  2  0  1  2  3
    ## [2929]  0  0  0  1  0  1  0  1  0  0  0  1  2  0  1  2  0  0  0  0  1  0  0  0
    ## [2953]  0  0  1  0  0  1  2  0  0  1  2  3  4  0  0  1  2  3  4  0  1  0  1  0
    ## [2977]  1  0  0  0  0  1  2  3  0  0  0  0  0  0  1  0  0  0  0  0  1  2  3  0
    ## [3001]  1  0  0  0  0  1  2  0  0  0  0  0  1  2  0  0  1  0  0  1  2  3  4  5
    ## [3025]  6  7  8  0  1  0  1  0  0  1  2  3  4  5  6  7  0  0  1  2  0  0  0  1
    ## [3049]  0  0  0  0  1  0  0  0  0  1  2  0  0  1  2  0  1  0  0  0  0  0  0  0
    ## [3073]  0  0  1  0  1  2  3  0  0  1  0  0  0  0  0  0  0  1  2  3  4  0  0  0
    ## [3097]  0  0  0  0  0  0  0  1  2  0  1  2  0  0  0  0  1  0  1  2  3  0  1  2
    ## [3121]  3  4  5  6  0  1  2  0  0  1  2  3  4  0  0  0  0  1  0  0  1  2  0  1
    ## [3145]  2  3  4  5  6  7  0  0  0  0  1  0  0  0  0  0  1  0  0  0  1  2  0  0
    ## [3169]  1  0  1  2  3  4  0  1  0  1  0  1  0  1  2  0  0  0  0  0  0  1  2  3
    ## [3193]  0  1  0  0  1  2  0  1  2  0  1  2  0  1  2  3  4  5  0  0  0  0  1  2
    ## [3217]  3  4  0  1  0  1  2  0  0  1  0  0  1  0  0  1  0  1  0  1  0  0  1  0
    ## [3241]  1  0  1  0  0  1  0  1  0  0  1  0  0  0  0  0  0  0  0  0  1  2  0  0
    ## [3265]  0  1  2  0  1  2  3  4  5  0  1  0  1  2  3  0  1  2  3  4  0  0  1  2
    ## [3289]  0  1  2  3  4  5  6  7  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
    ## [3313]  0  0  0  0  0  0  0  0  1  2  0  0  1  2  0  0  1  2  0  1  2  3  0  0
    ## [3337]  0  0  0  0  1  2  0  1  2  3  0  0  1  0  0  0  1  0  1  0  0  1  0  1
    ## [3361]  0  0  0  0  1  0  0  0  1  2  3  4  5  6  7  8  9 10 11 12  0  0  1  2
    ## [3385]  3  4  5  6  0  0  1  2  0  1  2  3  0  0  0  1  0  1  0  1  0  0  0  1
    ## [3409]  2  0  1  2  0  0  0  0  1  0  0  0  0  0  1  0  0  1  2  0  0  1  2  3
    ## [3433]  4  0  0  1  2  3  4  0  1  0  1  0  1  0  0  0  0  1  2  3  0  0  0  0
    ## [3457]  0  0  1  0  0  0  0  0  1  2  3  0  1  0  0  0  0  1  2  0  0  0  0  0
    ## [3481]  1  2  0  0  1  0  0  1  2  3  4  5  6  7  8  0  1  0  1  0  0  1  2  3
    ## [3505]  4  5  6  7  0  0  1  2  0  0  0  1  0  0  0  0  1  0  0  0  0  1  2  0
    ## [3529]  0  1  2  0  1  0  0  0  0  0  0  0  0  0  1  0  1  2  3  0  0  1  0  0
    ## [3553]  0  0  0  0  0  1  2  3  4  0  0  0  0  0  0  0  0  0  0  1  2  0  1  2
    ## [3577]  0  0  0  0  1  0  1  2  3  0  1  2  3  4  5  6  0  1  2  0  0  1  2  3
    ## [3601]  4  0  0  0  0  1  0  0  1  2  0  1  2  3  4  5  6  7  0  0  0  0  1  0
    ## [3625]  0  0  0  0  1  0  0  0  1  2  0  0  1  0  1  2  3  4  0  1  0  1  0  1
    ## [3649]  0  1  2  0  0  0  0  0  0  1  2  3  0  1  0  0  1  2  0  1  2  0  1  2
    ## [3673]  0  1  2  3  4  5  0  0  0  0  1  2  3  4  0  1  0  1  2  0  0  1  0  0
    ## [3697]  1  0  0  1  0  1  0  1  0  0  1  0  1  0  1  0  0  1  0  1  0  0  1  0
    ## [3721]  0  0  0  0  0  0  0  0  1  2  0  0  0  1  2  0  1  2  3  4  5  0  1  0
    ## [3745]  1  2  3  0  1  2  3  4  0  0  1  2  0  1  2  3  4  5  6  7  0  0  0  0
    ## [3769]  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  2  0  0
    ## [3793]  1  2  0  0  1  2  0  1  2  3  0  0  0  0  0  0  1  2  0  1  2  3  0  0
    ## [3817]  1  0  0  0  1  0  1  0  0  1  0  1  0  0  0  0  1  0  0  0  1  2  3  4
    ## [3841]  5  6  7  8  9 10 11 12  0  0  1  2  3  4  5  6  0  0  1  2  0  1  2  3
    ## [3865]  0  0  0  1  0  1  0  1  0  0  0  1  2  0  1  2  0  0  0  0  1  0  0  0
    ## [3889]  0  0  1  0  0  1  2  0  0  1  2  3  4  0  0  1  2  3  4  0  1  0  1  0
    ## [3913]  1  0  0  0  0  1  2  3  0  0  0  0  0  0  1  0  0  0  0  0  1  2  3  0
    ## [3937]  1  0  0  0  0  1  2  0  0  0  0  0  1  2  0  0  1  0  0  1  2  3  4  5
    ## [3961]  6  7  8  0  1  0  1  0  0  1  2  3  4  5  6  7  0  0  1  2  0  0  0  1
    ## [3985]  0  0  0  0  1  0  0  0  0  1  2  0  0  1  2  0  1  0  0  0  0  0  0  0
    ## [4009]  0  0  1  0  1  2  3  0  0  1  0  0  0  0  0  0  0  1  2  3  4  0  0  0
    ## [4033]  0  0  0  0  0  0  0  1  2  0  1  2  0  0  0  0  1  0  1  2  3  0  1  2
    ## [4057]  3  4  5  6  0  1  2  0  0  1  2  3  4  0  0  0  0  1  0  0  1  2  0  1
    ## [4081]  2  3  4  5  6  7  0  0  0  0  1  0  0  0  0  0  1  0  0  0  1  2  0  0
    ## [4105]  1  0  1  2  3  4  0  1  0  1  0  1  0  1  2  0  0  0  0  0  0  1  2  3
    ## [4129]  0  1  0  0  1  2  0  1  2  0  1  2  0  1  2  3  4  5  0  0  0  0  1  2
    ## [4153]  3  4  0  1  0  1  2  0  0  1  0  0  1  0  0  1  0  1  0  1  0  0  1  0
    ## [4177]  1  0  1  0  0  1  0  1  0  0  1  0  0  0  0  0  0  0  0  0  1  2  0  0
    ## [4201]  0  1  2  0  1  2  3  4  5  0  1  0

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
    ## [2329]  3  0  0  1  0  0  0  0  0  1  0  1  0  0  0  1  0  0  0  0  1  2  0  0
    ## [2353]  1  0  0  0  0  0  0  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
    ## [2377] 17 18 19 20 21 22 23 24  0  0  1  2  0  0  1  2  0  0  1  0  0  0  1  2
    ## [2401]  3  4  5  6  0  0  1  0  0  0  1  2  0  1  2  3  0  1  0  1  2  0  1  0
    ## [2425]  1  2  3  4  0  1  2  3  0  0  0  0  0  0  0  0  0  0  0  0  1  2  0  0
    ## [2449]  0  0  0  0  1  2  0  0  1  0  0  0  1  2  3  0  1  0  1  0  1  2  3  0
    ## [2473]  0  1  0  0  1  2  3  4  0  1  2  3  4  5  0  1  2  0  0  1  2  0  0  0
    ## [2497]  0  1  2  0  0  0  0  1  0  1  0  1  0  1  2  3  4  0  0  0  1  2  3  4
    ## [2521]  5  6  0  1  2  3  4  5  0  0  0  1  0  1  2  3  4  0  0  1  2  3  4  5
    ## [2545]  0  0  1  2  0  1  2  0  0  0  0  0  0  0  0  1  0  1  0  1  2  0  0  0
    ## [2569]  0  0  0  0  1  2  0  0  1  2  3  0  1  2  3  4  0  1  2  3  4  0  0  1
    ## [2593]  2  0  0  1  0  1  2  3  4  5  6  7  8  9  0  1  0  0  0  1  2  0  1  2
    ## [2617]  3  4  5  6  7  0  0  0  0  1  2  3  4  5  6  7  8  9 10  0  0  1  0  0
    ## [2641]  1  2  3  4  0  1  0  0  0  1  0  0  0  0  0  0  1  0  0  1  2  0  0  0
    ## [2665]  0  1  2  3  4  0  1  2  0  0  1  0  0  0  0  0  0  0  1  2  3  4  0  1
    ## [2689]  2  3  4  5  0  1  2  3  0  0  1  2  0  1  0  0  0  0  1  0  1  0  1  0
    ## [2713]  1  0  0  1  2  3  4  5  6  0  0  0  1  0  1  2  0  0  1  0  0  1  0  0
    ## [2737]  1  0  0  0  0  0  1  2  3  4  0  0  0  0  1  0  1  0  0  1  2  0  1  2
    ## [2761]  0  1  2  0  1  0  1  0  1  2  0  1  0  1  0  1  2  0  1  0  1  2  0  1
    ## [2785]  2  3  4  5  6  7  8  9  0  0  1  2  3  0  0  1  0  0  0  0  0  1  0  1
    ## [2809]  0  0  0  1  0  0  0  0  1  2  0  0  1  0  0  0  0  0  0  0  1  2  3  4
    ## [2833]  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24  0  0  1  2
    ## [2857]  0  0  1  2  0  0  1  0  0  0  1  2  3  4  5  6  0  0  1  0  0  0  1  2
    ## [2881]  0  1  2  3  0  1  0  1  2  0  1  0  1  2  3  4  0  1  2  3  0  0  0  0
    ## [2905]  0  0  0  0  0  0  0  0  1  2  0  0  0  0  0  0  1  2  0  0  1  0  0  0
    ## [2929]  1  2  3  0  1  0  1  0  1  2  3  0  0  1  0  0  1  2  3  4  0  1  2  3
    ## [2953]  4  5  0  1  2  0  0  1  2  0  0  0  0  1  2  0  0  0  0  1  0  1  0  1
    ## [2977]  0  1  2  3  4  0  0  0  1  2  3  4  5  6  0  1  2  3  4  5  0  0  0  1
    ## [3001]  0  1  2  3  4  0  0  1  2  3  4  5  0  0  1  2  0  1  2  0  0  0  0  0
    ## [3025]  0  0  0  1  0  1  0  1  2  0  0  0  0  0  0  0  1  2  0  0  1  2  3  0
    ## [3049]  1  2  3  4  0  1  2  3  4  0  0  1  2  0  0  1  0  1  2  3  4  5  6  7
    ## [3073]  8  9  0  1  0  0  0  1  2  0  1  2  3  4  5  6  7  0  0  0  0  1  2  3
    ## [3097]  4  5  6  7  8  9 10  0  0  1  0  0  1  2  3  4  0  1  0  0  0  1  0  0
    ## [3121]  0  0  0  0  1  0  0  1  2  0  0  0  0  1  2  3  4  0  1  2  0  0  1  0
    ## [3145]  0  0  0  0  0  0  1  2  3  4  0  1  2  3  4  5  0  1  2  3  0  0  1  2
    ## [3169]  0  1  0  0  0  0  1  0  1  0  1  0  1  0  0  1  2  3  4  5  6  0  0  0
    ## [3193]  1  0  1  2  0  0  1  0  0  1  0  0  1  0  0  0  0  0  1  2  3  4  0  0
    ## [3217]  0  0  1  0  1  0  0  1  2  0  1  2  0  1  2  0  1  0  1  0  1  2  0  1
    ## [3241]  0  1  0  1  2  0  1  0  1  2  0  1  2  3  4  5  6  7  8  9  0  0  1  2
    ## [3265]  3  0  0  1  0  0  0  0  0  1  0  1  0  0  0  1  0  0  0  0  1  2  0  0
    ## [3289]  1  0  0  0  0  0  0  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
    ## [3313] 17 18 19 20 21 22 23 24  0  0  1  2  0  0  1  2  0  0  1  0  0  0  1  2
    ## [3337]  3  4  5  6  0  0  1  0  0  0  1  2  0  1  2  3  0  1  0  1  2  0  1  0
    ## [3361]  1  2  3  4  0  1  2  3  0  0  0  0  0  0  0  0  0  0  0  0  1  2  0  0
    ## [3385]  0  0  0  0  1  2  0  0  1  0  0  0  1  2  3  0  1  0  1  0  1  2  3  0
    ## [3409]  0  1  0  0  1  2  3  4  0  1  2  3  4  5  0  1  2  0  0  1  2  0  0  0
    ## [3433]  0  1  2  0  0  0  0  1  0  1  0  1  0  1  2  3  4  0  0  0  1  2  3  4
    ## [3457]  5  6  0  1  2  3  4  5  0  0  0  1  0  1  2  3  4  0  0  1  2  3  4  5
    ## [3481]  0  0  1  2  0  1  2  0  0  0  0  0  0  0  0  1  0  1  0  1  2  0  0  0
    ## [3505]  0  0  0  0  1  2  0  0  1  2  3  0  1  2  3  4  0  1  2  3  4  0  0  1
    ## [3529]  2  0  0  1  0  1  2  3  4  5  6  7  8  9  0  1  0  0  0  1  2  0  1  2
    ## [3553]  3  4  5  6  7  0  0  0  0  1  2  3  4  5  6  7  8  9 10  0  0  1  0  0
    ## [3577]  1  2  3  4  0  1  0  0  0  1  0  0  0  0  0  0  1  0  0  1  2  0  0  0
    ## [3601]  0  1  2  3  4  0  1  2  0  0  1  0  0  0  0  0  0  0  1  2  3  4  0  1
    ## [3625]  2  3  4  5  0  1  2  3  0  0  1  2  0  1  0  0  0  0  1  0  1  0  1  0
    ## [3649]  1  0  0  1  2  3  4  5  6  0  0  0  1  0  1  2  0  0  1  0  0  1  0  0
    ## [3673]  1  0  0  0  0  0  1  2  3  4  0  0  0  0  1  0  1  0  0  1  2  0  1  2
    ## [3697]  0  1  2  0  1  0  1  0  1  2  0  1  0  1  0  1  2  0  1  0  1  2  0  1
    ## [3721]  2  3  4  5  6  7  8  9  0  0  1  2  3  0  0  1  0  0  0  0  0  1  0  1
    ## [3745]  0  0  0  1  0  0  0  0  1  2  0  0  1  0  0  0  0  0  0  0  1  2  3  4
    ## [3769]  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24  0  0  1  2
    ## [3793]  0  0  1  2  0  0  1  0  0  0  1  2  3  4  5  6  0  0  1  0  0  0  1  2
    ## [3817]  0  1  2  3  0  1  0  1  2  0  1  0  1  2  3  4  0  1  2  3  0  0  0  0
    ## [3841]  0  0  0  0  0  0  0  0  1  2  0  0  0  0  0  0  1  2  0  0  1  0  0  0
    ## [3865]  1  2  3  0  1  0  1  0  1  2  3  0  0  1  0  0  1  2  3  4  0  1  2  3
    ## [3889]  4  5  0  1  2  0  0  1  2  0  0  0  0  1  2  0  0  0  0  1  0  1  0  1
    ## [3913]  0  1  2  3  4  0  0  0  1  2  3  4  5  6  0  1  2  3  4  5  0  0  0  1
    ## [3937]  0  1  2  3  4  0  0  1  2  3  4  5  0  0  1  2  0  1  2  0  0  0  0  0
    ## [3961]  0  0  0  1  0  1  0  1  2  0  0  0  0  0  0  0  1  2  0  0  1  2  3  0
    ## [3985]  1  2  3  4  0  1  2  3  4  0  0  1  2  0  0  1  0  1  2  3  4  5  6  7
    ## [4009]  8  9  0  1  0  0  0  1  2  0  1  2  3  4  5  6  7  0  0  0  0  1  2  3
    ## [4033]  4  5  6  7  8  9 10  0  0  1  0  0  1  2  3  4  0  1  0  0  0  1  0  0
    ## [4057]  0  0  0  0  1  0  0  1  2  0  0  0  0  1  2  3  4  0  1  2  0  0  1  0
    ## [4081]  0  0  0  0  0  0  1  2  3  4  0  1  2  3  4  5  0  1  2  3  0  0  1  2
    ## [4105]  0  1  0  0  0  0  1  0  1  0  1  0  1  0  0  1  2  3  4  5  6  0  0  0
    ## [4129]  1  0  1  2  0  0  1  0  0  1  0  0  1  0  0  0  0  0  1  2  3  4  0  0
    ## [4153]  0  0  1  0  1  0  0  1  2  0  1  2  0  1  2  0  1  0  1  0  1  2  0  1
    ## [4177]  0  1  0  1  2  0  1  0  1  2  0  1  2  3  4  5  6  7  8  9  0  0  1  2
    ## [4201]  3  0  0  1  0  0  0  0  0  1  0  1

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
    ## [481] 2315 2326 2331 2337 2339 2343 2348 2352 2360 2386 2390 2394 2398 2406 2410
    ## [496] 2413 2417 2419 2422 2424 2429 2444 2452 2456 2460 2464 2466 2468 2473 2476
    ## [511] 2481 2487 2491 2497 2503 2505 2507 2509 2516 2523 2531 2533 2539 2546 2549
    ## [526] 2559 2561 2563 2572 2576 2580 2585 2591 2595 2597 2607 2611 2614 2625 2637
    ## [541] 2640 2645 2649 2656 2659 2665 2670 2674 2682 2687 2693 2698 2701 2706 2708
    ## [556] 2710 2712 2715 2724 2726 2730 2733 2736 2742 2750 2752 2755 2758 2761 2764
    ## [571] 2766 2768 2771 2773 2775 2778 2780 2783 2794 2799 2805 2807 2811 2816 2820
    ## [586] 2828 2854 2858 2862 2866 2874 2878 2881 2885 2887 2890 2892 2897 2912 2920
    ## [601] 2924 2928 2932 2934 2936 2941 2944 2949 2955 2959 2965 2971 2973 2975 2977
    ## [616] 2984 2991 2999 3001 3007 3014 3017 3027 3029 3031 3040 3044 3048 3053 3059
    ## [631] 3063 3065 3075 3079 3082 3093 3105 3108 3113 3117 3124 3127 3133 3138 3142
    ## [646] 3150 3155 3161 3166 3169 3174 3176 3178 3180 3183 3192 3194 3198 3201 3204
    ## [661] 3210 3218 3220 3223 3226 3229 3232 3234 3236 3239 3241 3243 3246 3248 3251
    ## [676] 3262 3267 3273 3275 3279 3284 3288 3296 3322 3326 3330 3334 3342 3346 3349
    ## [691] 3353 3355 3358 3360 3365 3380 3388 3392 3396 3400 3402 3404 3409 3412 3417
    ## [706] 3423 3427 3433 3439 3441 3443 3445 3452 3459 3467 3469 3475 3482 3485 3495
    ## [721] 3497 3499 3508 3512 3516 3521 3527 3531 3533 3543 3547 3550 3561 3573 3576
    ## [736] 3581 3585 3592 3595 3601 3606 3610 3618 3623 3629 3634 3637 3642 3644 3646
    ## [751] 3648 3651 3660 3662 3666 3669 3672 3678 3686 3688 3691 3694 3697 3700 3702
    ## [766] 3704 3707 3709 3711 3714 3716 3719 3730 3735 3741 3743 3747 3752 3756 3764
    ## [781] 3790 3794 3798 3802 3810 3814 3817 3821 3823 3826 3828 3833 3848 3856 3860
    ## [796] 3864 3868 3870 3872 3877 3880 3885 3891 3895 3901 3907 3909 3911 3913 3920
    ## [811] 3927 3935 3937 3943 3950 3953 3963 3965 3967 3976 3980 3984 3989 3995 3999
    ## [826] 4001 4011 4015 4018 4029 4041 4044 4049 4053 4060 4063 4069 4074 4078 4086
    ## [841] 4091 4097 4102 4105 4110 4112 4114 4116 4119 4128 4130 4134 4137 4140 4146
    ## [856] 4154 4156 4159 4162 4165 4168 4170 4172 4175 4177 4179 4182 4184 4187 4198
    ## [871] 4203 4209 4211

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
    ## [481] 2316 2327 2332 2338 2340 2344 2349 2353 2361 2387 2391 2395 2399 2407 2411
    ## [496] 2414 2418 2420 2423 2425 2430 2445 2453 2457 2461 2465 2467 2469 2474 2477
    ## [511] 2482 2488 2492 2498 2504 2506 2508 2510 2517 2524 2532 2534 2540 2547 2550
    ## [526] 2560 2562 2564 2573 2577 2581 2586 2592 2596 2598 2608 2612 2615 2626 2638
    ## [541] 2641 2646 2650 2657 2660 2666 2671 2675 2683 2688 2694 2699 2702 2707 2709
    ## [556] 2711 2713 2716 2725 2727 2731 2734 2737 2743 2751 2753 2756 2759 2762 2765
    ## [571] 2767 2769 2772 2774 2776 2779 2781 2784 2795 2800 2806 2808 2812 2817 2821
    ## [586] 2829 2855 2859 2863 2867 2875 2879 2882 2886 2888 2891 2893 2898 2913 2921
    ## [601] 2925 2929 2933 2935 2937 2942 2945 2950 2956 2960 2966 2972 2974 2976 2978
    ## [616] 2985 2992 3000 3002 3008 3015 3018 3028 3030 3032 3041 3045 3049 3054 3060
    ## [631] 3064 3066 3076 3080 3083 3094 3106 3109 3114 3118 3125 3128 3134 3139 3143
    ## [646] 3151 3156 3162 3167 3170 3175 3177 3179 3181 3184 3193 3195 3199 3202 3205
    ## [661] 3211 3219 3221 3224 3227 3230 3233 3235 3237 3240 3242 3244 3247 3249 3252
    ## [676] 3263 3268 3274 3276 3280 3285 3289 3297 3323 3327 3331 3335 3343 3347 3350
    ## [691] 3354 3356 3359 3361 3366 3381 3389 3393 3397 3401 3403 3405 3410 3413 3418
    ## [706] 3424 3428 3434 3440 3442 3444 3446 3453 3460 3468 3470 3476 3483 3486 3496
    ## [721] 3498 3500 3509 3513 3517 3522 3528 3532 3534 3544 3548 3551 3562 3574 3577
    ## [736] 3582 3586 3593 3596 3602 3607 3611 3619 3624 3630 3635 3638 3643 3645 3647
    ## [751] 3649 3652 3661 3663 3667 3670 3673 3679 3687 3689 3692 3695 3698 3701 3703
    ## [766] 3705 3708 3710 3712 3715 3717 3720 3731 3736 3742 3744 3748 3753 3757 3765
    ## [781] 3791 3795 3799 3803 3811 3815 3818 3822 3824 3827 3829 3834 3849 3857 3861
    ## [796] 3865 3869 3871 3873 3878 3881 3886 3892 3896 3902 3908 3910 3912 3914 3921
    ## [811] 3928 3936 3938 3944 3951 3954 3964 3966 3968 3977 3981 3985 3990 3996 4000
    ## [826] 4002 4012 4016 4019 4030 4042 4045 4050 4054 4061 4064 4070 4075 4079 4087
    ## [841] 4092 4098 4103 4106 4111 4113 4115 4117 4120 4129 4131 4135 4138 4141 4147
    ## [856] 4155 4157 4160 4163 4166 4169 4171 4173 4176 4178 4180 4183 4185 4188 4199
    ## [871] 4204 4210 4212

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
    ## [481] 2314 2324 2329 2332 2338 2340 2344 2350 2353 2384 2388 2392 2395 2404 2407
    ## [496] 2412 2416 2418 2421 2423 2428 2432 2446 2454 2457 2463 2465 2467 2471 2474
    ## [511] 2480 2486 2489 2493 2499 2504 2506 2508 2513 2522 2528 2532 2537 2544 2548
    ## [526] 2551 2560 2562 2565 2574 2579 2584 2589 2593 2596 2606 2608 2613 2621 2635
    ## [541] 2638 2644 2646 2650 2657 2661 2669 2672 2675 2686 2692 2696 2700 2702 2707
    ## [556] 2709 2711 2713 2721 2725 2728 2731 2734 2737 2746 2751 2753 2757 2760 2763
    ## [571] 2765 2767 2770 2772 2774 2777 2779 2782 2792 2797 2800 2806 2808 2812 2818
    ## [586] 2821 2852 2856 2860 2863 2872 2875 2880 2884 2886 2889 2891 2896 2900 2914
    ## [601] 2922 2925 2931 2933 2935 2939 2942 2948 2954 2957 2961 2967 2972 2974 2976
    ## [616] 2981 2990 2996 3000 3005 3012 3016 3019 3028 3030 3033 3042 3047 3052 3057
    ## [631] 3061 3064 3074 3076 3081 3089 3103 3106 3112 3114 3118 3125 3129 3137 3140
    ## [646] 3143 3154 3160 3164 3168 3170 3175 3177 3179 3181 3189 3193 3196 3199 3202
    ## [661] 3205 3214 3219 3221 3225 3228 3231 3233 3235 3238 3240 3242 3245 3247 3250
    ## [676] 3260 3265 3268 3274 3276 3280 3286 3289 3320 3324 3328 3331 3340 3343 3348
    ## [691] 3352 3354 3357 3359 3364 3368 3382 3390 3393 3399 3401 3403 3407 3410 3416
    ## [706] 3422 3425 3429 3435 3440 3442 3444 3449 3458 3464 3468 3473 3480 3484 3487
    ## [721] 3496 3498 3501 3510 3515 3520 3525 3529 3532 3542 3544 3549 3557 3571 3574
    ## [736] 3580 3582 3586 3593 3597 3605 3608 3611 3622 3628 3632 3636 3638 3643 3645
    ## [751] 3647 3649 3657 3661 3664 3667 3670 3673 3682 3687 3689 3693 3696 3699 3701
    ## [766] 3703 3706 3708 3710 3713 3715 3718 3728 3733 3736 3742 3744 3748 3754 3757
    ## [781] 3788 3792 3796 3799 3808 3811 3816 3820 3822 3825 3827 3832 3836 3850 3858
    ## [796] 3861 3867 3869 3871 3875 3878 3884 3890 3893 3897 3903 3908 3910 3912 3917
    ## [811] 3926 3932 3936 3941 3948 3952 3955 3964 3966 3969 3978 3983 3988 3993 3997
    ## [826] 4000 4010 4012 4017 4025 4039 4042 4048 4050 4054 4061 4065 4073 4076 4079
    ## [841] 4090 4096 4100 4104 4106 4111 4113 4115 4117 4125 4129 4132 4135 4138 4141
    ## [856] 4150 4155 4157 4161 4164 4167 4169 4171 4174 4176 4178 4181 4183 4186 4196
    ## [871] 4201 4204 4210

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
    ## [481] 2315 2326 2331 2337 2339 2343 2348 2352 2360 2386 2390 2394 2398 2406 2410
    ## [496] 2413 2417 2419 2422 2424 2429 2444 2452 2456 2460 2464 2466 2468 2473 2476
    ## [511] 2481 2487 2491 2497 2503 2505 2507 2509 2516 2523 2531 2533 2539 2546 2549
    ## [526] 2559 2561 2563 2572 2576 2580 2585 2591 2595 2597 2607 2611 2614 2625 2637
    ## [541] 2640 2645 2649 2656 2659 2665 2670 2674 2682 2687 2693 2698 2701 2706 2708
    ## [556] 2710 2712 2715 2724 2726 2730 2733 2736 2742 2750 2752 2755 2758 2761 2764
    ## [571] 2766 2768 2771 2773 2775 2778 2780 2783 2794 2799 2805 2807 2811 2816 2820
    ## [586] 2828 2854 2858 2862 2866 2874 2878 2881 2885 2887 2890 2892 2897 2912 2920
    ## [601] 2924 2928 2932 2934 2936 2941 2944 2949 2955 2959 2965 2971 2973 2975 2977
    ## [616] 2984 2991 2999 3001 3007 3014 3017 3027 3029 3031 3040 3044 3048 3053 3059
    ## [631] 3063 3065 3075 3079 3082 3093 3105 3108 3113 3117 3124 3127 3133 3138 3142
    ## [646] 3150 3155 3161 3166 3169 3174 3176 3178 3180 3183 3192 3194 3198 3201 3204
    ## [661] 3210 3218 3220 3223 3226 3229 3232 3234 3236 3239 3241 3243 3246 3248 3251
    ## [676] 3262 3267 3273 3275 3279 3284 3288 3296 3322 3326 3330 3334 3342 3346 3349
    ## [691] 3353 3355 3358 3360 3365 3380 3388 3392 3396 3400 3402 3404 3409 3412 3417
    ## [706] 3423 3427 3433 3439 3441 3443 3445 3452 3459 3467 3469 3475 3482 3485 3495
    ## [721] 3497 3499 3508 3512 3516 3521 3527 3531 3533 3543 3547 3550 3561 3573 3576
    ## [736] 3581 3585 3592 3595 3601 3606 3610 3618 3623 3629 3634 3637 3642 3644 3646
    ## [751] 3648 3651 3660 3662 3666 3669 3672 3678 3686 3688 3691 3694 3697 3700 3702
    ## [766] 3704 3707 3709 3711 3714 3716 3719 3730 3735 3741 3743 3747 3752 3756 3764
    ## [781] 3790 3794 3798 3802 3810 3814 3817 3821 3823 3826 3828 3833 3848 3856 3860
    ## [796] 3864 3868 3870 3872 3877 3880 3885 3891 3895 3901 3907 3909 3911 3913 3920
    ## [811] 3927 3935 3937 3943 3950 3953 3963 3965 3967 3976 3980 3984 3989 3995 3999
    ## [826] 4001 4011 4015 4018 4029 4041 4044 4049 4053 4060 4063 4069 4074 4078 4086
    ## [841] 4091 4097 4102 4105 4110 4112 4114 4116 4119 4128 4130 4134 4137 4140 4146
    ## [856] 4154 4156 4159 4162 4165 4168 4170 4172 4175 4177 4179 4182 4184 4187 4198
    ## [871] 4203 4209 4211

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
    ## [481] 2316 2327 2332 2338 2340 2344 2349 2353 2361 2387 2391 2395 2399 2407 2411
    ## [496] 2414 2418 2420 2423 2425 2430 2445 2453 2457 2461 2465 2467 2469 2474 2477
    ## [511] 2482 2488 2492 2498 2504 2506 2508 2510 2517 2524 2532 2534 2540 2547 2550
    ## [526] 2560 2562 2564 2573 2577 2581 2586 2592 2596 2598 2608 2612 2615 2626 2638
    ## [541] 2641 2646 2650 2657 2660 2666 2671 2675 2683 2688 2694 2699 2702 2707 2709
    ## [556] 2711 2713 2716 2725 2727 2731 2734 2737 2743 2751 2753 2756 2759 2762 2765
    ## [571] 2767 2769 2772 2774 2776 2779 2781 2784 2795 2800 2806 2808 2812 2817 2821
    ## [586] 2829 2855 2859 2863 2867 2875 2879 2882 2886 2888 2891 2893 2898 2913 2921
    ## [601] 2925 2929 2933 2935 2937 2942 2945 2950 2956 2960 2966 2972 2974 2976 2978
    ## [616] 2985 2992 3000 3002 3008 3015 3018 3028 3030 3032 3041 3045 3049 3054 3060
    ## [631] 3064 3066 3076 3080 3083 3094 3106 3109 3114 3118 3125 3128 3134 3139 3143
    ## [646] 3151 3156 3162 3167 3170 3175 3177 3179 3181 3184 3193 3195 3199 3202 3205
    ## [661] 3211 3219 3221 3224 3227 3230 3233 3235 3237 3240 3242 3244 3247 3249 3252
    ## [676] 3263 3268 3274 3276 3280 3285 3289 3297 3323 3327 3331 3335 3343 3347 3350
    ## [691] 3354 3356 3359 3361 3366 3381 3389 3393 3397 3401 3403 3405 3410 3413 3418
    ## [706] 3424 3428 3434 3440 3442 3444 3446 3453 3460 3468 3470 3476 3483 3486 3496
    ## [721] 3498 3500 3509 3513 3517 3522 3528 3532 3534 3544 3548 3551 3562 3574 3577
    ## [736] 3582 3586 3593 3596 3602 3607 3611 3619 3624 3630 3635 3638 3643 3645 3647
    ## [751] 3649 3652 3661 3663 3667 3670 3673 3679 3687 3689 3692 3695 3698 3701 3703
    ## [766] 3705 3708 3710 3712 3715 3717 3720 3731 3736 3742 3744 3748 3753 3757 3765
    ## [781] 3791 3795 3799 3803 3811 3815 3818 3822 3824 3827 3829 3834 3849 3857 3861
    ## [796] 3865 3869 3871 3873 3878 3881 3886 3892 3896 3902 3908 3910 3912 3914 3921
    ## [811] 3928 3936 3938 3944 3951 3954 3964 3966 3968 3977 3981 3985 3990 3996 4000
    ## [826] 4002 4012 4016 4019 4030 4042 4045 4050 4054 4061 4064 4070 4075 4079 4087
    ## [841] 4092 4098 4103 4106 4111 4113 4115 4117 4120 4129 4131 4135 4138 4141 4147
    ## [856] 4155 4157 4160 4163 4166 4169 4171 4173 4176 4178 4180 4183 4185 4188 4199
    ## [871] 4204 4210 4212

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
    ## [481] 2324 2329 2332 2338 2340 2344 2350 2353 2384 2388 2392 2395 2404 2407 2412
    ## [496] 2416 2418 2421 2423 2428 2432 2446 2454 2457 2463 2465 2467 2471 2474 2480
    ## [511] 2486 2489 2493 2499 2504 2506 2508 2513 2522 2528 2532 2537 2544 2548 2551
    ## [526] 2560 2562 2565 2574 2579 2584 2589 2593 2596 2606 2608 2613 2621 2635 2638
    ## [541] 2644 2646 2650 2657 2661 2669 2672 2675 2686 2692 2696 2700 2702 2707 2709
    ## [556] 2711 2713 2721 2725 2728 2731 2734 2737 2746 2751 2753 2757 2760 2763 2765
    ## [571] 2767 2770 2772 2774 2777 2779 2782 2792 2797 2800 2806 2808 2812 2818 2821
    ## [586] 2852 2856 2860 2863 2872 2875 2880 2884 2886 2889 2891 2896 2900 2914 2922
    ## [601] 2925 2931 2933 2935 2939 2942 2948 2954 2957 2961 2967 2972 2974 2976 2981
    ## [616] 2990 2996 3000 3005 3012 3016 3019 3028 3030 3033 3042 3047 3052 3057 3061
    ## [631] 3064 3074 3076 3081 3089 3103 3106 3112 3114 3118 3125 3129 3137 3140 3143
    ## [646] 3154 3160 3164 3168 3170 3175 3177 3179 3181 3189 3193 3196 3199 3202 3205
    ## [661] 3214 3219 3221 3225 3228 3231 3233 3235 3238 3240 3242 3245 3247 3250 3260
    ## [676] 3265 3268 3274 3276 3280 3286 3289 3320 3324 3328 3331 3340 3343 3348 3352
    ## [691] 3354 3357 3359 3364 3368 3382 3390 3393 3399 3401 3403 3407 3410 3416 3422
    ## [706] 3425 3429 3435 3440 3442 3444 3449 3458 3464 3468 3473 3480 3484 3487 3496
    ## [721] 3498 3501 3510 3515 3520 3525 3529 3532 3542 3544 3549 3557 3571 3574 3580
    ## [736] 3582 3586 3593 3597 3605 3608 3611 3622 3628 3632 3636 3638 3643 3645 3647
    ## [751] 3649 3657 3661 3664 3667 3670 3673 3682 3687 3689 3693 3696 3699 3701 3703
    ## [766] 3706 3708 3710 3713 3715 3718 3728 3733 3736 3742 3744 3748 3754 3757 3788
    ## [781] 3792 3796 3799 3808 3811 3816 3820 3822 3825 3827 3832 3836 3850 3858 3861
    ## [796] 3867 3869 3871 3875 3878 3884 3890 3893 3897 3903 3908 3910 3912 3917 3926
    ## [811] 3932 3936 3941 3948 3952 3955 3964 3966 3969 3978 3983 3988 3993 3997 4000
    ## [826] 4010 4012 4017 4025 4039 4042 4048 4050 4054 4061 4065 4073 4076 4079 4090
    ## [841] 4096 4100 4104 4106 4111 4113 4115 4117 4125 4129 4132 4135 4138 4141 4150
    ## [856] 4155 4157 4161 4164 4167 4169 4171 4174 4176 4178 4181 4183 4186 4196 4201
    ## [871] 4204 4210 4210

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

    ## [1] 873   4

    # Explore the output data frame
    glimpse(trfnumTable)

    ## Rows: 873
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
