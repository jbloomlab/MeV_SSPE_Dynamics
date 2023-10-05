## ---------------------------
##
## Script name: `rule add_metadata`
##
## Purpose of script: read in file with given accession and add in the requisite metadata
##
## Author: Will Hannon
##
## Date Created: 2020-06-26
##
## Copyright (c) Will Hannon, 2020
## Email: wwh22@uw.edu
##
## ---------------------------

require(tidyverse)

## ---------------------------

## ==== Get file paths from snakemake object ==== ##

vcf.filepath = snakemake@input[[1]] 
metadata.filepath = snakemake@params[[1]] 

## ==== Boolean to see if file is empty and there are no columns ==== ##
info = file.info(vcf.filepath)
empty = (info$size == 0)

## ==== Import Metadata, same table used to import samples from `Run Selector w/ minimal changes` ==== ##
#
# Add in any constrains on metadata here.
#
# Download and format metedata. But first, check if the vcf is empty..
if (!empty){
  
  metadata.df = read_csv(metadata.filepath) %>% rename(Accession = "Run")
  
  ## ==== Import and process the vcf -> txt files ==== ##
  
  # cols to split `ANN` vcf field into
  cols_into = c("ALT_Allele", 
                "Effect",
                "Putative_Impact",
                "Gene_Name",
                "Gene_ID",
                "Feature_Type",
                "Feature_ID",
                "Transcript_Biotype",
                "Rank",
                "HGVS.c",
                "HGVS.p",
                "cDNA_position",
                "CDS_position", 
                "Protein_position",
                "Distance_to_feature",
                "Errors")
  
  # extract metadata from filepath 
  sample.info = strsplit(basename(vcf.filepath), "[.]")[[1]]
  sample.name = sample.info[1]

  # If the the second element is either 'fwd' or 'rev', then sample.caller is the third element
  if (sample.info[2] == "Positive" | sample.info[2] == "Negative"){
    sample.direction = sample.info[2]
    sample.caller = sample.info[3]
  } else{
    sample.direction = NULL
    sample.caller = sample.info[2]
  }
    

  # different extraction depending on columns in the original VCF file
  if (sample.caller == "lofreq"){
    # read in vcf table
    sample.df = read_tsv(vcf.filepath) %>% 
      # split into cols with regex
      separate(col = ANN, sep = "\\|", remove = F, into = cols_into) %>% 
      select(POS, REF,  ALT,
             AF,  DP,
             Effect, Gene_Name, AA_Change="HGVS.p") %>% 
      # clean amino acid change column
      mutate(AA_Change = gsub("^p\\.", "", AA_Change)) %>% 
      # simplify the effect column
      mutate(Effect = case_when(Effect  == "missense_variant" ~ "Missense",
                                Effect  == "synonymous_variant" ~ "Synonymous",
                                Effect  ==  "upstream_gene_variant"  ~ "Synonymous",
                                Effect  ==  "downstream_gene_variant"  ~ "Synonymous",
                                Effect  ==  "stop_retained_variant"  ~ "Synonymous",
                                Effect  ==  "stop_lost" ~ "Nonsense",
                                Effect  ==  "start_lost"  ~ "Nonsense",
                                Effect  ==  "stop_gained" ~ "Nonsense")) %>% 
      # add sample ID
      mutate(Accession = sample.name) %>%
      mutate(Caller = sample.caller) 

  } else if (sample.caller == "varscan"){
    
    # read in vcf table
    sample.df = read_tsv(vcf.filepath) %>%
      # split into cols with regex
      separate(col = ANN, sep = "\\|", remove = F, into = cols_into) %>% 
      select(POS, REF,  ALT,
             AF="Sample1.FREQ",  DP="Sample1.DP",
             Effect, Gene_Name, AA_Change="HGVS.p") %>% 
      # clean amino acid change column
      mutate(AA_Change = gsub("^p\\.", "", AA_Change)) %>% 
      # clean allele freq column
      mutate(AF = gsub("%", "", AF)) %>% 
      # simplify the effect column
      mutate(Effect = case_when(Effect  == "missense_variant" ~ "Missense",
                                Effect  == "synonymous_variant" ~ "Synonymous",
                                Effect  ==  "upstream_gene_variant"  ~ "Synonymous",
                                Effect  ==  "downstream_gene_variant"  ~ "Synonymous",
                                Effect  ==  "stop_retained_variant"  ~ "Synonymous",
                                Effect  ==  "stop_lost" ~ "Nonsense",
                                Effect  ==  "start_lost"  ~ "Nonsense",
                                Effect  ==  "stop_gained" ~ "Nonsense")) %>% 
      # add sample ID
      mutate(Accession = sample.name) %>%
      mutate(Caller = sample.caller) 
    
    sample.df$AF = as.numeric(sample.df$AF) / 100
    
  } else if (sample.caller == "ivar"){
    
    # read in vcf table
    sample.df = read_tsv(vcf.filepath) %>%
      mutate(Effect = if_else(REF_AA == ALT_AA, "Synonymous", "Missense"),
             AA_Change = paste0(REF_AA, "-", ALT_AA),
             Gene_Name = NA) %>% 
      select(POS, REF,  ALT,
             AF="ALT_FREQ",  DP="TOTAL_DP",
             Effect, Gene_Name, AA_Change) %>% 
      # add sample ID
      mutate(Accession = sample.name) %>%
      mutate(Caller = sample.caller) %>% 
      distinct()
    
  }
  
  ## ==== Join the variant data with the metadata by 'Run' id ==== ##
  
  # Join with metadata
  sample.df = left_join(sample.df, metadata.df, by = "Accession")

  # If sample.direction isn't NULL, add it to the table
  if (!is.null(sample.direction)){
    sample.df = sample.df %>% mutate(Direction = sample.direction)
  }
    
  # Write out to a file
  write.csv(sample.df, file = snakemake@output[[1]], row.names=FALSE, sep="\t")
  
} else{
  # Write an empty output if file is empty
  cat(NULL,file=snakemake@output[[1]])
}


## ==== END ==== ##











