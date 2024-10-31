
#script takes a table of PCR genotypes and transforms into psuedoSNPs in the parent.call format used by LepMAP 3
#Designed for Y-markers used as part of a mass PCR genotyping project
#Requires .ped file used for linkage map construction and the table of results of the genotyping PCR for each sample
#input genotype file must be formatted as two columns, one with sample names (same as in .ped file) and the other as genotype call ("m" for male, "f" for female, "u" for unkown)
#default sex-linked region will be labelled as "Yregion" but this can be changed

### First, lets load in all the required files, and set up the required variables ###

args_in <- commandArgs(trailingOnly = TRUE)

genotype_file <- args_in[1]
ped_file <- args_in[2]
output_file <-  args_in[3]


if(length(args_in) < 4) {
  sex_linked_code <- "Yregion"
} else sex_linked_code <- (args_in[4])


sex_genotypes <- read.table(genotype_file, col.names = c("sample","genotype"))
ped_table <- read.table(ped_file)

XX_code <- "1.0 0 0 0 0 0 0 0 0 0"
XY_code <- "0 1.0 0 0 0 0 0 0 0 0"
Unkown_code <- "0.5 0.5 0 0 0 0 0 0 0 0"

output_table <- ped_table[2,]

#print(output_table)
print(sex_genotypes)

output_table[1,1] <- sex_linked_code
output_table[1,2] <- 100

### now loop over each column in the second row of the ped-table (which has the sample names), find the sample in the genotyping table, and add the relavent call ###

for(i in 3:ncol(ped_table)) {

    sample_index <- match(ped_table[2,i], sex_genotypes$sample)
    
    if(is.na(sample_index) == TRUE) {
    
        print(paste(ped_table[2,i], "not found in genotype file, treating as unkown"))
        output_table[1,i] <- Unkown_code
    }
    
    else if (sex_genotypes[sample_index,2] == "m") {
      
        print(paste(ped_table[2,i], "genotyped as male"))
        output_table[1,i] <- XY_code
    }
    
    else if (sex_genotypes[sample_index,2] == "f") {
      
        print(paste(ped_table[2,i], "genotyped as female"))
        output_table[1,i] <- XX_code
    }
    
    else if (sex_genotypes[sample_index,2] == "u") {
      
        print(paste(ped_table[2,i], "genotyped as unkown"))
        output_table[1,i] <- Unkown_code
    }

    else {
      
        print(paste(ped_table[2,i], "genotyped not parsed, treating as unkown"))
        output_table[1,i] <- Unkown_code
    }    

}

write.table(output_table,output_file,sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
