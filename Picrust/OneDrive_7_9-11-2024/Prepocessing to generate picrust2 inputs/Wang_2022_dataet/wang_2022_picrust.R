### Change the input files #################################################### wang 2022
seqtab2 = readRDS("/Users/if-mcs-1247c/Desktop/Osteoporosis/wang\ et\ al._2022/seqtabnochim.rds")
taxa2 = readRDS("/Users/if-mcs-1247c/Desktop/Osteoporosis/wang\ et\ al._2022/taxa_rdp.rds")          
map2 <- read.csv("/Users/if-mcs-1247c/Desktop/Osteoporosis/wang\ et\ al._2022/OP_wang_2022_map.csv")
View(map2)

library(phyloseq)
library(microbiome)

###Make a phyloseq Object
ps2 <- phyloseq(otu_table(seqtab2, taxa_are_rows=FALSE),
                tax_table(taxa2))
row.names(map2) = map2$X.Sample_ID
sample_data(ps2) = sample_data(map2) 
ps2



###Run pime to carry out prevalence filtering
pime.oob.error(ps2, "Groups")
per_variable_obj= pime.split.by.variable(ps2, "Groups")
prevalences=pime.prevalence(per_variable_obj)
set.seed(42)
best.prev=pime.best.prevalence(prevalences, "Groups")
imp50=best.prev$Importance$`Prevalence 50`
View(imp50)
prevalence.50 = prevalences$`50`



library(Biostrings)

#######convert Asv  sequence to fasta format
sequences <- Biostrings::DNAStringSet(taxa_names(prevalence.50))
names(sequences) <- taxa_names(prevalence.50)

#####merge the phyloseq object with sequence
prevalence.50 <- merge_phyloseq(prevalence.50, sequences)
prevalence.50
View(otu_table(prevalence.50))

#####Replace sequence colum with Asv
taxa_names(prevalence.50) <- paste0("ASV", seq(ntaxa(prevalence.50)))
View(tax_table(prevalence.50))

rs = refseq(prevalence.50)

### Export the sequence in fasta file
Biostrings::writeXStringSet(rs, "/Users/if-mcs-1247c/Desktop/Osteoporosis/wang\ et\ al._2022/Phyloseq_picrust\tax_wang22.fasta")

##transpose the otu_table
Otu_table_picrust2 = t(otu_table(prevalence.50))
View(Otu_table_picrust2)

###Covert the ASV row to a colum by adding a Column name
Otu_table_picrust2= as.data.frame(Otu_table_picrust2)
Otu_table_picrust2<- tibble::rownames_to_column(Otu_table_picrust2, "ASV_ID")
View(Otu_table_picrust2)

###Export the OTU_table
test= Otu_table_picrust2
View(test)
write.table(test, "/Users/if-mcs-1247c/Desktop/Osteoporosis/wang\ et\ al._2022/Phyloseq_picrust\test1.tsv", quote = FALSE, sep ="\t")
