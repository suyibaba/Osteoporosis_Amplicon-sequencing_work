#### Read Input files
seqtab3 = readRDS("/Users/if-mcs-1247c/Desktop/Osteoporosis/Korean_Osteoporosis_dataset/seqtab.nochim3")
taxa3 = readRDS("/Users/if-mcs-1247c/Desktop/Osteoporosis/Korean_Osteoporosis_dataset/taxa_rdp_korea.rds")          
map3 <- read.csv("/Users/if-mcs-1247c/Desktop/Osteoporosis/Korean_Osteoporosis_dataset/Osteoporosis_mapfile_new_f.csv") 


View(taxa3)
View(seqtab3)

#### Transpose the ASV_table so that it can be merged by row names with the Tax_table
seqtab3= t(seqtab3)
View(seqtab3)

###Make a phyloseq Object
ps3 <- phyloseq(otu_table(seqtab3, taxa_are_rows=TRUE),
                tax_table(taxa3))
row.names(map3) = map3$X.SampleID
sample_data(ps3) = sample_data(map3) 
ps3



###Run pime to carry out prevalence filtering
library(pime)
pime.oob.error(ps3_i, "Groups")

per_variable_obj= pime.split.by.variable(ps3_i, "Groups")
per_variable_obj

prevalences =pime.prevalence(per_variable_obj)
head(prevalences)

set.seed(42)
best.prev=pime.best.prevalence(prevalences, "Groups")

imp30=best.prev$Importance$`Prevalence 30`
View(imp30)
prevalence.30 = prevalences$`30`


#######convert Asv  sequence to fasta format
library(Biostrings)
sequences <- Biostrings::DNAStringSet(taxa_names(prevalence.30))
names(sequences) <- taxa_names(prevalence.30)

#####merge the phyloseq object with sequence
prevalence.30<- merge_phyloseq(prevalence.30, sequences)
prevalence.30
View(tax_table(prevalence.30))

#####Replace sequence colum with Asv
taxa_names(prevalence.30) <- paste0("ASV", seq(ntaxa(prevalence.30)))
View(tax_table(prevalence.30))

rs = refseq(prevalence.30)
rs

### Export the sequence in fasta file
Biostrings::writeXStringSet(rs, "/Users/if-mcs-1247c/Desktop/Osteoporosis/Korean_Osteoporosis_dataset/Phyloseq-picrust/korea.fasta")

##transpose the otu_table
Otu_table_picrust2 = t(otu_table(prevalence.30))
View(Otu_table_picrust2)

###Covert the ASV row to a column  by adding a Column name
Otu_table_picrust2= as.data.frame(Otu_table_picrust2)
Otu_table_picrust2<- tibble::rownames_to_column(Otu_table_picrust2, "ASV_ID")
View(Otu_table_picrust2)


###Export the OTU_table
test= Otu_table_picrust2

##### transpose the otu_table
test = t(test)
View(test)
write.table(test, "/Users/if-mcs-1247c/Desktop/Osteoporosis/Korean_Osteoporosis_dataset/Phyloseq-picrust/test12.tsv", quote = FALSE, sep = "\t")
