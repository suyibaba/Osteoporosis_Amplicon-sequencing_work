library(phyloseq)
seqtab = readRDS("/Users/if-mcs-1247c/Desktop/Osteoporosis/wang_et\ al_2017_dataset/seqtab.nochim2.rds")
taxa = readRDS("/Users/if-mcs-1247c/Desktop/Osteoporosis/wang_et\ al_2017_dataset/taxa_rdp.rds")          
map=  read.csv("/Users/if-mcs-1247c/Desktop/Osteoporosis/wang_et\ al_2017_dataset/wang_metadata_2017.csv")
View(map)


###Make a phyloseq Object
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
               tax_table(taxa))
row.names(map) = map$X.Sample_ID
sample_data(ps) = sample_data(map) 
ps



###Run pime to carry out prevalence filtering
library(pime)
pime.oob.error(ps, "Groups")
per_variable_obj= pime.split.by.variable(ps, "Groups")
prevalences=pime.prevalence(per_variable_obj)
set.seed(42)
best.prev=pime.best.prevalence(prevalences, "Groups")
imp50=best.prev$Importance$`Prevalence 50`
View(imp50)
prevalence.50 = prevalences$`50`
prevalence.50




library(Biostrings)

#######convert Asv  sequence to fasta format
sequences <- Biostrings::DNAStringSet(taxa_names(prevalence.50))
names(sequences) <- taxa_names(prevalence.50)

#####merge the phyloseq object with sequence
prevalence.50 <- merge_phyloseq(prevalence.50, sequences)
prevalence.50

#####Replace sequence colum with Asv
taxa_names(prevalence.50) <- paste0("ASV", seq(ntaxa(prevalence.50)))
View(tax_table(prevalence.50))

rs = refseq(prevalence.50)
rs

####write the refseq object in fasta format
 Biostrings::writeXStringSet(rs, "/Users/if-mcs-1247c/Desktop/Osteoporosis/wang_et\ al_2017_dataset/Phyloseq_picrust2/tax_wang17.fasta")


####transpose otu_table
Otu_table_picrust2 = t(otu_table(prevalence.50))
View(otu_table(prevalence.50))

###Covert the ASV row to a colum by adding a Column name
Otu_table_picrust2= as.data.frame(Otu_table_picrust2)
Otu_table_picrust2<- tibble::rownames_to_column(Otu_table_picrust2, "ASV_ID")
View(Otu_table_picrust2)

###Export the OTU_table
test = Otu_table_picrust2
View(test)
write.table(test, "/Users/if-mcs-1247c/Desktop/Osteoporosis/wang_et\ al_2017_dataset/Phyloseq_picrust2/test_new1.tsv", quote = FALSE, sep ="\t")
