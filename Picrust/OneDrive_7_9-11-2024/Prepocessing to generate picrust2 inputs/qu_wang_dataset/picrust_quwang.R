### Change the input files #################################################### qu wang
seqtab4 = readRDS("/Users/if-mcs-1247c/Desktop/Osteoporosis/Quwang\ _dataset/seqtabnochim_qi.rds")
taxa4 = readRDS("/Users/if-mcs-1247c/Desktop/Osteoporosis/Quwang\ _dataset/taxa_qirdb.rds")          
map4 <- read.csv("/Users/if-mcs-1247c/Desktop/Osteoporosis/Quwang\ _dataset/xu_metadata1.csv")              
###############################################################################

###Make a phyloseq Object
ps4 <- phyloseq(otu_table(seqtab4, taxa_are_rows=FALSE),
                tax_table(taxa4))
row.names(map4) = map4$X.Sample_ID
sample_data(ps4) = sample_data(map4) 
ps4


##remove cyanobacteria and mitochondria reads
ps4= subset_taxa(ps4, Phylum!="Cyanobacteria")
ps4= subset_taxa(ps4, Family!="Mitochondria")
ps4

###Run pime to carry out prevalence filtering
library(pime)
pime.oob.error(ps4, "Groups")
per_variable_obj= pime.split.by.variable(ps4, "Groups")
prevalences=pime.prevalence(per_variable_obj)
set.seed(42)
best.prev=pime.best.prevalence(prevalences, "Groups")
imp40=best.prev$Importance$`Prevalence 40`
View(imp40)
prevalence.40 = prevalences$`40`


library(Biostrings)
#######convert Asv  sequence to fasta format
sequences <- Biostrings::DNAStringSet(taxa_names(prevalence.40))
names(sequences) <- taxa_names(prevalence.40)

#####merge the phyloseq object with sequenc
prevalence.40<- merge_phyloseq(prevalence.40, sequences)
prevalence.40


#####Replace sequence column with Asv
taxa_names(prevalence.40) <- paste0("ASV", seq(ntaxa(prevalence.40)))
View(tax_table(prevalence.40))
rs = refseq(prevalence.40)

### Export the sequence in fasta file
Biostrings::writeXStringSet(rs, "/Users/if-mcs-1247c/Desktop/Osteoporosis/Quwang\ _dataset/Phyloseq_picrust/quwang.fasta")

##transpose the otu_table
Otu_table_picrust2 = t(otu_table(prevalence.40))
View(Otu_table_picrust2)

###Covert the ASV row to a colum by adding a Column name
Otu_table_picrust2= as.data.frame(Otu_table_picrust2)
Otu_table_picrust2<- tibble::rownames_to_column(Otu_table_picrust2, "ASV_ID")
View(Otu_table_picrust2)

###Export the OTU_table
test= Otu_table_picrust2
View(test)
write.table(test, "/Users/if-mcs-1247c/Desktop/Osteoporosis/Quwang\ _dataset/Phyloseq_picrust/test10.tsv", quote = FALSE, sep ="\t")
