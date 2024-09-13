###Load wang_2017 data
wang_2017 = read.table("/Users/if-mcs-1247c/Desktop/Osteoporosis/wang_et\ al_2017_dataset/PICRUST/picrust2_wang_2017_out_pipeline/pathways_out/path_abun_unstrat_descrip.tsv",sep = "\t", header = T) 
View(wang_2017)

## Add a row to the wang_2017 df so that the pathway column can be removed
rownames(wang_2017) <- paste0("ENZ", 1:nrow(wang_2017))
View(wang_2017)

######## Remove the pathway column 
wang_2017_n = wang_2017[,-1]
View(wang_2017_n)

######### Make the rows column, so it can be removed and the dataframe can be merged with others by the description column
#### Make the rows column
wang_2017_n = as.data.frame(wang_2017_n)
wang_2017_n<- tibble::rownames_to_column(wang_2017_n, "row_names")
View(wang_2017_n)

#Remove the row_names column
wang_2017_n = wang_2017_n[, -1]
View(wang_2017_n)

#######Remove pathways that are joint together
fi = c(165,240,248)
wang_2017_n_f = wang_2017_n[-fi,]
View(wang_2017_n_f)


###Load wang_2022 data
wang_2022 = read.table("/Users/if-mcs-1247c/Desktop/Osteoporosis/wang\ et\ al._2022/PICRUST/picrust2_wang2022_out_pipeline/pathways_out/path_abun_unstrat_descrip.tsv",sep = "\t", header = T) 
View(wang_2022)

## Add a row to the wang_2022 df so that the pathway column can be removed
rownames(wang_2022) <- paste0("ENZ", 1:nrow(wang_2022))
View(wang_2022)

######## Remove the pathway  column 
wang_2022_n = wang_2022[,-1]
View(wang_2022_n)

######### Make the rows column, so it can be removed and the dataframe can be merged with others by the description column
#### Make the rows column
wang_2022_n = as.data.frame(wang_2022_n)
wang_2022_n<- tibble::rownames_to_column(wang_2022_n, "row_names")
View(wang_2022_n)
#Remove the row_names column
wang_2022_n = wang_2022_n[, -1]
View(wang_2022_n)

#######Remove pathways that are joint together
fil = c(157, 231,240)
wang_2022_n_f = wang_2022_n[-fil,]
View(wang_2022_n_f)


#####################################
###Load korean_dataset
Korean = read.table("/Users/if-mcs-1247c/Desktop/Osteoporosis/Korean_Osteoporosis_dataset/PICRUST/picrust2_korean_out_pipeline/pathways_out/path_abun_unstrat_descrip.tsv",sep = "\t", header = T) 
View(Korean)

### Add a row to the korean df so that the pathway column can be removed
rownames(Korean) <- paste0("ENZ", 1:nrow(Korean))
View(Korean)

######## Remove the pathway  column 
Korean_n = Korean[,-1]
View(Korean_n)
dim(Korean_n)

######### Make the rows column, so it can be removed and the dataframe can be merged with others by the description column
#### Make the rows column
Korean_n = as.data.frame(Korean_n)
Korean_n<- tibble::rownames_to_column(Korean_n, "row_names")
View(Korean_n)

#Remove the row_names column
Korean_n = Korean_n[, -1]
View(Korean_n)

#######Filter out enzymes that are joint together
filt = c(136,204,212)

Korean_n_f = Korean_n[-filt,]
View(Korean_n_f)


###################### 
####load qu_wang dataset
quwang = read.table("/Users/if-mcs-1247c/Desktop/Osteoporosis/Quwang\ _dataset/PICRUST/picrust2_quwang_out_pipeline/pathways_out/path_abun_unstrat_descrip.tsv",sep = "\t", header = T) 
View(quwang)

### Add a row to th quwang dataset so that the pathway column can be removed

rownames(quwang) <- paste0("ENZ", 1:nrow(quwang)) 
View(quwang)

######## Remove the pathway column 
quwang_n = quwang[,-1]
View(quwang_n)
dim(quwang_n)
######### Make the rows column, so it can be removed and the dataframe can be merged with others by the description column
#### Make the rows column
quwang_n = as.data.frame(quwang_n)
quwang_n<- tibble::rownames_to_column(quwang_n, "row_names")
View(quwang_n)

#Remove the row_names column
quwang_n = quwang_n[, -1]
View(quwang_n)

#######Filter out enzymes that are joint together
fils = c(147,218,226)
quwang_n_f = quwang_n[-fils,]
View(quwang_n_f)

##load chein et al., 2021  dataset
chein_2021 = read.table("/Users/if-mcs-1247c/Desktop/Osteoporosis/Chein_2021/PICRUST/picrust2_chein_2021_out_pipeline/pathways_out/path_abun_unstrat_descrip.tsv", sep = "\t", header = T) 
View(chein_2021)
### Add a row to th chein_2021 dataset so that the Pathway column can be removed

rownames(chein_2021) <- paste0("ENZ", 1:nrow(chein_2021)) 
View(chein_2021)

######## Remove the pathway  column 
chein_2021_n = chein_2021[,-1]
View(chein_2021_n)

######### Make the rows column, so it can be removed and the dataframe can be merged with others by the description column
#### Make the rows column
chein_2021_n = as.data.frame(chein_2021_n)
chein_2021_n<- tibble::rownames_to_column(chein_2021_n, "row_names")
View(chein_2021_n)

#Remove the row_names column
chein_2021_n= chein_2021_n[, -1]
View(chein_2021_n)

filts = c(117,185,177)
chein_2021_n_f = chein_2021_n[-filts,]
View(chein_2021_n_f)


############### Merge wang_2017 and korea dataset by the  description column
dim(Korean_n_f)
dim(wang_2017_n_f)
merge1 = merge(wang_2017_n_f, Korean_n_f,  by = ("description"), all =TRUE)
dim(merge1)
View(merge1)


############### Merge wang_2022 and quwang dataset by the  description column
dim(wang_2022_n_f)
dim(quwang_n_f)
merge2 = merge(wang_2022_n_f, quwang_n_f,  by = ("description"), all =TRUE)
dim(merge2)
View(merge2)


############### Merge chein_2021 dataset  and merge2  by the  description column
dim(merge1)
dim(chein_2021_n_f)
merge3 = merge(chein_2021_n_f, merge1, by = ("description"), all =TRUE)
dim(merge3)
View(merge3)

########## merge merge2 and merge 3 by description column
dim(merge2)
dim(merge3)
merge4 = merge(merge3, merge2, by = ("description"), all =TRUE)
dim(merge4)
View(merge4)

#####Convert NAs to zero
merge4[is.na(merge4)] <-0
View(merge4)

#######Add a rownames to  merge 4 so that the TAX table and OTU table can be extracted and merged by common row names
merge4 = as.data.frame(merge4)
rownames(merge4) <- paste0("ENZ", 1:nrow(merge4)) 
View(merge4)

## Add an index column to the dataframe so that the tax table can be extracted
merge4$index <- paste0("Enzymes", merge4$index)
View(merge4)

### Extract TAX table
TAX_table = merge4[,c(353,1)]

View(TAX_table)
dim(TAX_table)

##Extract Otu_table
OTU_TABLE = merge4[,-c(353,1)]
View(OTU_TABLE)
dim(OTU_TABLE)

#####  Make NAs equal to zero
OTU_TABLE[is.na(OTU_TABLE)] <-0

### Make the OTU and TAX table matrix 
Tax_table =as.matrix(TAX_table)

Otu_table = as.matrix(OTU_TABLE)
class(Otu_table) <- "numeric"

###Make Phyloseq object
library(phyloseq)
ps_all <- phyloseq(otu_table(Otu_table, taxa_are_rows=TRUE),
                   tax_table(Tax_table))
ps_all

##Load metadata
sd_new = read.csv("/Users/if-mcs-1247c/Desktop/line/lol/sd_new1.csv")
dim(sd_new)

row.names(sd_new) = sd_new$X.Sample_ID
sample_data(ps_all) = sample_data(sd_new) 
ps_all



###View the reads found in each sample
library(data.table)
all.reads = data.table(as(sample_data(ps_all), "data.frame"),
                       TotalReads = sample_sums(ps_all), keep.rownames = TRUE)
View(all.reads)

######Rarify the dataset
set.seed(2125)
ps_all_r= rarefy_even_depth(ps_all, sample.size = 3074, replace = T)

### View the rarified phyloseq object
all.reads = data.table(as(sample_data(ps_all_r), "data.frame"),
                       TotalReads = sample_sums(ps_all_r), keep.rownames = TRUE)
View(all.reads)

####View number of samples after rarefaction
data.frame = as(sample_data(ps_all_r), "data.frame")
plyr::count(data.frame, "Groups")
ps_all_r

#######Extract otu_ and tax table
otu = otu_table(ps_all_r)
tax = tax_table(ps_all_r)

#####Merge otu_table and tax_table by column
otu_tax = cbind(tax,otu)
View(otu_tax)

###Remove the first column
otu_tax_new = otu_tax[,-1]
View(otu_tax_new)

########## convert the rows to column and remove it
otu_tax_new= as.data.frame(otu_tax_new)
otu_tax_new<- tibble::rownames_to_column(otu_tax_new, "row_names")
View(otu_tax_new)

#Remove the row_names column
otu_tax_new= otu_tax_new[, -1]

View(otu_tax_new)


#####Export stamp output
write.csv(otu_tax_new, "/Users/if-mcs-1247c/Desktop/Osteoporosis/stamp_input_path.csv")

