###Load wang_2017 data
wang_2017 = read.table("/Users/if-mcs-1247c/Desktop/Osteoporosis/wang_et\ al_2017_dataset/PICRUST/second_picrust/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv",sep = "\t", header = T) 
View(wang_2017)

### Add a row to the wang_2017 df so that the function column can be removed
rownames(wang_2017) <- paste0("ENZ", 1:nrow(wang_2017))
View(wang_2017)

######## Remove the function column 
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
dim(wang_2017_n)

#######Remove trimereric enzymes
fi = c(233,307,322,323,338,394,298, 299, 301,463,476,477,481,644,695,903)
wang_2017_n_f = wang_2017_n[-fi,]
View(wang_2017_n_f)


#####################################
###Load korean_dataset
Korean = read.table("/Users/if-mcs-1247c/Desktop/Osteoporosis/Korean_Osteoporosis_dataset/PICRUST/second\ picrust/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv",sep = "\t", header = T) 
View(Korean)

### Add a row to the korean df so that the function column can be removed
rownames(Korean) <- paste0("ENZ", 1:nrow(Korean))
View(Korean)

######## Remove the function column 
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
filt = c(249,271,273,274,348,366,422,487,496,498,504,601,701,905)

Korean_n_f = Korean_n[-filt,]
View(Korean_n_f)


###################### 
####load qu_wang dataset
quwang = read.table("/Users/if-mcs-1247c/Desktop/Osteoporosis/Quwang\ _dataset/PICRUST/second_picrust/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv",sep = "\t", header = T) 
View(quwang)

### Add a row to th quwang dataset so that the function column can be removed
rownames(quwang) <- paste0("ENZ", 1:nrow(quwang)) 
View(quwang)

######## Remove the function column 
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
fil = c(112,113,121,137,138,168,169,176,176,307,349,350,351)
quwang_n_f = quwang_n[-fil,]
View(quwang_n_f)

########### 

#load Wang 2022_dataset
wang_2022 = read.table("/Users/if-mcs-1247c/Desktop/Osteoporosis/wang\ et\ al._2022/PICRUST/second_wang-2022_picrust/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv" ,sep = "\t", header = T) 
View(wang_2022)

### Add a row to th wang 2022 dataset so that the function column can be removed

rownames(wang_2022) <- paste0("ENZ", 1:nrow(wang_2022)) 
View(wang_2022)

######## Remove the function column 
wang_2022 = wang_2022[,-1]
View(wang_2022)

######### Make the rows column, so it can be removed and the dataframe can be merged with others by the description column
#### Make the rows column
wang_2022 = as.data.frame(wang_2022)
wang_2022<- tibble::rownames_to_column(wang_2022, "row_names")
View(wang_2022)
dim(wang_2022)
#Remove the row_names column
wang_2022_n = wang_2022[, -1]
View(wang_2022_n)

#######Filter out enzymes that are joint together
f= c(120,147,148,155,181,278,316,319,320,322,327,526,477,722)
wang_2022_n_f = wang_2022_n[-f,]
View(wang_2022_n_f)



#############3
##load chein et al., 2021  dataset
chein_2021 = read.table("/Users/if-mcs-1247c/Desktop/Osteoporosis/Chein_2021/PICRUST/second_picrust/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv",sep = "\t", header = T) 

### Add a row to th chein_2021 dataset so that the function column can be removed

rownames(chein_2021) <- paste0("ENZ", 1:nrow(chein_2021)) 
View(chein_2021)

######## Remove the function column 
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
## Remove enzymes joint together 
filts = c(221,245,344,349,350,356,380,508,509,511,518,613,663,712,906)
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
View(merge2)
dim(merge2)

############### Merge chein_2021 dataset  and merge2  by the  description column
dim(merge2)
dim(chein_2021_n_f)
merge3 = merge(chein_2021_n_f, merge2, by = ("description"), all =TRUE)
dim(merge3)
View(merge3)

########## merge merge3and merge1 by description column
dim(merge1)
dim(merge3)
merge4 = merge(merge3, merge1, by = ("description"), all =TRUE)
dim(merge4)
View(merge4)

## Remove enzymes joint together 
merge4_new = merge4[-c(47,48,101,92,106,107,108,114,115,119,199,200,209,212,328,329,330,618,624,328,1496,1507,1508,471,94,93,9497,98,171,340,412,621,213,639,1249,1331,1332,1352,1436,1491,1495,1498,1317,1321,1323,1328,1497,1503,1508,1528,1856,1862,1873,1875,1964,1864,1871,1893,1965,1973,1974,1977,1958,1975),]
View(merge4_new)
dim(merge4_new)
write.csv(merge4_new, "/Users/if-mcs-1247c/Desktop/Osteoporosis/stamp_input_enzy1.33.csv")


#######Add a rownames to  merge 4 so that the TAX table and OTU table can be xtracted and merged by common row names
merge4_new = as.data.frame(merge4_new)
rownames(merge4_new) <- paste0("ENZ", 1:nrow(merge4_new)) 
View(merge4_new)

## Add an index column to the dataframe so that the tax table can be extracted
merge4_new$index <- paste0("Enzymes", merge4_new$index)
View(merge4_new)

### Extract TAX table
TAX_table = merge4_new[,c(353,1)]

View(TAX_table)
dim(TAX_table)

##Extract Otu_table
OTU_TABLE = merge4_new[,-c(353,1)]
View(OTU_TABLE)
dim(OTU_TABLE)

#####  Make NAs equal to zero
OTU_TABLE[is.na(OTU_TABLE)] <-0

### Make theOTU and TAX table matrix 
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

set.seed(2125)
ps_all_r= rarefy_even_depth(ps_all, sample.size = 10510, replace = T)
### View the rarified phyloseq object
all.reads = data.table(as(sample_data(ps_all_r), "data.frame"),
                       TotalReads = sample_sums(ps_all_r), keep.rownames = TRUE)
View(all.reads)

data.frame = as(sample_data(ps_all_r), "data.frame")
plyr::count(data.frame, "Groups")
ps_all_r

#######Extract otu and tax table
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
write.csv(otu_tax_new, "/Users/if-mcs-1247c/Desktop/Osteoporosis/stamp_input_enzy100.csv")



