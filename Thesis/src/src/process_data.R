library(PharmacoGx)
library(dplyr)
library(tidyr)


############ Parameters #########

# get rid of liquid tumors
# drop_tissues <- c("Myeloid","Lymphoid")
drop_tissues <- c()
# define a threshold for the number of cell lines a tissue must have to be included
count_threshold <- 5
#drugs



# file sourcing
output.dir <- "../data/processed/"
PSet.baseDir <- "../data/raw/"
PSet.names <- c("gCSI","GDSC2","Broad")
PSet.fileNames <- c("PSet_gCSI2019.rds", "GDSC2.rds")
stems <- c("EGAN","EGAR","SRR")





for( idx in seq(3)) {
	
	ps.name <- PSet.names[idx]
	print(ps.name)
	out.path <- paste0(output.dir, ps.name,"/")
	ccl_stem <- stems[idx]
	
	ifelse(!dir.exists(out.path),dir.create(out.path),FALSE)
	
	if(ps.name %in% c("gCSI","GDSC2")){
		
		ps.file <- PSet.fileNames[idx]
		ps <- readRDS(paste0(PSet.baseDir,ps.file))
		ps <- updateObject(ps)
		exp.pset <- ps 
		resp.pset <- ps

	} else {
		
		resp.pset <- readRDS("../data/raw/PSet_CTRPv2.rds")
		exp.pset <- readRDS("../data/raw/CCLE.rds")
		resp.pset <-  updateObject(resp.pset)
		exp.pset <- updateObject(exp.pset)

	}


	drug.data <- as.data.frame(drugInfo(resp.pset))%>% select(treatmentid,smiles)
	rownames(drug.data) <- NULL
	
	ccl_info <- as.data.frame(sampleInfo(exp.pset))
	keep_tissues <- names(which(table(ccl_info$tissueid)>count_threshold))

	keep_sample_ids <- ccl_info %>% 
		filter(tissueid %in% keep_tissues & !tissueid %in% drop_tissues) %>% 
		select(sampleid) %>% 
		pull(sampleid)
	
	ccl_info <- ccl_info %>% 
		filter(sampleid %in% keep_sample_ids) %>%
		select(sampleid, tissueid)

	write.csv(ccl_info,paste0(out.path,"ccl_data.csv"),row.names=FALSE)

	write.csv(drug.data,paste0(out.path,"drug_data.csv"),row.names=FALSE)


	experiment_info <- as.data.frame(treatmentResponse(resp.pset)$info)
	experiment_info <- experiment_info %>% filter(sampleid %in% keep_sample_ids)

	response_vars <- treatmentResponse(resp.pset)$profiles %>% select(aac_recomputed) %>% drop_na()

	targets <- merge(experiment_info,response_vars, by=0) %>% 
		select(sampleid, treatmentid, aac_recomputed) %>% 
		group_by(sampleid,treatmentid) %>%
		summarize(aac = mean(aac_recomputed))

	colnames(targets)<-c("Cell_line",'Drug','AAC')
		
	write.csv(as.data.frame(targets),paste0(out.path,"responses.csv"),row.names=FALSE)
	
	expression <- molecularProfiles(exp.pset)$Kallisto_0.46.1.rnaseq

	gene_info <- as.data.frame(rowData(expression)) %>% 
		select(gene_name) 
	
	all_genes <- as.vector(gene_info$gene_name)%>% unique()
	

	protein_coding_genes <- as.data.frame(rowData(expression)) %>%
		filter(gene_type == 'protein_coding') %>%
		select(gene_name)
	
	pc_genes <- as.vector(protein_coding_genes$gene_name) %>% unique()

	
	write(pc_genes,paste0(out.path,"pc_genes.txt"))
	write(all_genes,paste0(out.path,"all_genes.txt"))
	expression.data <- merge(protein_coding_genes, as.data.frame(assay(expression)), by=0)%>%
		select(-Row.names) %>%
		group_by(gene_name) %>% 
		summarise_at(vars(matches(ccl_stem)),median)

	
	keep.ccls <- as.data.frame(colData(expression)) %>%	
			filter(sampleid %in% keep_sample_ids) %>% 
			select(sampleid)
	
	colnames(keep.ccls)<-c("Cell_line")
	
	expression.data <- expression.data %>% 
		as.data.frame() %>% 
		tibble::column_to_rownames("gene_name") %>% t()


	expression.data <- keep.ccls %>% 
		select(Cell_line) %>% 
		merge(expression.data, by=0) %>% 
		select(-Row.names) %>% 
		group_by(Cell_line) %>% 
		summarize_all(median) %>%
		as.data.frame()

	write.csv(expression.data,paste0(out.path,"expression.csv"),row.names=FALSE)
}

broad_pc <- read.delim(paste0(output.dir,"Broad","/","pc_genes.txt"),header = FALSE)
broad_genes <- read.delim(paste0(output.dir,"Broad","/","all_genes.txt"),header = FALSE)
gdsc_pc <- read.delim(paste0(output.dir,"GDSC2","/","pc_genes.txt"),header = FALSE)
gdsc_genes <- read.delim(paste0(output.dir,"GDSC2","/","all_genes.txt"),header= FALSE)
gCSI_pc <- read.delim(paste0(output.dir,"gCSI","/","pc_genes.txt"), header = FALSE)
gCSI_genes <- read.delim(paste0(output.dir,"gCSI","/","all_genes.txt"), header = FALSE)

all_pcs <- as.vector(intersect(broad_pc$V1,intersect(gdsc_pc$V1,gCSI_pc$V1)))
all_genes <-  as.vector(intersect(broad_genes$V1,intersect(gdsc_genes$V1,gCSI_genes$V1)))

write(all_pcs, "../data/all_pcs.txt")
write(all_genes,"../data/all_genes.txt")
all_pcs <- gsub("-",".",all_pcs)

for(ds.name in c("gCSI","Broad","GDSC2")){
	print(ds.name)
	expr <- read.csv(paste0(output.dir,ds.name,"/expression.csv"))
	expr <- expr %>% select(all_of(c('Cell_line',all_pcs)))
	write.csv(expr,paste0(output.dir,ds.name,"/expression.csv"),row.names=FALSE)
}