## Run this line by line
# check each output is as expected for the antismash run
# import list of clusters from antismash
# each cluster is currently in a different folder because I ran antismash using the online server

# working directory
dir <- "G:/Secmet_TSA_induction_round2"
dir

# antismash directory
anti_dir <- paste0(dir,"/Anidulans_genome_files/Antismash_nocassis/")

# read in list of genes in each cluster
# and also find backbone names
# devtools::install_github("nvelden/geneviewer")
library(geneviewer)

gb_filenames <- list.files(path=anti_dir,recursive=TRUE,pattern="region",full.names = TRUE)
# only keep .gbk files
gb_filenames <- gb_filenames[grep("\\.gbk", gb_filenames)]

# also get shorter filename list
all_clusters <- list.files(path=anti_dir,recursive=TRUE,pattern="region",full.names = FALSE)
all_clusters <- all_clusters[grep("\\.gbk", all_clusters)]

# write a function to extract genes from each cluster
get_cluster_genes <- function(gbk_file) {
  # only extract features we need
  gb <- read_gbk(gbk_file,features=c("gene",
                                     "aSModule",
                                     "PFAM_domain"))  
  # get feature details for each gene
  names(gb)
  # generalize it
  genelist <- gb[[1]]$FEATURES$gene
  
  # extract only the locus tags
  
  old_locus_tags <- lapply(genelist, function(df) {
    df[["old_locus_tag"]]  # Extract the locus_tag column
  })
  old_locus_tags <- unlist(old_locus_tags)
  # remove everything after the .
  
  old_locus_tags <- sub("\\..*", "", old_locus_tags)
  
  # also get new locus tags
  locus_tags <- lapply(genelist, function(df) {
    df[["locus_tag"]]  # Extract the locus_tag column
  })
  locus_tags <- unlist(locus_tags)
  
  # match them into a df
  genes <- as.data.frame(cbind(locus_tags,old_locus_tags))
  genes$type <- NA
  
  # extract backbone list from aSModule
  core_list <- gb[[1]]$FEATURES$aSModule
  backbones <- lapply(core_list, function(df) {
    df[["locus_tags"]]
  })
  
  backbones <- unlist(backbones)
  
  # indicate which genes are backbone genes in the gene list
  head(genes)
  genes$type <- ifelse(genes$locus_tags %in% backbones, "backbone", genes$type)
  print(gbk_file)
  return(genes)
}




# now run it on every cluster and append the results together

# initialize the first entry
clusters <- get_cluster_genes(gb_filenames[1])
clusters$Cluster <- all_clusters[1]

for (m in c(1:35,37:length(gb_filenames))) {
  temp <- get_cluster_genes(gb_filenames[m])
  
  # add cluster name to temp
  temp$Cluster <- all_clusters[m]
  
  # rbind to result dataframe
  clusters <- rbind(clusters,temp)
}

# these are the entries that did not work:
# 36,58

# neither of these are importing for some reason. manually extract values into csv, then import here:
last_clusters <- read.csv(paste0(anti_dir,"entry36_58_antismash_genes.csv"))

clusters <- rbind(clusters,last_clusters)


# remove duplicate rows
clusters <- distinct(clusters)

# check number of clusters
length(unique(clusters$Cluster))

# answer is 58 - looks good!

# save resulting gene list
# write.csv(clusters,paste0(anti_dir,"antismash_clusters_v7.1.0.csv"),row.names=FALSE)



# extract PFAM domains
get_pfam_list <- function(gbk_file) {
  gb <- read_gbk(gbk_file,
                 features=c("gene",
                            "aSModule",
                            "PFAM_domain"),
                 origin=FALSE)  
  # extract pfam info
  pfam_list <- gb[[1]]$FEATURES$PFAM_domain
  pfam_list2 <- unlist(lapply(pfam_list, function(df) {df[["db_xref"]]}))
  pfam_list3 <- unlist(lapply(pfam_list, function(df) {df[["locus_tag"]]}))
  pfam_list4 <- unlist(lapply(pfam_list, function(df) {df[["description"]]}))
  
  # paste them together into a df
  pfam <- as.data.frame(cbind(pfam_list3,pfam_list2,pfam_list4))
  return(pfam)
}


# run get_pfam_list on all clusters
cluster_pfam <- get_pfam_list(gb_filenames[1])
cluster_pfam$Cluster <- all_clusters[1]
colnames(cluster_pfam)[1:3] <- c("locus_tags",
                                 "pfam_db_xref",
                                 "pfam_description")


for (n in c(1:35,37:length(gb_filenames))) {
  temp <- get_pfam_list(gb_filenames[n])
  temp$Cluster <- all_clusters[n]
  colnames(temp)[1:3] <- c("locus_tags",
                                   "pfam_db_xref",
                                   "pfam_description")
  
  # rbind to result dataframe
  cluster_pfam <- rbind(cluster_pfam,temp)
}

# these are the entries that did not work:
# 36,58
# compile these manually
# then append to the list
last_clusters_pfam <- read.csv(paste0(anti_dir,"entry36_58_antismash_pfam.csv"))
cluster_pfam <- rbind(cluster_pfam,last_clusters_pfam)

# remove duplicate entries
# for proteins with more than one of the same domain there can be multiple entries
cluster_pfam <-  distinct(cluster_pfam)

# fix cluster name column
cluster_pfam$Cluster <- sub(".*/", "", cluster_pfam$Cluster)
cluster_pfam$Cluster <- sub("\\.gbk.*", "", cluster_pfam$Cluster)

# export pfam file
write.csv(cluster_pfam, paste0(anti_dir,"all_clusters_pfam_domains_antismashv7.1.0.csv"),row.names=FALSE)


# find json files to extract backbone cluster type info
json_files <- list.files(path=anti_dir,recursive=TRUE,pattern=".json",full.names = TRUE)
json_files2 <- list.files(path=anti_dir,recursive=TRUE,pattern=".json",full.names = FALSE)

# write a function to extract cluster type from a json file
get_cluster_type <- function(json_file) {
  
  library(rjson)
  
  # import json file
  chr_data <- fromJSON(paste(readLines(json_file), collapse=""))
  
  # get cluster details
  chr_data$records[[1]]$areas
  
  # extract only the "products"
  products <- lapply(chr_data$records[[1]]$areas, function(df) {
    df[["products"]]
  })
  
  # collapse with a comma
  products <- sapply(products, function(x) paste(x, collapse = ", "))
  names(products) <- "cluster_type"
  return(products)
}



# get cluster type for each
# start with 1
# test the function
cluster_types <- as.data.frame(get_cluster_type(json_files[1]))
cluster_types$Chromosome <- json_files2[1]
colnames(cluster_types)[1] <- "cluster_type"

for (i in 2:length(json_files)) {
  temp <- as.data.frame(get_cluster_type(json_files[i]))
  temp$Chromosome <- json_files2[i]
  colnames(temp)[1] <- "cluster_type"
  cluster_types <- rbind(cluster_types,temp)
}

# add region numbers to each
cluster_types$Cluster <- all_clusters

# add backbone gene names to each
cluster_bb <- subset(clusters, type == "backbone")

# Comma-separate locus_tags for each unique Cluster using dplyr
cluster_bb2 <- cluster_bb %>%
  group_by(Cluster) %>%
  summarize(old_locus_tags = paste(old_locus_tags, collapse = ", "),
            locus_tags = paste(locus_tags, collapse = ", "))


# find list of entries that are missing backbone label
missing_clusters <- setdiff(all_clusters, cluster_bb2$Cluster)

# Print the missing clusters
print(missing_clusters)

# manually specify these in a csv using excel and then import
write.csv(missing_clusters,"missing_cluster_bb.csv",row.names=FALSE)
# merge these onto backbone cluster dataframe
missing_clusters_bb <- read.csv(paste0(anti_dir,"missing_cluster_bb.csv"))
clusters_bb_all <- rbind(cluster_bb2,missing_clusters_bb)


# also export list of all clusters
# then manually annotate with known cluster similarity from antismash index.html
write.csv(all_clusters,"all_clusters.csv",row.names=FALSE)
known_cluster_sum <- read.csv(paste0(anti_dir,"all_clusters_known_similarity.csv"))

# merge together final cluster summary
merged_cluster_sum <- merge(clusters_bb_all,cluster_types,by="Cluster")
merged_cluster_sum <- merge(merged_cluster_sum,known_cluster_sum,by="Cluster")

# clean up the Cluster names
merged_cluster_sum$Cluster <- sub(".*/", "", merged_cluster_sum$Cluster)
merged_cluster_sum$Cluster <- sub("\\.gbk.*", "", merged_cluster_sum$Cluster)

# simplify chromosome names
merged_cluster_sum$Chromosome <- sub("\\/.*", "", merged_cluster_sum$Chromosome)

# get list of chromosomes
# unique(merged_cluster_sum$Chromosome)

# define what to replace them with
old_chr <- c("NC_066257",
  "NC_066258",
  "NC_066259",
  "NC_066260",
  "NC_066261",
  "NC_066262",
  "NC_066263",
  "NC_066264")
new_chr <- c("Chr_I",
             "Chr_II",
             "Chr_III",
             "Chr_IV",
             "Chr_V",
             "Chr_VI",
             "Chr_VII",
             "Chr_VIII")

# make sure all indexes are found in both lists
indices <- which(merged_cluster_sum$Chromosome %in% old_chr)

# replace old NCBI ID values with clarified chromosome names
merged_cluster_sum$Chromosome[indices] <- new_chr[match(merged_cluster_sum$Chromosome[indices], old_chr)]

# save cluster summary output
write.csv(merged_cluster_sum, paste0(anti_dir,"all_clusters_summarized_antismashv7.1.0.csv"),row.names=FALSE)

# finish annotating gene types in the clusters df

# clean up cluster column
clusters$Cluster <- sub(".*/", "", clusters$Cluster)
clusters$Cluster <- sub("\\.gbk.*", "", clusters$Cluster)

# update backbone labels
# first separate entries with commas
all_bbs <- merged_cluster_sum$locus_tags
all_bbs_split <- unlist(strsplit(all_bbs, ", "))

# add label
clusters$type <- ifelse(clusters$locus_tags %in% all_bbs_split, "backbone", clusters$type)

# check summary
table(clusters$type)

# add transporter label
transporters <- read.csv(paste0(anti_dir,"../710transporters_anid_transportDB2.csv"))
clusters$type <- ifelse(clusters$old_locus_tags %in% transporters$old_locus_tags, "transporter", clusters$type)
table(clusters$type)

# add transcription factor label
tfs <- read.csv(paste0(anti_dir,"../Etxebeste_2021_TF_list.csv"))
clusters$type <- ifelse(clusters$old_locus_tags %in% tfs$old_locus_tags, "transcription factor", clusters$type)
table(clusters$type)

# if blank add NA
clusters$type[clusters$type == ""] <- NA
table(clusters$type)

# save output
write.csv(clusters,paste0(anti_dir,"all_clusters_genelist_antismashv7.1.0.csv"),row.names=FALSE)

# finally - add same labels to known cluster list
# this was manually curated based on Caesar et al. 2020, pmid 33035657
known_clusters <- read.csv(paste0(anti_dir,"../anid_known_clusters.csv"))

# add labels like above
table(known_clusters$type)

# first add additional backbone labels if necessary
all_bbs <- merged_cluster_sum$old_locus_tags
all_bbs_split <- unlist(strsplit(all_bbs, ", "))

known_clusters$type <- ifelse(known_clusters$GeneID %in% all_bbs_split, "backbone", known_clusters$type)
table(known_clusters$type)

# now add additional transporter labels
known_clusters$type <- ifelse(known_clusters$GeneID %in% transporters$old_locus_tags, "transporter", known_clusters$type)
table(known_clusters$type)

# add TF labels
known_clusters$type <- ifelse(known_clusters$GeneID %in% tfs$old_locus_tags, "transcription factor", known_clusters$type)
table(known_clusters$type)

# save as new file
write.csv(known_clusters,paste0(anti_dir,"../anid_known_clusters_manually_labeled.csv"))
