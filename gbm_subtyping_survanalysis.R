rm(list=ls())
# # Load the library
library("DCEM")
library("CancerSubtypes")

# Load GBM data
GBM_mRNA <- readRDS(file="GBM_mRNA_data.rds")
GBM_DNA <- readRDS(file="GBM_DNA_data.rds")
GBM_miRNA <- readRDS(file="GBM_miRNA_data.rds")
GBM_clinical <- readRDS(file="GBM_clinical_data.rds")

GBM_mRNA_FS=FSbyCox(GBM_mRNA,GBM_clinical$time,GBM_clinical$status,cutoff=0.00001)
GBM_DNA_FS=FSbyCox(GBM_DNA,GBM_clinical$time,GBM_clinical$status,cutoff=0.00001)
GBM_miRNA_FS=FSbyCox(GBM_miRNA,GBM_clinical$time,GBM_clinical$status,cutoff=0.00001)

hstackedGBM <- rbind(GBM_DNA_FS, GBM_mRNA_FS, GBM_miRNA_FS)
hstackedGBM=scale(hstackedGBM,center=T,scale=T)



start_time <- Sys.time()
# Train model
dcem_out = dcem_star_train(data = t(hstackedGBM), iteration_count = 100, num_clusters = 4, seeding = 'improved')
end_time <- Sys.time()

#dcem_out$prob         # estimated posterior probabilities
#dcem_out$meu          # estimated mean of the clusters
#dcem_out$sigma        # estimated covariance matrices
#dcem_out$priors       # estimated priors
#dcem_out$membership
clusters <- dcem_out$membership  # membership of data points based on maximum liklihood (posterior probabilities)

group = integer()

for (x in 1:276) {
  group[x] <- clusters[[x]]
} 


p_value=survAnalysis(mainTitle="Survival Analysis", GBM_clinical$time, GBM_clinical$status, group, 
                     distanceMatrix=NULL, similarity=TRUE)
print(end_time - start_time)




