rm(list=ls())
# # Load the library
library("DCEM")
library("CancerSubtypes")

#Load BRCA data
BRCA_mRNA <- readRDS(file="BRCA_mRNA_data.rds")
BRCA_DNA <- readRDS(file="BRCA_DNA_data.rds")
BRCA_miRNA <- readRDS(file="BRCA_miRNA_data.rds")
BRCA_clinical <- readRDS(file="BRCA_clinical_data.rds")

BRCA_mRNA_FS=FSbyCox(BRCA_mRNA,BRCA_clinical$time,BRCA_clinical$status,cutoff=0.00001)
BRCA_DNA_FS=FSbyCox(BRCA_DNA,BRCA_clinical$time,BRCA_clinical$status,cutoff=0.00001)
BRCA_miRNA_FS=FSbyCox(BRCA_miRNA,BRCA_clinical$time,BRCA_clinical$status,cutoff=0.00001)

hstackedBRCA <- rbind(BRCA_DNA_FS, BRCA_mRNA_FS, BRCA_miRNA_FS)
hstackedBRCA=scale(hstackedBRCA,center=T,scale=T)

start_time <- Sys.time()
# Train model
dcem_out = dcem_star_train(data = t(hstackedBRCA), iteration_count = 100, num_clusters = 4, seeding = 'improved')
end_time <- Sys.time()

#dcem_out$prob         # estimated posterior probabilities
#dcem_out$meu          # estimated mean of the clusters
#dcem_out$sigma        # estimated covariance matrices
#dcem_out$priors       # estimated priors
#dcem_out$membership

clusters <- dcem_out$membership  # membership of data points based on maximum liklihood (posterior probabilities)

group = integer()

for (x in 1:302) {
  group[x] <- clusters[[x]]
} 


p_value=survAnalysis(mainTitle="Survival Analysis", BRCA_clinical$time, BRCA_clinical$status, group, 
                     distanceMatrix=NULL, similarity=TRUE)
print(end_time - start_time)




