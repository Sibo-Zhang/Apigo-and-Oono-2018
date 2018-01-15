# Appendix C. Apigo and Oono (2018) Host Specificity in Foliar Fungal Endophytes. 
# Email: austen.apigo@lifesci.ucsb.edu

#############################
##### Table of Contents ##### 
#############################

# Appendix C1: Multivariate Cutoff Level Analysis
# Appendix C2: Calculate Structural, Network and Phylogenetic Specificity
# Appendix C3: Comparisons to a Random Community Assembly Model

# The code below functions to reproduce analyses from  Apigo and Oono (2018) 
# Host Specificity in Foliar Fungal Endophytes. You may need to amend the code 
# for your purposes.

# There are two separate OTU tables for (1) structural and network specificity & 
# (2) phylogenetic specificity because not all plant host species had available matk 
# and rbcl sequence data in GenBank (see methodology). 

#######################################################################
##### APPENDIX C1: Multivariate Cutoff Level Analysis (MultiCoLA) #####
#######################################################################
  
  # Set your working directory to the folder downloaded from Github.
  setwd("~/Desktop/Apigo-and-Oono-2018-master")  
  
  # Read in Appendix B1 for structural and network specificity and B2 for phylogenetic specificity.
  A.B1 <- read.csv("Appendix_B1_Structural_Network_OTUtable.csv", header = TRUE, row.names = 1) 
  A.B2 <- read.csv("Appendix_B2_Phylogenetic_OTUtable.csv", header = TRUE, row.names = 1) 

  # Add a new column with taxonomic information to each OTU table necessary
  # for the MultiCoLA (Gobet et al. 2010) code to run. The input file can be 
  # further modified to assess the impact within various taxonomic levels 
  # (e.g., Phylum, Class, Order, etc.), if desired. 
  A.B1$Kingdom <- "Fungi"
  A.B2$Kingdom <- "Fungi"

  # Move into MultiCoLA folder and load the taxa.pooler function. 
  setwd("~/Desktop/Apigo-and-Oono-2018-master/MultiCoLA_Gobet_et_al_2010")
  source("taxa.pooler.1.4.r")

  # Run the all_taxa_pooled function.
  all_taxa_pooled.A.B1 <- taxa.pooler(A.B1)
  all_taxa_pooled.A.B2 <- taxa.pooler(A.B2)
  # You will be prompted to answer questions --  
  # Number of samples? (e.g. 16)... 
    # List the number of samples [# of columns - 1 taxonomic column, in this scenario]. 
    # There are 37 samples for Appendix B1 and 30 samples for Appendix B2. 
  # Number of taxonomic levels? (e.g. phylum+class+order+family+genus=5)...
    # In this case, there is one taxonomic level.
  # Presence/absence tables as output? (y/n)...
    # This is optional and not needed for this analysis.
  # Output as text files? (y/n)...
    # This is optional and not needed for this analysis. 
  
  # Load the COtables function to truncate the OTU tables. 
  source("COtables.1.4.r")
  # The "dominant" option truncates the rare OTUs,  
  # whereas the "rare" option truncates the most abundant OTUs. 
  truncated.B1 <- COtables(all_taxa_pooled.A.B1[[2]], Type="ADS",typem="dominant")
  truncated.B2 <- COtables(all_taxa_pooled.A.B2[[2]], Type="ADS",typem="dominant")
  
  # Move out of MultiCoLA folder.  
  setwd("~/Desktop/Apigo-and-Oono-2018-master")
  # Make separate data frames for each OTU table.
  # '0.99' denotes a dataset with only 1% of sequences removed. 
  # Each data frame is transposed so they are in the same layout as original OTU table. 
  df.B1.100 <- read.csv("Appendix_B1_Structural_Network_OTUtable.csv", header = TRUE, row.names = 1)
  df.B1.99 <- t(as.data.frame(truncated.B1$`0.99`))
  df.B1.95 <- t(as.data.frame(truncated.B1$`0.95`))
  df.B1.90 <- t(as.data.frame(truncated.B1$`0.9`))
  df.B1.85 <- t(as.data.frame(truncated.B1$`0.85`))
  df.B1.80 <- t(as.data.frame(truncated.B1$`0.8`))
  # Make a list of the data frames
  df.B1.list <- list(df.B1.100, df.B1.99, df.B1.95, 
                          df.B1.90, df.B1.85, df.B1.80)
  # Repeat the code for Appendix B2.
  df.B2.100 <- read.csv("Appendix_B2_Phylogenetic_OTUtable.csv", header = TRUE, row.names = 1)
  df.B2.99 <- t(as.data.frame(truncated.B2$`0.99`))
  df.B2.95 <- t(as.data.frame(truncated.B2$`0.95`))
  df.B2.90 <- t(as.data.frame(truncated.B2$`0.9`))
  df.B2.85 <- t(as.data.frame(truncated.B2$`0.85`))
  df.B2.80 <- t(as.data.frame(truncated.B2$`0.8`))
  # Make a list of the data frames.
  df.B2.list <- list(df.B2.100, df.B2.99, df.B2.95, 
                         df.B2.90, df.B2.85, df.B2.80)
  # Citation
  # Gobet A, Quince C, Ramette A (2010) Multivariate Cutoff Level Analysis (MultiCoLA) 
  # of large community data sets. Nucleic Acids Res 38: e155
  
###################################################################################
##### APPENDIX C2: Calculate Structural, Network and Phylogenetic Specificity #####
###################################################################################
  
  ##########################
  # Structural Specificity #
  ##########################
  
  # The OTUs must populate the rows and samples must populate the columns. 
  library(vegan)
  # Calculate host richness across the list of data frames.  
  host.richness <- lapply(df.B1.list, specnumber)
  # Calculate Shannon's H across the list of data frames.
  shannons.h <- lapply(df.B1.list, diversity)
  # Citation for vegan package
  citation("vegan")
  # Construct final data frames for host richness and Shannon's H
  library(plyr)
  library(dplyr)
  library(data.table)
  # Read in the metadata file. 
  setwd("~/Desktop/Apigo-and-Oono-2018-master")
  meta <- read.csv("Appendix_C2_Metadata.csv", header = T)
  # Add a new column to each OTU table in the list with the percent of sequences removed. 
  pct <- list(0, 1, 5, 10, 15, 20)
  # Add the columns to each OTU table in the data frame list.
  percent.removed.hr <- Map(cbind, host.richness, Percent.Removed = pct) 
  percent.removed.sh <- Map(cbind, shannons.h, Percent.Removed = pct) 
  # For structural and phylogenetic specificity ONLY, negate the values and add them
  # to a new column so more positive values correspond with higher host specificity.
  neg.columns.hr <- list(percent.removed.hr[[1]][,1]*-1, percent.removed.hr[[2]][,1]*-1, 
                      percent.removed.hr[[3]][,1]*-1, percent.removed.hr[[4]][,1]*-1,
                      percent.removed.hr[[5]][,1]*-1, percent.removed.hr[[6]][,1]*-1
  )
  neg.columns.sh <- list(percent.removed.sh[[1]][,1]*-1, percent.removed.sh[[2]][,1]*-1, 
                      percent.removed.sh[[3]][,1]*-1, percent.removed.sh[[4]][,1]*-1,
                      percent.removed.sh[[5]][,1]*-1, percent.removed.sh[[6]][,1]*-1
  )
  # Add these columns to each OTU table in the list. 
  negated.hr <- Map(cbind, percent.removed.hr, Host.Richness = neg.columns.hr)
  negated.sh <- Map(cbind, percent.removed.sh, Shannons.H = neg.columns.sh)
  # Make new columns with the OTUID for each OTU table to match with the metadata. 
  rows.hr <- list(row.names(percent.removed.hr[[1]]), row.names(percent.removed.hr[[2]]), 
               row.names(percent.removed.hr[[3]]), row.names(percent.removed.hr[[4]]),
               row.names(percent.removed.hr[[5]]), row.names(percent.removed.hr[[6]])
  )
  rows.sh <- list(row.names(percent.removed.sh[[1]]), row.names(percent.removed.sh[[2]]), 
               row.names(percent.removed.sh[[3]]), row.names(percent.removed.sh[[4]]),
               row.names(percent.removed.sh[[5]]), row.names(percent.removed.sh[[6]])
  )
  # Add the new column to each OTU table in the list. 
  otu.row.hr <- Map(cbind, negated.hr, OTUs = rows.hr) 
  otu.row.sh <- Map(cbind, negated.sh, OTUs = rows.sh) 
  # Change the order of the columns. 
  otu.row.change.order.hr <- lapply(otu.row.hr, function(x) {x[, c(4, 2, 1, 3)]})
  otu.row.change.order.sh <- lapply(otu.row.sh, function(x) {x[, c(4, 2, 1, 3)]})
  # Join the metadata and specificity data within each file. 
  data.list.hr <- lapply(otu.row.change.order.hr, function(x) {join(as.data.frame(x), meta, by = "OTUs")})
  data.list.sh <- lapply(otu.row.change.order.sh, function(x) {join(as.data.frame(x), meta, by = "OTUs")})
  # Bind the list together into one dataframe (ignore warnings).
  host_richness_empirical <- bind_rows(data.list.hr) 
  shannons_h_empirical <- bind_rows(data.list.sh) 
  # Exclude the column with the positive values (only needed for structural and phylogenetic specificity).
  host_richness_empirical <- subset(host_richness_empirical, select = c(1,2,4:6))
  shannons_h_empirical <- subset(shannons_h_empirical, select = c(1,2,4:6))
  # Convert numeric colmuns to numeric. 
  host_richness_empirical$Percent.Removed <- as.numeric(as.character(host_richness_empirical$Percent.Removed))
  host_richness_empirical$Host.Richness <- as.numeric(as.character(host_richness_empirical$Host.Richness))
  shannons_h_empirical$Percent.Removed <- as.numeric(as.character(shannons_h_empirical$Percent.Removed))
  shannons_h_empirical$Shannons.H <- as.numeric(as.character(shannons_h_empirical$Shannons.H))
  # View the final datasets.
  View(host_richness_empirical)
  View(shannons_h_empirical)

  #######################
  # Network Specificity #
  #######################
  
  library(bipartite)
  # The OTUs must populate the columns and the samples must populate the rows. 
  # Transpose the list of data frames. 
  df.B1.list.transposed <- lapply(df.B1.list, function(x) {t(x)})
  # Make tables Presence-absence  
  df.B1.list.pa <- lapply(df.B1.list.transposed, function(x) {as.data.frame((x > 0) + 0)})
  # Calculate the Resource Range Index across the list of data frames.
  Resource.Range.Index <- lapply(df.B1.list.pa, PDI)
  # Calculate the Paired Difference Index across the list of data frames. 
  Paired.Difference.Index <- lapply(df.B1.list.transposed, PDI)
  # Citation for bipartite package
  citation("bipartite")
  # Construct final data frames for host richness and Shannon's H
  library(plyr)
  library(dplyr)
  library(data.table)
  # Read in the metadata file. 
  setwd("~/Desktop/Apigo-and-Oono-2018-master")
  meta <- read.csv("Appendix_C2_Metadata.csv", header = T)
  # Add a new column to each OTU table in the list with the percent of sequences removed. 
  pct <- list(0, 1, 5, 10, 15, 20)
  # Add columns to each OTU table in the data frame list 
  percent.removed.rri <- Map(cbind, Resource.Range.Index, Percent.Removed = pct) 
  percent.removed.pdi <- Map(cbind, Paired.Difference.Index, Percent.Removed = pct) 
  # Make new columns with the OTUID for each OTU table to match with the metadata. 
  rows.rri <- list(row.names(percent.removed.rri[[1]]), row.names(percent.removed.rri[[2]]), 
                  row.names(percent.removed.rri[[3]]), row.names(percent.removed.rri[[4]]),
                  row.names(percent.removed.rri[[5]]), row.names(percent.removed.rri[[6]])
  )
  rows.pdi <- list(row.names(percent.removed.pdi[[1]]), row.names(percent.removed.pdi[[2]]), 
                  row.names(percent.removed.pdi[[3]]), row.names(percent.removed.pdi[[4]]),
                  row.names(percent.removed.pdi[[5]]), row.names(percent.removed.pdi[[6]])
  )
  # Add the new column to each OTU table in the list. 
  otu.row.rri <- Map(cbind, percent.removed.rri, OTUs = rows.rri) 
  otu.row.pdi <- Map(cbind, percent.removed.pdi, OTUs = rows.pdi) 
  # Change the order and rename columns.
  otu.row.change.order.rri <- lapply(otu.row.rri, function(x) {x[, c(3, 2, 1)]})
  otu.row.change.order.pdi <- lapply(otu.row.pdi, function(x) {x[, c(3, 2, 1)]})
  # Join the metadata and specificity data within each file. 
  data.list.rri <- lapply(otu.row.change.order.rri, function(x) {join(as.data.frame(x), meta, by = "OTUs")})
  data.list.pdi <- lapply(otu.row.change.order.pdi, function(x) {join(as.data.frame(x), meta, by = "OTUs")})
  # Bind the list together into one dataframe (ignore warnings). 
  resource_range_empirical <- bind_rows(data.list.rri)
  paired_difference_empirical <- bind_rows(data.list.pdi) 
  # Convert numeric colmuns to numeric. 
  resource_range_empirical$Percent.Removed <- as.numeric(as.character(resource_range_empirical$Percent.Removed))
  resource_range_empirical$V3 <- as.numeric(as.character(resource_range_empirical$V3))
  paired_difference_empirical$Percent.Removed <- as.numeric(as.character(paired_difference_empirical$Percent.Removed))
  paired_difference_empirical$V3 <- as.numeric(as.character(paired_difference_empirical$V3))
  # Change column names to a more informative identifier.
  colnames(resource_range_empirical)[3] <- "Resource.Range"
  colnames(paired_difference_empirical)[3] <- "Paired.Difference" 
  # View the full host richness dataset
  View(resource_range_empirical)
  View(paired_difference_empirical)
  
  ############################
  # Phylogenetic Specificity #
  ############################
  
  # The OTUs must populate the rows and samples must populate the columns.  
  library(picante)
  # Read in ultrametric tree
  utree <- read.tree("Appendix_D_ultrametrictree.txt")
  
  # Calculate presence-absence mean pairwise phylogenetic distance (MPD)
  # across the list of data frames. I use the ses.mpd() function because it
  # prints the OTU IDs with the output, whereas the mpd() function does not. 
  # This code can take awhile. 
  mpd.pa <- lapply(df.B2.list, function(x) {ses.mpd(x, cophenetic(utree), 
                  null.model = "taxa.labels", runs = 1000, 
                  iterations = 1000, abundance.weighted = FALSE)}
  )
  # Select first two columns to keep 
  mpd.pa.subset <- lapply(mpd.pa, function(x) {subset(x, select = c(2))})
  # Replace NAs as zeros.
  mpd.pa.subset <- rapply(mpd.pa.subset, function(x) ifelse(is.na(x),0,x), how="replace")
  # Calculate abundance-weighted MPD across the list of data frames. 
  mpd.ab <- lapply(df.B2.list, function(x) {ses.mpd(x, cophenetic(utree), 
                                            null.model = "frequency", runs = 1000, 
                                            iterations = 1000, abundance.weighted = TRUE)}
  )
  # Select first two columns to keep 
  mpd.ab.subset <- lapply(mpd.ab, function(x) {subset(x, select = c(2))})
  # Replace NAs as zeros.
  mpd.ab.subset <- rapply(mpd.ab.subset, function(x) ifelse(is.na(x),0,x), how="replace")
  # Citation for picante package
  citation("picante")
  # Construct final data frames for presence-absence and abundance-weighted MPD
  library(plyr)
  library(dplyr)
  library(data.table)
  # Read in the metadata file. 
  setwd("~/Desktop/Apigo-and-Oono-2018-master")
  meta <- read.csv("Appendix_C2_Metadata.csv", header = T)
  # Add a new column to each OTU table in the list with the percent of sequences removed. 
  pct <- list(0, 1, 5, 10, 15, 20)
  # Add columns to each OTU table in the data frame list 
  percent.removedmpd.pa <- Map(cbind, mpd.pa.subset, Percent.Removed = pct) 
  percent.removedmpd.ab <- Map(cbind, mpd.ab.subset, Percent.Removed = pct) 
  # For structural and phylogenetic specificity ONLY, negate the values and add them
  # to a new column so more positive values correspond with higher host specificity.
  neg.columnsmpd.pa <- list(percent.removedmpd.pa[[1]][,1]*-1, percent.removedmpd.pa[[2]][,1]*-1, 
                            percent.removedmpd.pa[[3]][,1]*-1, percent.removedmpd.pa[[4]][,1]*-1,
                            percent.removedmpd.pa[[5]][,1]*-1, percent.removedmpd.pa[[6]][,1]*-1
  )
  neg.columnsmpd.ab <- list(percent.removedmpd.ab[[1]][,1]*-1, percent.removedmpd.ab[[2]][,1]*-1, 
                            percent.removedmpd.ab[[3]][,1]*-1, percent.removedmpd.ab[[4]][,1]*-1,
                            percent.removedmpd.ab[[5]][,1]*-1, percent.removedmpd.ab[[6]][,1]*-1
  )
  # Add these columns to each OTU table in the list. 
  negatedmpd.pa <- Map(cbind, percent.removedmpd.pa, MPD.Pres.Abs = neg.columnsmpd.pa)
  negatedmpd.ab <- Map(cbind, percent.removedmpd.ab, MPD.Abun.Weigh = neg.columnsmpd.ab)
  # Make new columns with the OTUID for each OTU table to match with the metadata. 
  rowsmpd.pa <- list(row.names(percent.removedmpd.pa[[1]]), row.names(percent.removedmpd.pa[[2]]), 
                     row.names(percent.removedmpd.pa[[3]]), row.names(percent.removedmpd.pa[[4]]),
                     row.names(percent.removedmpd.pa[[5]]), row.names(percent.removedmpd.pa[[6]])
  )
  rowsmpd.ab <- list(row.names(percent.removedmpd.ab[[1]]), row.names(percent.removedmpd.ab[[2]]), 
                     row.names(percent.removedmpd.ab[[3]]), row.names(percent.removedmpd.ab[[4]]),
                     row.names(percent.removedmpd.ab[[5]]), row.names(percent.removedmpd.ab[[6]])
  )
  # Add the new column to each OTU table in the list. 
  otu.rowmpd.pa <- Map(cbind, negatedmpd.pa, OTUs = rowsmpd.pa) 
  otu.rowmpd.ab <- Map(cbind, negatedmpd.ab, OTUs = rowsmpd.ab) 
  # Change the order of the columns. 
  otu.row.change.ordermpd.pa <- lapply(otu.rowmpd.pa, function(x) {x[, c(4, 2, 1, 3)]})
  otu.row.change.ordermpd.ab <- lapply(otu.rowmpd.ab, function(x) {x[, c(4, 2, 1, 3)]})
  # Join the metadata and specificity data within each file. 
  data.listmpd.pa <- lapply(otu.row.change.ordermpd.pa, function(x) {join(as.data.frame(x), meta, by = "OTUs")})
  data.listmpd.ab <- lapply(otu.row.change.ordermpd.ab, function(x) {join(as.data.frame(x), meta, by = "OTUs")})
  # Bind the list together into one dataframe. 
  mpd_pa_empirical <- bind_rows(data.listmpd.pa) 
  mpd_ab_empirical <- bind_rows(data.listmpd.ab) 
  # Exclude the column with the positive values (only needed for structural and phylogenetic specificity).
  mpd_pa_empirical <- subset(mpd_pa_empirical, select = c(1,2,4:6))
  mpd_ab_empirical <- subset(mpd_ab_empirical, select = c(1,2,4:6))
  # Convert numeric colmuns to numeric. 
  mpd_pa_empirical$Percent.Removed <- as.numeric(as.character(mpd_pa_empirical$Percent.Removed))
  mpd_pa_empirical$MPD.Pres.Abs <- as.numeric(as.character(mpd_pa_empirical$MPD.Pres.Abs))
  mpd_ab_empirical$Percent.Removed <- as.numeric(as.character(mpd_ab_empirical$Percent.Removed))
  mpd_ab_empirical$MPD.Abun.Weigh <- as.numeric(as.character(mpd_ab_empirical$MPD.Abun.Weigh))
  # View the full host richness dataset
  View(mpd_pa_empirical)
  View(mpd_ab_empirical)

#########################################################################
##### APPENDIX C3: Comparisons to a Random Community Assembly Model #####
#########################################################################
  
  # Make presence-absence randomized tables
  df.B1.list.pa <- lapply(df.B1.list, function(x) {as.data.frame((x > 0) + 0)})
  # Make transposed table for network specificity 
  df.B1.list.transposed <- lapply(df.B1.list, function(x) {t(x)})
  df.B1.list.transposed.pa <- lapply(df.B1.list.transposed, function(x) {as.data.frame((x > 0) + 0)})
  # Randomize presence-absence tables
  library(bipartite)
  null.otu.tables.B1.pa <- lapply(df.B1.list.pa, function(x) {shuffle.web(x, N=1000, legacy = FALSE)})
  # Create transposed null tables for network specificity
  null.otu.tables.B1.transposed.pa <- lapply(df.B1.list.transposed.pa, function(x) {shuffle.web(x, N=1000, legacy = FALSE)})
  # Randomize abundance-weighted tables
  null.otu.tables.B1.ab <- lapply(df.B1.list, function(x) {nullmodel(x, N=1000, method = 4)})
  # Create transposed null tables for network specificity
  null.otu.tables.B1.transposed.ab <- lapply(df.B1.list.transposed, function(x) {nullmodel(x, N=1000, method = 4)})
  
  ############################
  # Randomized Host Richness #
  ############################
  null.host.richness <-  lapply(null.otu.tables.B1.ab, sapply, specnumber)
  # Average host specificity measurements across 1000 randomized tables per OTU
  null.otu.means.host.richness <- lapply(null.host.richness, rowMeans)
  # Convert into a list of data.frames
  null.otu.means.host.richness <- lapply(null.otu.means.host.richness, as.data.frame)
  # Concatentate the list
  null.otu.means.host.richness <- bind_rows(null.otu.means.host.richness)
  # Combine with empirical dataset
  null.otu.means.host.richness <- bind_cols(null.otu.means.host.richness, host_richness_empirical)
  # Negate null averages
  null.otu.means.host.richness$Null.Host.Richness <- null.otu.means.host.richness$`X[[i]]` * -1
  # Subset dataset
  host_richness_full <- subset(null.otu.means.host.richness, select = c(2:7))
  # Re-order columns
  host_richness_full <- host_richness_full[, c(1, 2, 3, 6, 4, 5)]
  # Convert numeric columns to numeric
  host_richness_full$Percent.Removed <- as.numeric(as.character(host_richness_full$Percent.Removed))
  host_richness_full$Host.Richness <- as.numeric(as.character(host_richness_full$Host.Richness))
  host_richness_full$Null.Host.Richness <- as.numeric(as.character(host_richness_full$Null.Host.Richness))
  # View full dataset.
  View(host_richness_full)
  
  ##########################
  # Randomized Shannon's H #
  ##########################
  null.shannons.h <- lapply(null.otu.tables.B1.ab, sapply, diversity)
  # Average host specificity measurements across 1000 randomized tables per OTU
  null.otu.means.shannons.h <- lapply(null.shannons.h, rowMeans)
  # Convert into a list of data.frames
  null.otu.means.shannons.h <- lapply(null.otu.means.shannons.h, as.data.frame)
  # Concatentate the list
  null.otu.means.shannons.h <- bind_rows(null.otu.means.shannons.h)
  # Combine with empirical dataset
  null.otu.means.shannons.h <- bind_cols(null.otu.means.shannons.h, shannons_h_empirical)
  # Negate null averages
  null.otu.means.shannons.h$Null.Shannons.H <- null.otu.means.shannons.h$`X[[i]]` * -1
  # Subset dataset
  shannons_h_full <- subset(null.otu.means.shannons.h, select = c(2:7))
  # Re-order columns
  shannons_h_full <- shannons_h_full[, c(1, 2, 3, 6, 4, 5)]
  # Convert numeric columns to numeric
  shannons_h_full$Percent.Removed <- as.numeric(as.character(shannons_h_full$Percent.Removed))
  shannons_h_full$Shannons.H <- as.numeric(as.character(shannons_h_full$Shannons.H))
  shannons_h_full$Null.Shannons.H <- as.numeric(as.character(shannons_h_full$Null.Shannons.H))
  View(shannons_h_full)

  ###################################
  # Randomized Resource Range Index #
  ###################################
  null.resource.range <- lapply(null.otu.tables.B1.transposed.pa, sapply, PDI)
  # Average host specificity measurements across 1000 randomized tables per OTU
  null.otu.means.resource.range <- lapply(null.resource.range, rowMeans)
  # Convert into a list of data.frames
  null.otu.means.resource.range <- lapply(null.otu.means.resource.range, as.data.frame)
  # Concatentate the list
  null.otu.means.resource.range <- bind_rows(null.otu.means.resource.range)
  # Combine with empirical dataset
  resource_range_full <- bind_cols(null.otu.means.resource.range, resource_range_empirical)
  # Re-order columns
  resource_range_full <- resource_range_full[, c(2, 3, 4, 1, 5, 6)]
  # Convert numeric columns to numeric
  resource_range_full$Percent.Removed <- as.numeric(as.character(resource_range_full$Percent.Removed))
  resource_range_full$Resource.Range <- as.numeric(as.character(resource_range_full$Resource.Range))
  resource_range_full$'X[[i]]' <- as.numeric(as.character(resource_range_full$'X[[i]]'))
  # Change column names to a more informative identifier
  colnames(resource_range_full)[4] <- "Null.Resource.Range"
  # View full dataset.
  View(resource_range_full)
  
  ######################################
  # Randomized Paired Difference Index #
  ######################################
  null.paired.difference <- lapply(null.otu.tables.B1.transposed.ab, sapply, PDI)
  # Average host specificity measurements across 1000 randomized tables per OTU
  null.otu.means.paired.difference <- lapply(null.paired.difference, rowMeans)
  # Convert into a list of data.frames
  null.otu.means.paired.difference <- lapply(null.otu.means.paired.difference, as.data.frame)
  # Concatentate the list
  null.otu.means.paired.difference <- bind_rows(null.otu.means.paired.difference)
  # Combine with empirical dataset
  paired_difference_full <- bind_cols(null.otu.means.paired.difference, paired_difference_empirical)
  # Re-order columns
  paired_difference_full <- paired_difference_full[, c(2, 3, 4, 1, 5, 6)]
  # Convert numeric columns to numeric
  paired_difference_full$Percent.Removed <- as.numeric(as.character(paired_difference_full$Percent.Removed))
  paired_difference_full$Paired.Difference <- as.numeric(as.character(paired_difference_full$Paired.Difference))
  paired_difference_full$'X[[i]]' <- as.numeric(as.character(paired_difference_full$'X[[i]]'))
  # Change column names to a more informative identifier
  colnames(paired_difference_full)[4] <- "Null.Paired.Difference"
  # View full dataset.
  View(paired_difference_full)

  ###################################################################
  # Randomized Presence-Absence Mean Pairwise Phylogenetic Distance #
  ###################################################################
  # Rename the ses.mpd() output.
  null.mpd.pa <- mpd.pa
  # Average host specificity measurements across 1000 randomized tables per OTU
  null.otu.means.mpd.pa <- lapply(null.mpd.pa, rowMeans)
  # Convert into a list of data.frames
  null.otu.means.mpd.pa <- lapply(null.otu.means.mpd.pa, as.data.frame)
  # Concatentate the list by row
  null.otu.means.mpd.pa <- bind_rows(null.otu.means.mpd.pa)
  # Combine with empirical dataset
  null.otu.means.mpd.pa <- bind_cols(null.otu.means.mpd.pa, mpd_pa_empirical)
  is.nan.data.frame <- function(x)
  do.call(cbind, lapply(null.otu.means.mpd.pa, is.nan))
  null.otu.means.mpd.pa[is.nan(null.otu.means.mpd.pa)] <- 0
  # Negate null averages
  null.otu.means.mpd.pa$Null.MPD.Pres.Ab <- null.otu.means.mpd.pa$`X[[i]]` * -1
  # Subset dataset
  mpd_pa_full <- subset(null.otu.means.mpd.pa, select = c(2:7))
  # Re-order columns
  mpd_pa_full <- mpd_pa_full[, c(1, 2, 3, 6, 4, 5)]
  # Convert numeric columns to numeric
  mpd_pa_full$Percent.Removed <- as.numeric(as.character(mpd_pa_full$Percent.Removed))
  mpd_pa_full$MPD.Pres.Abs <- as.numeric(as.character(mpd_pa_full$MPD.Pres.Abs))
  mpd_pa_full$Null.MPD.Pres.Ab <- as.numeric(as.character(mpd_pa_full$Null.MPD.Pres.Ab))
  # View full dataset.
  View(mpd_pa_full)
  
  #####################################################################
  # Randomized Abundance-weighted Mean Pairwise Phylogenetic Distance #
  #####################################################################
  # Rename ses.mpd() output. 
  null.mpd.ab <- mpd.ab
  # Average host specificity measurements across 1000 randomized tables per OTU
  null.otu.means.mpd.ab <- lapply(null.mpd.ab, rowMeans)
  
  # Convert into a list of data.frames
  null.otu.means.mpd.ab <- lapply(null.otu.means.mpd.ab, as.data.frame)
  # Concatentate the list by row
  null.otu.means.mpd.ab <- bind_rows(null.otu.means.mpd.ab)
  # Combine with empirical dataset
  null.otu.means.mpd.ab <- bind_cols(null.otu.means.mpd.ab, mpd_ab_empirical)
  is.nan.data.frame <- function(x)
    do.call(cbind, lapply(null.otu.means.mpd.ab, is.na))
  null.otu.means.mpd.ab[is.nan(null.otu.means.mpd.ab)] <- 0
  # Negate null averages
  null.otu.means.mpd.ab$Null.MPD.Abun.Weigh <- null.otu.means.mpd.ab$`X[[i]]` * -1
  # Subset dataset
  mpd_ab_full <- subset(null.otu.means.mpd.ab, select = c(2:7))
  # Re-order columns
  mpd_ab_full <- mpd_ab_full[, c(1, 2, 3, 6, 4, 5)]
  # Convert numeric columns to numeric
  mpd_ab_full$Percent.Removed <- as.numeric(as.character(mpd_ab_full$Percent.Removed))
  mpd_ab_full$MPD.Abun.Weigh <- as.numeric(as.character(mpd_ab_full$MPD.Abun.Weigh))
  mpd_ab_full$Null.MPD.Abun.Weigh <- as.numeric(as.character(mpd_ab_full$Null.MPD.Abun.Weigh))
  # View full dataset.
  View(mpd_ab_full)
 