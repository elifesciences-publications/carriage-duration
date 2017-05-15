require(dplyr)
require(msm)
require(foreach)
require(minqa)
require(data.table)

setwd("~/Documents/PhD/Maela carriage length/")

#*********************************#
# Functions                       #
#*********************************#

# Removes NTs from latex sweeps; combines 15B and 15C, and 6A and 6C from WHO data
combine_sero <- function(sero_list, latex=F)
{
  sero_list <- replace(sero_list, sero_list=="15B", "15B/C")
  sero_list <- replace(sero_list, sero_list=="15C", "15B/C")
  
  if (latex==T)
  {
    vec_ret <- sero_list[sero_list != "NT" & sero_list !=""]
  }
  else
  {
    sero_list <- replace(sero_list, sero_list=="6A", "6A/C")
    sero_list <- replace(sero_list, sero_list=="6C", "6A/C")
    vec_ret <- sero_list[sero_list !=""]
  }
  return(vec_ret)
}

# Combines WHO and latex serotypes
routine_parse_sero <- function(sero_list)
{
  sero_list <- unique(unlist(c(combine_sero(c(sero_list[1], sero_list[2], sero_list[3], sero_list[4]), latex=T),
               combine_sero(c(sero_list[5], sero_list[6], sero_list[7])))))
  sero_list <- sero_list[!is.na(sero_list)]
  
  return(sero_list)
}

# Create dummy variable for observed state. Various versions depending on HMM used
createDummy <- function(serotype, carried, carrying_serotypes)
{
  dummy <- 1
  
  carrying_serotypes <- unlist(strsplit(carrying_serotypes, split=","))
  if (serotype %in% carrying_serotypes)
  {
    dummy <- 2
  }
  else if (carried == 1)
  {
    dummy <-3
  }
  
  return(dummy)
}

createDummyTwoState <- function(serotype, carried, carrying_serotypes)
{
  dummy <- 1
  
  carrying_serotypes <- unlist(strsplit(carrying_serotypes, split=","))
  if (serotype %in% carrying_serotypes)
  {
    dummy <- 2
  }
  
  return(dummy)
}

createDummyThreeState <- function(serotype, carried, carrying_serotypes)
{
  dummy <- 1
  
  carrying_serotypes <- unlist(strsplit(carrying_serotypes, split=","))
  if (serotype %in% carrying_serotypes)
  {
    if (length(carrying_serotypes) == 1)
    {
      dummy <- 2
    }
    else
    {
      dummy <- 3
    }
  }
  
  return(dummy)
}

# iterate through path to get length
getDateBoundary <- function(obs_row, codenum, path, observations, direction=-1)
{
  boundary <- pathTraverse(obs_row, codenum, path, direction)
  if (!is.na(path$subject[boundary+direction]) & path$subject[boundary+direction] == codenum)
  {
    date_boundary <- observations$specdate[boundary+direction] + 
        (observations$specdate[boundary] - observations$specdate[boundary+direction])/2
  }
  else
  {
    date_boundary <- observations$specdate[boundary]
  }
  
  return(date_boundary)
}

# similar to above function, but gets whether this is the first carriage episode or a subsequent one
getCarried <- function(obs_row, codenum, path, observations)
{
  boundary <- pathTraverse(obs_row, codenum, path)
  
  if (path$subject[boundary-1] == codenum)
  {
    # As will always be carrying at the time of a positive swab
    carried <- observations$carried[boundary-1]
  }
  else
  {
    # The first episode of carriage
    carried <- 0
  }
  
  return(carried)
}

# Path traversing utility for above two functions
pathTraverse <- function(obs_row, codenum, path, direction=-1)
{
  boundary <- obs_row
  while(!is.na(path$subject[boundary+direction]) & path$subject[boundary+direction] == codenum)
  {
    # Always presume a positive swab is carriage (sometimes Viterbi will output no carriage
    # due to low init and/or transition probs)
    if(path$fitted[boundary+direction] == 2 || path$observed[boundary+direction] == 2)
    {
      boundary <- boundary + direction
    }
    else
    {
      break
    }
  }
  return(boundary)
}

# Average resistance measure over a time period
averageResistance <- function(res_vector)
{
  order <- c("SENSITIVE", "INTERMEDIATE", "RESISTANT")
  ordered_vec <- ordered(res_vector, order)
  levs <- levels(ordered_vec)
  
  m <- median(as.integer(ordered_vec), na.rm = T)
  if (is.na(m))
  {
    median <- NA
  }
  else
  {
    if (floor(m) != m)
    {
      m <- floor(m)
    }
    median <- levs[m]
  }
  return(median)
}

#*********************************#
# Main                            #
#*********************************#

# Read sequence metadata in, do some basic formatting
sequence_metadata <- dplyr::as_data_frame(read.delim("sequence_metadata.txt", 
                                                     header=T, stringsAsFactors = F))

sequence_metadata %<>%
  mutate(specdate=as.Date(sequence_metadata$specdate,"%d-%b-%y")) %>%
  filter(category=="INFANT") %>%
  dplyr::select(lane, specdate, codenum, serotype) %>%
  arrange(codenum, specdate)

saveRDS(sequence_metadata, "sequence_metadata.Rdata")

# Read observation data in
routine_data <- dplyr::as_data_frame(read.csv("20151026_maela_rout_sweep_data.csv",
                                            header = T, stringsAsFactors = F))
immunology_data <- dplyr::as_data_frame(read.csv("20151026_maela_imm_who_data.csv",
                                                 header = T, stringsAsFactors = F))

# 6A/6C are combined (as cannot be separated by sweep). NTs are treated separately
# 15B/C spontaneously intraconvert, so are combined
observed_serotypes <- unique(unlist(c(immunology_data[,10:12], routine_data[,c(10:13,15:17)])))
observed_serotypes <- observed_serotypes[observed_serotypes != "" & !is.na(observed_serotypes) 
                                         & observed_serotypes != "6A" & observed_serotypes != "6C"
                                         & observed_serotypes != "15B" & observed_serotypes != "15C"
                                         & observed_serotypes != "NT"]
observed_serotypes <- c(observed_serotypes, "15B/C")

# Filter and format. Immunology first
# Dates into dates from beginning, kids only, sort by person not date
immunology_data %<>%
  mutate(sampleday=(as.numeric(as.Date(immunology_data$specdate,"%d/%m/%Y") - 
                                 min(as.Date(immunology_data$specdate,"%d/%m/%Y"))))) %>%
  mutate(specdate=as.Date(specdate,"%d/%m/%Y")) %>%
  filter(category=="Infant" & !(codenum=="ARI-0218" & specnum=="08B09098")) %>% # Kids only; duplicate observation to remove
  group_by(codenum) %>%
  mutate(sampleday=sampleday-min(sampleday), collection="immunology") %>%
  rowwise() %>%
  mutate(serotype=paste(combine_sero(c(whoserotype1, whoserotype2, whoserotype3)), collapse = ",")) %>%
  ungroup() %>%
  arrange(codenum, sampleday) %>%
  dplyr::select(codenum, collection, age_d, specdate, sampleday, whopnc, serotype) %>%
  distinct()

# Routine data
# NTs from latex sweeps not detected
routine_data %<>%
  mutate(sampleday=(as.numeric(as.Date(routine_data$specdate,"%d/%m/%Y") - 
                                 min(as.Date(routine_data$specdate,"%d/%m/%Y"))))) %>%
  mutate(specdate=as.Date(specdate,"%d/%m/%Y")) %>%
  filter(!(specnum=="08B03382" | specnum=="09B02164")) %>% # Remove double observations: see problem_rows.txt
  group_by(codenum) %>%
  mutate(sampleday=sampleday-min(sampleday), collection="routine") %>%
  rowwise() %>%
  mutate(serotype=paste(routine_parse_sero(c(sweepserotype1, sweepserotype2, sweepserotype3, sweepserotype4,
                                           whoserotype1, whoserotype2, whoserotype3)),
                          collapse=",")) %>%
  ungroup() %>%
  group_by(codenum, sampleday) %>%
  filter(n() == 1 | whopnc != "") %>% # Identifies cases with two observations, and takes better of two
  ungroup() %>%
  arrange(codenum, sampleday) %>%
  dplyr::select(codenum, collection, age_d, specdate, sampleday, whopnc, serotype) %>%
  distinct()

all_observations <- dplyr::rbind_list(immunology_data, routine_data)

# Normalise times
time_dev <- as.numeric(ungroup(all_observations) %>% summarise(dev=sd(sampleday)))
all_observations <- mutate(all_observations, normday=sampleday/time_dev)
                                 
# Add a variable of whether carriage has been observed before
# Not vectorised
prev_id <- ""
carry_vec <- rep(0, nrow(all_observations))
for (i in 1:nrow(all_observations))
{
  if (prev_id != as.character(all_observations[i,1]))
  {
    prev_id <- as.character(all_observations[i,1])
    carrying <- 0
  }
  
  if(!carrying && nrow(all_observations %>% 
     slice(i) %>% 
     filter(serotype != "")) > 0)
  {
    carrying <- 1
  }
    
  carry_vec[i] <- carrying
}

# Add the carried vector, and remove kids with only one observation
all_observations %<>%
  mutate(carried=carry_vec) %>%
  group_by(codenum) %>%
  filter(n() > 1)

# Compare three below models with AIC

# Set up HMM (three state: clear, carrying, carried)
Qm <- rbind(c(0, 0.1, 0.1),
            c(0, 0, 0.1),
            c(0, 0.1, 0))
ematrix <- rbind(c(0, 0, 0),
                 c(0.1, 0, 0.1),
                 c(0.1, 0, 0))

with_dummy <- rowwise(all_observations) %>%
                mutate(state=createDummy("19F", carried, serotype))
  
nineteenf.msm <- msm(state ~ normday, subject = codenum,
                     data = with_dummy, qmatrix = Qm,
                     est.initprobs = T, ematrix = ematrix, 
                     opt.method = "bobyqa", control=list(maxfun=2000))

# (19F - -2*log-L = 8935.592; k = 8; AIC = 8951.592)
# W/ covar (19F - -2*log-L = 9419.8547; k = 13; AIC = 9445.8547)

# Set up HMM (two state: clear, carrying)
Qm <- rbind(c(0, 0.1),
            c(0.1, 0))
ematrix <- rbind(c(0, 0),
                 c(0.1, 0))

with_dummy <- rowwise(all_observations) %>%
mutate(state=createDummyTwoState("19F", carried, serotype))

nineteenfTwoState.msm <- msm(state ~ normday, subject = codenum,
                             data = with_dummy, qmatrix = Qm,
                             ematrix = ematrix, est.initprobs = T, opt.method = "bobyqa",
                             control=list(maxfun=2000, iprint=3))

# (19F - -2*log-L = 7087.664; k = 5; AIC = 7097.664)

# With covar
nineteenfTwoState.msm <- msm(state ~ normday, subject = codenum,
                             data = with_dummy, qmatrix = Qm,
                             ematrix = ematrix, est.initprobs = T, opt.method = "bobyqa",
                             control=list(maxfun=2000, iprint=3), covariates = ~ carried)

# (19F - -2*log-L = 7081.576; k = 6; AIC = 7093.576)

# Set up HMM (three state: clear, carriage, multiple carriage)
Qm <- rbind(c(0, 0.1, 0),
            c(0.1, 0, 0.1),
            c(0, 0.1, 0))
ematrix <- rbind(c(0, 0, 0),
                 c(0.1, 0, 0),
                 c(0.1, 0.1, 0))

with_dummy <- rowwise(all_observations) %>%
  mutate(state=createDummyThreeState("19F", carried, serotype))

nineteenfThreeState.msm <- msm(state ~ normday, subject = codenum,
                             data = with_dummy, qmatrix = Qm,
                             ematrix = ematrix, est.initprobs = T, opt.method = "bobyqa",
                             control=list(maxfun=2000, iprint=3))

# (19F - -2*log-L = 8933.4596; k = 9; AIC = 8951.4596)

# With covar
nineteenfThreeState.msm <- msm(state ~ normday, subject = codenum,
                             data = with_dummy, qmatrix = Qm,
                             ematrix = ematrix, est.initprobs = T, opt.method = "bobyqa",
                             control=list(maxfun=2000, iprint=3), covariates = ~ carried)

# (19F - -2*log-L = 8921.6363; k = 13; AIC = 8947.636)

# So two state with covariate has best AIC. Run on all serotypes

statetables <- list()
models <- list()
paths <- list()

Qm <- rbind(c(0, 0.1),
            c(0.1, 0))
ematrix <- rbind(c(0, 0),
                 c(0.1, 0))

observed_serotypes <- observed_serotypes[-(observed_serotypes == "15B" | observed_serotypes == "15C")]

foreach(i=1:length(observed_serotypes)) %do% {
  # Dummy variables for each serotype: 
  # 1 clear, 2 carrying, 3 colonised before but currently clear
  with_dummy <- rowwise(all_observations) %>%
    mutate(state=createDummyTwoState(observed_serotypes[i], carried, serotype))
  
  statetables[[observed_serotypes[i]]] <- statetable.msm(state, codenum, with_dummy)
  
  models[[observed_serotypes[i]]] <- msm(state ~ normday, subject = codenum,
                                         data = with_dummy, qmatrix = Qm,
                                         est.initprobs = T, ematrix = ematrix, 
                                         covariates = ~ carried,
                                         opt.method = "bobyqa", control=list(maxfun=4000))
  
  paths[[observed_serotypes[i]]] <- viterbi.msm(models[[observed_serotypes[i]]])
}

# NT can only be detected where WHO serotyping has been done
NT_obs <- filter(all_observations, whopnc != "") %>%
  group_by(codenum) %>%
  filter(n() > 1) %>%
  rowwise() %>%
  mutate(state=createDummyTwoState("NT", carried, serotype))

statetables[['NT']] <- statetable.msm(state, codenum, NT_obs)

models[["NT"]] <- msm(state ~ normday, subject = codenum,
                    data = NT_obs, qmatrix = Qm,
                    est.initprobs = T, ematrix = ematrix,
                    covariates = ~ carried,
                    opt.method = "bobyqa", control=list(maxfun=4000))

paths[["NT"]] <- viterbi.msm(models[["NT"]])

# to fix the misclassification probability at 12%
ematrix <- rbind(c(0, 0),
                 c(0.12, 0))
models[['NT_fixed']] <- msm(state ~ normday, subject = codenum, 
                            data = NT_obs, qmatrix = Qm, 
                            est.initprobs = T, ematrix = ematrix,
                            opt.method = "bobyqa", 
                            control=list(maxfun=4000, iprint=3), fixedpars = c(3))

saveRDS(statetables, file="statetables.Rdata")
saveRDS(models, file="models.Rdata")
saveRDS(paths, file="paths.Rdata")

sojourn.msm(models[['NT_fixed']])[2,]*time_dev
envisits.msm(models[['NT_fixed']], tot = 3.32, ci="normal")

# Ascertained from state tables and model fit
good_models <- c("19F", "23F", "6B", "14", "6A/C", "NT")
bad_models <- observed_serotypes[!(observed_serotypes %in% good_models)]

# Get paths for models which don't have enough data to converge, 
# or error rates overestimate/transition rates underestimated
# Assume 19F, for which we have the most data, is a good HMM for other serotypes, and use the viterbi
# algorithm on them
new_model <- models[['19F']]
foreach(i=1:length(bad_models)) %do% {
  new_model$data <- models[[bad_models[i]]]$data
  paths[[bad_models[i]]] <- viterbi.msm(new_model)
}

saveRDS(models, file="good_models.Rdata")
saveRDS(paths, file="good_paths.Rdata")

# Convert these paths into a pheno file for the metadata
# Not vectorised
sequence_metadata <- filter(sequence_metadata, 
                            codenum %in% all_observations$codenum)

pheno <- rep(0, nrow(sequence_metadata))
age <- rep(0, nrow(sequence_metadata))
carried_out <- rep(0, nrow(sequence_metadata))
foreach(i=1:nrow(sequence_metadata)) %do% {
  # Get the corresponding row in the swab data
  serotype <- sequence_metadata$serotype[i]
  codenum <- sequence_metadata$codenum[i]
  obs_row <- 0
  
  if (serotype == "15B" || serotype == "15C")
  {
    serotype <- "15B/C"
  }
  else if (serotype == "6A" || serotype == "6C")
  {
    serotype <- "6A/C"
  }
  
  # NT uses a different data frame
  if (serotype == "NT")
  {
    obs_row <- which(NT_obs$codenum==codenum 
                     & NT_obs$specdate==sequence_metadata$specdate[i])
    # Check the sequence serotype and swab serotype match
    if(!(serotype %in% unlist(strsplit(NT_obs$serotype[obs_row], split=","))))
    {
      # Check warnings to make sure there are no multiple carriage mismatches
      print(paste0("WARNING: Sample ", i, ",", obs_row, " NT mismatch ", serotype, "_", all_observations$serotype[obs_row]))
      serotype <- strsplit(NT_obs$serotype[obs_row], split=",")[[1]]
    }
  }
  else
  {
    obs_row <- which(all_observations$codenum==codenum
                     & all_observations$specdate==sequence_metadata$specdate[i])
    # Check the sequence serotype and swab serotype match
    if(!(serotype %in% unlist(strsplit(all_observations$serotype[obs_row], split=","))))
    {
      # Check warnings to make sure there are no multiple carriage mismatches
      print(paste0("WARNING: Sample ", i, ",", obs_row, " serotype mismatch ", serotype, "_", all_observations$serotype[obs_row]))
      serotype <- strsplit(all_observations$serotype[obs_row], split=",")[[1]]
      
      if(serotype == "NT")
      {
        obs_row <- which(NT_obs$codenum==codenum 
                         & NT_obs$specdate==sequence_metadata$specdate[i])
      }
    }
  }
  
  if(length(obs_row) > 1)
  {
    print(paste0("WARNING: Sample ", i, " multiple obs: ", obs_row))
    obs_row<-obs_row[1]
  }
  
  if(length(serotype) > 1)
  {
    print(paste0("WARNING: Sample ", i, " multiple serotypes ", serotype))
  }
  
  # Get date at the start of carriage
  if (serotype != "NT")
  {
     age[i] <- as.numeric(all_observations[obs_row, "age_d"])
     carried_out[i] <- getCarried(obs_row, codenum, paths[[serotype]], all_observations)
     start <- getDateBoundary(obs_row, codenum, paths[[serotype]], all_observations, direction=-1)
     end <- getDateBoundary(obs_row, codenum, paths[[serotype]], all_observations, direction=1)
  }
  else
  {
    age[i] <- as.numeric(NT_obs[obs_row, "age_d"])
    carried_out[i] <- getCarried(obs_row, codenum, paths[[serotype]], NT_obs)
    start <- getDateBoundary(obs_row, codenum, paths[[serotype]], NT_obs, direction=-1)
    end <- getDateBoundary(obs_row, codenum, paths[[serotype]], NT_obs, direction=1)
  }
  
  pheno[i] <- as.numeric(end-start)
}

# Log transform carriage length
pheno <- log(pheno)
age <- log(age)

metadata_out <- sequence_metadata %>%
                  mutate(phenotype=pheno, fid=lane) %>%
                  dplyr::select(lane, fid, phenotype)

covariates_out <- sequence_metadata %>%
                    mutate(age=age, carried=carried_out, fid=lane) %>%
                    dplyr::select(lane, fid, age, carried)

write.table(metadata_out, file="length.pheno", 
            quote = F, sep="\t", row.names=F, col.names=F)
saveRDS(metadata_out, file = "metadata_out.Rdata")

write.table(covariates_out, file="covariates.txt", 
            quote = F, sep="\t", row.names=F, col.names=F)

# All carriage episode lengths
pathcopy <- paths
observed_serotypes <- c(observed_serotypes, "NT")

# Read in all resistance swabs
resistance_metadata <- dplyr::as_data_frame(read.delim("ARI_AMR_data.csv", sep = ",",
                                                     header=T, stringsAsFactors = F))
resistance_metadata %<>%
  filter(Category=="Infant") %>%
  mutate(specdate=as.Date(Specimen.date,"%d/%m/%Y"))

X <- list()
y <- list()
for (serotype in names(pathcopy))
{
  # clean up paths by removing false positives - force fitted to be observed 
  pathcopy[[serotype]][pathcopy[[serotype]][,3] == 2, 4] = 2
  
  # Add in fitted state from above or below, to find start and end of carriage
  pathcopy_new <- as.data.table(pathcopy[[serotype]])
  pathcopy_new[, ("lag") := shift(fitted, 1), by=subject]
  pathcopy_new[, ("lead") := shift(fitted, 1, type = "lead"), by=subject]

  start_points <- which(pathcopy_new$fitted == 2 & (pathcopy_new$lag != 2 | is.na(pathcopy_new$lag)))
  end_points <- which(pathcopy_new$fitted == 2 & (pathcopy_new$lead != 2 | is.na(pathcopy_new$lead)))
  
  if (length(start_points) > 0)
  {
    # Now construct the y and X
    carried <- rep(0, length(start_points))
    length <- rep(0, length(start_points))
    age <- all_observations[start_points, "age_d"]
    resistant <- array(data = NA, dim=c(length(start_points),7)) 
    colnames(resistant) = c("Ceftriaxone", "Chloramphenicol", "Clindamycin", "Erythromycin", "Penicillin", "Sulpha.trimethoprim", "Tetracycline")
    for (i in 1:length(start_points))
    {
      if (serotype != "NT")
      {  
        obsframe <- all_observations
      }
      else
      {
        obsframe <- NT_obs
      }
      codenum <- as.character(obsframe[start_points[i],"codenum"])
      carried[i] <- getCarried(start_points[i], codenum, pathcopy[[serotype]], obsframe)
      start <- getDateBoundary(start_points[i], codenum, pathcopy[[serotype]], obsframe, direction=-1)
      end <- getDateBoundary(start_points[i], codenum, pathcopy[[serotype]], obsframe, direction=1)
      length[i] <- as.numeric(end-start)
      
      # Get the AB phenotype
      if (serotype == "15B/C")
      {
        swabs <- resistance_metadata %>% filter(specdate <= obsframe[end_points[i], "specdate"][[1]] & 
                                                  specdate >= obsframe[start_points[i], "specdate"][[1]] &
                                                  Study.code == codenum &
                                                  (serotype == "15B" | serotype == "15C" | serotype == "15B/C")) %>%
          dplyr::select(Ceftriaxone, Chloramphenicol, Clindamycin, Erythromycin, Penicillin, Sulpha.trimethoprim, Tetracycline)
      }
      else
      {
        swabs <- resistance_metadata %>% filter(specdate <= obsframe[end_points[i], "specdate"][[1]] & 
                                                  specdate >= obsframe[start_points[i], "specdate"][[1]] &
                                                  Study.code == codenum &
                                                  serotype == serotype) %>%
          dplyr::select(Ceftriaxone, Chloramphenicol, Clindamycin, Erythromycin, Penicillin, Sulpha.trimethoprim, Tetracycline)
      }
      resistant[i,] <- apply(swabs, 2, averageResistance)
      
    }
    X[[serotype]] <- data.frame(age=age, carried=carried, serotype=serotype, resistant)
    y[[serotype]] <- data.frame(length=length)
  }
}

all_X <- do.call(rbind, X)
all_y <- do.call(rbind, y)

saveRDS(all_X, file = "all_X.Rdata")
saveRDS(all_y, file = "all_y.Rdata")
