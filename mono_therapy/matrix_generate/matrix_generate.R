require(MASS)
library(devtools)
library(readxl)
library(openxlsx)
library(strex)
library(stringi)
library(strex)
library(PharmacoGx)
library(Cairo)
options(scipen = 999)

#compounds used

compounds <- c("GEM", "TAX", "5FU", "OXA", "IRI", "AFB", "OLA", "CPT-11")

#get excel sheet names (samples)

excel <- "~/Desktop/compass/resource/COMPASS Annotated Raw Data 20171201 - 20190432.xlsx"
patients_sheet <- excel_sheets(path = excel)
patients_sheet <- patients_sheet[2:43]

#concatenate samples with patient ID to form unique ID

identifiers <- read_xlsx("~/Desktop/compass/resource/List of COMPASS PDO Identifiers.xlsx", sheet=1)
identifiers$`Study ID` <- sub("^(.{4})", "\\1_", identifiers$`Study ID`)
patient_identifiers <- paste0(sub("\\-.*", "", identifiers$`Study ID`), "(", identifiers$`Model ID`, ")")

#create matrix + header for each sample in each sheet:

for (aa in patients_sheet) {
  raw <- read.xlsx("~/Desktop/compass/resource/COMPASS Annotated Raw Data 20171201 - 20190432.xlsx", sheet=aa, skipEmptyRows = T, skipEmptyCols = T)
  raw[1,] <- NA
  raw[,1] <- NA
  raw <- raw[, colSums(is.na(raw)) != nrow(raw)]
  # untreated <- grep("AVE", raw)
  # untreated <- as.numeric(rownames(raw[grep("AVE", raw[,untreated[1]]),]))
  # raw <- raw[1:(untreated -2),]
  row_remove <-vector()
  for (row in 1:nrow(raw)){
    
    if (length(grep("uM", raw[row,])) > 1){
      row_remove <- c(row_remove, row)
    }
    
    if (length(grep("nM", raw[row,])) > 1){
      row_remove <- c(row_remove, row)
    }
    
  }
  
  raw <- raw[1:(row_remove[2]-1),]
  if (length(grep("PPTO", raw)) > 0) {
    extra_sample_col <- grep("PPTO", raw)
    extra_sample_row <- as.numeric(rownames(raw[grep("PPTO", raw[,extra_sample_col[1]]),]))
    if (length(extra_sample_row) > 1) {
      raw[tail(extra_sample_row, n=1) ,extra_sample_col] <- NA
      
    }
  }
  
  
  if (length(grep("repeat", raw)) > 0) {
    notes_col <- grep("repeat", raw)
    notes_row <- as.numeric(rownames(raw[grep("repeat", raw[,notes_col[1]]),]))
    raw[notes_row, notes_col] <- NA
  }
  #Gather metadata for each drug + sample and convert to matrix:
  
  drug_found <- grep("M", raw)
  drugs_all <- vector()
  for (d in 1:length(drug_found)) {
    d_name <- raw[,drug_found[d]][grep("M", raw[,drug_found[d]])]
    drugs_all <- c(drugs_all, d_name)
    
    if (grepl("^M", drugs_all) == TRUE){
      drugs_all <- drugs_all[-grep("^M", drugs_all)]
    }
    
    if (grepl("^2019", drugs_all) == TRUE){
      drugs_all <- drugs_all[-grep("^2019", drugs_all)]
    }
    
  }
  
  
  
  drugs_dose <- sub(" .*", "", drugs_all)
  
  
  if (length(grep("PPTO", raw)) > 0) {
    plate_found <- grep("PPTO", raw)
    plates_all <- vector()
    for (p in 1:length(plate_found)) {
      p_name <- raw[,plate_found[p]][grep("PPTO", raw[,plate_found[p]])]
      plates_all <- c(p_name, plates_all)
    }
  } else {
    
    
    plate_found <- grep("XDO", raw)
    plates_all <- vector()
    for (p in 1:length(plate_found)) {
      p_name <- raw[,plate_found[p]][grep("XDO", raw[,plate_found[p]])]
      plates_all <- c(p_name, plates_all)
    }
    
  }
  
  drug_columns <- grep("M", raw)
  dd <- 1
  
  for (col in 1:length(drug_columns)) {
    
    if (col == length(drug_columns))
    { raw_data <- raw[,drug_columns[col]:ncol(raw)]
    print(raw_data)
    } else {
      raw_data <- raw[,drug_columns[col]:(drug_columns[(col + 1)]-1)]
      print(raw_data)
    }
    
    if (length(plates_all[which(plates_all %in% raw_data[,1])]) == 0){
      plate <- plate_used
    } else {
      plate <- plates_all[which(plates_all %in% raw_data[,1])]
      plate_used <- vector()
      plate_used <- c(plate, plate_used)
    }
    
    if (length(grep("Untreated", raw_data)) == 1){
      
      untreated_col <- grep("Untreated", raw_data)
      untreated_row <- as.numeric(rownames(raw[grep("Untreated", raw_data[,untreated_col[1]]),]))
      untreated_values <- vector()
      for (unt in 1:length(untreated_col)){
        untreated_values <- c(untreated_values,raw_data[untreated_row,(untreated_col[unt]+1)] )
        raw_data[untreated_row,untreated_col[unt]] <- NA
        raw_data[untreated_row,(untreated_col[unt]+1)] <- NA
      }
      
      untreated_ave_id <- untreated_values[1]
    } else{
      
      untreated_ave_id <- tail(untreated_values, n=1)
    }
    
    
    drug <- drugs_dose[dd]
    drugs_all2 <- sub("5-FU", "FU", drugs_all)
    drugs_all2 <- sub("5FU", "FU", drugs_all2)
    drugs_all2 <- sub("CPT-11", "CPT", drugs_all2)
    drugs_all2 <- sub("CPT11", "CPT", drugs_all2)
    dose <- str_nth_number(drugs_all2[dd], n = 1, decimals = TRUE)
    dose_temp <- dose
    fold <- str_nth_number(drugs_all2[dd], n = 2, decimals = TRUE)
    
    raw_final <- raw_data[(which(raw_data == drugs_all[dd])+1):nrow(raw_data),]
    raw_final <- Filter(function(x) !all(is.na(x)), raw_final)
    raw_final <- raw_final[complete.cases(raw_final),]
    
    
    dose_scale <- dose
    
    for (s in 1:(nrow(raw_final)-1)) {
      
      x <- dose_temp
      y <- x/fold
      dose_scale <- c(dose_scale,y)
      dose_temp <- y
    }
    
    
    #mrx1 <- matrix(c(dose_scale), nrow=length(dose_scale), ncol = 1)
    #mrx2 <- cbind(mrx1, raw_final)
    
    date <- str_nth_number(plate, n = 1)
    if(length(grep("PPTO",raw)) > 0){
      sample_id <- str_nth_number(plate, n = 2)
      sample_id <- print(paste0("PPTO.", sample_id))
    } else {
      sample_id <- "XDO-COMP0037"
      
    }
    plate <- str_nth_number(plate, n = 3)
    
    if (length(grep("uM", drugs_all[1])) == 1) {
      unit <- "uM"
    } else {
      unit <- "nM"
    }
    
    drug <- sub("-", "", drug)
    setwd("~/Desktop/compass/resource/matrices")
    plate <- paste0("p", plate)
    identifier <- identifiers$`Study ID`[which(identifiers$`Model ID` == sample_id)]
    filename <- paste(drug, sample_id, plate, identifier, sep = "__", collapse = NULL)
    colnames(raw_final) <- NULL
    
    date_id <- print(paste0("#date:", date))
    sample_id <- print(paste0("#sampleID:", sample_id))
    plate_id <- print(paste0("#phase:", plate))
    drug_id <- print(paste0("#drug:", drug))
    unit_id <-  print(paste0("#unit:", unit))
    untreat_id <- print(paste0("#untreated_control:", untreated_ave_id))
    start_dose_id <- print(paste0("#start_dose:", dose_scale[1]))
    fold_id <- print(paste0("#fold:", fold))
    
    
    if (grepl("+", drugs_all[dd], fixed=TRUE) == TRUE){
      print("skip drug combo")
      
    }else{
      
      write.matrix(raw_final, file=paste(filename, ".txt", sep=""), sep="\t")
      
      fConn <- file(paste(filename, ".txt", sep=""), 'r+') 
      Lines <- readLines(fConn) 
      writeLines(c(date_id, sample_id, drug_id, plate_id, unit_id, untreat_id, start_dose_id, fold_id, Lines), con = fConn) 
      close(fConn) 
      
      
    }
    
    dose <- 0
    fold <- 0
    dd <- dd + 1
  }
  
  dose <- 0
  fold <- 0
  dd <- 0
  date_id <- ""
  sample_id <- ""
  plate_id <- ""
  drug_id <- ""
  unit_id <-  ""
  start_dose_id <- ""
  fold_id <- ""
  
}