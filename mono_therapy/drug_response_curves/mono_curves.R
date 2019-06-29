home_dir <- "~/Desktop/compass/resource/matrices/"
out_dir <- "~/Desktop/compass/resource/plots/"
plotting_flag <- TRUE

library(PharmacoGx)
library(Cairo)
library(devtools)
library(readxl)
library(openxlsx)

options(stringsAsFactors = FALSE)

identifiers <- read_xlsx("~/Desktop/compass/resource/List of COMPASS PDO Identifiers.xlsx", sheet=1)
identifiers$`Study ID` <- sub("^(.{4})", "\\1_", identifiers$`Study ID`)
ids <- identifiers$`Study ID`
names(ids) <- identifiers$`Model ID`

drugDoseResponseCurve2 <- function(drug, cellline, conc, viability, aac, nreplicates, nconcentration, plate) {    
  doses <- list(); responses <- list(); legend.values <- list(); j <- 0; pSetIndex <- list()
  doses[[1]] <- conc
  responses[[1]] <- viability
  
  dose.range <- c(10^100 , 0)
  viability.range <- c(0 , 10)
  for(i in 1:length(doses)) {
    dose.range <- c(min(dose.range[1], min(doses[[i]], na.rm=TRUE), na.rm=TRUE), 
                    max(dose.range[2], max(doses[[i]], na.rm=TRUE), na.rm=TRUE))
    viability.range <- c(0, max(viability.range[2], max(responses[[i]], na.rm=TRUE), na.rm=TRUE))
  }
  x1 <- 10 ^ 10; x2 <- 0
  
  plot(NA, xlab="Concentration (uM)", ylab="% Viability", axes =FALSE,
       main=sprintf("%s, %s (%s) treated with %s", cellline, plate, ids[cellline], drug), 
       log="x", ylim=c(0,115), xlim=dose.range, cex=1, cex.main=1, 
       sub=paste0("#doses = ", nconcentration, ", #replicates = ", nreplicates))
  
  magicaxis::magaxis(side=1:2, frame.plot=TRUE, tcl=-.3, majorn=c(5,3), minorn=c(5,2))
  
  if (length(doses) > 1) {
    rect(xleft=x1, xright=x2, ybottom=viability.range[1] , ytop=viability.range[2] , 
         col=rgb(240, 240, 240, maxColorValue = 255), border=FALSE)
  }
  
  points(doses[[1]],responses[[1]], col="black", pch=20)
  
  log_logistic_params <- PharmacoGx::logLogisticRegression(conc = doses[[i]], viability = responses[[i]])
  log10_x_vals <- PharmacoGx:::.GetSupportVec(log10(doses[[i]]))
  lines(10 ^ log10_x_vals, PharmacoGx:::.Hill(log10_x_vals, pars=c(log_logistic_params$HS, 
                                                                   log_logistic_params$E_inf/100, log10(log_logistic_params$EC50))) * 100 ,lty=1, lwd=2, col="red")
  
  legend("bottomleft", legend=paste0("AAC = ", format(round(aac, 4), nsmall = 4)), bty="n")
  
  abline(h=50, lty=2, lwd=0.5, col="grey")
  abline(h=100, lty=2, lwd=0.5, col="grey")
}


options(stringsAsFactors = FALSE)

all_files <- list.files(home_dir)
all_files <- sub('.txt', '', all_files)

aac_output <- c()
for(idx in 1:length(all_files)) {
  buff <- unlist(strsplit(all_files[idx], "__"))
  drug_name <- buff[1]
  cell_line <- buff[2]
  plate_number <- buff[3]
  
  header_scan <- scan(paste0(home_dir, all_files[idx], ".txt"), nlines = 8, what = character(), sep="\n")
  headers <- data.frame(strsplit(header_scan, ":"))[2,]
  names(headers) <- gsub("#", "", data.frame(strsplit(header_scan, ":"))[1,])
  print(headers$sampleID)
  
  # modify
  untreated_control <- as.numeric(headers[6])
  dose <- as.numeric(headers[7]) 
  fold <- as.numeric(headers[8])
  
  matrix <- read.table(paste0(home_dir, all_files[idx], ".txt"), sep="\t")
  no_of_replicates <- ncol(matrix)
  no_of_concentrations <- nrow(matrix)
  
  replicate <- (matrix/untreated_control) * 100
  
  dose_scale <- dose
  dose_temp <- dose
  for (s in 1:(no_of_concentrations-1)) {
    x <- dose_temp
    y <- x/fold
    dose_scale <- c(dose_scale, y)
    dose_temp <- y
  }
  
  input <- data.frame(concentration = rep(dose_scale, times=3), viability = unlist(replicate))
  
  aac <- computeAUC(concentration = input$concentration, viability = input$viability, viability_as_pct = TRUE, verbose = F) / 100
  aac_output <- rbind(aac_output, c(drug_name, cell_line, plate_number, no_of_concentrations, no_of_replicates, ids[cell_line], aac))
  
  if(plotting_flag) {
    png_file_name = paste0(out_dir, drug_name, "__", cell_line, "__", plate_number, "__", ids[cell_line], ".png")
    Cairo(width = 1000, height = 800, file = png_file_name, type = "png", bg = "white", canvas = "white", units = "px", dpi = 200)
    drugDoseResponseCurve2(drug = drug_name, cellline = cell_line, conc = input$concentration, viability = input$viability,
                           aac = aac, nreplicates = no_of_replicates, nconcentration = no_of_concentrations, plate = plate_number)
    dev.off()
  }
}
colnames(aac_output) <- c("Drug", "Sample", "Plate", "Doses", "Replicates", "Patient", "AAC")
write.table(aac_output, file=paste0(out_dir, "AAC_profiles.txt"), row.names=F, col.names=T, sep="\t", quote=F)

q("no")