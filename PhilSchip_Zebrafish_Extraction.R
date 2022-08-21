## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see https://www.gnu.org/licenses/.


## Changing the working directory to the folder of the evaluation script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#1: Pre-analysis parameters-----

## If TRUE uses previously processed raw data (located in the same folder and named as "Preprocessed_Data")
Preprocess.done <- FALSE

## If TRUE uses batch correction on the picked features using pooled QC samples
batchcor <- TRUE

## Sample injection order, arranged based on alphabetical order of classes and samples
order <- c(2,8,10,4,16, 32,33,35,36,3, 12,23,11,34,28, 26,30,22,24,6,
           21,18,14,27,29, 5,15,17,9,20, 1,7,13,19,25,31,37)

## Marks every sample group with a specific colour, can be adapted freely
colours_classes <- c("Extraction I" = "#0926B7", "Extraction II" = "#FF217B",
                     "Extraction III" = "#05992C", "Extraction III_B" = "#7DF212",
                     "Extraction IV" = "#FF1705",  "Extraction IV_B" = "#FA9106", "QC" = "#000000")

column <- "PhenHex"      ## used for some result annotation (PhenHex and HILIC supported)
polarity <- "negative"   ## used for some result annotation (positive and negative supported)

QC.groupsize <- 7        ## amount of measured QC samples

## Determines if internal standards (IS), or a curated list of molecules of interest (MoInt) are evaluated
Target_eval <- "MoInt" ## allows "MoInt" and "IS"

#2: Loading of needed packages-----

library("tidyverse")
library("xcms")
library("CAMERA")
library("ggrepel")
library("Rtsne")
library("gplots")
library("MASS")
library("caret")
library("matrixStats")
library("remotes")
library("devtools")
library("ggpubr") 
library("writexl")
library("readxl")


## Loads metaboLib database from Manier 2020 (https://github.com/saskema/metaboLib)
source("metaboLib.R")

#3: Folder setup and loading of raw data files-----

## Name and create results folder
dir.create("Statistics_Results", showWarnings = FALSE)

## Loading of the datafiles (stored groupwise in subfolders)
files <- list.files(path = getwd(), recursive = TRUE, full.names = TRUE, pattern = ".mzXML")

## Loading of the optimized centwave parameters
## using a modified version of the CentwaveOpt algorithm from Manier 2020 (https://github.com/saskema/centWaveOpt)
## which can be found at https://github.com/philschip/centWabeOpt_modified
parameter <- scan("centWaveOptResults/parameter.csv", sep = ",")

#4: Preprocessing of raw data-----

## Loads already preprocessed data if available
if(Preprocess.done ==  TRUE) {
  load("Preprocessed_Data")
}else{ ##or performs peak picking, grouping, retention time correction and peak filling
  set <- xcmsSet(files, method = "centWave", peakwidth = c(parameter[1], parameter[2]),
                 ppm = parameter[3], snthresh = parameter[4], mzdiff = parameter[5],
                 prefilter = c(parameter[6], parameter[7]))
  set <- group.density(set, bw = parameter[8])
  set <- retcor(set, method = "obiwarp", plottype = "deviation")
  set <- group.density(set, bw = parameter[8])
  set <- retcor(set, method = "obiwarp", plottype = "deviation")
  set <- group.density(set, bw = parameter[8])
  set <- fillPeaks(set)
}

## Saving the preprocessed data as shortcut for future analysis
save(set, file = "Preprocessed_Data")

## Annotate features using the CAMERA package and save the results as csv file
set.CAMERA <- CAMERA::annotate(set, polarity  = polarity)
write.csv(getPeaklist(set.CAMERA), file = "peaklist_CAMERA.csv")

#5: Bacis processing-----

## Setup of classes and class levels based on the raw data
Classes <- set$class
class.levels <- levels(as.factor(set$class))

## Generating a peaklist out of the preprocessed data
peaklist <- peakTable(set)
rownames(peaklist) <- groupnames(set)

peak.matrix <- peaklist[,(8 + length(class.levels)):length(peaklist)]

## Log10 Transformation and feature abundance imputation, acc. to Wehrens 2016
matrix.transformed <- log10.matrix(peak.matrix)

## Batch correction if marked as TRUE in pre-analysis parameters
if(batchcor == TRUE) {
  matrix.batchcor <- norm.single.batch(matrix = matrix.transformed, class = set$class, order = order,
                                           type = "pooled", standard = (peaklist$QC == 7), poly = 3, plot = 50)
  
  dev.off()
  
} else {
  matrix.batchcor <- as.matrix(matrix.transformed)
}

#6: Statistical analysis-----

pdf("Statistics_Results/Evaluation.pdf", width = 8, height = 5)

## Generates bonferroni correction factor based on amount of detected features
Bonf.Corr <- nrow(matrix.batchcor)

## Calculates significance by ANOVA, including Bonferroni Correction of the initial p-value
significance.table <- t(matrix.batchcor)[set$class != "QC",] 
significance.result <- evalANOVA(features = significance.table,
                                 names = rownames(matrix.batchcor),
                                 classes = set$class[set$class != "QC"],
                                 p.value = (0.001/Bonf.Corr),
                                 plot = TRUE,
                                 labels = TRUE)

## Filters all detected features by ANOVA significance
matrix.sig <- column_to_rownames(filter(rownames_to_column(as.data.frame(matrix.batchcor)),
                                                    significance.result$Significant == "TRUE"))
matrix.sig <- t(matrix.sig)

## Perform multivariate statistics on significant features
if (length(matrix.sig) > 1) {
  
  ## Plot heatmap with dendrogram
  heatmap.2(x = t(matrix.sig), scale = "row", distfun = dist, margins = c(7,8),
            trace = "none", cexRow=0.9, cexCol=0.9)
  
  
  ## PCA and tSNE Parameters
  PCA.center <- TRUE
  PCA.scale <- FALSE
  
  tSNE.center <- TRUE
  tSNE.scale <- TRUE
  tSNE.perplexity <- 4
  
  ## Perform t-SNE analysis and plot results
  tsne <- Rtsne(matrix.sig, dim = 2, perplexity = tSNE.perplexity, epoch = 100,
                verbose = FALSE, pca_center = tSNE.center, pca_scale = tSNE.scale)
  print(
    ggplot(data = as.data.frame(tsne$Y), aes(x = tsne$Y[,1], y = tsne$Y[,2], col = Classes)) +
      geom_point(size = 3) +
      scale_color_manual(values = colours_classes) +
      ggtitle("t-SNE") +
      xlab("X") +
      ylab("Y") +
      theme_bw()
  )
  
  ## Perform PCA and plot results
  pca <- evalPCA(matrix.sig, plot = TRUE, annotate.scores = TRUE, annotate.loadings = 15,
                 classes = set$class, scale = PCA.scale, center = PCA.center)
  
  ## Perform PC-DFA
  n.pc <- length(which(pca$sdev^2 > mean(pca$sdev^2)))
  if(n.pc < 3) {
    
    n.pc <- 3
    
  }
  pc.ldam <- lda(as.matrix(pca$x[,1:n.pc]), set$class)
  pc.lda.loadings <- pca$rotation[,1:n.pc] %*% pc.ldam$scaling
  
  ## Calculate prediction accuracy and cohens kappa
  monte.carlo <- train(pca$x[,1:n.pc], set$class,
                       method = "lda",
                       trControl = trainControl(method = "LGOCV"))
  cvaccuracy <- as.numeric(round(monte.carlo$results[2]*100, 0))
  cvkappa <- as.numeric(round(monte.carlo$results[3]*100, 0))
  
  ## Predict model for test data
  p <- predict(pc.ldam, as.data.frame(pca$x[,1:n.pc]))
  
  ## Plot PC-DFA Scores without labels
  print(
    ggplot(data = as.data.frame(p$x), aes(x = LD1, y = LD2)) +
      geom_point(size = 3, aes(col = Classes)) +
      scale_color_manual(values = colours_classes) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      ggtitle(substitute(paste("PCs used: ", n_pc, ", ", "Accuracy: ", cvacc, "%, ",
                               kappa, " = ", cvka, "%"),
                         list(n_pc = n.pc, cvacc = cvaccuracy, cvka = cvkappa))) +
      theme_bw()
  )
  
  classes <- Classes  
  ## Plot PC-DFA Loadings
  print(
    ggplot(data = as.data.frame(pc.lda.loadings), aes(x = LD1, y = LD2)) +
      geom_point(size = 3, aes(col = "red"), show.legend = FALSE) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      geom_text_repel(point.padding = 0.2, data = as.data.frame(pc.lda.loadings),
                      aes(label = rownames(pc.lda.loadings))) +
      ggtitle(substitute(paste("PCs used: ", n_pc, ", ", "Accuracy: ", cvacc, "%, ",
                               kappa, " = ", cvka, "%"),
                         list(n_pc = n.pc, cvacc = cvaccuracy, cvka = cvkappa))) +
      theme_bw()
)
}

dev.off()


peaklist.Camera <- getPeaklist(set.CAMERA)
rownames(peaklist.Camera) <- groupnames(set)

## Generating a csv file of all significant features
if(batchcor == TRUE) {
  peaklist.Camera_a <- peaklist.Camera[which(peaklist.Camera$QC == QC.groupsize), ]
  peaklist.Camera_b <- peaklist.Camera[which(peaklist.Camera$QC != QC.groupsize), ]
  
  
  significant.features <- dplyr::select(significance.result, feature, pvalue)
  
  significant.features$column <- column
  significant.features$polarity <- polarity
  
  significant.features$mz <- rbind(peaklist[which(peaklist$QC == QC.groupsize),],
                                   peaklist[which(peaklist$QC != QC.groupsize),])$mz
  significant.features$rt <- rbind(peaklist[which(peaklist$QC == QC.groupsize),],
                                   peaklist[which(peaklist$QC != QC.groupsize),])$rt
  
  significant.features$adduct <- rbind(peaklist.Camera_a, peaklist.Camera_b)$adduct
  significant.features$isotopes <- rbind(peaklist.Camera_a, peaklist.Camera_b)$isotopes
  significant.features$pcgroup <- rbind(peaklist.Camera_a, peaklist.Camera_b)$pcgroup
  
  significant.features <- significant.features[which(significance.result$Significant == TRUE),]
  significant.features$abs.loading <- sqrt(abs(pc.lda.loadings[,1])^2 + abs(pc.lda.loadings[,2])^2)
  significant.features <- arrange(significant.features, desc(abs.loading))

} else {
  
  significant.features <- dplyr::select(significance.result, feature, pvalue)
  
  significant.features$column <- column
  significant.features$polarity <- polarity
  
  significant.features$mz <- peaklist$mz
  significant.features$rt <- peaklist$rt
  
  significant.features$adduct <- peaklist.Camera$adduct
  significant.features$isotopes <- peaklist.Camera$isotopes
  significant.features$pcgroup <- peaklist.Camera$pcgroup
  
  significant.features <- significant.features[which(significance.result$Significant == TRUE),]
  significant.features$abs.loading <- sqrt(abs(pc.lda.loadings[,1])^2 + abs(pc.lda.loadings[,2])^2)
  significant.features <- arrange(significant.features, desc(abs.loading))
}

write.csv(significant.features, file = "Statistics_Results/Significant_features.csv")

#7: Analysis of feature count-----

## Plot feature count (peaknumber)
peaknumber <- as.data.frame(t(peak.matrix))
peaknumber[peaknumber == 0] <- NA
peaknumber[is.na(peaknumber) == FALSE] <- 1
peaknumber[is.na(peaknumber) == TRUE] <- 0
peaknumber$sum <- apply(peaknumber, 1, sum)
peaknumber$extraction <- set$class
peaknumber_order <- c("Extraction_I", "Extraction_II", "Extraction_III", "Extraction_IV",
                      "Extraction_III_B", "Extraction_IV_B", "QC")

## Plot peaknumber of extractions without bead homogenization and test significance against extraction I
p.peaknumber <- ggplot(peaknumber, aes(x = factor(extraction, level = peaknumber_order), y = sum)) +
  geom_boxplot(aes(group = extraction)) +
  stat_compare_means(method = "anova", size = 5, label.x = 0.75)+
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "Extraction_I", size = 6) +
  xlab("Extraction Method") +
  ylab("Feature Count [AU]") +
  scale_x_discrete(limits = c("Extraction_I", "Extraction_II", "Extraction_III", "Extraction_IV")) +
  theme_bw(base_size = 16)


## Plot peaknumber of all extractions
## and calculate significance between extractions with and without beads (III, IV and III_B, IV_B)
p.peaknumber_homogenization <- ggplot(peaknumber, aes(x = factor(extraction, level = peaknumber_order), y = sum)) +
  geom_boxplot(aes(group = extraction)) +
  stat_compare_means(comparisons = list(c("Extraction_III", "Extraction_III_B")), label = "p.signif", method = "t.test", ref.group = "Extraction_III", color="#FF0000", size = 6, label.y = 6750, bracket.size = 1.2) +
  stat_compare_means(comparisons = list(c("Extraction_IV", "Extraction_IV_B")), label = "p.signif", method = "t.test", ref.group = "Extraction_IV", color="#006400", size = 6, label.y = 6900, bracket.size = 1.2) +
  xlab("Extraction Method") +
  ylab("Feature Count [AU]") +
  theme_bw(base_size = 16)


ggsave("Peaknumbers_procedure.pdf", p.peaknumber, width = 8, height = 6, device = "pdf")
ggsave("Peaknumbers_homogenization.pdf", p.peaknumber_homogenization, width = 8, height = 6, device = "pdf")

Peaknumber_Results <- as.data.frame(peaknumber$sum)
colnames(Peaknumber_Results) <- "Feature_count"
Peaknumber_Results$Sample <- rownames(peaknumber)
write_xlsx(Peaknumber_Results, "Feature_Count.xlsx")

#8: Evaluation of Internal standards or MoInts

## Generating of result folder name and reading of inclusion list for targeted evaluation
if (Target_eval == "MoInt") {
  Target_Folder <- "MoInt_Evaluation"
  if (column == "PhenHex" & polarity == "positive"){
    exceldata = read_excel("MoInt_Inclusion_PhenHex_positive.xlsx")
  }
  if (column == "PhenHex" & polarity == "negative"){
    exceldata = read_excel("MoInt_Inclusion_PhenHex_negative.xlsx")
  }
  if (column == "HILIC" & polarity == "positive"){
    exceldata = read_excel("MoInt_Inclusion_HILIC_positive.xlsx")
  }
  if (column == "HILIC" & polarity == "negative"){
    exceldata = read_excel("MoInt_Inclusion_HILIC_negative.xlsx")
  }
}

if (Target_eval == "IS") {
  Target_Folder <- "IS_Evaluation"
  if (column == "PhenHex" & polarity == "positive"){
    exceldata = read_excel("IS_Inclusion_PhenHex_positive.xlsx")
  }
  if (column == "PhenHex" & polarity == "negative"){
    exceldata = read_excel("IS_Inclusion_PhenHex_negative.xlsx")
  }
  if (column == "HILIC" & polarity == "positive"){
    exceldata = read_excel("IS_Inclusion_HILIC_positive.xlsx")
  }
  if (column == "HILIC" & polarity == "negative"){
    exceldata = read_excel("IS_Inclusion_HILIC_negative.xlsx")
  }
}

## Creating the folder for upcoming results
dir.create(Target_Folder)

Targetnames <- exceldata[ ,1]
Target.mz <- exceldata[ ,2]
Target.rt <- exceldata[ ,3]


## Preparing dataframes for upcoming analysis
Targetlength <- length(t(Targetnames))
Sample.amount <- nrow(peaknumber)

Target_data_Area <- data.frame(Extraction_I =numeric(), Extraction_II =numeric(), Extraction_III =numeric(),
                               Extraction_III_B =numeric(), Extraction_IV =numeric(), Extraction_IV_B =numeric(),
                               QC =numeric(), stringsAsFactors=FALSE)
Target_data_Area[Targetlength,1] <- NA
row.names(Target_data_Area) <- Targetnames$Name

Target_data_Stabw <- Target_data_Area[1:Targetlength, ]
Target_data_CV <- Target_data_Area

Target_data_Stabw[Targetlength,1] <- NA
Target_data_CV[Targetlength,1] <- NA

complete_data <- data.frame(matrix(vector(), length(t(Targetnames)), Sample.amount))

## Looping through every loaded target, extracting features by retention time and ppm deviation
## and generating groupwise Mean
for(i in 1:Targetlength){
  mz.Target <- as.numeric(Target.mz[i,1])
  rt.Target <- as.numeric(Target.rt[i,1])
  name.Target <- (Targetnames[i, 1])
  
  peaklist.annotated <- annotate.compound(peaklist, "Target", mz.Target, rt.Target, rtlim = 20, ppmlim = 5)
  is.Target <- rownames(peaklist.annotated)[which(peaklist.annotated$Compound == "Target")]
  
  Target_data <- data.frame(matrix(vector(), 7, 3,
                                   dimnames=list(c( "Extraction_I", "Extraction_II", "Extraction_III",
                                                    "Extraction_III_B", "Extraction_IV", "Extraction_IV_B",
                                                    "QC"),
                                                 c("Average Area", "Stabw", "VarCoeff"))), stringsAsFactors=F)
  
  
  
  
  
  
  
  Target.Boxplot <- as.data.frame(matrix(ncol = 3, nrow = 37, dimnames = list(NULL, c("area", "extraction", "experiment"))))
  PeakTarget.Boxplot <- peaklist.annotated[which(peaklist.annotated$Compound == "Target"),(8 + length(class.levels)):length(peaklist)]
  
  if(i == 18 && column == "PhenHex" && polarity == "positive"){ ##needed due to Arginine on PH pos being complicated
    Target.Boxplot$area <- as.data.frame(t(PeakTarget.Boxplot[2, ]))
  }else{
    Target.Boxplot$area <- as.data.frame(t(PeakTarget.Boxplot[1, ]))
  }
  colnames(Target.Boxplot[,1]) <- "area"
  Target.Boxplot$extraction <- set$class
  Boxplot_order <- c("Extraction_I", "Extraction_II", "Extraction_III", "Extraction_III_B", "Extraction_IV", "Extraction_IV_B", "QC")
  
  extraction.order <- factor(Target.Boxplot$extraction, level=peaknumber_order)
  p.Target <- ggplot(Target.Boxplot, aes(x = extraction.order, y = Target.Boxplot[,1] )) +
    geom_boxplot(aes(group = extraction)) +
    stat_compare_means(method = "anova", label.y = Target.Boxplot[24,1]*1.5)+
    stat_compare_means(label = "p.signif", method = "t.test", ref.group = "Extraction_I", label.y = 0) +
    xlab("Extraction Method") +
    ylab("Area") +
    theme_bw()
  
  
  
  
  
  
  
  ## added step for multiple detection of the same peak
  ## manually select the better integrated peak, which was the first one except Arginine on PhenHex positive
  if(length(is.Target) >0){
      if(name.Target == "Arginine" && column == "PhenHex" && polarity == "positive"){
      Target <- peak.matrix[is.Target[2], ]
    }else{
      Target <- peak.matrix[is.Target[1], ]
    }
    
    
    ## Plots a Boxplot of the peak area of each target and calculates t.test against Extraction_I
    Target.Boxplot <- as.data.frame(matrix(ncol = 3, nrow = nrow(peaknumber),
                                           dimnames = list(NULL, c("area", "extraction", "experiment"))))
    Target.Boxplot$area <- as.data.frame(t(Target))
    
    colnames(Target.Boxplot[,1]) <- "area"
    Target.Boxplot$extraction <- set$class
    Boxplot_order <- c("Extraction_I", "Extraction_II", "Extraction_III", "Extraction_III_B", "Extraction_IV", "Extraction_IV_B", "QC")
    
    extraction.order <- factor(Target.Boxplot$extraction, level=peaknumber_order)
    p.Target <- ggplot(Target.Boxplot, aes(x = extraction.order, y = Target.Boxplot[,1] )) +
      geom_boxplot(aes(group = extraction)) +
      stat_compare_means(method = "anova", size = 3.5)+
      stat_compare_means(label = "p.signif", method = "t.test", ref.group = "Extraction_I", size = 6) +
      xlab("Extraction Method") +
      ylab("Peak area") +
      theme_bw(base_size = 16)
    
    
    ## adds the results of each target to the complete_data object
    complete_data[i, ] <- Target
    
    Targetname <- Targetnames[i,1]
    Target_header <- paste(Targetname, "pdf", sep = ".")
    Target_header <- paste(Target_Folder, "/", Target_header, sep = "")
    pdf(Target_header, width = 12, height = 8)
    
    ## Plots the EIC of each target and the Boxplot in a separate PDF datafile
    plot(p.Target)
    groupid <- c(is.Target)
    groupid <- groupid[groupid != 0]
    eic.Target <- getEIC(set, groupidx = groupid, rt = "corrected")
    plot(eic.Target, set, groupidx = groupnames(eic.Target))
    
    
    dev.off()
   
    
  }
}

colnames(complete_data) <- set$class
complete_data <- as.data.frame(t(complete_data))
colnames(complete_data) <- Targetnames$Name
complete_data_excel <- cbind(" "=rownames(complete_data), complete_data)

write_xlsx(complete_data_excel, "Target_peakareas.xlsx")
