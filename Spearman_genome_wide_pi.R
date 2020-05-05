rm(list = ls())

# libraries   

library(Hmisc)
library(argparser)
library(factoextra)
library("cluster")
library(corrplot)
library(dplyr)
library("ggpubr")
library(ggmap)
library(rjson)

# arguments
p <- arg_parser("This program finds potential correlation between genome-wide pi and environmental factors")
p <- add_argument(p, "-D", help="input file, metadata", default="BS_MAGv2_sample_metadata.tsv")
p <- add_argument(p, "-P", help="Main directory with subdirectory containing genomic pi parameter data", default="./POGENOM_OUTPUT")
p <- add_argument(p, "-C", help="cutoff in reduced number of samples (all environmental factors), max.p.adj.value", default=0.05)
p <- add_argument(p, "-T", help="cutoff in all samples (only environmental factors measured in all samples), max.p.adj.value", default=0.0125)
p <- add_argument(p, "-O", help="Output directory", default="./Correlations")
p <- add_argument(p, "-N", help="minimum number of samples in the analysis", default=6)

argv <- parse_args(p)

# Define Functions

Selecting_factors <- function(DB, var_vector, cutoff) {

  lv <- length(var_vector)
  p_values <- rep(0,lv)
  rho <- rep(0, lv)
  selection <- c()
  message <- c()
  best = NULL
  best_factor = NULL
  for ( i in (seq(2,lv+1,1)) ) {
    f <-cor.test(scale(DB$Norm_intra_pi), scale(DB[[i]]), method=c("spearman")) # Spearman correlation Normilsed pi vs environmental factor 
    rho[[i-1]] <- f$estimate
    p_values[[i-1]] <- f$p.value
  }
  p_values_adj <- p.adjust(p_values, method="fdr")

  p <- list()
  r <- list()
  for (v in (seq(1,lv,1))){
    if (p_values_adj[[v]] < cutoff) {  # P.adj.value is below user-defined cutoff
      message <- c(message,paste("      Selected factor:", var_vector[v], "- P.adj.value:", p_values_adj[[v]], "- rho:", rho[[v]],"\n"))
      selection <-c(selection, var_vector[v])  # list of environmental factors with signicant correlation with genome-wide pi
      r[var_vector[v]] <-rho[[v]]
      p[var_vector[v]] <-p_values_adj[[v]]
    }
  }
  if (length(selection) > 1) { # If there are more than one environmental factor significally correlated 
    maximo <- 0
    for (o in seq(1,length(selection),1)) { 
      if (abs(as.numeric(r[o][1])) > maximo) { 
        maximo <- abs(as.numeric(r[o][1]))
        best <- paste("      The best factor:", selection[o][1], "\n", "     rho:", r[o][1], "\n", "     p.adj.value:", p[o][1],"\n", "\n")  
        best_factor <- selection[o][1] # Environmental factor with the highest rho value
      }
    }
  }
  
  return(list(best_factor, best, selection, message))
}

pdf_figures_selected_factors <- function(DB, vec_factors, genome, num_samples) { # printing out correlation figure for significant environmental factors
for (env in vec_factors) {
  dir.create(argv$O)
  directorio <- paste(argv$O, mag, sep = "/")
  dir.create(directorio)
  name2 = paste(env, genome, "Correlation.pdf", sep = "_")
  filename3 = paste(directorio, name2, sep = "/")
  pdf(filename3)
  print(
    ggscatter(DB, x = env, y = "Norm_intra_pi", 
              add = "reg.line", conf.int = TRUE, 
              #        cor.coef = TRUE, cor.method = "spearman",
              xlab = env, ylab = "Norm_pi",  title = paste(genome, "- Num_Samples_analysed:", num_samples))
  )
  dev.off()
  }
}

printing_map <- function(DB, Genome){
  
  transect_samples <- c('P1994_122', 'P1994_125', 'P1994_107', 'P1994_119', 'P1994_104', 'P1994_110', 'P1994_128', 'P1994_116', 'P1994_101', 'P1994_113')
  lmo_samples <- c('P4201_102', 'P4201_110', 'P4201_116', 'P4201_121', 'P4201_124', 'P4201_112', 'P4201_103', 'P4201_120', 'P4201_123', 'P4201_108', 'P4201_101', 'P4201_122', 'P4201_118', 'P4201_106', 'P4201_109', 'P4201_107', 'P4201_105', 'P4201_119', 'P4201_111', 'P4201_104', 'P4201_114')
  costal_samples <- c('P6071_523', 'P6071_534', 'P6071_515', 'P6071_503', 'P6071_529', 'P6071_510', 'P6071_513', 'P6071_530', 'P6071_508', 'P6071_527', 'P6071_502', 'P6071_528', 'P6071_531', 'P6071_506', 'P6071_501', 'P6071_518', 'P6071_504', 'P6071_532', 'P6071_533', 'P6071_509', 'P6071_517', 'P6071_507', 'P6071_521', 'P6071_524', 'P6071_522', 'P6071_520', 'P6071_519', 'P6071_505', 'P6071_525', 'P6071_511', 'P6071_512', 'P6071_514', 'P6071_516', 'P6071_526')
  
  # Dataset with selected samples
  cost <- DB %>% filter(samples %in% costal_samples )
  t <- DB %>% filter(samples %in% transect_samples )
  lmo <- DB %>% filter(samples %in% lmo_samples )
  
  sbbox <- make_bbox(lon = c(5, 32), lat = c(52, 66), f = .05)
  # get map
  brisbane = get_map(location=sbbox, zoom=5,
                     maptype="terrain")
  # create map
  brisbanemap = ggmap(brisbane)
  
  
  directorio <- paste(argv$O, mag, sep = "/")
  dir.create(directorio)
  name4 = paste("MAP", Genome, "samples.pdf", sep = "_")
  filename4 = paste(directorio, name4, sep = "/")
  pdf(filename4)
  
  # display map
  
  print(
    brisbanemap +
      geom_point(data = cost, mapping = aes(x = Lon, y = Lat), 
                 color = "darkblue", size = cost$Sal/2, alpha = .3) +
      geom_point(data = t, mapping = aes(x = Lon, y = Lat), 
                 color = "red", size = t$Sal/2, alpha = .3) +
      geom_point(data = lmo, mapping = aes(x = Lon, y = Lat), 
                 color = "black", size = lmo$Sal/2, alpha = .3) +
      scale_x_continuous(name="Longitude", limits=c(5, 32)) +
      scale_y_continuous(name="Latitude", limits=c(52, 66)) 
  )
  
  dev.off()
  
}

all_samples_correlation<-function(m_factors, m_pi, cutoff1, GENOME){ # Using all samples but only some environmental factors
####USING ALL samples when data is available   
factors_in_all_samples <- c("Sal", "Temp", "PO4", "Chla", "Silicate", "Lat", "Lon") # list of environmental factors measured in all samples
allDPI <- data.frame(m_factors, m_pi)
all_datasN <-data.frame(allDPI$Norm_intra_pi, allDPI$Sal, allDPI$Temp, allDPI$PO4, allDPI$Chla, allDPI$Silicate, allDPI$Lat, allDPI$Lon, check.names = TRUE, check.rows = TRUE, row.names= allDPI$samples)

names(all_datasN) <- c("Norm_intra_pi", "Sal", "Temp", "PO4", "Chla", "Silicate", "Lat", "Lon")
all_best_factor <- Selecting_factors(all_datasN, factors_in_all_samples, cutoff1)

cat("      Running analysis with:\n")
cat("      Environmental factors:", factors_in_all_samples, "\n")
cat("      Total number of samples analysed:", length(allDPI$samples),"- P.adj.value cut-off:", cutoff1,"\n")

if (is.null(all_best_factor[[1]])) {
  cat("      With these parameters, there is not environmental factor that significally correlates with genome-wide pi\n", "\n")
} else {
  cat(all_best_factor[[4]])
  cat(all_best_factor[[2]])
  dir.create(argv$O)
  directorio <- paste(argv$O, mag, sep = "/")
  dir.create(directorio)
  name1 = paste(GENOME, "Selected_factors.txt", sep = "_")
  filename1 = paste(directorio, name1, sep = "/")
  cat("***Running analysis with parameters***\n","Environmental factors:\n", factors_in_all_samples,
                "\nTotal number of samples analysed:",length(allDPI$samples),"\nP adj value cut-off:",cutoff1, "\n", 
                "\n", all_best_factor[[4]], "\n", all_best_factor[[2]], file=filename1)
 
  all_sN <- scale(all_datasN)
  fig_name <- paste(GENOME, "N_samples_analysed", length(allDPI$samples), "scaled", "spearman", "correlations.pdf", sep = "_")
  pdf(paste(directorio,fig_name, sep = "/"))
# printing out correlation matrix - all vs all
    corrplot(cor(all_sN), method="number", col = c("red", "blue"), bg = "white", tl.col= "black", title = GENOME)
 # corrplot(cor(all_datasN), method="number", col = c("red", "blue"), bg = "white", tl.col= "black", title = GENOME)
  
  dev.off()
  pdf_figures_selected_factors(allDPI, all_best_factor[[3]], GENOME, length(allDPI$samples))
  printing_map(allDPI, GENOME)
  }
return(all_best_factor[[1]])
}

analysis_per_file <- function(finaldf, pidata, GENOME, cutoff1, cutoff2) {
  
pi <- read.csv(pidata, sep="\t") # Reading genome-wide pi data per sample
#selecting relevant samples for this MAG
sel <- data.frame() # Only relevant samples 
for (sample in pi$Sample){
  if (sample != "All_samples_combined") {
    sel <- rbind(sel, subset(finaldf, finaldf$samples == sample))
  }
}

#Ordering according to sample name
sel<-sel[order(sel$samples),] 
pi<-pi[order(pi$Sample),]
pi <- pi[!startsWith(as.character(pi$Sample), 'SRX'),] #removing Caspian sea samples

DPI <- na.omit(data.frame(sel, pi[-1,])) #Pi without "all_samples_combined" and removing Na values
DPI <- within(DPI, rm(Sample)) # Keepin only one colum with sample names
DPI <- DPI[-c(19)] # removing Average_read_depth column (there is no variation) 

if (length(DPI$samples) > (argv$N-1)) { # if the number of samples is higher than user-defined threshold in the "reduced dataset" (samples when all environmental factors where measured)
  # correlation Pi vs Env factors - Using samples with data available for all of the environmental factors 
  variables <- names(sel)
  variables <- variables[-1] # removing "samples" from this list
  best_factor <- Selecting_factors(DPI, variables, cutoff1)

    if (is.null(best_factor[[1]])) { # if no correlation was found using the "reduced dataset" (samples when all environmental factors where measured)
      cat("      Environmental factors:", variables, "\n")
      cat("      Total number of samples analysed:", length(DPI$samples),"- P.adj.value cut-off:", cutoff1,"\n")
      cat("      With these parameters, there is not environmental factor that significally correlates with genome-wide pi\n", "\n")
      cat("      Analysis will be run with all samples, excluding NA values in environmental factors","\n")
      ####Trying with ALL samples but less environmental factors 
      validatingA <- all_samples_correlation(sel, pi[-1,], cutoff2, GENOME)
      if (is.null(validatingA)) { validatingB <- all_samples_correlation(sel, pi[-1,], cutoff1, GENOME)}
      if (!is.null(validatingB)){ anytoprint <- FALSE
      } else {anytoprint <- FALSE
              cat("***** There is not environmental factor that significally correlates with genome-wide pi\n", "\n")
              }
  
    } else { # if a correlation was found using the "reduced dataset"
      anytoprint <- TRUE
      factors_in_all_samples <- c("Sal", "Temp", "PO4", "Chla", "Silicate", "Lat", "Lon")
      for (fac in factors_in_all_samples){
        if (best_factor[[1]] == fac) {  # if the environmental factor with the highest significant correlation with pi is one of the environmental factors measured in all samples
          ####USING ALL samples when data is available
          validating <- all_samples_correlation(sel, pi[-1,], cutoff2, GENOME) # trying the lowest p.value cutoff
          if (is.null(validating)) { 
              validating2 <- all_samples_correlation(sel, pi[-1,], cutoff1, GENOME) # if no correlation found, then it tries the highest p.value cutoff
              if (!is.null(validating2)){anytoprint <- FALSE} 
          } else { anytoprint <- FALSE }
          break
        }  
      }
    }
    
    if (anytoprint == TRUE) { # when correlation is significant only using the "reduced dataset"
       cat("      Running analysis with parameters:\n")
       cat("      Environmental factors:", variables, "\n")
       cat("      Total number of samples analysed:", length(DPI$samples),"\n")
       cat("      P adj value cut-off:", cutoff1, "\n")
       cat(best_factor[[4]])
       cat(best_factor[[2]])
      
       dir.create(argv$O)
       directorio <- paste(argv$O, mag, sep = "/")
       dir.create(directorio)
       name2 = paste(GENOME, "Selected_factors.txt", sep = "_")
       filename2 = paste(directorio, name2, sep = "/")
       cat("***Running analysis with parameters***\n","Environmental factors:\n", variables,
           "\nTotal number of samples analysed:",length(DPI$samples),"\nP adj value cut-off:",cutoff1, "\n", 
           "\n", best_factor[[4]], "\n", best_factor[[2]], file=filename2)

        factors <- scale(data.frame(select(DPI, Norm_intra_pi, Sal, Depth, O2, Temp, NH4, NO3, PO4, Chla, DOC, Silicate, Lat, Lon), row.names = DPI$samples))
#       factors <- data.frame(select(DPI, Norm_intra_pi, Sal, Depth, O2, Temp, NH4, NO3, PO4, Chla, DOC, Silicate, Lat, Lon), row.names = DPI$samples)
        # printing out correlation matrix - all vs all - using the "reduced dataset"
        fig_name1 <- paste(GENOME, "N_samples-analysed", length(DPI$samples), "scaled", "spearman", "correlations.pdf", sep = "_")
        pdf(paste(directorio,fig_name1, sep = "/"))
        corrplot(cor(factors), method="number", col = c("red", "blue"), bg = "white", tl.col= "black", title = GENOME)
        dev.off()
        pdf_figures_selected_factors(DPI, best_factor[[3]], GENOME, length(DPI$samples))
        printing_map(DPI, GENOME)
       }
} else { # If the number of samples in the "reduced dataset" (samples when all environmental factors where measured) is too small (below user-defined threshold)
  ####USING ALL samples when data is available
  cat("      warning: Samples with data for all of the environmental factors are:", length(DPI$samples),"\n")
  cat("      Analysis will be run with all samples, excluding environmental factors with NA values","\n")
  validating3 <- all_samples_correlation(sel, pi[-1,], cutoff2, GENOME) # trying the lowest p.value cutoff
  if (is.null(validating3)) { validating4 <- all_samples_correlation(sel, pi[-1,], cutoff1, GENOME) # if no correlation found, then it tries the highest p.value cutoff
      if (is.null(validating4)) { cat("***** There is not environmental factor that significally correlates with genome-wide pi\n", "\n")}
      }
  }  

}

# Reading Metadata

df <- read.csv(argv$D, sep="\t", header = FALSE) # Reading sample metadata
df2 <- as.data.frame(t(df))  # row to columns
names(df2)<-df$V1 # Column names
d <- df2[-1,] # removing unnecessary lines (Row of names twice)

# converting to character instead of factors
d$samples <- as.character(d$samples)

# converting to numeric instead of factors
for (inp in c("Sal", "Depth", "O2", "Temp", "NH4", "NO3", "PO4", "Chla", "DOC", "Silicate", "Lat", "Lon")) {
  d[[inp]]<-as.numeric(as.character(d[[inp]]))
}
#selecting relevant data: name of the samples, and Environmental factor measurments 
finaldf <- select(d, samples, Sal, Depth, O2, Temp, NH4, NO3, PO4, Chla, DOC, Silicate, Lat, Lon)

# List of Genomes
MAGs <- list.files(path = argv$P, pattern = "bin", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

# Correlation for each Genome
for (mag in MAGs) {
  prefix = paste("Pogenom", mag, sep = "_")
  file_name = paste(prefix, "intradiv.txt", sep = ".")
  file = paste(argv$P, mag, file_name, sep = "/")
  cat("INFO: Running analysis on Genome", mag, "\n")
  analysis_per_file(finaldf, file, mag, argv$C, argv$T)
}


