rm(list = ls())

# libraries   

library(Hmisc)
library(argparser)
library(dplyr)
library(ggmap)
library(rjson)

# arguments
p <- arg_parser("This program generates the Map with sampling locations. Each marker colour corresponds to a sample set, Red = Transect 2014, blue = Coastal 2015, black = LMO 2013-2014. The marker size indicates the salinity of the water sample, ranging from 2.44 - 28.05%. All samples were collected from surface waters.")
p <- add_argument(p, "-D", help="input file, metadata", default="BS_MAGv2_sample_metadata.tsv")

argv <- parse_args(p)


# Reading Metadata
df <- read.csv(argv$D, sep="\t", header = FALSE) # Reading sample metadata
df2 <- as.data.frame(t(df))  # row to columns
names(df2)<-df$V1 # Column names
d <- df2[-1,] # removing unnecessary lines (Row of names twice)
# converting to character instead of factors
d$samples <- as.character(d$samples)
var <- c("Sal", "Depth", "O2", "Temp", "NH4", "NO3", "PO4", "Chla", "DOC", "Silicate", "Lat", "Lon")
for (inp in var) {
  d[[inp]]<-as.numeric(as.character(d[[inp]]))
}
#selecting relevant data

finaldf <- select(d, samples, Sal, Depth, O2, Temp, NH4, NO3, PO4, Chla, DOC, Silicate, Lat, Lon)

#LMO <- finaldf[startsWith(as.character(finaldf$samples), 'P4201'),]
#Costal <- finaldf[startsWith(as.character(finaldf$samples), 'P6071'),]
#Ask <- finaldf[startsWith(as.character(finaldf$samples), 'SR'),]
#Redoxcline <- finaldf[startsWith(as.character(finaldf$samples), 'P2'),]
#Transect <- finaldf[startsWith(as.character(finaldf$samples), 'P1'),]

# selected samples
transect_samples <- c('P1994_122', 'P1994_125', 'P1994_107', 'P1994_119', 'P1994_104', 'P1994_110', 'P1994_128', 'P1994_116', 'P1994_101', 'P1994_113')
lmo_samples <- c('P4201_102', 'P4201_110', 'P4201_116', 'P4201_121', 'P4201_124', 'P4201_112', 'P4201_103', 'P4201_120', 'P4201_123', 'P4201_108', 'P4201_101', 'P4201_122', 'P4201_118', 'P4201_106', 'P4201_109', 'P4201_107', 'P4201_105', 'P4201_119', 'P4201_111', 'P4201_104', 'P4201_114')
costal_samples <- c('P6071_523', 'P6071_534', 'P6071_515', 'P6071_503', 'P6071_529', 'P6071_510', 'P6071_513', 'P6071_530', 'P6071_508', 'P6071_527', 'P6071_502', 'P6071_528', 'P6071_531', 'P6071_506', 'P6071_501', 'P6071_518', 'P6071_504', 'P6071_532', 'P6071_533', 'P6071_509', 'P6071_517', 'P6071_507', 'P6071_521', 'P6071_524', 'P6071_522', 'P6071_520', 'P6071_519', 'P6071_505', 'P6071_525', 'P6071_511', 'P6071_512', 'P6071_514', 'P6071_516', 'P6071_526')

# Dataset with selected samples
cost <- finaldf %>% filter(samples %in% costal_samples )
cost$Dataset <- rep("Costal", length(cost$samples))
t <- finaldf %>% filter(samples %in% transect_samples )
t$Dataset <- rep("Transect", length(t$samples))
lmo <- finaldf %>% filter(samples %in% lmo_samples )
lmo$Dataset <- rep("LMO", length(lmo$samples))


all <-rbind(cost, lmo, t)

sbbox <- make_bbox(lon = c(6, 31), lat = c(53, 66), f = .05)
# get map
baltic = get_map(location=sbbox, zoom=6, maptype="terrain")
# create map
balticmap = ggmap(baltic)

# display map
N <- balticmap + geom_point(data = all, mapping = aes(x = Lon, y = Lat, fill=Dataset, color=Dataset, size=Sal), #shape=Dataset), 
            alpha = .4)+
  scale_color_manual(values=c("red", "blue", "black"))+
  theme(legend.title = element_text(colour="black", size=10, face="bold"))+
  theme(legend.text = element_text(colour="black", size=10, face="bold"))+
  scale_size(name = "Salinity",  breaks = c(3,8,14,22,28))+
  scale_x_continuous(name="Longitude", limits=c(6, 31)) +
  scale_y_continuous(name="Latitude", limits=c(53, 66)) 

ggsave("input_pogenom_validation_samples.pdf", N, width = 20, height = 20, units = "cm")  


# MAX, MEDIAN, MEAN; MIN of environmental factors in dataset

for (f in seq(2,length(var)+1,1)){
  J <-na.omit(all[[f]])
  cat(paste(var[f-1],round(max(J), 3), round(median(J), 3), round(mean(J), 3), round(min(J),3), " ", sep="\n"))
}

