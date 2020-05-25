x <- 5
y <- 10
number <- x+y

glengths <- c(4.6, 3000, 50000)

species <- c("ecoli", "human", "corn")

combined <- c(glengths, species)

expression <- c("low", "high", "medium", "high", "low", "medium", "high")
expression <- factor(expression)

samplegroup <- c("CTL", "CTL", "CTL", "KO", "KO", "KO", "OE", "OE", "OE")
samplegroup <- factor(samplegroup)

df <- data.frame(species, glengths)

list1 <- list(species, df, number)
list2 <- list(species,glengths,number)

mean(glengths)

square_it <- function(x) {
  square <- x*x
  return(square)
}
square_it(5)

metadata <- read.csv(file="data/mouse_exp_design.csv")

age <- c(15, 22, 45, 52, 73, 81)
alphabets <- c("C","D","X","L","F")

expression <- factor(expression, levels=c("low", "medium", "high"))

samplegroup <- factor(samplegroup, levels=c("KO","CTL","OE"))

comp2 <- list1[[2]]
random <- list(metadata,age,list1,samplegroup,number)

names(list1) <- c("species","df","number")
names(random) <- c("metadata","age","list1","samplegroup","number")
sub_meta <- metadata[which(metadata$replicate>1),]
write.csv(sub_meta, file="data/subset_meta.csv")
write(glengths, file="data/genome_lengths.txt", ncolumns=1)

rpkm_data <- read.csv("data/counts.rpkm.csv")

A <- c(1,3,5,7,9,11)   # odd numbers
B <- c(2,4,6,8,10,12)  # even numbers
B <- c(2,4,6,8,1,5)  # add some odd numbers in 

A <- c(10,20,30,40,50)
B <- c(50,40,30,20,10)

x <- rownames(metadata)
y <- colnames(rpkm_data)
important_genes <- c("ENSMUSG00000083700", "ENSMUSG00000080990", "ENSMUSG00000065619", "ENSMUSG00000047945", "ENSMUSG00000081010", "ENSMUSG00000030970")

rpkm_data[rownames(rpkm_data) %in% important_genes,]
rpkm_data[important_genes,]

teaching_team <- c("Jihe", "Mary", "Meeta", "Radhika")
reorder_teach <- teaching_team[c(4, 2, 1, 3)]

first <- c("A","B","C","D","E")
second <- c("B","D","E","A","C")
second_reordered <- second[match(first,second)]
genomic_idx <- match(rownames(metadata), colnames(rpkm_data))
rpkm_ordered  <- rpkm_data[,genomic_idx]

library(tidyverse)
library(purrr)
samplemeans <- map_dbl(rpkm_ordered, mean)

list_purrr <- list(c(0:10), c(20:30), c(40:50))
list_purrr
map(list_purrr, median)
map_dbl(list_purrr, median)
map_chr(list_purrr, median)

age_in_days <- c(40, 32, 38, 35, 41, 32, 34, 26, 28, 28, 30, 32)
new_metadata <- data.frame(metadata, samplemeans, age_in_days) 

library(ggplot2)
ggplot(new_metadata) +
  geom_point(aes(x=age_in_days, y=samplemeans, color=genotype, shape=celltype), size=3.0) +
  theme_bw() +
  theme(axis.title=element_text(size=rel(1.5)), plot.title=element_text(size=rel(1.5), hjust=0.5)) +
  xlab("Age (days)") + ylab("Mean expression") +
  ggtitle("Sample plot")

ggplot(new_metadata) +
  geom_boxplot(aes(x = genotype, y=samplemeans, fill = celltype)) +
  ggtitle("Genotype differences in average gene expression") +
  xlab("Genotype") +
  ylab("Mean expression") +
  theme(axis.title=element_text(size=rel(1.5))) +
  theme(plot.title=element_text(size=rel(1.5))) +
  scale_fill_manual(values=c("purple","orange")) 
  
## Open device for writing
pdf("figures/scatterplot.pdf")
## Make a plot which will be written to the open device, in this case the temp file created by pdf()/png()
ggplot(new_metadata) +
  geom_point(aes(x = age_in_days, y= samplemeans, color = genotype,
                 shape=celltype), size=rel(3.0)) 
## Closing the device is essential to save the temporary file created by pdf()/png()
dev.off()

hist(new_metadata$samplemeans, xlab="Mean expression level", main="", col="darkgrey", border=FALSE) 

rep_number <- metadata$replicate
head(factor(rep_number))
factor(rep_number) %>% head()

functional_GO_results <- read_delim(file = "data/gprofiler_results_Mov10oe.csv", delim = "\t" )

bp_oe <- functional_GO_results %>%
  filter(domain == "BP")

bp_oe <- bp_oe %>%
  filter(relative.depth > 4)
bp_oe <- bp_oe %>%
  select(term.id, term.name, p.value, query.size, term.size, overlap.size, intersection)
bp_oe <- bp_oe %>%
  arrange(p.value)
bp_oe <- bp_oe %>% 
  dplyr::rename(GO_id = term.id, 
                GO_term = term.name)
bp_oe <- bp_oe %>% 
  dplyr::rename(genes = intersection)

bp_oe <- bp_oe %>%
  mutate(gene_ratio = overlap.size / query.size)

bp_oe <- bp_oe %>%
  mutate(term_percent = overlap.size / term.size)
bp_plot <- bp_oe[1:30, ]

library(RColorBrewer)
# Testing the palette with three colors
display.brewer.pal(3, "YlOrRd")

# Define a palette
mypalette <- brewer.pal(3, "YlOrRd")

# how are the colors represented in the mypalette vector?
mypalette

library(ggplot2)
ggplot(bp_plot) +
  geom_point(aes(x = gene_ratio, y = GO_term, color = -log10(p.value)), 
             size = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(size=rel(1.15)),
        axis.title = element_text(size=rel(1.15))) +
  xlab("Gene ratios") +
  ylab("Top 30 significant GO terms") +
  ggtitle("Dotplot of top 30 significant GO terms") +
  theme(plot.title = element_text(hjust=0.5, 
                                  face = "bold")) +
  scale_color_gradientn(colors = mypalette)

ggplot(bp_plot) +
  geom_point(aes(x = gene_ratio, y = GO_term, color = -log10(p.value)), 
             size = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(size=rel(1.15)),
        axis.title = element_text(size=rel(1.15))) +
  xlab("Gene ratios") +
  ylab("Top 30 significant GO terms") +
  ggtitle("Dotplot of top 30 significant GO terms") +
  theme(plot.title = element_text(hjust=0.5, 
                                  face = "bold")) +
  scale_color_gradientn(name = "Significance \n (-log10(padj))", colors = mypalette) +
  theme(legend.title = element_text(size=rel(1.15),
                                    hjust=0.5, 
                                    face="bold"))


ggplot(bp_plot) +
  geom_col(aes(x = GO_term, y = overlap.size),
           fill = "royalblue",
           color = "black") +
  theme(axis.text.x = element_text(size=rel(1.15)),
        axis.title = element_text(size=rel(1.15))) +
  theme(plot.title = element_text(hjust=0.5, 
                                  face = "bold")) +
  labs(title = "DE genes per GO process", x = NULL, y =  "# DE genes") +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1)) + 
  theme(plot.margin = unit(c(1,1,1,3), "cm")) +
  scale_y_continuous(expand = c(0, 0)) 


y <- c(1:10)
z <- rep(y, length(y))

##Based on the number of genes associated with each GO term (“term.size” column) we can categorize them into “small”, “large” or “medium” categories. Once we have done that, we want to determine what the spread of p-values is for each category; we can do this by drawing a boxplot.
#Use the following code to create a new column in bp_oe tibble for the new categories
x <- bp_oe$term.size
sizes <- rep(NA, length(x) )
sizes[which(x > 3000)] <- "large"
sizes[which(x <= 3000 & x > 500 )] <- "medium"
sizes[which(x < 500)] <- "small"
bp_oe$term_cat <- factor(sizes, levels = c("small","medium","large"))
