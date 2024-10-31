
library(ggplot2)
library(grid)
library(gridExtra)
library(gridtext)
library(gtable)
library(grDevices)


Linkage_map_file_1 <- "Trit_Y_map_1_mapped.txt"
Linkage_map_file_2 <- "Lissotriton_Y_map_7_mapped_SA_reorder.txt"
Blast_result_file_1 <- 'Trit_RADmap_blast_1_All_filtered.txt'
Blast_result_file_2 <- "Liss_lowMAF_blast_filtered_with_Y_translated.txt"
Genome_struct_file <- "Pluro_genome_struct_chr.txt"
output_name <- output_string

# Import the linkage map, filtered blast results and genome structure

Linkage_map_table_1 <- read.table(Linkage_map_file_1)
Linkage_map_table_2 <- read.table(Linkage_map_file_2)
Blast_result_table_1 <- read.table(Blast_result_file_1)
Blast_result_table_2 <- read.table(Blast_result_file_2)
Genome_struct_table <- read.table(Genome_struct_file)

# Get the start positions of the B chromosome segments

Genome_struct_table[,"start_pos"] <- 0

for(i in 1:nrow(Genome_struct_table)) {
  if( Genome_struct_table[i,4] == 2) {Genome_struct_table[i,5] <- Genome_struct_table[i - 1,2] }
}

Genome_struct_table[,"extra_length"] <- 0

for(i in 2:nrow(Genome_struct_table)) {
  if( Genome_struct_table[i,4] == 2) {Genome_struct_table[i - 1,"extra_length"] <- Genome_struct_table[i,2] }
}

Genome_struct_table[,"total_length"] <- Genome_struct_table$V2 + Genome_struct_table$extra_length

Chr_lengths <- subset(Genome_struct_table, V4 == 1, select = c("V3","total_length"))

row.names(Chr_lengths) <- Chr_lengths$V3
Chr_lengths[,"Start_pos"] <- 0

for(chr in 2:nrow(Chr_lengths)) {
  
  Chr_lengths[chr,"Start_pos"] <- Chr_lengths[chr - 1,2] + Chr_lengths[chr - 1,3]
}

Chr_lengths[,"End_pos"] <- Chr_lengths[,"Start_pos"] + Chr_lengths[,"total_length"]
Chr_lengths[,"Mid_point"] <- Chr_lengths$Start_pos + (0.5 * Chr_lengths$total_length)

# Reconstruct blast coordinates then add blast data to map

Add_blast_coords <- function(linkage_map_in, blast_table_in) {
  
  Blast_result_table <- blast_table_in
  Linkage_map_table <- linkage_map_in
  
  Blast_result_table[,"Raw_scaff"] <- match(Blast_result_table$V2,Genome_struct_table$V1)
  
  for(i in 1:nrow(Blast_result_table)) {
    
    raw_scaff_no <- Blast_result_table[i,"Raw_scaff"]
    
    Blast_result_table[i,9] <- Blast_result_table[i,9] + Genome_struct_table[raw_scaff_no,5]
    Blast_result_table[i,10] <- Blast_result_table[i,10] + Genome_struct_table[raw_scaff_no,5]
    Blast_result_table[i,2] <- Genome_struct_table[raw_scaff_no,3]
  }
  
  Blast_result_table[,"Map_no"] <- match(Blast_result_table$V1,Linkage_map_table$V1)
  
  Linkage_map_table[,"Genome_Chr"] <- NA
  Linkage_map_table[,"Chr_pos"] <- NA
  Linkage_map_table[,"Genome_pos"] <- NA
  
  for(i in 1:nrow(Blast_result_table)) {
    
    Map_no <- Blast_result_table[i,"Map_no"]
    
    if(Map_no %in% row.names(Linkage_map_table)){
      Linkage_map_table[Map_no,"Genome_Chr"] <- Blast_result_table[i,2]
      Linkage_map_table[Map_no,"Chr_pos"] <- Blast_result_table[i,9]
      Linkage_map_table[Map_no,"Genome_pos"] <- Linkage_map_table[Map_no,"Chr_pos"] + Chr_lengths[Linkage_map_table[Map_no,"Genome_Chr"],"Start_pos"]
    }
  }
  
  return(Linkage_map_table)
}

Liss_linkage_with_blast <- Add_blast_coords(Linkage_map_table_2, Blast_result_table_2)

# Generate stats with chromosomes

match_chromosomes <- function(Linkage_map_table){
  
  nGroups <- max(Linkage_map_table$V3)
  
  stats_table <- data.frame(matrix(NA, nrow = nGroups, ncol = nrow(Chr_lengths)))
  colnames(stats_table) <- Chr_lengths$V3
  
  Hits_table <- na.omit(Linkage_map_table)
  
  for(i in 1:nGroups){
    
    group_table <- subset(Hits_table, V3==i)
    
    for(j in 1:nrow(Chr_lengths)){
      
      stats_table[i,j] <- sum(group_table$Genome_Chr == Chr_lengths$V3[j], na.rm=T)
    }
  }
  
  stats_table[,"Total_hits"] <- rowSums(stats_table)
  stats_table[,"Best_match"] <- NA
  stats_table[,"recipricol"] <- NA
  
  for(i in 1:nGroups){
    
    best_match_index <- which.max(stats_table[i,1:nrow(Chr_lengths)])
    stats_table[i,"Best_match"] <- colnames(stats_table)[best_match_index]
    stats_table[i,"recipricol"] <- which.max(stats_table[,best_match_index])
    stats_table[i,"Best_match_hits"] <- stats_table[i,best_match_index]
  }
  
  stats_table[,"correlation"] <- NA
  stats_table[,"inverted"] <- FALSE
  stats_table[,"length"] <- FALSE
  
  for(i in 1:nGroups){
    
    group_table <- subset(Hits_table, V3==i & Genome_Chr == stats_table[i,"Best_match"])
    
    coefficient <- cor.test(group_table$V2, group_table$Chr_pos, method = "spearman")
    stats_table[i,"correlation"] <- coefficient$estimate
    
    if(coefficient$estimate < 0) { stats_table[i,"inverted"] <- TRUE}
    
    group_end <- findInterval(i,Linkage_map_table$V3)
    stats_table[i,"length"] <- as.numeric(Linkage_map_table[group_end,2])
  }
  
  ordered_stats <- stats_table[order(as.numeric(stats_table$Best_match)),]
  
  ordered_stats[,"Start"] <- 0
  
  for(group in 2:nGroups) {
    
    ordered_stats[group,"Start"] <- ordered_stats[group - 1,"length"] + ordered_stats[group - 1,"Start"]
      Chr_lengths[chr - 1,2] + Chr_lengths[chr - 1,3]
  }

  ordered_stats[,"End"] <- ordered_stats[,"Start"] + ordered_stats[,"length"]
  ordered_stats[,"Mid"] <- ordered_stats[,"Start"] + ordered_stats[,"length"] * 0.5

  return(ordered_stats)
}

Liss_stats_2 <- match_chromosomes(Liss_linkage_with_blast)

# Re-order map

order_map <- function(Linkage_map_table, stats_table){
  
  Hits_table <- na.omit(Linkage_map_table)
  Hits_table[,"map_pos"] <- NA
  
  print(Hits_table[1,])
  
  output_table <- Hits_table[0,]
  
  nGroups <- max(Hits_table$V3)
  
  for(group in 1:nGroups){
    group_table <- subset(Hits_table, V3==group)
    
    if(stats_table[as.character(group),"inverted"] == TRUE){
      group_table$V2 <- (group_table$V2 * -1) + stats_table[as.character(group),"length"]
    }
    group_table$map_pos <- group_table$V2 + stats_table[as.character(group),"Start"]
    print(group_table[1,])
    
    output_table <- rbind(output_table,group_table)

  }
  return(output_table)
}

Liss_map_ordered <- order_map(Liss_linkage_with_blast, Liss_stats_2)

# Make table with map characteristics

Liss_map_groups <- Liss_stats_2[,c("recipricol","Start","End","Mid")]

Trit_linkage_with_blast <- Add_blast_coords(Linkage_map_table_1, Blast_result_table_1)
Trit_stats_2 <- match_chromosomes(Trit_linkage_with_blast)
Trit_map_ordered <- order_map(Trit_linkage_with_blast, Trit_stats_2)

Trit_map_ordered$Colour <- "black"
Trit_map_ordered$Size <- 1

for(i in 1:nrow(Trit_map_ordered)){
  if(grepl("Y",Trit_map_ordered[i,"V1"])==TRUE) {Trit_map_ordered[i,"Colour"] <- "red"}
  if(grepl("Y",Trit_map_ordered[i,"V1"])==TRUE) {Trit_map_ordered[i,"Size"] <- 3}
}

Liss_map_ordered$Colour <- "black"
Liss_map_ordered$Size <- 1

for(i in 1:nrow(Liss_map_ordered)){
  if(grepl("Y",Liss_map_ordered[i,"V1"])==TRUE) {Liss_map_ordered[i,"Colour"] <- "blue"}
  if(grepl("Y",Liss_map_ordered[i,"V1"])==TRUE) {Liss_map_ordered[i,"Size"] <- 3}
}



Trit_map_groups <- Trit_stats_2[,c("recipricol","Start","End","Mid")]

# Make this plot:

Pick_scale <- function(total_length,fraction = 0.25) {
  optimal_scale_length <- total_length * fraction
  scale_mag <- 10^(floor(log10(total_length)) - 1)
  Scale_factors <- c(1,2,5,10)
  Scale_options <- Scale_factors * scale_mag
  Scale_results <- abs(Scale_options - optimal_scale_length)
  scale_best <- Scale_options[which.min(Scale_results)]
  return(scale_best)
}

X_text_1 <- "_Triturus ivanbureschi_ Linkage Group"
X_text_2 <- "_Lissotriton vulgaris_ Linkage Group"
X_length_1 <- max(Trit_map_ordered$map_pos)
X_length_2 <- max(Liss_map_ordered$map_pos)
X_axis_offset <- -2.0e+09
X_scale <- Pick_scale(X_length)
X_scale_text <- paste(X_scale,"cM")
x_text_offset <- 0.75
x_midpoint_1 <- X_length_1 / 2
x_midpoint_2 <- X_length_2 / 2

y_text <- "_Pleurodeles waltl_ Chromosome"
y_axis_offset <- -60
y_scale <- 5e+09
y_scale_text <- "5 Gbp"
y_text_offset <- 0.75
y_midpoint <- max(Trit_map_ordered$Genome_pos) / 2

axis_par <- gpar(fontsize = 10)

Trit_plot <- ggplot(Trit_map_ordered, aes(x=map_pos, y=Genome_pos)) + geom_point(colour = Trit_map_ordered$Colour, size = Trit_map_ordered$Size) +
  scale_x_continuous(breaks = Trit_map_groups$Mid, minor_breaks = Trit_map_groups$Start, expand = c(0,0), labels=Trit_map_groups$recipricol) +
  scale_y_continuous(breaks = Chr_lengths$Mid_point, minor_breaks = Chr_lengths$Start_pos, expand = c(0,0), labels=Chr_lengths$V3, limits = c(0,sum(Chr_lengths$total_length))) +
  theme(panel.grid.minor = element_line(colour="black", linewidth=0.5), panel.grid.major = element_blank(), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA)) +
  theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank()) +
  
  annotation_custom(grid::linesGrob(gp = gpar(lwd = unit(2, "pt"))), xmin = 0, xmax = X_scale, ymin = X_axis_offset, ymax = X_axis_offset)  +
  annotation_custom(richtext_grob(gp = axis_par, X_scale_text, vjust = -0.5), xmin = 0, xmax = X_scale, ymin = X_axis_offset, ymax = X_axis_offset) +
  annotation_custom(richtext_grob(gp = axis_par, X_text_1, vjust = -0.5), xmin = x_midpoint_1, xmax = x_midpoint, ymin = X_axis_offset, ymax = X_axis_offset) +
  
  annotation_custom(grid::linesGrob(gp = gpar(lwd = unit(2, "pt"))), xmin = y_axis_offset, xmax = y_axis_offset, ymin = 0, ymax = y_scale)  +
  annotation_custom(richtext_grob(gp = axis_par, y_scale_text, rot = 90, vjust = -0.5), xmin = y_axis_offset, xmax = y_axis_offset, ymin = 0, ymax = y_scale) +
  annotation_custom(richtext_grob(gp = axis_par, y_text, rot = 90, vjust = -0.5), xmin = y_axis_offset, xmax = y_axis_offset, ymin = y_midpoint, ymax = y_midpoint) +
  coord_cartesian(clip = 'off') + theme(plot.margin=unit(c(0.1,0.1,0.4,0.35), 'in'))
  
  

Liss_plot <- ggplot(Liss_map_ordered, aes(x=map_pos, y=Genome_pos)) + geom_point(colour = Liss_map_ordered$Colour, size = Liss_map_ordered$Size) +
  scale_x_continuous(breaks = Liss_map_groups$Mid, minor_breaks = Liss_map_groups$Start, expand = c(0,0), labels=Liss_map_groups$recipricol) +
  scale_y_continuous(breaks = NULL, minor_breaks = Chr_lengths$Start_pos, expand = c(0,0), labels = NULL, name = NULL, limits = c(0,sum(Chr_lengths$total_length))) +
  theme(panel.grid.minor = element_line(colour="black", linewidth=0.5), panel.grid.major = element_blank(), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA)) +
  theme(axis.title.x = element_blank()) +
  
  annotation_custom(grid::linesGrob(gp = gpar(lwd = unit(2, "pt"))), xmin = 0, xmax = X_scale, ymin = X_axis_offset, ymax = X_axis_offset)  +
  annotation_custom(richtext_grob(gp = axis_par, X_scale_text, vjust = -0.5), xmin = 0, xmax = X_scale, ymin = X_axis_offset, ymax = X_axis_offset) +
  annotation_custom(richtext_grob(gp = axis_par, X_text_2, vjust = -0.5), xmin = x_midpoint_2, xmax = x_midpoint, ymin = X_axis_offset, ymax = X_axis_offset) +
  coord_cartesian(clip = 'off') + theme(plot.margin=unit(c(0.1,0.1,0.4,0.1), 'in')) 

Trit_Y_chr_len_rel <- Chr_lengths$total_length[2] / sum(Chr_lengths$total_length)
Trit_Y_chr_mid_rel <- Chr_lengths$Mid_point[2] / sum(Chr_lengths$total_length)

Liss_Y_chr_len_rel <- Chr_lengths$total_length[5] / sum(Chr_lengths$total_length)
Liss_Y_chr_mid_rel <- Chr_lengths$Mid_point[5] / sum(Chr_lengths$total_length)

Trit_grob <- ggplotGrob(Trit_plot)
Liss_grob <- ggplotGrob(Liss_plot)
combined_grob <- cbind(Trit_grob,Liss_grob)
combined_grob <- gtable_add_grob(combined_grob,rectGrob(gp=gpar(fill="red", alpha = 0.2), height = Trit_Y_chr_len_rel, y = Trit_Y_chr_mid_rel), t = 3, b = 9, l = 7, r = 24)
combined_grob <- gtable_add_grob(combined_grob,rectGrob(gp=gpar(fill="blue", alpha = 0.2), height = Liss_Y_chr_len_rel, y = Liss_Y_chr_mid_rel), t = 3, b = 9, l = 7, r = 24, name = "bluey")
grid.newpage()
grid.draw(combined_grob)

ggplot_build()

Trit_plot_grob3 <- GGplotGrob(Liss_plot)

Trit_plot_grob3 <- ggplot_gtable(ggplot_build(Trit_plot))

Trit_plot_grob3 <- ggplot_gtable(ggplot_build(Liss_plot))

grid.arrange(Trit_plot, Liss_plot, nrow = 1)
# bodge the two together


# Get X co-ordinates

for(i in 1:nrow(Hits_table)) {
  
  marker_chr <- as.numeric(Hits_table[i,"CHR_no"])
  Hits_table[i,"Y_coord"] <- Hits_table[i,5] + Chr_lengths[marker_chr,"Start_pos"]
  if(grepl("Y",Hits_table[i,"V1"])==TRUE) {Hits_table[i,"Colour"] <- "red"}
  if(grepl("Y",Hits_table[i,"V1"])==TRUE) {Hits_table[i,"Size"] <- 3}
  
}

Hits_table <- Hits_table[order(Hits_table$Size),]

# Make plot

Pick_scale <- function(total_length,fraction = 0.25) {
  optimal_scale_length <- total_length * fraction
  scale_mag <- 10^(floor(log10(total_length)) - 1)
  Scale_factors <- c(1,2,5,10)
  Scale_options <- Scale_factors * scale_mag
  Scale_results <- abs(Scale_options - optimal_scale_length)
  scale_best <- Scale_options[which.min(Scale_results)]
  return(scale_best)
}

Group_lenghts[,"Mid_point"] <- Group_lenghts$Start_pos + (0.5 * Group_lenghts$Length)
Chr_lengths[,"Mid_point"] <- Chr_lengths$Start_pos + (0.5 * Chr_lengths$total_length)

Group_font <- 20
Axis_font <- 20

X_text <- "_Lissotriton vulgaris_ Linkage Group"
X_length <- max(Hits_table$X_coord)
X_axis_offset <- -2.0e+09
X_scale <- Pick_scale(X_length)
X_scale_text <- paste(X_scale,"cM")
x_text_offset <- 0.75
x_midpoint <- X_length / 2
  
y_text <- "<i>Pleurodeles waltl</i> Chromosome"
y_axis_offset <- -50
y_scale <- 5e+09
y_scale_text <- "5 Gbp"
y_text_offset <- 0.75
y_midpoint <- max(Hits_table$Y_coord) / 2

#axis_par <- gpar(fontsize = Axis_font, fontfamily = "TT Arial")
axis_par <- gpar(fontsize = Axis_font)

final_plot <- ggplot(Hits_table, aes(x=X_coord, y=Y_coord)) + geom_point(colour = Hits_table$Colour, size = Hits_table$Size) + 
  scale_x_continuous(breaks = Group_lenghts$Mid_point, minor_breaks = Group_lenghts$Start_pos, expand = c(0,0), labels=Group_lenghts$X1.nGroups) +
  scale_y_continuous(breaks = Chr_lengths$Mid_point, minor_breaks = Chr_lengths$Start_pos, expand = c(0,0), labels=Chr_lengths$V3) +
  theme(panel.grid.minor = element_line(colour="black", linewidth=0.5), panel.grid.major = element_blank(), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(color= "black", size = Group_font)) +
  theme(axis.text.y = element_text(color= "black", size = Group_font)) +
  
  annotation_custom(grid::linesGrob(gp = gpar(lwd = unit(2, "pt"))), xmin = 0, xmax = X_scale, ymin = X_axis_offset, ymax = X_axis_offset)  +
  annotation_custom(richtext_grob(gp = axis_par, X_scale_text, vjust = -0.5), xmin = 0, xmax = X_scale, ymin = X_axis_offset, ymax = X_axis_offset) +
  annotation_custom(richtext_grob(gp = axis_par, X_text, vjust = -0.5), xmin = x_midpoint, xmax = x_midpoint, ymin = X_axis_offset, ymax = X_axis_offset) +
  
  annotation_custom(grid::linesGrob(gp = gpar(lwd = unit(2, "pt"))), xmin = y_axis_offset, xmax = y_axis_offset, ymin = 0, ymax = y_scale)  +
  annotation_custom(richtext_grob(gp = axis_par, y_scale_text, rot = 90, vjust = -0.5), xmin = y_axis_offset, xmax = y_axis_offset, ymin = 0, ymax = y_scale) +
  annotation_custom(richtext_grob(gp = axis_par, y_text, rot = 90, vjust = -0.5), xmin = y_axis_offset, xmax = y_axis_offset, ymin = y_midpoint, ymax = y_midpoint) +
  
  coord_cartesian(clip = 'off') +
  theme(plot.margin=unit(c(0.25,0.25,1,1), 'in')) +
  labs(title=" ") +
  theme(plot.title = element_text(size = 25))
  
ggsave(paste(output_name,".pdf", sep = ""), plot = final_plot, width = 14, height = 12, units = "in")

write.table(stats_table, paste(output_name,"_stats.txt", sep = ""), sep = "\t", quote = FALSE)

