library("wesanderson")

violin_plot <- function(data, my_title, my_ylab = "log10(Intensities)", my_filename = my_title, color_palette = wes_palette("Darjeeling1", 40, type = "continuous") ){
  
  violin_plot <- melt(log10(data))
  violin_plot$fraction <- as.factor(sort(unlist(sapply(as.character(violin_plot$variable), function(x) substr(x, 1, nchar(x)-1)))))
  if(length(which(is.na(violin_plot$value))) > 0){
    violin_plot <- violin_plot[-is.na(violin_plot$value), ]
  }
  
  
  ggplot(violin_plot, aes(x=variable, y = value , fill = as.numeric(fraction))) + geom_violin(show.legend = FALSE, na.rm = TRUE) + theme_classic(base_size = 18) + xlab("") + ylab(my_ylab) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_gradientn(colours = color_palette) + ggtitle(my_title) + stat_summary(fun =median, geom="point", show.legend = FALSE, na.rm = TRUE)
  
  ggsave(paste("./Plots/violin_", my_filename,".jpg", sep =""), device = "jpg", width = 16 , height = 14, units = "cm")
  
}
