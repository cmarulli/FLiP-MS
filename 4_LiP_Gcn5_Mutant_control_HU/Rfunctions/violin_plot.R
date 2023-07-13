violin_plot <- function(data, my_title, my_ylab = "log10(Intensities)", my_filename = my_title ){
  violin_plot <- melt(log10(data))
  violin_plot$fraction <- as.factor(sort(unlist(sapply(as.character(violin_plot$variable), function(x) substr(x, 1, nchar(x)-1)))))
  
  ggplot(violin_plot, aes(x=variable, y = value , fill = as.numeric(fraction))) + geom_violin(show.legend = FALSE) + theme_classic(base_size = 18) + xlab("") + ylab(my_ylab) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_gradientn(colours = pal) + ggtitle(my_title) + stat_summary(fun.y=median, geom="point", show.legend = FALSE)
  
  ggsave(paste("./Plots/violin_", my_filename,".jpg", sep =""), device = "jpg", width = 16 , height = 14, units = "cm")
  
}
