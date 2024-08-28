library(ggplot2)
library("wesanderson")
pal <- wes_palette("Darjeeling1", 40, type = "continuous")



barplot_numbers_ggplot <- function(conditions, number, file_location, ylab = "", geom_text = TRUE, y_lim, rotate_axis = 0){
  
  df <- data.frame(conditions = conditions, number = number)
  df$conditions <- factor(df$conditions, levels = df$conditions)
  
  if(geom_text == TRUE){
    ggplot(df, aes(x = conditions, y = number, fill = c(1:length(conditions)))) + geom_bar( stat="identity", show.legend = FALSE) +
      scale_fill_gradientn(colours = pal) + theme_classic(base_size = 26) + xlab("") + ylab (ylab)  +  
      geom_text(aes(label = round(number, 1)), vjust = -0.5, size = 8) + ylim(0, y_lim) +  theme(axis.text.x = element_text(angle = rotate_axis))
  }else{
    ggplot(df, aes(x = conditions, y = number, fill = c(1:length(conditions)))) + geom_bar( stat="identity", show.legend = FALSE) +
      scale_fill_gradientn(colours = pal) + theme_classic(base_size = 26) + xlab("") + ylab (ylab) + ylim(0,1)
  }
  
  
  ggsave(paste(file_location, ".jpg", sep =""), device = "jpg", width = 16 , height = 14, units = "cm")
  
}
  
