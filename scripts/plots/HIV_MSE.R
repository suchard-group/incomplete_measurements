library(ggplot2)
library(stringr)
library(wesanderson)
source("custom_boxplot.R")
options(device="windows")


mse.pathspvl <- "hiv_prediction.csv"
mse.dataspvl <- read.csv(mse.pathspvl)
mse.outfilespvl <- "HIV_MSE_SPVL.pdf"

# mse.data$wrappedProcess = str_wrap(mse.data$Process, width=12)

make_mse_plot <- function(data, out_file){
  
  my_theme <- function(){
    theme_bw(base_size=12) %+replace%
      theme(
        
      )
  }
  
  quantiles_95 <- function(x) {
    r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }
    
  p <- ggplot(data, aes(x = Process, y = MSE, fill = Dimension)) +
        stat_boxplot_custom(qs = c(0.025, 0.25, 0.5, 0.75, 0.975), outlier.shape=NA) +
        ylab("log MSE") +
        # geom_boxplot(outlier.size=-1) +
        # stat_summary(fun.data = quantiles_95, geom="boxplot") +
        # scale_y_continuous(trans='log2') + 
        scale_fill_manual(values=wes_palette(n=2, name="Darjeeling2")) +
        my_theme() + 
        theme(panel.grid.major.x = element_blank(),
              axis.title.x = element_blank(),
              legend.title = element_blank(),
              axis.text.x = element_text(hjust = 0.5))
  
  ggsave(out_file, 
         height = 8, width = 14.7, units = "cm", dpi = 320)
}


mse.dataspvl.sub <- mse.dataspvl
make_mse_plot(mse.dataspvl.sub, mse.outfilespvl)