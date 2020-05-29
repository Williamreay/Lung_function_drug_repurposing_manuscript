###############################
## LCV forest plot
## William Reay -  June 2020
###############################

## Load dependencies

library(readxl)
library(easyGgplot2)

LCV_FP_input <- read_excel("Desktop/Pneumonia_cytokine_lung_function/LCV/LCV_FP_input.xlsx")

## Calculate pm.GCP +/- its SE

LCV_FP_input$Min_GCP <- LCV_FP_input$GCP - LCV_FP_input$SE
LCV_FP_input$Max_GCP <- LCV_FP_input$GCP + LCV_FP_input$SE

## Forest plot code

FP <- ggplot(data = LCV_FP_input, aes(x=Metabolite_hormone, y=GCP, ymin=Min_GCP, ymax=Max_GCP, colour=Spirometry_measure)) +
  geom_pointrange() +
  geom_hline(yintercept=0, lty=2) +
  facet_wrap(~Spirometry_measure,strip.position="left",nrow=3,scales = "free_y") +
  coord_flip() +
  xlab("Metabolite or Hormone trait") + ylab(expression("GCP > 0 "%->% "Partial genetic causality on lung function" )) +
  theme_bw() +
  labs(colour="Spirometry measure")
  
 
