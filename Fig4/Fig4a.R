# Draw Volcano plot
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(readxl)
library(EnhancedVolcano)
library(dplyr)

################################################################################

data <- read_excel('Fig4_data.xlsx', sheet = 'BIN1v1_Network')
keyvals <- as.character(data$color)
names(keyvals)[keyvals == '#4DAF4A'] <- 'APC/C'
names(keyvals)[keyvals == '#E78AC3'] <- 'CLT'
names(keyvals)[keyvals == '#FDB462'] <- 'RIN'
names(keyvals)[keyvals == '#FFD92F'] <- 'PIK3'
names(keyvals)[keyvals == '#FCCDE5'] <- 'TIM'
names(keyvals)[keyvals == '#BEBADA'] <- 'AP2'
names(keyvals)[keyvals == '#7DAEE0'] <- 'UP'
names(keyvals)[keyvals == '#ff7f00'] <- 'BIN1'
names(keyvals)[keyvals == 'grey'] <- 'NS'

filtered_data <- data %>% filter(color %in% c("#4DAF4A",'#FDB462'))
label_list = c(filtered_data$Genes,'BIN1')

EnhancedVolcano(data,
                lab = data$Genes,
                x = 'log2FC',
                y = 'adjpval',
                selectLab = data$Genes[which(names(keyvals) %in% c('BIN1','APC/C',
                                                                   'RIN','TIM','AP2'))],
                title = 'BIN1v1',
                subtitle = NULL,
                pCutoff = 0.1,
                FCcutoff = 1.29,
                pointSize = c(ifelse(data$Genes %in% label_list,4,3)),
                labSize = 10,
                colCustom = keyvals,
                colAlpha = c(ifelse(data$Genes %in% label_list,1,0.7)),
                drawConnectors = TRUE,
                legendPosition = '',
                widthConnectors = 1.0,
                arrowheads = FALSE,
                caption = '',
                axisLabSize = 30,
                ylab = bquote(~-log[10]~ '(Adjusted P)'),
                xlab = bquote(~log[2]~ '(Fold Change)')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  scale_x_continuous(limits = c(-1,6.5), breaks = seq(-1,6,1)) +
  scale_y_continuous(limits = c(0,5), breaks = seq(0,5,1))

#ggsave('BIN1_v1.svg',height = 10, width = 12)

