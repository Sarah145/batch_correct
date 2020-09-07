#! /usr/bin/env Rscript

sample_name <- commandArgs(trailingOnly = TRUE)

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(R2HTML)

methods <- c('raw', 'Seurat3','MNNCorrect', 'limma', 'Harmony', 'Combat', 'Scanorama', 'BBKNN')

pngs <- c()
for(i in 1:length(methods)){
  f <- paste0('data/', sample_name, '/asw/umap_', methods[i], '.', sample_name, '.png')
  if(file.exists(f)){
    system(paste('base64', f, '>', str_replace(f, '.png', '.txt')))
  }
  p <- ifelse(file.exists(f), readChar(str_replace(f, '.png', '.txt'), file.info(str_replace(f, '.png', '.txt'))$size), 'NA')
  pngs <- c(pngs, p)
}

batch_asws <- c()
for(i in 1:length(methods)){
  f <- ifelse(file.exists(paste0('data/', sample_name, '/asw/asw_', methods[i], '.', sample_name, '.csv')), read.csv(paste0('data/', sample_name, '/asw/asw_', methods[i], '.', sample_name, '.csv'), sep = '\t')$Batch_ASW, 'NA')
  batch_asws <- c(batch_asws, f)
  batch_asws <- unlist(batch_asws)
}

cell_asws <- c()
for(i in 1:length(methods)){
  f <- ifelse(file.exists(paste0('data/', sample_name, '/asw/asw_', methods[i], '.', sample_name, '.csv')), read.csv(paste0('data/', sample_name, '/asw/asw_', methods[i], '.', sample_name, '.csv'), sep = '\t')$Cell_type_ASW, 'NA')
  cell_asws <- c(cell_asws, f)
  cell_asws <- unlist(cell_asws)
}


batch_entropys <- c()
for(i in 1:length(methods)){
  f <- ifelse(file.exists(paste0('data/', sample_name, '/entropy/entropy_', methods[i], '.', sample_name, '.csv')), mean(read.csv(paste0('data/', sample_name, '/entropy/entropy_', methods[i], '.', sample_name, '.csv'), sep = '\t')[,1]), 'NA')
  batch_entropys <- c(batch_entropys, f)
  batch_entropys <- unlist(batch_entropys)
}

cell_entropys <- c()
for(i in 1:length(methods)){
  f <- ifelse(file.exists(paste0('data/', sample_name, '/entropy/entropy_', methods[i], '.', sample_name, '.csv')), mean(read.csv(paste0('data/', sample_name, '/entropy/entropy_', methods[i], '.', sample_name, '.csv'), sep = '\t')[,2]), 'NA')
  cell_entropys <- c(cell_entropys, f)
  cell_entropys <- unlist(cell_entropys)
}

jacc_sims <- c()
for(i in 1:length(methods)){
  f <- ifelse(file.exists(paste0('data/', sample_name, '/jaccard_sim/jacc_', methods[i], '.', sample_name, '.csv')), read.csv(paste0('data/', sample_name, '/jaccard_sim/jacc_', methods[i], '.', sample_name, '.csv'), sep = '\t')$Jacc_index, 'NA')
  jacc_sims <- c(jacc_sims, f)
  jacc_sims <- unlist(jacc_sims)
}

nf_report <- read.csv('trace.txt', sep = '\t')
nf_report$duration <- as.character(nf_report$duration)
process_methods <- c('raw', 'Seurat_3','MNNCorrect', 'limma', 'Harmony', 'Combat', 'Scanorama', 'BBKNN')
times <- c()
for(i in 1:length(methods)){
  r <- ifelse(length(unique(str_detect(nf_report$name, paste0('^', process_methods[i])))) > 1,  subset(nf_report, str_detect(nf_report$name, paste0('^', process_methods[i], '.*', sample_name, '\\)')) == T)$duration, 'NA')
  times <- c(times, r)
  }

df <- data.frame(method = methods, batch_asw = batch_asws, cell_asw = cell_asws, batch_entropy = batch_entropys, cell_entropy = cell_entropys, jacc_sim = jacc_sims)
for(i in 2:ncol(df)){df[,i] <- as.character(df[,i])}
df <- df %>% pivot_longer(c(batch_asw, cell_asw, batch_entropy, cell_entropy, jacc_sim), names_to = 'metric')
df$value <- as.character(df$value)
df <- df %>% mutate(value1 = case_when(value != 'NA' ~ value))
df$value1 <- as.numeric(df$value1)
df <- df %>% mutate(value1 = case_when(metric != 'batch_asw' & metric != 'cell_entropy' ~ value1,
                                 metric == 'batch_asw' ~ (1 - as.numeric(value)),
                                 metric == 'cell_entropy' ~ (1 - as.numeric(value))))
df$value1 <- as.numeric(df$value1)
p1 <- ggplot(df, aes(x = metric, y = method, fill = value1)) +
  geom_tile(col = 'white', size = 1) +
  scale_x_discrete(labels = c('1 - batch \nASW', 'Batch entropy',
                                                'Cell type ASW', '1 - cell type \nentropy',
                                                'Similarity of \nmarker genes')) +
  scale_y_discrete(limits = rev(unique(df$method))) +
  labs(x = NULL, y = NULL, fill = NULL) +  
  theme_void(base_size = 18) +
  theme(axis.text.x = element_text(angle = 70, hjust = 0.5, vjust = 0.5),
        axis.text = element_text(colour = 'black'))
ggsave(paste0(sample_name, '_heatmap.png'), plot = p1, height = 10, width = 8)

system(paste0('base64 ', sample_name, '_heatmap.png > ', sample_name, '_heatmap.txt'))
heatmap <- readChar(paste0(sample_name, '_heatmap.txt'), file.info(paste0(sample_name, '_heatmap.txt'))$size)

target <- HTMLInitFile(outdir = '.', filename=paste0(sample_name, '_summary'))
HTML('<link rel="stylesheet" href="https://fonts.googleapis.com/css?family=Roboto">', file=target)
HTML("<div class='title'>", file=target)
HTML.title(' Batch Correction Summary', HR=1, file = target)
HTML("</div>", file = target)
HTML.title(sample_name, HR=2, file = target)
for(i in 1:length(methods)){
  HTML.title(methods[i], HR=3, file = target)
  HTML(paste0("<img src='data:image/png;base64,", pngs[i], "' width=100%>"), file=target)
  HTML("<div class='boxed' align='center'>", file=target)
  HTML('<table style="width:100%">', file=target)
  HTML(paste('<tr> <td>', 'Runtime', '</td> <td align="right">', times[i], '</td> </tr>'), file=target)
  HTML(paste('<tr> <td>', 'Batch ASW', '</td> <td align="right">', batch_asws[i], '</td> </tr>'), file=target)
  HTML(paste('<tr> <td>', 'Cell type ASW', '</td> <td align="right">', cell_asws[i], '</td> </tr>'), file=target)
  HTML(paste('<tr> <td>', 'Batch entropy', '</td> <td align="right">', batch_entropys[i], '</td> </tr>'), file=target)
  HTML(paste('<tr> <td>', 'Cell type entropy', '</td> <td align="right">', cell_entropys[i], '</td> </tr>'), file=target)
  HTML(paste('<tr> <td>', 'Marker gene recovery (jaccard index)', '</td> <td align="right">', jacc_sims[i], '</td> </tr>'), file=target)
  HTML('</table>', file=target)
  HTML('</div>', file=target)
}
HTML(paste0("<img src='data:image/png;base64,", heatmap, "' width=100%>"), file=target)
HTML('<style type="text/css">
		.title {
    			background-color: #0972D5;
    			padding: 8px;
    			color: white;
    			position: fixed;
    			top: 0;
    			left: 0;
    			z-index: 999;
    			width: 100%;
		}
		.boxed {
  			border: 1px solid #868D96;
  			padding: 10px;
  			margin: 20px;
		}
		h1 {
			font-family: "Roboto";
			font-size: 33px;
		}
		h2 {
			font-family: "Roboto";
			font-size: 26px;
		}
		h3 {
			font-family: "Roboto";
			font-size: 18px;
		}
		table tr:nth-child(even) {
  			background-color: #eee;
		}
		table tr:nth-child(odd) {
  			background-color: #fff;
		}
		table {
  			font-family: "Roboto";
			font-size: 20px;
			border: 1px solid #868D96;
		}
		#mathplayer{
  			height: 80px;
		}
		</style> </head>', file=target)
HTMLEndFile()
  