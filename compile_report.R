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
  p <- ifelse(file.exists(paste0('data/', sample_name, '/asw/umap_', methods[i], '.', sample_name, '.png')), paste0('data/', sample_name, '/asw/umap_', methods[i], '.', sample_name, '.png'), 'NA')
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
  f <- ifelse(file.exists(paste0('data/', sample_name, '/entropy/entropy_', methods[i], '.', sample_name, '.csv')), mean(read.csv(paste0('data/', sample_name, '/entropy/entropy_', methods[i], '.', sample_name, '.csv'), sep = '\t')$Batch_entropy), 'NA')
  batch_entropys <- c(batch_entropys, f)
  batch_entropys <- unlist(batch_entropys)
}

cell_entropys <- c()
for(i in 1:length(methods)){
  f <- ifelse(file.exists(paste0('data/', sample_name, '/entropy/entropy_', methods[i], '.', sample_name, '.csv')), mean(read.csv(paste0('data/', sample_name, '/entropy/entropy_', methods[i], '.', sample_name, '.csv'), sep = '\t')$Cell_type_entropy), 'NA')
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
  r <- ifelse(length(unique(str_detect(nf_report$name, paste0('^', process_methods[i])))) > 1,  subset(nf_report, str_detect(nf_report$name, paste0('^', process_methods[i])) == T)$duration, 'NA')
  times <- c(times, r)
  }

df <- data.frame(method = methods, batch_asw = batch_asws, cell_asw = cell_asws, batch_entropy = batch_entropys, cell_entropy = cell_entropys, jacc_sim = jacc_sims)
df <- df %>% pivot_longer(c(batch_asw, cell_asw, batch_entropy, cell_entropy, jacc_sim), names_to = 'metric')
df$value <- as.character(df$value)
df <- df %>% mutate(value1 = case_when(value != 'NA' ~ value))
df$value1 <- as.numeric(df$value1)
df <- df %>% mutate(value1 = case_when(TRUE ~ value1,
                                 metric == 'cell_asw' ~ (1 - value1),
                                 metric == 'cell_entropy' ~ (1 - value1)))
df$method <- factor(df$method, levels = unique(df$method))
ggplot(df, aes(x = metric, y = method, fill = value1)) +
  geom_tile(col = 'white', size = 1) +
  scale_x_discrete(position = "top", labels = c('Batch ASW', 'Batch entropy',
                                                '1 - cell type \nASW', '1 - cell type \nentropy',
                                                'Similarity of \nmarker genes')) +
  scale_y_discrete(limits = rev(unique(df$method))) +
  labs(x = NULL, y = NULL, fill = NULL) +
  scale_fill_distiller(palette = 'PuRd', direction = 1, limits = c(0,1)) +
  theme_minimal(base_size = 22) +
  theme(axis.text.x = element_text(angle = 70, hjust = 0.5, vjust = 0.5))
ggsave('heatmap.png', height = 16, width  = 10)


target <- HTMLInitFile(outdir = '.', filename=paste0(sample_name, '_summary'))
HTML('<link rel="stylesheet" href="https://fonts.googleapis.com/css?family=Roboto">', file=target)
HTML("<div class='title'>", file=target)
HTML.title(' Batch Correction Summary', HR=1, file = target)
HTML("</div>", file = target)
HTML.title(sample_name, HR=2, file = target)
for(i in 1:length(methods)){
  HTML.title(methods[i], HR=3, file = target)
  HTML(paste0("<img src='", pngs[i], "' width=100%>"), file=target)
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
HTML("<img src='heatmap.png' width=100%>", file=target)
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
  