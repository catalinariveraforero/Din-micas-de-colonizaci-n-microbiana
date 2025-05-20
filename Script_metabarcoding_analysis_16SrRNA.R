# Author: Catalina Rivera-Forero
# email: c.riveraf2@uniandes.edu.co | catalina.riveraf@utadeo.edu.co | catalinariveraforero.bm@gmail.com
# El siguiente script es una modificación del script desarrollado previamente en un estudio del microbioma de Madracis auretenra realizado por Ruiz-Toquica at al. (2025).
#Ruiz-Toquica J, Franco Herrera A, Medina M (2025) Endozoicomonas dominance and Vibrionaceae stability underpin resilience in urban coral Madracis auretenra. PeerJ 13:e19226. https://doi.org/10.7717/peerj.19226
## Licencia
# Este proyecto está licenciado bajo la Licencia BSD-3 Clause. Consulte el archivo [LICENSE] para más detalles.

###################################################### INSPECCIÓN, FILTRADO Y ASIGNACIÓN TAXONÓMICA DE DATOS ################################################
library(dada2)
setwd("/Users/ktari/OneDrive/Escritorio/tesis16S/")
path <- "/Users/ktari/OneDrive/Escritorio/tesis16S/" 
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)
plotQualityProfile(fnFs[1:4])
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen = c(250, 250), 
                     trimLeft = 10, 
                     maxN = 0, 
                     maxEE = c(2, 2), 
                     truncQ = 2, 
                     rm.phix = TRUE, 
                     compress = TRUE, 
                     verbose = TRUE,
                     multithread = FALSE)
out
(out[,1]-out[,2])/out[,1]*100
errF <- learnErrors(filtFs, multithread=FALSE, verbose = TRUE)
errR <- learnErrors(filtRs, multithread=FALSE, verbose = TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)
dadaFs[[57]]
dadaRs[[1]]
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/ktari/OneDrive/Escritorio/tesis16S/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=FALSE, verbose=TRUE)
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track
write.table(taxa, "taxa.txt", sep = "\t", quote = F)
write.table(seqtab.nochim, "asv.txt", sep = "\t", quote = F)
write.table(track, "track_datos_filtrados.txt", sep = "\t", quote = F)

#################################################### CREACIÓN OBJETO PHYLOSEQ, RAREFACCION Y FILTRADO ###############################################################################################
library(phyloseq)
library(tidyverse)
library(ggplot2)

asv  <- read_delim("/Users/ktari/OneDrive/Escritorio/tesis16S/asv_table.txt")
taxa   <- read_delim("/Users/ktari/OneDrive/Escritorio/tesis16S/taxa_table.txt")
metadata <- read_delim("/Users/ktari/OneDrive/Escritorio/tesis16S/metadata_16S_blocks.txt")

sampledata <- metadata |>
  column_to_rownames(var = "sampleID") |> 
  sample_data()

Year <- as.factor(sampledata$Year)
Location <- as.factor(sampledata$Location)
Type <- as.factor(sampledata$Type)
Month <- as.factor(sampledata$Month)
Depth <- as.factor(sampledata$Depth)

asvtable  <- asv |>
  column_to_rownames(var = "asvID") |> 
  as.matrix() |> 
  otu_table(taxa_are_rows = TRUE)

taxatable  <- taxa |>
  column_to_rownames(var = "asvID") |> 
  as.matrix() |> 
  tax_table()

phso <- phyloseq(asvtable, taxatable, sampledata)
phso

library(dplyr)
library(tidyr)
library(data.table)
library(writexl)

# NORMALIZACIÓN DE LAS ABUNDANCIAS DE ASVs POR LA MEDIANA DEL NÚMERO TOTAL DE MUESTRAS
total = median(sample_sums(phso))
standf = function(x, t=total) round(t * (x / sum(x)))
phso = transform_sample_counts(phso, standf)

readcount = data.table(as(sample_data(phso), "data.frame"),
                       TotalReads = sample_sums(phso), 
                       keep.rownames = TRUE)
setnames(readcount, "rn", "SampleID")
ggplot(readcount, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")

tab <- otu_table(phso)
class(tab) <- "matrix"
tab <- t(tab)
sum_seq <- rowSums(tab)

source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/load-extra-functions.R")
rarefaccion <- ggrare(phso, step = 100, color = "Location", label = "Sample", plot = TRUE, parallel = TRUE, se = FALSE) 

rarefac <- rarefaccion + theme_classic() +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        axis.ticks.length = unit(.25, "cm")) +
  theme(plot.margin = margin(2, 2, 2, 2, "cm")) +
  scale_x_continuous(labels = scales::comma) +
  scale_x_continuous(limits = c(0, 18800), breaks = seq(0, 20000, by = 2000)) +
  scale_color_manual(values = c("cadetblue3", "darkseagreen")) +
  theme(plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"))

rarefac

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Curva_rarefaccion.png", 
       rarefac, width = 9, height = 8, dpi = 600)

#ELIMINACIÓN DE ASVs SIN ASIGNACIÓN Y SECUENCIAS DE CLOROPLASTO Y MITROCONDRIA
table(tax_table(phso) [, "Phylum"], exclude = NULL)
phso <- subset_taxa(phso, !is.na(Phylum) & !Phylum %in% c("", "<NA>"))
table(tax_table(phso) [, "Class"], exclude = NULL)
phso <- subset_taxa(phso, !is.na(Order) & !Order %in% c("", "Chloroplast"))
table(tax_table(phso) [, "Family"], exclude = NULL)
phso <- subset_taxa(phso, !is.na(Family) & !Family %in% c("", "Mitochondria"))
phso
#SE VUELVE A NORMALIZAR DESPUÉS DEL FILTRADO
total = median(sample_sums(phso))
standf = function(x, t=total) round(t * (x / sum(x)))
phso = transform_sample_counts(phso, standf)

#PREVALENCIA
prevdf <- apply(X = otu_table(phso),
                MARGIN = ifelse(taxa_are_rows(phso), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(phso),
                     tax_table(phso))
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence), sum(df1$Prevalence))})

prevdf1 <- subset(prevdf, Phylum %in% get_taxa_unique(phso, "Phylum"))

#Dividir los filos en tres grupos de 13
unique_phyla <- unique(prevdf1$Phylum)
phylum_groups <- split(unique_phyla, ceiling(seq_along(unique_phyla) / 13))
plots <- list()

for (i in seq_along(phylum_groups)) {
  group_phylum <- phylum_groups[[i]]
  subset_prevdf <- subset(prevdf1, Phylum %in% group_phylum)
  plots[[i]] <- ggplot(subset_prevdf, aes(TotalAbundance, Prevalence / nsamples(phso), color=Phylum)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +
    scale_x_log10() +
    xlab("Total Abundance") + 
    ylab("Prevalence [Frac. Samples]") +
    facet_wrap(~Phylum) + 
    theme(legend.position="none") +
    theme_classic() +
    ggtitle(paste("Group", i, "of Phyla"))

  ggsave(paste0("/Users/ktari/OneDrive/Escritorio/tesis16S/Prevalencia_phylum", i, ".png"), 
         plots[[i]], width = 8, height = 6, dpi = 600)
}

#ELIMACIÓN DE ASVs CON PREVALENCIA < 1 % (POSIBLEMENTE ARTEFACTOS)
prevalenceThreshold <- 0.01 * nsamples(phso)
keepTaxa <- rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
phso1 <- prune_taxa(keepTaxa, phso)
keepTaxa <- rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]

#SE VUELVE A NORMALIZAR DESPUÉS DEL FILTRADO DE ASVs CON PREVALENCIA < 1%
total = median(sample_sums(phso1))
standf = function(x, t=total) round(t * (x / sum(x)))
phso1 = transform_sample_counts(phso1, standf)

write.csv(cbind(data.frame(otu_table(phso1)),
                tax_table(phso1)), 
          file="/Users/ktari/OneDrive/Escritorio/tesis16S/ASvs_filtradas.csv")
phso

############################################### ABUNDANCIA RELATIVA ###############################################################
library(dplyr)
library(RColorBrewer)
library(paletteer)

phso1_phylum <- tax_glom(phso, "Phylum", NArm = TRUE)
phso1_phylum_relabun <- transform_sample_counts(phso1_phylum, function(OTU) OTU/sum(OTU) * 100)
taxa_abundance_table_phylum <- psmelt(phso1_phylum_relabun)
phso1.rel = transform_sample_counts(phso1, function(x) x/sum(x)*100)

#BARPLOT A NIVEL DE PHYLUM
glom <- tax_glom(phso1.rel, taxrank = 'Phylum', NArm = FALSE)
phso1.melt <- psmelt(glom)
phso1.melt$Phylum <- as.character(phso1.melt$Phylum)
phso1.melt <- phso1.melt %>%
  group_by(Sample, Month, Location, Phylum) %>%
  mutate(median=median(Abundance))
keep <- unique(phso1.melt$Phylum[phso1.melt$median > 0])
phso1.melt$Phylum[!(phso1.melt$Phylum %in% keep)] <- "< 0%"
phso1.melt_sum <- phso1.melt %>%
  group_by(Sample, Month, Location, Phylum) %>%
  summarise(Abundance=sum(Abundance))

#BARPLOT PHYLUM
phy <- ggplot(phso1.melt_sum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", aes(fill=Phylum), alpha = 0.9) + 
  labs(x="", y="Relative abundance (%)") +
  facet_wrap(~Location + Month, scales= "free_x", nrow=1) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x.bottom = element_text(angle = 90, size = 7),
        axis.title=element_text(size=18),
        axis.text=element_text(size=14),
        legend.text=element_text(size=14),
        legend.title = element_text(size = 18),
        strip.text = element_text(size = 16),
        axis.ticks.length = unit(.3, "cm")) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(100))

phy

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Barplot_phylum.png", 
       phy, width = 12, height = 6, dpi = 600)


#BARPLOT A NIVEL DE FAMILIA
glom <- tax_glom(phso1.rel, taxrank = 'Family', NArm = FALSE)
phso1.melt <- psmelt(glom)
phso1.melt$Family <- as.character(phso1.melt$Family)
phso1.melt <- phso1.melt %>%
  group_by(Sample, Month, Location, Family) %>%
  mutate(median=median(Abundance))
keep <- unique(phso1.melt$Family[phso1.melt$median > 3])
phso1.melt$Family[!(phso1.melt$Family %in% keep)] <- "< 3%"
phso1.melt_sum <- phso1.melt %>%
  group_by(Sample, Month, Location, Family) %>%
  summarise(Abundance=sum(Abundance))


fam <- ggplot(phso1.melt_sum, aes(x = Sample, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity", aes(fill=Family), alpha = 0.9) + 
  labs(x="", y="Relative abundance (%)") +
  facet_wrap(~Location + Month, scales= "free_x", nrow=1) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x.bottom = element_text(angle = 90, size = 7),  
        axis.title=element_text(size=18), 
        axis.text=element_text(size=14),    
        legend.text=element_text(size=14), 
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 16),   
        axis.ticks.length = unit(.3, "cm")) +   
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(100))

fam

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Barplot_familia.png", 
       fam, width = 13.5, height = 5.7, dpi = 600)


#BARPLOT A NIVEL DE GENERO
glom <- tax_glom(phso1.rel, taxrank = 'Genus', NArm = FALSE)
phso1.melt <- psmelt(glom)
phso1.melt$Genus <- as.character(phso1.melt$Genus)
phso1.melt <- phso1.melt %>%
  group_by(Sample, Month, Location, Genus) %>%
  mutate(median=median(Abundance))
keep <- unique(phso1.melt$Genus[phso1.melt$median > 2.5])
phso1.melt$Genus[!(phso1.melt$Genus %in% keep)] <- "< 2.5%"
phso1.melt_sum <- phso1.melt %>%
  group_by(Sample, Month, Location, Genus) %>%
  summarise(Abundance=sum(Abundance))

gen <- ggplot(phso1.melt_sum, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", aes(fill=Genus), alpha = 0.9) + 
  labs(x="", y="Relative abundance (%)") +
  facet_wrap(~Location + Month, scales= "free_x", nrow=1) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x.bottom = element_text(angle = 90, size = 7)) +  
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=14), 
        legend.text=element_text(size=14),
        legend.title = element_text(size = 18),
        strip.text = element_text(size = 16),
        axis.ticks.length = unit(.3, "cm")) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(100))

gen

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Barplot_genero.png", 
       gen, width = 13.5, height = 5.7, dpi = 600)


# TABLAS DE LAS ABUNDANCIAS RELATIVAS DE CADA TAXÓN
library(dplyr)
library(writexl)

taxrank <- c("Phylum", "Family", "Genus")
# Función para calcular el porcentaje de abundancia relativa por ubicación y nivel taxonómico
asv_relative_percentage <- function(phso, taxrank, location, output_path) {
  # Filtrar las muestras para la ubicación específica
  phso_location <- subset_samples(phso, Location == location)
  # Normalizar las abundancias relativas
  phso_rel <- transform_sample_counts(phso_location, function(x) x / sum(x) * 100)
  # Extraer la tabla taxonómica
  taxa_table <- as.data.frame(tax_table(phso_rel))
  # Extraer la tabla de abundancia
  abundance_table <- as.data.frame(otu_table(phso_rel))
  combined_table <- cbind(taxa_table, abundance_table)
  percentage_table <- combined_table %>%
    group_by(.data[[taxrank]]) %>%
    summarise(Relative_Abundance = sum(across(where(is.numeric)))) %>%
    mutate(Percentage = Relative_Abundance / sum(Relative_Abundance) * 100)
  file_path <- paste0(output_path, "/ASV_relative_percentage_", taxrank, "_", location, ".xlsx")
  write_xlsx(percentage_table, file_path)
  print(paste("Resultados para", location, "-", taxrank))
  print(percentage_table)
}

output_path <- "/Users/ktari/OneDrive/Escritorio/tesis16S"
locations <- c("GOR", "SAI")
taxranks <- c("Phylum", "Family", "Genus")

for (location in locations) {
  for (taxrank in taxranks) {
    asv_relative_percentage(phso1, taxrank, location, output_path)
  }
}

# ANÁLISIS ESTADÍSTICO COMPARACIÓN DE CADA TAXÓN ENTRE CONDICIONES
library(openxlsx)
library(dplyr)
library(FSA)

melted_data <- psmelt(phso1)
melted_data$Depth <- as.factor(melted_data$Depth)  
melted_data$Month <- as.factor(melted_data$Month)  
melted_data$Location <- as.factor(melted_data$Location)

locations <- unique(melted_data$Location)
taxranks <- c("Phylum", "Family", "Genus")
months <- unique(melted_data$Month)
depths <- unique(melted_data$Depth) 

output_path <- "/Users/ktari/OneDrive/Escritorio/tesis16S/stats_abundancia_relativa"
if (!dir.exists(output_path)) {
  dir.create(output_path)
}
results <- list()
# Kruskal-Wallis y Dunn Test
add_results <- function(data, group_var, taxa, taxa_name, comparison_results) {
  if (length(unique(data[[group_var]])) > 1) {
    # Kruskal-Wallis Test
    kruskal_result <- tryCatch(
      kruskal.test(Abundance ~ as.factor(data[[group_var]]), data = data),
      error = function(e) NULL
    )
    
    if (!is.null(kruskal_result)) {
      # Agregar resultados de Kruskal-Wallis
      comparison_results[[taxa_name]] <- list(Kruskal_Wallis_p = kruskal_result$p.value)
      
      # Test de Dunn
      dunn_result <- tryCatch(
        dunnTest(Abundance ~ as.factor(data[[group_var]]), data = data, method = "bh"),
        error = function(e) NULL
      )
      
      if (!is.null(dunn_result)) {
        comparison_results[[taxa_name]]$Dunn_Test <- dunn_result$res
      } else {
        comparison_results[[taxa_name]]$Dunn_Test <- "Dunn Test no pudo ejecutarse"
      }
    } else {
      comparison_results[[taxa_name]] <- "Kruskal-Wallis no pudo ejecutarse"
    }
  } else {
    comparison_results[[taxa_name]] <- "Menos de dos grupos con datos distintos"
  }
  return(comparison_results)
}

# 1. Comparaciones entre meses en cada ubicación
for (location in locations) {
  for (taxrank in taxranks) {
    comparison_name <- paste("Comparaciones_meses", location, taxrank, sep = "_")
    comparison_results <- list()
    
    data_filtered <- melted_data %>%
      filter(Location == location) %>%
      group_by(Sample, Month, !!sym(taxrank)) %>%
      summarise(Abundance = sum(Abundance), .groups = "drop")
    
    for (taxa in unique(data_filtered[[taxrank]])) {
      taxa_data <- data_filtered %>% filter(!!sym(taxrank) == taxa)
      comparison_results <- add_results(taxa_data, "Month", taxa, taxa, comparison_results)
    }
    
    results[[comparison_name]] <- comparison_results
  }
}

# 2. Comparaciones entre profundidades en cada ubicación
for (location in locations) {
  for (taxrank in taxranks) {
    comparison_name <- paste("Comparaciones_profundidades", location, taxrank, sep = "_")
    comparison_results <- list()
    
    data_filtered <- melted_data %>%
      filter(Location == location) %>%
      group_by(Sample, Depth, !!sym(taxrank)) %>%
      summarise(Abundance = sum(Abundance), .groups = "drop")
    
    for (taxa in unique(data_filtered[[taxrank]])) {
      taxa_data <- data_filtered %>% filter(!!sym(taxrank) == taxa)
      comparison_results <- add_results(taxa_data, "Depth", taxa, taxa, comparison_results)
    }
    
    results[[comparison_name]] <- comparison_results
  }
}

# 3. Comparaciones entre ubicaciones
for (taxrank in taxranks) {
  comparison_name <- paste("Comparaciones_ubicaciones", taxrank, sep = "_")
  comparison_results <- list()
  
  data_filtered <- melted_data %>%
    group_by(Sample, Location, !!sym(taxrank)) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")
  
  for (taxa in unique(data_filtered[[taxrank]])) {
    taxa_data <- data_filtered %>% filter(!!sym(taxrank) == taxa)
    comparison_results <- add_results(taxa_data, "Location", taxa, taxa, comparison_results)
  }
  
  results[[comparison_name]] <- comparison_results
}

# Guardar resultados
for (comparison_name in names(results)) {
  comparison_data <- results[[comparison_name]]
  wb <- createWorkbook()
  addWorksheet(wb, "Resultados")
  row <- 1
  for (taxa in names(comparison_data)) {
    writeData(wb, "Resultados", paste("Taxa:", taxa), startRow = row, startCol = 1)
    row <- row + 1  # Avanzar a la siguiente fila
    if (is.list(comparison_data[[taxa]])) {
      if (!is.null(comparison_data[[taxa]]$Kruskal_Wallis_p)) {
        writeData(wb, "Resultados", paste("Kruskal-Wallis p-value:", comparison_data[[taxa]]$Kruskal_Wallis_p), startRow = row, startCol = 1)
        row <- row + 1
      }
      if (!is.null(comparison_data[[taxa]]$Dunn_Test)) {
        if (is.data.frame(comparison_data[[taxa]]$Dunn_Test)) {
          writeData(wb, "Resultados", comparison_data[[taxa]]$Dunn_Test, startRow = row, startCol = 1)
          row <- row + nrow(comparison_data[[taxa]]$Dunn_Test) + 2
        } else {
          writeData(wb, "Resultados", data.frame(Message = comparison_data[[taxa]]$Dunn_Test), startRow = row, startCol = 1)
          row <- row + 2
        }
      }
    } else {
      writeData(wb, "Resultados", data.frame(Message = comparison_data[[taxa]]), startRow = row, startCol = 1)
      row <- row + 2
    }
  }
  saveWorkbook(wb, file.path(output_path, paste0(comparison_name, ".xlsx")), overwrite = TRUE)
}


# Análisis que busca diferencias en las abundancias relativas de las comunidades microbiológicas entre diferentes condiciones
library(dplyr)
library(FSA)

# 1. Comparaciones entre meses en cada ubicación
cat("### Comparaciones generales entre meses por ubicación ###\n")
for (location in locations) {
  cat("\nUbicación:", location, "\n")
  data_filtered <- melted_data %>%
    filter(Location == location) %>%
    group_by(Sample, Month) %>%
    summarise(Total_Abundance = sum(Abundance), .groups = "drop")
  
  if (length(unique(data_filtered$Month)) > 1) {
    # Kruskal-Wallis Test
    kruskal_result <- tryCatch(
      kruskal.test(Total_Abundance ~ as.factor(Month), data = data_filtered),
      error = function(e) NULL
    )
    
    if (!is.null(kruskal_result)) {
      cat("Kruskal-Wallis p-value:", kruskal_result$p.value, "\n")
    } else {
      cat("Kruskal-Wallis no pudo ejecutarse\n")
    }
    
    # Test de Dunn (independientemente del resultado de Kruskal-Wallis)
    dunn_result <- tryCatch(
      dunnTest(Total_Abundance ~ as.factor(Month), data = data_filtered, method = "bh"),
      error = function(e) NULL
    )
    
    if (!is.null(dunn_result)) {
      cat("Dunn Test Resultados:\n")
      print(dunn_result$res)
    } else {
      cat("Dunn Test no pudo ejecutarse\n")
    }
  } else {
    cat("Menos de dos meses con datos disponibles.\n")
  }
}

# 2. Comparaciones entre profundidades en cada ubicación
cat("\n### Comparaciones generales entre profundidades por ubicación ###\n")
for (location in locations) {
  cat("\nUbicación:", location, "\n")
  data_filtered <- melted_data %>%
    filter(Location == location) %>%
    group_by(Sample, Depth) %>%
    summarise(Total_Abundance = sum(Abundance), .groups = "drop")
  
  if (length(unique(data_filtered$Depth)) > 1) {
    # Kruskal-Wallis Test
    kruskal_result <- tryCatch(
      kruskal.test(Total_Abundance ~ as.factor(Depth), data = data_filtered),
      error = function(e) NULL
    )
    
    if (!is.null(kruskal_result)) {
      cat("Kruskal-Wallis p-value:", kruskal_result$p.value, "\n")
    } else {
      cat("Kruskal-Wallis no pudo ejecutarse\n")
    }
    
    # Test de Dunn (independientemente del resultado de Kruskal-Wallis)
    dunn_result <- tryCatch(
      dunnTest(Total_Abundance ~ as.factor(Depth), data = data_filtered, method = "bh"),
      error = function(e) NULL
    )
    
    if (!is.null(dunn_result)) {
      cat("Dunn Test Resultados:\n")
      print(dunn_result$res)
    } else {
      cat("Dunn Test no pudo ejecutarse\n")
    }
  } else {
    cat("Menos de dos profundidades con datos disponibles.\n")
  }
}

# 3. Comparaciones entre ubicaciones
cat("\n### Comparaciones generales entre ubicaciones ###\n")
data_filtered <- melted_data %>%
  group_by(Sample, Location) %>%
  summarise(Total_Abundance = sum(Abundance), .groups = "drop")

if (length(unique(data_filtered$Location)) > 1) {
  # Kruskal-Wallis Test
  kruskal_result <- tryCatch(
    kruskal.test(Total_Abundance ~ as.factor(Location), data = data_filtered),
    error = function(e) NULL
  )
  
  if (!is.null(kruskal_result)) {
    cat("Kruskal-Wallis p-value:", kruskal_result$p.value, "\n")
  } else {
    cat("Kruskal-Wallis no pudo ejecutarse\n")
  }
  
  # Test de Dunn (independientemente del resultado de Kruskal-Wallis)
  dunn_result <- tryCatch(
    dunnTest(Total_Abundance ~ as.factor(Location), data = data_filtered, method = "bh"),
    error = function(e) NULL
  )
  
  if (!is.null(dunn_result)) {
    cat("Dunn Test Resultados:\n")
    print(dunn_result$res)
  } else {
    cat("Dunn Test no pudo ejecutarse\n")
  }
} else {
  cat("Menos de dos ubicaciones con datos disponibles.\n")
}


phso1_SAI <- subset_samples(phso1, Location == "SAI")

phso1_GOR <- subset_samples(phso1, Location == "GOR")

phso1_SAI
phso1_GOR 

########################################################################### ANALISIS PREVALENCIA Y ABUNDANCIA ###############################################################################
library(phyloseq)
library(ggplot2)
library(dplyr)
taxa <- c("Phylum", "Family", "Genus")

# CALCULO DE PREVALENCIA (número de muestras en las que cada taxón está presente) Y ABUINDANCIA TOTAL  (suma de lecturas de cada taxón en todas las muestras)
calculate_prevalence_abundance <- function(physeq) {
  prevalence <- apply(otu_table(physeq), 1, function(x) sum(x > 0))
  abundance_total <- rowSums(otu_table(physeq))
  prev_abund <- data.frame(
    Taxa = rownames(otu_table(physeq)),
    Prevalence = prevalence,
    TotalAbundance = abundance_total,
    tax_table(physeq)
  )
  return(prev_abund)
}

prev_abund_all <- calculate_prevalence_abundance(phso)
prev_abund_SAI <- calculate_prevalence_abundance(phso1_SAI)
prev_abund_GOR <- calculate_prevalence_abundance(phso1_GOR)

prev_abund_by_month <- list()
for (location in c("SAI", "GOR")) {
  for (month in unique(sample_data(phso1)$Month)) {
    subset_physeq <- subset_samples(phso1, Location == location & Month == month)
    prev_abund_by_month[[paste(location, month, sep = "_")]] <- calculate_prevalence_abundance(subset_physeq)
  }
}

# GRAFICO PREVALENCIA VS. ABUNDANCIA 
plot_prevalence_abundance <- function(prev_abund, title = "") {
  ggplot(prev_abund, aes(x = Prevalence, y = TotalAbundance, color = Phylum)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_y_log10() +
    labs(
      x = "Prevalence (Number of Samples)",
      y = "Total Abundance (Log Scale)",
      title = title
    ) +
    theme_classic() +
    theme(legend.position = "bottom")
}

library(RColorBrewer)

plot_prevalence_abundance <- function(prev_abund, title = "") {
  ggplot(prev_abund, aes(x = Prevalence, y = TotalAbundance, color = Phylum)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_y_log10() +  
    scale_color_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(length(unique(prev_abund$Phylum)))) +  
    labs(
      x = "Prevalence (Number of Samples)",
      y = "Total Abundance (Log Scale)",
      title = title
    ) +
    theme_classic() +
    theme(legend.position = "bottom")
}

library(viridis)

plot_prevalence_abundance <- function(prev_abund, title = "") {
  ggplot(prev_abund, aes(x = Prevalence, y = TotalAbundance, color = Phylum)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_y_log10() + 
    scale_color_viridis(discrete = TRUE, option = "D") + 
    labs(
      x = "Prevalence (Number of Samples)",
      y = "Total Abundance (Log Scale)",
      title = title
    ) +
    theme_classic() +
    theme(legend.position = "bottom")
}



plot_prevalence_abundance(prev_abund_all, title = "All Samples")
plot_prevalence_abundance(prev_abund_SAI, title = "SAI Samples")
plot_prevalence_abundance(prev_abund_GOR, title = "GOR Samples")


for (name in names(prev_abund_by_month)) {
  print(plot_prevalence_abundance(prev_abund_by_month[[name]], title = name))
}


# ID TAXONES CON ALTA PREVALENCIA Y ABUNDANCIA
high_prev_abund <- prev_abund_all %>%
  filter(Prevalence > 0.5 * nsamples(phso1), TotalAbundance > median(TotalAbundance))
head(high_prev_abund)

output_dir <- "/Users/ktari/OneDrive/Escritorio/tesis16S/Prevalencia_abundancia_graficos"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

save_plot <- function(plot, name, output_dir) {
  file_path <- file.path(output_dir, paste0(name, ".png"))
  ggsave(file_path, plot, width = 10, height = 8, dpi = 600)
}


plot_all <- plot_prevalence_abundance(prev_abund_all, title = "All Samples")
save_plot(plot_all, "all_samples", output_dir)

plot_SAI <- plot_prevalence_abundance(prev_abund_SAI, title = "SAI Samples")
save_plot(plot_SAI, "SAI_samples", output_dir)

plot_GOR <- plot_prevalence_abundance(prev_abund_GOR, title = "GOR Samples")
save_plot(plot_GOR, "GOR_samples", output_dir)

for (name in names(prev_abund_by_month)) {
  plot <- plot_prevalence_abundance(prev_abund_by_month[[name]], title = name)
  save_plot(plot, name, output_dir)
}


################################################# ANÁLISIS DE CORRELQCIÓN ENTRE TIEMPO Y ABUNBDANCIA #########################################################################
library(phyloseq)
library(dplyr)
library(ggplot2)
library(openxlsx)

# CORRELACIÓN ENRE TIEMPO Y ABUNDANCIA PARA UN TAXÓN EN UNA UBICACIÓN
calculate_correlation <- function(taxa_data) {
  if (nrow(taxa_data) > 1) {
    correlation_result <- cor.test(taxa_data$Month, taxa_data$Abundance, method = "spearman")
    return(correlation_result)
  } else {
    return(NULL)
  }
}

# TASA DE CAMBIO EN LA ABUND. DE UN TAXÓN EN UNA UBICACIÓN
calculate_taxa_change <- function(taxa_data) {
  if (nrow(taxa_data) > 1) {
    taxa_data <- taxa_data %>%
      group_by(Month) %>%
      summarise(MeanAbundance = mean(Abundance), .groups = "drop") %>%
      mutate(Change = c(NA, diff(MeanAbundance)))
    return(taxa_data)
  } else {
    return(NULL)
  }
}

melted_data <- psmelt(phso1)
tax_levels <- c("Phylum", "Family", "Genus")
locations <- unique(melted_data$Location)
correlation_results <- list()
change_results <- list()

# CALCULO CORRELACIÓN TASA DE CMBIO PARA C/D NICEL TAXONÓMICO EN C/D UBICACIÓN
for (tax_level in tax_levels) {
  unique_taxa <- unique(melted_data[[tax_level]])
  for (location in locations) {
    for (taxa in unique_taxa) {
      taxa_data <- melted_data %>%
        filter(!!sym(tax_level) == taxa, Location == location)
      correlation_results[[paste(location, tax_level, taxa, sep = "_")]] <- calculate_correlation(taxa_data)
      change_results[[paste(location, tax_level, taxa, sep = "_")]] <- calculate_taxa_change(taxa_data)
    }
  }
}

extract_correlation_results <- function(correlation_result) {
  if (is.null(correlation_result) || is.na(correlation_result$estimate)) {
    return(NULL)
  }
  data.frame(
    Estimate = correlation_result$estimate,
    P_value = correlation_result$p.value,
    Method = correlation_result$method,
    Alternative = correlation_result$alternative
  )
}


correlation_df <- data.frame()
for (name in names(correlation_results)) {
  result <- correlation_results[[name]]
  result_df <- extract_correlation_results(result)
  if (!is.null(result_df)) {
    result_df$Taxa_Location <- name  # Agregar el nombre del taxón y ubicación
    correlation_df <- rbind(correlation_df, result_df)
  }
}

output_dir <- "/Users/ktari/OneDrive/Escritorio/tesis16S/analisis_sucesion"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
write.xlsx(correlation_df, file.path(output_dir, "resultados_correlacion.xlsx"), rowNames = FALSE)
shorten_name <- function(name) {
  if (nchar(name) > 31) {
    return(substr(name, 1, 31))  # Truncar el nombre a 31 caracteres
  } else {
    return(name)
  }
}
wb <- createWorkbook()

for (name in names(change_results)) {
  result <- change_results[[name]]
  if (!is.null(result)) {
    sheet_name <- shorten_name(name)
    addWorksheet(wb, sheetName = sheet_name)
    writeData(wb, sheet = sheet_name, x = result)
  }
}

saveWorkbook(wb, file.path(output_dir, "resultados_tasa_cambio.xlsx"), overwrite = TRUE)

#################################################### CREACCIÓN DE SUBSETS ##################################################################################################################################################################################################################################################################################################################################################################################################################
phso1_SAI_2 <- subset_samples(phso1_SAI, Month == "2")

phso1_SAI_5 <- subset_samples(phso1_SAI, Month == "6")

phso1_SAI_7 <- subset_samples(phso1_SAI, Month == "8")

phso1_SAI_9 <- subset_samples(phso1_SAI, Month == "10")

phso1_GOR_2 <- subset_samples(phso1_GOR, Month == "2")

phso1_GOR_5 <- subset_samples(phso1_GOR, Month == "6")

phso1_GOR_7 <- subset_samples(phso1_GOR, Month == "8")

phso1_GOR_9 <- subset_samples(phso1_GOR, Month == "10")

################################################## ANALISIS DE CONJUNTOS ##########################################################################################################
library(phyloseq)
library(MicEco)
library(ggplot2)
library(ggpubr)

phso1_phylum <- tax_glom(phso1, taxrank = "Phylum")
phso1_family <- tax_glom(phso1, taxrank = "Family")
phso1_genus <- tax_glom(phso1, taxrank = "Genus")


venn_phylum <- ps_venn(phso1_phylum, group = "Location", fill = c("cadetblue3", "darkseagreen"), 
                     alpha = 0.5, quantities = list(type = c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))

venn_familia <- ps_venn(phso1_family, group = "Location", fill = c("cadetblue3", "darkseagreen"), 
                     alpha = 0.5, quantities = list(type = c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))

venn_genero <- ps_venn(phso1_genus, group = "Location", fill = c("cadetblue3", "darkseagreen"), 
                    alpha = 0.5, quantities = list(type = c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))

venn <- ggarrange(venn_phylum, venn_familia, venn_genero,
                         ncol = 3, nrow = 1)
venn

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Diagrama_venn_Ubicacion_phylum_familia_genero.png", 
       venn, width = 8, height = 6, dpi = 600)

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/plot_venn_Location_genero_2024.png", 
       venn_genero, width = 7, height = 6, dpi = 600)


phso1_SAI_phylum <- tax_glom(phso1_SAI, taxrank = "Phylum")
phso1_SAI_family <- tax_glom(phso1_SAI, taxrank = "Family")
phso1_SAI_genus <- tax_glom(phso1_SAI, taxrank = "Genus")

venn_phylum_SAI <- ps_venn(phso1_SAI_phylum, group = "Month", fill = c("pink2", "palevioletred", "plum3", "paleturquoise3"), 
                         alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.01, weight = TRUE, col = "grey16"))

venn_familia_SAI <- ps_venn(phso1_SAI_family, group = "Month", fill = c("pink2", "palevioletred", "plum3", "paleturquoise3"), 
                         alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.01, weight = TRUE, col = "grey16"))

venn_genero_SAI <- ps_venn(phso1_SAI_genus, group = "Month", fill = c("pink2", "palevioletred", "plum3", "paleturquoise3"), 
                        alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.01, weight = TRUE, col = "grey16"))

venn_SAI <- ggarrange(venn_phylum_SAI, venn_familia_SAI,venn_genero_SAI,
                         ncol = 3, nrow = 1)
venn_SAI

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Diagrama_venn_SAI_Meses_phylum_familia_genero.png", 
       venn_SAI, width = 20, height = 6, dpi = 600)

venn_SAI <- ggarrange(venn_familia_SAI,venn_genero_SAI,
                         ncol = 2, nrow = 1)
venn_SAI

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Diagrama_venn_SAI_Meses_familia_genero.png", 
       venn_SAI, width = 10, height = 4, dpi = 600)


phso1_GOR_phylum <- tax_glom(phso1_GOR, taxrank = "Phylum")
phso1_GOR_family <- tax_glom(phso1_GOR, taxrank = "Family")
phso1_GOR_genus <- tax_glom(phso1_GOR, taxrank = "Genus")

venn_phylum_GOR <- ps_venn(phso1_GOR_phylum, group = "Month", fill = c("pink2", "palevioletred", "plum3", "paleturquoise3"), 
                         alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.01, weight = TRUE, col = "grey16"))

venn_familia_GOR <- ps_venn(phso1_GOR_family, group = "Month", fill = c("pink2", "palevioletred", "plum3", "paleturquoise3"), 
                         alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.01, weight = TRUE, col = "grey16"))

venn_genero_GOR <- ps_venn(phso1_GOR_genus, group = "Month", fill = c("pink2", "palevioletred", "plum3", "paleturquoise3"), 
                        alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.01, weight = TRUE, col = "grey16"))

venn_GOR <- ggarrange(venn_phylum_GOR, venn_familia_GOR, venn_genero_GOR,
                         ncol = 3, nrow = 1)
venn_GOR
ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Diagrama_venn_GOR_Meses_phylum_familia_genero.png", 
       venn_GOR, width = 20, height = 6, dpi = 600)

venn_GOR <- ggarrange(venn_familia_GOR,venn_genero_GOR,
                         ncol = 2, nrow = 1)
venn_GOR
ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Diagrama_venn_GOR_Meses_familia_genero.png", 
       venn_GOR, width = 10, height = 4, dpi = 600)

phso1_SAI_phylum_2 <- tax_glom(phso1_SAI_2, taxrank = "Phylum")
phso1_SAI_family_2 <- tax_glom(phso1_SAI_2, taxrank = "Family")
phso1_SAI_genus_2 <- tax_glom(phso1_SAI_2, taxrank = "Genus")


venn_phylum_SAI_2_prof <- ps_venn(phso1_SAI_phylum_2 , group = "Depth", fill = c("6" = "beige", "12" = "bisque3"), 
                                alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))

venn_familia_SAI_2_prof <- ps_venn(phso1_SAI_family_2 , group = "Depth", fill = c("6" = "beige", "12" = "bisque3"), 
                                alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))

venn_genero_SAI_2_prof <- ps_venn(phso1_SAI_genus_2 , group = "Depth", fill = c("6" = "beige", "12" = "bisque3"), 
                               alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))

phso1_SAI_phylum_5 <- tax_glom(phso1_SAI_5, taxrank = "Phylum")
phso1_SAI_family_5 <- tax_glom(phso1_SAI_5, taxrank = "Family")
phso1_SAI_genus_5 <- tax_glom(phso1_SAI_5, taxrank = "Genus")

venn_phylum_SAI_5_prof <- ps_venn(phso1_SAI_phylum_5 , group = "Depth", fill = c("6" = "beige", "12" = "bisque3"), 
                                alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))

venn_familia_SAI_5_prof <- ps_venn(phso1_SAI_family_5 , group = "Depth", fill = c("6" = "beige", "12" = "bisque3"), 
                                alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))

venn_genero_SAI_5_prof <- ps_venn(phso1_SAI_genus_5 , group = "Depth", fill = c("6" = "beige", "12" = "bisque3"), 
                               alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))

phso1_SAI_phylum_7 <- tax_glom(phso1_SAI_7, taxrank = "Phylum")
phso1_SAI_family_7 <- tax_glom(phso1_SAI_7, taxrank = "Family")
phso1_SAI_genus_7 <- tax_glom(phso1_SAI_7, taxrank = "Genus")


venn_phylum_SAI_7_prof <- ps_venn(phso1_SAI_phylum_7 , group = "Depth", fill = c("6" = "beige", "12" = "bisque3"), 
                                alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))

venn_familia_SAI_7_prof <- ps_venn(phso1_SAI_family_7 , group = "Depth", fill = c("6" = "beige", "12" = "bisque3"), 
                                alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))

venn_genero_SAI_7_prof <- ps_venn(phso1_SAI_genus_7 , group = "Depth", fill = c("6" = "beige", "12" = "bisque3"), 
                               alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))


phso1_SAI_phylum_9 <- tax_glom(phso1_SAI_9, taxrank = "Phylum")
phso1_SAI_family_9 <- tax_glom(phso1_SAI_9, taxrank = "Family")
phso1_SAI_genus_9 <- tax_glom(phso1_SAI_9, taxrank = "Genus")


venn_phylum_SAI_9_prof <- ps_venn(phso1_SAI_phylum_9 , group = "Depth", fill = c("6" = "beige", "12" = "bisque3"), 
                                alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))

venn_familia_SAI_9_prof <- ps_venn(phso1_SAI_family_9 , group = "Depth", fill = c("6" = "beige", "12" = "bisque3"), 
                                alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))

venn_genero_SAI_9_prof <- ps_venn(phso1_SAI_genus_9 , group = "Depth", fill = c("6" = "beige", "12" = "bisque3"), 
                               alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))

venn_SAI_prof <- ggarrange(venn_genero_SAI_2_prof, venn_genero_SAI_5_prof, venn_genero_SAI_7_prof, venn_genero_SAI_9_prof,
                           ncol = 4, nrow = 1)
venn_SAI_prof

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Diagarama_venn_SAI_Genero_Profundidad.png", 
       venn_SAI_prof, width = 20, height = 6, dpi = 600)

venn_SAI_prof <- ggarrange(venn_familia_SAI_2_prof, venn_familia_SAI_5_prof, venn_familia_SAI_7_prof, venn_familia_SAI_9_prof,
                           ncol = 4, nrow = 1)
venn_SAI_prof

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Diagarama_venn_SAI_Familia_Profundidad.png", 
       venn_SAI_prof, width = 20, height = 6, dpi = 600)

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Diagrama_venn_SAI_Genero_Profundidad_Mes_2.png", 
       venn_genero_SAI_2_prof, width = 2.5, height = 2.5, dpi = 600)

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Diagrama_venn_SAI_Genero_Profundidad_Mes_6.png", 
       venn_genero_SAI_5_prof, width = 2.5, height = 2.5, dpi = 600)

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Diagrama_venn_SAI_Genero_Profundidad_Mes_8.png", 
       venn_genero_SAI_7_prof, width = 2.5, height = 2.5, dpi = 600)

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Diagrama_venn_SAI_Genero_Profundidad_Mes_10.png", 
       venn_genero_SAI_9_prof, width = 2.5, height = 2.5, dpi = 600)


#No hay muestras de 12 metros a los 2 meses de GOR

phso1_GOR_phylum_5 <- tax_glom(phso1_GOR_5, taxrank = "Phylum")
phso1_GOR_family_5 <- tax_glom(phso1_GOR_5, taxrank = "Family")
phso1_GOR_genus_5 <- tax_glom(phso1_GOR_5, taxrank = "Genus")


venn_phylum_GOR_5_prof <- ps_venn(phso1_GOR_phylum_5 , group = "Depth", fill = c("6" = "beige", "12" = "bisque3"), 
                                alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))

venn_familia_GOR_5_prof <- ps_venn(phso1_GOR_family_5 , group = "Depth", fill = c("6" = "beige", "12" = "bisque3"), 
                                alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))

venn_genero_GOR_5_prof <- ps_venn(phso1_GOR_genus_5 , group = "Depth", fill = c("6" = "beige", "12" = "bisque3"), 
                               alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))


phso1_GOR_phylum_7 <- tax_glom(phso1_GOR_7, taxrank = "Phylum")
phso1_GOR_family_7 <- tax_glom(phso1_GOR_7, taxrank = "Family")
phso1_GOR_genus_7 <- tax_glom(phso1_GOR_7, taxrank = "Genus")


venn_phylum_GOR_7_prof <- ps_venn(phso1_GOR_phylum_7 , group = "Depth", fill = c("6" = "beige", "12" = "bisque3"), 
                                alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))

venn_familia_GOR_7_prof <- ps_venn(phso1_GOR_family_7 , group = "Depth", fill = c("6" = "beige", "12" = "bisque3"), 
                                alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))

venn_genero_GOR_7_prof <- ps_venn(phso1_GOR_genus_7 , group = "Depth", fill = c("6" = "beige", "12" = "bisque3"), 
                               alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))


phso1_GOR_phylum_9 <- tax_glom(phso1_GOR_9, taxrank = "Phylum")
phso1_GOR_family_9 <- tax_glom(phso1_GOR_9, taxrank = "Family")
phso1_GOR_genus_9 <- tax_glom(phso1_GOR_9, taxrank = "Genus")


venn_phylum_GOR_9_prof <- ps_venn(phso1_GOR_phylum_9 , group = "Depth", fill = c("6" = "beige", "12" = "bisque3"), 
                                alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))

venn_familia_GOR_9_prof <- ps_venn(phso1_GOR_family_9 , group = "Depth", fill = c("6" = "beige", "12" = "bisque3"), 
                                alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))

venn_genero_GOR_9_prof <- ps_venn(phso1_GOR_genus_9 , group = "Depth", fill = c("6" = "beige", "12" = "bisque3"), 
                               alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))

venn_GOR_prof <- ggarrange(venn_genero_GOR_5_prof, venn_genero_GOR_7_prof, venn_genero_GOR_9_prof,
                           ncol = 3, nrow = 1)
venn_GOR_prof

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Diagarama_venn_GOR_Genero_Profundidad.png", 
       venn_GOR_prof, width = 20, height = 6, dpi = 600)

venn_GOR_prof <- ggarrange(venn_familia_SAI_5_prof, venn_familia_SAI_7_prof, venn_familia_SAI_9_prof,
                           ncol = 3, nrow = 1)
venn_GOR_prof

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Diagarama_venn_GOR_Familia_Profundidad.png", 
       venn_GOR_prof, width = 20, height = 6, dpi = 600)


ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Diagrama_venn_GOR_Genero_Profundidad_Mes_6.png", 
       venn_genero_GOR_5_prof, width = 2.5, height = 2.5, dpi = 600)

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Diagrama_venn_GOR_Genero_Profundidad_Mes_8.png", 
       venn_genero_GOR_7_prof, width = 2.5, height = 2.5, dpi = 600)

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Diagrama_venn_GOR_Genero_Profundidad_Mes_10.png", 
       venn_genero_GOR_9_prof, width = 2.5, height = 2.5, dpi = 600)


##################################################### DIVERSIDAD ALFA (SHANNON Y DOMINANCIA DE SIMPSON) #########################################################################
alpha_diversity <- estimate_richness(phso1, measures = c("Shannon", "Simpson", "Observed"))
colnames(metadata)
colnames(sample_data(phso1))
sample_data(phso1)$sampleID <- rownames(sample_data(phso1))

alpha_diversity <- alpha_diversity %>%
  rownames_to_column(var = "sampleID") %>%
  left_join(as.data.frame(sample_data(phso1)), by = "sampleID")

ggboxplot(alpha_diversity, x = "Location", y = "Shannon",
          color = "Location", add = "jitter", palette = "jco") +
  stat_compare_means(method = "kruskal.test", label = "p.signif") +
  labs(title = "Comparación de Diversidad Alfa (Shannon) por Location")

ggboxplot(alpha_diversity, x = "Location", y = "Simpson",
          color = "Location", add = "jitter", palette = "jco") +
  stat_compare_means(method = "kruskal.test", label = "p.signif") +
  labs(title = "Comparación de Diversidad Alfa (Simpson) por Location")

ggboxplot(subset(alpha_diversity, Location == "SAI"), x = "Month", y = "Shannon",
          color = "Month", add = "jitter", palette = "jco") +
  stat_compare_means(method = "kruskal.test", label = "p.signif") +
  labs(title = "Comparación de Diversidad Alfa (Shannon) en SAI por Mes")

ggboxplot(subset(alpha_diversity, Location == "GOR"), x = "Month", y = "Shannon",
          color = "Month", add = "jitter", palette = "jco") +
  stat_compare_means(method = "kruskal.test", label = "p.signif") +
  labs(title = "Comparación de Diversidad Alfa (Shannon) en GOR por Mes")

# VIOLINPLOTS
library(ggplot2)
library(dplyr)

alpha_diversity <- alpha_diversity %>%
  mutate(Month = factor(Month, levels = c("2", "6", "8", "10")))
alpha_diversity <- estimate_richness(phso1, measures = c("Shannon", "Simpson", "Observed"))

plot_alpha_diversity <- function(data, x_var, y_var, group_filter = NULL, title, palette_colors) {
  if (!is.null(group_filter)) {
    data <- data %>% filter(!!sym(names(group_filter)) == group_filter[[1]])
  }
  ggplot(data, aes_string(x = x_var, y = y_var, fill = x_var)) +
    geom_violin(trim = FALSE, alpha = 0.7, color = "black") +
    geom_boxplot(width = 0.1, fill = "white", position = position_dodge(width = 0.9)) +  # Boxplot
    scale_fill_manual(values = palette_colors) +
    labs(title = title, x = x_var, y = y_var) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 18),
      legend.position = "none",
      axis.line = element_line(color = "black")
    )
}

custom_palette_location <- c("SAI" = "darkseagreen", "GOR" = "cadetblue3")
alpha_diversity <- alpha_diversity %>%
  mutate(Month = factor(Month, levels = c("2", "6", "8", "10")))
custom_palette_month <- c("2" = "pink2", "6" = "palevioletred", "8" = "plum3", "10" = "paleturquoise3")

shannon <- plot_alpha_diversity(alpha_diversity, 
                             x_var = "Location", 
                             y_var = "Shannon", 
                             title = "Shannon", 
                             palette_colors = custom_palette_location)
shannon

simpson <- plot_alpha_diversity(alpha_diversity, 
                               x_var = "Location", 
                               y_var = "Simpson", 
                               title = "Simpson", 
                               palette_colors = custom_palette_location)
simpson

#MESES
shannon_SAI <- plot_alpha_diversity(alpha_diversity, 
                             x_var = "Month", 
                             y_var = "Shannon", 
                             group_filter = list("Location" = "SAI"),
                             title = "Shannon - Meses SAI") +
  scale_fill_manual(values = custom_palette_month)

shannon_SAI

simpson_SAI <- plot_alpha_diversity(alpha_diversity, 
                               x_var = "Month", 
                               y_var = "Simpson", 
                               group_filter = list("Location" = "SAI"),
                               title = "Simpson - Meses SAI", 
                               palette_colors = custom_palette_month)
simpson_SAI


shannon_GOR <- plot_alpha_diversity(alpha_diversity, 
                             x_var = "Month", 
                             y_var = "Shannon", 
                             group_filter = list("Location" = "GOR"),
                             title = "Shannon - Meses GOR", 
                             palette_colors = custom_palette_month)

shannon_GOR

simpson_GOR <- plot_alpha_diversity(alpha_diversity, 
                               x_var = "Month", 
                               y_var = "Simpson", 
                               group_filter = list("Location" = "GOR"),
                               title = "Simpson - Meses GOR", 
                               palette_colors = custom_palette_month)

simpson_GOR

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Violin_plot_Shannon.png", 
       shannon, width = 5, height = 4, dpi = 600)

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Violin_plot_Shannon_Meses_SAI.png", 
       shannon_SAI, width = 5, height = 4, dpi = 600)

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Violin_plot_Shannon_Meses_GOR.png", 
       shannon_GOR, width = 5, height = 4, dpi = 600)

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Violin_plot_Simpson.png", 
       simpson, width = 5, height = 4, dpi = 600)

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Violin_plot_Simpson_Meses_SAI.png", 
       simpson_SAI, width = 5, height = 4, dpi = 600)

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Violin_plot_Simpson_Meses_GOR.png", 
       simpson_GOR, width = 5, height = 4, dpi = 600)


# PROFUNDIDAD
alpha_diversity <- alpha_diversity %>%
  mutate(Depth = factor(Depth, levels = c("6", "12")))

plot_alpha_diversity_depth <- function(data, x_var, y_var, facet_var, group_filter = NULL, title, palette_colors) {
  if (!is.null(group_filter)) {
    data <- data %>% filter(!!sym(names(group_filter)) == group_filter[[1]])
  }
  
  ggplot(data, aes_string(x = x_var, y = y_var, fill = x_var)) +
    geom_violin(trim = FALSE, alpha = 0.7, color = "black") +  
    geom_boxplot(width = 0.1, fill = "white", position = position_dodge(width = 0.9)) +  
    scale_fill_manual(values = palette_colors) +
    facet_wrap(~get(facet_var)) + 
    labs(title = title, x = x_var, y = y_var) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 18),
      legend.position = "none",
      axis.line = element_line(color = "black")
    )
}

custom_palette_depth <- c("6" = "beige", "12" = "bisque3")

shannon_SAI_prof <- plot_alpha_diversity_depth(alpha_diversity, "Depth", "Shannon", "Month", 
                                          group_filter = list("Location" = "SAI"),
                                          title = "Shannon - SAI por Mes y Profundidad", 
                                          palette_colors = custom_palette_depth)

shannon_SAI_prof

simpson_SAI_prof <- plot_alpha_diversity_depth(alpha_diversity, "Depth", "Simpson", "Month", 
                                          group_filter = list("Location" = "SAI"),
                                          title = "Simpson - SAI por Mes y Profundidad", 
                                          palette_colors = custom_palette_depth)

simpson_SAI_prof


alpha_diversity_gor <- alpha_diversity %>%
  filter(!(Location == "GOR" & Month == "2" & Depth == "12"))

shannon_GOR_prof <- plot_alpha_diversity_depth(alpha_diversity_gor, "Depth", "Shannon", "Month", 
                                          group_filter = list("Location" = "GOR"),
                                          title = "Shannon - GOR por Mes y Profundidad", 
                                          palette_colors = custom_palette_depth)

shannon_GOR_prof

simpson_GOR_prof <- plot_alpha_diversity_depth(alpha_diversity_gor, "Depth", "Simpson", "Month", 
                                          group_filter = list("Location" = "GOR"),
                                          title = "Simpson - GOR por Mes y Profundidad", 
                                          palette_colors = custom_palette_depth)

simpson_GOR_prof


ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Violin_plot_Shannon_Profundidad_SAI.png", 
       shannon_SAI_prof, width = 5, height = 4, dpi = 600)

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Violin_plot_Simpson_Profundidad_SAI.png", 
       simpson_SAI_prof, width = 5, height = 4, dpi = 600)

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Violin_plot_Shannon_Profundidad_GOR.png", 
       shannon_GOR_prof, width = 5, height = 4, dpi = 600)

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/Violin_plot_Simpson_Profundidad_GOR.png", 
       simpson_GOR_prof, width = 5, height = 4, dpi = 600)

################################################## COMPARACIONES ESTADÍSTICAS PARA DIVERSIDAD ALFA ##################################################################
library(FSA)
# 1. Comparaciones para Location == "SAI"
alpha_sai <- filter(alpha_diversity, Location == "SAI")

#Test de Kruskal-Wallis para Month dentro de SAI
kruskal_results_sai <- kruskal.test(Shannon ~ Month, data = alpha_sai)
cat("Resultados de Kruskal-Wallis para Diversidad Alfa en SAI por Mes:\n")
print(kruskal_results_sai)

#Comparaciones post-hoc entre los meses en SAI
dunn_results_sai <- dunnTest(Shannon ~ Month, data = alpha_sai, method = "bonferroni")
print(dunn_results_sai)

#Test de Kruskal-Wallis para Month dentro de SAI
kruskal_results_sai <- kruskal.test(Simpson ~ Month, data = alpha_sai)
cat("Resultados de Kruskal-Wallis para Diversidad Alfa en SAI por Mes:\n")
print(kruskal_results_sai)

#Comparaciones post-hoc entre los meses en SAI
dunn_results_sai <- dunnTest(Simpson ~ Month, data = alpha_sai, method = "bonferroni")
print(dunn_results_sai)

# 2. Comparaciones para Location == "GOR"
alpha_gor <- filter(alpha_diversity, Location == "GOR")

#Test de Kruskal-Wallis para Month dentro de GOR
kruskal_results_gor <- kruskal.test(Shannon ~ Month, data = alpha_gor)
print(kruskal_results_gor)

#Comparaciones post-hoc entre los meses en GOR
dunn_results_gor <- dunnTest(Shannon ~ Month, data = alpha_gor, method = "bonferroni")
print(dunn_results_gor)


#Test de Kruskal-Wallis para Month dentro de GOR
kruskal_results_gor <- kruskal.test(Simpson ~ Month, data = alpha_gor)
print(kruskal_results_gor)

#Comparaciones post-hoc entre los meses en GOR
dunn_results_gor <- dunnTest(Simpson ~ Month, data = alpha_gor, method = "bonferroni")
print(dunn_results_gor)


# Comparar diversidad alfa por Ubicacion
wilcox_results_location <- wilcox.test(Shannon ~ Location, data = alpha_diversity)
cat("Resultados de Wilcoxon para Diversidad Alfa por Ubicacion:\n")
print(wilcox_results_location)

wilcox_results_location <- wilcox.test(Simpson ~ Location, data = alpha_diversity)
cat("Resultados de Wilcoxon para Diversidad Alfa por Ubicacion:\n")
print(wilcox_results_location)


# Comparar diversidad alfa por "Depth" dentro de Location == "SAI"
wilcox_results_sai_depth <- wilcox.test(Shannon ~ Depth, data = subset(alpha_diversity, Location == "SAI"))
cat("Resultados de Wilcoxon para Diversidad Alfa en SAI por Profundidad:\n")
print(wilcox_results_sai_depth)

wilcox_results_sai_depth <- wilcox.test(Simpson ~ Depth, data = subset(alpha_diversity, Location == "SAI"))
cat("Resultados de Wilcoxon para Diversidad Alfa en SAI por Profundidad:\n")
print(wilcox_results_sai_depth)

# Comparar diversidad alfa por "Depth" dentro de Location == "GOR"
wilcox_results_gor_depth <- wilcox.test(Shannon ~ Depth, data = subset(alpha_diversity, Location == "GOR"))
cat("Resultados de Wilcoxon para Diversidad Alfa en GOR por Profundidad:\n")
print(wilcox_results_gor_depth)

wilcox_results_gor_depth <- wilcox.test(Simpson ~ Depth, data = subset(alpha_diversity, Location == "GOR"))
cat("Resultados de Wilcoxon para Diversidad Alfa en GOR por Profundidad:\n")
print(wilcox_results_gor_depth)


compare_depths <- function(data, location, month, index) {
  filtered_data <- data %>%
    filter(Location == location & Month == month)
  wilcox_test <- wilcox.test(filtered_data[[index]] ~ Depth, data = filtered_data)
  cat("Comparación de", index, "entre profundidades en", location, "para el mes", month, ":\n")
  print(wilcox_test)
  cat("\n")
}

months <- c("2", "6", "8", "10")
locations <- c("SAI", "GOR")

for (location in locations) {
  for (month in months) {
    if (location == "GOR" & month == "2") {
      cat("No hay datos de 12m para GOR en el mes 2. Se omite la comparación.\n")
      next
    }
    compare_depths(alpha_diversity, location, month, "Shannon")
    compare_depths(alpha_diversity, location, month, "Simpson")
  }
}

################################################## ANALISIS DE DIVERSIDAD BETA ######################################################################
library(phyloseq)
library(ggplot2)
library(vegan)

sample_data(phso1)$Month <- factor(sample_data(phso1)$Month, levels = c("2", "6", "8", "10"))                            
phso1_SAI <- subset_samples(phso1, Location == "SAI")
phso1_GOR <- subset_samples(phso1, Location == "GOR")
sample_data(phso1_SAI)$Month <- factor(sample_data(phso1_SAI)$Month, levels = c("2", "6", "8", "10"))
sample_data(phso1_GOR)$Month <- factor(sample_data(phso1_GOR)$Month, levels = c("2", "6", "8", "10"))
dist_bray <- phyloseq::distance(phso1, method = "bray")
pcoa <- ordinate(phso1, method = "PCoA", distance = dist_bray)

pcoa_location <- plot_ordination(phso1, pcoa, color = "Location", shape = "Location") +
  stat_ellipse(aes(group = Location), linetype = 2) +
  geom_point(size = 3, alpha = 0.8) + 
  geom_text(aes(label = Location), hjust = 1.5, vjust = 1.5, size = 3.5) +
  scale_color_manual(values = c("SAI" = "darkseagreen", "GOR" = "cadetblue3")) +  
  labs(title = "Diversidad Beta por Ubicación") + 
  labs(x = "PCo1 [11.5%]", y = "PCo2 [7.8%]") +
  ggtitle("ADONIS: p= 0.001, R2= 0.1054") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0, size = 16),  
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18), 
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 18), 
    panel.grid = element_blank() 
  )

pcoa_location

pcoa_sai_month <- plot_ordination(phso1_SAI, pcoa, color = "Month") +
  stat_ellipse(aes(group = Month), linetype = 2) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(aes(label = Month), hjust = 1.5, vjust = 1.5, size = 3) + 
  scale_color_manual(values = c("2" = "pink2", "6" = "palevioletred", "8" = "plum3", "10" = "paleturquoise3")) +  
  labs(title = "Diversidad Beta en SAI por Mes") +
  labs(x = "PCo1 [11.5%]", y = "PCo2 [7.8%]") +
  ggtitle("ADONIS: p= 0.001, R2= 0.14969") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0, size = 16),
    axis.title = element_text(size = 18), 
    axis.text = element_text(size = 18), 
    legend.text = element_text(size = 14), 
    legend.title = element_text(size = 18), 
    panel.grid = element_blank()
  )
pcoa_sai_month

pcoa_gor_month <- plot_ordination(phso1_GOR, pcoa, color = "Month") +
  stat_ellipse(aes(group = Month), linetype = 2) +
  geom_text(aes(label = Month), hjust = 1.5, vjust = 1.5, size = 3) +  
  scale_color_manual(values = c("2" = "pink2", "6" = "palevioletred", "8" = "plum3", "10" = "paleturquoise3")) +
  labs(title = "Diversidad Beta en GOR por Mes") +
  labs(x = "PCo1 [11.5%]", y = "PCo2 [7.8%]") +
  ggtitle("ADONIS: p= 0.001, R2= 0.23031") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0, size = 16),
    axis.title = element_text(size = 18), 
    axis.text = element_text(size = 18), 
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 18), 
    panel.grid = element_blank() 
  )
pcoa_gor_month

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/PCoA_Ubicacion.png", 
       pcoa_location, width = 8, height = 6, dpi = 600)

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/PCoA_SAI_Meses.png", 
       pcoa_sai_month, width = 8, height = 6, dpi = 600)

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/PCoA_GOR_Meses.png", 
       pcoa_gor_month, width = 8, height = 6, dpi = 600)


sample_data(phso1)$Depth <- factor(sample_data(phso1)$Depth, levels = c("6", "12"))
phso1_SAI <- subset_samples(phso1, Location == "SAI")
phso1_GOR <- subset_samples(phso1, Location == "GOR")
sample_data(phso1_SAI)$Depth <- factor(sample_data(phso1_SAI)$Depth, levels = c("6", "12"))
sample_data(phso1_GOR)$Depth <- factor(sample_data(phso1_GOR)$Depth, levels = c("6", "12"))
dist_bray <- phyloseq::distance(phso1, method = "bray")
pcoa <- ordinate(phso1, method = "PCoA", distance = dist_bray)

pcoa_sai_depth <- plot_ordination(phso1_SAI, pcoa, color = "Depth") +
  stat_ellipse(aes(group = Depth), linetype = 2) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(aes(label = Depth), hjust = 1.5, vjust = 1.5, size = 3) +
  scale_color_manual(values = c("6" = "bisque2", "12" = "bisque3")) +
  labs(title = "Diversidad Beta en SAI por Profundidad") +
  labs(x = "PCo1 [11.5%]", y = "PCo2 [7.8%]") +
  ggtitle("ADONIS: p= 0.001, R2= 0.09325") +
  theme_classic()
pcoa_sai_depth

pcoa_gor_depth <- plot_ordination(phso1_GOR, pcoa, color = "Depth") +
  stat_ellipse(aes(group = Depth), linetype = 2) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(aes(label = Depth), hjust = 1.5, vjust = 1.5, size = 3) +
  scale_color_manual(values = c("6" = "bisque2", "12" = "bisque3")) +
  labs(title = "Diversidad Beta en GOR por Profundidad") +
  labs(x = "PCo1 [11.5%]", y = "PCo2 [7.8%]") +
  ggtitle("ADONIS: p= 0.007, R2= 0.07649") +
  theme_classic()
pcoa_gor_depth

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/PCoA_SAI_Profundidad.png", 
       pcoa_sai_depth, width = 8, height = 6, dpi = 600)

ggsave("/Users/ktari/OneDrive/Escritorio/tesis16S/PCoA_GOR_Profundidad.png", 
       pcoa_gor_depth, width = 8, height = 6, dpi = 600)


################################################## COMPARACIONES ESTADISTICAS DIVERSIDAD BETA #####################################################################################
# Comparar diversidad beta entre ubicaciones usando PERMANOVA
adonis_location <- adonis2(dist_bray ~ Location, data = as(sample_data(phso1), "data.frame"))
cat("Resultados de PERMANOVA para Diversidad Beta por Ubicación:\n")
print(adonis_location)

# Comparar diversidad beta entre meses en SAI usando PERMANOVA
adonis_sai_month <- adonis2(phyloseq::distance(phso1_SAI, method = "bray") ~ Month, 
                            data = as(sample_data(phso1_SAI), "data.frame"))
cat("Resultados de PERMANOVA para Diversidad Beta en SAI por Mes:\n")
print(adonis_sai_month)

# Comparar diversidad beta entre meses en GOR usando PERMANOVA
adonis_gor_month <- adonis2(phyloseq::distance(phso1_GOR, method = "bray") ~ Month, 
                            data = as(sample_data(phso1_GOR), "data.frame"))
cat("Resultados de PERMANOVA para Diversidad Beta en GOR por Mes:\n")
print(adonis_gor_month)

# Comparar diversidad beta entre profundidades en SAI usando PERMANOVA
adonis_sai_depth <- adonis2(phyloseq::distance(phso1_SAI, method = "bray") ~ Depth, 
                            data = as(sample_data(phso1_SAI), "data.frame"))
cat("Resultados de PERMANOVA para Diversidad Beta en SAI por Profundidad:\n")
print(adonis_sai_depth)

# Comparar diversidad beta entre profundidades en GOR usando PERMANOVA
adonis_gor_depth <- adonis2(phyloseq::distance(phso1_GOR, method = "bray") ~ Depth, 
                            data = as(sample_data(phso1_GOR), "data.frame"))
cat("Resultados de PERMANOVA para Diversidad Beta en GOR por Profundidad:\n")
print(adonis_gor_depth)


# Comparar diversidad beta entre profundidades en SAI mes 2
adonis_sai_2 <- adonis2(phyloseq::distance(subset_samples(phso1_SAI, Month == 2), method = "bray") ~ Depth, 
                        data = as(sample_data(subset_samples(phso1_SAI, Month == 2)), "data.frame"))
cat("Resultados de PERMANOVA para Diversidad Beta en SAI  mes 2 por Profundidad:\n")
print(adonis_sai_2)

# Comparar diversidad beta entre profundidades en SAI  mes 6
adonis_sai_6 <- adonis2(phyloseq::distance(subset_samples(phso1_SAI, Month == 6), method = "bray") ~ Depth, 
                        data = as(sample_data(subset_samples(phso1_SAI, Month == 6)), "data.frame"))
cat("Resultados de PERMANOVA para Diversidad Beta en SAI mes 6 por Profundidad:\n")
print(adonis_sai_6)

# Comparar diversidad beta entre profundidades en SAI mes 8
adonis_sai_8 <- adonis2(phyloseq::distance(subset_samples(phso1_SAI, Month == 8), method = "bray") ~ Depth, 
                        data = as(sample_data(subset_samples(phso1_SAI, Month == 8)), "data.frame"))
cat("Resultados de PERMANOVA para Diversidad Beta en SAI mes 8 por Profundidad:\n")
print(adonis_sai_8)

# Comparar diversidad beta entre profundidades en SAI mes 10
adonis_sai_10 <- adonis2(phyloseq::distance(subset_samples(phso1_SAI, Month == 10), method = "bray") ~ Depth, 
                         data = as(sample_data(subset_samples(phso1_SAI, Month == 10)), "data.frame"))
cat("Resultados de PERMANOVA para Diversidad Beta en SAI mes 10 por Profundidad:\n")
print(adonis_sai_10)



# Comparar diversidad beta entre profundidades en GOR  mes 6
adonis_gor_6 <- adonis2(phyloseq::distance(subset_samples(phso1_GOR, Month == 6), method = "bray") ~ Depth, 
                        data = as(sample_data(subset_samples(phso1_GOR, Month == 6)), "data.frame"))
cat("Resultados de PERMANOVA para Diversidad Beta en GOR  mes 6 por Profundidad:\n")
print(adonis_gor_6)

# Comparar diversidad beta entre profundidades en GOR  mes 8
adonis_gor_8 <- adonis2(phyloseq::distance(subset_samples(phso1_GOR, Month == 8), method = "bray") ~ Depth, 
                        data = as(sample_data(subset_samples(phso1_GOR, Month == 8)), "data.frame"))
cat("Resultados de PERMANOVA para Diversidad Beta en GOR mes 8 por Profundidad:\n")
print(adonis_gor_8)

# Comparar diversidad beta entre profundidades en GOR mes 10
adonis_gor_10 <- adonis2(phyloseq::distance(subset_samples(phso1_GOR, Month == 10), method = "bray") ~ Depth, 
                         data = as(sample_data(subset_samples(phso1_GOR, Month == 10)), "data.frame"))
cat("Resultados de PERMANOVA para Diversidad Beta en GOR mes 10 por Profundidad:\n")
print(adonis_gor_10)

######################################################## ANCOM-BC ####################################################################333
library(ANCOMBC)
library(phyloseq)
library(ggplot2)
library(openxlsx)
library(dplyr)
library(stringr)
library(vegan)
library(ggplot2)
library(RColorBrewer)

# ANCOM-BC A NIVEL DE FAMILIA
phso1.family <- tax_glom(phso1, taxrank = 'Family', NArm = TRUE)
phso1.taxa.family <- subset_samples(phso1.family, Location %in% c("SAI", "GOR"))

out <- ancombc(data = phso1.taxa.family, 
               assay_name = "counts", 
               tax_level = "Family", 
               formula = "Location", 
               p_adj_method = "holm", 
               prv_cut = 0.10, 
               lib_cut = 1000, 
               group = "Location", 
               struc_zero = TRUE, 
               neg_lb = TRUE, 
               tol = 1e-5, 
               max_iter = 100, 
               conserve = TRUE, 
               alpha = 0.05, 
               global = TRUE,
               n_cl = 1, 
               verbose = TRUE)

res <- out$res
res_global = out$res_global

tab_lfc = res$lfc
col_name = colnames(tab_lfc)
colnames(tab_lfc) = col_name
tab_lfc %>% 
  datatable(caption = "Log Fold Changes from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_se = res$se
colnames(tab_se) = col_name
tab_se %>% 
  datatable(caption = "SEs from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_w = res$W
colnames(tab_w) = col_name
tab_w %>% 
  datatable(caption = "Test Statistics from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_p = res$p_val
colnames(tab_p) = col_name
tab_p %>% 
  datatable(caption = "P-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_q = res$q
colnames(tab_q) = col_name
tab_q %>% 
  datatable(caption = "Adjusted p-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_diff = res$diff_abn
colnames(tab_diff) = col_name
tab_diff %>% 
  datatable(caption = "Differentially Abundant Taxa from the Primary Result")

# Que bacterias contribuyen mas a las diferencias en la composicion de las comunidades entre los grupos (beta-diversidad)
#TODO: hacer un analisis PERMANOVA
pseq <- phso1
pseq.rel <- microbiome::transform(pseq, "compositional")
otu <- abundances(pseq.rel)
meta <- meta(pseq.rel)
permanova <- adonis(t(otu) ~ Location,
                    data = meta, permutations=999, method = "bray")
print(as.data.frame(permanova$aov.tab)["Location", "Pr(>F)"])
dist <- vegdist(t(otu))
anova(betadisper(dist, meta$Location))
permutest(betadisper(dist, meta$Location), pairwise = TRUE)
coef <- coefficients(permanova)["Location1",]
top.coef <- coef[rev(order(abs(coef)))[1:40]]
par(mar = c(3, 6, 3, 6))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")
tab_p$LocationSAI <- as.numeric(tab_p$LocationSAI[, "structural_zero (Location = SAI)"])
tab_lfc$taxon <- sapply(strsplit(tab_lfc$taxon, "_"), function(x) x[length(x)])
tab_lfc

combined_data <- data.frame(
  Family = tab_lfc$taxon,  # Usar la columna taxon de tab_lfc
  LFC_LocationSAI = tab_lfc$LocationSAI,  # Log Fold Change para LocationSAI
  p_value_LocationSAI = tab_p$LocationSAI  # Valores p para LocationSAI
)

head(tab_lfc)
head(tab_p)
str(tab_lfc)
str(tab_p)

significant_taxa <- combined_data[combined_data$p_value_LocationSAI < 0.05, ]
significant_taxa

p_lfc_significant <- ggplot(significant_taxa, aes(x = reorder(Family, LFC_LocationSAI), y = LFC_LocationSAI, fill = Family)) +
  geom_bar(stat = "identity", width = 0.99) +  # Aumentar el grosor de las barras
  coord_flip() +
  labs(
    x = "Familia",
    y = "Log Fold Change (LFC)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 11),
    legend.position = "none",
    plot.margin = margin(t = 40, r = 20, b = 20, l = 20),
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5),
  ) + 
  scale_x_discrete(expand = expansion(mult = c(0, 0.2))) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(9, "Set3"))(length(significant_taxa$Family))) +
  annotate("segment", y = -max(abs(significant_taxa$LFC_LocationSAI)), yend = max(abs(significant_taxa$LFC_LocationSAI)), x = nrow(significant_taxa) + 1, xend = nrow(significant_taxa) + 1, color = "black", size = 0.5, arrow = arrow(ends = "both", length = unit(0.2, "cm"))) +  # Flecha
  annotate("text", x = nrow(significant_taxa) + 2, y = -max(abs(significant_taxa$LFC_LocationSAI)) * 0.9, label = "GOR", color = "black", size = 3) +  # Etiqueta SAI
  annotate("text", x = nrow(significant_taxa) + 2, y = max(abs(significant_taxa$LFC_LocationSAI)) * 0.9, label = "SAI", color = "black", size = 3)  # Etiqueta GOR

print(p_lfc_significant)

output_dir <- "/Users/ktari/OneDrive/Escritorio/tesis16S/ANCOMBC"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

ggsave(file.path(output_dir, "Top_Tax_LogFC_a_Family_estrecho.png"), p_lfc_significant, width = 7, height = 12, dpi = 300)
ggsave(file.path(output_dir, "Top_Tax_LogFC_a_Family.png"), p_lfc_significant, width = 11, height = 10, dpi = 300)


# ANCOM-BC A NIVEL DE GENERO
phso1.genus <- tax_glom(phso1, taxrank = 'Genus', NArm = TRUE)
phso1.taxa.genus <- subset_samples(phso1.genus, Location %in% c("SAI", "GOR"))

out <- ancombc(data = phso1.taxa.genus, 
               assay_name = "counts", 
               tax_level = "Genus", 
               formula = "Location", 
               p_adj_method = "holm", 
               prv_cut = 0.10, 
               lib_cut = 1000, 
               group = "Location", 
               struc_zero = TRUE, 
               neg_lb = TRUE, 
               tol = 1e-5, 
               max_iter = 100, 
               conserve = TRUE, 
               alpha = 0.05, 
               global = TRUE,
               n_cl = 1, 
               verbose = TRUE)

res <- out$res
res_global = out$res_global

tab_lfc = res$lfc
col_name = colnames(tab_lfc)
colnames(tab_lfc) = col_name
tab_lfc %>% 
  datatable(caption = "Log Fold Changes from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_se = res$se
colnames(tab_se) = col_name
tab_se %>% 
  datatable(caption = "SEs from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_w = res$W
colnames(tab_w) = col_name
tab_w %>% 
  datatable(caption = "Test Statistics from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_p = res$p_val
colnames(tab_p) = col_name
tab_p %>% 
  datatable(caption = "P-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_q = res$q
colnames(tab_q) = col_name
tab_q %>% 
  datatable(caption = "Adjusted p-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_diff = res$diff_abn
colnames(tab_diff) = col_name
tab_diff %>% 
  datatable(caption = "Differentially Abundant Taxa from the Primary Result")
pseq <- phso1
pseq.rel <- microbiome::transform(pseq, "compositional")
otu <- abundances(pseq.rel)
meta <- meta(pseq.rel)
permanova <- adonis(t(otu) ~ Location,
                    data = meta, permutations=999, method = "bray")

print(as.data.frame(permanova$aov.tab)["Location", "Pr(>F)"])
dist <- vegdist(t(otu))
anova(betadisper(dist, meta$Location))
permutest(betadisper(dist, meta$Location), pairwise = TRUE)
coef <- coefficients(permanova)["Location1",]
top.coef <- coef[rev(order(abs(coef)))[1:40]]
par(mar = c(3, 6, 3, 6))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")
tab_p$LocationSAI <- as.numeric(tab_p$LocationSAI[, "structural_zero (Location = SAI)"])
tab_lfc$taxon <- sapply(strsplit(tab_lfc$taxon, "_"), function(x) x[length(x)])
tab_lfc
combined_data <- data.frame(
  Genus = tab_lfc$taxon, 
  LFC_LocationSAI = tab_lfc$LocationSAI,  
  p_value_LocationSAI = tab_p$LocationSAI  
)

head(tab_lfc)
head(tab_p)
str(tab_lfc)
str(tab_p)

significant_taxa <- combined_data[combined_data$p_value_LocationSAI < 0.05, ]
significant_taxa

p_lfc_significant_genus <- ggplot(significant_taxa, aes(x = reorder(Genus, LFC_LocationSAI), y = LFC_LocationSAI, fill = Genus)) +
  geom_bar(stat = "identity", width = 0.99) +
  coord_flip() +
  labs(
    x = "Genero",
    y = "Log Fold Change (LFC)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9),
    legend.position = "none",
    plot.margin = margin(t = 40, r = 20, b = 20, l = 20),
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5),
  ) + 
  scale_x_discrete(expand = expansion(mult = c(0, 0.2))) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(9, "Set3"))(length(significant_taxa$Genus))) +
  annotate("segment", y = -max(abs(significant_taxa$LFC_LocationSAI)), yend = max(abs(significant_taxa$LFC_LocationSAI)), x = nrow(significant_taxa) + 1, xend = nrow(significant_taxa) + 1, color = "black", size = 0.5, arrow = arrow(ends = "both", length = unit(0.2, "cm"))) +  # Flecha
  annotate("text", x = nrow(significant_taxa) + 2, y = -max(abs(significant_taxa$LFC_LocationSAI)) * 0.9, label = "GOR", color = "black", size = 3) +  # Etiqueta SAI
  annotate("text", x = nrow(significant_taxa) + 2, y = max(abs(significant_taxa$LFC_LocationSAI)) * 0.9, label = "SAI", color = "black", size = 3)  # Etiqueta GOR

print(p_lfc_significant_genus)

output_dir <- "/Users/ktari/OneDrive/Escritorio/tesis16S/ANCOMBC"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

ggsave(file.path(output_dir, "Top_Tax_LogFC_a_Genus_estrecho.png"), p_lfc_significant_genus, width = 8, height = 12, dpi = 300)

ggsave(file.path(output_dir, "Top_Tax_LogFC_a_Genus.png"), p_lfc_significant_genus, width = 11, height = 8, dpi = 300)

