## ----echo=FALSE, results='hide', warning=FALSE, message=FALSE-----------------
library(scTEI)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(Seurat))
suppressWarnings(suppressPackageStartupMessages(library(SeuratData)))
suppressPackageStartupMessages(library(monocle3))

## ----warning=FALSE, message=FALSE---------------------------------------------
# load Packer and Zhu et al (2019) data set
SeuratData::InstallData("celegans.embryo.SeuratData")
library(celegans.embryo.SeuratData)
data("celegans.embryo")
celegans <- UpdateSeuratObject(celegans.embryo)
#celegans <- SeuratData::LoadData("celegans.embryo")
dim(celegans)
head(rownames(celegans))

## ----warning=FALSE, message=FALSE---------------------------------------------
# preprocess scRNA
all.genes <- rownames(celegans)
celegans <- Seurat::NormalizeData(
    celegans,
    normalization.method = "LogNormalize",
    scale.factor = 10000) |>
    Seurat::FindVariableFeatures(selection.method = "vst",
    nfeatures = 2000) |>
    Seurat::ScaleData(features = all.genes) |>
    Seurat::RunPCA(dims=50) |>
    Seurat::RunUMAP(dims = 1:10)

## ----fig.width=7, fig.align="center"------------------------------------------
Seurat::Idents(celegans) <- "embryo.time.bin"
p1 <- DimPlot(celegans)
print(p1)

## ----warning=FALSE, message=FALSE---------------------------------------------
# overwrite NA in cell.type + add embryo.time.bin x cell.type
celegans@meta.data$cell.type[is.na(celegans@meta.data$cell.type)] <- "notClassified"
celegans@meta.data["embryo.time.bin.cell.type"] <- paste0(
    unlist(celegans@meta.data["embryo.time.bin"]),
    "-",
    unlist(celegans@meta.data["cell.type"]))

# load Caenorhabditis elegans gene age estimation
celegans_ps <- readr::read_tsv(file = system.file("extdata",
    "Sun2021_Orthomap.tsv", package = "scTEI"))
table(celegans_ps$Phylostratum)

# define Phylostratum
ps_vec <- setNames(as.numeric(celegans_ps$Phylostratum),
    celegans_ps$GeneID)

## -----------------------------------------------------------------------------
# add TEI values
celegans@meta.data["TEI"] <- TEI(
    ExpressionSet = GetAssayData(celegans, assay="RNA", layer="counts"),
    Phylostratum = ps_vec
)

## ----fig.width=7, fig.align="center"------------------------------------------
# make FeaturePlot
p2 <- FeaturePlot(
    object = celegans,
    features = "TEI",
    min.cutoff='q05',
    max.cutoff='q95',
    cols = viridis(3))
print(p2)

## ----fig.width=7, fig.align="center"------------------------------------------
# make RidgePlot
p3 <- RidgePlot(object = celegans,
    features = "TEI",
    group.by = "embryo.time.bin")
print(p3)

## ----fig.width=7, fig.align="center"------------------------------------------
# make RidgePlot by cell type
Seurat::Idents(celegans) <- "cell.type"
p4 <- RidgePlot(object = celegans,
    features = "TEI",
    group.by = "cell.type") +
    Seurat::NoLegend()
print(p4)

## ----fig.width=7, fig.align="center"------------------------------------------
# subset to specific cell type - ADF + notClassified
ADF <- subset(celegans,
    cells = c(grep("ADF",
        celegans@meta.data$cell.type),
    grep("notClassified",
         celegans@meta.data$cell.type)))
p5 <- DimPlot(object = ADF)
Seurat::Idents(ADF) <- "embryo.time.bin"
p6 <- DimPlot(object = ADF)
p7 <- RidgePlot(ADF, "TEI")
Seurat::Idents(ADF) <- "embryo.time.bin.cell.type"
p8 <- RidgePlot(object = ADF, features = "TEI") +
    Seurat::NoLegend()

# make grid plot
print(plot_grid(p5, p6))
print(plot_grid(p7, p8))

## -----------------------------------------------------------------------------
# use pMatrix as data to cluster

# Seurat v4

#celegans.TEI <- Seurat::CreateSeuratObject(counts = celegans@assays$RNA@counts,
#    meta.data = celegans@meta.data, row.names = rownames(celegans@assays$RNA@counts))
#celegans.TEI@assays$RNA@data <- pMatrixTEI(
#    ExpressionSet = celegans.TEI@assays$RNA@counts,
#    Phylostratum = ps_vec
#)

# Seurat v5

celegans.TEI <- Seurat::CreateSeuratObject(
    counts = GetAssayData(celegans, assay="RNA", layer="counts"),
    meta.data = celegans@meta.data)
celegans.TEI <- SetAssayData(
    object = celegans.TEI,
    layer = "data",
    new.data = pMatrixTEI(
        ExpressionSet = GetAssayData(celegans.TEI, assay="RNA", layer="counts"),
        Phylostratum = ps_vec
    )
)
all.genes <- rownames(GetAssayData(celegans.TEI, assay="RNA", layer="data"))
celegans.TEI <- Seurat::FindVariableFeatures(
    celegans.TEI,
    selection.method = "vst",
    nfeatures = 2000) %>%
    Seurat::ScaleData(do.scale = FALSE, do.center = FALSE,
    features = all.genes) %>%
    Seurat::RunPCA(dims=50) %>%
    Seurat::RunUMAP(dims = 1:20)

## ----fig.width=7, fig.align="center"------------------------------------------
Seurat::Idents(celegans.TEI) <- "embryo.time.bin"
p9 <- DimPlot(celegans.TEI)
print(p9)

## ----fig.width=7, fig.align="center"------------------------------------------
# make FeaturePlot
p10 <- FeaturePlot(
    object = celegans.TEI,
    features = "TEI",
    min.cutoff='q05',
    max.cutoff='q95',
    cols = viridis(3))
print(p10)

## ----fig.width=7, fig.align="center"------------------------------------------
Seurat::Idents(celegans) <- "cell.type"
p11 <- DimPlot(celegans)
Seurat::Idents(celegans.TEI) <- "cell.type"
p12 <- DimPlot(celegans.TEI)
# make grid plot
print(plot_grid(p2, p11))
print(plot_grid(p10, p12))

## -----------------------------------------------------------------------------
# get TEI per strata
pS <- pStrataTEI(
    ExpressionSet = GetAssayData(celegans, assay="RNA", layer="counts"),
    Phylostratum = ps_vec
)

## -----------------------------------------------------------------------------
# get permutations

# Seurat v4

#bM <- bootTEI(
#    ExpressionSet = celegans@assays$RNA@counts,
#    Phylostratum = ps_vec,
#    permutations = 100
#)

# Seurat v5

bM <- bootTEI(
    ExpressionSet = GetAssayData(celegans, assay="RNA", layer="counts"),
    Phylostratum = ps_vec,
    permutations = 100
)

## -----------------------------------------------------------------------------
# get mean expression matrix

# Seurat v4

#meanMatrix <- REMatrix(
#    ExpressionSet = celegans@assays$RNA@data,
#    Phylostratum = ps_vec
#)

# Seurat v5

meanMatrix <- REMatrix(
    ExpressionSet = GetAssayData(celegans, assay="RNA", layer="data"),
    Phylostratum = ps_vec
)

## ----fig.width=7, fig.align="center"------------------------------------------
# get mean expression matrix with groups

cell_groups <- setNames(
    lapply(names(table(celegans@meta.data$cell.type)),
        function(x){which(celegans@meta.data$cell.type==x)}),
    names(table(celegans@meta.data$cell.type))
)

# Seurat v4

#meanMatrix_by_cell.type <- REMatrix(
#    ExpressionSet = celegans@assays$RNA@scale.data,
#    Phylostratum = ps_vec,
#    groups = cell_groups
#)

# Seurat v5

meanMatrix_by_cell.type <- REMatrix(
    ExpressionSet = GetAssayData(celegans, assay="RNA", layer="scale.data"),
    Phylostratum = ps_vec,
    groups = cell_groups
)

ComplexHeatmap::Heatmap(meanMatrix_by_cell.type,
    cluster_rows = FALSE, cluster_columns = FALSE,
    col = viridis::viridis(3))

## ----fig.width=7, fig.align="center"------------------------------------------
# computing relative expression profile over cell types

cell_groups <- setNames(
    lapply(names(table(celegans@meta.data$cell.type)),
        function(x){which(celegans@meta.data$cell.type==x)}),
    names(table(celegans@meta.data$cell.type))
)

# Seurat v4

#reMatrix_by_cell.type <- REMatrix(
#    ExpressionSet = celegans@assays$RNA@scale.data,
#    Phylostratum = ps_vec,
#    groups = cell_groups,
#    by = "row"
#)

# Seurat v5

reMatrix_by_cell.type <- REMatrix(
    ExpressionSet = GetAssayData(celegans, assay="RNA", layer="scale.data"),
    Phylostratum = ps_vec,
    groups = cell_groups,
    by = "row"
)

ComplexHeatmap::Heatmap(reMatrix_by_cell.type,
    cluster_rows = FALSE)

## -----------------------------------------------------------------------------
# load Packer and Zhu et al (2019) data set

#expression_matrix <- readRDS(
#    url(
#    paste0("http://staff.washington.edu/hpliner/data/",
#    "packer_embryo_expression.rds")
#    )
#)

#cell_metadata <- readRDS(
#    url(
#    paste0("http://staff.washington.edu/hpliner/data/",
#    "packer_embryo_colData.rds")
#    )
#)

#gene_annotation <- readRDS(
#    url(
#    paste0("http://staff.washington.edu/hpliner/data/",
#    "packer_embryo_rowData.rds")
#    )
#)

packer_embryo_expression_path <- system.file("extdata",
    "packer_embryo_expression.rds",
    package = "scTEI")
expression_matrix <- readRDS(packer_embryo_expression_path)

packer_embryo_colData_path <- system.file("extdata",
    "packer_embryo_colData.rds",
    package = "scTEI")
cell_metadata <- readRDS(packer_embryo_colData_path)

packer_embryo_rowData_path <- system.file("extdata",
    "packer_embryo_rowData.rds",
    package = "scTEI")
gene_annotation <- readRDS(packer_embryo_rowData_path)

cds <- new_cell_data_set(
    expression_data = expression_matrix,
    cell_metadata = cell_metadata,
    gene_metadata = gene_annotation
)

## ----warning=FALSE, message=FALSE---------------------------------------------
# preprocess scRNA
cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, alignment_group = "batch",
    residual_model_formula_str = "~ bg.300.loading +
        bg.400.loading + bg.500.1.loading + bg.500.2.loading +
        bg.r17.loading + bg.b01.loading + bg.b02.loading")
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

## ----warning=FALSE, message=FALSE, fig.width=7, fig.align="center"------------
p1 <- plot_cells(cds,
    label_groups_by_cluster=FALSE,
    color_cells_by = "embryo.time.bin",
    group_label_size = 5,
    label_cell_groups = FALSE,
    label_leaves = TRUE,
    label_branch_points = TRUE,
    graph_label_size=1.5)
print(p1)

## ----fig.width=7, fig.align="center"------------------------------------------
# order cells - select youngest time point according to embryo.time
# - center of the plot
cds <- order_cells(cds,
    root_cells = colnames(cds)[which(
    colData(cds)$embryo.time == min(colData(cds)$embryo.time))])
colData(cds)["pseudotime"] <- pseudotime(cds)

# plot by pseudotime
p2 <- plot_cells(cds,
    label_groups_by_cluster=FALSE,
    color_cells_by = "pseudotime",
    group_label_size = 5,
    label_cell_groups = FALSE,
    label_leaves = TRUE,
    label_branch_points = TRUE,
    graph_label_size=1.5)
print(p2)

## -----------------------------------------------------------------------------
# load Caenorhabditis elegans gene age estimation
celegans_ps <- readr::read_tsv(file = system.file("extdata",
    "Sun2021_Orthomap.tsv", package = "scTEI"))
table(celegans_ps$Phylostratum)

# define Phylostratum
ps_vec <- setNames(as.numeric(celegans_ps$Phylostratum),
    celegans_ps$GeneID)

## -----------------------------------------------------------------------------
# add TEI values
colData(cds)["TEI"] <- TEI(
    ExpressionSet = counts(cds),
    Phylostratum = ps_vec
)

## ----fig.width=7, fig.align="center"------------------------------------------
# make FeaturePlot
p3 <- plot_cells(cds,
    label_groups_by_cluster=FALSE,
    color_cells_by = "TEI",
    group_label_size = 5,
    label_cell_groups = FALSE,
    label_leaves = TRUE,
    label_branch_points = TRUE,
    graph_label_size=1.5,
    cell_size = 1)
print(p3)

## ----fig.width=7, fig.align="center"------------------------------------------
# make Boxplot
p4 <- ggplot2::ggplot(data.frame(colData(cds)),
    aes(x=embryo.time.bin, y=TEI, fill=embryo.time.bin)) +
    geom_violin() +
    geom_boxplot(width=0.1)
print(p4)

## ----fig.width=7, fig.align="center"------------------------------------------
# make scatter plot - TEI vs pseudotime
p5 <- ggplot(data.frame(colData(cds)),
    aes(x=pseudotime, y=TEI, col=pseudotime)) +
    geom_point()
print(p5)
p6 <- ggplot(data.frame(colData(cds)),
    aes(x=pseudotime, y=TEI, col=embryo.time.bin)) +
    geom_point()
print(p6)

## ----fig.width=7, fig.height=7, fig.align="center"----------------------------
# make grid plot
print(plot_grid(p1, p2, p3, p4))

## ----sessionInfo, echo=TRUE---------------------------------------------------
sessionInfo()

