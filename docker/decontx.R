suppressMessages(suppressWarnings(library("singleCellTK")))
suppressMessages(suppressWarnings(library("Matrix")))
suppressMessages(suppressWarnings(library("optparse")))

# args from command line:
args <- commandArgs(TRUE)

option_list <- list(
    make_option(
        c('-f', '--input_file'),
        help='PATH to count matrix input'
    )
)

opt <- parse_args(OptionParser(option_list=option_list))


# change the working directory to co-locate with the counts file:
working_dir <- dirname(opt$input_file)
setwd(working_dir)

# Import counts as a data.frame
# expects the data frame write output as:
# write.table(
#     counts(sce),
#     file = filename,
#     col.names = TRUE,
#     row.names = TRUE,
#     sep = "\t",
#     quote = FALSE
# )
cnts <- read.table(
    file = opt$input_file,
    sep = "\t",
    row.names = 1,
    header=T,
    check.names = FALSE
)

# above, we set check.names=F to prevent the mangling of the sample names.
# Now, we stash those original sample names and run make.names, so that any 
# downstream functions, etc. don't run into trouble. In the end, we convert 
# back to the original names.
orig_cols = colnames(cnts)
new_colnames = make.names(orig_cols)
colnames(cnts) = new_colnames

colname_mapping = data.frame(
    orig_names = orig_cols,
    adjusted_names=new_colnames,
    stringsAsFactors=F
)

# Convert full matrix to sparse so SCTK can work with matrix
cnts <- as(as.matrix(cnts), "sparseMatrix")

# Create an SCE object from the counts
sce <- SingleCellExperiment(
    assays=list(counts=cnts)
)

# Run DecontX
sce.flt <- runDecontX(
    sce,
    useAssay = "counts",
    estimateDelta = T,
    maxIter = 500,
    convergence = 0.001,
    varGenes = 5000,
    dbscanEps = 1
)

# Create a data frame from the DecontX factors and the SCE object
# cell barcodes
df.decontx <- data.frame(
    cell_barcode = as.vector(colnames(sce.flt)),
    decontx_contamination = as.vector(sce.flt$decontX_contamination),
    decontx_class = as.vector(sce.flt$decontX_clusters)
)

# map back to the original colnames
df.final = merge(df.decontx, colname_mapping, by.x='cell_barcode', by.y='adjusted_names')
rownames(df.final) = df.final[, 'orig_names']
df.final = df.final[, c('decontx_contamination', 'decontx_class')]

output_filename <- paste(working_dir, 'sctk_decontx_output.tsv', sep='/')
write.table(
    df.final,
    output_filename,
    sep = "\t",
    quote = F
)

# for WebMEV compatability, need to create an outputs.json file.
json_str = paste0(
       '{"decontaminate_output":"', output_filename, '"}'
)
output_json <- paste(working_dir, 'outputs.json', sep='/')
write(json_str, output_json)