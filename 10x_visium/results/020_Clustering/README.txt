AB1-AB4 = Non-diabetic
AB5-AB8 = Diabetic

.rds files are Seurat objects that contain SCT_snn clusters run at different resolutions:
 - 0.6
 - 0.8
 - 1.0

The resolutions performed can be checked by typing:
colnames(seurat_obj@meta.data)

To perform an action on specifically 1 resolution, set the resolution as follows:

Idents(ab1) <- "RNA_snn_res.0.8"  # Change to resolution 0.8
DimPlot(ab1, reduction = "umap", label = TRUE)
