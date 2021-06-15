## Switch to MSI

The function "multiSetMEs" was not working on Sapelo2 (even on example data). I uploaded the following files to MSI:
 - consensusTom_signed.RData 
 - Consensus-dataInput.RData
 - multiExpr.RData

#### checksums
````bash
md5sum ConsensusTOM_signed.RData # bbdeaeb153008c79b83c10a9e909e44f  ConsensusTOM_signed.RData
md5sum Consensus-dataInput.RData # b796038fa0d5f93e63568e3f050079b5  Consensus-dataInput.RData
md5sum multiExpr.RData # a3548031c642345d4d93137cbc6f3272  multiExpr.RData
````

#### copied to google drive
````bash
rclone copy -v /scratch/eld72413/Salty_Nut/CultivatedOnly/DE_Analyses_Inbred/Consensus_Network/ConsensusTOM_signed.RData google_drive_emily:Sunflower_data
rclone copy -v /lustre2/scratch/eld72413/Salty_Nut/CultivatedOnly/DE_Analyses_Inbred/Consensus-dataInput.RData google_drive_emily:Sunflower_data
rclone copy -v /lustre2/scratch/eld72413/Salty_Nut/CultivatedOnly/DE_Analyses_Inbred/multiExpr.RData google_drive_emily:Sunflower_data
````
- need to be in xfer node to use rclone

#### Upload to MSI
Logged into MSI

```bash
cd /scratch.global/edittmar/WGCNA
module load rclone/1.54

rclone copy -v google_drive_emily:Sunflower_data/ConsensusTOM_signed.RData .
md5sum ConsensusTOM_signed.RData # bbdeaeb153008c79b83c10a9e909e44f  ConsensusTOM_signed.RData

rclone copy -v google_drive_emily:Sunflower_data/Consensus-dataInput.RData .
md5sum Consensus-dataInput.RData # b796038fa0d5f93e63568e3f050079b5  Consensus-dataInput.RData

rclone copy -v google_drive_emily:Sunflower_data/multiExpr.RData .
md5sum multiExpr.RData # a3548031c642345d4d93137cbc6f3272  multiExpr.RData
```

#### Install WGCNA

```bash
module load R/3.6.0
R
```

```R
library(BiocManager)

r <- getOption("repos");
r["CRAN"] <- "http://cran.rstudio.com/";
options(repos=r);
#install.packages("packagename");
BiocManager::install("WGCNA") # followed 'yes' at prompts to install a personal library
```
