# SaltxNut


## Diagnostics

- Called SNPs to identify outliers
	- missing 90% of SNP data: 37, 45, 63, 186, 254, 261, 336
	- Mis-ID'd: 37, 323, 335
	- high pairwise error: 39, 67
- Library sizes
	- smaller than 5 M: 45, 63, 254, 261, 336; 186 is smaller than 7 M

- PCA
	- outliers: 37, 63, 254, 261, 323, 335, 336 (slight outliers: 39, 45, 67)
		- all 5 that have small library sizes (below 5 M.)
		- all 3 that are Mis'ID'd
		- 2 with high pairwise error

- Dendrogram
	- If cut height is 300, outliers: 37, 45, 63, 186, 254, 261, 323, 336
		- other possible outlier: 56 RHA 373 highsalt
		- does not include #39, 67, 335 (though 67 still a slight outlier)

### Summary:
- Removing 3 mislabelled samples, 6 with library sizes smaller than 7 M - 8 of these are outliers on dendrogram
	- 37 mislabelled, dendrogram outlier; also PCA outlier
	- 45 library size smaller than 5 M, dendrogram oulier
	- 63 
