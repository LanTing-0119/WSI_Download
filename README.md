# WSI_Download

Utilities for downloading public TCGA-COAD diagnostic HE whole-slide images from the GDC and pairing them with open masked somatic mutation MAF files.

This repository stores code and lightweight documentation only. Downloaded SVS and MAF files should remain local and should not be committed.

## Quick start

```bash
python3 scripts/download_tcga_coad_paired.py --n 100 --out tcga_coad_paired_he_mutation
```

Outputs:

- `slides_svs/`: downloaded SVS diagnostic HE slides
- `maf/`: paired masked somatic mutation MAF files
- `metadata/cohort_manifest.tsv`: case, slide, MAF, size, checksum, and path mapping
- `metadata/mutation_labels.tsv`: sample-level mutation flags for common CRC genes

Source: NCI Genomic Data Commons, project `TCGA-COAD`.
