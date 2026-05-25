# WSI_Download

Utilities for downloading public TCGA-COAD diagnostic HE whole-slide images from the GDC and pairing them with open masked somatic mutation MAF files.

This repository stores code and lightweight documentation only. Downloaded SVS and MAF files should remain local and should not be committed.

## Quick start

```bash
python3 scripts/download_tcga_coad_paired.py --n 457 --file-retries 100 --out tcga_coad_paired_he_mutation
```

Outputs:

- `slides_svs/`: downloaded SVS diagnostic HE slides
- `maf/`: paired masked somatic mutation MAF files
- `metadata/cohort_manifest.tsv`: case, slide, MAF, size, checksum, and path mapping
- `metadata/mutation_labels.tsv`: sample-level mutation flags for common CRC genes

Source: NCI Genomic Data Commons, project `TCGA-COAD`.

## Manual download commands

Run a small smoke test first:

```bash
python3 -u scripts/download_tcga_coad_paired.py --n 5 --file-retries 100 --out tcga_coad_paired_he_mutation
```

Run all currently paired TCGA-COAD cases interactively:

```bash
python3 -u scripts/download_tcga_coad_paired.py --n 457 --file-retries 100 --out tcga_coad_paired_he_mutation
```

Run all currently paired TCGA-COAD cases in the background on macOS/Linux:

```bash
nohup python3 -u scripts/download_tcga_coad_paired.py --n 457 --file-retries 100 --out tcga_coad_paired_he_mutation > tcga_coad_download_all.log 2>&1 &
tail -f tcga_coad_download_all.log
```

If `screen` is available, this is easier to monitor and reconnect:

```bash
screen -dmS tcga_coad_all bash -lc 'python3 -u scripts/download_tcga_coad_paired.py --n 457 --file-retries 100 --out tcga_coad_paired_he_mutation >> tcga_coad_download_all_screen.log 2>&1'
screen -ls
tail -f tcga_coad_download_all_screen.log
```

The downloader is restart-safe for interrupted files: it deletes any stale `.partial` file and any final output whose size does not match GDC metadata before downloading that file again. Use `--file-retries` to control how many times each individual GDC file is retried before the batch exits.
