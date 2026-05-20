#!/usr/bin/env python3
import argparse
import csv
import gzip
import json
import os
import sys
import time
import urllib.parse
import urllib.request


GDC = "https://api.gdc.cancer.gov"
GENES = ["APC", "TP53", "KRAS", "NRAS", "BRAF", "PIK3CA", "SMAD4", "FBXW7"]
CHUNK_SIZE = 1024 * 1024


def post_json(endpoint, payload):
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        f"{GDC}/{endpoint}",
        data=data,
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                return json.load(resp)
        except Exception:
            if attempt == 4:
                raise
            time.sleep(2 * (attempt + 1))


def query_files(filters, fields, page_size=250):
    hits = []
    offset = 0
    while True:
        payload = {
            "filters": filters,
            "fields": ",".join(fields),
            "format": "JSON",
            "size": str(page_size),
            "from": str(offset),
        }
        data = post_json("files", payload)["data"]
        hits.extend(data["hits"])
        total = data["pagination"]["total"]
        if len(hits) >= total or not data["hits"]:
            return hits
        offset += page_size


def filt(*items):
    return {"op": "and", "content": list(items)}


def eq(field, values):
    return {"op": "=", "content": {"field": field, "value": values}}


def first_case_id(hit):
    return hit["cases"][0]["submitter_id"]


def sample_id(hit):
    return hit["cases"][0]["samples"][0]["submitter_id"]


def download_file(file_id, file_name, out_dir, expected_size):
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, file_name)
    tmp_path = out_path + ".partial"
    if os.path.exists(out_path) and os.path.getsize(out_path) == expected_size:
        print(f"exists: {out_path}")
        return out_path

    if os.path.exists(out_path):
        actual_size = os.path.getsize(out_path)
        print(f"removing incomplete file: {out_path} ({actual_size} != {expected_size} bytes)")
        os.remove(out_path)
    if os.path.exists(tmp_path):
        print(f"removing interrupted partial download: {tmp_path}")
        os.remove(tmp_path)

    url = f"{GDC}/data/{urllib.parse.quote(file_id)}"
    for attempt in range(1, 6):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "codex-gdc-downloader"})
            with urllib.request.urlopen(req, timeout=120) as resp, open(tmp_path, "wb") as out:
                total = 0
                last_report = time.time()
                while True:
                    chunk = resp.read(CHUNK_SIZE)
                    if not chunk:
                        break
                    out.write(chunk)
                    total += len(chunk)
                    now = time.time()
                    if now - last_report > 10:
                        pct = 100 * total / expected_size if expected_size else 0
                        print(f"downloading {file_name}: {total / 1e6:.1f} MB ({pct:.1f}%)")
                        last_report = now
            actual_size = os.path.getsize(tmp_path)
            if actual_size != expected_size:
                raise IOError(f"incomplete download for {file_name}: {actual_size} != {expected_size} bytes")
            os.replace(tmp_path, out_path)
            return out_path
        except (OSError, IOError, ConnectionError) as exc:
            if os.path.exists(tmp_path):
                try:
                    os.remove(tmp_path)
                    print(f"deleted interrupted partial download: {tmp_path}")
                except OSError as cleanup_exc:
                    print(f"warning: could not delete partial file {tmp_path}: {cleanup_exc}", file=sys.stderr)
            if attempt == 5:
                raise
            wait_seconds = 10 * attempt
            print(f"download interrupted for {file_name}: {exc}; retrying in {wait_seconds}s ({attempt}/5)")
            time.sleep(wait_seconds)
    return out_path


def parse_maf(maf_path):
    genes = {gene: 0 for gene in GENES}
    counts = {gene: 0 for gene in GENES}
    opener = gzip.open if maf_path.endswith(".gz") else open
    with opener(maf_path, "rt", encoding="utf-8", errors="replace") as handle:
        reader = None
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            reader = csv.DictReader([line] + list(handle), delimiter="\t")
            break
        if reader is None:
            return genes, counts, 0
        total = 0
        for row in reader:
            total += 1
            gene = row.get("Hugo_Symbol", "")
            if gene in genes:
                genes[gene] = 1
                counts[gene] += 1
    return genes, counts, total


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=3, help="number of paired cases to download")
    parser.add_argument("--out", default="tcga_coad_paired_he_mutation")
    args = parser.parse_args()

    slide_filters = filt(
        eq("cases.project.project_id", ["TCGA-COAD"]),
        eq("files.data_type", ["Slide Image"]),
        eq("files.data_format", ["SVS"]),
        eq("files.access", ["open"]),
        eq("files.experimental_strategy", ["Diagnostic Slide"]),
        eq("cases.samples.sample_type", ["Primary Tumor"]),
    )
    maf_filters = filt(
        eq("cases.project.project_id", ["TCGA-COAD"]),
        eq("files.data_type", ["Masked Somatic Mutation"]),
        eq("files.data_format", ["MAF"]),
        eq("files.access", ["open"]),
        eq("files.analysis.workflow_type", ["Aliquot Ensemble Somatic Variant Merging and Masking"]),
    )

    slide_fields = [
        "file_id",
        "file_name",
        "file_size",
        "md5sum",
        "experimental_strategy",
        "cases.submitter_id",
        "cases.samples.submitter_id",
        "cases.samples.sample_type",
    ]
    maf_fields = [
        "file_id",
        "file_name",
        "file_size",
        "md5sum",
        "analysis.workflow_type",
        "cases.submitter_id",
    ]

    print("Querying GDC slides...")
    slides = query_files(slide_filters, slide_fields)
    print(f"Found {len(slides)} primary tumor diagnostic SVS files.")

    print("Querying GDC MAF files...")
    mafs = query_files(maf_filters, maf_fields)
    print(f"Found {len(mafs)} open masked somatic MAF files.")

    maf_by_case = {first_case_id(hit): hit for hit in mafs}
    best_slide_by_case = {}
    for hit in slides:
        case = first_case_id(hit)
        if case not in maf_by_case:
            continue
        current = best_slide_by_case.get(case)
        if current is None or hit["file_size"] < current["file_size"]:
            best_slide_by_case[case] = hit

    selected = sorted(best_slide_by_case.values(), key=lambda hit: hit["file_size"])[: args.n]
    if len(selected) < args.n:
        print(f"Only found {len(selected)} paired cases.", file=sys.stderr)

    root = os.path.abspath(args.out)
    slide_dir = os.path.join(root, "slides_svs")
    maf_dir = os.path.join(root, "maf")
    meta_dir = os.path.join(root, "metadata")
    os.makedirs(meta_dir, exist_ok=True)

    cohort_rows = []
    label_rows = []
    for slide in selected:
        case = first_case_id(slide)
        maf = maf_by_case[case]
        print(f"Downloading paired case {case}")
        slide_path = download_file(slide["file_id"], slide["file_name"], slide_dir, slide["file_size"])
        maf_path = download_file(maf["file_id"], maf["file_name"], maf_dir, maf["file_size"])
        flags, counts, total_variants = parse_maf(maf_path)
        cohort_rows.append(
            {
                "case_submitter_id": case,
                "sample_submitter_id": sample_id(slide),
                "slide_file_id": slide["file_id"],
                "slide_file_name": slide["file_name"],
                "slide_file_size": slide["file_size"],
                "slide_md5sum": slide.get("md5sum", ""),
                "slide_path": slide_path,
                "maf_file_id": maf["file_id"],
                "maf_file_name": maf["file_name"],
                "maf_file_size": maf["file_size"],
                "maf_md5sum": maf.get("md5sum", ""),
                "maf_path": maf_path,
                "maf_total_variants": total_variants,
            }
        )
        row = {"case_submitter_id": case, "maf_total_variants": total_variants}
        for gene in GENES:
            row[f"{gene}_mut"] = flags[gene]
            row[f"{gene}_variant_count"] = counts[gene]
        label_rows.append(row)

    cohort_path = os.path.join(meta_dir, "cohort_manifest.tsv")
    label_path = os.path.join(meta_dir, "mutation_labels.tsv")
    with open(cohort_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(cohort_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(cohort_rows)
    with open(label_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(label_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(label_rows)

    summary_path = os.path.join(root, "README.md")
    with open(summary_path, "w", encoding="utf-8") as handle:
        handle.write("# TCGA-COAD paired HE WSI and mutation data\n\n")
        handle.write("Source: GDC API, project TCGA-COAD. Images are public SVS diagnostic slides; mutations are open masked somatic MAF files.\n\n")
        handle.write(f"Downloaded paired cases: {len(cohort_rows)}\n\n")
        handle.write("Metadata files:\n")
        handle.write("- metadata/cohort_manifest.tsv\n")
        handle.write("- metadata/mutation_labels.tsv\n")

    print(f"Wrote {cohort_path}")
    print(f"Wrote {label_path}")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
