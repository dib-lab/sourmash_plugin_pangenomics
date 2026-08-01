#! /usr/bin/env python

import json
import csv
import gzip
import hashlib
import os

outdir = "sample-genomes"
sigdir = os.path.join(outdir, "signatures")
os.makedirs(sigdir, exist_ok=True)

hashes = [100,200,300,400,500,600,700,800,900,1000]

manifest_rows = []
lineage_rows = []

for genome in range(1, 11):

    ident = f"GCA_{genome}.1"
    name = f"Genome {genome}"
    filename = f"sample-{genome}.sig.gz"
    internal_location = f"signatures/{filename}"
    sig_path = os.path.join(sigdir, filename)

    abundances = [(i + 1) * genome for i in range(10)]

    sig = [{
        "class": "sourmash_signature",
        "email": "",
        "hash_function": "0.murmur64",
        "filename": sig_path,
        "name": ident + " " + name,
        "license": "CC0",
        "signatures": [{
            "num": 0,
            "ksize": 31,
            "seed": 42,
            "max_hash": 18446744073709552,
            "mins": hashes,
            "md5sum": hashlib.md5(
                "".join(map(str, hashes)).encode()
            ).hexdigest(),
            "abundances": abundances,
            "molecule": "DNA"
        }],
        "version": 0.4
    }]

    with gzip.open(sig_path, "wt") as fp:
        json.dump(sig, fp)

    md5 = hashlib.md5(
        "".join(map(str, hashes)).encode()
    ).hexdigest()

    manifest_rows.append([
        internal_location,
        md5,
        md5[:8],
        31,
        "DNA",
        0,
        1000,
        len(hashes),
        True,
        name,
        sig_path
    ])

    lineage_rows.append([
        ident,
        "t",
        "d__Test",
        "p__Test",
        "c__Test",
        "o__Test",
        "f__Test",
        "g__Test",
        "s__Test"
    ])

manifest = os.path.join(outdir, "SOURMASH-MANIFEST.csv")

with open(manifest, "w", newline="") as fp:
    fp.write("# SOURMASH-MANIFEST-VERSION: 1.0\n")

    writer = csv.writer(fp)
    writer.writerow([
        "internal_location",
        "md5",
        "md5short",
        "ksize",
        "moltype",
        "num",
        "scaled",
        "n_hashes",
        "with_abundance",
        "name",
        "filename"
    ])

    writer.writerows(manifest_rows)

with open("lineages.sample-genomes.csv", "w") as fp:
    writer = csv.writer(fp)

    writer.writerow([
        "ident",
        "gtdb_representative",
        "superkingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species"
    ])

    writer.writerows(lineage_rows)
