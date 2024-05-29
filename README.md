# sourmash_plugin_pangenomics: tools for sourmash-based pangenome analyses

## Installation

```
pip install sourmash_plugin_pangenomics
```

## Quickstart

You can run all of these commands in the `test_workflow` directory of the git repository.

### Build a pangenome database using lineages

(CTB: explain contents!)

The following command builds a pangenome database for the species present in the lineages file `gtdb-rs214-agatha.lineages.csv.gz ` (currently only `s__Agathobacter faecis`), using the sketches present in the `gtdb-rs214-agatha-k21.zip`.

```
sourmash scripts pangenome_createdb \
    gtdb-rs214-agatha-k21.zip \
    -t gtdb-rs214-agatha.lineages.csv.gz \
    -o agatha-merged.sig.zip --abund -k 21
```

The output file is `agatha-merged.sig.zip` and contains the following:
```
% sourmash sig summarize agatha-merged.sig.zip

...
num signatures: 1
** examining manifest...
total hashes: 27398
summary of sketches:
   1 sketches with DNA, k=21, scaled=1000, abund      27398 total hashes
```

Note: the command `pangenome_merge` (see below) will construct a pangenome
sketch by merging all provided signatures.

### Build a pangenome "ranktable"

A "ranktable" is our name for a database that assigns hashes a pangenomic "rank" - central core, external core, shell, inner cloud, or surface cloud.

The following command builds a ranktable for the species `s__Agathobacter faecis`, selected from the pangenome database created above:
```
sourmash scripts pangenome_ranktable \
    agatha-merged.sig.zip \
    -o test_output/agathobacter_faecis.csv \
    -k 21 -l 'GCF_020557615 s__Agathobacter faecis'
```

The output file is `test_output/agathobacter_faecis.csv`, and it contains two columns:

```
hashval,pangenome_classification
96834755571756,1
119187685848053,1
129679169912030,1
...
18440589591308259,4
18443409651295626,4
18446214016691046,4
```
where the first column is the hash value, and the second column is the pangenome rank for that hash.

### Summarize the ranks of the hashes in a sketch

We can now use our ranktable to summarize _any_ sketch, including a metagenome. Here we use a human gut metagenome, `SRR5650070`:

```
sourmash scripts pangenome_classify \
    SRR5650070.trim.sig.zip \
    test_output/agathobacter_faecis.csv \
    -k 21
```

This will yield the following output:

```
For 'test_output/agathobacter_faecis.csv', signature 'SRR5650070' contains:
         497 (12.5%) hashes are classified as central core
         427 (10.8%) hashes are classified as external core
         1791 (45.2%) hashes are classified as shell
         1251 (31.5%) hashes are classified as inner cloud
         0 (0.0%) hashes are classified as surface cloud
         ...and 262716 hashes are NOT IN the csv file
```

### Build a pangenome sketch without using lineages

(CTB: explain contents!)

The following command builds a pangenome sketch by combining all provided sketches. Here we use the sketches present in the `gtdb-rs214-agatha-k21.zip` file:

```
sourmash scripts pangenome_merge \
    gtdb-rs214-agatha-k21.zip \
    -o agatha-merged-2.sig.zip-k 21
```

The output file is `agatha-merged-2.sig.zip` and is identical
(via e.g. `sourmash compare`) to the `agatha-merged.sig.zip` file.

## Support

We suggest filing issues in [the main sourmash issue tracker](https://github.com/dib-lab/sourmash/issues) as that receives more attention (and is monitored by the same people anyway)!

## Dev docs

`sourmash_plugin_pangenomics` is developed at https://github.com/sourmash-bio/sourmash_plugin_pangenomics.

### Testing

The current tests are implemented as Snakemake workflow in `test_workflow/`. To run them, execute the following command in the main directory:
```
make cleanrun
```

### Generating a release

Bump version number in `pyproject.toml` and push.

Make a new release on github.

Then pull, and:

```
python -m build
```

followed by `twine upload dist/...`.
