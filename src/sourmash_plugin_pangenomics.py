"""pangenomics plugin"""

epilog = """
See https://github.com/dib-lab/sourmash_plugin_pangenomics for more examples.

Need help? Have questions? Ask at http://github.com/sourmash-bio/sourmash/issues!
"""

import argparse
import sys
from collections import Counter, defaultdict
import csv
import os
import re
import pprint

import sourmash
import sourmash_utils
from sourmash import sourmash_args
from sourmash.tax import tax_utils
from sourmash.logging import debug_literal
from sourmash.plugins import CommandLinePlugin
from sourmash.save_load import SaveSignaturesToLocation


CENTRAL_CORE = 1
EXTERNAL_CORE = 2
SHELL = 3
INNER_CLOUD = 4
SURFACE_CLOUD = 5

NAMES = {
    CENTRAL_CORE: "central core",
    EXTERNAL_CORE: "external core",
    SHELL: "shell",
    INNER_CLOUD: "inner cloud",
    SURFACE_CLOUD: "surface cloud",
}


###

#
# CLI plugins - supports 'sourmash scripts <commands>'
#


class Command_CreateDB(CommandLinePlugin):
    command = "pangenome_createdb"  # 'scripts <command>'
    description = "create a database of sketches merged at given rank"  # output with -h
    usage = "pangenome_createdb <db> -t <lineagedb> -o <merged>.zip"  # output with no args/bad args as well as -h
    epilog = epilog  # output with -h
    formatter_class = argparse.RawTextHelpFormatter  # do not reformat multiline

    def __init__(self, subparser):
        super().__init__(subparser)
        p = subparser

        p.add_argument(
            "-t",
            "--taxonomy-file",
            "--taxonomy",
            metavar="FILE",
            action="extend",
            nargs="+",
            required=True,
            help="database lineages file",
        )
        p.add_argument("sketches", nargs="+", help="sketches to combine")
        p.add_argument(
            "-o",
            "--output",
            required=True,
            help="Define a filename for the pangenome signatures (.zip preferred).",
        )
        p.add_argument(
            "--csv",
            help="A CSV file generated to contain the lineage rank, genome name, hash count, and genome count.",
        )
        p.add_argument("-r", "--rank", default="species")
        p.add_argument(
            "-a",
            "--abund",
            action="store_true",
            help="Enable abundance tracking of hashes across rank selection.",
        )
        sourmash_utils.add_standard_minhash_args(p)

    def main(self, args):
        super().main(args)
        return pangenome_createdb_main(args)


class Command_Merge(CommandLinePlugin):
    command = "pangenome_merge"  # 'scripts <command>'
    description = "merge sketches into a pangenome database a la createdb"  # output with -h
    usage = "pangenome_merge <sketches> -o <merged>.zip"  # output with no args/bad args as well as -h
    epilog = epilog  # output with -h
    formatter_class = argparse.RawTextHelpFormatter  # do not reformat multiline

    def __init__(self, subparser):
        super().__init__(subparser)
        p = subparser

        p.add_argument("sketches", nargs="+", help="sketches to combine")
        p.add_argument(
            "-o",
            "--output",
            required=True,
            help="Define a filename for the pangenome signatures (.zip preferred).",
        )
        sourmash_utils.add_standard_minhash_args(p)

    def main(self, args):
        super().main(args)
        return pangenome_merge_main(args)


class Command_RankTable(CommandLinePlugin):
    command = "pangenome_ranktable"  # 'scripts <command>'
    description = "create a CSV ranktable that annotates hashes with pangenome characters"  # output with -h
    usage = "pangenome_ranktable <merged.zip> -o <ranktable>.csv -l <lineage>"  # output with no args/bad args as well as -h
    epilog = epilog  # output with -h
    formatter_class = argparse.RawTextHelpFormatter  # do not reformat multiline

    def __init__(self, subparser):
        super().__init__(subparser)
        p = subparser
        p.add_argument(
            "data",
            metavar="SOURMASH_DATABASE",
            help="The sourmash dictionary created from 'pangenome_creatdb --abund'",
        )
        p.add_argument(
            "-l",
            "--lineage",
            help='The specific lineage to extract from the sourmash pangenome database (e.g. "s__Escherichia coli")',
        )
        p.add_argument(
            "-i",
            "--ignore-case",
            action="store_true",
            help="Ignore the casing of search terms",
        )
        p.add_argument(
            "-o",
            "--output-hash-classification",
            required=False,
            help="CSV file containing classification of each hash",
        )
        sourmash_utils.add_standard_minhash_args(p)

    def main(self, args):
        super().main(args)
        return pangenome_ranktable_main(args)


class Command_Classify(CommandLinePlugin):
    command = "pangenome_classify"  # 'scripts <command>'
    description = "classify the hashes in a sketch based on given ranktable characters"  # output with -h
    usage = "pangenome_classify <sketch> <ranktable1> [<ranktable2> ...]"  # output with no args/bad args as well as -h
    epilog = epilog  # output with -h
    formatter_class = argparse.RawTextHelpFormatter  # do not reformat multiline

    def __init__(self, subparser):
        super().__init__(subparser)
        p = subparser
        p.add_argument("metagenome_sig")
        p.add_argument("ranktable_csv_files", nargs="+",
                       help="rank tables produced by pangenome_ranktable")
        sourmash_utils.add_standard_minhash_args(p)

    def main(self, args):
        super().main(args)
        return classify_hashes_main(args)


#
# pangenome_createdb
#

def pangenome_createdb_main(args):
    print(f"loading taxonomies from {args.taxonomy_file}")
    taxdb = sourmash.tax.tax_utils.MultiLineageDB.load(args.taxonomy_file)
    print(f"found {len(taxdb)} identifiers in taxdb.")

    ident_d = {}
    revtax_d = {}
    accum = defaultdict(dict)
    if args.abund:
        counts = {}
    if args.csv:
        csv_file = check_csv(args.csv)

    select_mh = sourmash_utils.create_minhash_from_args(args)
    print(f"selecting sketches: {select_mh}")

    # Load the database
    for filename in args.sketches:
        print(f"loading sketches from file {filename}")
        db = sourmash_utils.load_index_and_select(filename, select_mh)

        if args.csv:
            chunk = []

        # Work on a single signature at a time across the database
        for n, ss in enumerate(db.signatures()):
            if n and n % 1000 == 0:
                print(f"...{n} - loading")

            name = ss.name
            ident = tax_utils.get_ident(name)

            # grab relevant lineage name
            try:
                lineage_tup = taxdb[ident]
            except KeyError:  # older versions of genbank are not named with version!
                try:
                    if "." in ident:
                        short_ident = ident.split(".")[0]
                        lineage_tup = taxdb[short_ident]
                    else:
                        for i in range(1, 10):
                            try:
                                new_ident = f"{ident}.{i}"
                                lineage_tup = taxdb[new_ident]
                                break
                            except KeyError:
                                continue
                except:
                    print("Wow, that sucks!")

            lineage_tup = tax_utils.RankLineageInfo(lineage=lineage_tup)
            lineage_pair = lineage_tup.lineage_at_rank(args.rank)
            lineage_name = lineage_pair[-1].name

            ident_d[lineage_name] = (
                ident  # pick an ident to represent this set of pangenome sketches
            )

            # Accumulate the count within lineage names if `--abund` in cli
            if args.abund:
                c = Counter(ss.minhash.hashes)
                if lineage_name in counts:
                    counts[lineage_name].update(c)
                else:
                    counts[lineage_name] = c

            # track merged sketches
            mh = revtax_d.get(lineage_name)

            if mh is None:
                mh = ss.minhash.to_mutable()
                revtax_d[lineage_name] = mh
            else:
                mh += ss.minhash

            ## Add {name, hash_count} to a lineage key then
            ## create a simpler dict for writing the csv
            if args.csv:
                # Accumulated counts of hashes in lineage by genome
                hash_count = len(mh.hashes)

                accum[lineage_name][name] = (
                    accum[lineage_name].get(name, 0) + hash_count
                )
                chunk.append(
                    {
                        "lineage": lineage_name,
                        "sig_name": name,
                        "hash_count": hash_count,
                        "genome_count": n,
                    }
                )

            if args.csv and len(chunk) >= 1000:  # args.chunk_size?
                write_chunk(chunk, csv_file)  # args.outputfilenameforcsv?
                accum = defaultdict(dict)
                chunk = []

        # Write remaining data
        if args.csv and len(chunk) > 0:
            write_chunk(chunk, csv_file)
            accum = defaultdict(dict)
            chunk = []

    # save!
    print(f"Writing output sketches to '{args.output}'")
    with sourmash_args.SaveSignaturesToLocation(args.output) as save_sigs:
        for n, (lineage_name, ident) in enumerate(ident_d.items()):
            if n and n % 1000 == 0:
                print(f"...{n} - saving")

            sig_name = f"{ident} {lineage_name}"

            # retrieve merged MinHash
            mh = revtax_d[lineage_name]

            # Add abundance to signature if `--abund` in cli
            if args.abund:
                abund_d = dict(counts[lineage_name])

                abund_mh = mh.copy_and_clear()
                abund_mh.track_abundance = True
                abund_mh.set_abundances(abund_d)

                ss = sourmash.SourmashSignature(abund_mh, name=sig_name)
            else:
                ss = sourmash.SourmashSignature(mh, name=sig_name)

            save_sigs.add(ss)


# Chunk function to limit the memory used by the hash_count dict and list
def write_chunk(chunk, output_file):
    with open(output_file, "a", newline="") as csvfile:
        fieldnames = ["lineage", "sig_name", "hash_count", "genome_count"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writerows(chunk)


# @CTB do we need this?
def check_csv(csv_file):
    if csv_file is None:
        return

    if os.path.exists(csv_file):
        raise argparse.ArgumentTypeError("\n%s already exists" % csv_file)

    count_csv = os.path.splitext(csv_file)[0] + ".csv"
    return count_csv


#
# pangenome_merge
#

def pangenome_merge_main(args):
    select_mh = sourmash_utils.create_minhash_from_args(args)
    print(f"selecting sketches: {select_mh}")

    # Load the database
    for filename in args.sketches:
        print(f"loading sketches from file {filename}")
        db = sourmash_utils.load_index_and_select(filename, select_mh)

        c = Counter()
        mh = None

        # work across the entire database
        for n, ss in enumerate(db.signatures()):
            if n and n % 1000 == 0:
                print(f"...{n} - loading")

            c.update(ss.minhash.hashes)

            if mh is None:
                mh = ss.minhash.to_mutable()
            else:
                mh += ss.minhash

    # save!
    print(f"Writing output sketches to '{args.output}'")

    sig_name = "merged" # @CTB update with --name?

    abund_mh = mh.copy_and_clear() # hmm, don't need mh tracked above?
    abund_mh.track_abundance = True
    abund_mh.set_abundances(dict(c))

    with sourmash_args.SaveSignaturesToLocation(args.output) as save_sigs:
        print(f"saving to '{args.output}'")

        ss = sourmash.SourmashSignature(abund_mh, name=sig_name)
        save_sigs.add(ss)


def db_process(
        filename,
        ignore_case,
        select_mh=None,
        lineage_name="None",
):
    print(f"selecting sketches: {select_mh}")

    ss_dict = {}
    print(f"loading sketches from file '{filename}'")
    db = sourmash_utils.load_index_and_select(filename, select_mh)
    print(f"'{filename}' contains {len(db)} signatures")

    if lineage_name:
        print(f"Looking for {lineage_name} signature\n")

        # build the search pattern
        pattern = lineage_name
        if ignore_case:
            pattern = re.compile(pattern, re.IGNORECASE)
        else:
            pattern = re.compile(pattern)

        def search_pattern(vals):
            return any(pattern.search(val) for val in vals)

        # find all matching rows.
        mf = db.manifest
        sub_mf = mf.filter_on_columns(search_pattern, ["name", "filename", "md5"])

        selected_sigs = []
        print(f"Found {len(sub_mf)} signatures in '{filename}':")

        for n, row in enumerate(sub_mf.rows, start=1):
            print(f'{n:<15} \033[0;31m{row.get("name")}\033[0m')
            selected_sigs.append(row.get("name"))

        sub_picklist = sub_mf.to_picklist()

        try:
            db = db.select(picklist=sub_picklist)
        except ValueError:
            error("Chosen lineage name input not supported")

        for ss in db.signatures():
            name = ss.name

            print(f"Found \033[0;31m{name}\033[0m in '{filename}'")

            mh = ss.minhash
            hashes = mh.hashes
            ss_dict[name] = hashes

        total_rows_examined = 0
        total_rows_examined += len(mf)

        print(
            f"\nloaded {total_rows_examined} total that matched ksize & molecule type"
        )

        if ss_dict:
            print(f"extracted {len(ss_dict)} signatures from {len(db)} file(s)\n")

        else:
            print("no matching signatures found!\n")
            sys.exit(-1)

    else:
        # process the entire database

        for n, ss in enumerate(db.signatures()):
            if n % 10 == 0:
                print(f"...Processing {n} of {len(mf)}", end="\r", flush=True)

                name = ss.name

                mh = ss.minhash
                hashes = mh.hashes
                ss_dict[name] = hashes

            print(f"...Processed {n} of {len(mf)} \n")

    return ss_dict


def pangenome_elements(data):
    # get the pangenome elements of the dicts for each rank pangenome
    for i, (key, nested_dict) in enumerate(data.items()):
        max_value = max(nested_dict.values())

        central_core_threshold = 0.95  # 0.95 is core , 90% is technically soft core
        external_core_threshold = 0.90
        shell_threshold = 0.10  # 0.10
        inner_cloud_threshold = 0.01  # 0.0 is the full cloud, but trimming (0.001?) may be necessary to create the viz...?
        surface_cloud_threshold = 0.00

        central_core = []
        external_core = []
        shell = []
        inner_cloud = []
        surface_cloud = []

        for nested_key, nested_value in nested_dict.items():
            if nested_value >= max_value * central_core_threshold:
                central_core.append((nested_key, nested_value))
            elif nested_value >= max_value * external_core_threshold:
                external_core.append((nested_key, nested_value))
            elif nested_value >= max_value * shell_threshold:
                shell.append((nested_key, nested_value))
            elif nested_value >= max_value * inner_cloud_threshold:
                inner_cloud.append((nested_key, nested_value))
            elif nested_value >= max_value * surface_cloud_threshold:
                surface_cloud.append((nested_key, nested_value))
        return central_core, external_core, shell, inner_cloud, surface_cloud

#
# pangenome_ranktable
#

def pangenome_ranktable_main(args):
    select_mh = sourmash_utils.create_minhash_from_args(args)

    ss_dict = db_process(
        filename=args.data,
        select_mh=select_mh,
        lineage_name=args.lineage,
        ignore_case=args.ignore_case,
    )
    results = pangenome_elements(ss_dict)

    if args.output_hash_classification:
        print(
            f"Writing hash classification to CSV file '{args.output_hash_classification}'"
        )
        with open(args.output_hash_classification, "w", newline="") as fp:
            w = csv.writer(fp)
            w.writerow(["hashval", "pangenome_classification"])
            central_core, external_core, shell, inner_cloud, surface_cloud = results

            for xx, classify_code in (
                (central_core, CENTRAL_CORE),
                (external_core, EXTERNAL_CORE),
                (shell, SHELL),
                (inner_cloud, INNER_CLOUD),
                (surface_cloud, SURFACE_CLOUD),
            ):
                for hashval, _ in xx:
                    w.writerow([hashval, classify_code])


#
# pangenome_classify
#

def classify_hashes_main(args):
    select_mh = sourmash_utils.create_minhash_from_args(args)
    print(f"selecting sketches: {select_mh}")

    db = sourmash_utils.load_index_and_select(args.metagenome_sig, select_mh)
    sketches = list(db.signatures())
    assert len(sketches) == 1
    sketch = sketches[0]
    minhash = sketch.minhash
    hashes = minhash.hashes

    central_core_mh = minhash.copy_and_clear()
    shell_mh = minhash.copy_and_clear()

    for csv_file in args.ranktable_csv_files:
        with open(csv_file, "r", newline="") as fp:
            r = csv.DictReader(fp)

            classify_d = {}
            for row in r:
                hashval = int(row["hashval"])
                classify_as = int(row["pangenome_classification"])
                classify_d[hashval] = classify_as

        # now, classify_d has the classifications we care about. Let's
        # do terrible things to the sketch now.

        counter_d = defaultdict(int)
        total_classified = 0
        for hashval in hashes:
            classify = classify_d.get(hashval, -1)
            counter_d[classify] += 1
            if classify >= 0:
                total_classified += 1
            if classify == CENTRAL_CORE:
                central_core_mh.add_hash(hashval)
            elif classify == SHELL:
                shell_mh.add_hash(hashval)

        print(f"For '{csv_file}', signature '{sketch.name}' contains:")
        for int_id in sorted(NAMES):
            name = NAMES[int_id]
            count = counter_d.get(int_id, 0)
            percent = count / total_classified * 100
            print(f"\t {count} ({percent:.1f}%) hashes are classified as {name}")

        count = counter_d.get(-1, 0)
        print(f"\t ...and {count} hashes are NOT IN the csv file")
