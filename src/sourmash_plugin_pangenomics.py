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
from difflib import get_close_matches

import sourmash
import sourmash_utils
from sourmash import sourmash_args
from sourmash.tax import tax_utils
from sourmash.logging import debug_literal
from sourmash.plugins import CommandLinePlugin
from sourmash.save_load import SaveSignaturesToLocation

DEFAULT_THRESHOLDS = dict(
    CENTRAL_CORE = 0.95,
    EXTERNAL_CORE = 0.90,
    SHELL = 0.10,
    INNER_CLOUD = 0.01,
    SURFACE_CLOUD = 0.00,
)

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
            help="Enable abundance tracking of hashes across rank selection. I.e. Bastardize the abundance metric to count frequency across genomes instead of abundance of kmer across genomes.",
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
            help="The sourmash dictionary created from 'pangenome_creatdb --abund'. Must be a single signature.",
        )
        p.add_argument(
            "-t",
            "--taxonomy-file",
            help="The taxonomy file that was used to create the sourmash pangenome database the signature is from..."
        )
        p.add_argument(
            "-r",
            "--rank",
            help="The rank of the lineage being used as the header of the taxonomy file. I.e. species, genus, family so on"
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
            required=True,
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
        p.add_argument("--thresholds",
                       help="colon-separated thresholds for central core, external core, shell, inner cloud, surface cloud, e.g. 95:90:10:01:00 (which is the default")
        sourmash_utils.add_standard_minhash_args(p)

    def main(self, args):
        super().main(args)

        return classify_hashes_main(args)


#
# pangenome_createdb
#

def pangenome_createdb_main(args):
    print(f"Loading taxonomies from {args.taxonomy_file}")
    taxdb = tax_utils.MultiLineageDB.load(args.taxonomy_file)
    print(f"Found {len(taxdb)} identifiers in taxdb.")

    ident_d = {}
    revtax_d = {}
    do_abund = args.abund
    do_csv = args.csv
    target_rank = args.rank

    counts = defaultdict(Counter) if args.abund else None

    if args.csv:
        csv_file = check_csv(args.csv)
        chunk = []

    select_mh = sourmash_utils.create_minhash_from_args(args)
    print(f"Selecting sketches: {select_mh}")

    for filename in args.sketches:
        print(f"Loading sketches from file {filename}")
        db = sourmash_utils.load_index_and_select(filename, select_mh)

        for n, ss in enumerate(db.signatures()):
            if n > 0 and n % 1000 == 0:
                print(f"...{n} - loading")

            sig_name = ss.name
            ss_mh = ss.minhash
            ident = tax_utils.get_ident(sig_name)

            lineage_tup = taxdb.get(ident)

            if lineage_tup is None:
                if "." in ident:
                    lineage_tup = taxdb.get(ident.split(".")[0])
                else:
                    for i in range(1, 10):
                        lineage_tup = taxdb.get(f"{ident}.{i}")
                        if lineage_tup is not None:
                            break

            if lineage_tup is None:
                print(f"Cannot find ident '{ident}' in the provided taxonomy file.")
                print(f"The three closest matches to '{ident}' are:")
                for k in get_close_matches(ident, taxdb.keys(), n=3):
                    print(f"* '{k}'")
                sys.exit(-1)

            lineage_tup = tax_utils.RankLineageInfo(lineage=lineage_tup)
            lineage_pair = lineage_tup.lineage_at_rank(target_rank)
            lineage_name = lineage_pair[-1].name

            ident_d[lineage_name] = ident

            if lineage_name not in revtax_d:
                revtax_d[lineage_name] = ss.minhash.to_mutable()
            else:
                revtax_d[lineage_name] += ss.minhash

            if do_abund:
                counts[lineage_name].update(set(ss.minhash.hashes))

            if do_csv:
                hash_count = len(revtax_d[lineage_name].hashes)
                chunk.append({
                    "lineage": lineage_name,
                    "sig_name": name,
                    "hash_count": hash_count,
                    "genome_count": n,
                })

                # Write chunk if it reaches capacity
                if len(chunk) >= 1000:
                    write_chunk(chunk, csv_file)
                    chunk = []

    if do_csv and chunk:
        write_chunk(chunk, csv_file)

    print(f"Writing output sketches to '{args.output}'")
    with sourmash_args.SaveSignaturesToLocation(args.output) as save_sigs:
        for n, (lineage_name, ident) in enumerate(ident_d.items()):
            if n > 0 and n % 1000 == 0:
                print(f"...{n} - saving")

            sig_name = f"{ident} {lineage_name}"
            mh = revtax_d[lineage_name]

            if do_abund:
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

    # save!
    print(f"Writing output sketches to '{args.output}'")

    sig_name = "merged" # @CTB update with --name?

    abund_mh = select_mh.copy_and_clear() # hmm, don't need mh tracked above?
    abund_mh.track_abundance = True
    abund_mh.set_abundances(dict(c))

    assert not os.path.exists(args.output) # @CTB
    with sourmash_args.SaveSignaturesToLocation(args.output) as save_sigs:
        print(f"saving to '{args.output}'")

        ss = sourmash.SourmashSignature(abund_mh, name=sig_name)
        save_sigs.add(ss)


def load_all_sketches(
        filename,
        *,
        select_mh=None,
):
    """
    Process a database of sketches into a single merged pangenome sketch.
    """
    print(f"selecting sketches: {select_mh}")

    ss_dict = {}
    print(f"loading sketches from file '{filename}'")
    db = sourmash_utils.load_index_and_select(filename, select_mh)
    print(f"'{filename}' contains {len(db)} signatures")

    assert len(db) == 1, len(db)

    # load the entire database
    for n, ss in enumerate(db.signatures()):
        if n % 10 == 0:
            print(f"...Processing {n} of {len(db)}", end="\r", flush=True)

            name = ss.name

            mh = ss.minhash
            hashes = mh.hashes
            ss_dict[name] = hashes

        print(f"...Processed {n+1} of {len(db)} \n")

    return ss_dict


def load_sketches_by_lineage(filename,
                             lineage_name,
                             *,
                             ignore_case=True,
                             select_mh=None,
):
    """
    Process a database of sketches into a set of merged pangenome sketches,
    based on matches to lineage.
    """
    print(f"selecting sketches: {select_mh}")

    ss_dict = {}
    print(f"loading sketches from file '{filename}'")
    db = sourmash_utils.load_index_and_select(filename, select_mh)
    print(f"'{filename}' contains {len(db)} signatures")

    print(f"Looking for {lineage_name} signature\n")

    # build the search pattern
    if ignore_case:
        pattern = re.compile(lineage_name, re.IGNORECASE)
    else:
        pattern = re.compile(lineage_name)

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
        raise

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

    return ss_dict

def lineage_count(taxonomy, name, rank="species"):
    count = 0
    with open(taxonomy, newline="") as fp:
        reader = csv.DictReader(fp)
        for row in reader:
            if row[rank] == name:
                count += 1

    return count

def calc_pangenome_element_frequency(data, taxonomy, rank):

    # get the pangenome elements of the dicts for each rank pangenome
    for name, hash_dict in data.items():
        lineage_name = " ".join(name.split(" ")[1:])
        if lineage_name.endswith(" singlehash"):
            lineage_name = lineage_name.removesuffix(' singlehash') #https://stackoverflow.com/questions/3663450/remove-substring-only-at-the-end-of-string#comment110394307_61432766
        print(lineage_name)
        tax_max_value = lineage_count(taxonomy, lineage_name, rank="species")
        # get max abundance in genome
        max_value = max(hash_dict.values())
        # number of genomes with this lineage
        print(max_value, tax_max_value)

        if max_value <= tax_max_value:
            max_value = tax_max_value
        # return all hashvals / hash_abunds, along with associated max value
        items = hash_dict.items()
        # sort by abund, highest first
        items = sorted(items, key=lambda x: -x[1])
        for hashval, hash_abund in items:
            freq = hash_abund / max_value
            freq = round(freq, 4)
            yield hashval, freq, hash_abund, max_value


def classify_pangenome_element(freq, *, thresholds=DEFAULT_THRESHOLDS):
    # validate thresholds
    min_threshold = 1.
    for k, v in thresholds.items():
        assert v == float(v)    # must be number
        min_threshold = min(min_threshold, v)
    assert min_threshold == 0   # must catch all hashes :)

    if freq >= thresholds['CENTRAL_CORE']:
        return CENTRAL_CORE
    if freq >= thresholds['EXTERNAL_CORE']:
        return EXTERNAL_CORE
    if freq >= thresholds['SHELL']:
        return SHELL
    if freq >= thresholds['INNER_CLOUD']:
        return INNER_CLOUD
    if freq >= thresholds['SURFACE_CLOUD']:
        return SURFACE_CLOUD

    raise Exception("a hash slipped through the cracks")

#
# pangenome_ranktable
#

# @CTB - rename to count_hashes or something?
def pangenome_ranktable_main(args):
    select_mh = sourmash_utils.create_minhash_from_args(args)

    # load a pre-existing merged/etc database, calc hash info.
    if args.lineage:
        ss_dict = load_sketches_by_lineage(args.data,
                                           args.lineage,
                                           ignore_case=args.ignore_case,
                                           select_mh=select_mh)
    else:
        ss_dict = load_all_sketches(args.data, select_mh=select_mh)

    frequencies = calc_pangenome_element_frequency(ss_dict, args.taxonomy_file, args.rank)

    if args.output_hash_classification:
        print(
            f"Writing hash classification to CSV file '{args.output_hash_classification}'"
        )
        with open(args.output_hash_classification, "w", newline="") as fp:
            w = csv.writer(fp)
            w.writerow(["hashval", "freq", "abund", "max_abund"])

            for hashval, freq, hash_abund, max_value in frequencies:
                w.writerow([hashval, freq, hash_abund, max_value])


#
# pangenome_classify
#

def classify_hashes_main(args):
    # @CTB print out the thresholds or something
    thresholds_str = args.thresholds
    thresholds = dict(DEFAULT_THRESHOLDS)
    if thresholds_str:
        threshold_vals = thresholds_str.split(':')
        assert len(threshold_vals) == 5
        threshold_vals = [ int(t)/100 for t in threshold_vals ]
        assert max(threshold_vals) <= 1
        assert min(threshold_vals) >= 0

        thresholds['CENTRAL_CORE'] = threshold_vals[0]
        thresholds['EXTERNAL_CORE'] = threshold_vals[1]
        thresholds['SHELL'] = threshold_vals[2]
        thresholds['INNER_CLOUD'] = threshold_vals[3]
        thresholds['SURFACE_CLOUD'] = threshold_vals[4]
    
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

    # load in all the frequencies etc, and classfy
    for csv_file in args.ranktable_csv_files:
        with open(csv_file, "r", newline="") as fp:
            r = csv.DictReader(fp)

            classify_d = {}
            for row in r:
                hashval = int(row["hashval"])
                abund = int(row["abund"])
                max_abund = int(row["max_abund"])
                freq = abund / max_abund
                classify_as = classify_pangenome_element(freq,
                                                         thresholds=thresholds)

                assert hashval not in classify_d, "hashval already encountered"
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
