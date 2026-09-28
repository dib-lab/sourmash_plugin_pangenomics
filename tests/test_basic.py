"""
Tests for sourmash_plugin_pangenomics.
"""
import os
import pytest
import collections
from dataclasses import dataclass

import sourmash
import sourmash_tst_utils as utils
from sourmash_tst_utils import SourmashCommandFailed


from sourmash_plugin_pangenomics import pangenome_createdb_main


def get_workflow_data(filename):
    return os.path.join(os.path.dirname(__file__),
                        '../test_workflow',
                        filename)

# basic test, does nothing except check the test harness
# (sourmash should fail to run ;)
def test_run_sourmash(runtmp):
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('', fail_ok=True)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert runtmp.last_result.status != 0                    # no args provided, ok ;)


# make a fake 'args' object to mimic argparse
@dataclass
class FakeCreateDBArgs:
    taxonomy_file: list
    sketches: list
    output: str

    ksize: int = 21
    rank: str = 'species'
    abund: bool = False
    scaled: int = 1000
    dna: bool = True
    protein: bool = False
    dayhoff: bool = False
    hp: bool = False
    skipm1n3: bool = False
    skipm2n3: bool = False
    csv: object = None


# can we get createdb to run?
def test_createdb_1_basic(runtmp):
    taxdb_path = get_workflow_data('gtdb-rs214-agatha.lineages.csv.gz')
    print(taxdb_path)

    db_path = get_workflow_data('gtdb-rs214-agatha-k21.zip')
    print(db_path)

    outfile = runtmp.output('output.sig.zip')

    args = FakeCreateDBArgs(sketches=[db_path], taxonomy_file=[taxdb_path],
                            output=outfile)

    pangenome_createdb_main(args)

    assert os.path.exists(outfile)

    db = sourmash.load_file_as_index(outfile)
    assert len(db) == 1

    siglist = list(db.signatures())
    assert siglist[0].name == 'GCF_020557615 s__Agathobacter faecis'


# test multiple species
def test_createdb_2_multiple_species(runtmp):
    taxdb_path = get_workflow_data('gtdb-rs214-agatha.lineages.csv.gz')
    print(taxdb_path)

    db_path = get_workflow_data('gtdb-rs214-agatha-k21.zip')
    print(db_path)

    # add something from s__Agathobacter sp900317585 per GTDB
    extra_sig = utils.get_test_data('GCA_945878335.1.sig.zip')

    outfile = runtmp.output('output.sig.zip')

    args = FakeCreateDBArgs(sketches=[db_path, extra_sig],
                            taxonomy_file=[taxdb_path],
                            output=outfile)

    pangenome_createdb_main(args)

    assert os.path.exists(outfile)

    db = sourmash.load_file_as_index(outfile)
    assert len(db) == 2

    siglist = list(db.signatures())
    names = [ ss.name for ss in siglist ]
    assert 'GCF_020557615 s__Agathobacter faecis' in names
    assert 'GCA_945878335 s__Agathobacter sp900317585' in names


# switch to using a smaller, simpler data set that we can muck around with
# first test: does merge do what we expect it to do with the minhashes?
def test_createdb_3_xeno(runtmp):
    taxdb_path = utils.get_test_data('xeno.lineages.csv')

    db_path1 = utils.get_test_data('xeno-1.sig')
    db_path2 = utils.get_test_data('xeno-2.sig')

    outfile = runtmp.output('output.sig.zip')

    args = FakeCreateDBArgs(sketches=[db_path1, db_path2],
                            taxonomy_file=[taxdb_path],
                            output=outfile, ksize=31)

    pangenome_createdb_main(args)

    assert os.path.exists(outfile)

    db = sourmash.load_file_as_index(outfile)
    assert len(db) == 1

    siglist = list(db.signatures())
    names = [ ss.name for ss in siglist ]
    assert 'GCF_900116635 s__Xenorhabdus koppenhoeferi' in names

    # check minhash by doing a manual merge
    ss1 = next(iter(sourmash.load_file_as_signatures(db_path1)))
    ss2 = next(iter(sourmash.load_file_as_signatures(db_path2)))
    mh1 = ss1.minhash
    mh2 = ss2.minhash
    mh3 = mh1 + mh2

    species_mh = siglist[0].minhash
    assert species_mh == mh3


# does merge do what we expect it to do with the minhashes, with just one mh?
def test_createdb_4_xeno(runtmp):
    taxdb_path = utils.get_test_data('xeno.lineages.csv')

    db_path1 = utils.get_test_data('xeno-2.sig')

    outfile = runtmp.output('output.sig.zip')

    args = FakeCreateDBArgs(sketches=[db_path1],
                            taxonomy_file=[taxdb_path],
                            output=outfile, ksize=31)

    pangenome_createdb_main(args)

    assert os.path.exists(outfile)

    db = sourmash.load_file_as_index(outfile)
    assert len(db) == 1

    siglist = list(db.signatures())
    names = [ ss.name for ss in siglist ]
    assert 'GCF_900116635 s__Xenorhabdus koppenhoeferi' in names

    # check minhash is same
    ss1 = next(iter(sourmash.load_file_as_signatures(db_path1)))
    mh1 = ss1.minhash

    species_mh = siglist[0].minhash
    assert species_mh == mh1


# mess around with version numbers in lineages file
def test_createdb_5_xeno(runtmp):
    taxdb_path = utils.get_test_data('xeno-version.lineages.csv')

    db_path1 = utils.get_test_data('xeno-1.sig')
    db_path2 = utils.get_test_data('xeno-2.sig')
    outfile = runtmp.output('output.sig.zip')

    args = FakeCreateDBArgs(sketches=[db_path1, db_path2],
                            taxonomy_file=[taxdb_path],
                            output=outfile, ksize=31)

    pangenome_createdb_main(args)

    assert os.path.exists(outfile)

    db = sourmash.load_file_as_index(outfile)
    assert len(db) == 1

    siglist = list(db.signatures())
    names = [ ss.name for ss in siglist ]
    assert 'GCF_900116635 s__Xenorhabdus koppenhoeferi' in names

    # check minhash by doing a manual merge
    ss1 = next(iter(sourmash.load_file_as_signatures(db_path1)))
    ss2 = next(iter(sourmash.load_file_as_signatures(db_path2)))
    mh1 = ss1.minhash
    mh2 = ss2.minhash
    mh3 = mh1 + mh2

    species_mh = siglist[0].minhash
    assert species_mh == mh3


# remove version numbers completely
def test_createdb_5_xeno(runtmp):
    taxdb_path = utils.get_test_data('xeno-noversion.lineages.csv')

    db_path1 = utils.get_test_data('xeno-1.sig')
    db_path2 = utils.get_test_data('xeno-2.sig')
    outfile = runtmp.output('output.sig.zip')

    args = FakeCreateDBArgs(sketches=[db_path1, db_path2],
                            taxonomy_file=[taxdb_path],
                            output=outfile, ksize=31)

    pangenome_createdb_main(args)

    assert os.path.exists(outfile)

    db = sourmash.load_file_as_index(outfile)
    assert len(db) == 1

    siglist = list(db.signatures())
    names = [ ss.name for ss in siglist ]
    assert 'GCF_900116635 s__Xenorhabdus koppenhoeferi' in names

    # check minhash by doing a manual merge
    ss1 = next(iter(sourmash.load_file_as_signatures(db_path1)))
    ss2 = next(iter(sourmash.load_file_as_signatures(db_path2)))
    mh1 = ss1.minhash
    mh2 = ss2.minhash
    mh3 = mh1 + mh2

    species_mh = siglist[0].minhash
    assert species_mh == mh3


# remove a needed accession and check for failure
def test_createdb_6_xeno(runtmp):
    taxdb_path = utils.get_test_data('xeno-missing.lineages.csv')

    db_path1 = utils.get_test_data('xeno-1.sig')
    db_path2 = utils.get_test_data('xeno-2.sig')
    outfile = runtmp.output('output.sig.zip')

    args = FakeCreateDBArgs(sketches=[db_path1, db_path2],
                            taxonomy_file=[taxdb_path],
                            output=outfile, ksize=31)

    with pytest.raises(SystemExit) as e:
        pangenome_createdb_main(args)

    assert e.value.code == -1
