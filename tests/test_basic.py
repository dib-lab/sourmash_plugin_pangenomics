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
def test_createdb_1(runtmp):
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
def test_createdb_2(runtmp):
    taxdb_path = get_workflow_data('gtdb-rs214-agatha.lineages.csv.gz')
    print(taxdb_path)

    db_path = get_workflow_data('gtdb-rs214-agatha-k21.zip')
    print(db_path)

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
