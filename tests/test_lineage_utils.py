#! /usr/bin/env python

from sourmash.tax.tax_utils import LineagePair, LINLineageInfo #BaseLineageInfo


def build_test_lineage():
    lineage = [
        LineagePair("superkingdom", "Bacteria"),
        LineagePair("phylum", "Firmicutes"),
        LineagePair("class", "Bacilli"),
        LineagePair("order", "Lactobacillales"),
        LineagePair("family", "Lactobacillaceae"),
        LineagePair("genus", "Lactobacillus"),
        LineagePair("species", "Lactobacillus_acidophilus"),
        LineagePair("sub_species", "Lactobacillus_acidophilus_NCFM"),
    ]
    return LINLineageInfo(lineage)


def test_lineage_methods():
    lin = build_test_lineage()

    print("\nlineage_at_rank('superkingdom'):", lin.lineage_at_rank("superkingdom"))
    print("\nlineage_above_rank('superkingdom'):", lin.lineage_above_rank("superkingdom"))
    print("\nlineage_below_rank('superkingdom'):", lin.lineage_below_rank("superkingdom"))


    print("\nlineage_at_rank('class'):", lin.lineage_at_rank("class"))
    print("\nlineage_above_rank('class'):", lin.lineage_above_rank("class"))
    print("\nlineage_below_rank('class'):", lin.lineage_below_rank("class"))

    print("\nlineage_at_rank('species'):", lin.lineage_at_rank("species"))
    print("\nlineage_above_rank('species'):", lin.lineage_above_rank("species"))
    print("\nlineage_below_rank('species'):", lin.lineage_below_rank("species"))

    print("\nlineage_at_rank('sub_species'):", lin.lineage_at_rank("sub_species"))
    print("\nlineage_above_rank('sub_species'):", lin.lineage_above_rank("sub_species"))
    print("\nlineage_below_rank('sub_species'):", lin.lineage_below_rank("sub_species"))

    try:
        lin.lineage_at_rank("flavor")
    except ValueError as e:
        print("\nCaught expected error:", e)


if __name__ == "__main__":
    test_lineage_methods()
