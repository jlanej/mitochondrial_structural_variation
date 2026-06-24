#!/usr/bin/env python3
"""Unit tests for the caller output parsers.

These run WITHOUT the Docker image or any caller installed: they feed each
parser a literal sample of real caller output and assert the normalised record
is correct (in particular that the ~4977 bp common deletion is recognised).

Run:  python3 test/test_parsers.py        (no pytest needed)
  or: pytest test/test_parsers.py
"""
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline", "lib"))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
import parsers  # noqa: E402


def _tmp(content):
    fd, path = tempfile.mkstemp()
    with os.fdopen(fd, "w") as fh:
        fh.write(content)
    return path


# --------------------------------------------------------------------------- #
def test_is_common_deletion():
    assert parsers.is_common_deletion(8469, 13447)
    assert parsers.is_common_deletion(8482, 13446)      # breakpoints in repeat
    assert not parsers.is_common_deletion(5999, 10999)  # del6000
    assert not parsers.is_common_deletion(None, 13447)


def test_eklipse():
    # ";"-delimited, double-quoted, comma decimals (European locale)
    content = (
        '"Title";"5\' breakpoint";"3\' breakpoint";"Freq";"Freq For";"Freq Rev";'
        '"5\' Blast";"3\' Blast";"5\' Depth";"3\' Depth";"Repetition"\n'
        '"sampleA";"8469";"13447";"26,8";"24,1";"29,5";"2";"23";"1393";"412";'
        '"8470-GA-8471 | 13447-CT-13448"\n'
    )
    path = _tmp(content)
    recs = parsers.parse_eklipse(path, sample="sampleA")
    os.unlink(path)
    assert len(recs) == 1
    r = recs[0]
    assert r["caller"] == "eklipse" and r["sv_type"] == "deletion"
    assert r["bp5"] == 8469 and r["bp3"] == 13447
    assert r["svlen"] == 4978
    assert abs(r["het"] - 0.268) < 1e-6        # 26.8% -> 0.268
    assert parsers.is_common_deletion(r["bp5"], r["bp3"])


def test_mitosalt():
    # Real delplot.R header (delplot.R:454-460): cluster.id, del.size, event
    # values "del"/"dup", and a "seq" direct-repeat column (no "direct.repeat").
    content = (
        "sample\tcluster.id\talt.reads\tref.reads\theteroplasmy\tdel.start.range\t"
        "del.end.range\tdel.size\tfinal.event\tfinal.start\tfinal.end\tfinal.size\t"
        "seq1\tseq2\tseq\n"
        "S1\t1\t120\t300\t0.286\t8469-8470\t13446-13447\t4977\tdel\t8469\t"
        "13447\t4978\tACGTACGT\tTGCATGCA\tCTCATTCACC\n"
        "S1\t2\t40\t260\t0.133\t5999-6000\t10999-11000\t4999\tdup\t6000\t"
        "7000\t1000\tAAAA\tTTTT\t-\n"
    )
    path = _tmp(content)
    recs = parsers.parse_mitosalt(path, sample="S1")
    os.unlink(path)
    assert len(recs) == 2
    d = recs[0]
    assert d["bp5"] == 8469 and d["bp3"] == 13447 and d["svlen"] == 4978
    assert abs(d["het"] - 0.286) < 1e-6 and d["support"] == 120
    assert d["sv_type"] == "deletion"
    assert d["extra"] == "repeat=CTCATTCACC"      # read from the real "seq" column
    assert parsers.is_common_deletion(d["bp5"], d["bp3"])
    assert recs[1]["sv_type"] == "duplication"


def test_splicebreak2():
    content = (
        "Sample_ID  Reference_Genome  MapSplice_Breakpoint  5'_Break  3'_Break  "
        "Deletion_Size_bp  Deletion_Reads  Benchmark_Coverage  Deletion_Read_%  "
        "Annotation  Left_Overhang  Right_Overhang\n"
        "S1  NC_012920.1_rCRS  8469-13447  8469  13447  4978  473  5000  9.46  "
        "high_frequency_deletion  AAAA  TTTT\n"
    )
    path = _tmp(content)
    recs = parsers.parse_splicebreak2(path)
    os.unlink(path)
    assert len(recs) == 1
    r = recs[0]
    assert r["sample"] == "S1"
    assert r["bp5"] == 8469 and r["bp3"] == 13447 and r["svlen"] == 4978
    assert r["support"] == 473 and abs(r["het"] - 0.0946) < 1e-6
    assert parsers.is_common_deletion(r["bp5"], r["bp3"])


def test_mitomut():
    content = (
        "Total Reads\tS1 Reads\tS2 Reads\tStart\tEnd\tHeteroplasmy Level\n"
        "133\t64\t69\t8469\t13447\t0.0751442329808187\n"
    )
    path = _tmp(content)
    recs = parsers.parse_mitomut(path, sample="S1")
    os.unlink(path)
    assert len(recs) == 1
    r = recs[0]
    assert r["bp5"] == 8469 and r["bp3"] == 13447 and r["svlen"] == 4978
    assert r["support"] == 133 and abs(r["het"] - 0.0751442329808187) < 1e-9
    assert parsers.is_common_deletion(r["bp5"], r["bp3"])


def test_mitohpc():
    # MitoHPC sv.tab (the reference caller); a real row from example/sv.tab.
    content = (
        "sample\tchrom\tpos_bp5\tend_bp3\tsvlen\tsvclaim\tjr\tsr\taf_junction\t"
        "af_coverage\tafdiff\tcvgr\tflank_dp\thomlen\thomseq\tdelclass\tcommon\t"
        "ngene\tgene_list\thgvs\tfilter\tflags\tsrcons\tsrsb\tjsup\tsvconf\t"
        "svimpact\tsvimpact_band\n"
        "sv_del4977_h30\tchrM\t8482\t13446\t4964\tDJ\t154\t401\t0.277\t0.268\t"
        "0.010\t0.732\t575\t13\tACCTCCCTCACCA\tI\t1\t12\tATP8:P\t"
        "NC_012920.1:m.8483_13446del\tPASS\tREPEAT\t1.000\t0.496\tHIGH\t83\t85\t"
        "SEVERE\n")
    path = _tmp(content)
    recs = parsers.parse_mitohpc(path, sample="sv_del4977_h30")
    os.unlink(path)
    assert len(recs) == 1
    r = recs[0]
    assert r["caller"] == "mitohpc" and r["sv_type"] == "deletion"
    assert r["bp5"] == 8482 and r["bp3"] == 13446 and r["svlen"] == 4964
    assert abs(r["het"] - 0.268) < 1e-6 and r["support"] == 154
    assert "filter=PASS" in r["extra"] and "flags=REPEAT" in r["extra"]
    assert parsers.is_common_deletion(r["bp5"], r["bp3"])


def test_mitohpc_first_in_callers():
    assert parsers.CALLERS[0] == "mitohpc" and len(parsers.CALLERS) == 6


def test_mitoseek_discordant():
    content = (
        "#MitoChr\tMitoPos\tVarChr\tVarPos\tSupportedReads\tMeanMappingQuality\t"
        "Mito.gene\tMito.genedetail\n"
        "MT\t12491\t18\t57409214\t4\t37.000\tND5\tMTND5:mRNA:12338-14149:(+)\n"
    )
    path = _tmp(content)
    recs = parsers.parse_mitoseek_discordant(path, sample="S1")
    os.unlink(path)
    assert len(recs) == 1
    assert recs[0]["sv_type"] == "breakpoint" and recs[0]["bp5"] == 12491
    assert recs[0]["support"] == 4


def test_mitoseek_large_deletion():
    # Three read pairs spanning the common deletion: left mate POS 8369-8371
    # (100M -> 5' break ~8469), mate (PNEXT) at 13447, large TLEN.
    lines = []
    for i, pos in enumerate((8369, 8370, 8371)):
        lines.append(
            "r%d\t99\tchrM\t%d\t60\t100M\t=\t13447\t5178\t%s\t%s"
            % (i, pos, "A" * 100, "I" * 100))
    content = "\n".join(lines) + "\n"
    path = _tmp(content)
    recs = parsers.parse_mitoseek_large_deletion(path, sample="S1",
                                                 min_support=3, bin_size=25)
    os.unlink(path)
    assert len(recs) == 1, recs
    r = recs[0]
    assert r["sv_type"] == "deletion" and r["support"] == 3.0
    assert parsers.is_common_deletion(r["bp5"], r["bp3"]), (r["bp5"], r["bp3"])


def test_parse_sample_dir_and_postprocess(tmp_path=None):
    import shutil
    root = tempfile.mkdtemp()
    try:
        sd = os.path.join(root, "S1")
        for c in ("eklipse", "mitomut"):
            os.makedirs(os.path.join(sd, c))
        with open(os.path.join(sd, "status.tsv"), "w") as fh:
            fh.write("caller\tstatus\tseconds\n")
        with open(os.path.join(sd, "eklipse", "eKLIPse_deletions.csv"), "w") as fh:
            fh.write('"Title";"5\' breakpoint";"3\' breakpoint";"Freq";"Freq For";'
                     '"Freq Rev";"5\' Blast";"3\' Blast";"5\' Depth";"3\' Depth";'
                     '"Repetition"\n'
                     '"S1";"8469";"13447";"26,8";"";"";"2";"23";"1393";"412";"x"\n')
        with open(os.path.join(sd, "mitomut", "mitomut_results.txt"), "w") as fh:
            fh.write("Total Reads\tS1 Reads\tS2 Reads\tStart\tEnd\tHeteroplasmy Level\n"
                     "133\t64\t69\t8469\t13447\t0.075\n")
        recs = parsers.parse_sample_dir(sd, "S1")
        assert len(recs) == 2

        # Exercise the postprocess CLI end-to-end.
        sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
        import postprocess
        rc = postprocess.main(["--root", root, "--out-dir", root])
        assert rc == 0
        for f in ("cohort_sv_calls.tsv", "cohort_common_deletion.tsv",
                  "cohort_caller_matrix.tsv", "cohort_summary.txt"):
            assert os.path.isfile(os.path.join(root, f)), f
        with open(os.path.join(root, "cohort_common_deletion.tsv")) as fh:
            body = fh.read()
        assert "S1\teklipse\t1" in body and "S1\tmitomut\t1" in body
    finally:
        shutil.rmtree(root, ignore_errors=True)


def _main():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in fns:
        try:
            fn()
            print("PASS", fn.__name__)
        except AssertionError as e:
            failed += 1
            print("FAIL", fn.__name__, "->", e)
        except Exception as e:  # noqa: BLE001
            failed += 1
            print("ERROR", fn.__name__, "->", repr(e))
    print("\n%d/%d passed" % (len(fns) - failed, len(fns)))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(_main())
