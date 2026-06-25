#!/usr/bin/env bash
###############################################################################
# 20_eklipse_splicebreak.sh — Python 2 / Java 8 toolbox
#
# eKLIPse (Python 2.7, BLAST+ >=2.3, circos) and Splice-Break2 (bash driver
# that shells out to /usr/bin/python [must be py2] and Java 8, with MapSplice2
# and bbmap bundled as ELF binaries) share one Python-2 conda environment.
###############################################################################
set -euxo pipefail

# Pin the last Python-2.7-compatible biopython/numpy/tqdm so the py2.7 solve is
# deterministic. These pins matter: recent tqdm noarch builds still declare
# `python >=2.7` but their code uses py3 f-strings (import SyntaxError under
# py2.7), and modern biopython has no py2.7 build. tqdm 4.62.x is the last line
# that still supports Python 2.7 (4.63.0 dropped it).
# blast/circos/samtools/openjdk are python-independent.
micromamba create -y -n py2tools \
    python=2.7 \
    biopython=1.76 \
    numpy=1.16 \
    tqdm=4.62 \
    'blast>=2.9' \
    circos \
    'samtools>=1.9,<1.16' \
    openjdk=8

# --- eKLIPse ---------------------------------------------------------------
git clone https://github.com/dooguypapua/eKLIPse /opt/eKLIPse
git -C /opt/eKLIPse checkout "${EKLIPSE_SHA}"
# The bundled rCRS GenBank annotation we rely on:
test -f /opt/eKLIPse/data/NC_012920.1.gb
# eKLIPse imports local modules (pybam/spinner/tabulate) -> run from its dir.
# Assert the Python deps AND the external binaries eKLIPse shells out to are
# present, so a missing tool fails the build rather than the first analysis run.
micromamba run -n py2tools python -c "from Bio import SeqIO; import tqdm; print('biopython ok')"
micromamba run -n py2tools bash -c 'command -v blastn makeblastdb circos samtools'

# --- Splice-Break2 ---------------------------------------------------------
git clone https://github.com/brookehjelm/Splice-Break2 /tmp/Splice-Break2-src
git -C /tmp/Splice-Break2-src checkout "${SPLICEBREAK_SHA}"
mkdir -p /opt/Splice-Break2
tar -xzf /tmp/Splice-Break2-src/Paired-End_Download/Splice-Break2-v3.0.2_PAIRED-END.tar.gz \
    -C /opt/Splice-Break2
rm -rf /tmp/Splice-Break2-src
# Verify the extracted layout the wrapper depends on.
test -f /opt/Splice-Break2/Splice-Break2-v3.0.2_PAIRED-END/Splice-Break2_paired-end.sh
test -f /opt/Splice-Break2/Splice-Break2-v3.0.2_PAIRED-END/NC_012920.1/NC.fa
chmod -R u+rwX /opt/Splice-Break2

# --- Patch CountBases.py (the inner mpileup base-counter) -------------------
# Upstream reads only a SINGLE digit of an indel's +/-<len><seq> length and
# indexes the read-base column one char at a time, so a multi-digit indel or an
# indel token at the END of a WGS-depth pileup column throws IndexError. That
# uncaught Python-2 exception aborts CountBases.py, which trips the inner
# Splice-Break2_0725.sh `if [ $? -gt 0 ]; then exit` guard BEFORE it writes any
# deletion table — yielding a 0-byte *_LargeMTDeletions_*.txt on every sample
# (the del4977 junction is in MapSplice's junctions.txt but never processed).
# The replacement reads the full multi-digit length and slices (bounds-safe);
# all non-indel counting is byte-for-byte the upstream behaviour.
SB_PE=/opt/Splice-Break2/Splice-Break2-v3.0.2_PAIRED-END
for cb in "$SB_PE/Splice-Break_v3.0.1/reference/CountBases.py" \
          "$SB_PE/Splice-Break_v3.0.1/ref_Nsub/CountBases.py"; do
    [[ -f "$cb" ]] || continue
    cat > "$cb" <<'PYEOF'
import sys
# NOTE: vendored Splice-Break2 helper, patched at image build
# (docker/install/20_eklipse_splicebreak.sh). Only the +/- indel parse changed:
# it now reads the FULL (multi-digit) indel length and slices the inserted bases
# (bounds-safe), so WGS-depth pileups no longer crash it with IndexError. All
# other counting is identical to upstream.
inFile = open(sys.argv[1],'r')

print 'bp\tA\tG\tC\tT\tdel\tins\tinserted\tambiguous'
for line in inFile:
        data = line.strip().split('\t')
        if len(data) < 5:
                continue
        bp = data[1]
        bases = data[4].upper()
        ref = data[2].upper()

        types = {'A':0,'G':0,'C':0,'T':0,'-':0,'+':[],'X':[]}

        i = 0
        while i < len(bases):
                base = bases[i]
                if base == '^' or base == '$':
                        i += 1
                elif base == '+' or base == '-':
                        i += 1
                        num = ''
                        while i < len(bases) and bases[i].isdigit():
                                num += bases[i]
                                i += 1
                        addNum = int(num) if num else 0
                        addSeq = bases[i:i+addNum]
                        i += addNum - 1
                        if base == '+':
                                types['+'].append(addSeq)
                elif base == '*':
                        types['-'] += 1
                elif base == '.' or base == ',':
                        types[ref] += 1
                else:
                        if types.has_key(base):
                                types[base] += 1
                        else:
                                types['X'].append(base)

                i += 1

        adds = '.'
        if len(types['+']) > 0:
                adds = ','.join(types['+'])

        amb = '.'
        if len(types['X']) > 0:
                amb = ','.join(types['X'])

        out = [bp,types['A'],types['G'],types['C'],types['T'],types['-'],len(types['+']),adds,amb]
        print '\t'.join([str(x) for x in out])
PYEOF
    # Assert the patch actually survives the crashing input class (a terminal
    # indel and a multi-digit indel length) and still emits rows — a single
    # /dev/null check would pass even if the patch silently regressed.
    printf 'chrM\t100\tA\t3\t..+1\tIII\nchrM\t101\tC\t4\t,,,+12ACGTACGTACGT\tIIII\n' > /tmp/_cb_test.pileup
    cb_lines="$(micromamba run -n py2tools python "$cb" /tmp/_cb_test.pileup | wc -l || echo 0)"
    [[ "$cb_lines" -ge 3 ]] || { echo "ERROR: patched CountBases.py crashed/empty on indel input ($cb)"; exit 1; }
    echo "patched + crash-tested $cb ($cb_lines lines)"
done
rm -f /tmp/_cb_test.pileup

# The bundled samtools v1.8 was built on CentOS and needs libcrypto.so.10
# (OpenSSL 1.0), absent on the Debian base, so it fails at exec — used by the
# driver (`samtools sort`) and the inner Splice-Break2_0725.sh (`samtools
# mpileup`), both via ${SB_Path}/samtools/v1.8/bin/samtools. Replace it with a
# wrapper to the py2tools env's samtools (1.9-1.15; sort/mpileup are
# interface-compatible with the classic pileup format Splice-Break parses).
SB_SAMTOOLS=/opt/Splice-Break2/Splice-Break2-v3.0.2_PAIRED-END/samtools/v1.8/bin/samtools
cat > "$SB_SAMTOOLS" <<'EOF'
#!/usr/bin/env bash
exec /opt/conda/envs/py2tools/bin/samtools "$@"
EOF
chmod +x "$SB_SAMTOOLS"
# Sanity: the replacement runs and is a 1.x samtools.
micromamba run -n py2tools "$SB_SAMTOOLS" --version | head -1

micromamba run -n py2tools java -version 2>&1 | head -1
micromamba clean -a -y
