#!/usr/bin/env bash
# Check that the .h5 rustar-aligner writes is a valid HDF5 file, using HDF5's
# own C library rather than our own reader.
#
# The in-tree writer (`src/io/hdf5.rs`) has no dependency on libhdf5, which is
# the point of it; the flip side is that `cargo test` cannot tell whether what
# it produced is a conformant file. Only libhdf5 can say that, so this check
# lives here rather than in the test suite, alongside the differential tests
# that also need an external oracle.
#
#   test/h5_conformance.sh <path/to/some.h5>
#
# Needs the HDF5 command-line tools (`brew install hdf5`, `apt install
# hdf5-tools`) and, for the reader checks, python3 with h5py and scipy.
set -uo pipefail

F=${1:?usage: h5_conformance.sh <file.h5>}
fail=0
step() { printf '\n== %s\n' "$1"; }
ok()   { printf '   OK   %s\n' "$1"; }
bad()  { printf '   FAIL %s\n' "$1"; fail=1; }

command -v h5dump >/dev/null || { echo "h5dump not found; install the HDF5 tools" >&2; exit 2; }

step "libhdf5 opens it and walks the whole tree"
h5ls -r "$F" && ok "h5ls -r" || bad "h5ls -r"
h5dump -H "$F" >/dev/null && ok "h5dump -H" || bad "h5dump -H"
h5stat "$F" >/dev/null && ok "h5stat" || bad "h5stat"

step "libhdf5 rewrites it from scratch and finds no difference"
# h5repack reads the file with the C library and writes a fresh one; h5diff then
# compares them object by object. Passing means our bytes carry exactly the
# content libhdf5 would have written itself.
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
if h5repack "$F" "$tmp/repacked.h5"; then
  ok "h5repack"
  if h5diff "$tmp/repacked.h5" "$F"; then ok "h5diff (no difference)"; else bad "h5diff"; fi
else
  bad "h5repack"
fi

step "the readers this file exists for"
# Optional: these need h5py and scipy. Point PYTHON at an interpreter that has
# them if the default one does not, rather than treating their absence as a
# failure of the file.
PYTHON=${PYTHON:-python3}
if ! "$PYTHON" -c "import h5py, scipy" 2>/dev/null; then
  printf '   SKIP %s\n' "$PYTHON has no h5py/scipy; set PYTHON to an interpreter that does"
else
"$PYTHON" - "$F" <<'PY' && ok "python readers" || bad "python readers"
import sys
import h5py
import scipy.sparse as sp
import numpy as np

path = sys.argv[1]

# These two readers only apply to a feature-barcode matrix. The writer is also
# used for other trees (its own golden test file, for one), and those are
# checked by the libhdf5 steps above.
with h5py.File(path, 'r') as f:
    if '/matrix' not in f or 'barcodes' not in f['matrix']:
        print("   not a feature-barcode matrix; skipping the 10x readers")
        sys.exit(0)

# scanpy's _read_v3_10x_h5, transcribed.
with h5py.File(path, 'r') as f:
    assert '/matrix' in f, "no /matrix group: scanpy would treat this as v2"
    d = {}
    f['matrix'].visititems(
        lambda name, obj: d.__setitem__(name.split('/')[-1], obj[...])
        if isinstance(obj, h5py.Dataset) else None)
    M, N = d['shape']
    csr = sp.csr_matrix((d['data'].astype('float32'), d['indices'], d['indptr']),
                        shape=(N, M))
    barcodes = d['barcodes'].astype(str)
    ids = d['id'].astype(str)
    print(f"   scanpy    : {csr.shape} cells x genes, nnz={csr.nnz}, "
          f"sum={int(csr.sum())}, first barcode {barcodes[0]}, first gene {ids[0]}")

# CellBender's get_matrix_from_cellranger_h5, transcribed.
with h5py.File(path, 'r') as f:
    assert 'matrix' in f.keys(), "CellBender would detect this as CellRanger v2"
    g = f['matrix']
    d = {}
    for k in g.keys():
        if k == 'features':
            for k2 in g['features'].keys():
                d[k2] = np.array(g['features'][k2])
        else:
            d[k] = np.array(g[k])
    csc = sp.csc_matrix((d['data'], d['indices'], d['indptr']), shape=d['shape'])
    print(f"   cellbender: {csc.transpose().tocsr().shape} cells x genes, "
          f"nnz={csc.nnz}, sum={int(csc.sum())}")
PY
fi

printf '\n'
if [ "$fail" = 0 ]; then echo "PASS: $F"; else echo "FAIL: $F"; fi
exit "$fail"
