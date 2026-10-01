#!/usr/bin/env bash
#
# refine_tetmesh.sh -- offline uniform 1->8 tet refinement via the Gmsh CLI.
#
# Usage:
#   refine_tetmesh.sh [--verify] <in.msh> <L> [<outdir>]
#
# Applies Gmsh's uniform red refinement (each tet split into 8, each boundary
# triangle into 4) to an ASCII MSH 4.1 linear-tet mesh, L times, writing
# <stem>_L1.msh ... <stem>_L<L>.msh (stem = <in.msh> basename without the
# extension). Level k is produced from level k-1 (level 1 from the input);
# an L=2 run therefore leaves both _L1 and _L2 on disk.
#
# ONLY the .msh is ever touched -- the .geo is never re-read or re-meshed.
# Re-meshing from the .geo would put new boundary nodes back on the CAD
# surface (e.g. the sphere), changing the polyhedron a curved case encloses;
# Gmsh's "-refine" on a discrete .msh instead places new nodes at edge
# midpoints, which is what keeps every refined level exactly enclosing the
# original polyhedron the P2 mesh-fidelity gate compares against.
#
# Measured with Gmsh 4.12.1 on tetmesh_sphere.msh (the discriminating case,
# since a curved boundary is the only one a midpoint-snap would visibly move):
# the bare command
#   gmsh <in>.msh -refine -format msh41 -setnumber Mesh.Binary 0 -o <out>.msh
# runs in batch mode and exits 0 (no GUI, no prompt), places new nodes at
# edge midpoints (not reparametrised onto the sphere) and leaves every
# original node bit-identical -- sphere L=1: tets 6526->52208 (x8 exactly),
# boundary triangles 1504->6016 (x4 exactly), V_relerr=0.000e+00 and
# A_relerr=0.000e+00 (both within 1e-13), orig_nodes_missing=0, element
# types restricted to {1,2,4,15} (linear only). No extra option (e.g.
# Mesh.SecondOrderLinear) was needed.
#
# Output naming: <outdir>/<stem>_L<k>.msh plus a <outdir>/<stem>_L<k>.gmsh.log
# per level. Default outdir is <build_dir>/tetmesh_refined, where build_dir is
# read from <repo>/.claude/gismo-dev.local.json; an explicit <outdir> argument
# overrides it. The refined files are deliberately never committed -- they
# live under the (gitignored) build directory, or an outdir the caller chose,
# and refine_tetmesh.sh refuses any outdir that resolves under
# $REPO_ROOT/filedata to keep that true by construction.
#
# Size note: tet count scales as 8^L. The sphere mesh has ~6.5k tets, so
# L=2 is ~4e5 tets and L=3 is ~8x that again (~3.3e6) -- expect the L=3 gmsh
# pass and, under --verify, the Python summation over it, to take noticeably
# longer than L<=2.
#
# --verify parses the original and every produced level with an embedded
# Python reader that mirrors examples/gsTetMeshClip.h:readMsh41 one check at
# a time (cited inline below), then checks level k against the ORIGINAL:
#   - exact tet count ratio 8^k and boundary-triangle count ratio 4^k;
#   - mesh volume and boundary area (fsum-summed, the same quantities as
#     meshVolumeExact/unclippedBoundaryArea, gsTetMeshClip.h:925-1022)
#     agree to <= 1e-13 relative;
#   - every original node coordinate triple is present, unmoved, among the
#     level's node coordinates (orig_nodes_missing == 0);
#   - every element type is in readMsh41's allowed set {15,1,2,4}.
# It exits non-zero if any level fails.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: refine_tetmesh.sh [--verify] <in.msh> <L> [<outdir>]

  <in.msh>   ASCII MSH 4.1 linear-tet mesh (e.g. filedata/volumes/tetmesh_sphere.msh)
  <L>        number of uniform 1->8 refinement levels to apply (positive integer)
  <outdir>   optional output directory; default <build_dir>/tetmesh_refined,
             build_dir read from <repo>/.claude/gismo-dev.local.json

  --verify   after writing all levels, verify each against the original
             (tet/boundary-triangle count ratios, volume/area, unmoved nodes,
             allowed element types) and exit non-zero on any failure.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

VERIFY=0
if [[ "${1:-}" == "--verify" ]]; then
  VERIFY=1
  shift
fi

if [[ $# -lt 2 || $# -gt 3 ]]; then
  usage
  exit 1
fi

IN_MSH="$1"
L="$2"
OUTDIR_ARG="${3:-}"

if [[ ! -f "$IN_MSH" ]]; then
  echo "refine_tetmesh.sh: input file not found: $IN_MSH" >&2
  exit 1
fi
case "$IN_MSH" in
  *.msh) ;;
  *) echo "refine_tetmesh.sh: input file must end in .msh: $IN_MSH" >&2; exit 1 ;;
esac
if ! [[ "$L" =~ ^[1-9][0-9]*$ ]]; then
  echo "refine_tetmesh.sh: L must be a positive integer, got: $L" >&2
  exit 1
fi
if ! command -v python3 >/dev/null 2>&1; then
  echo "refine_tetmesh.sh: python3 not found on PATH" >&2
  exit 1
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"

if [[ -n "$OUTDIR_ARG" ]]; then
  OUTDIR="$OUTDIR_ARG"
else
  CONFIG="$REPO_ROOT/.claude/gismo-dev.local.json"
  if [[ ! -f "$CONFIG" ]]; then
    echo "refine_tetmesh.sh: no output dir given and config file not found: $CONFIG" >&2
    exit 1
  fi
  BUILD_DIR="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["build_dir"])' "$CONFIG")"
  OUTDIR="$BUILD_DIR/tetmesh_refined"
fi

# Resolve (without requiring the directory to exist yet) and refuse anything
# under filedata/ BEFORE creating it, so a refused call never touches the
# committed tree.
OUTDIR_RESOLVED="$(realpath -m -- "$OUTDIR")"
case "$OUTDIR_RESOLVED" in
  "$REPO_ROOT/filedata"|"$REPO_ROOT/filedata"/*)
    echo "refine_tetmesh.sh: refusing output dir under \$REPO_ROOT/filedata: $OUTDIR_RESOLVED" >&2
    exit 1
    ;;
esac
mkdir -p "$OUTDIR_RESOLVED"
OUTDIR="$OUTDIR_RESOLVED"

GMSH="${GMSH:-/usr/bin/gmsh}"

BASENAME="$(basename "$IN_MSH")"
STEM="${BASENAME%.msh}"

LEVEL_FILES=()
PREV="$IN_MSH"
for ((k = 1; k <= L; ++k)); do
  OUT="$OUTDIR/${STEM}_L${k}.msh"
  LOG="$OUTDIR/${STEM}_L${k}.gmsh.log"
  START=$SECONDS
  rm -f -- "$OUT"
  if ! timeout 600 "$GMSH" "$PREV" -refine -format msh41 -setnumber Mesh.Binary 0 -o "$OUT" >"$LOG" 2>&1; then
    rm -f -- "$OUT"
    echo "refine_tetmesh.sh: gmsh -refine failed at level $k (see $LOG)" >&2
    exit 1
  fi
  if [[ ! -s "$OUT" ]]; then
    echo "refine_tetmesh.sh: gmsh produced a missing or empty output at level $k: $OUT" >&2
    exit 1
  fi
  ELAPSED=$((SECONDS - START))
  echo "refine_tetmesh.sh: L=$k file=$OUT time=${ELAPSED}s"
  LEVEL_FILES+=("$OUT")
  PREV="$OUT"
done

if [[ "$VERIFY" -eq 1 ]]; then
  python3 - "$IN_MSH" "${LEVEL_FILES[@]}" <<'PY'
# Mirrors examples/gsTetMeshClip.h:readMsh41 (:578-736) check by check, then
# compares each refined level against the original mesh it was refined from.
import sys, os, math
from collections import defaultdict

ALLOWED_ETYPES = {15, 1, 2, 4}


def sub3(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def cross3(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def dot3(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def det3(a, b, c):
    return dot3(a, cross3(b, c))


def norm3(a):
    return math.sqrt(dot3(a, a))


def parse_msh41(path):
    with open(path, 'r') as f:
        text = f.read().replace('\r\n', '\n')
    lines = text.split('\n')
    n = len(lines)
    pos = [0]

    def next_line():
        if pos[0] >= n:
            raise RuntimeError('%s: unexpected end of file' % path)
        line = lines[pos[0]]
        pos[0] += 1
        return line

    # readMsh41:584-585 -- scan forward for the literal $MeshFormat line.
    while pos[0] < n and lines[pos[0]] != '$MeshFormat':
        pos[0] += 1
    if pos[0] >= n:
        raise RuntimeError('%s: no $MeshFormat section' % path)
    pos[0] += 1

    # readMsh41:589-595 -- version must be 4.1, filetype 0 (ASCII).
    verline = next_line().split()
    version, filetype = float(verline[0]), int(verline[1])
    if version != 4.1:
        raise RuntimeError('%s: MeshFormat version %s (need 4.1)' % (path, version))
    if filetype != 0:
        raise RuntimeError('%s: binary MSH file (filetype %d)' % (path, filetype))
    next_line()  # $EndMeshFormat

    nodes = []
    tag2idx = {}
    tets = []
    etypes = set()
    type2_count = 0

    while pos[0] < n:
        line = next_line()
        if line == '$Nodes':
            # readMsh41:606-634
            hdr = next_line().split()
            num_blocks, num_nodes = int(hdr[0]), int(hdr[1])
            nodes = [None] * num_nodes
            next_idx = 0
            for _ in range(num_blocks):
                bh = next_line().split()
                param, cnt = int(bh[2]), int(bh[3])
                if param != 0:
                    raise RuntimeError('%s: parametric node block is not supported' % path)
                tags = [int(next_line()) for _ in range(cnt)]
                for kk in range(cnt):
                    x, y, z = next_line().split()
                    nodes[next_idx] = (float(x), float(y), float(z))
                    tag2idx[tags[kk]] = next_idx
                    next_idx += 1
            next_line()  # $EndNodes
            continue
        if line == '$Elements':
            # readMsh41:641-669
            hdr = next_line().split()
            num_blocks = int(hdr[0])
            for _ in range(num_blocks):
                bh = next_line().split()
                dim, etype, cnt = int(bh[0]), int(bh[2]), int(bh[3])
                etypes.add(etype)
                if etype not in ALLOWED_ETYPES:
                    raise RuntimeError('%s: disallowed element type %d' % (path, etype))
                if dim == 3 and etype != 4:
                    raise RuntimeError('%s: dim==3 block has non-tet type %d' % (path, etype))
                if etype == 2:
                    type2_count += cnt
                if etype != 4:
                    for _ in range(cnt):
                        next_line()
                    continue
                for _ in range(cnt):
                    parts = next_line().split()
                    etag = int(parts[0])
                    node_tags = (int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4]))
                    idxs = []
                    for nt in node_tags:
                        if nt not in tag2idx:
                            raise RuntimeError('%s: element %d references unknown node tag %d' % (path, etag, nt))
                        idxs.append(tag2idx[nt])
                    tets.append(tuple(idxs))
            next_line()  # $EndElements
            continue
        if line.startswith('$') and not line.startswith('$End'):
            end_tag = '$End' + line[1:]
            while pos[0] < n and lines[pos[0]] != end_tag:
                pos[0] += 1
            if pos[0] < n:
                pos[0] += 1
            continue

    # readMsh41:679
    if not tets:
        raise RuntimeError('%s: no type-4 tetrahedra found' % path)

    # readMsh41:681-686 -- zero-volume tet check.
    for t in tets:
        p0, p1, p2, p3 = nodes[t[0]], nodes[t[1]], nodes[t[2]], nodes[t[3]]
        d = det3(sub3(p1, p0), sub3(p2, p0), sub3(p3, p0))
        if d == 0.0:
            raise RuntimeError('%s: tet with zero volume' % path)

    # readMsh41:694-733 -- boundary = faces owned by exactly one tet;
    # multiplicity > 2 is non-manifold.
    faces = defaultdict(list)
    for ti, t in enumerate(tets):
        for i_excl in range(4):
            f = tuple(t[j] for j in range(4) if j != i_excl)
            faces[tuple(sorted(f))].append((f, t[i_excl]))
    bdr_tris = []
    for recs in faces.values():
        if len(recs) > 2:
            raise RuntimeError('%s: non-manifold face' % path)
        if len(recs) == 1:
            (a, b, c), opp = recs[0]
            pa, pb, pc, po = nodes[a], nodes[b], nodes[c], nodes[opp]
            orient = dot3(cross3(sub3(pb, pa), sub3(pc, pa)), sub3(po, pa))
            if orient == 0.0:
                raise RuntimeError('%s: degenerate boundary-face orientation' % path)
            if orient > 0.0:
                b, c = c, b
            bdr_tris.append((a, b, c))

    return {
        'nodes': nodes, 'tag2idx': tag2idx, 'tets': tets, 'bdr_tris': bdr_tris,
        'etypes': etypes, 'type2_count': type2_count,
    }


def mesh_volume(data):
    # gsTetMeshClip.h:925-931 -- fsum(|det(P1-P0,P2-P0,P3-P0)|)/6.
    nodes = data['nodes']
    return math.fsum(
        abs(det3(sub3(nodes[t[1]], nodes[t[0]]), sub3(nodes[t[2]], nodes[t[0]]), sub3(nodes[t[3]], nodes[t[0]])))
        for t in data['tets']
    ) / 6.0


def mesh_area(data):
    # gsTetMeshClip.h:933-942 -- fsum(0.5*|cross(B-A,C-A)|).
    nodes = data['nodes']
    return math.fsum(
        0.5 * norm3(cross3(sub3(nodes[t[1]], nodes[t[0]]), sub3(nodes[t[2]], nodes[t[0]])))
        for t in data['bdr_tris']
    )


def fmt_ratio(n, n0):
    if n0 and n % n0 == 0:
        return str(n // n0)
    return '%.6f' % (n / n0 if n0 else float('nan'))


def main():
    argv = sys.argv[1:]
    if len(argv) < 2:
        print('refine_verify: FAIL usage: python3 - <orig> <L1> ... <LL>')
        return 1

    orig_path, level_paths = argv[0], argv[1:]
    orig = parse_msh41(orig_path)
    n0, m0 = len(orig['tets']), len(orig['bdr_tris'])
    v0, a0 = mesh_volume(orig), mesh_area(orig)
    orig_coords = set(orig['nodes'])

    all_pass = True

    def emit(k, path, data):
        nonlocal all_pass
        base = os.path.basename(path)
        n, m = len(data['tets']), len(data['bdr_tris'])
        expect_tet, expect_bdr = 8 ** k, 4 ** k
        v, a = mesh_volume(data), mesh_area(data)
        v_relerr = abs(v - v0) / abs(v0) if v0 != 0 else abs(v - v0)
        a_relerr = abs(a - a0) / abs(a0) if a0 != 0 else abs(a - a0)
        level_coords = set(data['nodes'])
        missing = sum(1 for c in orig_coords if c not in level_coords)
        etypes_ok = data['etypes'] <= ALLOWED_ETYPES
        etypes_str = ','.join(str(x) for x in sorted(data['etypes']))

        ok = (n == expect_tet * n0 and m == expect_bdr * m0 and
              v_relerr <= 1e-13 and a_relerr <= 1e-13 and
              missing == 0 and etypes_ok)
        if not ok:
            all_pass = False
        status = 'PASS' if ok else 'FAIL'
        print('refine_verify: L=%d file=%s tets=%d tet_ratio=%s expect=%d '
              'bdr_tris=%d bdr_ratio=%s expect=%d type2_tris=%d '
              'V=%.17g V_relerr=%.3e A=%.17g A_relerr=%.3e '
              'orig_nodes_missing=%d etypes=%s %s'
              % (k, base, n, fmt_ratio(n, n0), expect_tet,
                 m, fmt_ratio(m, m0), expect_bdr, data['type2_count'],
                 v, v_relerr, a, a_relerr, missing, etypes_str, status))
        if missing > 0:
            common = set(orig['tag2idx']) & set(data['tag2idx'])
            max_shift = 0.0
            for tag in common:
                d = norm3(sub3(orig['nodes'][orig['tag2idx'][tag]], data['nodes'][data['tag2idx'][tag]]))
                max_shift = max(max_shift, d)
            print('refine_verify: L=%d file=%s diagnostic max_shift_common_tags=%.6e' % (k, base, max_shift))

    emit(0, orig_path, orig)
    for k, path in enumerate(level_paths, start=1):
        emit(k, path, parse_msh41(path))

    print('refine_verify: ALL PASS' if all_pass else 'refine_verify: FAIL')
    return 0 if all_pass else 1


sys.exit(main())
PY
fi
