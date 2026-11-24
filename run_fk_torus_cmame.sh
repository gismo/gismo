#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=48G
#SBATCH --account=compmech
#SBATCH --job-name=fk_torus
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE

# Fisher-Kolmogorov on the solid torus: body-fitted (quadratic, C1) vs
# immersed (Algoim cut cells), CMAME comparison.
# Parameters: Lorenzo et al. (2024) median D and rho, dt = 1 day (as in the
# paper), T = 730 days (2 years), seed R = 2.5 mm, sigma = 1 mm.
# Two lanes run AT THE SAME TIME, each as its own srun step with half the CPUs:
#   lane BF: body-fitted cases,  lane IM: immersed cases.
# Within a lane the cases run one after the other, LIGHT FIRST, HEAVY LAST.
# A failing case is reported and its lane moves on to the next case.
#
# This run (2026-11-24): the cases again EXCEPT bf_r5 (running in its own job),
# bf_r3_dt5, im_r5 (quadDepth 2) and the --cond cases (already done), now with
# --saveAt: every case writes its
# solution spline at 0, 365 and 730 days (<case>/N_t<time>.xml). After both
# lanes, the body-fitted code compares the saved fields (DSC, RMSE, e_L2)
# against bf_r4 -> comparison.csv. Needs the executables with --saveAt and
# --cmpA: rebuild on the server first.
#
# Linear solver: the codes' default --solver cgdiag (CG + Jacobi); add
# --solver direct to a case for the old SimplicialLDLT.
# Memory: immersed -r 5 needs ~2.5 GB (--quadDepth 2) or ~1.3 GB (--quadDepth 0);
# -r 6 with --quadDepth 0: 1.4 GB quadrature cache (measured).

set -eo pipefail

source /opt/intel/oneapi/setvars.sh --force
set -u

ROOT_DIR=/home/lucasventa/gismo_tumor
BF_EXE="$ROOT_DIR/build/bin/fischerkolmogrov_torus_example_bodyfitted"
IM_EXE="$ROOT_DIR/build/bin/fischerkolmogrov_torus_example_immersed"
# New folder: run_case deletes <case> before running, so the old results in
# 00_results_fk_torus_cmame stay untouched.
RESULTS_DIR="$ROOT_DIR/00_results_fk_torus_cmame_save"

NCPU=40              # whole job (= --cpus-per-task above)
LANE_CPU=20          # CPUs per lane (2 lanes)
BF_MEM=16G           # memory per srun step; the two must fit in --mem
IM_MEM=32G           # immersed -r 6 is the memory-hungry one

# Physics, identical for every case
PHYS=(-D 1.26e-3 -p 1.86e-3
      -R 2.5 -s 1 -a 0.5 -L 50
      --xc 7.88 --yc 25 --zc 27.12)
# Time stepping, passed per case (an option may be given only once)
DT1=(-t 1 -N 730)    # dt = 1 day, as in the paper, T = 730 days
DT5=(-t 5 -N 146)    # dt = 5 days, same T (time-step check)
DTC=(-t 1 -N 1)      # one step only: condition-number runs (--cond)
# Saved solutions for the field comparison (multiples of both dt = 1 and 5)
SAVE=(--saveAt 0,365,730)
SAVE_TIMES=(0 365 730)

for exe in "$BF_EXE" "$IM_EXE"; do
  if [[ ! -x "$exe" ]]; then
    echo "Error: executable not found or not executable: $exe" >&2
    exit 1
  fi
  if ! "$exe" --help 2>&1 | grep -q -- "--saveAt"; then
    echo "Error: $exe has no --saveAt option: rebuild it first." >&2
    exit 1
  fi
done

# run_case <name> <mem> <executable> [extra arguments...]
#   output folder: $RESULTS_DIR/<name>, log: <name>.log / <name>.err
#   --exact and an explicit --mem let two steps share the allocation
#   (without them the first step takes all CPUs/memory and the second waits).
run_case() {
  local name=$1
  local mem=$2
  local exe=$3
  shift 3

  echo "[$(date '+%F %T')] Starting $name with $LANE_CPU CPUs, $mem: $(basename "$exe") $*"
  rm -rf "$RESULTS_DIR/$name"
  if OMP_NUM_THREADS=$LANE_CPU srun --exact --ntasks=1 --cpus-per-task=$LANE_CPU \
       --mem="$mem" "$exe" \
       "${PHYS[@]}" "$@" \
       -o "$RESULTS_DIR/$name" \
       > "$RESULTS_DIR/$name.log" 2> "$RESULTS_DIR/$name.err"; then
    echo "[$(date '+%F %T')] Finished $name"
  else
    echo "[$(date '+%F %T')] FAILED $name (exit $?), see $RESULTS_DIR/$name.err" >&2
  fi
}

mkdir -p "$RESULTS_DIR"

# ---------------------------------------------------------------------------
# Lane BF: body-fitted, light -> heavy
# ---------------------------------------------------------------------------
lane_bf() {
  # condition numbers: already done (00_results_fk_torus_cmame/cond_*)
  # run_case cond_bf_r2 "$BF_MEM" "$BF_EXE" "${DTC[@]}" -r 2 --refTor 1 --cond
  # run_case cond_bf_r3 "$BF_MEM" "$BF_EXE" "${DTC[@]}" -r 3 --refTor 1 --cond
  # run_case cond_bf_r4 "$BF_MEM" "$BF_EXE" "${DTC[@]}" -r 4 --refTor 1 --cond

  run_case bf_r2      "$BF_MEM" "$BF_EXE" "${DT1[@]}" -r 2 --refTor 1 "${SAVE[@]}"
  run_case bf_r3      "$BF_MEM" "$BF_EXE" "${DT1[@]}" -r 3 --refTor 1 "${SAVE[@]}"
  # time-step check (dt = 5): already done, not rerun
  # run_case bf_r3_dt5  "$BF_MEM" "$BF_EXE" "${DT5[@]}" -r 3 --refTor 1 "${SAVE[@]}"
  run_case bf_r4     "$BF_MEM" "$BF_EXE" "${DT1[@]}" -r 4 --refTor 1 --plot "${SAVE[@]}"
  # bf_r5: NOT here, it runs in its own job
}

# ---------------------------------------------------------------------------
# Lane IM: immersed, light -> heavy
# ---------------------------------------------------------------------------
lane_im() {
  # condition numbers: already done (00_results_fk_torus_cmame/cond_*)
  # run_case cond_im_r4 "$IM_MEM" "$IM_EXE" "${DTC[@]}" -r 4 --cond
  # run_case cond_im_r5 "$IM_MEM" "$IM_EXE" "${DTC[@]}" -r 5 --quadDepth 0 --cond

  run_case im_r4      "$IM_MEM" "$IM_EXE" "${DT1[@]}" -r 4 "${SAVE[@]}"
  run_case im_r5_q0   "$IM_MEM" "$IM_EXE" "${DT1[@]}" -r 5 --quadDepth 0 "${SAVE[@]}"
  # im_r5 (quadDepth 2): not rerun, the paper uses im_r5_q0
  # run_case im_r5      "$IM_MEM" "$IM_EXE" "${DT1[@]}" -r 5 --plot "${SAVE[@]}"
  run_case im_r6_q0   "$IM_MEM" "$IM_EXE" "${DT1[@]}" -r 6 --quadDepth 0 "${SAVE[@]}"
}

# ---------------------------------------------------------------------------
# Field comparison (after both lanes): DSC, RMSE, e_L2 of <A> against <B>,
# one row per saved time in comparison.csv. Quadrature on the bf r4 mesh (the
# finest mesh involved) with 4 Gauss points per direction; takes seconds.
# ---------------------------------------------------------------------------
compare() {   # compare <A case> <A type bf|im> <B case> <B type bf|im>
  local a=$1 ta=$2 b=$3 tb=$4
  for t in "${SAVE_TIMES[@]}"; do
    local fa="$RESULTS_DIR/$a/N_t$t.xml" fb="$RESULTS_DIR/$b/N_t$t.xml"
    if [[ -f $fa && -f $fb ]]; then
      echo "[$(date '+%F %T')] Comparing $a vs $b at t = $t"
      OMP_NUM_THREADS=$NCPU "$BF_EXE" -L 50 -r 4 --refTor 1 --cmpQuB 2 \
        --cmpA "$fa" --cmpAType "$ta" --cmpB "$fb" --cmpBType "$tb" \
        --cmpOut "$RESULTS_DIR/comparison.csv" \
        || echo "Comparison FAILED: $a vs $b at t = $t" >&2
    else
      echo "Missing $fa or $fb: comparison skipped" >&2
    fi
  done
}

lane_cmp() {
  rm -f "$RESULTS_DIR/comparison.csv"
  compare bf_r4     bf bf_r4 bf   # sanity check: DSC = 1, RMSE = 0
  # convergence, both methods, against the reference bf_r4
  compare bf_r2     bf bf_r4 bf
  compare bf_r3     bf bf_r4 bf
  compare im_r4     im bf_r4 bf
  compare im_r5_q0  im bf_r4 bf
  compare im_r6_q0  im bf_r4 bf
  # matched-resolution pairs: immersed vs body-fitted at similar h
  compare im_r5_q0  im bf_r3 bf
  compare im_r4     im bf_r2 bf
}

lane_bf > "$RESULTS_DIR/lane_bf.txt" 2>&1 &
pid_bf=$!
lane_im > "$RESULTS_DIR/lane_im.txt" 2>&1 &
pid_im=$!
echo "[$(date '+%F %T')] Lanes started (BF pid $pid_bf, IM pid $pid_im); progress in lane_bf.txt / lane_im.txt"

wait "$pid_bf" || true
echo "[$(date '+%F %T')] Lane BF done"
wait "$pid_im" || true
echo "[$(date '+%F %T')] Lane IM done"

lane_cmp > "$RESULTS_DIR/lane_cmp.txt" 2>&1
echo "[$(date '+%F %T')] Comparisons done: $RESULTS_DIR/comparison.csv"
