#!/usr/bin/env bash
# Paper Fig. 7: curvilinear SGBC coated cylinder, then the digitized-curve check.
set -o pipefail
cd /home/alejandro/workspace/dgtd

export OMP_NUM_THREADS=1
LOG=exports/paper_fig7_g2.log
CASE=testData/maxwellInputs/2D_RCS_SGBC_Circle_G2_1m_LowFreq/2D_RCS_SGBC_Circle_G2_1m_LowFreq.json
BIN=build/gnu-release-mpi/bin

echo "=== paper fig7 G2 start $(date -Is) ===" | tee "$LOG"
stdbuf -oL -eL mpirun -np 12 "$BIN/opensemba_dgtd" --restart -i "$CASE" \
  2>&1 | tee -a "$LOG"
dgtd_rc=${PIPESTATUS[0]}
echo "DGTD_EXIT:${dgtd_rc}" | tee -a "$LOG"
if [[ ${dgtd_rc} -ne 0 ]]; then
  echo "PAPER_FIG7_DONE dgtd_failed" | tee -a "$LOG"
  exit "${dgtd_rc}"
fi

echo "=== RCS start $(date -Is) ===" | tee -a "$LOG"
stdbuf -oL -eL "$BIN/opensemba_rcs" \
  -i testData/rcsInputs/2D_RCS_SGBC_Circle_G2_1m_LowFreq.rcs.json \
  2>&1 | tee -a "$LOG"
rcs_rc=${PIPESTATUS[0]}
echo "RCS_EXIT:${rcs_rc}" | tee -a "$LOG"
if [[ ${rcs_rc} -ne 0 ]]; then
  echo "PAPER_FIG7_DONE rcs_failed" | tee -a "$LOG"
  exit "${rcs_rc}"
fi

DAT=exports/SimulationData/mpi-12/2D_RCS_SGBC_Circle_G2_1m_LowFreq/RCSSurface/cylinder_rcs/rcs/rcsData_Th_1.570796_Phi_3.141593_dgtd.dat
echo "=== figure check $(date -Is) ===" | tee -a "$LOG"
python3 testData/maxwellInputs/2D_RCS_SGBC_Circle_G2_1m_LowFreq/paper_fig7_check.py "$DAT" \
  2>&1 | tee -a "$LOG"
check_rc=${PIPESTATUS[0]}
echo "CHECK_EXIT:${check_rc}" | tee -a "$LOG"
echo "PAPER_FIG7_DONE check_rc=${check_rc}" | tee -a "$LOG"
exit "${check_rc}"
