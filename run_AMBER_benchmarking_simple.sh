#!/bin/bash
set -e

# Usage: ./run_AMBER_benchmarking_simple.sh <molecule_name>
# Example: ./run_AMBER_benchmarking_simple.sh ethane
#          ./run_AMBER_benchmarking_simple.sh propane
MOLECULE="${1:-ethane}"

# AMBER setup. Override AMBERHOME in your shell to point at a local install;
# otherwise default to the ambertools conda env.
: "${AMBERHOME:=$HOME/miniconda3/envs/ambertools}"
export AMBERHOME
export PATH=$AMBERHOME/bin:$PATH
source $AMBERHOME/amber.sh

# Directories
INPUT="input_molecules_copy/${MOLECULE}.xyz"
WORKDIR="amber_runs/${MOLECULE}"

if [[ ! -f "$INPUT" ]]; then
  echo "Error: input file '$INPUT' not found."
  exit 1
fi

rm -rf $WORKDIR
mkdir -p $WORKDIR

echo "=== $(echo "$MOLECULE" | tr '[:lower:]' '[:upper:]') PIPELINE START ==="

# Step 1: xyz -> mol2
echo "Step 1: Open Babel"
obabel $INPUT -O $WORKDIR/${MOLECULE}.raw.mol2

# Step 2: antechamber
echo "Step 2: antechamber"
antechamber -i $WORKDIR/${MOLECULE}.raw.mol2 \
            -fi mol2 \
            -o $WORKDIR/${MOLECULE}.mol2 \
            -fo mol2 \
            -c bcc \
            -at gaff2 \
            -pf y \
            > $WORKDIR/antechamber.log 2>&1

# Step 3: parmchk2
echo "Step 3: parmchk2"
parmchk2 -i $WORKDIR/${MOLECULE}.mol2 \
         -f mol2 \
         -o $WORKDIR/${MOLECULE}.frcmod

# Step 4: tleap
echo "Step 4: tleap"
cat <<EOF > $WORKDIR/tleap.in
source leaprc.gaff2
mol = loadmol2 ${MOLECULE}.mol2
loadamberparams ${MOLECULE}.frcmod
saveamberparm mol ${MOLECULE}.prmtop ${MOLECULE}.inpcrd
quit
EOF

(cd $WORKDIR && tleap -f tleap.in > tleap.log 2>&1)

# Step 5: single point energy
echo "Step 5: sander"
cat <<EOF > $WORKDIR/sp.in
Single point
&cntrl
  imin=1,
  maxcyc=0,
  ntb=0,
  cut=999.0,
/
EOF

(cd $WORKDIR && sander -O -i sp.in -p ${MOLECULE}.prmtop -c ${MOLECULE}.inpcrd -o sp.out -r sp.rst)

# Step 6: energy extraction
echo "Step 6: energy"
grep -A 9 "FINAL RESULTS" $WORKDIR/sp.out || echo "Energy not found (check $WORKDIR/sp.out)"

echo "=== DONE ==="
