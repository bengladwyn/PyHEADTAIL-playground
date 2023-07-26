###############################################################################
# Script to launch IW2D simulations
# Also checks the results of the test cases against results obtained on LXPLUS at CERN
###############################################################################

set -e

#######
# Paths
#######

SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IW2D_PATH="$(cd "$SCRIPT_PATH/.." && pwd)"
BIN_PATH="$IW2D_PATH/IW2D/cpp"
TEST_FILES_PATH="$IW2D_PATH/examples/outputs"
INPUT_FILES_PATH="$IW2D_PATH/examples/input_files"

####################
# Round chamber case
####################

printf "\n\n**********************\nTesting round chamber\n"

cd $TEST_FILES_PATH

$BIN_PATH/roundchamber.x < $INPUT_FILES_PATH/RoundChamberInputFile.txt > out.txt
ROOTNAME="WLHC_1layers2.00mm_test_coll_round.dat"

for component in 'long' 'xdip' 'ydip' 'xquad' 'yquad'; do
    cmp --silent Z${component}${ROOTNAME} $TEST_FILES_PATH/round_chamber/Z${component}${ROOTNAME} || echo "Round chamber: Z${component} files are different"
done

tail -n 1 out.txt
rm -f out.txt *${ROOTNAME}

printf "Done\n**********************\n"

##############################
# Round chamber case with SiO2
##############################

printf "\n\n**********************\nTesting round chamber - SiO2\n"

cd $TEST_FILES_PATH

$BIN_PATH/roundchamber.x < $INPUT_FILES_PATH/RoundChamberSiO2InputFile.txt > out.txt
ROOTNAME="WLHC_2layers1.00mm_test_coll_round_sio2.dat"

for component in 'long' 'xdip' 'ydip' 'xquad' 'yquad'; do
    cmp --silent Z${component}${ROOTNAME} $TEST_FILES_PATH/round_chamber_sio2/Z${component}${ROOTNAME} || echo "Round chamber SiO2: Z${component} files are different"
done

tail -n 1 out.txt
rm -f out.txt *${ROOTNAME}

printf "Done\n**********************\n"

###################
# Flat chamber case
###################

printf "\n\n**********************\nTesting flat chamber\n"

cd $TEST_FILES_PATH

$BIN_PATH/flatchamber.x < $INPUT_FILES_PATH/FlatChamberInputFile.txt > out.txt
ROOTNAME="WLHC_1layersup_0layersdown4.00mm_test_TCDQ_Cu.dat"

for component in 'long' 'xdip' 'ydip' 'xquad' 'yquad' 'ycst'; do
    cmp --silent Z${component}${ROOTNAME} $TEST_FILES_PATH/flat_chamber/Z${component}${ROOTNAME} || echo "Flat chamber: Z${component} files are different"
done

tail -n 1 out.txt
rm -f out.txt *${ROOTNAME}

printf "Done\n**********************\n"

###############################
# Flat chamber, non-linear case
###############################

printf "\n\n**********************\nTesting flat chamber - non-linear case\n"

cd $TEST_FILES_PATH

$BIN_PATH/flatchamber.x < $INPUT_FILES_PATH/FlatChamberInputFileNonlinear.txt > out.txt
ROOTNAME="WLHC_1layersup_0layersdown4.00mm_test_TCDQ_nonlinear.dat"

for component in 'l0001' 'l0002' 'l0020' 'l0100' 'l0101' 'l0200' 'l1010' 'l2000' 'long' 'x0011' 'x0110' 'x1001' 'x1100' 'xdip' 'xquad' 'y0002' 'y0020' 'y0101' 'y0200' 'y1010' 'y2000' 'ycst' 'ydip' 'yquad'; do
	cmp --silent Z${component}${ROOTNAME} $TEST_FILES_PATH/nonlinear_flat_chamber/Z${component}${ROOTNAME} || echo "Flat chamber (non-linear): Z${component} files are different"
done

tail -n 1 out.txt
rm -f out.txt *${ROOTNAME}

printf "Done\n**********************\n"

#########################
# Wake round chamber case
#########################

printf "\n\n**********************\nTesting wake round chamber\n"

cd $TEST_FILES_PATH

$BIN_PATH/wake_roundchamber.x < $INPUT_FILES_PATH/WakeRoundChamberInputFile.txt > out.txt
ROOTNAME="WLHC_1layers4.00mm_test_long_wake_in_front_ss_gamma_1e4"

for component in 'long' 'xdip' 'ydip' 'xquad' 'yquad'; do
    cmp --silent W${component}${ROOTNAME}_precise.dat $TEST_FILES_PATH/round_chamber_wake/W${component}${ROOTNAME}_precise.dat\
     || (echo "Wake round chamber: W${component} files are different"; \
     pr -m -t -s\  W${component}${ROOTNAME}_precise.dat $TEST_FILES_PATH/round_chamber_wake/W${component}${ROOTNAME}_precise.dat\
     | awk -v "sum=0" -v "comp=${component}" 'function abs(x){return ((x < 0.0) ? -x : x)}\
     {if ( (NR>1) && ($1!=0) ) {max=((abs($2-$4)>max)?abs($2-$4):max)}; if ( (NR>1) && ($1!=0) && ( (abs($2-$4)>(1e7*((comp=="long")?100:1))) ) ) {sum+=1} }\
     END {print sum " out of " NR-2 " values are off by more than 1e" 7+((comp=="long")?2:0) " (max diff.=" max ")"}')
done

tail -n 1 out.txt
rm -f out.txt *${ROOTNAME}*.dat

printf "Done\n**********************\n"

########################
# Wake flat chamber case
########################

printf "\n\n**********************\nTesting wake flat chamber\n"

cd $TEST_FILES_PATH

$BIN_PATH/wake_flatchamber.x < $INPUT_FILES_PATH/WakeFlatChamberInputFile.txt > out.txt
ROOTNAME="WLHC_2layersup_0layersdown6.50mm"

for component in 'long' 'xdip' 'ydip' 'xquad' 'yquad' 'ycst'; do
    cmp --silent W${component}${ROOTNAME}_precise.dat $TEST_FILES_PATH/flat_chamber_wake/W${component}${ROOTNAME}_precise.dat\
     || (echo "Wake flat chamber: W${component} files are different" ; \
     pr -m -t -s\  W${component}${ROOTNAME}_precise.dat $TEST_FILES_PATH/flat_chamber_wake/W${component}${ROOTNAME}_precise.dat\
     | awk -v "sum=0" -v "comp=${component}" 'function abs(x){return ((x < 0.0) ? -x : x)}\
     {if ( (NR>1) && ($1!=0) ) {max=((abs($2-$4)>max)?abs($2-$4):max)}; if ( (NR>1) && ($1!=0) && ( (abs($2-$4)>(1e7*((comp=="long")?100:1))) ) ) {sum+=1} }\
     END {print sum " out of " NR-2 " values are off by more than 1e" 7+((comp=="long")?2:0) " (max diff.=" max ")"}')
done

tail -n 1 out.txt
rm -f out.txt *${ROOTNAME}*.dat

printf "Done\n**********************\n"
