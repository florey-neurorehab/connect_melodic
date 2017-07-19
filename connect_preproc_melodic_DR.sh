#Takes NON-filtered cleaned, smoothed output from CoNNECT preproc script,
#filters, runs melodic to 70 dims. 

hdir=/home/peter/Desktop/Connect/rest/output/
holddir=/home/peter/Desktop/melodic/connect/subj_filtered_raw/
outdir=/home/peter/Desktop/melodic/connect/

mask=$FSL_DIR'/data/standard/MNI152_T1_3mm_brain_mask.nii.gz'

mkdir $holddir

###Select only baseline files (note - the -v in grep inverts the pattern search), writes to file###
subj_files=$(find $hdir/*/preproc_epis -name 'sresiduals_global_trans.nii' | grep -v P2)
printf "$subj_files" > $outdir'baseline_subs_test.txt' #Output to text file for troubleshooting.

#Get subj info, loop:
for subj in $subj_files;
    do subj_path=`dirname $subj`
    subj_id=${subj_path/$hdir'/'/''} #Remove $hdir path retaining subj id
    subj_id=${subj_id%/*}  #Retain subj_id, dump preproc path.
    
    printf $subj_id

###FILTER###
    #Extract temporal mean (Tmean) image.
    fslmaths $subj -Tmean $holddir'/tempmean.nii.gz'#

    #Filter (16.66...7 = hp in volumes, -1 = no LP), re-add mean
    fslmaths $subj -bptf 16.666666667 -1 -add $holddir'/tempmean.nii.gz' $holddir/$subj_id'_epi_hpfilt.nii'
    
    #Remove temporal mean file.

done




###Exclude participants with bad data (see .../Connect/connect_bl_qc.xlsx for list)###
declare -a ex_list=("D_H06" "D_H16" "D_H20" "D_H24" "D_S21" "D_S29" "D_S30" "D_S05" "D_S11" "D_S17" "D_S33" "D_S37")

subj_files=$(find $holddir -name *hpfilt* | sort)
printf "$subj_files" > $holddir/subj_files.txt

culled_subj_files=$subj_files

for i in "${ex_list[@]}";
    do printf Removing $i
    culled_subj_files=$(printf "$culled_subj_files" | grep -v $i)
done

#Output culled list - NOTE: Output is NOT ORDERED by group!
printf "$culled_subj_files" > $holddir/subj_files_culled.txt





###Create variables and output with correct members of each group:###
culled_subj_files=$(cat $holddir/subj_files_culled.txt)

#Control
declare -a con_arr=('H02' 'H07' 'H08' 'H09' 'H12' 'H14' 'H15' 'H17' 'H19' 'H25' 'H26' 'H27' 'H28')

first='H01'

#Initialise variable
con=""

subj=$(printf "$culled_subj_files" | grep $first)
con=$subj

for c in "${con_arr[@]}";
do
    subj=$(printf "$culled_subj_files" | grep $c)
    con="${con}\n${subj}"
done

#Left
declare -a left_arr=('S10' 'S12' 'S14' 'S15' 'S19' 'S22' 'S24' 'S27' 'S31' 'S32' 'S35' 'S36' 'S38')

first='S04'

left=""
subj=$(printf "$culled_subj_files" | grep $first)
left=$subj

for l in "${left_arr[@]}";
    do subj=$(printf "$culled_subj_files" | grep $l)
    left="${left}\n${subj}" 
done


#Right
declare -a right_arr=('S02' 'S03' 'S06' 'S08' 'S09' 'S13' 'S16' 'S18' 'S20' 'S23' 'S25' 'S28' 'S34')

first='S01'

right=""
subj=$(printf "$culled_subj_files" | grep $first)
right=$subj

for r in "${right_arr[@]}";
    do subj=$(printf "$culled_subj_files" | grep $r)
    right="${right}\n${subj}" 
done


subj_files_culled_ordered="$con\n$left\n$right"

printf "$subj_files_culled_ordered" > $holddir/subj_files_culled_ordered.txt



###Run melodic using a tensor ICA approach (instead of a space x time decomposition, data is decomposed into
#space x time x subject matricies that allow for modellling of between subject variance.
#NOTE: CONCAT IS USED BECAUSE SIMILAR SPATIAL COMPONENTS ARE EXPECTED BUT SIMILAR TIME SERIES (TICA) ARE NOT. ###
melodic -i $holddir/subj_files_culled_ordered.txt -o $outdir'baseline' -v --nobet --bgthreshold=1 --tr=3.000 --report --mmthresh=0.5 --Ostats --approach=concat --dim=20

#Run FSLview to select components for removal
fslview -m ortho $FSLDIR/data/standard/MNI152_T1_3mm_brain.nii.gz $outdir'baseline/melodic_IC.nii.gz' -l 'Hot' -b 1.9,10 &
 
#Remove artifactual data from ICA
mkdir $outdir'baseline/split/'
fslsplit $outdir'baseline/melodic_IC.nii.gz' $outdir'baseline/split/'

printf '\nPlease enter artifactual components (starting from 0), with a space between each component number.\n'

read rem_comps

for x in $rem_comps; 
    do printf Removing component $x 
    if [ "$x" -gt "9" ] 
    then 
        rm -f $outdir'baseline/split/00'$x'.nii.gz'
    else
        rm -f $outdir'baseline/split/000'$x'.nii.gz'
    fi
done

comp_files=$(ls $outdir'baseline/split/'*.nii.gz)

fslmerge -a $outdir'baseline/clean_ica.nii' $comp_files

rm -f -R $outdir'baseline/split/'


###Running dual regression###

printf '\n***RUNNING DUAL REGRESSION***\n'
infiles=$(cat $holddir/subj_files_culled_ordered.txt)

#Run dual_regression to get component ts / images (BUT NOT RANDOMISE - comes later, explicit control of mask).

mkdir -p $outdir'DR/baseline-test'
dual_regression_thr $outdir'baseline/melodic_IC.nii.gz' 1 $outdir'design_all_culled/design_all_culled.mat' $outdir'design_all_culled/design_all_culled.con' 0 $outdir'DR/baseline-test' $infiles #Used to run Stage 2 dual regression (DR without randomise).

#Run randomise (more control over masking)
dr_infiles=$(find $outdir'DR' -name '*stage2*Z.nii.gz' | sort) #Keep dual regression in subject order.
fslmerge -a $outdir'baseline/dr_files.nii.gz' $dr_infiles

randomise -i $outdir'baseline/dr_files.nii.gz' -o $outdir'randomise/' -m $mask -d $outdir'design_all_culled/design_all_culled.mat' -t $outdir'design_all_culled/design_all_culled.con' -T -n 5000



