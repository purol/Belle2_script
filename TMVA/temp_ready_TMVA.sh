mkdir -p ./root_file_for_TMVA

cp ./CHG_analysis/test_v001/ver0_signal/final_output/Mxs_large/final_output_merge_Mxs_larger.root ./root_file_for_TMVA/final_output_merge_Mxs_larger_CHG_test.root
cp ./CHG_analysis/test_v001/ver0_signal/final_output/Mxs_small/final_output_merge_Mxs_smaller.root ./root_file_for_TMVA/final_output_merge_Mxs_smaller_CHG_test.root
cp ./CHG_analysis/train_v001/ver0_signal/final_output/Mxs_large/final_output_merge_Mxs_larger.root ./root_file_for_TMVA/final_output_merge_Mxs_larger_CHG_train.root
cp ./CHG_analysis/train_v001/ver0_signal/final_output/Mxs_small/final_output_merge_Mxs_smaller.root ./root_file_for_TMVA/final_output_merge_Mxs_smaller_CHG_train.root

cp ./MIX_analysis/test_v001/ver0_signal/final_output/Mxs_large/final_output_merge_Mxs_larger.root ./root_file_for_TMVA/final_output_merge_Mxs_larger_MIX_test.root
cp ./MIX_analysis/test_v001/ver0_signal/final_output/Mxs_small/final_output_merge_Mxs_smaller.root ./root_file_for_TMVA/final_output_merge_Mxs_smaller_MIX_test.root
cp ./MIX_analysis/train_v001/ver0_signal/final_output/Mxs_large/final_output_merge_Mxs_larger.root ./root_file_for_TMVA/final_output_merge_Mxs_larger_MIX_train.root
cp ./MIX_analysis/train_v001/ver0_signal/final_output/Mxs_small/final_output_merge_Mxs_smaller.root ./root_file_for_TMVA/final_output_merge_Mxs_smaller_MIX_train.root

cp ./UUBAR_analysis/test_v001/ver0_signal/final_output/Mxs_large/final_output_merge_Mxs_larger.root ./root_file_for_TMVA/final_output_merge_Mxs_larger_UUBAR_test.root
cp ./UUBAR_analysis/test_v001/ver0_signal/final_output/Mxs_small/final_output_merge_Mxs_smaller.root ./root_file_for_TMVA/final_output_merge_Mxs_smaller_UUBAR_test.root
cp ./UUBAR_analysis/train_v001/ver0_signal/final_output/Mxs_large/final_output_merge_Mxs_larger.root ./root_file_for_TMVA/final_output_merge_Mxs_larger_UUBAR_train.root
cp ./UUBAR_analysis/train_v001/ver0_signal/final_output/Mxs_small/final_output_merge_Mxs_smaller.root ./root_file_for_TMVA/final_output_merge_Mxs_smaller_UUBAR_train.root

cp ./DDBAR_analysis/test_v001/ver0_signal/final_output/Mxs_large/final_output_merge_Mxs_larger.root ./root_file_for_TMVA/final_output_merge_Mxs_larger_DDBAR_test.root
cp ./DDBAR_analysis/test_v001/ver0_signal/final_output/Mxs_small/final_output_merge_Mxs_smaller.root ./root_file_for_TMVA/final_output_merge_Mxs_smaller_DDBAR_test.root
cp ./DDBAR_analysis/train_v001/ver0_signal/final_output/Mxs_large/final_output_merge_Mxs_larger.root ./root_file_for_TMVA/final_output_merge_Mxs_larger_DDBAR_train.root
cp ./DDBAR_analysis/train_v001/ver0_signal/final_output/Mxs_small/final_output_merge_Mxs_smaller.root ./root_file_for_TMVA/final_output_merge_Mxs_smaller_DDBAR_train.root

cp ./SSBAR_analysis/test_v001/ver0_signal/final_output/Mxs_large/final_output_merge_Mxs_larger.root ./root_file_for_TMVA/final_output_merge_Mxs_larger_SSBAR_test.root
cp ./SSBAR_analysis/test_v001/ver0_signal/final_output/Mxs_small/final_output_merge_Mxs_smaller.root ./root_file_for_TMVA/final_output_merge_Mxs_smaller_SSBAR_test.root
cp ./SSBAR_analysis/train_v001/ver0_signal/final_output/Mxs_large/final_output_merge_Mxs_larger.root ./root_file_for_TMVA/final_output_merge_Mxs_larger_SSBAR_train.root
cp ./SSBAR_analysis/train_v001/ver0_signal/final_output/Mxs_small/final_output_merge_Mxs_smaller.root ./root_file_for_TMVA/final_output_merge_Mxs_smaller_SSBAR_train.root

cp ./CHARM_analysis/test_v001/ver0_signal/final_output/Mxs_large/final_output_merge_Mxs_larger.root ./root_file_for_TMVA/final_output_merge_Mxs_larger_CHARM_test.root
cp ./CHARM_analysis/test_v001/ver0_signal/final_output/Mxs_small/final_output_merge_Mxs_smaller.root ./root_file_for_TMVA/final_output_merge_Mxs_smaller_CHARM_test.root
cp ./CHARM_analysis/train_v001/ver0_signal/final_output/Mxs_large/final_output_merge_Mxs_larger.root ./root_file_for_TMVA/final_output_merge_Mxs_larger_CHARM_train.root
cp ./CHARM_analysis/train_v001/ver0_signal/final_output/Mxs_small/final_output_merge_Mxs_smaller.root ./root_file_for_TMVA/final_output_merge_Mxs_smaller_CHARM_train.root

cd ./root_file_for_TMVA
cp /media/sf_virtualbox_folder/20211109/root_file_for_TMVA/ReadRootFiles_r_sp.c ./
.q | root ./ReadRootFiles_r_sp.c
cd -

mkdir -p ./For_TMVA_BKG/test/Mxs_large
mkdir -p ./For_TMVA_BKG/test/Mxs_small
mkdir -p ./For_TMVA_BKG/train/Mxs_large
mkdir -p ./For_TMVA_BKG/train/Mxs_small

mv ./root_file_for_TMVA/*larger*test_data.root ./For_TMVA_BKG/test/Mxs_large
mv ./root_file_for_TMVA/*smaller*test_data.root ./For_TMVA_BKG/test/Mxs_small
mv ./root_file_for_TMVA/*larger*train_data.root ./For_TMVA_BKG/train/Mxs_large
mv ./root_file_for_TMVA/*smaller*train_data.root ./For_TMVA_BKG/train/Mxs_small
