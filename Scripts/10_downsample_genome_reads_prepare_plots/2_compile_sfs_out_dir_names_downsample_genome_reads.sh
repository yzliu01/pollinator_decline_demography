
## all: checked methods
input_sfs_dir=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/data/bee_proj_data/SFS_data
cd $input_sfs_dir

## downsample reads
#out_sfs_dir=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/sofwtare/stairway_plot_v2/stairway_plot_v2.1.2/stairway_plot_blueprint/pool_downsample_reads/templates

## downsample genome
out_sfs_dir=/home/yzliu/eDNA/faststorage/yzliu/DK_proj/sofwtare/stairway_plot_v2/stairway_plot_v2.1.2/stairway_plot_blueprint/pool_shuf_downsample_genome/templates

for species in {Hae,Mar,Pas,Vet};do ## downsample genome

#for species in {Hae,Mar};do ## downsample reads
    ## different sort orders
    ## downsample genome
    sfs1=`find -maxdepth 1 -print | grep "$species" | egrep "P_[0-1][0-9]" | grep "equal_self.sfs$" | sort -V`
  
    ## downsample reads
    #sfs1=`find -maxdepth 1 -print | grep "$species" | egrep '_[0-9]P|_[0-1][0-9]P' | grep "equal_self.sfs$" | perl -pe 'm/_(\d+)P_downsample_reads\.(\d+)x/; ($p, $xn) = ($1, $2); $_ = sprintf("%05d_%05d_%s", $xn, $p, $_);' | sort -t '_' -k1,1n -k2,2n | cut -d'_' -f3-`
    sfs2=`find -maxdepth 1 -print | grep "xxxxx_old_script" | sort -V`
# species: Hae
    if [ "$species" = "Hae" ];then
        for file_name in $sfs1;do
            #    echo -e "\n"
            new_file_name=${file_name//}
            awk -v var=$file_name 'BEGIN{FS=OFS=" "}NR>=1,NR<=1{print $0 var OFS var}' $file_name | cut -d ' ' -f 1-41 >> $out_sfs_dir/"$species"_39_dipS.1_3_5_7_10x.P_01_02_03_04_05_06_07_08_09_10.sfs;

        done
        
        # no need to change
        for file_name in $sfs2;do
            awk -v var=$file_name 'BEGIN{FS=OFS=" "}NR>=1,NR<=1{print $0 var OFS var}' $file_name | cut -d ' ' -f 1-41 >> $out_sfs_dir/"$species"_39_dipS.1_3_5_7_10x.P_01_02_03_04_05_06_07_08_09_10.sfs;
        done

# species: Mar
    elif [ "$species" = "Mar" ];then
        for file_name in $sfs1;do
            # empty sfs data for 10x
            awk -v var=$file_name 'BEGIN{FS=OFS=" "}NR>=1,NR<=1{print $0 var OFS var}' $file_name | cut -d ' ' -f 1-42 >> $out_sfs_dir/"$species"_40_dipS.1_3_5_7_10x.P_01_02_03_04_05_06_07_08_09_10.sfs;
        done
        # no need to change
        for file_name in $sfs2;do
            awk -v var=$file_name 'BEGIN{FS=OFS=" "}NR>=1,NR<=1{print $0 var OFS var}' $file_name | cut -d ' ' -f 1-42 >> $out_sfs_dir/"$species"_40_dipS.1_3_5_7_10x.P_01_02_03_04_05_06_07_08_09_10.sfs;
        done

# species: Pas
    elif [ "$species" = "Pas" ];then
        for file_name in $sfs1;do
            awk -v var=$file_name 'BEGIN{FS=OFS=" "}NR>=1,NR<=1{print $0 var OFS var}' $file_name | cut -d ' ' -f 1-36 >> $out_sfs_dir/"$species"_34_dipS.1_3_5_7_10x.P_01_02_03_04_05_06_07_08_09_10.sfs;
        done
        # no need to change
        for file_name in $sfs2;do
            awk -v var=$file_name 'BEGIN{FS=OFS=" "}NR>=1,NR<=1{print $0 var OFS var}' $file_name | cut -d ' ' -f 1-36 >> $out_sfs_dir/"$species"_34_dipS.1_3_5_7_10x.P_01_02_03_04_05_06_07_08_09_10.sfs;
        done
    else

# species: Vet
        for file_name in $sfs1;do
            awk -v var=$file_name 'BEGIN{FS=OFS=" "}NR>=1,NR<=1{print $0 var OFS var}' $file_name | cut -d ' ' -f 1-31 >> $out_sfs_dir/"$species"_29_dipS.1_3_5_7_10x.P_01_02_03_04_05_06_07_08_09_10.sfs;

        done
        # no need to change
        for file_name in $sfs2;do
            awk -v var=$file_name 'BEGIN{FS=OFS=" "}NR>=1,NR<=1{print $0 var OFS var}' $file_name | cut -d ' ' -f 1-31 >> $out_sfs_dir/"$species"_29_dipS.1_3_5_7_10x.P_01_02_03_04_05_06_07_08_09_10.sfs;
        done

    fi
done
