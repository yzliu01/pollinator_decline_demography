
cd /home/yzliu/eDNA/faststorage/yzliu/DK_proj/sofwtare/stairway_plot_v2/stairway_plot_v2.1.2/stairway_plot_blueprint/pool_shuf_downsample_genome/templates

# Hae_39_dipS.1_3_5_7_10x.P_01_02_03_04_05_06_07_08_09_10.sfs
for species in {Hae_39,Mar_40,Pas_34,Vet_29};do
    if [ "$species" = "Hae_39" ];then
        cat "$species"_dipS.1_3_5_7_10x.P_01_02_03_04_05_06_07_08_09_10.sfs |\
        # cut: the delimiter must be a single character
        #cut -d "/" -f 1,3,5 | \
        awk -F "./" '{print $1 $2 $3}' |\
        sed -e 's/concated\.//g' -e 's/_New//g' -e 's/100kb_g1500x_regions.GQ_issue_solved.SNP_softmask_genic_bi_FMT_//g' -e 's/_noMS.shuf//g' -e 's/.equal_self.sfs//g' -e 's/_1500x//g' > ./"$species"_dipS.1_3_5_7_10x.P_01_02_03_04_05_06_07_08_09_10.sfs.blueprint.txt

    elif [ "$species" = "Mar_40" ];then
        cat "$species"_dipS.1_3_5_7_10x.P_01_02_03_04_05_06_07_08_09_10.sfs |\
        awk -F "./" '{print $1 $2 $3}' |\
        sed -e 's/concated\.//g' -e 's/_New//g' -e 's/100kb_g1500x_regions.all_chr.sorted.GQ_issue_solved.SNP_softmask_genic_bi_FMT_//g' -e 's/_noMS.shuf//g' -e 's/.equal_self.sfs//g' -e 's/_1500x//g' > ./"$species"_dipS.1_3_5_7_10x.P_01_02_03_04_05_06_07_08_09_10.sfs.blueprint.txt

   
    elif [ "$species" = "Pas_34" ];then
        cat "$species"_dipS.1_3_5_7_10x.P_01_02_03_04_05_06_07_08_09_10.sfs |\
        awk -F "./" '{print $1 $2 $3}' |\
        sed -e 's/concated\.//g' -e 's/_New//g' -e 's/100kb_g1500x_regions.all_chr.sorted.GQ_issue_solved.SNP_softmask_genic_bi_FMT_//g' -e 's/_noMS.shuf//g' -e 's/.equal_self.sfs//g' -e 's/_1500x//g' > ./"$species"_dipS.1_3_5_7_10x.P_01_02_03_04_05_06_07_08_09_10.sfs.blueprint.txt

    else
        cat "$species"_dipS.1_3_5_7_10x.P_01_02_03_04_05_06_07_08_09_10.sfs |\
        awk -F "./" '{print $1 $2 $3}' |\
        sed -e 's/concated\.//g' -e 's/_New//g' -e 's/100kb_g1500x_regions.SNP_softmask_genic_bi_FMT_//g' -e 's/_noMS.shuf//g' -e 's/.equal_self.sfs//g' -e 's/_1500x//g' > ./"$species"_dipS.1_3_5_7_10x.P_01_02_03_04_05_06_07_08_09_10.sfs.blueprint.txt

    fi
done

