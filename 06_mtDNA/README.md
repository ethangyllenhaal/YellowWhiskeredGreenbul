### Analysis of mitogenome with the program MitoFinder

### _Install MitoFinder_
```{r install MitoFinder, warning=FALSE, message=FALSE}
wget https://github.com/RemiAllio/MitoFinder/archive/master.zip
unzip master.zip
mv MitoFinder-master MitoFinder
cd MitoFinder
./install.sh

PATH/TO/MITOFINDER/mitofinder -h 
```
### _Make a .txt file "mitofinder-cdm.txt" with each sample to run in MitoFinder_
```{r run MitoFinder, warning=FALSE, message=FALSE}
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__COD_Luki__FMNH490025 -1 E_latirostris__COD_Luki__FMNH490025__R1.fastq.gz -2 E_latirostris__COD_Luki__FMNH490025__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__COD_Luki__FMNH490026 -1 E_latirostris__COD_Luki__FMNH490026__R1.fastq.gz -2 E_latirostris__COD_Luki__FMNH490026__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__COD_Tschambi__FMNH473400 -1 E_latirostris__COD_Tschambi__FMNH473400__R1.fastq.gz -2 E_latirostris__COD_Tschambi__FMNH473400__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__COD_Tschambi__FMNH473401 -1 E_latirostris__COD_Tschambi__FMNH473401__R1.fastq.gz -2 E_latirostris__COD_Tschambi__FMNH473401__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__COD_Tschambi__FMNH473402 -1 E_latirostris__COD_Tschambi__FMNH473402__R1.fastq.gz -2 E_latirostris__COD_Tschambi__FMNH473402__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__GAB_Doussala__FMNH389319 -1 E_latirostris__GAB_Doussala__FMNH389319__R1.fastq.gz -2 E_latirostris__GAB_Doussala__FMNH389319__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__GAB_Doussala__FMNH396330 -1 E_latirostris__GAB_Doussala__FMNH396330__R1.fastq.gz -2 E_latirostris__GAB_Doussala__FMNH396330__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__GAB_Doussala__FMNH396331 -1 E_latirostris__GAB_Doussala__FMNH396331__R1.fastq.gz -2 E_latirostris__GAB_Doussala__FMNH396331__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__GHA_AssinFoso__FMNH396575 -1 E_latirostris__GHA_AssinFoso__FMNH396575__R1.fastq.gz -2 E_latirostris__GHA_AssinFoso__FMNH396575__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__GHA_AssinFoso__FMNH396576 -1 E_latirostris__GHA_AssinFoso__FMNH396576__R1.fastq.gz -2 E_latirostris__GHA_AssinFoso__FMNH396576__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__GHA_AssinFoso__FMNH396577 -1 E_latirostris__GHA_AssinFoso__FMNH396577__R1.fastq.gz -2 E_latirostris__GHA_AssinFoso__FMNH396577__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__GNQ_Moka__KU32866 -1 E_latirostris__GNQ_Moka__KU32866__R1.fastq.gz -2 E_latirostris__GNQ_Moka__KU32866__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__GNQ_Moka__KU32908 -1 E_latirostris__GNQ_Moka__KU32908__R1.fastq.gz -2 E_latirostris__GNQ_Moka__KU32908__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__GNQ_Moka__KU32961 -1 E_latirostris__GNQ_Moka__KU32961__R1.fastq.gz -2 E_latirostris__GNQ_Moka__KU32961__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__GNQ_MonteAlen__KU8325 -1 E_latirostris__GNQ_MonteAlen__KU8325__R1.fastq.gz -2 E_latirostris__GNQ_MonteAlen__KU8325__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__GNQ_MonteAlen__KU8326 -1 E_latirostris__GNQ_MonteAlen__KU8326__R1.fastq.gz -2 E_latirostris__GNQ_MonteAlen__KU8326__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__GNQ_MonteAlen__KU8332 -1 E_latirostris__GNQ_MonteAlen__KU8332__R1.fastq.gz -2 E_latirostris__GNQ_MonteAlen__KU8332__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__SLE_KambuiHills__KU19685 -1 E_latirostris__SLE_KambuiHills__KU19685__R1.fastq.gz -2 E_latirostris__SLE_KambuiHills__KU19685__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__SLE_KambuiHills__KU19697 -1 E_latirostris__SLE_KambuiHills__KU19697__R1.fastq.gz -2 E_latirostris__SLE_KambuiHills__KU19697__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__SLE_KambuiHills__KU19776 -1 E_latirostris__SLE_KambuiHills__KU19776__R1.fastq.gz -2 E_latirostris__SLE_KambuiHills__KU19776__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__UGA_AgoroAgu__FMNH489311 -1 E_latirostris__UGA_AgoroAgu__FMNH489311__R1.fastq.gz -2 E_latirostris__UGA_AgoroAgu__FMNH489311__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__UGA_AgoroAgu__FMNH489312 -1 E_latirostris__UGA_AgoroAgu__FMNH489312__R1.fastq.gz -2 E_latirostris__UGA_AgoroAgu__FMNH489312__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__UGA_AgoroAgu__FMNH489313 -1 E_latirostris__UGA_AgoroAgu__FMNH489313__R1.fastq.gz -2 E_latirostris__UGA_AgoroAgu__FMNH489313__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__UGA_Echuya__FMNH384894 -1 E_latirostris__UGA_Echuya__FMNH384894__R1.fastq.gz -2 E_latirostris__UGA_Echuya__FMNH384894__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__UGA_Echuya__FMNH384895 -1 E_latirostris__UGA_Echuya__FMNH384895__R1.fastq.gz -2 E_latirostris__UGA_Echuya__FMNH384895__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__UGA_Echuya__FMNH384896 -1 E_latirostris__UGA_Echuya__FMNH384896__R1.fastq.gz -2 E_latirostris__UGA_Echuya__FMNH384896__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__UGA_Nteko__FMNH384909 -1 E_latirostris__UGA_Nteko__FMNH384909__R1.fastq.gz -2 E_latirostris__UGA_Nteko__FMNH384909__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__UGA_Nteko__FMNH384910 -1 E_latirostris__UGA_Nteko__FMNH384910__R1.fastq.gz -2 E_latirostris__UGA_Nteko__FMNH384910__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__UGA_Nteko__FMNH384911 -1 E_latirostris__UGA_Nteko__FMNH384911__R1.fastq.gz -2 E_latirostris__UGA_Nteko__FMNH384911__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__UGA_Rwenzori__FMNH355307 -1 E_latirostris__UGA_Rwenzori__FMNH355307__R1.fastq.gz -2 E_latirostris__UGA_Rwenzori__FMNH355307__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__UGA_Rwenzori__FMNH355308 -1 E_latirostris__UGA_Rwenzori__FMNH355308__R1.fastq.gz -2 E_latirostris__UGA_Rwenzori__FMNH355308__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__UGA_Rwenzori__FMNH355309 -1 E_latirostris__UGA_Rwenzori__FMNH355309__R1.fastq.gz -2 E_latirostris__UGA_Rwenzori__FMNH355309__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__UGA_Rwenzori__FMNH355310 -1 E_latirostris__UGA_Rwenzori__FMNH355310__R1.fastq.gz -2 E_latirostris__UGA_Rwenzori__FMNH355310__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_latirostris__UGA_Rwenzori__FMNH355311 -1 E_latirostris__UGA_Rwenzori__FMNH355311__R1.fastq.gz -2 E_latirostris__UGA_Rwenzori__FMNH355311__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_virens__CMR_Nzimdipnenkum__FMNH511205 -1 E_virens__CMR_Nzimdipnenkum__FMNH511205__R1.fastq.gz -2 E_virens__CMR_Nzimdipnenkum__FMNH511205__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
singularity run /symbiont/apps/mitofinder_v1.4.2.sif -j E_virens__COD_Uma__FMNH501571 -1 E_virens__COD_Uma__FMNH501571__R1.fastq.gz -2 E_virens__COD_Uma__FMNH501571__R2.fastq.gz -r E_virens_MG762099.1.gb -o 2 -p 20 -m 10
```

### _Run MitoFinder_
```{r run MitoFinder, warning=FALSE, message=FALSE}
nohup  parallel -j10 < mitofinder-cdm.txt
```
