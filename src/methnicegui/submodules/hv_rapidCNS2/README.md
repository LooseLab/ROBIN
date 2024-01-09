### Nextflow wf-human-variation and rapidCNS2 combined workflow

It takes an aligned bam of ONT data along with the .bai index and runs epi2me-labs wf-human-variation to generate SNP, SV, and CNV variants and also aggregate methylation data (the --snp, --sv, --cnv, and --methyl options).

Outputs from wf-human-variation are then analysed using Rapid-CNS2 to generate the final report. They are also optionally analysed using the Sturgeon classifier (see optional extra parameters below).

Input BAM file **MUST** contain methylation data (MM:Z tags). The workflow should ERROR if MM:Z tags are not present.
The .bai index must be in the same location as the input bam.

The provided reference should be the same reference used to align reads in the input BAM.

At present, it relies upon the *outdir* being within the current working directory. Setting a full path to some other location will cause errors.

### Latest updates in v.0.8.1:
- wf-human-variation "--mod" module uses the BAM file with the merged mods
- wf-human-variation "--snp" uses the raw BAM subsetted to just the target regions
- wf-human-variation "--sv" use the subsetted BAM, which has all associated supplementary alignments retained
- wf-human-variation "--cnv" uses the raw input BAM
- two versions of the output report are created; a lite v with a simple mutations table and a full v with interactive igvreport
- hv_rapidCNS2 workflow version number now included in the final report
- sequencer information derived from input BAM and included in report (can be specified with --seq)
### Update in v0.8.2:
- MGMT coverage, methylation status, and methylation plot now correctly displayed in lite and full reports

### Requirements:
```
docker
nextflow
```

### Pull latest docker image:
```
docker pull graefox/rapid_cns2:latest
```

### pull/clone latest version of workflow, and required dependencies
```
git clone https://github.com/graemefox/hv_rapidCNS2.git
nextflow pull epi2me-labs/wf-human-variation
```


### Example command:
```
## define sample name, ID and output directory, input BAM and reference genome:

SAMPLE=sample_01
PATIENT=JohnDoe
OUTPUT_DIR=${SAMPLE}_output
BAM=my_data.bam
REFERENCE=my_reference.fa.gz
ANNOTATIONS=my_annotation_set.gtf

## run the pipeline

nextflow run hv_rapidCNS2/main.nf \
-with-docker graefox/rapid_cns2:latest \
-with-report ${OUTPUT_DIR}/${SAMPLE}_nextflow_report.html \
--sample ${SAMPLE} \
--patient ${PATIENT} \
--bam ${BAM} \
--outdir ${OUTPUT_DIR} \
--reference ${REFERENCE} \
--annotations ${ANNOTATIONS}

```

### Optional extra parameters (with their default values)
These a have default values specified in the nextflow.config file, but you may override them on the CLI.
```
--threads 16 (CPUs to use [default: 64]) 
--bam_min_coverage (minimum coverage required to run the epi2melabs/wf-human-variation stages [ default: 5]) 
--minimum_mgmt_cov (minimum avg coverage at the mgmt promoter. Coverage must be greater than this to run the analysis of mgmt methylation)
--sturgeon (the nextflow will ALSO run the sturgeon (https://github.com/marcpaga/sturgeon) classifier if the --sturgeon flag is passed [Defualt behaviour is to NOT run sturgeon])
--nanoplot (nextflow will ALSO run NanoPlot to generate a QC report[ Default behaviour is to NOT run nanoplot])

```

### To run with slurm
Add `-process.executor='slurm'` to your nextflow command, then run as normal. You do not need to submit a script with SBATCH, just run the nextflow command as normal and nextflow knows
to submit each process into SLURM.

### Troubleshooting tips
If the run seems to hang forever at the cnvpytor step, it may be that you have not indexed your input bam. This is also just quite a long process.

If you get the Docker Error: "docker: permission denied while trying to connect to the docker daemon socket".... on Ubuntu (based) systems, you need to add your user to the docker group. 
Follow the instructions here: (https://www.digitalocean.com/community/questions/how-to-fix-docker-got-permission-denied-while-trying-to-connect-to-the-docker-daemon-socket)
