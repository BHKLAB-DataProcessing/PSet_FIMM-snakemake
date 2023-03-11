from os import path
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(
    access_key_id=config["key"],
    secret_access_key=config["secret"],
    host=config["host"],
    stay_on_remote=False
)

prefix = config["prefix"]
rna_tool = 'Kallisto-0.46.1'
rna_ref = 'Gencode_v33'
basePath = "https://orcestradata.blob.core.windows.net/gcsi/gCSI/2018"

rna_tool_dir = rna_tool.replace('-', '_')
rnaseq_dir = path.join(prefix, "processed",
                       rna_tool_dir, rna_tool_dir + '_' + rna_ref)
rna_ref_file = rna_ref.replace('_', '.') + '.annotation.RData'

rule get_fimm:
    input:
        prefix + "processed/profiles.RData",
        prefix + "processed/drug.info.rds",
        prefix + "processed/cell.info.rds",
        prefix + "processed/curationCell.rds",
        prefix + "processed/curationDrug.rds",
        prefix + "processed/curationTissue.rds",
        prefix + "processed/sens.info.rds",
        prefix + "processed/sens.prof.rds",
        prefix + "processed/sens.raw.rds",
        prefix + "download/drugs_with_ids.csv",
        prefix + "download/cell_annotation_all.csv",
    output:
        prefix + "FIMM.rds"
    shell:
        """
        Rscript scripts/getFIMM.R {prefix}
        """

rule recalculate_and_assemble:
    input:
        prefix + "processed/raw_sense_slices.zip",
    output:
        prefix + "processed/profiles.RData"
    shell:
        """
        Rscript scripts/recalculateAndAssemble.R {prefix}
        """

rule process_fimm:
    input:
        prefix + "download/drugs_with_ids.csv",
        prefix + "download/cell_annotation_all.csv",
        prefix + "download/nature20171-s1.xls",
        prefix + "download/nature20171-s2.xlsx"
    output:
        prefix + "processed/drug.info.rds",
        prefix + "processed/cell.info.rds",
        prefix + "processed/curationCell.rds",
        prefix + "processed/curationDrug.rds",
        prefix + "processed/curationTissue.rds",
        prefix + "processed/sens.info.rds",
        prefix + "processed/sens.prof.rds",
        prefix + "processed/sens.raw.rds",
        prefix + "processed/raw_sense_slices.zip",
    shell:
        """
        Rscript scripts/processFIMM.R {prefix}
        """

rule download_annotation:
    output:
        prefix + "download/drugs_with_ids.csv",
        prefix + "download/cell_annotation_all.csv",
        prefix + 'download/' + rna_ref_file
    shell:
        """
        wget 'https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/drugs_with_ids.csv' \
            -O {prefix}download/drugs_with_ids.csv
        wget 'https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/cell_annotation_all.csv' \
            -O {prefix}download/cell_annotation_all.csv
        wget 'https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/{rna_ref_file}' \
            -O {prefix}download/{rna_ref_file}
        """

rule download_data:
    output:
        prefix + "download/nature20171-s1.xls",
        prefix + "download/nature20171-s2.xlsx"
    shell:
        """
        wget https://media.nature.com/original/nature-assets/nature/journal/v540/n7631/extref/nature20171-s1.xls \
            -O {prefix}download/nature20171-s1.xls
        wget https://media.nature.com/original/nature-assets/nature/journal/v540/n7631/extref/nature20171-s2.xlsx \
            -O {prefix}download/nature20171-s2.xlsx
        """
