#!/usr/bin/env nextflow




Channel.fromPath(params.dataset_list)
  .flatMap{ it.readLines() }
  .view()
  .set { DATALIST }

// fetch data
process get_datasets {
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 2.GB + 10.GB * (task.attempt - 1) }
 	tag "$datasetname"

    	input:
    	val datasetname from DATALIST
    	output:
    	set val(datasetname), file('*.Rds') into SCE_2_H5AD, RUN_SEU, RUN_COM, RUN_LIM, RUN_MNN, RUN_HAR
	set val(datasetname), val('raw'), file("*.Rds") into RAW_ENT    	
    
		
    	shell:
    	'''
    	Rfile="!{params.data_dir}/!{datasetname}.Rds"
    	
    	if [[ ! -e $Rfile ]]; then
    	  echo "Please check existence of $Rfile"
    	  false
    	fi
    	ln -s $Rfile .
    	'''
}

	 
process sce_2_h5ad{
        publishDir "./data"
        errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
        memory = { 10.GB + 20.GB * (task.attempt - 1) }
        tag "sce_2_h5ad $datasetname"

        input:
        set val(datasetname), file(datain) from SCE_2_H5AD


        output:
        set val(datasetname), file('*.h5ad') into RUN_BBK, RUN_SCA
	set val(datasetname), val('raw'), file("*.h5ad") into RAW_ASW

        """
        sce2h5ad.R\
                --input_object ${datain}\
                --output ${datasetname}.h5ad
        """
}

process Seurat_3{
 	publishDir "./data/${datasetname}/Corrected_objects" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "Seurat 3 $datasetname"
    	
	input:
    	set val(datasetname), file(datain) from RUN_SEU
    	
	output:
    	set val(datasetname), val('Seurat3'), file('Seurat3.*.rds') into SEU_ENT, SEU_2_H5, SEU_MARK
	
    	"""
    	seurat_integrate.R\
		--input_object ${datain}\
		--batch_key ${params.batch_key}\
		--output_object Seurat3.${datasetname}.rds
    	"""
}

process Combat{
        publishDir "./data/${datasetname}/Corrected_objects" 
        errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
        memory = { 10.GB + 20.GB * (task.attempt - 1) }
        tag "Combat $datasetname"

        input:
        set val(datasetname), file(datain) from RUN_COM
          
        output:
        set val(datasetname), val('Combat'), file('Combat.*.rds') into COM_ENT, COM_2_H5, COM_MARK

        """
        combat.R\
                --input_object ${datain}\
                --batch_key ${params.batch_key}\
                --output_object Combat.${datasetname}.rds
        """
}

process limma{
        publishDir "./data/${datasetname}/Corrected_objects"
        errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
        memory = { 10.GB + 20.GB * (task.attempt - 1) }
        tag "limma $datasetname"

        input:
        set val(datasetname), file(datain) from RUN_LIM

        output:
        set val(datasetname), val('limma'), file('limma.*.rds') into LIM_ENT, LIM_2_H5, LIM_MARK

        """
        limma.R\
                --input_object ${datain}\
                --batch_key ${params.batch_key}\
                --output_object limma.${datasetname}.rds
        """
}

process MNNCorrect{
        publishDir "./data/${datasetname}/Corrected_objects"
        errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
        memory = { 10.GB + 20.GB * (task.attempt - 1) }
        tag "MNNCorrect $datasetname"
        
        input:
        set val(datasetname), file(datain) from RUN_MNN 
          
        output:
        set val(datasetname), val('MNNCorrect'), file('MNNCorrect.*.rds') into MNN_ENT, MNN_2_H5, MNN_MARK

        """
        mnnCorrect.R\
                --input_object ${datain}\
                --batch_key ${params.batch_key}\
                --output_object MNNCorrect.${datasetname}.rds
        """
}

process Harmony{
        publishDir "./data/${datasetname}/Corrected_objects"
        errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
        memory = { 10.GB + 20.GB * (task.attempt - 1) }
        tag "Harmony $datasetname"

        input:
        set val(datasetname), file(datain) from RUN_HAR

        output:
        set val(datasetname), val('Harmony'), file('Harmony.*.rds') into HAR_ENT, HAR_2_H5

        """
        harmony.R\
                --input_object ${datain}\
                --batch_key ${params.batch_key}\
                --output_object Harmony.${datasetname}.rds
        """
}

process BBKNN{
	publishDir "./data/${datasetname}/Corrected_objects"
        errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
        memory = { 10.GB + 20.GB * (task.attempt - 1) }
        tag "BBKNN $datasetname"

        input:
        set val(datasetname), file(datain) from RUN_BBK

        output:
        set val(datasetname), val('BBKNN'), file('BBKNN.*.h5ad') into BBK_ENT, BBK_ASW

        """
        BBKNN.py\
		--input_object ${datain}\
                --batch_key ${params.batch_key}\
                --output_object BBKNN.${datasetname}.h5ad
        """
}

process Scanorama{
        publishDir "./data/${datasetname}/Corrected_objects"
        errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
        memory = { 10.GB + 20.GB * (task.attempt - 1) }
        tag "Scanorama $datasetname"

        input:
        set val(datasetname), file(datain) from RUN_SCA

        output:
        set val(datasetname), val('Scanorama'), file('Scanorama.*.h5ad') into SCA_ENT, SCA_ASW, SCA_2_SEU

        """
        Scanorama.py\
                --input_object ${datain}\
                --batch_key ${params.batch_key}\
                --output_object Scanorama.${datasetname}.h5ad
        """
}


RAW_ENT.mix(SEU_ENT, COM_ENT, LIM_ENT, MNN_ENT, HAR_ENT).set{R_ENTROPY}

process r_entropy{
        publishDir "./data/${datasetname}/entropy"
        errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
        memory = { 10.GB + 20.GB * (task.attempt - 1) }
        tag "Entropy $datain $method $datasetname"

        input:
        set val(datasetname), val(method), file(datain) from R_ENTROPY

        output:
        file('*.csv')

        """
        R_entropy.R\
	         --input_object ${datain}\
                --batch_key ${params.batch_key}\
                --celltype_key ${params.celltype_key}\
                --output_entropy entropy_${method}.${datasetname}.csv
        """
}

BBK_ENT.mix(SCA_ENT).set{PY_ENTROPY}

process py_entropy {
        publishDir "./data/${datasetname}/entropy"
        errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
        memory = { 10.GB + 20.GB * (task.attempt - 1) }
        tag "entropy (python) $datain $method $datasetname"

        input:
        set val(datasetname), val(method), file(datain) from PY_ENTROPY
        output:
        file('*.csv')

        """
        py_entropy.py\
                --input ${datain}\
                --batch_key ${params.batch_key}\
                --celltype_key ${params.celltype_key}\
                --output_entropy  entropy_${method}.${datasetname}.csv
        """
}

SEU_2_H5.mix(COM_2_H5, LIM_2_H5, MNN_2_H5, HAR_2_H5).set{RES_2_H5}


process Results2H5AD{
        publishDir "./data/${datasetname}/Corrected_objects"
        errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
        memory = { 10.GB + 20.GB * (task.attempt - 1) }
        tag "Convert to h5ad $method $datasetname"

        input:
        set val(datasetname), val(method), file(datain) from RES_2_H5

        output:
        set val(datasetname), val(method), file('*.h5ad') into R_ASW

        """
        results2h5ad.R\
                --input_object ${datain}\
                --output_object ${method}.${datasetname}.h5ad
        """
}


RAW_ASW.mix(BBK_ASW, SCA_ASW, R_ASW).set{ASW}

process asw {
        publishDir "./data/${datasetname}/asw"
        errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
        memory = { 10.GB + 20.GB * (task.attempt - 1) }
        tag "ASW $datain $method $datasetname"

        input:
        set val(datasetname), val(method), file(datain) from ASW
        
        output:
        file('*.csv')
	file('*.png')

        """
        asw.py\
                --input ${datain}\
                --batch_key ${params.batch_key}\
                --celltype_key ${params.celltype_key}\
                --output_asw  asw_${method}.${datasetname}.csv\
                --output_png _${method}.${datasetname}.png
        """
}

process Sca_2_seu{
        publishDir "./data/${datasetname}/Corrected_objects"
        errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
        memory = { 10.GB + 20.GB * (task.attempt - 1) }
        tag "Scanorama to Seurat $datasetname"

        input:
        set val(datasetname), val(method), file(datain) from SCA_2_SEU

        output:
        set val(datasetname), val('Scanorama'), file('Scanorama.*.rds') into SCA_MARK

        """
        h5ad2seurat.R\
                --input_object ${datain}\
                --output_object Scanorama.${datasetname}.rds
        """
}

SEU_MARK.mix(COM_MARK, LIM_MARK, MNN_MARK, SCA_MARK).set{FIND_MARKERS}
process markers {
        publishDir "./data/${datasetname}/jaccard_sim"
        errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
        memory = { 10.GB + 20.GB * (task.attempt - 1) }
        tag "Markers similarity $datain $method $datasetname"

        input:
        set val(datasetname), val(method), file(datain) from FIND_MARKERS

        output:
        file('*.csv')

        """
        find_markers.R\
                --input_object ${datain}\
                --batch_key ${params.batch_key}\
                --celltype_key ${params.celltype_key}\
                --output_file  jacc_${method}.${datasetname}.csv
        """
}

