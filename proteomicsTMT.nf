#!/usr/bin/env nextflow
/*
=======================================================================================


ProteomicsTMT


 #### Homepage / Documentation
 +++++++++++++++++++++++++++++
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/proteomicslfq --spectra '*.mzML' --database '*.fasta' -profile docker

    Main arguments:
      --input                       Path/URI to PRIDE Sample to data relation format file (SDRF) OR path to input spectra as mzML or Thermo Raw

      For SDRF:
      --root_folder                 (Optional) If given, looks for the filenames in the SDRF in this folder, locally
      --local_input_type            (Optional) If given and 'root_folder' was specified, it overwrites the filetype in the SDRF for local lookup and matches only the basename.

      For mzML/raw files:
      --expdesign                   (Optional) Path to an experimental design file (if not given, it assumes unfractionated, unrelated samples)

      And:
      --database                    Path to input protein database as fasta

    Decoy database:
      --add_decoys                  Add decoys to the given fasta
      --decoy_affix                 The decoy prefix or suffix used or to be used (default: DECOY_)
      --affix_type                  Prefix (default) or suffix (WARNING: Percolator only supports prefices)

    Isobaric analyze:
      --min_precursor_intensity     Minimum intensity of the precursor to be extracted. (Default 1.0)
      --normalization 				Enable normalization of channel intensities with respect to the reference channel.(default false)
      --reference_channel			Number of the reference channel
      --isotope_correction			Enable isotope correction (highly recommended,default true)


    Database Search:
      --search_engines               Which search engine: "comet" (default) or "msgf"
      --enzyme                      Enzymatic cleavage (e.g. 'unspecific cleavage' or 'Trypsin' [default], see OpenMS enzymes)
      --num_enzyme_termini          Specify the termini where the cleavage rule has to match (default:
                                         'fully' valid: 'semi', 'fully')
      --num_hits                    Number of peptide hits per spectrum (PSMs) in output file (default: '1')
      --fixed_mods                  Fixed modifications ('Carbamidomethyl (C)', see OpenMS modifications)
      --variable_mods               Variable modifications ('Oxidation (M)', see OpenMS modifications)
      --precursor_mass_tolerance    Mass tolerance of precursor mass (default: 5)
      --precursor_mass_tolerance_unit Da or ppm (default: ppm)
      --fragment_mass_tolerance     Mass tolerance for fragment masses (currently only controls Comets fragment_bin_tol) (default: 0.03)
      --fragment_mass_tolerance_unit Da or ppm (default: Da)
      --allowed_missed_cleavages    Allowed missed cleavages (default: 2)
      --min_precursor_charge        Minimum precursor ion charge (default: 2)
      --max_precursor_charge        Maximum precursor ion charge (default: 4)
      --min_peptide_length          Minimum peptide length to consider (default: 6)
      --max_peptide_length          Maximum peptide length to consider (default: 40)
      --instrument                  Type of instrument that generated the data (currently only 'high_res' [default] and 'low_res' supported)
      --protocol                    Used labeling or enrichment protocol (if any, phospho,iTRAQ,iTRAQ_phospho,TMT)
      --fragment_method             Used fragmentation method (currently unused since we let the search engines consider all MS2 spectra and let them determine from the spectrum metadata)
      --max_mods                    Maximum number of modifications per peptide. If this value is large, the search may take very long
            --db_debug                    Debug level during database search

      //TODO probably also still some options missing. Try to consolidate them whenever the two search engines share them

    Peak picking:
      --openms_peakpicking          Use the OpenMS PeakPicker to ADDITIONALLY pick the spectra before the search. This is usually done
                                    during conversion already. Only activate if something goes wrong.
      --peakpicking_inmemory        Perform OpenMS peakpicking in-memory. Needs at least the size of the mzML file as RAM but is faster. default: false
      --peakpicking_ms_levels       Which MS levels to pick. default: [] which means auto-convert all non-centroided

    Peptide Re-indexing:
      --IL_equivalent               Should isoleucine and leucine be treated interchangeably? Default: true
      --allow_unmatched             Ignore unmatched peptides (Default: false; only activate if you double-checked all other settings)

    PSM Rescoring:
      --posterior_probabilities     How to calculate posterior probabilities for PSMs:
                                    "percolator" = Re-score based on PSM-feature-based SVM and transform distance
                                        to hyperplane for posteriors
                                    "fit_distributions" = Fit positive and negative distributions to scores
                                        (similar to PeptideProphet)
      --rescoring_debug             Debug level during PSM rescoring
      --psm_pep_fdr_cutoff          FDR cutoff on PSM level (or potential peptide level; see Percolator options) before going into
                                    feature finding, map alignment and inference.

      Percolator specific:
      --train_FDR                   False discovery rate threshold to define positive examples in training. Set to testFDR if 0
      --test_FDR                    False discovery rate threshold for evaluating best cross validation result and reported end result
      --percolator_fdr_level        Level of FDR calculation ('peptide-level-fdrs' or 'psm-level-fdrs')
      --description_correct_features Description of correct features for Percolator (0, 1, 2, 4, 8, see Percolator retention time and calibration)
      --generic_feature_set         Use only generic (i.e. not search engine specific) features. Generating search engine specific
                                    features for common search engines by PSMFeatureExtractor will typically boost the identification rate significantly.
      --subset_max_train            Only train an SVM on a subset of PSMs, and use the resulting score vector to evaluate the other
                                    PSMs. Recommended when analyzing huge numbers (>1 million) of PSMs. When set to 0, all PSMs are used for training as normal.
      --klammer                     Retention time features are calculated as in Klammer et al. instead of with Elude

      Distribution specific:
      --outlier_handling            How to handle outliers during fitting:
                                    - ignore_iqr_outliers (default): ignore outliers outside of 3*IQR from Q1/Q3 for fitting
                                    - set_iqr_to_closest_valid: set IQR-based outliers to the last valid value for fitting
                                    - ignore_extreme_percentiles: ignore everything outside 99th and 1st percentile (also removes equal values like potential censored max values in XTandem)
                                    - none: do nothing
      --top_hits_only               Use only the top hits for fitting

      //TODO add more options for rescoring part

    IDFilter:
      --protgroup_score_cutoff      The score which should be reached by a protein group to be kept.Use in combination with 'delete_unreferenced_peptide_
                                    hits' to remove affected peptides.  (Default: 0.0)


    IDMapper:
      --rt_tolerance				RT tolerance (in seconds) for the matching of peptide identifications and consensus features (Default 5.0).
      --mz_tolerance				m/z tolerance (in ppm or Da) for the matching of peptide identifications and consensus features(Default 20.0).
      --mz_measure					Unit of 'mz_tolerance'. ('ppm', 'Da',Default 'ppm')
      --mz_reference				Source of m/z values for peptide identifications. If 'precursor', the precursor-m/z from the idXML is used. 
									If 'peptide',masses are computed from the sequences of peptide hits.(Defalut 'peptide')

	FileMerger:
	  --annotate_file_origin		Store the original filename in each feature using meta value "file_origin".(Default false).
	  --append_method				Append consensusMaps rowise or colwise.(Default append_rows).




    Inference:
      --protein_fdr            		Additionally calculate the target-decoy FDR on protein-level based on the posteriors(Default false).
      --greedy_group_resolution		Default none.
      --top_PSMs					Consider only top X PSMs per spectrum. 0 considers all.(Default 1).

    IDConflictResolver:
      --resolve_between_features	A map may contain multiple features with both identical (possibly modified i.e. not stripped) sequence and charge state. 							   The feature with the 'highest intensity' is very likely the most reliable one.Default(off),highest_intensity

    ProteinQuantifier:
      --top 						Calculate protein abundance from this number of proteotypic peptides (most abundant first; '0' for all, Default 3)
      --average						Averaging method used to compute protein abundances from peptide abundances.
      								(median,mean,weighted_mean,sum Default median).
      --best_charge_and_fraction	Distinguish between fraction and charge states of a peptide.(Default false).
      --ratios						Add the log2 ratios of the abundance values to the output.(Default false)
      --normalize					Scale peptide abundances so that medians of all samples are equal.(Default false)
      --fix_peptides				Use the same peptides for protein quantification across all samples.(Default false)



    Other options:
      --outdir [file]                 The output directory where the results will be saved
      --publish_dir_mode [str]        Mode for publishing results in the output directory. Available: symlink, rellink, link, copy, copyNoFollow, move (Default: copy)
      --email [email]                 Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail [email]         Same as --email, except only send mail if the workflow is not successful
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch options:
      --awsqueue [str]                The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]               The AWS Region for your AWS Batch job to run on
      --awscli [str]                  Path to the AWS CLI tool

    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name

if (!(workflow.runName == ~/[a-z]+_[a-z]+/)) {
	custom_runName = workflow.runName
}

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}


// Stage config files
//ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)
//ch_output_docs_images = file("$baseDir/docs/images/", checkIfExists: true)

// Validate input
if (isCollectionOrArray(params.input))
{
  tocheck = params.input[0]
} else {
  tocheck = params.input
}
tocheck = params.input

sdrf_file = null

if (tocheck.toLowerCase().endsWith("sdrf") || tocheck.toLowerCase().endsWith("tsv")) {
  sdrf_file = params.input
} else if (tocheck.toLowerCase().endsWith("mzml") || tocheck.toLowerCase().endsWith("raw")) {
  spectra_files = params.input
} else {
	log.error "EITHER spectra data (mzml/raw) OR an SDRF needs to be  provided as input."; exit 1
}

params.database = params.database ?: { log.error "No protein database provided. Make sure you have used the '--database' option."; exit 1 }()
params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()


/*
 * Create a channel for input files
 */

 //Filename        FixedModifications      VariableModifications   Label   PrecursorMassTolerance  PrecursorMassToleranceUnit      FragmentMassTolerance   DissociationMethod      Enzyme


if (!sdrf_file) {
  ch_spectra = Channel.fromPath(spectra_files, checkIfExists: true)
  ch_spectra
  .multiMap{ it -> id = it.toString().md5()
                    comet_settings: msgf_settings: tuple(id,
                                    params.fixed_mods,
                                    params.variable_mods,
                                    "", //labelling modifications currently not supported
                                    params.precursor_mass_tolerance,
                                    params.precursor_mass_tolerance_unit,
                                    params.fragment_mass_tolerance,
                                    params.fragment_mass_tolerance_unit,
                                    params.fragment_method,
                                    params.enzyme)
                    idx_settings: tuple(id,
                                    params.enzyme)
                    luciphor_settings:
                                  tuple(id,
                                    params.fragment_method)
                    mzmls: tuple(id,it)}
  .set{ch_sdrf_config}
} else {
  ch_sdrf = Channel.fromPath(sdrf_file, checkIfExists: true)
  /*
   * STEP 0 - SDRF parsing
   */
  process sdrf_parsing {

      publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

      input:
       file sdrf from ch_sdrf

      output:
       file "experimental_design.tsv" into ch_expdesign, ch_pro_quant_exp
       file "openms.tsv" into ch_sdrf_config_file

      when:
        sdrf_file

      script:
       """
       ## -t2 since the one-table format parser is broken in OpenMS2.5
       ## -l for legacy behavior to always add sample columns
       parse_sdrf convert-openms -t1 -l -s ${sdrf} > sdrf_parsing.log
       """
  }

  //TODO use header and reference by col name instead of index
  ch_sdrf_config_file
  .splitCsv(skip: 1, sep: '\t')
  .multiMap{ row -> id = row.toString().md5()
  					isobaricanalyzer_settings: tuple(id,
  									row[4],
  									row[9])
                    msgf_settings: tuple(id,
                                    row[2],
                                    row[3],
                                    row[4],
                                    row[5],
                                    row[6],
                                    row[7],
                                    row[8],
                                    row[9],
                                    row[10])
                    idx_settings: tuple(id,
                                    row[10])
                    luciphor_settings:
                                  tuple(id,
                                    row[9])
                    mzmls: tuple(id, !params.root_folder ?
                                    row[0] :
                                    params.root_folder + "/" + (params.local_input_type ?
                                        row[1].take(row[1].lastIndexOf('.')) + '.' + params.local_input_type :
                                        row[1]))}
  .set{ch_sdrf_config}
}



ch_db_for_decoy_creation = Channel.fromPath(params.database)

// overwrite experimental design if given additionally to SDRF
//TODO think about that
if (params.expdesign)
{
    Channel
        .fromPath(params.expdesign)
        .set { ch_expdesign }
        .set { ch_pro_quant_exp }
}


ch_sdrf_config.mzmls
.branch {
        raw: hasExtension(it[1], 'raw')
        mzML: hasExtension(it[1], 'mzML')
}
.set {branched_input}



//TODO we could also check for outdated mzML versions and try to update them
branched_input.mzML
.branch {
    nonIndexedMzML: file(it[1]).withReader {
                        f = it;
                        1.upto(5) {
                            if (f.readLine().contains("indexedmzML")) return false;
                        }
                        return true;
                    }
    inputIndexedMzML: file(it[1]).withReader {
                        f = it;
                        1.upto(5) {
                            if (f.readLine().contains("indexedmzML")) return true;
                        }
                        return false;
                    }
}
.set {branched_input_mzMLs}

//Push raw files through process that does the conversion, everything else directly to downstream Channel with mzMLs

//This piece only runs on data that is a.) raw and b.) needs conversion
//mzML files will be mixed after this step to provide output for downstream processing - allowing you to even specify mzMLs and RAW files in a mixed mode as input :-)

/*
 * STEP 0.1 - Raw file conversion
 */
process raw_file_conversion {

    label 'process_low'
    label 'process_single_thread'


    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     tuple mzml_id, path(rawfile) from branched_input.raw

    output:
     tuple mzml_id, file("*.mzML") into mzmls_converted

    script:
     """
     ThermoRawFileParser.sh -i=${rawfile} -f=2 -o=./ > ${rawfile}_conversion.log
     """
}

/*
 * STEP 0.2 - MzML indexing
 */
process mzml_indexing {

    label 'process_low'

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     tuple mzml_id, path(mzmlfile) from branched_input_mzMLs.nonIndexedMzML

    output:
     tuple mzml_id, file("out/*.mzML") into mzmls_indexed
     file "*.log"

    script:
     """
     mkdir out
     FileConverter -in ${mzmlfile} -out out/${mzmlfile.baseName}.mzML > ${mzmlfile.baseName}_mzmlindexing.log
     """
}

//Mix the converted raw data with the already supplied mzMLs and push these to the same channels as before

if (params.openms_peakpicking)
{
  branched_input_mzMLs.inputIndexedMzML.mix(mzmls_converted).mix(mzmls_indexed).set{mzmls_pp}
  (mzmls_isobaric_analyzer, mzmls_msgf, mzmls_luciphor, mzmls_plfq) = [Channel.empty(), Channel.empty(), Channel.empty(), Channel.empty()]
}
else
{
  branched_input_mzMLs.inputIndexedMzML.mix(mzmls_converted).mix(mzmls_indexed).into{mzmls_isobaric_analyzer;mzmls_msgf; mzmls_luciphor; mzmls_plfq}
  mzmls_pp = Channel.empty()
}

//Fill the channels with empty Channels in case that we want to add decoys. Otherwise fill with output from database.
(searchengine_in_db_msgf, pepidx_in_db, plfq_in_db) = ( params.add_decoys
                    ? [ Channel.empty(), Channel.empty(), Channel.empty(), Channel.empty() ]
                    : [ Channel.fromPath(params.database), Channel.fromPath(params.database), Channel.fromPath(params.database), Channel.fromPath(params.database) ] )

//Add decoys if params.add_decoys is set appropriately
process generate_decoy_database {

    label 'process_very_low'
    label 'process_single_thread'

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     file(mydatabase) from ch_db_for_decoy_creation

    output:
     file "${mydatabase.baseName}_decoy.fasta" into searchengine_in_db_decoy_msgf, pepidx_in_db_decoy, plfq_in_db_decoy
     file "*.log"

    when:
     params.add_decoys

    script:
     """
     DecoyDatabase  -in ${mydatabase} \\
                 -out ${mydatabase.baseName}_decoy.fasta \\
                 -decoy_string ${params.decoy_affix} \\
                 -decoy_string_position ${params.affix_type} \\
                 > ${mydatabase.baseName}_decoy_database.log
     """
}


process openms_peakpicker {

    label 'process_low'

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     tuple mzml_id, path(mzml_file) from mzmls_pp

    when:
      params.openms_peakpicking

    output:
     tuple mzml_id, file("out/${mzml_file.baseName}.mzML") into isobaric_analyzer_picked, mzmls_msgf_picked, mzmls_plfq_picked
     file "*.log"

    script:
     in_mem = params.peakpicking_inmemory ? "inmemory" : "lowmemory"
     lvls = params.peakpicking_ms_levels ? "-algorithm:ms_levels ${params.peakpicking_ms_levels}" : ""
     """
     mkdir out
     PeakPickerHiRes -in ${mzml_file} \\
                     -out out/${mzml_file.baseName}.mzML \\
                     -threads ${task.cpus} \\
                     -debug ${params.pp_debug} \\
                     -processOption ${in_mem} \\
                     ${lvls} \\
                     > ${mzml_file.baseName}_pp.log
     """
}

if (params.enzyme == "unspecific cleavage")
{
  params.num_enzyme_termini == "none"
}

pepidx_num_enzyme_termini = params.num_enzyme_termini
if (params.num_enzyme_termini == "fully")
{
  pepidx_num_enzyme_termini = "full"
}


process isobaric_analyzer {
	label 'process_medium'
	
	publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

	input:
	  tuple mzml_id, path(mzml_file), label,  diss_meth from mzmls_isobaric_analyzer.mix(isobaric_analyzer_picked).join(ch_sdrf_config.isobaricanalyzer_settings)

	when:
	  params.isobaric_analyzer


	output:
	  tuple mzml_id, file("${mzml_file.baseName}_iso.consensusXML") into id_files_consensusXML
      file "*.log"

	script:
	  if (diss_meth == 'HCD') diss_meth = 'High-energy collision-induced dissociation'
	  else if (diss_meth == 'CID') diss_meth = 'Collision-induced dissociation'
	  else if (diss_meth == 'ETD') diss_meth = 'Electron transfer dissociation'
	  else if (diss_meth == 'ECD') diss_meth = 'Electron capture dissociation'
	  """
	  IsobaricAnalyzer -type ${label} \\
	  				   -in ${mzml_file} \\
	  				   -threads ${task.cpus} \\
	  				   -quantification:normalization \\
	  				   -extraction:select_activation "${diss_meth}" \\
	  				   -extraction:min_reporter_intensity ${params.min_reporter_intensity} \\
	  				   -extraction:min_precursor_purity ${params.min_precursor_purity} \\
	  				   -extraction:precursor_isotope_deviation ${params.precursor_isotope_deviation} \\
	  				   -debug ${params.iso_debug} \\
	  				   -out ${mzml_file.baseName}_iso.consensusXML \\
	  				   > ${mzml_file.baseName}_isob.log
	  """


}
process search_engine_msgf {

    label 'process_medium'

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    // ---------------------------------------------------------------------------------------------------------------------
    // ------------- WARNING: THIS IS A HACK. IT JUST DOES NOT WORK IF THIS PROCESS IS RETRIED -----------------------------
    // ---------------------------------------------------------------------------------------------------------------------
    // I actually dont know, where else this would be needed.
    errorStrategy 'terminate'

    input:
     tuple file(database), mzml_id, path(mzml_file), fixed, variable, label, prec_tol, prec_tol_unit, frag_tol, frag_tol_unit, diss_meth, enzyme from searchengine_in_db_msgf.mix(searchengine_in_db_decoy_msgf).combine(mzmls_msgf.mix(mzmls_msgf_picked).join(ch_sdrf_config.msgf_settings))

     // This was another way of handling the combination
     //file database from searchengine_in_db.mix(searchengine_in_db_decoy)
     //each file(mzml_file) from mzmls
    when:
      params.search_engines.contains("msgf")

    output:
     tuple mzml_id, file("${mzml_file.baseName}_msgf.idXML") into id_files_msgf
     file "*.log"

    script:
      if (enzyme == 'Trypsin') enzyme = 'Trypsin/P'
      else if (enzyme == 'Arg-C') enzyme = 'Arg-C/P'
      else if (enzyme == 'Asp-N') enzyme = 'Asp-N/B'
      else if (enzyme == 'Chymotrypsin') enzyme = 'Chymotrypsin/P'
      else if (enzyme == 'Lys-C') enzyme = 'Lys-C/P'

      if ((frag_tol.toDouble() < 50 && frag_tol_unit == "ppm") || (frag_tol.toDouble() < 0.1 && frag_tol_unit == "Da"))
      {
        inst = params.instrument ?: "high_res"
      } else {
        inst = params.instrument ?: "low_res"
      }
     """
     MSGFPlusAdapter -in ${mzml_file} \\
                     -out ${mzml_file.baseName}_msgf.idXML \\
                     -threads ${task.cpus} \\
                     -java_memory ${task.memory.toMega()} \\
                     -database "${database}" \\
                     -instrument ${inst} \\
                     -protocol "${params.protocol}" \\
                     -matches_per_spec ${params.num_hits} \\
                     -min_precursor_charge ${params.min_precursor_charge} \\
                     -max_precursor_charge ${params.max_precursor_charge} \\
                     -min_peptide_length ${params.min_peptide_length} \\
                     -max_peptide_length ${params.max_peptide_length} \\
                     -enzyme "${enzyme}" \\
                     -tryptic ${params.num_enzyme_termini} \\
                     -precursor_mass_tolerance ${prec_tol} \\
                     -precursor_error_units ${prec_tol_unit} \\
                     -fixed_modifications ${fixed.tokenize(',').collect { "'${it}'" }.join(" ") } \\
                     -variable_modifications ${variable.tokenize(',').collect { "'${it}'" }.join(" ") } \\
                     -max_mods ${params.max_mods} \\
                     -debug ${params.db_debug} \\
                     > ${mzml_file.baseName}_msgf.log
     """
}


process index_peptides {

    label 'process_low'

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     tuple mzml_id, file(id_file), val(enzyme), file(database) from id_files_msgf.combine(ch_sdrf_config.idx_settings, by: 0).combine(pepidx_in_db.mix(pepidx_in_db_decoy))

    output:
     tuple mzml_id, file("${id_file.baseName}_idx.idXML") into id_files_idx_ForPerc
     file "*.log"

    script:
     def il = params.IL_equivalent ? '-IL_equivalent' : ''
     def allow_um = params.allow_unmatched ? '-allow_unmatched' : ''
     // see comment in CometAdapter. Alternative here in PeptideIndexer is to let it auto-detect the enzyme by not specifying.
     if (params.search_engines.contains("msgf"))
     {
        if (enzyme == 'Trypsin') enzyme = 'Trypsin/P'
        else if (enzyme == 'Arg-C') enzyme = 'Arg-C/P'
        else if (enzyme == 'Asp-N') enzyme = 'Asp-N/B'
        else if (enzyme == 'Chymotrypsin') enzyme = 'Chymotrypsin/P'
        else if (enzyme == 'Lys-C') enzyme = 'Lys-C/P'
     }
     """
     PeptideIndexer -in ${id_file} \\
                    -out ${id_file.baseName}_idx.idXML \\
                    -threads ${task.cpus} \\
                    -fasta ${database} \\
                    -enzyme:name "${enzyme}" \\
                    -enzyme:specificity ${pepidx_num_enzyme_termini} \\
                    ${il} \\
                    ${allow_um} \\
                    > ${id_file.baseName}_index_peptides.log
     """
}



process extract_percolator_features {

    label 'process_very_low'
    label 'process_single_thread'

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     tuple mzml_id, file(id_file) from id_files_idx_ForPerc

    output:
     tuple mzml_id, file("${id_file.baseName}_feat.idXML") into id_files_idx_feat
     file "*.log"

    when:
     params.posterior_probabilities == "percolator"

    script:
     """
     PSMFeatureExtractor -in ${id_file} \\
                         -out ${id_file.baseName}_feat.idXML \\
                         -threads ${task.cpus} \\
                         > ${id_file.baseName}_extract_percolator_features.log
     """
}


//Note: from here, we do not need any settings anymore. so we can skip adding the mzml_id to the channels
//TODO find a way to run across all runs merged
process percolator {

    //TODO Actually it heavily depends on the subset_max_train option and the number of IDs
    // would be cool to get an estimate by parsing the number of IDs from previous tools.
    label 'process_medium'
    //Since percolator 3.5 it allows for 27 parallel tasks
    cpus { check_max( 27, 'cpus' ) }

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'
    publishDir "${params.outdir}/raw_ids", mode: 'copy', pattern: '*.idXML'

    input:
     tuple mzml_id, file(id_file) from id_files_idx_feat

    output:
     tuple mzml_id, file("${id_file.baseName}_perc.idXML"), val("MS:1001491") into id_files_perc
     file "*.log"

    when:
     params.posterior_probabilities == "percolator"

    // NICE-TO-HAVE: the decoy-pattern is automatically detected from PeptideIndexer.
    // Parse its output and put the correct one here.
    script:
      if (params.klammer && params.description_correct_features == 0) {
          log.warn('Klammer was specified, but description of correct features was still 0. Please provide a description of correct features greater than 0.')
          log.warn('Klammer will be implicitly off!')
      }

      // currently post-processing-tdc is always set since we do not support separate TD databases
      """
      ## Percolator does not have a threads parameter. Set it via OpenMP env variable,
      ## to honor threads on clusters
      OMP_NUM_THREADS=${task.cpus} PercolatorAdapter \\
                          -in ${id_file} \\
                          -out ${id_file.baseName}_perc.idXML \\
                          -threads ${task.cpus} \\
                          -subset_max_train ${params.subset_max_train} \\
                          -decoy_pattern ${params.decoy_affix} \\
                          -post_processing_tdc \\
                          -score_type pep \\
                          > ${id_file.baseName}_percolator.log
      """
}


process idfilter {

    label 'process_very_low'
    label 'process_single_thread'

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'
    publishDir "${params.outdir}/ids", mode: 'copy', pattern: '*.idXML'

    input:
     tuple mzml_id, file(id_file), val(q_val) from id_files_perc

    output:
     tuple mzml_id, file("${id_file.baseName}_filter.idXML") into id_filtered
     file "*.log"

    script:
     """
     IDFilter -in ${id_file} \\
              -out ${id_file.baseName}_filter.idXML \\
              -threads ${task.cpus} \\
              -score:pep ${params.psm_pep_fdr_cutoff} \\
              > ${id_file.baseName}_idfilter.log
     """
}

process idmapper{
	label 'process_medium'

	publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

	input:
	 tuple mzml_id, file(id_file_filter), file(consensusXML) from id_filtered.combine(id_files_consensusXML, by: 0)
	 
	output:
	 file("${id_file_filter.baseName}_map.consensusXML") into id_map_to_merger

	script:
	 """
	 IDMapper -id ${id_file_filter} \\
	 		  -in ${consensusXML} \\
	 		  -threads ${task.cpus} \\
	 		  -rt_tolerance ${params.rt_tolerance} \\
	 		  -mz_tolerance ${params.mz_tolerance} \\
	 		  -mz_measure ${params.mz_measure} \\
	 		  -mz_reference ${params.mz_reference} \\
	 		  -debug 1 \\
	 		  -out ${id_file_filter.baseName}_map.consensusXML \\
	 		  > ${id_file_filter.baseName}_map.log
	 """
}

process file_merge{
	label 'process_medium'
	label 'process_single_thread'

	publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

	input:
	 file(id_map) from id_map_to_merger.collect()

	output:
	 file("ID_mapper_merge.consensusXML") into id_merge_to_epi


	script:
	 """
	 FileMerger -in ${(id_map as List).join(' ')} \\
	 			-annotate_file_origin  \\
	 			-append_method ${params.append_method} \\
	 			-threads ${task.cpus} \\
	 			-debug 1 \\
	 			-out ID_mapper_merge.consensusXML \\
	 			> ID_mapper_merge.log
	 """
}

process epifany{
	label 'process_medium'
    
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

	input:
	 file(consus_file) from id_merge_to_epi
	 file expdes from ch_expdesign


	output:
	 file("${consus_file.baseName}_epi.consensusXML") into epi_idfilter

	 // expdes currently unused
	script:
	 """
	 Epifany -in ${consus_file} \\
	 		 -protein_fdr ${params.protein_fdr} \\
	 		 -threads ${task.cpus} \\
	 		 -debug 1 \\
	 		 -greedy_group_resolution ${params.greedy_group_resolution} \\
			 -algorithm:top_PSMs ${params.top_PSMs} \\
			 -out ${consus_file.baseName}_epi.consensusXML \\
			 > ${consus_file.baseName}_epi.log
	 """                                       
}


process epi_filter{
	label 'process_very_low'
    label 'process_single_thread'

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     file(consus_epi) from epi_idfilter


    output:
     file("${consus_epi.baseName}_filt.consensusXML") into confict_res


    script:
    """
    IDFilter -in ${consus_epi} \\
             -out ${consus_epi.baseName}_filt.consensusXML \\
             -threads ${task.cpus} \\
             -score:protgroup ${params.protgroup_score_cutoff} 
             > ${consus_epi.baseName}_idfilter.log
    """
}


process resolve_conflict{
	label 'process_medium'

	publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

	input:
	 file(consus_epi_filt) from confict_res


	output:
	 file "${consus_epi_filt.baseName}_resconf.consensusXML" into pro_quant, ch_mztabexport


	script:
	"""IDConflictResolver -in ${consus_epi_filt} \\
						  -threads ${task.cpus} \\
						  -debug 1 \\
						  -resolve_between_features ${params.res_between_fet} \\
						  -out ${consus_epi_filt.baseName}_resconf.consensusXML \\
						  > ${consus_epi_filt.baseName}_resconf.log

	"""
}


process mztab_export{
	label 'process_medium'

	publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'
    publishDir "${params.outdir}/mztab", mode: 'copy'

	input:
	 file(cons_epi_filt_resconf) from ch_mztabexport


	output:
	 file "out.mzTab"

	script:
	"""
	MzTabExporter -in ${cons_epi_filt_resconf} \\
				  -out out.mztab \\
				  -debug 1 \\
				  -threads ${task.cpus} \\
				  > mztabexporter.log
	"""
}



process pro_quant{
	label 'process_medium'

	publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'
    publishDir "${params.outdir}/pro_quant", mode: 'copy'


	input:
	 file(epi_filt_resolve) from pro_quant
	 file pro_quant_exp from ch_pro_quant_exp


	output:
	 file "protein_out.csv" optional true into downstreams_tool_A
	 file "peptide_out.csv" into downstreams_tool_B
	 file "out.mzTab" into downstreams_tool_C


	 script:
	 """
	 ProteinQuantifier -in ${epi_filt_resolve} \\
	 				   -out protein_out.csv \\
	 				   -peptide_out peptide_out.csv \\
	 				   -mztab out.mztab \\
	 				   -top ${params.top} \\
	 				   -average ${params.average} \\
	 				   -best_charge_and_fraction \\
	 				   -ratios \\
	 				   -threads ${task.cpus} \\
	 				   -consensus:normalize \\
	 				   -consensus:fix_peptides \\
	 				   > pro_quant.log
	 """
}


//--------------------------------------------------------------- //
//---------------------- Utility functions  --------------------- //
//--------------------------------------------------------------- //

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Check class of an Object for "List" type
boolean isCollectionOrArray(object) {
    [Collection, Object[]].any { it.isAssignableFrom(object.getClass()) }
}