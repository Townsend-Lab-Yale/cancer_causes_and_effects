#' @import cancereffectsizeR
#' 
#' Get signature contributions to variants
#' 
#' Builds a data table containing the contributions of each SNV signature to each variant in CES selection output
#' for all tumors with each variant.
#' 
#' For each selected variant and each tumor with the variant, set the contribution of each SNV signature to be the 
#' weight of the signature in the tumor multipled by the relative rate of the trinuc context within the signature, 
#' divided by the relative rate of the trinuc SNV in the tumor. All required information (SNV signature definitions,
#' tumor signature weights, trinuc-specific mutation rates, variant annotations) is taken from the input CESAnalysis.
#' 
#' @param ces_output CESAnalysis with selection intensities calculated for variants of interest
#' @param provided_signature_weights For use when get_signature_weights() will not work on CESAnalysis
#' because you provided the rates directly.
#' @param min_variant_freq Minimum number of affected tumors for a variant to be included in output
#' @param remove_nonsplice_silent Exclude synonymous coding mutations unless they are near splice sites
#' @return A data table covering all tumors and variants

population_scaled_effect_per_tumor <- function(ces_output, provided_signature_weights = NULL,min_variant_freq = 2, remove_nonsplice_silent = F) {
  if(! is(ces_output, "CESAnalysis")) {
    stop("ces_ouput should be a CESAnalysis object", call. = F)
  }
  
  # Get signature definition matrix (rows = signature names, columns = relative rates of trinuc-context-specific SNVs)
  # This works with all versions of CES that have a reference_data accessor (earliest versions put
  # definitions directly at $snv_signatures, newer have a list of information there)
  signature_defs = ces_output$reference_data$snv_signatures
  if (! is(signature_defs, "data.frame")) {
    signature_defs = signature_defs$signatures
  }
  signature_names = rownames(signature_defs)
  
  # Get tumor-specific signature weights
  if(is.null(provided_signature_weights)){
    signature_weights = get_signature_weights(ces_output, include_tumors_without_data = TRUE)
  }else{
    signature_weights = provided_signature_weights 
  }
  
  # For efficiency, drop signatures where weights are always 0
  not_zero = signature_weights[, (sapply(.SD, function(x) any(x != 0))), .SDcols = signature_names]
  not_zero = names(which(not_zero == T))
  signature_weights = signature_weights[, c("Unique_Patient_Identifier", ..not_zero)]
  setkey(signature_weights, "Unique_Patient_Identifier")
  
  signature_names = signature_names[signature_names %in% not_zero]
  signature_defs = signature_defs[not_zero, ]
  
  flux_colnames = paste(signature_names, "flux", sep = "_")
  nrsi_colnames = paste(signature_names, "nrsi", sep = "_")
  
  
  # Get SIs of variants of interest
  selection_data = ces_output$selection[tumors_with_variant >= min_variant_freq]
  
  
  # if requested, remove silent variants that are near a splice site
  if (remove_nonsplice_silent) {
    aac_ids = selection_data[variant_type == "aac", variant]
    silent_nonsplice = ces_output$mutations$amino_acid_change[aac_id %in% aac_ids & next_to_splice == F & aa_ref == aa_alt, aac_id]
    selection_data= selection_data[! variant %in% silent_nonsplice]
  }
  
  
  trinuc_rates = ces_output$trinuc_rates
  setkey(trinuc_rates, "Unique_Patient_Identifier")
  
  
  # Ensure that we're only handling AACs and SNVs
  other_variant_type_index = selection_data[! variant_type %in% c("snv", "aac"), which = T]
  if (length(other_variant_type_index) > 0) {
    selection_data = selection_data[! other_variant_type_index, ]
    warning("Some variants are neither coding mutations nor noncoding SNVs; these were skipped.")
  }
  
  
  # Colllect SNV IDs: just the variant ID for SNVs; for AACs, pull information from mutation annotations
  selection_data[variant_type == "snv", snv_ids := list(list(variant)), by = "variant"]
  
  selection_data[variant_type == "aac", snv_ids := ces_output$mutations$amino_acid_change[variant, all_snv_ids]]
  missing_aac_index = selection_data[, bad := is.null(snv_ids[[1]]), by = "variant"][bad == T, which = T]
  if (length(missing_aac_index) > 0) {
    missing = selection_data[missing_aac_index, variant]
    selection_data = selection_data[! missing_aac_index]
    for (variant in missing) {
      warning(sprintf("Variant %s has incomplete annotation, so it was skipped.", variant))
    }
  }
  
  variant_snv_pairings = selection_data[, .(snv = unlist(snv_ids), selection_intensity), by = "variant"]
  variant_snv_pairings[, trinuc_context := ces_output$mutations$snv[snv, trinuc_mut]]
  
  # Ensure that all SNVs have trinuc context annotation
  bad_snv_index = variant_snv_pairings[is.na(trinuc_context), which = T]
  if (length(bad_snv_index) > 0) {
    missing = variant_snv_pairings[bad_snv_index, variant]
    for (variant in missing) {
      warning(sprintf("Variant %s has incomplete SNV annotations, so it was skipped.", variant))
    }
    variant_snv_pairings = variant_snv_pairings[! variant %in% missing]
  }
  
  # Subset MAF to just the SNVs we need for quicker access
  maf = ces_output$maf[snv_id %in% variant_snv_pairings$snv]
  setkey(maf, "snv_id")
  
  # For each SNV that constitutes a coding mutation (or just the one mutation for nonconding SNVs), 
  # get the contribution of each signature to each tumor that has the mutation.
  # Note that for AACs, cancereffectsizeR's exclusion of di- and tri-nucleotide substitutions (as well as indels)
  # makes it so that each tumor with the AAC has exactly one of the constituent SNVs.
  process_snv = function(variant_id, snv, trinuc_context, selection_intensity) {
    
    tumors = maf[snv, Unique_Patient_Identifier, nomatch = NULL]
    if (length(tumors) == 0) {
      return(NULL)
    }
    # relative rates of the given trinuc SNV for each signature
    sig_trinuc_rate = as.numeric(signature_defs[, trinuc_context])
    snv_output = data.table(variant = variant_id, Unique_Patient_Identifier = tumors)
    
    tumor_trinuc_rates = trinuc_rates[tumors, ..trinuc_context][[1]] # get vector, not 1-column DT
    weights = signature_weights[tumors, ..signature_names]
    
    # For each row (tumor) in weights table, do (weights * sig_trinuc_rate) / tumor_trinuc_rate
    # Opaque syntax for efficiency
    contrib = setDT(Map(`*`, weights, sig_trinuc_rate))
    if(contrib[, .N] == 1) {
      contrib = as.data.table(t(apply(contrib, 2, function(x) x/tumor_trinuc_rates)))
    } else {
      contrib = as.data.table(apply(contrib, 2, function(x) x / tumor_trinuc_rates))
    }
    
    # Multiply contributions by the scalar selection intensity
    # Absurdly, the obscure syntax runs 5x faster than just using "nrsi = contrib * selection_intensity"
    nrsi = setDT(Map(`*`, contrib, selection_intensity))
    snv_output[, (flux_colnames) := contrib]
    snv_output[, (nrsi_colnames) := nrsi]
    return(snv_output)
  }
  
  return(rbindlist(mapply(process_snv, variant_snv_pairings$variant, variant_snv_pairings$snv, 
                          variant_snv_pairings$trinuc_context, variant_snv_pairings$selection_intensity)))
}



