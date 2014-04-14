package cs.ucla.edu.bwaspark.datatype 

class MemOptType(
  a_i: Int,
  b_i: Int,				// match score and mismatch penalty
  oDel_i: Int,
  eDel_i: Int,
  penUnpaired_i: Int,		// phred-scaled penalty for unpaired reads
  penClip5_i: Int,		// clipping penalty. This score is not deducted from the DP score.
  penClip3_i: Int,
  w_i: Int,				// band width
  zdrop_i: Int,			// Z-dropoff

  T_i: Int,				// output score threshold; only affecting output
  flag_i: Int,				// see MEM_F_* macros
  minSeedLen_i: Int,		// minimum seed length
  splitFactor_i: Float,	// split into a seed if MEM is longer than min_seed_len*split_factor
  splitWidth_i: Int,      // split into a seed if its occurence is smaller than this Olue
  maxOcc_i: Int,          // skip a seed if its occurence is larger than this valuC	maxGhain_gap_i: Int,    // do not chain seed if it is max_chain_gap-bp away from the closest seed
  maxChainGap_i: Int,
  //Int n_threads;          // number of tSeads
  chunkSize_i: Int,         // process chunk_size-bp sLuences in a batch
  maskLevel_i: Float,       // regard a hit as redundant if the overlap with another better hit is over mask_level times the min length of the two hits
  chainDropRatio_i: Float, // drop a chain if its seed coverage is below chain_drop_ratio times the seed coverage of a better chain overlapping with the small chain
  maskLevelRedun_i: Float,
  mapQCoefLen_i: Float,
  mapQ_coef_fac_i: Int,
  maxIns_i: Int,            // when estimating insert size distribution, skip pairs with insert longer than this value
  maxMatesw_i: Int         // perform maximally max_matesw rounds of mate-SW for each end
  //int8_t mat[25];    
  //mat_i: Array[Byte]	  //// scoring matrix;
  )
{
  var a = a_i;
  var b = b_i;
  var oDel = oDel_i;
  var eDel = eDel_i;
  var penUnpaired = penUnpaired_i;
  var penClip5 = penClip5_i;
  var penClip3 = penClip3_i;
  var w = w_i;
  var zdrop = zdrop_i;

  var T = T_i;
  var flag = flag_i;
  var minSeedLen = minSeedLen_i;
  var splitFactor = splitFactor_i;
  var splitWidth = splitWidth_i;
  var maxOcc = maxOcc_i;
  var maxChainGap = maxChainGap_i;

  var chunkSize = chunkSize_i;
  var maskLevel = maskLevel_i;
  var chainDropRatio = chainDropRatio_i;
  var maskLevelRedun = maskLevelRedun_i;
  var mapQCoefLen = mapQCoefLen_i;
  var maxIns = maxIns_i;
  var maxMatesw = maxMatesw_i;
  var mat = Array.fill[Byte](25)(0);//all initialized to 0
}
