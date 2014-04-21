package cs.ucla.edu.bwaspark.datatype

class MemAlnType {
  var pos: Long = _          // forward strand 5'-end mapping position
  var rid: Int = _           // reference sequence index in bntseq_t; <0 for unmapped
  var flag: Int = _          // extra flag
  var isRev: Boolean = _     // is_rev: whether on the reverse strand
  var mapq: Short = _        // mapq: mapping quality
  var NM: Int = _            // NM: edit distance
  var nCigar: Int = _        // number of CIGAR operations
  var cigar: Array[Int] = _  // CIGAR in the BAM encoding: opLen<<4|op; op to integer mapping: MIDSH=>01234
  var score: Int = _
  var sub: Int = _
}

