package cs.ucla.edu.bwaspark.datatype

class MemAlnRegType {
  var rBeg: Long = _       // [rBeg,rEnd): reference sequence in the alignment
  var rEnd: Long = _       
  var qBeg: Int = _        // [qBeg,qEnd): query sequence in the alignment
  var qEnd: Int = _
  var score: Int = _       // best local SW score
  var trueScore: Int = _   // actual score corresponding to the aligned region; possibly smaller than $score
  var sub: Int = _         // 2nd best SW score
  var csub: Int = _        // SW score of a tandem hit
  var subNum: Int = _      // approximate number of suboptimal hits
  var width: Int = _       // actual band width used in extension
  var seedCov: Int = _     // length of regions coverged by seeds
  var secondary: Int = _   // index of the parent hit shadowing the current hit; <0 if primary
  var hash: Long = _
}

