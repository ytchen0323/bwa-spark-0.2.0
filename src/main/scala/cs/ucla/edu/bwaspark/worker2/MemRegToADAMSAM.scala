package cs.ucla.edu.bwaspark.worker2

import scala.util.control.Breaks._
import scala.List
import scala.collection.mutable.MutableList
import scala.math.log
import scala.math.abs

import cs.ucla.edu.bwaspark.datatype._
import cs.ucla.edu.bwaspark.util.BNTSeqUtil._

object MemRegToADAMSAM {
  val MEM_F_ALL = 0x8
  val MEM_F_NO_MULTI = 0x10
  val MEM_MAPQ_COEF = 30.0

  /**
    *  Transform the alignment registers to SAM format
    *  
    *  @param opt the input MemOptType object
    *  @param bns the input BNTSeqType object
    *  @param pac the PAC array
    *  @param seq the read (NOTE: currently we use Array[Byte] first. may need to be changed back to Array[Char]!!!)
    *  @param regs the alignment registers to be transformed
    *  @param extraFlag
    *  @param alns 
    */
  def memRegToSAMSe(opt: MemOptType, bns: BNTSeqType, pac: Array[Byte], seq: Array[Byte], regs: MutableList[MemAlnRegType], extraFlag: Int, alnsIn: MutableList[MemAlnType]) {
    var alns: MutableList[MemAlnType] = new MutableList[MemAlnType]

    // NOTE: set opt.flag manually here!!! This should be modified from the logger!!!
    opt.flag = 24
/*   
    println("[opt object] flag: " + opt.flag + " T: " + opt.T + " minSeedLen: " + opt.minSeedLen + " a: " + opt.a + " b: " + opt.b + " mapQCoefLen: " + opt.mapQCoefLen + " mapQCoefFac: " + opt.mapQCoefFac)

    var j = 0
    regs.foreach(r => {
      print("Reg " + j + "(")
      print(r.rBeg + ", " + r.rEnd + ", " + r.qBeg + ", " + r.qEnd + ", " + r.score + ", " + r.trueScore + ", ")
      println(r.sub + ", "  + r.csub + ", " + r.subNum + ", " + r.width + ", " + r.seedCov + ", " + r.secondary + ")")
      j += 1
      } )
*/
    for(i <- 0 to (regs.length - 1)) {
      if(regs(i).score >= opt.T) {
        if(regs(i).secondary < 0 || ((opt.flag & MEM_F_ALL) > 0)) {
          if(regs(i).secondary < 0 || regs(i).score >= regs(regs(i).secondary).score * 0.5) {
            // debugging
            print("i=" + i + " ")
            var aln = memRegToAln(opt, bns, pac, 101, seq, regs(i))   // NOTE: current data structure has not been obtained from RDD. We assume the length to be 101 here
            alns += aln
            aln.flag |= extraFlag   // flag secondary
            if(regs(i).secondary >= 0) aln.sub = -1   // don't output sub-optimal score
            if(i > 0 && regs(i).secondary < 0)   // if supplementary
              if((opt.flag & MEM_F_NO_MULTI) > 0) aln.flag |= 0x10000
              else aln.flag |= 0x800

            if(i > 0 && aln.mapq > alns(0).mapq) aln.mapq = alns(0).mapq            
          }
        }
      }
    }

    // no alignments good enough; then write an unaligned record
    if(alns.length == 0) { 
      var aln = memRegToAln(opt, bns, pac, 101, seq, null)
      aln.flag |= extraFlag
      //memAlnToSAM
    }
    else {
      //alns.foreach(memAlnToSAM)
    }
  }


  //basic hit->SAM conversion
  /**
    *
    @param l1
    @param l2
    @param score
    @param a
    @param q
    @param r
    @param w, output
    */

   def inferBw( l1: Int, l2: Int, score: Int, a: Int, q: Int, r: Int) : Int = {
   //private def inferBw( l1: Int, l2: Int, score: Int, a: Int, q: Int, r: Int) : Int = {
     var w: Int = 0;
     if(l1 == l2 && (((l1 * a) - score) < ((q + r - a) << 1))) {
     }
     else {
       var wDouble: Double = 0
       wDouble = (((if( l1 < l2 ) l1 else l2) * a - score - q).toDouble / r.toDouble + 2.toDouble)
       var absDifference: Int = if( l1 < l2 ) ( l2 - l1) else ( l1 - l2)
       if( wDouble < absDifference ) w = absDifference
       else w = wDouble.toInt
       //w is the returned Integer
     }
     w
   }

  // wrapper implementation only for now
  /**
    *  Transform the alignment registers to alignment type
    *
    *  @param opt the input MemOptType object
    *  @param bns the input BNTSeqType object
    *  @param pac the PAC array
    *  @param seqLen the length of the input sequence
    *  @param seq the input sequence (NOTE: currently we use Array[Byte] first. may need to be changed back to Array[Char]!!!)
    *  @param reg the input alignment register
    */
  private def memRegToAln(opt: MemOptType, bns: BNTSeqType, pac: Array[Byte], seqLen: Int, seq: Array[Byte], reg: MemAlnRegType): MemAlnType = {
    val aln = new MemAlnType
    if(reg != null) {
      println(memApproxMapqSe(opt, reg))
      bwaFixXref2(opt.mat, opt.oDel, opt.eDel, opt.oIns, opt.eIns, opt.w, bns, pac, seq, reg.qBeg, reg.qEnd, reg.rBeg, reg.rEnd)
    }
    aln
  }

  /**
    *  Calculate the approximate mapq value
    *
    *  @param opt the input MemOptType object
    *  @param reg the input alignment register
    */
  private def memApproxMapqSe(opt: MemOptType, reg: MemAlnRegType): Int = {
    var sub = 0
    
    if(reg.sub > 0) sub = reg.sub
    else sub = opt.minSeedLen * opt.a

    if(reg.csub > sub) sub = reg.csub
   
    if(sub >= reg.score) 0
    else {
      var len = 0
      var mapq = 0
      
      if(reg.qEnd - reg.qBeg > reg.rEnd - reg.rBeg) len = reg.qEnd - reg.qBeg
      else len = (reg.rEnd - reg.rBeg).toInt
      
      val identity = 1.0 - (len * opt.a - reg.score).toDouble / (opt.a + opt.b) / len
      if(reg.score == 0) mapq = 0
      else if(opt.mapQCoefLen > 0) {
        var tmp: Double = 1.0
        if(len > opt.mapQCoefLen) tmp = opt.mapQCoefFac / log(len)
        tmp *= identity * identity
        mapq = (6.02 * (reg.score - sub) / opt.a * tmp * tmp + 0.499).toInt
      }
      else {
        mapq = (MEM_MAPQ_COEF * (1.0 - sub.toDouble / reg.score) * log(reg.seedCov) + .499).toInt
        if(identity < 0.95) mapq = (mapq * identity * identity + .499).toInt
      }
   
      if(reg.subNum > 0) mapq -= (4.343 * log(reg.subNum + 1) + .499).toInt
      if(mapq > 60) mapq = 60
      if(mapq < 0) mapq = 0

      mapq
    }

  }  

  // wrapper implementation only for now
  private def bwaFixXref2(mat: Array[Byte], oDel: Int, eDel: Int, oIns: Int, eIns: Int, w: Int, bns: BNTSeqType, 
    pac: Array[Byte], query: Array[Byte], qBeg: Int, qEnd: Int, rBeg: Long, rEnd: Long): Array[Int] = {
    var retArray = new Array[Int](5)
    bwaGenCigar2(mat, oDel, eDel, oIns, eIns, w, bns.l_pac, pac, qEnd - qBeg, query, rBeg, rEnd)
    retArray
  }

  private def bwaGenCigar2(mat: Array[Byte], oDel: Int, eDel: Int, oIns: Int, eIns: Int, w: Int, pacLen: Long, pac: Array[Byte], 
    queryLen: Int, query_i: Array[Byte], rBeg: Long, rEnd: Long): (Int, Int, Int, Array[Int]) = {

    var numCigar = 0
    var NM = -1
    var score = 0
    var cigar: Array[Int] = null

    // reject if negative length or bridging the forward and reverse strand
    if(queryLen <= 0 || rBeg >= rEnd || (rBeg < pacLen && rEnd > pacLen)) (0, 0, 0, null)
    else {
      val ret = bnsGetSeq(pacLen, pac, rBeg, rEnd)
      var rseq = ret._1
      val rlen = ret._2

      // possible if out of range
      if(rEnd - rBeg != rlen) (0, 0, 0, null)
      else {
        var query = query_i

        // then reverse both query and rseq; this is to ensure indels to be placed at the leftmost position
        if(rBeg >= pacLen) {
          for(i <- 0 to ((queryLen >> 1) - 1)) {
            var tmp = query(i)
            query(i) = query(queryLen - 1 - i)
            query(queryLen - 1 - i) = tmp
          }
            
          for(i <- 0 to ((rlen >> 1) - 1).toInt) {
            var tmp = rseq(i)
            rseq(i) = rseq((rlen - 1 - i).toInt)
            rseq((rlen - 1 - i).toInt) = tmp
          }
        }        
        // no gap; no need to do DP
        if(queryLen == (rEnd - rBeg) && w == 0) {
          // FIXME: due to an issue in mem_reg2aln(), we never come to this block. This does not affect accuracy, but it hurts performance. (in original C implementation)
          cigar = new Array[Int](1)
          cigar(0) = queryLen << 4 | 0
          numCigar = 1

          for(i <- 0 to (queryLen - 1)) 
            score += mat(rseq(i) * 5 + query(i))
        }
        else {
          val maxIns = ((((queryLen + 1) >> 1) * mat(0) - oIns).toDouble / eIns + 1.0).toInt
          val maxDel = ((((queryLen + 1) >> 1) * mat(0) - oDel).toDouble / eDel + 1.0).toInt
          var maxGap = maxIns
          if(maxIns < maxDel) maxGap = maxDel

          var width = (maxGap + abs((rlen - queryLen) + 1)) >> 1
          if(width > w) width = w
          val minWidth = abs(rlen - queryLen) + 3
          if(width < minWidth) width = minWidth
          // NW alignment
          //score = ksw_global2(l_query, query, rlen, rseq, 5, mat, o_del, e_del, o_ins, e_ins, w, n_cigar, &cigar);
        }
       
        // compute NM and MD
        // str.l = str.m = *n_cigar * 4; str.s = (char*)cigar; // append MD to CIGAR   
        var int2base = ""
        var x = 0
        var y = 0
        var u = 0
        var n_mm = 0
        var nGap = 0

        if(rBeg < pacLen) int2base = "ACGTN"
        else int2base = "TGCAN"        

        for(k <- 0 to (numCigar - 1)) {
          // cigar = (uint32_t*)str.s;
          var op = cigar(k) & 0xf
          var len = cigar(k) >> 4
          
          // match
          if(op == 0) {

          }
        }
        
        (0, 0, 0, cigar)
      }
    }

  }
}

