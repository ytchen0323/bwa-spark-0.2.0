package cs.ucla.edu.bwaspark.worker2

import scala.util.control.Breaks._
import scala.List
import scala.collection.mutable.MutableList
import scala.math.log

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


  private def bwaFixXref2(mat: Array[Byte], oDel: Int, eDel: Int, oIns: Int, eIns: Int, w: Int, bns: BNTSeqType, 
    pac: Array[Byte], query: Array[Byte], qBeg: Int, qEnd: Int, rBeg: Long, rEnd: Long): Array[Int] = {
    var retArray = new Array[Int](5)
    bwaGenCigar2(mat, oDel, eDel, oIns, eIns, w, bns.l_pac, pac, qEnd - qBeg, query, rBeg, rEnd)
    retArray
  }

  private def bwaGenCigar2(mat: Array[Byte], oDel: Int, eDel: Int, oIns: Int, eIns: Int, w: Int, pacLen: Long, pac: Array[Byte], 
    queryLen: Int, query: Array[Byte], rb: Long, re: Long): (Int, Int, Int, Array[Int]) = {
    var cigar: Array[Int] = new Array[Int](10)    

    (0, 0, 0, cigar)
  }
}

