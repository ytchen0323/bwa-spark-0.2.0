package cs.ucla.edu.bwaspark.worker1

import scala.util.control.Breaks._
import scala.List
import scala.collection.mutable.MutableList

import cs.ucla.edu.bwaspark.datatype._

object MemRegToADAMSAM {
  val MEM_F_ALL = 0x8
  /**
    *  Transform the alignment registers to SAM format
    *  
    *  @param opt the input MemOptType object
    *  @param bns the input BNSSeqType object
    *  @param pac the PAC array
    *  @param seq the read (NOTE: currently we use Array[Byte] first. may need to be changed!!!)
    *  @param regsIn the alignment registers to be transformed
    *  @param extraFlag
    *  @param alns 
    */
  def memRegToSAMSe(opt: MemOptType, bns: BNTSeqType, pac: Array[Byte], seq: Array[Byte], regsIn: MutableList[MemAlnRegType], extraFlag: Int, alns: List[MemAlnReg]) {
    var regs = regsIn.filter(_.score >= opt.T).filter(r => (r.secondary < 0) || ((opt.flag & MEM_F_ALL) > 0))
  }
}

