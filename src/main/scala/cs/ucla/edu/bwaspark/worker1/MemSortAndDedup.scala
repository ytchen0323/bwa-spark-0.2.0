package cs.ucla.edu.bwaspark.worker1

import scala.util.control.Breaks._
import scala.List
import scala.collection.mutable.MutableList

import cs.ucla.edu.bwaspark.datatype._

object MemSortAndDedup {
  /**
    *  Sort the MemAlnRegs according to the given order
    *  and remove the redundant MemAlnRegs
    *  
    *  @param regsIn alignment registers, which are the output of chain to alignment (after memChainToAln() is applied)
    *  @param maskLevelRedun mask level of redundant alignment registers (from MemOptType object)
    */
  def memSortAndDedup(regsIn: MutableList[MemAlnRegType], maskLevelRedun: Float): MutableList[MemAlnRegType] = {
    if(regsIn.length <= 1) {
      regsIn
    } 
    else {
      var regs = regsIn.sortBy(_.rEnd)
      
      for(i <- 1 to (regs.length - 1)) {
        if(regs(i).rBeg < regs(i-1).rEnd) {
          var j = i - 1

          breakable {
            while(j >= 0 && regs(i).rBeg < regs(j).rEnd) {
              // a[j] has been excluded
              if(regs(j).qEnd != regs(j).qBeg) { 
                var oq = 0
                var mr: Long = 0
                var mq = 0
                var or = regs(j).rEnd - regs(i).rBeg // overlap length on the reference
                // overlap length on the query
                if(regs(j).qBeg < regs(i).qBeg) oq = regs(j).qEnd - regs(i).qBeg
                else oq = regs(i).qEnd - regs(j).qBeg
                // min ref len in alignment
                if(regs(j).rEnd - regs(j).rBeg < regs(i).rEnd - regs(i).rBeg) mr = regs(j).rEnd - regs(j).rBeg
                else mr = regs(i).rEnd - regs(i).rBeg
                // min qry len in alignment
                if(regs(j).qEnd - regs(j).qBeg < regs(i).qEnd - regs(i).qBeg) mq = regs(j).qEnd - regs(j).qBeg
                else mq = regs(i).qEnd - regs(i).qBeg
                // one of the hits is redundant
                if(or.toFloat > maskLevelRedun * mr && oq.toFloat > maskLevelRedun * mq) {
                  if(regs(i).score < regs(j).score) {
                    regs(i).qEnd = regs(i).qBeg
                    break
                  }
                  else
                    regs(j).qEnd = regs(j).qBeg
                }
              }             
 
              j -= 1
            }
          }

        }
      }

      // exclude identical hits
      regs = regs.filter(r => (r.qEnd > r.qBeg))
      regs = regs.sortBy(r => (- r.score, r.rBeg, r.qBeg))
      
      for(i <- 1 to (regs.length - 1))
        if(regs(i).score == regs(i-1).score && regs(i).rBeg == regs(i-1).rBeg && regs(i).qBeg == regs(i-1).qBeg)
          regs(i).qEnd = regs(i).qBeg
        
      regs.filter(r => (r.qEnd > r.qBeg))
    }
  }
}

