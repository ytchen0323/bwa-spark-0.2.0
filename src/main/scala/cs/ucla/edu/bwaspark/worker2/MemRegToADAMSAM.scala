package cs.ucla.edu.bwaspark.worker2

import scala.collection.mutable.MutableList
import scala.math.log
import scala.math.abs

import cs.ucla.edu.bwaspark.datatype._
import cs.ucla.edu.bwaspark.sam.SAMHeader._
import cs.ucla.edu.bwaspark.util.BNTSeqUtil._
import cs.ucla.edu.bwaspark.util.SWUtil._

object MemRegToADAMSAM {
  val MEM_F_ALL = 0x8
  val MEM_F_NO_MULTI = 0x10
  val MEM_MAPQ_COEF = 30.0
  val int2op = Array('M', 'I', 'D', 'S', 'H')
  val int2forward = Array('A', 'C', 'G', 'T', 'N')
  val int2reverse = Array('T', 'G', 'C', 'A', 'N')
 

  class SAMString {
    var str: Array[Char] = new Array[Char](8192)
    var idx: Int = 0
    var size: Int = 8192

    def addCharArray(in: Array[Char]) {
      if((idx + in.size + 1) >= size) {
        size = size << 2
        val old = str
        str = new Array[Char](size)
        old.copyToArray(str, 0, idx + 1)
      }
 
      var i = 0
      while(i < in.size) {
        str(idx) = in(i)
        i += 1
        idx += 1
      }
    }
 
    def addChar(c: Char) {
      if((idx + 1) >= size) {
        size = size << 2
        val old = str
        str = new Array[Char](size)
        old.copyToArray(str, 0, idx + 1)
      }

      str(idx) = c
      idx += 1
    }
  }

  /**
    *  Transform the alignment registers to SAM format
    *  
    *  @param opt the input MemOptType object
    *  @param bns the input BNTSeqType object
    *  @param pac the PAC array
    *  @param seq the read (NOTE: currently we use Array[Byte] first. may need to be changed back to Array[Char]!!!)
    *  @param regs the alignment registers to be transformed
    *  @param extraFlag
    *  @param alnIn currently we skip this parameter
    */
  def memRegToSAMSe(opt: MemOptType, bns: BNTSeqType, pac: Array[Byte], seq: FASTQSingleNode, regs: Array[MemAlnRegType], extraFlag: Int, alnMate: MemAlnType) {
    var alns: MutableList[MemAlnType] = new MutableList[MemAlnType]

    // NOTE: set opt.flag manually here!!! This should be modified from the logger!!!
    opt.flag = 24

    //pre-process: transform A/C/G/T to 0,1,2,3
    def locusEncode(locus: Char): Byte = {
      //transforming from A/C/G/T to 0,1,2,3
      locus match {
        case 'A' => 0
        case 'a' => 0
        case 'C' => 1
        case 'c' => 1
        case 'G' => 2
        case 'g' => 2
        case 'T' => 3
        case 't' => 3
        case '-' => 5
        case _ => 4
      }
    }

    val seqTrans = seq.seq.toCharArray.map(ele => locusEncode(ele))

/*   
    println("[opt object] flag: " + opt.flag + " T: " + opt.T + " minSeedLen: " + opt.minSeedLen + " a: " + opt.a + " b: " + opt.b + " mapQCoefLen: " + opt.mapQCoefLen + " mapQCoefFac: " + opt.mapQCoefFac)
*/
/*
    var j = 0
    regs.foreach(r => {
      print("Reg " + j + "(")
      print(r.rBeg + ", " + r.rEnd + ", " + r.qBeg + ", " + r.qEnd + ", " + r.score + ", " + r.trueScore + ", ")
      println(r.sub + ", "  + r.csub + ", " + r.subNum + ", " + r.width + ", " + r.seedCov + ", " + r.secondary + ")")
      j += 1
      } )
*/
    if(regs != null) {
      var i = 0
      while(i < regs.length) {
        if(regs(i).score >= opt.T) {
          if(regs(i).secondary < 0 || ((opt.flag & MEM_F_ALL) > 0)) {
            if(regs(i).secondary < 0 || regs(i).score >= regs(regs(i).secondary).score * 0.5) {
              // debugging
              //print("Aln " + i + " " +  regs(i).score + " ")
              var aln = memRegToAln(opt, bns, pac, seq.seqLen, seqTrans, regs(i))   // NOTE: current data structure has not been obtained from RDD. We assume the length to be 101 here
              alns += aln
              aln.flag |= extraFlag   // flag secondary
              if(regs(i).secondary >= 0) aln.sub = -1   // don't output sub-optimal score
              if(i > 0 && regs(i).secondary < 0)   // if supplementary
                if((opt.flag & MEM_F_NO_MULTI) > 0) aln.flag |= 0x10000
                else aln.flag |= 0x800

              if(i > 0 && aln.mapq > alns.head.mapq) aln.mapq = alns.head.mapq            
            }
          }
        }

        i += 1
      }
    }

    //seqTrans.foreach(print(_))
    //println

    // no alignments good enough; then write an unaligned record

    var samStr = new SAMString

    if(alns.length == 0) { 
      var aln = memRegToAln(opt, bns, pac, seq.seqLen, seqTrans, null)
      aln.flag |= extraFlag
      var alnList = new Array[MemAlnType](1)
      alnList(0) = aln

      memAlnToSAM(bns, seq, seqTrans, alnList, 0, alnMate, samStr)
    }
    else {
      var k = 0
      val alnsArray = alns.toArray
      while(k < alns.size) {
        memAlnToSAM(bns, seq, seqTrans, alnsArray, k, alnMate, samStr)
        k += 1
      }
    }

    seq.sam = samStr.str.dropRight(samStr.size - samStr.idx).mkString
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

   //def inferBw( l1: Int, l2: Int, score: Int, a: Int, q: Int, r: Int) : Int = {
   private def inferBw( l1: Int, l2: Int, score: Int, a: Int, q: Int, r: Int) : Int = {
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

   /*
   get_rLen actually this function is barely called in the current flow
   @param n_cigar
   @param *cigar 
   @param l output


   */
  // this one is not tested because it is never been used in 20 reads dataset

  private def getRlen(cigarSegs: Vector[CigarSegType]) : Int = {
    var l: Int = 0

    var k = 0
    while(k < cigarSegs.size) {
      if(cigarSegs(k).op == 0 || cigarSegs(k).op == 2) l += cigarSegs(k).len

      k += 1
    }
    l
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

    if(reg == null || reg.rBeg < 0 || reg.rEnd < 0) {
      aln.rid = -1
      aln.pos = -1
      aln.flag |= 0x4
      aln
    }
    else {
      val qb = reg.qBeg
      val qe = reg.qEnd
      val rb = reg.rBeg
      val re = reg.rEnd

      if(reg.secondary < 0) 
        aln.mapq = memApproxMapqSe(opt, reg).toShort
      else
        aln.mapq = 0

      // secondary alignment
      if(reg.secondary >= 0) aln.flag |= 0x100 

      val ret = bwaFixXref2(opt.mat, opt.oDel, opt.eDel, opt.oIns, opt.eIns, opt.w, bns, pac, seq, reg.qBeg, reg.qEnd, reg.rBeg, reg.rEnd)
      val iden = ret._5
      if(iden < 0) {
        println("[Error] If you see this message, please let the developer know. Abort. Sorry.")
        assert(false, "bwaFixXref2() problem")
      }

      var tmp = inferBw(qe - qb, (re - rb).toInt, reg.trueScore, opt.a, opt.oDel, opt.eDel)
      var w2 = inferBw(qe - qb, (re - rb).toInt, reg.trueScore, opt.a, opt.oIns, opt.eIns)
      if(w2 < tmp) w2 = tmp
      if(w2 > opt.w) {
        if(w2 > reg.width) w2 = reg.width
      }
      //else w2 = opt.w  // TODO: check if we need this line on long reads. On 1-800bp reads, it does not matter and it should be. (In original C implementation)

      var i = 0
      aln.cigar = null
      var score = 0
      var lastScore = -(1 << 30)

      var isBreak = false
      do {
        // make a copy to pass into bwaGenCigar2
        // if there is performance issue later, we may modify the underlying implementation
        var query: Array[Byte] = new Array[Byte](qe - qb)

        var j = 0
        while(j < (qe - qb)) {
          query(j) = seq(qb + j)
          j += 1
        }
 
        var ret = bwaGenCigar2(opt.mat, opt.oDel, opt.eDel, opt.oIns, opt.eIns, w2, bns.l_pac, pac, qe - qb, query, rb, re)
        score = ret._1
        aln.nCigar = ret._2
        aln.NM = ret._3
        aln.cigar = ret._4

        if(score == lastScore) isBreak = true
    
        if(!isBreak) {
          lastScore = score
          w2 <<= 1
        }

        i += 1
      } while(i < 3 && score < reg.trueScore - opt.a && !isBreak)


      var pos: Long = 0
      var isRev: Int = 0
 
      if(rb < bns.l_pac) {
        val ret = bnsDepos(bns, rb)
        pos = ret._1
        isRev = ret._2
      }
      else {
        val ret = bnsDepos(bns, re - 1)
        pos = ret._1
        isRev = ret._2
      }

      aln.isRev = isRev.toByte

      // squeeze out leading or trailing deletions
      if(aln.nCigar > 0) {
        if(aln.cigar.cigarSegs(0).op == 2) {
          pos += aln.cigar.cigarSegs(0).len
          //aln.cigar.cigarSegs = aln.cigar.cigarSegs.drop(1)
          aln.cigar.cigarSegs = aln.cigar.cigarSegs.tail
          aln.nCigar -= 1
        }
        else if(aln.cigar.cigarSegs(aln.nCigar - 1).op == 2) {            
          aln.cigar.cigarSegs = aln.cigar.cigarSegs.dropRight(1)
          aln.nCigar -= 1
        }
      }

      // add clipping to CIGAR
      if(qb != 0 || qe != seqLen) {
        var clip5 = 0
        var clip3 = 0

        if(isRev > 0) clip5 = seqLen - qe
        else clip5 = qb

        if(isRev > 0) clip3 = qb
        else clip3 = seqLen - qe

        if(clip5 > 0) {
          val cigarSeg = new CigarSegType
          cigarSeg.op = 3
          cigarSeg.len = clip5
          aln.cigar.cigarSegs = aln.cigar.cigarSegs.+:(cigarSeg)  // prepend (Vector type)
          aln.nCigar += 1
        }
        
        if(clip3 > 0) {
          val cigarSeg = new CigarSegType
          cigarSeg.op = 3
          cigarSeg.len = clip3
          aln.cigar.cigarSegs = aln.cigar.cigarSegs :+ cigarSeg  // append (Vector type)
          aln.nCigar += 1
        }
      }

      aln.rid = bnsPosToRid(bns, pos)
      aln.pos = pos - bns.anns(aln.rid).offset
      aln.score = reg.score
      if(reg.sub > reg.csub) aln.sub = reg.sub
      else aln.sub = reg.csub
      
      aln
    }
  }


  //private def memAlnToSAM(bns: BNTSeqType, seq: FASTQSingleNode, seqTrans: Array[Byte], alnList: Array[MemAlnType], which: Int, alnMate: MemAlnType): String = {
  private def memAlnToSAM(bns: BNTSeqType, seq: FASTQSingleNode, seqTrans: Array[Byte], alnList: Array[MemAlnType], which: Int, alnMate: MemAlnType, samStr: SAMString) {
    var aln = alnList(which)
    var alnTmp = aln.copy
    var alnMateTmp: MemAlnType = null
    if(alnMate != null) alnMateTmp = alnMate.copy 

    // set flag
    if(alnMateTmp != null) alnTmp.flag |= 0x1 // is paired in sequencing
    if(alnTmp.rid < 0) alnTmp.flag |= 0x4 // is mapped
    if(alnMateTmp != null && alnMateTmp.rid < 0) alnTmp.flag |= 0x8 // is mate mapped
    if(alnTmp.rid < 0 && alnMateTmp != null && alnMateTmp.rid >= 0) { // copy mate to alignment
      alnTmp.rid = alnMateTmp.rid
      alnTmp.pos = alnMateTmp.pos
      alnTmp.isRev = alnMateTmp.isRev
      alnTmp.nCigar = 0
    }
    if(alnMateTmp != null && alnMateTmp.rid < 0 && alnTmp.rid >= 0) { // copy alignment to mate
      alnMateTmp.rid = alnTmp.rid
      alnMateTmp.pos = alnTmp.pos
      alnMateTmp.isRev = alnTmp.isRev
      alnMateTmp.nCigar = 0
    }
    if(alnTmp.isRev > 0) alnTmp.flag |= 0x10 // is on the reverse strand
    if(alnMateTmp != null && alnMateTmp.isRev > 0) alnTmp.flag |= 0x20 // is mate on the reverse strand
       
    // print up to CIGAR
    samStr.addCharArray(seq.name.toCharArray)   // QNAME
    samStr.addChar('\t')
    if((alnTmp.flag & 0x10000) > 0) alnTmp.flag = (alnTmp.flag & 0xffff) | 0x100   // FLAG
    else alnTmp.flag = (alnTmp.flag & 0xffff) | 0
    samStr.addCharArray(alnTmp.flag.toString.toCharArray)
    samStr.addChar('\t')
    if(alnTmp.rid >= 0) { // with coordinate
      samStr.addCharArray(bns.anns(alnTmp.rid).name.toCharArray)   // RNAME
      samStr.addChar('\t')
      samStr.addCharArray((alnTmp.pos + 1).toString.toCharArray)   // POS
      samStr.addChar('\t')
      samStr.addCharArray(alnTmp.mapq.toString.toCharArray)   // MAPQ
      samStr.addChar('\t')

      if(alnTmp.nCigar > 0) {   // aligned
        var i = 0
        while(i < alnTmp.nCigar) {
          var c = alnTmp.cigar.cigarSegs(i).op
          if(c == 3 || c == 4) 
            if(which > 0) c = 4   // use hard clipping for supplementary alignments
            else c = 3
          samStr.addCharArray(alnTmp.cigar.cigarSegs(i).len.toString.toCharArray)
          samStr.addChar(int2op(c))
          i += 1
        }
      }
      else samStr.addChar('*')
    }
    else samStr.addCharArray("*\t0\t0\t*".toCharArray)   // without coordinte
    samStr.addChar('\t')

    // print the mate position if applicable
    if(alnMateTmp != null && alnMateTmp.rid >= 0) {
      if(alnTmp.rid == alnMateTmp.rid) samStr.addChar('=')
      else samStr.addCharArray(bns.anns(alnMateTmp.rid).name.toString.toCharArray)
      samStr.addChar('\t')
      samStr.addCharArray((alnMateTmp.pos + 1).toString.toCharArray)
      samStr.addChar('\t')
      if(alnTmp.rid == alnMateTmp.rid) {
        var p0: Long = -1
        var p1: Long = -1
        if(alnTmp.isRev > 0) p0 = alnTmp.pos + getRlen(alnTmp.cigar.cigarSegs) - 1
        else p0 = alnTmp.pos
        if(alnMateTmp.isRev > 0) p1 = alnMateTmp.pos + getRlen(alnMateTmp.cigar.cigarSegs) - 1
        else p1 = alnMateTmp.pos
        if(alnMateTmp.nCigar == 0 || alnTmp.nCigar == 0) samStr.addChar('0')
        else {
          if(p0 > p1) samStr.addCharArray((-(p0 - p1 + 1)).toString.toCharArray)
          else if(p0 < p1) samStr.addCharArray((-(p0 - p1 - 1)).toString.toCharArray)
          else samStr.addCharArray((-(p0 - p1)).toString.toCharArray)
        }
      }
      else samStr.addChar('0')
    }
    else samStr.addCharArray("*\t0\t0".toCharArray)
    samStr.addChar('\t')
    
    // print SEQ and QUAL
    if((alnTmp.flag & 0x100) > 0) {   // for secondary alignments, don't write SEQ and QUAL
      samStr.addCharArray("*\t*".toCharArray)
    }
    else if(alnTmp.isRev == 0) {   // the forward strand
      //println("Forward")
      var qb = 0
      var qe = seq.seqLen

      if(alnTmp.nCigar > 0) {
        if(which > 0 && (alnTmp.cigar.cigarSegs(0).op == 4 || alnTmp.cigar.cigarSegs(0).op == 3)) qb += alnTmp.cigar.cigarSegs(0).len
        if(which > 0 && (alnTmp.cigar.cigarSegs(alnTmp.nCigar - 1).op == 4 || alnTmp.cigar.cigarSegs(alnTmp.nCigar - 1).op == 3)) qe -= alnTmp.cigar.cigarSegs(alnTmp.nCigar - 1).len
      }

      var i = qb
      while(i < qe) {
        samStr.addChar(int2forward(seqTrans(i)))
        i += 1
      }
      samStr.addChar('\t')

      if(seq.qual != "") {
        var i = qb
        val seqArray = seq.qual.toCharArray
        while(i < qe) {
          samStr.addChar(seqArray(i))
          i += 1
        }
      }
      else samStr.addChar('*')
    }
    else {   // the reverse strand
      //println("Reverse")
      var qb = 0
      var qe = seq.seqLen
      
      if(alnTmp.nCigar > 0) {
        if(which > 0 && (alnTmp.cigar.cigarSegs(0).op == 4 || alnTmp.cigar.cigarSegs(0).op == 3)) qe -= alnTmp.cigar.cigarSegs(0).len
        if(which > 0 && (alnTmp.cigar.cigarSegs(alnTmp.nCigar - 1).op == 4 || alnTmp.cigar.cigarSegs(alnTmp.nCigar - 1).op == 3)) qb += alnTmp.cigar.cigarSegs(alnTmp.nCigar - 1).len
      }

      var i = qe - 1
      while(i >= qb) {
        samStr.addChar(int2reverse(seqTrans(i)))
        i -= 1
      }
      samStr.addChar('\t')

      if(seq.qual != "") {
        var i = qe - 1
        val seqArray = seq.qual.toCharArray
        while(i >= qb) {
          samStr.addChar(seqArray(i))
          i -= 1
        }
      }
      else samStr.addChar('*')
    }

    // print optional tags
    if(alnTmp.nCigar > 0) {
      samStr.addCharArray("\tNM:i:".toCharArray)
      samStr.addCharArray(alnTmp.NM.toString.toCharArray)
      samStr.addCharArray("\tMD:Z:".toCharArray)
      samStr.addCharArray(alnTmp.cigar.cigarStr.toCharArray)
    }
    if(alnTmp.score >= 0) {
      samStr.addCharArray("\tAS:i:".toCharArray)
      samStr.addCharArray(alnTmp.score.toString.toCharArray)
    }
    if(alnTmp.sub >= 0) {
      samStr.addCharArray("\tXS:i:".toCharArray)
      samStr.addCharArray(alnTmp.sub.toString.toCharArray)
    }
    // Read group is read using SAMHeader class 
    if(bwaReadGroupID != "") {
      samStr.addCharArray("\tRG:Z:".toCharArray)
      samStr.addCharArray(bwaReadGroupID.toCharArray)
    }
    
    if((alnTmp.flag & 0x100) == 0) { // not multi-hit
      var i = 0
      var isBreak = false
      while(i < alnList.size && !isBreak) {
        if(i != which && (alnList(i).flag & 0x100) == 0) { 
          isBreak = true 
          i -= 1
        }
        i += 1
      }

      if(i < alnList.size) { // there are other primary hits; output them
        samStr.addCharArray("\tSA:Z:".toCharArray)
        var j = 0
        while(j < alnList.size) {
          if(j != which && (alnList(j).flag & 0x100) == 0) { // proceed if: 1) different from the current; 2) not shadowed multi hit
            samStr.addCharArray(bns.anns(alnList(j).rid).name.toCharArray)
            samStr.addChar(',')
            samStr.addCharArray((alnList(j).pos + 1).toString.toCharArray)
            samStr.addChar(',')
            if(alnList(j).isRev == 0) samStr.addChar('+')
            else samStr.addChar('-')
            samStr.addChar(',')
            
            var k = 0
            while(k < alnList(j).nCigar) {
              samStr.addCharArray(alnList(j).cigar.cigarSegs(k).len.toString.toCharArray)
              samStr.addChar(int2op(alnList(j).cigar.cigarSegs(k).op))
              k += 1
            }

            samStr.addChar(',')
            samStr.addCharArray(alnList(j).mapq.toString.toCharArray)
            samStr.addChar(',')
            samStr.addCharArray(alnList(j).NM.toString.toCharArray)
            samStr.addChar(';')

          }
          j += 1
        }
      } 
    }

    if(seq.comment != "") {
      samStr.addChar('\t')
      samStr.addCharArray(seq.comment.toCharArray)
    }
    samStr.addChar('\n')

    //samStr.str.dropRight(samStr.size - samStr.idx).mkString
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
    pac: Array[Byte], query: Array[Byte], qBeg: Int, qEnd: Int, rBeg: Long, rEnd: Long): (Int, Int, Long, Long, Int) = {
    var retArray = new Array[Int](5)

    // cross the for-rev boundary; actually with BWA-MEM, we should never come to here
    if(rBeg < bns.l_pac && rEnd > bns.l_pac) {
      println("[ERROR] In BWA-mem, we should never come to here")
      (-1, -1, -1, -1, -1)  // unable to fix
    }
    else {
      val ret = bnsDepos(bns, (rBeg + rEnd) >> 1)  // coordinate of the middle point on the forward strand
      val fm = ret._1
      val isRev = ret._2
      val ra = bns.anns(bnsPosToRid(bns, fm))  // annotation of chr corresponding to the middle point
      var cb = ra.offset
      if(isRev > 0) cb = (bns.l_pac << 1) - (ra.offset + ra.len)
      var ce = cb + ra.len   // chr end
      var qBegRet = qBeg
      var qEndRet = qEnd
      var rBegRet = rBeg
      var rEndRet = rEnd

      // fix is needed
      if(cb > rBeg || ce < rEnd) {
        if(cb < rBeg) cb = rBeg
        if(ce > rEnd) ce = rEnd

        var queryArr: Array[Byte] = new Array[Byte](qEnd - qBeg)
        // make a copy to pass into bwaGenCigar2 
        // if there is performance issue later, we may modify the underlying implementation
        var i = 0
        while(i < (qEnd - qBeg)) {
          queryArr(i) = query(qBeg + i)
          i += 1
        }

        val ret = bwaGenCigar2(mat, oDel, eDel, oIns, eIns, w, bns.l_pac, pac, qEnd - qBeg, queryArr, rBeg, rEnd)
        val numCigar = ret._2
        val cigar = ret._4

        var x = rBeg
        var y = qBeg

        var isBreak = false
        i = 0
        while(i < numCigar && !isBreak) {
          val op = cigar.cigarSegs(i).op
          val len = cigar.cigarSegs(i).len

          if(op == 0) {
            if(x <= cb && cb < x + len) {
              qBegRet = (y + (cb - x)).toInt
              rBegRet = cb
            }
              
            if(x < ce && ce <= x + len) {
              qEndRet = (y + (ce - x)).toInt
              rEndRet = ce
              isBreak = true
            }
            else {
              x += len
              y += len
            }
          }
          else if(op == 1) {
            y += len
          } 
          else if(op == 2) {
            if(x <= cb && cb < x + len) {
              qBegRet = y
              rBegRet = x + len
            }
              
            if(x < ce && ce <= x + len) {
              qEndRet = y
              rEndRet = x
            }
            else x += len
          }
          else {
            println("[Error] Should not be here!!!")
            assert(false, "in bwaFixXref2()")
          }

          i += 1
        }
        // NOTE!!!: Need to be implemented!!! temporarily skip this for loop
      }
    
      var iden = 0
      if(qBegRet == qEndRet || rBegRet == rEndRet) iden = -2
      (qBegRet, qEndRet, rBegRet, rEndRet, iden)
    }

  }

  private def bwaGenCigar2(mat: Array[Byte], oDel: Int, eDel: Int, oIns: Int, eIns: Int, w: Int, pacLen: Long, pac: Array[Byte], 
    queryLen: Int, query_i: Array[Byte], rBeg: Long, rEnd: Long): (Int, Int, Int, CigarType) = {

    var numCigar = 0
    var NM = -1
    var score = 0
    var cigar = new CigarType

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
          var i = 0
          while(i < (queryLen >> 1)) {
            var tmp = query(i)
            query(i) = query(queryLen - 1 - i)
            query(queryLen - 1 - i) = tmp

            i += 1
          }
            
          i = 0
          while(i < (rlen >> 1).toInt) {
            var tmp = rseq(i)
            rseq(i) = rseq((rlen - 1 - i).toInt)
            rseq((rlen - 1 - i).toInt) = tmp

            i += 1
          }
        }        
        // no gap; no need to do DP
        if(queryLen == (rEnd - rBeg) && w == 0) {
          // FIXME: due to an issue in mem_reg2aln(), we never come to this block. This does not affect accuracy, but it hurts performance. (in original C implementation)
          //println("ENTER!!!")
          var cigarSeg = new CigarSegType
          cigarSeg.len = queryLen 
          cigarSeg.op = 0
          numCigar = 1
          cigar.cigarSegs = cigar.cigarSegs :+ cigarSeg   // append (Vector type)

          var i = 0
          while(i < queryLen) {
            score += mat(rseq(i) * 5 + query(i))
            i += 1
          }
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
          //val ret = SWGlobal(queryLen, query, rlen.toInt, rseq, 5, mat, oDel, eDel, oIns, eIns, width.toInt, numCigar, cigar.cigarSegs)
          val ret = SWGlobal(queryLen, query, rlen.toInt, rseq, 5, mat, oDel, eDel, oIns, eIns, width.toInt, numCigar, cigar)
          score = ret._1
          numCigar = ret._2
          //println("maxGap " + maxGap + "; numCigar " + numCigar)
        }
       
        // compute NM and MD
        // str.l = str.m = *n_cigar * 4; str.s = (char*)cigar; // append MD to CIGAR   
        var int2base = new Array[Char](5)
        var x = 0
        var y = 0
        var u = 0
        var n_mm = 0
        var nGap = 0

        if(rBeg < pacLen) int2base = Array('A', 'C', 'G', 'T', 'N')
        else int2base = Array('T', 'G', 'C', 'A', 'N')        

        //println("numCigar: " + numCigar)

        var k = 0
        while(k < numCigar) {
          val op = cigar.cigarSegs(k).op
          val len = cigar.cigarSegs(k).len
        
          //println("op " + op + ", len " + len)
  
          // match
          if(op == 0) {
            var i = 0
            while(i < len) {
              if(query(x + i) != rseq(y + i)) {
                cigar.cigarStr += u.toString
                cigar.cigarStr += int2base(rseq(y + i))
                n_mm += 1
                u = 0
              }
              else u += 1

              i += 1
            }

            x += len
            y += len
          }
          // deletion
          else if(op == 2) {
            // don't do the following if D is the first or the last CIGAR
            if(k > 0 && k < numCigar - 1) {
              cigar.cigarStr += u.toString
              cigar.cigarStr += '^'
              
              var i = 0
              while(i < len) {
                cigar.cigarStr += int2base(rseq(y + i))
                i += 1
              }
              
              u = 0
              nGap += len
            }

            y += len
          }
          // insertion
          else if(op == 1) {
            x += len
            nGap += len
          }

          k += 1
        }
        
        cigar.cigarStr += u.toString
        var NM = n_mm + nGap

        // reverse back query 
        // This is done in original C implementation. However, in the current implementation, the query is a copy
        // from the original array. Therefore, we do not need to reverse this copy back
        //if(rBeg >= pacLen) 
          //for(i <- 0 to ((queryLen >> 1) - 1)) {
            //var tmp = query(i)
            //query(i) = query(queryLen - 1 - i)
            //query(queryLen - 1 - i) = tmp
          //}
        
        //println(cigar.cigarStr)

        (score, numCigar, NM, cigar)
      }
    }

  }
}

