package cs.ucla.edu.bwaspark.util

import scala.util.control.Breaks._
import scala.math.abs

object SWUtil {
  /**
    *  EHType: use for Smith Waterman Extension
    */
  class EHType(e_i: Int, h_i: Int) {
    var e: Int = e_i
    var h: Int = h_i
  }


  /**
    *  Smith-Waterman Extension
    *  The key function in both left and right extension
    *
    *  @param qLen the query length
    *  @param query the query (read)
    *  @param tLen the target (reference) length retrieved
    *  @param target the target (reference) retrieved
    *  @param m the (mat) array size in one dimension 
    *  @param mat the mat array
    *  @param oDel oDel in the input MemOptType object
    *  @param eDel eDel in the input MemOptType object
    *  @param oIns oIns in the input MemOptType object
    *  @param eIns eIns in the input MemOptType object
    *  @param w_i the input weight
    *  @param endBonus endBonus in the input MemOptType object
    *  @param zdrop zdrop in the input MemOptType object
    *  @param h0 initial S-W score
    */  
  def SWExtend(
    qLen: Int, query: Array[Byte], tLen: Int, target: Array[Byte], m: Int, mat: Array[Byte],
    oDel: Int, eDel: Int, oIns: Int, eIns: Int, w_i: Int, endBonus: Int, zdrop: Int, h0: Int): Array[Int] =  
  {
    var retArray: Array[Int] = new Array[Int](6)
    var eh: Array[EHType] = new Array[EHType](qLen + 1) // score array
    var qp: Array[Byte] = new Array[Byte](qLen * m) // query profile
    var oeDel = oDel + eDel
    var oeIns = oIns + eIns
    var i = 0 
    var j = 0 
    var k = 0
    var w = w_i

    for(i <- 0 until (qLen + 1))
      eh(i) = new EHType(0, 0)

    // generate the query profile
    i = 0
    for(k <- 0 to (m - 1)) {
      val p = k * m
      
      for(j <- 0 to (qLen - 1)) {
        qp(i) = mat(p + query(j))
        i += 1
      }
    }
    
    // fill the first row
    eh(0).h = h0
    if(h0 > oeIns) eh(1).h = h0 - oeIns
    else eh(1).h = 0
    j = 2
    while(j <= qLen && eh(j-1).h > eIns) {
      eh(j).h = eh(j-1).h - eIns
      j += 1
    }

    // adjust $w if it is too large
    k = m * m
    var max = 0
    for(i <- 0 to (k - 1)) // get the max score
      if(max < mat(i))
        max = mat(i)
    var maxIns = ((qLen * max + endBonus - oIns).toDouble / eIns + 1.0).toInt
    if(maxIns < 1) maxIns = 1
    var maxDel = ((qLen * max + endBonus - oDel).toDouble / eDel + 1.0).toInt
    if(maxDel < 1) maxDel = 1
    if(w > maxDel) w = maxDel  // TODO: is this necessary? (in original C implementation)

    // DP loop
    max = h0
    var max_i = -1
    var max_j = -1
    var max_ie = -1
    var gscore = -1
    var max_off = 0
    var beg = 0
    var end = qLen

    breakable {
      for(i <- 0 to (tLen - 1)) {
        var t = 0
        var f = 0
        var h1 = 0
        var m = 0
        var mj = -1
        var qPtr = target(i) * qLen
        // compute the first column
        h1 = h0 - (oDel + eDel * (i + 1))
        if(h1 < 0) h1 = 0
        // apply the band and the constraint (if provided)
        if (beg < i - w) beg = i - w
        if (end > i + w + 1) end = i + w + 1
        if (end > qLen) end = qLen

      
        for(j <- beg to (end - 1)) {
          // At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
          // Similar to SSE2-SW, cells are computed in the following order:
          //   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
          //   E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
          //   F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape
          var h = eh(j).h
          var e = eh(j).e   // get H(i-1,j-1) and E(i-1,j)
          eh(j).h = h1
          h += qp(qPtr + j)
          if(h < e) h = e
          if(h < f) h = f 
          h1 = h            // save H(i,j) to h1 for the next column
          if(m < h) { 
            mj = j          // record the position where max score is achieved
            m = h           // m is stored at eh[mj+1]
          }
          t = h - oeDel
          if(t < 0) t = 0
          e -= eDel
          if(e < t) e = t   // computed E(i+1,j)
          eh(j).e = e       // save E(i+1,j) for the next row
          t = h - oeIns
          if(t < 0) t = 0
          f -= eIns
          if(f < t) f = t
        }
      
        eh(end).h = h1
        eh(end).e = 0
        // end == j after the previous loop
        if(end == qLen) {
          if(gscore < h1) {
            max_ie = i
            gscore = h1
          }
        }

        if(m == 0) break

        if(m > max) {
          max = m
          max_i = i
          max_j = mj

          if(max_off < abs(mj - i)) max_off = abs(mj - i)
        }
        else if(zdrop > 0) {
          if((i - max_i) > (mj - max_j)) 
            if(max - m - ((i - max_i) - (mj - max_j)) * eDel > zdrop) break
          else
            if(max - m - ((mj - max_j) - (i - max_i)) * eIns > zdrop) break;
        }
        
        // update beg and end for the next round
        j = mj
        while(j >= beg && eh(j).h > 0) {
          j -= 1
        }
        beg = j + 1

        j = mj + 2
        while(j <= end && eh(j).h > 0) {
          j += 1
        }
        end = j
      }

    }

    retArray(0) = max
    retArray(1) = max_j + 1
    retArray(2) = max_i + 1
    retArray(3) = max_ie + 1
    retArray(4) = gscore
    retArray(5) = max_off
    
    retArray
  }

}
