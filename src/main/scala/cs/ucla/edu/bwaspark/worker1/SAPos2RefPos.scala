package cs.ucla.edu.bwaspark.worker1

import cs.ucla.edu.bwaspark.datatype._

//standalone object for transforming from position in suffix array
//to position in reference genome
object SAPos2RefPos {

  def occAux(y: Long, c: Int): Int = {

    var res1 = 0l
    var res2 = 0l

    var tmp1 = 0l
    var tmp2 = 0l
    var tmp3 = 0l

    if ( (c & 2) != 0 ) tmp1 = y else tmp1 = ~y

    tmp1 = tmp1 >>> 1

    if ( (c & 1) != 0 ) tmp2 = y else tmp2 = ~y

    tmp3 = 0x5555555555555555l

    res1 = tmp1 & tmp2 & tmp3

    tmp3 = 0x3333333333333333l

    tmp1 = (res1 & tmp3)

    tmp2 = (res1 >>> 2) & tmp3

    res2 = tmp1 + tmp2

    tmp3 = 0x0f0f0f0f0f0f0f0fl

    tmp1 = (res2 + (res2 >>> 4)) & tmp3

    tmp1 = tmp1 * 0x0101010101010101l

    (tmp1 >>> 56).toInt
  }



  def bwtOcc(bwt: BWTType, pos: Long, x: Long): Long = {
    var n: Long = 0l
    var k = pos

    //transform c from Long to ubyte_t?(Int instead)
    var c = (x & 0xffl).toInt

    if (k == bwt.seqLen) bwt.L2(c+1) - bwt.L2(c)
    else if (k == -1l) 0l //the original is (uint64_t)(-1)
    else {
      if (k >= bwt.primary) k = k - 1

      //calculate new pointer position
      var newStartPoint = (k >>> 7) << 4

      n = bwt.bwt(c)

      //jump to the start of the first bwt cell

      newStartPoint += 8 //size of Long: 8 bytes

      val occIntvMask = (1l << 7) - 1l

      var newEndPoint = newStartPoint + (((k >>> 5) - ((k & ~occIntvMask) >>> 5)) << 1)
      while (newStartPoint < newEndPoint) {
        n = n + occAux(bwt.bwt(newStartPoint.toInt) << 32 | bwt.bwt(newStartPoint.toInt + 1), c)
        newStartPoint += 2
      }
      n += occAux((bwt.bwt(newStartPoint.toInt) << 32 | bwt.bwt(newStartPoint.toInt + 1)) & ~((1l << ((~k & 31) << 1)) - 1), c)
      if (c == 0) n -= ~k & 31

      n

    }
  }

  //compute inverse CSA
  def bwtInvPsi(bwt:BWTType, k: Long): Long = {
    var x: Long = if (k > bwt.primary) k - 1 else k

    var bwtBwt = bwt.bwt(((k >>> 7 << 4) + 8 + ((k & 0x7f) >>> 4)).toInt)
    x = bwtBwt >> (((~k & 0xf) << 1) & 3)
    x = bwt.L2(x.toInt) + bwtOcc(bwt, k, x)
    if (k == bwt.primary) 0
    else x
  }

  def suffixArrayPos2ReferencePos(bwt: BWTType, k: Long /*uint64_t*/): Long /*uint64_t*/ = {
 
    //initialization
    var sa = 0l
    var mask = (bwt.saIntv - 1).toLong

    var pos = k

    while ( (pos & mask) != 0 ) {
      sa += 1l
      pos = bwtInvPsi(bwt, pos)
    }

    sa + bwt.sa((pos / bwt.saIntv).toInt)

  }
}
