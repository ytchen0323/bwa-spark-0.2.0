package cs.ucla.edu.bwaspark.worker1

import cs.ucla.edu.bwaspark.datatype._
import scala.util.control.Breaks._
import scala.collection.mutable.MutableList
//import org.scalatest.Assertions._

class BWTSMem {

  //constant value in the orignal algo
  val OCC_INTV_SHIFT = 7  //bwt.h
  val OCC_INTERVAL = (1 << OCC_INTV_SHIFT) //bwt.h, assume 0x80
  val OCC_INTV_MASK = OCC_INTERVAL - 1

  //parameter/local variable in the original bwt_smem1
  var m_ik : BWTIntvType = _
  var m_ok : Array[BWTIntvType] = new Array[BWTIntvType](4)
  var mBWT : BWTType = _

  //local variable in the original bwt_extend
  var m_tk : Array[Long] = new Array[Long](4)
  var m_tl : Array[Long] = new Array[Long](4)

  def bwtOccIntv(k: Long) : Long = {
    (k >>> 7) << 4
  }

  def occAux4(b: Int) : Long = {
    mBWT.cntTable(b & 0xff) + mBWT.cntTable((b >>> 8) & 0xff) + mBWT.cntTable((b >>> 16) & 0xff) + mBWT.cntTable(b >>> 24)
  }

  def bwt_occ4(k: Long) : Array[Long] = {
    var cnt : Array[Long] = Array(0, 0, 0, 0);
    if (k == -1) return cnt

    var _k = k
    if (k >= mBWT.primary) _k -= 1
    var index = bwtOccIntv(_k)
    var tmp_idx = index.toInt
    assert(tmp_idx >= 0) //assertion enabled?
    //mimic memcpy in the orignal program, test needed!!!
    for (i <- 0 to 3) {
      cnt(i) = mBWT.bwt(tmp_idx) << 32 | mBWT.bwt(tmp_idx + 1)
      tmp_idx += 2
    }
    index += 8
    var end_idx = index + ((_k >>> 4) - ((_k & ~OCC_INTV_MASK) >>> 4))
    var x : Long = 0
    while (index < end_idx) {
      x += occAux4(mBWT.bwt(index.toInt))
      index += 1
    }

    var tmp : Int = mBWT.bwt(index.toInt) & ~((1 << ((~_k & 15) << 1)) - 1)
    x += occAux4(tmp) - (~_k & 15)
    cnt(0) += x & 0xff
    cnt(1) += (x >>> 8) & 0xff
    cnt(2) += (x >>> 16) & 0xff
    cnt(3) += (x >>> 24)
    cnt
  }

  def bwt_2occ4(k: Long, l: Long) {
    var _k = k
    var _l = l
    if (k >= mBWT.primary) _k -= 1
    if (l >= mBWT.primary) _l -= 1
    if (_l >>> OCC_INTV_SHIFT != _k >>> OCC_INTV_SHIFT || k == -1 || l == -1) {
      m_tk = bwt_occ4(k)
      m_tl = bwt_occ4(l)
    }
    else {
      if (k >= mBWT.primary) _k = k - 1
      if (l >= mBWT.primary) _l = l - 1
      var index = bwtOccIntv(_k)
      var tmp_idx = index.toInt
      assert(tmp_idx >= 0) //assertion enabled?
      //mimic memcpy in the orignal program, test needed!!!
      for (i <- 0 to 3) {
        m_tk(i) = mBWT.bwt(tmp_idx) << 32 | mBWT.bwt(tmp_idx + 1)
	tmp_idx += 2
      }
      index += 8 //sizeof(bwtint_t) in the original algo
      var endk_idx = index + ((_k >>> 4) - ((_k & ~OCC_INTV_MASK) >>> 4))
      var endl_idx = index + ((_l >>> 4) - ((_l & ~OCC_INTV_MASK) >>> 4))
      var x : Long = 0
      var y : Long = 0
      while (index < endk_idx) {
	x += occAux4(mBWT.bwt(index.toInt))
	index += 1
      }
      y = x
      var tmp : Int = 0
      tmp = mBWT.bwt(index.toInt) & ~((1 << ((~_k & 15) << 1)) - 1)
      x += occAux4(tmp) - (~_k & 15)
      while (index < endl_idx) {
	y += occAux4(mBWT.bwt(index.toInt))
	index += 1
      }
      tmp = mBWT.bwt(index.toInt) & ~((1 << ((~_l & 15) << 1)) - 1)
      y += occAux4(tmp) - (~_l & 15)
      m_tk.copyToArray(m_tl)
      m_tk(0) += x & 0xff
      m_tk(1) += (x >>> 8) & 0xff
      m_tk(2) += (x >>> 16) & 0xff
      m_tk(3) += (x >>> 24)

      m_tl(0) += y & 0xff
      m_tl(1) += (y >>> 8) & 0xff
      m_tl(2) += (y >>> 16) & 0xff
      m_tl(3) += (y >>> 24)
    }
  }

  def bwtExtend(is_back: Boolean) {
    var cond : Boolean = false
    if (is_back) {
      bwt_2occ4(m_ik.k - 1, m_ik.k - 1 + m_ik.s)
      for (i <- 0 to 3) {
	m_ok(i).k = mBWT.L2(i) + 1 + m_tk(i)
	m_ok(i).s = m_tl(i) - m_tk(i)
      }
      cond = ((m_ik.k <= mBWT.primary) && (m_ik.k + m_ik.s - 1 >= mBWT.primary))
      m_ok(3).l = m_ik.l
      if (cond) m_ok(3).l += 1
      m_ok(2).l = m_ok(3).l + m_ok(3).s
      m_ok(1).l = m_ok(2).l + m_ok(2).s
      m_ok(0).l = m_ok(1).l + m_ok(1).s
    }
    else {
      bwt_2occ4(m_ik.l - 1, m_ik.l - 1 + m_ik.s)
      for (i <- 0 to 3) {
	m_ok(i).l = mBWT.L2(i) + 1 + m_tk(i)
	m_ok(i).s = m_tl(i) - m_tk(i)
      }
      cond = ((m_ik.l <= mBWT.primary) && (m_ik.l + m_ik.s - 1 >= mBWT.primary))
      m_ok(3).k = m_ik.k
      if (cond) m_ok(3).k += 1
      m_ok(2).k = m_ok(3).k + m_ok(3).s
      m_ok(1).k = m_ok(2).k + m_ok(2).s
      m_ok(0).k = m_ok(1).k + m_ok(1).s
    }
  }

  //len: q's length
  def bwtSMem1(bwt: BWTType, len: Int, q: Array[Byte], x: Int, min_intv: Int, mem: MutableList[BWTIntvType], tmpvec_0: MutableList[BWTIntvType], tmpvec_1: MutableList[BWTIntvType]) : Int = {
    if (q(x) > 3)
      return x + 1

    var min_intv_copy : Int = 0
    if (min_intv < 1)
      min_intv_copy = 1
    else
      min_intv_copy = min_intv

    var prev = tmpvec_0
    var curr = tmpvec_1
    var swap : MutableList[BWTIntvType] = new MutableList[BWTIntvType]()
    mBWT = bwt
    // start pos for ik = 0?
    m_ik = new BWTIntvType(0, x + 1, bwt.L2(q(x)) + 1, bwt.L2(3 - q(x)) + 1, bwt.L2(q(x) + 1) - bwt.L2(q(x)))
    //var ok : Array[BWTIntvType] = new Array[BWTIntvType](4)
    var c : Int = 0

    breakable { for (i <- x + 1 to len - 1) { //forward search
      if (q(i) < 4) {	//an A/C/G/T base
	c = 3 - q(i)
	bwtExtend(false) //bwt_forward_extend
	if (m_ok(c).s != m_ik.s) {
	  curr.+=:(m_ik)
	  if (m_ok(c).s < min_intv) break
	}
	m_ik = m_ok(c)
	m_ik.endPoint = i + 1
      }
      else {
	curr.+=:(m_ik) //prepend the item to the list -- no need to reverse later
	break
      }
    } }
    
    var ret : Int = curr(0).endPoint
    swap = curr
    curr = prev
    prev = swap

    breakable { for (i <- x - 1 to -1) { //backward extension
      if (i < 0)
	c = -1
      else {
	if (q(i) < 4) c = q(i) else c = -1
      }

      breakable { for (j <- 0 to prev.length - 1) {
	m_ik = prev(j)
	bwtExtend(true) //bwt_extend on p
	if (c < 0 || m_ok(c).s < min_intv) {
	  if (curr.isEmpty) {
	    if (mem.isEmpty || i + 1 < mem.last.startPoint) {
	      m_ik.startPoint = i + 1
	      mem.+=:(m_ik)
	    }
	  }
	}
	else if (curr.isEmpty || m_ok(c).s != curr.last.s) {
	  m_ok(c).startPoint = m_ik.startPoint
	  m_ok(c).endPoint = m_ik.endPoint
	  curr += m_ok(c)
	}
      } }
      if (curr.isEmpty) break
      swap = curr
      curr = prev
      prev = swap
    } }

    ret
  }
}