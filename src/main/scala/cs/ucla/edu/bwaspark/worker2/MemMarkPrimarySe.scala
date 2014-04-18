package cs.ucla.edu.bwaspark.worker2

import scala.util.control.Breaks._
import scala.List
import scala.math.abs
import scala.collection.mutable.MutableList

import cs.ucla.edu.bwaspark.datatype._

//MemMarkPrimarySe
object MemMarkPrimarySe {
  
  
  /**
   * The main function of memMarkPrimarySe class

   * @param opt the MemOptType object
   * @param n the Int object
   * @param a the MemAlnRegType object
   * @param id the Long object
   */
  //??opt a list or a object??
  def memMarkPrimarySe(opt: MemOptType, n: Int, a: List[MemAlnRegType], id: Long) {
    var i: Int = 0
    var j: Int = 0
    var tmp: Int = 0
    var k: Int = 0
    var z: MutableList[Int] = new MutableList

    if( n != 0) {
      for( i <- 0 until n ) {
	a(i).sub = 0
	a(i).secondary = -1
	a(i).hash = hash64((id + i.toLong))
      }

      //ks_introsort(mem_ars_hash, n, a)
      //#define alnreg_hlt(a, b) ((a).score > (b).score || ((a).score == (b).score && (a).hash < (b).hash))
      a.sortWith( (x, y) => ((x.score > y.score) || ( x.score == y.score && x.hash < y.hash)  ) )
      tmp = opt.a + opt.b
      if((opt.oDel + opt.eDel) > tmp) {
	tmp = opt.oDel + opt.eDel
      }
      if((opt.oIns + opt.eIns) > tmp) {
	tmp = opt.oIns + opt.eIns
      }
      //kv_push()
      z += 0
      for(i <- 1 until n) {
	breakable {
	  for(k <- 0 until z.size) {
	    j = z(k)
	    var bMax: Int = if(a(j).qBeg > a(i).qBeg) a(j).qBeg else a(i).qBeg
	    var eMin: Int = if(a(j).qEnd < a(i).qEnd) a(j).qEnd else a(i).qEnd
	    // have overlap
	    if( eMin > bMax ) {
	      var minL: Int = if ((a(i).qEnd - a(i).qBeg)<(a(j).qEnd - a(j).qBeg)) (a(i).qEnd - a(i).qBeg) else (a(j).qEnd - a(j).qBeg)
	      //have significant overlap
	      if((eMin - bMax)>= minL * opt.maskLevel) {
		if(a(j).sub == 0) a(j).sub = a(i).score
		if((a(j).score - a(i).score) <= tmp) a(j).subNum = a(j).subNum + 1
		break
	      }
	    }
	  }
	  if(k == z.size) z += i
	  else a(i).secondary = z(k)
	}
      }
    }
  }

  private def hash64( key: Long ) : Long = {
    var keyVar: Long = key
    keyVar += ~(keyVar << 32)
    keyVar ^= (keyVar >>> 22)
    keyVar += ~(keyVar << 13)
    keyVar ^= (keyVar >>> 8)
    keyVar += (keyVar << 3)
    keyVar ^= (keyVar >>> 15)
    keyVar += ~(keyVar <<27)
    keyVar ^= (keyVar >>> 31)
    keyVar
  }
}
