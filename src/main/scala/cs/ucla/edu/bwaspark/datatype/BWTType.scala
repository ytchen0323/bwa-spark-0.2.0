package cs.ucla.edu.bwaspark.datatype

import java.io.{FileInputStream, IOException}
import java.nio.{ByteBuffer}
import java.nio.file.{Files, Path, Paths}
import java.nio.channels.FileChannel
import java.nio.ByteOrder

import scala.util.control.Breaks._
import scala.List


class BWTType {
  // Data Structure
  var primary: Long = _
  var L2: Array[Long] = new Array[Long](5)
  var seqLen: Long = _
  var bwtSize: Long = _
  var bwt: Array[Int] = _
  var cntTable: Array[Int] = new Array[Int](250)
  var saIntv: Int = _
  var numSa: Long = _
  var sa: Array[Long] = _

  val readBufSize = 0x80000

  def BWTLoad(bwtFile: String) {
    val fs = new FileInputStream(bwtFile)
    val fc = fs.getChannel
    val fileSize = fc.size
    val buf = ByteBuffer.allocate(readBufSize)
    buf.order(ByteOrder.LITTLE_ENDIAN)

    // Read primary member variable in BWTType
    val smallBuf = ByteBuffer.allocate(8)
    smallBuf.order(ByteOrder.nativeOrder)
    fc.read(smallBuf)
    smallBuf.flip
    primary = smallBuf.getLong(0)
    smallBuf.rewind

    // Read L[1] - L[4] member variables in BWTType
    L2(0) = 0
    for(i <- 1 to 4) {
      fc.read(smallBuf)
      smallBuf.flip
      L2(i) = smallBuf.getLong(0)
      smallBuf.rewind
    }
    seqLen = L2(4)

    // Read the bwt array in BWTType
    var ret = 0
    var bytesRead: Long = 0
    var numIntRead = 0
    bwtSize = (fileSize - 8 * 5) >> 2

    println("File size: " + fileSize)

    // Potential array boundary problem here
    // Enforce the type of numElem to be Int
    // This is because the size of a scala array cannot exceed the Int upper bound
    bwt = new Array[Int](bwtSize.asInstanceOf[Int])

    while(ret >= 0) {
      ret = fc.read(buf)
      buf.flip
      bytesRead += ret
      var curBytesRead = 0

      // Fill the data from buf
      while(curBytesRead < readBufSize && buf.hasRemaining) {
        val bwtSeg = buf.getInt
        bwt(numIntRead) = bwtSeg
        curBytesRead += 4
        numIntRead += 1
      }

      buf.rewind
    } 

    fc.close

    // Debugging messages
    println("Bytes Read: " + (bytesRead + 1))
    println("Num Integers Read: " + numIntRead)
    println("primary: " + primary)
    println("bwtSize: " + bwtSize)
    println("seqLen: " + seqLen)
    for(i <- 0 to 4) 
      println("L2(" + i + "): " + L2(i))
  }
}
