package cs.ucla.edu.bwaspark.datatype

//BWAIdxType: maintaining all the information of BWA Index generated from FastA Reference
class BWAIdxType {  

  //1st: BWTType(".bwt", ".sa" files)
  var bwt: BWTType = _

  //2nd: BNTSeqType(".ann", ".amb" files)
  var bns: BNTSeqType = _

  //3rd: PACType(".pac" file)
  var pac: Array[Byte] = _  //uint8_t[]

  //loading files into fields
  //prefix: prefix of filenames
  //which: !!!to add!!!
  def load(prefix: String, which: Int) {
    //There is a function called "bwa_idx_infer_prefix" in original BWA,
    //but it seems to be useless

    //loading bwt
    if (which & BWA_IDX_BWT) {
      bwt = new BWTType
      bwt.load(prefix)
    }

    //loading bns
    if (which & BWA_IDX_BNS) {
      bns = new BNTSeqType
      bns.load(prefix)
      
      //loading pac
      if (which & BWA_IDX_PAC) {
        def pacLoader(filename: String, length: Long): Array[Byte] = {
          //to add: reading binary file
        }
        pac = pacLoader(prefix+".pac", bns.l_pac)
      }
    }
  }
}
