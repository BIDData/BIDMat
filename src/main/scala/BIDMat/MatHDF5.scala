package BIDMat
import ncsa.hdf.hdf5lib.structs._
import ncsa.hdf.hdf5lib.H5._
import ncsa.hdf.hdf5lib.HDF5Constants._

object MatHDF5 {
  var refcount:Long = -1

  def setCompressionPlist(dplist_id:Int, dims:Array[Long]) = {
	if (Mat.compressType > 0) {
	  if (dims.length == 1) {
		if (dims(0) > 1024) {
		  val cdims = new Array[Long](1)
		  cdims(0) = math.max(1, math.min(dims(0), Mat.chunkSize))
		  H5Pset_chunk(dplist_id, 1, cdims)
		  if (Mat.compressType == 1) {
		  	H5Pset_deflate(dplist_id, Mat.compressionLevel)
		  } else {
		  	H5Pset_szip(dplist_id, H5_SZIP_EC_OPTION_MASK, Mat.szipBlock)
		  }
		}
	  } else {
		if (dims(0)*dims(1) > 1024) {
		  val cdims = new Array[Long](2)
		  cdims(0) = math.max(1, math.min(dims(0), 1+Mat.chunkSize/dims(1)))
		  cdims(1) = math.max(1, dims(1))
		  if (Mat.compressType == 1) {
		  	H5Pset_deflate(dplist_id, Mat.compressionLevel)
		  } else {
		  	H5Pset_szip(dplist_id, H5_SZIP_EC_OPTION_MASK, Mat.szipBlock)
		  }
		}
	  }
	}
  }
  
  def getStringAttr(id:Int, obj_name:String, attr_name:String):String = { 
	val attr_id = H5Aopen_by_name(id, obj_name, attr_name, H5P_DEFAULT, H5P_DEFAULT)
	val attr_type_id = H5Aget_type(attr_id)
	val attr_type_size = H5Tget_size(attr_type_id)
	val sbuf = new Array[Byte](attr_type_size + 1)
	H5Aread(attr_id, attr_type_id, sbuf)
	H5Tclose(attr_type_id)
	H5Aclose(attr_id)
	new String(sbuf).trim()
  }

  def putStringAttr(id:Int, attr_name:String, attr_val:String) = { 
	val space_id = H5Screate(H5S_SCALAR)
	val memtype_id = H5Tcopy(H5T_FORTRAN_S1) 
	H5Tset_size(memtype_id, attr_val.length()) 
	val attr_id = H5Acreate(id, attr_name, memtype_id, space_id, H5P_DEFAULT, H5P_DEFAULT)
	H5Awrite(attr_id, memtype_id, attr_val.getBytes())
	H5Tclose(memtype_id)
	H5Aclose(attr_id)
	H5Sclose(space_id)
  }

  def getLongAttr(id:Int, obj_name:String, attr_name:String):Long = { 
	val attr_id = H5Aopen_by_name(id, obj_name, attr_name, H5P_DEFAULT, H5P_DEFAULT)
	val attr_type_id = H5Aget_type(attr_id)
	val attr_type_size = H5Tget_size(attr_type_id)
	val sbuf = new Array[Long](attr_type_size/8)
	H5Aread(attr_id, attr_type_id, sbuf)
	H5Tclose(attr_type_id)
	H5Aclose(attr_id)
	sbuf(0)
  }

  def putIntAttr(id:Int, attr_name:String, attr_val:Int) = { 
	val space_id = H5Screate(H5S_SCALAR)
	val attr_id = H5Acreate(id, attr_name, H5T_NATIVE_INT, space_id, H5P_DEFAULT, H5P_DEFAULT)
	val lbuf = Array[Int](1)
	lbuf(0) = attr_val
	H5Awrite(attr_id, H5T_NATIVE_INT, lbuf)
	H5Aclose(attr_id)
	H5Sclose(space_id)
  }

  def putLongAttr(id:Int, attr_name:String, attr_val:Long) = { 
	val space_id = H5Screate(H5S_SCALAR)
	val attr_id = H5Acreate(id, attr_name, H5T_NATIVE_LLONG, space_id, H5P_DEFAULT, H5P_DEFAULT)
	val lbuf = Array[Long](1)
	lbuf(0) = attr_val
	H5Awrite(attr_id, H5T_NATIVE_LLONG, lbuf)
	H5Aclose(attr_id)
	H5Sclose(space_id)
  }

  def putByteAttr(id:Int, attr_name:String, attr_val:Byte) = { 
	val space_id = H5Screate(H5S_SCALAR)
	val attr_id = H5Acreate(id, attr_name, H5T_NATIVE_UCHAR, space_id, H5P_DEFAULT, H5P_DEFAULT)
	val lbuf = Array[Byte](1)
	lbuf(0) = attr_val
	H5Awrite(attr_id, H5T_NATIVE_UCHAR, lbuf)
	H5Aclose(attr_id)
	H5Sclose(space_id)
  }

  def getMatDims(data_id:Int):Array[Long] = { 
	val space_id = H5Dget_space(data_id)
	val dims = new Array[Long](2)
	val ok = H5Sget_simple_extent_dims(space_id, dims, null)
	H5Sclose(space_id)
	dims
  }

  def readMatDims(fname:String, varname:String):(Long, Long) = { 
	val fid = H5Fopen(fname,H5F_ACC_RDONLY,H5P_DEFAULT)
	val data_id = H5Dopen(fid, varname, H5P_DEFAULT)
    val dims = getMatDims(data_id)
	H5Dclose(data_id)
	H5Fclose(fid)
    if (dims(1) == 0) { 
      (dims(0), dims(1))
    } else { 
      (dims(1), dims(0))
    }
  }

  def getDenseMat[T : ClassManifest](fid:Int, varname:String, h5class:Int, dsize:Int):DenseMat[T] = {
	val data_id = H5Dopen(fid, varname, H5P_DEFAULT)
	val data_type_id = H5Dget_type(data_id)
	val data_class = H5Tget_class(data_type_id)
	val data_size = H5Tget_size(data_type_id)
	val dims = getMatDims(data_id)
	var mdata:DenseMat[T] = null
	if (data_class == h5class && data_size == dsize) {
	  mdata = new DenseMat[T](dims(1).intValue, dims(0).intValue)
	  H5Dread(data_id, data_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, mdata.data)
	} else {
	  throw new RuntimeException("Bad class or data size "+data_class+" "+data_size)
	}
	H5Tclose(data_type_id)
	H5Dclose(data_id)
	mdata
  }

  def getCellMat(fid:Int, varname:String):CSMat = {
	val data_id = H5Dopen(fid, varname, H5P_DEFAULT)
	val data_type_id = H5Dget_type(data_id)
	val data_class = H5Tget_class(data_type_id)
	val data_size = H5Tget_size(data_type_id)
	val dims = getMatDims(data_id)
	var mdata:CSMat = null
	mdata = CSMat(dims(1).intValue, dims(0).intValue)
	val bdata = new Array[Array[Byte]]((dims(0)*dims(1)).intValue)
	for (i <- 0 until bdata.length) {
	  bdata(i) = new Array[Byte](data_size)
	}
	H5Dread(data_id, data_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, bdata)
	val obj_type_out = new Array[Int](1)
	obj_type_out(0) = 1
	for (i <- 0 until bdata.length) {
	  val item_id = H5Rdereference(fid, H5R_OBJECT, bdata(i))
	  mdata.data(i) = getMat(item_id, ".").asInstanceOf[String]
	  H5Oclose(item_id)
	}
	H5Tclose(data_type_id)
	H5Dclose(data_id)
	mdata
  }

  def getMatString(fid:Int, varname:String):String = {
	val data_id = H5Dopen(fid, varname, H5P_DEFAULT)
	val data_type_id = H5Dget_type(data_id)
	val data_class = H5Tget_class(data_type_id)
	val data_size = H5Tget_size(data_type_id)
	val dims = getMatDims(data_id)
	val nrows = dims(0).intValue
	val ncols = dims(1).intValue
	val sbuf = new Array[Byte](data_size*nrows*ncols)
	H5Dread(data_id, data_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, sbuf)
	H5Tclose(data_type_id)
	H5Dclose(data_id)
	new String(sbuf, "UTF_16LE").trim()
  }

  def getSparseMat[T](fid:Int, varname:String)(implicit manifest:Manifest[T], numeric:Numeric[T]):SparseMat[T] = {
	val nrows = getLongAttr(fid, varname, "MATLAB_sparse").intValue
	val jc_id = H5Dopen(fid, varname+"/jc", H5P_DEFAULT)
	val ncols = getMatDims(jc_id)(0).intValue - 1
	val data_id = H5Dopen(fid, varname+"/data", H5P_DEFAULT)
	val data_type_id = H5Dget_type(data_id)
	val nnz = getMatDims(data_id)(0).intValue
	var ir_id = -1
	try {
	  ir_id = H5Dopen(fid, varname+"/ir", H5P_DEFAULT)
	} catch {
	  case _ => {}
	}
	val sdata = if (ir_id >= 0) {
	  SparseMat(nrows, ncols, nnz) 
	} else {
	  SparseMat.noRows(nrows, ncols, nnz)
	}
	val convert_ints = H5Tcopy(H5T_NATIVE_INT)
	H5Dread_int(jc_id, convert_ints, H5S_ALL, H5S_ALL, H5P_DEFAULT, sdata.jc)
	addOne(sdata.jc)
	H5Dclose(jc_id)
	if (ir_id >= 0) {
	  H5Dread_int(ir_id, convert_ints, H5S_ALL, H5S_ALL, H5P_DEFAULT, sdata.ir)
	  addOne(sdata.ir)
	  H5Dclose(ir_id)
	}
	H5Tclose(convert_ints)
	H5Dread(data_id, data_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, sdata.data)
	H5Tclose(data_type_id)
	H5Dclose(data_id)
	sdata
  }

  def getMat(fid:Int, varname:String):AnyRef = {
	if (fid > 0 && H5Aexists_by_name(fid, varname, "MATLAB_class", H5P_DEFAULT)) {
	  val attr_class = getStringAttr(fid, varname, "MATLAB_class")
	  if (attr_class.equals("double")) {
		if (H5Aexists_by_name(fid, varname, "MATLAB_sparse", H5P_DEFAULT)) {
		  SDMat(getSparseMat[Double](fid, varname))
		} else {
		  DMat(getDenseMat[Double](fid, varname, H5T_FLOAT, 8))
		}
	  } else if (attr_class.equals("single")) {
		if (H5Aexists_by_name(fid, varname, "MATLAB_sparse", H5P_DEFAULT)) {
		  SMat(getSparseMat[Float](fid, varname))
		} else {
		  FMat(getDenseMat[Float](fid, varname, H5T_FLOAT, 4))
		}
	  } else if (attr_class.equals("int32")) {
		if (H5Aexists_by_name(fid, varname, "MATLAB_sparse", H5P_DEFAULT)) {
		  throw new RuntimeException("Sparse arrays of ints unsupported")
		} else {
		  IMat(getDenseMat[Int](fid, varname, H5T_INTEGER, 4))
		}
	  } else if (attr_class.equals("int8")) {
		if (H5Aexists_by_name(fid, varname, "MATLAB_sparse", H5P_DEFAULT)) {
		  BMat(getSparseMat[Byte](fid, varname))
		} else {
		  throw new RuntimeException("Dense arrays of bytes unsupported")
		}
	  } else if (attr_class.equals("char")) {
		if (H5Aexists_by_name(fid, varname, "MATLAB_sparse", H5P_DEFAULT)) {
		  throw new RuntimeException("Sparse arrays of char unsupported")
		} else {
		  getMatString(fid, varname)
		}
	  } else if (attr_class.equals("cell")) {
		if (H5Aexists_by_name(fid, varname, "MATLAB_sparse", H5P_DEFAULT)) {
		  throw new RuntimeException("Sparse cell arrays unsupported")
		} else {
		  getCellMat(fid, varname)
		}
	  } else throw new RuntimeException("Couldnt read storage class "+attr_class)
	} else throw new RuntimeException("Couldnt find matlab var named "+varname)
  }

  def writeMatHeader(fname:String) = {
	val ff = new java.io.RandomAccessFile(fname,"rws")
	val sp = new scala.sys.SystemProperties()
	val hstring = "MATLAB 7.3 MAT-file, Platform: "+sp.get("os.arch").get+" "+sp.get("os.name").get+" "+sp.get("os.version").get+ " "+
	"Created by BIDMat on "+(new java.text.SimpleDateFormat("EEE MMM d HH:mm:ss yyyy")).format(new java.util.Date())+
	" HDF5 Schema 1.0 ."
	val hb = hstring.getBytes()
	val hbytes = new Array[Byte](512)
	for (i <- 0 until 116) hbytes(i) = 32
	System.arraycopy(hb, 0, hbytes, 0, math.min(hstring.length(), 116))
	val version:Byte = 2
	hbytes(125) = version
	hbytes(126) = 0x49
	hbytes(127) = 0x4D
	ff.write(hbytes)
	//		ff.write(emptyHDF5file)
	ff.close()
  }

  def putDenseMat[T](fid:Int, a:DenseMat[T], aname:String, h5class:Int, matclass:String):Array[Byte] = {
	val dims = new Array[Long](2)
	dims(0) = a.ncols
	dims(1) = a.nrows
	val filespace_id = H5Screate_simple(2, dims, null)
	val dplist_id = H5Pcreate(H5P_DATASET_CREATE)
//	setCompressionPlist(dplist_id, dims)
	val dataset_id = H5Dcreate(fid, "/"+aname, h5class, filespace_id, H5P_DEFAULT, dplist_id, H5P_DEFAULT)
	H5Dwrite(dataset_id, h5class, H5S_ALL, H5S_ALL, H5P_DEFAULT, a.data)
	H5Pclose(dplist_id)
	putStringAttr(dataset_id, "MATLAB_class", matclass)
	val ref = H5Rcreate(dataset_id, ".", H5R_OBJECT, -1)
	H5Dclose(dataset_id)
	H5Sclose(filespace_id)
	ref
  }

  def putEmptyRef(id:Int):Array[Byte] = {
	val dims = new Array[Long](1)
	dims(0) = 2
	val tmp = Array[Long](2)
	val dmatspace_id = H5Screate_simple(1, dims, null)
	val dmat_id = H5Dcreate(id, "0", H5T_NATIVE_ULLONG, dmatspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)
	H5Dwrite(dmat_id, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp)
	putStringAttr(dmat_id, "MATLAB_class", "canonical empty")
	putByteAttr(dmat_id, "MATLAB_empty", 1)		
	val ref = H5Rcreate(dmat_id, ".", H5R_OBJECT, -1)
	H5Dclose(dmat_id)
	H5Sclose(dmatspace_id)
	ref
  }

  def putCellMat(fid:Int, varname:String, a:CSMat) = {
	var group_id = 0
	if (refcount < 0) { 
	  group_id = H5Gcreate(fid, "/#refs#", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)
	  putEmptyRef(group_id)
	  refcount = 1
	} else {
	  group_id = H5Gopen(fid, "/#refs#", H5P_DEFAULT)
	}
	val refIds = new Array[Array[Byte]](a.length)
	for (i <- 0 until a.length) {
	  val newname = "%x" format refcount
	  refcount += 1
	  refIds(i) = putMat(group_id, a.data(i), newname)
	}
	val dims = new Array[Long](2)
	dims(0) = a.ncols
	dims(1) = a.nrows
	val dplist_id = H5Pcreate(H5P_DATASET_CREATE)
	setCompressionPlist(dplist_id, dims)
	val refspace_id = H5Screate_simple(2, dims, null)
	val refs_id = H5Dcreate(fid, varname, H5T_STD_REF_OBJ, refspace_id, H5P_DEFAULT, dplist_id, H5P_DEFAULT)
	H5Dwrite(refs_id, H5T_STD_REF_OBJ, H5S_ALL, H5S_ALL, H5P_DEFAULT, refIds)
	putStringAttr(refs_id, "MATLAB_class", "cell")
	val ref = H5Rcreate(refs_id, ".", H5R_OBJECT, -1)
	H5Dclose(refs_id)
	H5Sclose(refspace_id)
	H5Pclose(dplist_id)
	H5Gclose(group_id)
	ref
  }

  def putSparseMat[T](fid:Int, a:SparseMat[T], varname:String, nativeClass:Int, className:String):Array[Byte] = {
	val dims = new Array[Long](1)
	val group_id = H5Gcreate(fid, "/"+varname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)
	putStringAttr(group_id, "MATLAB_class", className)
	putLongAttr(group_id, "MATLAB_sparse", a.nrows)
	val convert_ints = H5Tcopy(H5T_NATIVE_INT)
	dims(0) = a.ncols + 1
	var dplist_id = H5Pcreate(H5P_DATASET_CREATE)
	setCompressionPlist(dplist_id, dims)
	val jcs_id = H5Screate_simple(1, dims, null)
	val jc_id = H5Dcreate(group_id, "jc", H5T_NATIVE_LLONG, jcs_id, H5P_DEFAULT, dplist_id, H5P_DEFAULT)
	subOne(a.jc)
	try {
	  H5Dwrite(jc_id, convert_ints, H5S_ALL, H5S_ALL, H5P_DEFAULT, a.jc)
	} catch  {
	  case e => {
		addOne(a.jc)
		throw new RuntimeException("Error writing sparse mat "+e)
	  }
	}
	addOne(a.jc)
	H5Dclose(jc_id)
	H5Sclose(jcs_id)
	H5Pclose(dplist_id)

	dims(0) = a.nnz
	dplist_id = H5Pcreate(H5P_DATASET_CREATE)
	setCompressionPlist(dplist_id, dims)
	if (a.ir != null) {
	  val irs_id = H5Screate_simple(1, dims, null)
	  val ir_id = H5Dcreate(group_id, "ir", H5T_NATIVE_LLONG, irs_id, H5P_DEFAULT, dplist_id, H5P_DEFAULT)
	  subOne(a.ir)
	  try {
		H5Dwrite(ir_id, convert_ints, H5S_ALL, H5S_ALL, H5P_DEFAULT, a.ir)
	  } catch  {
		case e => {
		  addOne(a.ir)
		  throw new RuntimeException("Error writing sparse mat "+e)
		}
	  }
	  addOne(a.ir)
	  H5Dclose(ir_id)
	  H5Sclose(irs_id)
	}

	val dataspace_id = H5Screate_simple(1, dims, null)
	val data_id = H5Dcreate(group_id, "data", nativeClass, dataspace_id, H5P_DEFAULT, dplist_id, H5P_DEFAULT)
	H5Dwrite(data_id, nativeClass, H5S_ALL, H5S_ALL, H5P_DEFAULT, a.data)
	H5Dclose(data_id)
	H5Sclose(dataspace_id)
	H5Pclose(dplist_id)
	H5Tclose(convert_ints)
	val ref = H5Rcreate(group_id, ".", H5R_OBJECT, -1)
	H5Gclose(group_id)
	ref
  }

  def putMatString(id:Int, varname:String, str:String):Array[Byte] = { 
	val dims = new Array[Long](2)
	dims(0) = str.length
	dims(1) = 1
	val dplist_id = H5Pcreate(H5P_DATASET_CREATE)
	setCompressionPlist(dplist_id, dims)
	val sbytes = str.getBytes("UTF_16LE")
	val strspace_id = H5Screate_simple(2, dims, null) 
	val str_id = H5Dcreate(id, varname, H5T_NATIVE_USHORT, strspace_id, H5P_DEFAULT, dplist_id, H5P_DEFAULT)
	putStringAttr(str_id, "MATLAB_class", "char")
	putIntAttr(str_id, "MATLAB_int_decode", 2)
	H5Dwrite(str_id, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sbytes)
	val ref = H5Rcreate(str_id, ".", H5R_OBJECT, -1)
	H5Dclose(str_id)
	H5Sclose(strspace_id)
	H5Pclose(dplist_id)
	ref
  }

  def putMat(fid:Int, a:AnyRef, aname:String):Array[Byte] = { 
	a match { 
	  case aa:DMat => putDenseMat[Double](fid, aa, aname, H5T_NATIVE_DOUBLE, "double")
	  case aa:FMat => putDenseMat[Float](fid, aa, aname, H5T_NATIVE_FLOAT, "single")
	  case aa:IMat => putDenseMat[Int](fid, aa, aname, H5T_NATIVE_INT, "int32")
	  case aa:BMat => putSparseMat[Byte](fid, aa, aname, H5T_NATIVE_CHAR, "int8")
	  case aa:SMat => putSparseMat[Float](fid, aa, aname, H5T_NATIVE_FLOAT, "single")
	  case aa:SDMat => putSparseMat[Double](fid, aa, aname, H5T_NATIVE_DOUBLE, "double")
	  case aa:CSMat => putCellMat(fid, aname, aa)
	  case aa:String => putMatString(fid, aname, aa)
	  case _ => throw new RuntimeException("unsupported matrix type to save")
	}
  }

  def hload(fname:String, vname:String):AnyRef = {
  val fapl = H5Pcreate(H5P_FILE_ACCESS)
//  H5Pset_fapl_core(fapl, 16*1024*1024, false);  println("core driver")
  if (Mat.useStdio) H5Pset_fapl_stdio(fapl); //println("stdio driver")                     // Not working on windows
	val fid = H5Fopen(fname,H5F_ACC_RDONLY,fapl)
	H5Pclose(fapl)
	val mat = getMat(fid, vname)
	H5Fclose(fid)
	mat
  }

  def hload(fname:String, vnames:List[String]):List[AnyRef] = {
  val fapl = H5Pcreate(H5P_FILE_ACCESS)
//  H5Pset_fapl_core(fapl, 32*1024*1024, false);  println("core driver")
  if (Mat.useStdio)  H5Pset_fapl_stdio(fapl); //println("stdio driver")
	val fid = H5Fopen(fname,H5F_ACC_RDONLY,fapl)
	H5Pclose(fapl)
	val mats = vnames.map((vname) => getMat(fid, vname))
	H5Fclose(fid)
	mats
  }

  def hsaveAsHDF5(fname:String, args:List[AnyRef]) = {
	refcount = -1
	val fapl_id = H5Pcreate (H5P_FILE_ACCESS)
	if (Mat.useStdio) H5Pset_fapl_stdio(fapl_id)
	val fid = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id)
	H5Pclose(fapl_id)
	saveAsImpl(fid, args)
	H5Fclose(fid)
  }

  def hsaveAs(fname:String, args:List[AnyRef]) = {
	refcount = -1
	val fapl_id = H5Pcreate (H5P_FILE_ACCESS)
	if (Mat.useStdio) H5Pset_fapl_stdio(fapl_id)
	val fcplist_id = H5Pcreate(H5P_FILE_CREATE)
	H5Pset_userblock(fcplist_id, 512)
	val fid = H5Fcreate(fname, H5F_ACC_TRUNC, fcplist_id, fapl_id)
	H5Pclose(fcplist_id)
	H5Pclose(fapl_id)
	saveAsImpl(fid, args)
	H5Fclose(fid)
	writeMatHeader(fname)
  }

  def saveAsImpl(fid:Int, argList:List[AnyRef]) = {
	var i = 0
	while (i < argList.length) {
	  argList(i) match {
		case a:Mat => {
		  argList(i+1) match {
			case str:String => putMat(fid, a, str)
			case _ => throw new RuntimeException("odd numbered args must be String variable names")
		  }
		}
		case _ => throw new RuntimeException("even numbered args must be Mat variables")
	  }
	  i += 2
	}
  }

  def addOne(ii:Array[Int]) = {
	if (Mat.ioneBased == 1) {
	  var i = 0
	  while (i < ii.length) {
		ii(i) += 1
		i += 1
	  }
	}
  }

  def subOne(ii:Array[Int]) = {
	if (Mat.ioneBased == 1) {
	  var i = 0
	  while (i < ii.length) {
		ii(i) = ii(i) - 1
		i += 1
	  }
	}
  }
}
