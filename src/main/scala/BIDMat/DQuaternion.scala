//-*-coding:utf-8-*-
package BIDMat

class DQuaternion(data0:Array[Double]) extends DMat(Array(1,4), data0) {

    def c:DQuaternion = cc(null);

    def c(omat:DQuaternion):DQuaternion = cc(omat);

    def cc(omat:DQuaternion):DQuaternion = {
  	val out = DQuaternion.newOrCheckDQuaternion(omat, GUID, "DQuaternionc".##)
	val odata = out.data;
	odata(0) = data(0);
	odata(1) = -data(1);
	odata(2) = -data(2);
	odata(3) = -data(3);
	out
    }

    override def inv:DQuaternion = cc(null);

    def inv(omat:DQuaternion):DQuaternion = inv_(omat);

    def inv_(omat:DQuaternion):DQuaternion = {
  	val out = DQuaternion.newOrCheckDQuaternion(omat, GUID, "DQuaternioninv".##)
	val odata = out.data;
	val v = normsq;
	odata(0) = data(0)/v;
	odata(1) = -data(1)/v;
	odata(2) = -data(2)/v;
	odata(3) = -data(3)/v;
	out
    }

    def * (b:DQuaternion):DQuaternion = {
	val out = DQuaternion.newOrCheckDQuaternion(null, GUID, b.GUID, "DQuaternion*".##)
	val odata = out.data;
	val bdata = b.data;
	odata(0) = data(0) * bdata(0) - data(1) * bdata(1) - data(2) * bdata(2) - data(3) * bdata(3);
	odata(1) = data(0) * bdata(1) + data(1) * bdata(0) + data(2) * bdata(3) - data(3) * bdata(2);
	odata(2) = data(0) * bdata(2) + data(2) * bdata(0) - data(1) * bdata(3) + data(3) * bdata(1);
	odata(3) = data(0) * bdata(3) + data(3) * bdata(0) + data(1) * bdata(2) - data(2) * bdata(1);
	out;
    };

    override def * (m:DMat):DQuaternion = {
	this * DQuaternion(m);
    }

    def *^ (b:DQuaternion):DQuaternion = {
	val out = DQuaternion.newOrCheckDQuaternion(null, GUID, b.GUID, "DQuaternion*^".##)
	val odata = out.data;
	val bdata = b.data;
	odata(0) = data(0) * bdata(0) + data(1) * bdata(1) + data(2) * bdata(2) + data(3) * bdata(3);
	odata(1) = - data(0) * bdata(1) + data(1) * bdata(0) - data(2) * bdata(3) + data(3) * bdata(2);
	odata(2) = - data(0) * bdata(2) + data(2) * bdata(0) + data(1) * bdata(3) - data(3) * bdata(1);
	odata(3) = - data(0) * bdata(3) + data(3) * bdata(0) - data(1) * bdata(2) + data(2) * bdata(1);
	out;
    };

    override def *^ (m:DMat):DQuaternion = {
	this *^ DQuaternion(m);
    }

    def ^* (b:DQuaternion):DQuaternion = {
	val out = DQuaternion.newOrCheckDQuaternion(null, GUID, b.GUID, "DQuaternion^*".##)
	val odata = out.data;
	val bdata = b.data;
	odata(0) = data(0) * bdata(0) + data(1) * bdata(1) + data(2) * bdata(2) + data(3) * bdata(3);
	odata(1) = data(0) * bdata(1) - data(1) * bdata(0) - data(2) * bdata(3) + data(3) * bdata(2);
	odata(2) = data(0) * bdata(2) - data(2) * bdata(0) + data(1) * bdata(3) - data(3) * bdata(1);
	odata(3) = data(0) * bdata(3) - data(3) * bdata(0) - data(1) * bdata(2) + data(2) * bdata(1);
	out;
    };

    override def ^* (m:DMat):DQuaternion = {
	this ^* DQuaternion(m);
    };

    def + (b:DQuaternion):DQuaternion = {
	val out = DQuaternion.newOrCheckDQuaternion(null, GUID, b.GUID, "DQuaternion+".##)
	val odata = out.data;
	val bdata = b.data;
	odata(0) = data(0) + bdata(0);
	odata(1) = data(1) + bdata(1);
	odata(2) = data(2) + bdata(2);
	odata(3) = data(3) + bdata(3);
	out;
    };

    def - (b:DQuaternion):DQuaternion = {
	val out = DQuaternion.newOrCheckDQuaternion(null, GUID, b.GUID, "DQuaternion-".##)
	val odata = out.data;
	val bdata = b.data;
	odata(0) = data(0) - bdata(0);
	odata(1) = data(1) - bdata(1);
	odata(2) = data(2) - bdata(2);
	odata(3) = data(3) - bdata(3);
	out;
    };

    def / (q:DQuaternion):DQuaternion = {
	this *^ q / q.normsq;
    }

    override def * (v:Double):DQuaternion = {
	val out = DQuaternion.newOrCheckDQuaternion(null, GUID, v.hashCode, "DQuaternion*Scalar".##);
	val odata = out.data;
	odata(0) = data(0)*v;
	odata(1) = data(1)*v;
	odata(2) = data(2)*v;
	odata(3) = data(3)*v;
	out;
    }

    override def / (v:Double):DQuaternion = {
	val out = DQuaternion.newOrCheckDQuaternion(null, GUID, v.hashCode, "DQuaternion/Scalar".##);
	val odata = out.data;
	odata(0) = data(0)/v;
	odata(1) = data(1)/v;
	odata(2) = data(2)/v;
	odata(3) = data(3)/v;
	out;
    }

    override def + (v:Double):DQuaternion = {
	val out = DQuaternion.newOrCheckDQuaternion(null, GUID, v.hashCode, "DQuaternion+Scalar".##);
	val odata = out.data;
	odata(0) = data(0)+v;
	odata(1) = data(1)+v;
	odata(2) = data(2)+v;
	odata(3) = data(3)+v;
	out;
    }

    override def - (v:Double):DQuaternion = {
	val out = DQuaternion.newOrCheckDQuaternion(null, GUID, v.hashCode, "DQuaternion-Scalar".##);
	val odata = out.data;
	odata(0) = data(0)-v;
	odata(1) = data(1)-v;
	odata(2) = data(2)-v;
	odata(3) = data(3)-v;
	out;
    }


    def toMat_(omat:DMat):DMat = {
	val out = DMat.newOrCheckDMat(4, 4, omat, GUID, "DQuaternionToMat".##);
	val odata = out.data;
	odata(0) = data(0);
	odata(1) = data(1);
	odata(2) = data(2);
	odata(3) = data(3);

	odata(4) = -data(1);
	odata(5) = data(0);
	odata(6) = data(3);
	odata(7) = -data(2);

	odata(8) = -data(2);
	odata(9) = -data(3);
	odata(10) = data(0);
	odata(11) = data(1);

	odata(12) = -data(3);
	odata(13) = data(2);
	odata(14) = -data(1);
	odata(15) = data(0);

	out;
    };

    def toMat:DMat = {
	toMat_(null);
    };

    def toMat(omat:DMat):DMat = {
	toMat_(omat);
    };

    def toRotMat_(omat:DMat):DMat = {
	val out = DMat.newOrCheckDMat(3, 3, omat, GUID, "DQuaternionToRotMat".##);
	val odata = out.data;
	val q0 = data(0);
	val qx = data(1);
	val qy = data(2);
	val qz = data(3);
	
	odata(0) = q0*q0 + qx*qx - qy*qy - qz*qz;
	odata(1) = 2*(qx*qy + q0*qz);
	odata(2) = 2*(qx*qz - q0*qy);

	odata(3) = 2*(qx*qy - q0*qz);
	odata(4) = q0*q0 - qx*qx + qy*qy - qz*qz;
	odata(5) = 2*(qy*qz + q0*qx);

	odata(6) = 2*(qx*qz + q0*qy);
	odata(7) = 2*(qy*qz - q0*qx);
	odata(8) = q0*q0 - qx*qx - qy*qy + qz*qz;

	out;
    };

    def toRotMat:DMat = {
	toRotMat_(null);
    };

    def toRotMat(omat:DMat):DMat = {
	toRotMat_(omat);
    };

    def normsq:Double = {
	data(0)*data(0) + data(1)*data(1) + data(2)*data(2) + data(3)*data(3);
    }

    def norm:Double = {
	math.sqrt(data(0)*data(0) + data(1)*data(1) + data(2)*data(2) + data(3)*data(3));
    }
}

object DQuaternion {
    def apply():DQuaternion = {
	new DQuaternion(new Array[Double](4));
    }

    def apply(m:DMat):DQuaternion = {
	val out = newOrCheckDQuaternion(null, m.GUID, "DQuaternionFromDMat".##);
	if (m.length == 3) {
	    out.data(0) = 0;
	    out.data(1) = m.data(0);
	    out.data(2) = m.data(1);
	    out.data(3) = m.data(2);
	} else if (m.length == 4) {
	    out.data(0) = m.data(0);
	    out.data(1) = m.data(1);
	    out.data(2) = m.data(2);
	    out.data(3) = m.data(3);
	} else if (m.nrows == 4 && m.ncols == 4) {
	    out.data(0) = m.data(0);
	    out.data(1) = m.data(1);
	    out.data(2) = m.data(2);
	    out.data(3) = m.data(3);
	} else {
 	    throw new RuntimeException("DQuaternion constructor: DMat must be a 3 or 4 dimensional vector, or a 4x4 Matrix");
	}
	out
    }

    def fromRotMat(m:DMat) = {
	val out = newOrCheckDQuaternion(null, m.GUID, "DQuaternionFromRotMat".##);
	val m00 = m.data(0);
	val m10 = m.data(1);
	val m20 = m.data(2);
	val m01 = m.data(3);
	val m11 = m.data(4);
	val m21 = m.data(5);
	val m02 = m.data(6);
	val m12 = m.data(7);
	val m22 = m.data(8);
	
	val tr = m00 + m11 + m22;
	var q0 = 0.0;
	var qx = 0.0;
	var qy = 0.0;
	var qz = 0.0;
	if (tr > 0) { // q0^2 > 0.25;
	    q0 = math.sqrt(tr+1.0) * 0.5; 
	    val scale = 0.25f / q0;
	    qx = (m21 - m12) * scale;
	    qy = (m02 - m20) * scale; 
	    qz = (m10 - m01) * scale; 
	} else if (m00 - m11 - m22 > 0) { // qx^2 > 0.25;
	    qx = math.sqrt(1.0 + m00 - m11 - m22) * 0.5;
	    val scale = 0.25f / qx;
	    q0 = (m21 - m12) * scale;
	    qy = (m01 + m10) * scale; 
	    qz = (m02 + m20) * scale; 
	} else if (m11 - m00 - m22 > 0) { // qy^2 > 0.25;
	    qy = math.sqrt(1.0 + m11 - m00 - m22) * 0.5;
	    val scale = 0.25f / qy;
	    q0 = (m02 - m20) * scale;
	    qx = (m01 + m10) * scale; 
	    qz = (m12 + m21) * scale; 
	} else {  // qz^2 is at least 0.25
	    qz = math.sqrt(1.0 + m22 - m00 - m11) * 0.5;
	    val scale = 0.25f / qz;
	    q0 = (m10 - m01) * scale;
	    qx = (m02 + m20) * scale;
	    qy = (m12 + m21) * scale;
	}
	val odata = out.data;
	odata(0) = q0;
	odata(1) = qx;
	odata(2) = qy;
	odata(3) = qz;
	out;
    }

    def apply(q0:Double, qx:Double, qy:Double, qz:Double) = {
	val out = new DQuaternion(new Array[Double](4));
	out.data(0) = q0;
	out.data(1) = qx;
	out.data(2) = qy;
	out.data(3) = qz;
	out;
    }
	

    def identity():DQuaternion = {
	val out = new DQuaternion(new Array[Double](4));
	out.data(0) = 1;
	out.data(0) = 0;
	out.data(0) = 0;
	out.data(0) = 0;
	out;
    }

    def zero():DQuaternion = {
	val out = new DQuaternion(new Array[Double](4));
	out.data(0) = 0;
	out.data(0) = 0;
	out.data(0) = 0;
	out.data(0) = 0;
	out;
    }

    def newOrCheckDQuaternion(outmat:DQuaternion):DQuaternion = {
	if (outmat.asInstanceOf[AnyRef] == null) {
	    DQuaternion();
	} else {
	    outmat;
	}
    }
  
    def newOrCheckDQuaternion(outmat:DQuaternion, matGuid:Long, opHash:Int, forceCache:Boolean):DQuaternion = {
	if (outmat.asInstanceOf[AnyRef] != null || (!Mat.useCache && !forceCache)) {
	    newOrCheckDQuaternion(outmat);
	} else {
	    val key = (matGuid, opHash);
	    val res = Mat.cache2(key).asInstanceOf[DQuaternion];
	    if (res != null) {
		newOrCheckDQuaternion(res);
	    } else {
		val omat = newOrCheckDQuaternion(null);
		Mat.cache2put(key, omat);
		omat
	    }
	}
    }

    def newOrCheckDQuaternion(outmat:DQuaternion,guid1:Long, opHash:Int):DQuaternion =
	newOrCheckDQuaternion(outmat, guid1, opHash, false);

  
    def newOrCheckDQuaternion(outmat:DQuaternion, guid1:Long, guid2:Long, opHash:Int, forceCache:Boolean):DQuaternion = {
	if (outmat.asInstanceOf[AnyRef] != null || (!Mat.useCache && !forceCache)) {
	    newOrCheckDQuaternion(outmat);
	} else {
	    val key = (guid1, guid2, opHash);
	    val res = Mat.cache3(key).asInstanceOf[DQuaternion];
	    if (res != null) {
		newOrCheckDQuaternion(res);
	    } else {
		val omat = newOrCheckDQuaternion(null);
		Mat.cache3put(key, omat);
		omat;
	    }
	}
    }
  
    def newOrCheckDQuaternion(outmat:DQuaternion,guid1:Long, guid2:Long, opHash:Int):DQuaternion =
	newOrCheckDQuaternion(outmat, guid1, guid2, opHash, false);
  
    def newOrCheckDQuaternion(outmat:DQuaternion, guid1:Long, guid2:Long, guid3:Long, opHash:Int, forceCache:Boolean):DQuaternion = {
	if (outmat.asInstanceOf[AnyRef] != null || (!Mat.useCache && !forceCache)) {
	    newOrCheckDQuaternion(outmat);
	} else {
	    val key = (guid1, guid2, guid3, opHash);
	    val res = Mat.cache4(key).asInstanceOf[DQuaternion];
	    if (res != null) {
		newOrCheckDQuaternion(res);
	    } else {
		val omat = newOrCheckDQuaternion(null);
		Mat.cache4put(key, omat);
		omat
	    }
	}
    }
  
    def newOrCheckDQuaternion(outmat:DQuaternion,guid1:Long, guid2:Long, guid3:Long, opHash:Int):DQuaternion =
	newOrCheckDQuaternion(outmat, guid1, guid2, guid3, opHash, false);

}
