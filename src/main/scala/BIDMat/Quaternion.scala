//-*-coding:utf-8-*-
package BIDMat

class Quaternion(data0:Array[Float]) extends FMat(Array(1,4), data0) {

    def c:Quaternion = cc(null);

    def c(omat:Quaternion):Quaternion = cc(omat);

    def cc(omat:Quaternion):Quaternion = {
  	val out = Quaternion.newOrCheckQuaternion(omat, GUID, "Quaternionc".##)
	val odata = out.data;
	odata(0) = data(0);
	odata(1) = -data(1);
	odata(2) = -data(2);
	odata(3) = -data(3);
	out
    }

    override def inv:Quaternion = cc(null);

    def inv(omat:Quaternion):Quaternion = inv_(omat);

    def inv_(omat:Quaternion):Quaternion = {
  	val out = Quaternion.newOrCheckQuaternion(omat, GUID, "Quaternioninv".##)
	val odata = out.data;
	val v = normsq;
	odata(0) = data(0)/v;
	odata(1) = -data(1)/v;
	odata(2) = -data(2)/v;
	odata(3) = -data(3)/v;
	out
    }

    def * (b:Quaternion):Quaternion = {
	val out = Quaternion.newOrCheckQuaternion(null, GUID, b.GUID, "Quaternion*".##)
	val odata = out.data;
	val bdata = b.data;
	odata(0) = data(0) * bdata(0) - data(1) * bdata(1) - data(2) * bdata(2) - data(3) * bdata(3);
	odata(1) = data(0) * bdata(1) + data(1) * bdata(0) + data(2) * bdata(3) - data(3) * bdata(2);
	odata(2) = data(0) * bdata(2) + data(2) * bdata(0) - data(1) * bdata(3) + data(3) * bdata(1);
	odata(3) = data(0) * bdata(3) + data(3) * bdata(0) + data(1) * bdata(2) - data(2) * bdata(1);
	out;
    };

    override def * (m:FMat):Quaternion = {
	this * Quaternion(m);
    }

    def *^ (b:Quaternion):Quaternion = {
	val out = Quaternion.newOrCheckQuaternion(null, GUID, b.GUID, "Quaternion*^".##)
	val odata = out.data;
	val bdata = b.data;
	odata(0) = data(0) * bdata(0) + data(1) * bdata(1) + data(2) * bdata(2) + data(3) * bdata(3);
	odata(1) = - data(0) * bdata(1) + data(1) * bdata(0) - data(2) * bdata(3) + data(3) * bdata(2);
	odata(2) = - data(0) * bdata(2) + data(2) * bdata(0) + data(1) * bdata(3) - data(3) * bdata(1);
	odata(3) = - data(0) * bdata(3) + data(3) * bdata(0) - data(1) * bdata(2) + data(2) * bdata(1);
	out;
    };

    override def *^ (m:FMat):Quaternion = {
	this *^ Quaternion(m);
    }

    def ^* (b:Quaternion):Quaternion = {
	val out = Quaternion.newOrCheckQuaternion(null, GUID, b.GUID, "Quaternion^*".##)
	val odata = out.data;
	val bdata = b.data;
	odata(0) = data(0) * bdata(0) + data(1) * bdata(1) + data(2) * bdata(2) + data(3) * bdata(3);
	odata(1) = data(0) * bdata(1) - data(1) * bdata(0) - data(2) * bdata(3) + data(3) * bdata(2);
	odata(2) = data(0) * bdata(2) - data(2) * bdata(0) + data(1) * bdata(3) - data(3) * bdata(1);
	odata(3) = data(0) * bdata(3) - data(3) * bdata(0) - data(1) * bdata(2) + data(2) * bdata(1);
	out;
    };

    override def ^* (m:FMat):Quaternion = {
	this ^* Quaternion(m);
    };

    def + (b:Quaternion):Quaternion = {
	val out = Quaternion.newOrCheckQuaternion(null, GUID, b.GUID, "Quaternion+".##)
	val odata = out.data;
	val bdata = b.data;
	odata(0) = data(0) + bdata(0);
	odata(1) = data(1) + bdata(1);
	odata(2) = data(2) + bdata(2);
	odata(3) = data(3) + bdata(3);
	out;
    };

    def - (b:Quaternion):Quaternion = {
	val out = Quaternion.newOrCheckQuaternion(null, GUID, b.GUID, "Quaternion-".##)
	val odata = out.data;
	val bdata = b.data;
	odata(0) = data(0) - bdata(0);
	odata(1) = data(1) - bdata(1);
	odata(2) = data(2) - bdata(2);
	odata(3) = data(3) - bdata(3);
	out;
    };

    def / (q:Quaternion):Quaternion = {
	this *^ q / q.normsq;
    }

    override def * (v:Float):Quaternion = {
	val out = Quaternion.newOrCheckQuaternion(null, GUID, v.hashCode, "Quaternion*Scalar".##);
	val odata = out.data;
	odata(0) = data(0)*v;
	odata(1) = data(1)*v;
	odata(2) = data(2)*v;
	odata(3) = data(3)*v;
	out;
    }

    override def / (v:Float):Quaternion = {
	val out = Quaternion.newOrCheckQuaternion(null, GUID, v.hashCode, "Quaternion/Scalar".##);
	val odata = out.data;
	odata(0) = data(0)/v;
	odata(1) = data(1)/v;
	odata(2) = data(2)/v;
	odata(3) = data(3)/v;
	out;
    }

    override def + (v:Float):Quaternion = {
	val out = Quaternion.newOrCheckQuaternion(null, GUID, v.hashCode, "Quaternion+Scalar".##);
	val odata = out.data;
	odata(0) = data(0)+v;
	odata(1) = data(1)+v;
	odata(2) = data(2)+v;
	odata(3) = data(3)+v;
	out;
    }

    override def - (v:Float):Quaternion = {
	val out = Quaternion.newOrCheckQuaternion(null, GUID, v.hashCode, "Quaternion-Scalar".##);
	val odata = out.data;
	odata(0) = data(0)-v;
	odata(1) = data(1)-v;
	odata(2) = data(2)-v;
	odata(3) = data(3)-v;
	out;
    }


    def toMat_(omat:FMat):FMat = {
	val out = FMat.newOrCheckFMat(4, 4, omat, GUID, "QuaternionToMat".##);
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

    def toMat:FMat = {
	toMat_(null);
    };

    def toMat(omat:FMat):FMat = {
	toMat_(omat);
    };

    def toRotMat_(omat:FMat):FMat = {
	val out = FMat.newOrCheckFMat(3, 3, omat, GUID, "QuaternionToRotMat".##);
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

    def toRotMat:FMat = {
	toRotMat_(null);
    };

    def toRotMat(omat:FMat):FMat = {
	toRotMat_(omat);
    };

    def normsq:Float = {
	data(0)*data(0) + data(1)*data(1) + data(2)*data(2) + data(3)*data(3);
    }

    def norm:Float = {
	math.sqrt(data(0)*data(0) + data(1)*data(1) + data(2)*data(2) + data(3)*data(3)).toFloat;
    }
}

object Quaternion {
    def apply():Quaternion = {
	new Quaternion(new Array[Float](4));
    }

    def apply(m:FMat):Quaternion = {
	val out = newOrCheckQuaternion(null, m.GUID, "QuaternionFromFMat".##);
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
 	    throw new RuntimeException("Quaternion constructor: FMat must be a 3 or 4 dimensional vector, or a 4x4 Matrix");
	}
	out
    }

    def fromRotMat(m:FMat) = {
	val out = newOrCheckQuaternion(null, m.GUID, "QuaternionFromRotMat".##);
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
	var q0 = 0f;
	var qx = 0f;
	var qy = 0f;
	var qz = 0f;
	if (tr > 0) { // q0^2 > 0.25;
	    q0 = math.sqrt(tr+1.0).toFloat * 0.5f; 
	    val scale = 0.25f / q0;
	    qx = (m21 - m12) * scale;
	    qy = (m02 - m20) * scale; 
	    qz = (m10 - m01) * scale; 
	} else if (m00 - m11 - m22 > 0) { // qx^2 > 0.25;
	    qx = math.sqrt(1.0 + m00 - m11 - m22).toFloat * 0.5f;
	    val scale = 0.25f / qx;
	    q0 = (m21 - m12) * scale;
	    qy = (m01 + m10) * scale; 
	    qz = (m02 + m20) * scale; 
	} else if (m11 - m00 - m22 > 0) { // qy^2 > 0.25;
	    qy = math.sqrt(1.0 + m11 - m00 - m22).toFloat * 0.5f;
	    val scale = 0.25f / qy;
	    q0 = (m02 - m20) * scale;
	    qx = (m01 + m10) * scale; 
	    qz = (m12 + m21) * scale; 
	} else {  // qz^2 is at least 0.25
	    qz = math.sqrt(1.0 + m22 - m00 - m11).toFloat * 0.5f;
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

    def apply(q0:Float, qx:Float, qy:Float, qz:Float) = {
	val out = new Quaternion(new Array[Float](4));
	out.data(0) = q0;
	out.data(1) = qx;
	out.data(2) = qy;
	out.data(3) = qz;
	out;
    }
	

    def identity():Quaternion = {
	val out = new Quaternion(new Array[Float](4));
	out.data(0) = 1;
	out.data(0) = 0;
	out.data(0) = 0;
	out.data(0) = 0;
	out;
    }

    def zero():Quaternion = {
	val out = new Quaternion(new Array[Float](4));
	out.data(0) = 0;
	out.data(0) = 0;
	out.data(0) = 0;
	out.data(0) = 0;
	out;
    }

    def newOrCheckQuaternion(outmat:Quaternion):Quaternion = {
	if (outmat.asInstanceOf[AnyRef] == null) {
	    Quaternion();
	} else {
	    outmat;
	}
    }
  
    def newOrCheckQuaternion(outmat:Quaternion, matGuid:Long, opHash:Int, forceCache:Boolean):Quaternion = {
	if (outmat.asInstanceOf[AnyRef] != null || (!Mat.useCache && !forceCache)) {
	    newOrCheckQuaternion(outmat);
	} else {
	    val key = (matGuid, opHash);
	    val res = Mat.cache2(key).asInstanceOf[Quaternion];
	    if (res != null) {
		newOrCheckQuaternion(res);
	    } else {
		val omat = newOrCheckQuaternion(null);
		Mat.cache2put(key, omat);
		omat
	    }
	}
    }

    def newOrCheckQuaternion(outmat:Quaternion,guid1:Long, opHash:Int):Quaternion =
	newOrCheckQuaternion(outmat, guid1, opHash, false);

  
    def newOrCheckQuaternion(outmat:Quaternion, guid1:Long, guid2:Long, opHash:Int, forceCache:Boolean):Quaternion = {
	if (outmat.asInstanceOf[AnyRef] != null || (!Mat.useCache && !forceCache)) {
	    newOrCheckQuaternion(outmat);
	} else {
	    val key = (guid1, guid2, opHash);
	    val res = Mat.cache3(key).asInstanceOf[Quaternion];
	    if (res != null) {
		newOrCheckQuaternion(res);
	    } else {
		val omat = newOrCheckQuaternion(null);
		Mat.cache3put(key, omat);
		omat;
	    }
	}
    }
  
    def newOrCheckQuaternion(outmat:Quaternion,guid1:Long, guid2:Long, opHash:Int):Quaternion =
	newOrCheckQuaternion(outmat, guid1, guid2, opHash, false);
  
    def newOrCheckQuaternion(outmat:Quaternion, guid1:Long, guid2:Long, guid3:Long, opHash:Int, forceCache:Boolean):Quaternion = {
	if (outmat.asInstanceOf[AnyRef] != null || (!Mat.useCache && !forceCache)) {
	    newOrCheckQuaternion(outmat);
	} else {
	    val key = (guid1, guid2, guid3, opHash);
	    val res = Mat.cache4(key).asInstanceOf[Quaternion];
	    if (res != null) {
		newOrCheckQuaternion(res);
	    } else {
		val omat = newOrCheckQuaternion(null);
		Mat.cache4put(key, omat);
		omat
	    }
	}
    }
  
    def newOrCheckQuaternion(outmat:Quaternion,guid1:Long, guid2:Long, guid3:Long, opHash:Int):Quaternion =
	newOrCheckQuaternion(outmat, guid1, guid2, guid3, opHash, false);

}
